"""Frozen-input OW3D batch controller. No physics, email transport or local fetch."""
import argparse
import concurrent.futures
import datetime
import hashlib
import json
import os
import pathlib
import re
import signal
import subprocess
import threading
import time


def utc():
    return datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')


def atomic(path, value):
    path = pathlib.Path(path)
    temporary = path.with_name(path.name + f'.{threading.get_ident()}.tmp')
    temporary.write_text(value)
    os.replace(temporary, path)


def sha(path):
    with pathlib.Path(path).open('rb') as data:
        return hashlib.file_digest(data, 'sha256').hexdigest()


def native_completion(case, directory, code):
    if code != 0:
        raise RuntimeError(f'OW3D exited {code}')
    if 'JOB IS COMPLETE' not in (directory/'ow3d.log').read_text(errors='replace'):
        raise RuntimeError('Missing native completion marker (exit 0 alone is insufficient)')
    with (directory/'OceanWave3D.end').open() as stream:
        stream.readline()
        values = stream.readline().replace('D', 'E').split()
    if abs(float(values[4])-case['expected_final_time_s']) > 1e-8:
        raise RuntimeError('Native final time does not match the frozen input')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('manifest')
    args = parser.parse_args()
    manifest_path = pathlib.Path(args.manifest).resolve()
    root = manifest_path.parent
    m = json.loads(manifest_path.read_text())
    assert len(m['cases']) == 16 and m['concurrency'] == 16
    assert not (root/'started.utc').exists(), 'Refusing to rerun an already-started batch.'
    assert sha(m['solver']) == m['solver_sha256']
    for relative, digest in m['source_hashes'].items():
        assert sha(root/'source'/relative) == digest, f'Source changed: {relative}'
    for case in m['cases']:
        directory = root/'cases'/case['case_id']
        for filename, digest in case['input_hashes'].items():
            assert sha(directory/filename) == digest
    mem = pathlib.Path('/proc/meminfo').read_text()
    available = int(re.search(r'^MemAvailable:\s+(\d+)', mem, re.M).group(1))*1024
    assert available >= m['minimum_available_memory_bytes'], 'Insufficient fresh free memory for 16 jobs'
    atomic(root/'started.utc', utc()+'\n')
    atomic(root/'controller.pid', str(os.getpid())+'\n')
    atomic(root/'status.txt', 'running\n')
    lock = threading.Lock()
    validation_lock = threading.Lock()
    stopping = threading.Event()
    processes = {}
    states = {c['case_id']:dict(state='queued', pid=None, peak_observed_RSS_KiB=0) for c in m['cases']}

    def update(identifier, **values):
        with lock:
            states[identifier].update(values)

    def stop(signum, frame):
        stopping.set()
        with lock:
            children = list(processes.values())
        for child in children:
            if child.poll() is None:
                try:
                    os.killpg(child.pid, signal.SIGTERM)
                except ProcessLookupError:
                    pass

    signal.signal(signal.SIGTERM, stop)
    signal.signal(signal.SIGINT, stop)

    def worker(case):
        identifier = case['case_id']
        directory = root/'cases'/identifier
        try:
            if stopping.is_set():
                raise RuntimeError('Batch was cancelled before launch')
            env = os.environ.copy()
            env.update(GFORTRAN_UNBUFFERED_ALL='y', OMP_NUM_THREADS='1',
                       OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
            with (directory/'ow3d.log').open('xb') as output:
                child = subprocess.Popen([m['solver']], cwd=directory, env=env,
                                         stdout=output, stderr=subprocess.STDOUT,
                                         stdin=subprocess.DEVNULL, start_new_session=True)
                with lock:
                    processes[identifier] = child
                atomic(directory/'started.utc', utc()+'\n')
                atomic(directory/'pid.txt', str(child.pid)+'\n')
                atomic(directory/'status.txt', 'running\n')
                update(identifier, state='running', pid=child.pid, started_utc=utc())
                print(f'{utc()} START {identifier} pid={child.pid}', flush=True)
                code = child.wait()
            atomic(directory/'native-exit-code.txt', str(code)+'\n')
            update(identifier, native_exit_code=code)
            if stopping.is_set():
                raise RuntimeError('Batch cancelled')
            native_completion(case, directory, code)
            update(identifier, state='awaiting_validation')
            atomic(directory/'status.txt', 'awaiting_validation\n')
            # Keep expensive MATLAB startup/field extraction serial while other OW3D jobs run.
            with validation_lock:
                if stopping.is_set():
                    raise RuntimeError('Batch cancelled before validation')
                update(identifier, state='validating')
                atomic(directory/'status.txt', 'validating\n')
                source = root/'source/research/directional_wave_data'
                quote = lambda x: str(x).replace("'", "''")
                command = f"addpath('{quote(source)}'); validate_ow3d_campaign_case('{quote(directory)}');"
                with (directory/'validation.log').open('xb') as output:
                    validator = subprocess.Popen([m['matlab'], '-batch', command], cwd=root,
                        env=env, stdout=output, stderr=subprocess.STDOUT,
                        stdin=subprocess.DEVNULL, start_new_session=True)
                    with lock:
                        processes[identifier] = validator
                    update(identifier, validator_pid=validator.pid)
                    validation_code = validator.wait()
                if validation_code != 0:
                    raise RuntimeError(f'MATLAB output validator exited {validation_code}')
            v = json.loads((directory/'processed/validation.json').read_text())
            assert v['status'] == 'RAW_OUTPUT_INTEGRITY_VERIFIED'
            update(identifier, state='completed_raw_verified',
                   psi_checkpoint_warning=v['psi_checkpoint_relative_tolerance_warning'])
            atomic(directory/'status.txt', 'completed_raw_verified\n')
            atomic(directory/'exit-code.txt', '0\n')
            print(f'{utc()} VERIFIED {identifier}', flush=True)
        except Exception as error:
            update(identifier, state='cancelled' if stopping.is_set() else 'failed', error=str(error))
            atomic(directory/'status.txt', ('cancelled' if stopping.is_set() else 'failed')+'\n')
            atomic(directory/'exit-code.txt', '1\n')
            print(f'{utc()} FAILED {identifier}: {error}', flush=True)
        finally:
            atomic(directory/'finished.utc', utc()+'\n')

    def snapshot():
        with lock:
            for identifier, child in processes.items():
                if child.poll() is None and states[identifier]['state'] == 'running':
                    try:
                        status = pathlib.Path(f'/proc/{child.pid}/status').read_text()
                        peak = re.search(r'^VmHWM:\s+(\d+)', status, re.M)
                        rss = re.search(r'^VmRSS:\s+(\d+)', status, re.M)
                        states[identifier]['rss_KiB'] = int(rss.group(1)) if rss else 0
                        if peak:
                            states[identifier]['peak_observed_RSS_KiB'] = max(
                                states[identifier]['peak_observed_RSS_KiB'], int(peak.group(1)))
                        kin = list((root/'cases'/identifier).glob('Kinematics*.bin'))
                        states[identifier]['kinematics_bytes'] = sum(p.stat().st_size for p in kin)
                    except (FileNotFoundError, ProcessLookupError):
                        pass
            data = dict(updated_utc=utc(), run_id=m['run_id'], cases={k:dict(v) for k,v in states.items()})
        atomic(root/'status.json', json.dumps(data, indent=2))
        return data

    with concurrent.futures.ThreadPoolExecutor(max_workers=16) as pool:
        futures = [pool.submit(worker, case) for case in m['cases']]
        while not all(f.done() for f in futures):
            snapshot()
            stopping.wait(10) if not stopping.is_set() else time.sleep(1)
        for future in futures:
            future.result()
    final = snapshot()
    completed = sum(v['state']=='completed_raw_verified' for v in states.values())
    warnings = sum(bool(v.get('psi_checkpoint_warning')) for v in states.values())
    code = 0 if completed == 16 else 1
    state = 'completed_raw_verified_16_of_16' if code == 0 else f'finished_with_failures_{completed}_of_16_verified'
    if stopping.is_set():
        state = f'cancelled_{completed}_of_16_verified'
    final.update(completed=completed, total=16, psi_checkpoint_warnings=warnings,
                 physical_accuracy_certified=False, automatic_fetch=False)
    atomic(root/'summary.json', json.dumps(final, indent=2))
    atomic(root/'status.txt', state+'\n')
    atomic(root/'exit-code.txt', str(code)+'\n')
    print(f'{utc()} FINAL {state}; psi checkpoint warnings={warnings}; physical accuracy not certified', flush=True)
    print(f'Remote data: {root}/cases; no local raw-data fetch.', flush=True)
    atomic(root/'finished.utc', utc()+'\n')
    return code


if __name__ == '__main__':
    raise SystemExit(main())
