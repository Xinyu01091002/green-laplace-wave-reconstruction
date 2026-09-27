"""One new kpHs/2=.06 family from frozen order-two input, then existing GL."""
import csv
import hashlib
import json
import math
import os
import pathlib
import re
import shutil
import signal
import subprocess
import time

ROOT = pathlib.Path(__file__).resolve().parent
SOURCE = pathlib.Path('/home/lxy/green-laplace-unidirectional-time-series-runs/hos-jonswap-r4gl-80tp-20260926-v2')
PRIOR = pathlib.Path('/home/lxy/green-laplace-unidirectional-time-series-runs/hos-directional-akp002-20260926-v1')
PINNED = {
    'high/inputs/initial_fields.mat': 'dd16b1676443e5254bfbeb42f4a94e9b0260afd9a36bdac33df68cb208a1916a',
    'high/settings.json': 'b75359acc1a818e0a6f64041297d240fc2b83f206ec9410d370c9c9b290dfde0',
    'low/inputs/initial_fields.mat': '4cc3bd3797782a9967849402204a1d5e81ba26bd0ececa11220a2783b8a03237',
    'compare_jonswap_time.m': '88023c5af9a97cd5529acfdc614bc6d4df0ef5b9392d21b338e7d62e4d237dc9',
    'run_directional_hos.py': 'fdb57ecbbb345fe98917fee3c48471c71c5fc3366a9bc135288108a3fb67c1ff',
    'export_jonswap_hos.m': '66a7bd5c5d26d9cef472ed874e7bd5b863e80949bc6cca2395e291189ce5c3f4',
}
ENV = os.environ.copy()
ENV.update({k: '1' for k in ['OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS']})
STATE = {'stage': 'starting', 'outcomes': {}, 'kpHs_over_2': .06, 'source_campaign': str(SOURCE)}


def utc():
    return time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def atomic(path, value):
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2))
    temporary.replace(path)


def status(stage, **extra):
    STATE.update(stage=stage, updated_utc=utc(), **extra)
    atomic(ROOT / 'pipeline-status.json', STATE)


def copy_once(source, destination):
    assert source.is_file() and not destination.exists(), (source, destination)
    shutil.copy2(source, destination)
    assert sha(source) == sha(destination)


def rss_tree(pid):
    info = {}
    for path in pathlib.Path('/proc').iterdir():
        if not path.name.isdigit():
            continue
        try:
            text = (path / 'status').read_text()
            parent = int(re.search(r'^PPid:\s*(\d+)', text, re.M)[1])
            rss = re.search(r'^VmRSS:\s*(\d+)', text, re.M)
            info[int(path.name)] = (parent, int(rss[1]) if rss else 0)
        except (OSError, TypeError, ValueError):
            pass
    selected = {pid}
    while extra := {p for p, (parent, _) in info.items() if parent in selected} - selected:
        selected |= extra
    return sum(info[p][1] for p in selected if p in info)


def matlab(folder, command, stem, budget):
    peak = 0
    start = time.monotonic()
    with (folder / (stem + '.log')).open('x') as log:
        process = subprocess.Popen([
            '/usr/bin/taskset', '-c', '40', str(PRIOR / 'bin/time'), '-v', '-o',
            str(folder / (stem + '-resources.txt')),
            '/home/lxy/Desktop/matlabr2026a/bin/matlab', '-singleCompThread', '-batch',
            "addpath('%s');%s" % (ROOT, command),
        ], cwd=folder, env=ENV, stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
        try:
            while process.poll() is None:
                rss = rss_tree(process.pid)
                peak = max(peak, rss)
                available = int(re.search(r'MemAvailable:\s*(\d+)', pathlib.Path('/proc/meminfo').read_text())[1])
                atomic(folder / (stem + '-monitor.json'), dict(pid=process.pid,
                       wall_seconds=time.monotonic() - start, peak_RSS_KiB=peak, current_RSS_KiB=rss))
                assert rss < budget * 1024**2 and available > 64 * 1024**2, 'MATLAB own-job memory guard'
                time.sleep(3)
            assert process.returncode == 0, (folder / (stem + '.log')).read_text()[-3000:]
        finally:
            if process.poll() is None:
                os.killpg(process.pid, signal.SIGTERM)
                try:
                    process.wait(timeout=15)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGKILL)
                    process.wait()


def prepare_runtime():
    for name, expected in PINNED.items():
        assert sha(SOURCE / name) == expected, 'Frozen source mismatch: ' + name
    manifest = json.loads((ROOT / 'launch-source-manifest.json').read_text())
    for name, expected in manifest['files_sha256'].items():
        assert sha(ROOT / name) == expected, 'Launch source mismatch: ' + name
    dependencies = dict(PINNED)
    for name in ['run_directional_hos.py', 'export_jonswap_hos.m', 'compare_jonswap_time.m']:
        copy_once(SOURCE / name, ROOT / name)
    copy_once(SOURCE / 'completion_email.sh', ROOT / 'completion_email.sh')
    mail = (ROOT / 'completion_email.sh').read_text()
    replacements = {
        'JONSWAP gamma3.3 low then high kpHs/2=.02/.12: four phases each, approx80Tp, R4-GL initialization.':
        'JONSWAP gamma3.3 kpHs/2=.06: four phases, approx80Tp, same R4-GL input with linear x1/2 and quadratic x1/4.',
        'Same99-percent-energy support and random phases; Hs normalized separately for both levels.':
        'Same11128 parents and random phases as the .12 source; no raw-data download.',
        'Completion means both HOS families and GL time-series processing completed; inspect summary.json for failures.':
        'Completion means this four-phase HOS family and GL time-series processing completed; inspect summary.json for failures.',
        'low/gl-time-comparison-v1 and high/gl-time-comparison-v1': 'medium/gl-time-comparison-v1',
    }
    for old, new in replacements.items():
        assert mail.count(old) == 1
        mail = mail.replace(old, new)
    (ROOT / 'completion_email.sh').write_text(mail, newline='\n')
    subprocess.run(['sh', '-n', str(ROOT / 'completion_email.sh')], check=True)
    # Validate existing notification configuration without printing its contents.
    config = pathlib.Path('/usr/local/bin/run.sh').read_text()
    assert re.search(r'^EMAIL="\$\{EMAIL:-[^}]+\}"$', config, re.M)
    assert pathlib.Path('/usr/bin/mail').is_file()
    for name in ['initialization-audit.json', 'source-manifest.json', 'postprocess-runtime-patch.json',
                 'postprocess-rank-manifest.json', 'postprocess-display-manifest.json', 'completion_email.sh']:
        dependencies[name] = sha(SOURCE / name)
    for path in [PRIOR / 'bin/time', PRIOR / 'directional-preview-source.tar.gz']:
        dependencies[str(path)] = sha(path)
    atomic(ROOT / 'inherited-source-manifest.json', dict(source_campaign=str(SOURCE),
           dependencies_sha256=dependencies,
           deployed_sha256={name: sha(ROOT / name) for name in [
               'run_directional_hos.py', 'export_jonswap_hos.m', 'compare_jonswap_time.m', 'completion_email.sh']},
           note='Only email wording adapted; HOS binary, controller and GL postprocessor unchanged.'))


def prepare_family(run):
    (run / 'bin').mkdir()
    for name in ['HOS-Ocean', 'time']:
        copy_once(SOURCE / 'high/bin' / name, run / 'bin' / name)
    for name in ['environment.json', 'build-provenance.json', 'benchmark-validation.json']:
        copy_once(SOURCE / 'high' / name, run / name)
    for name in ['run_directional_hos.py', 'export_jonswap_hos.m']:
        copy_once(ROOT / name, run / name)
    atomic(run / 'selection.json', dict(chosen_ranks=8,
           reason='Same validated binary and MPI export; amplitude-only rerun, no new MPI benchmark.'))


def compare(run):
    out = run / 'gl-time-comparison-v1'
    out.mkdir()
    settings = json.loads((run / 'settings.json').read_text())
    hashes = {}
    for phase in [0, 90, 180, 270]:
        lines = (run / 'cases' / ('phi%03d' % phase) / 'Results/probes.dat').read_text().splitlines()
        index = next(i for i, line in enumerate(lines) if line.startswith('VARIABLES'))
        rows = [[float(v) for v in line.split()] for line in lines[index + 1:] if line.strip()]
        assert len(rows) == settings['expected_samples']
        assert all(len(row) == 6 and all(math.isfinite(v) for v in row) for row in rows)
        assert max(abs(row[0] - settings['output_dt_s'] * i) for i, row in enumerate(rows)) < 1e-8
        path = out / ('phi%03d.csv' % phase)
        with path.open('x', newline='') as stream:
            csv.writer(stream).writerows(rows)
        hashes[path.name] = sha(path)
    atomic(out / 'snapshot.json', dict(sample_count=settings['expected_samples'],
           end_time_s=settings['duration_s'], partial=False, files_sha256=hashes, source_run=str(run)))
    (out / 'source').mkdir()
    archive = PRIOR / 'directional-preview-source.tar.gz'
    inherited = json.loads((ROOT / 'inherited-source-manifest.json').read_text())
    assert sha(archive) == inherited['dependencies_sha256'][str(archive)]
    subprocess.run(['tar', '-xzf', str(archive), '-C', str(out / 'source')], check=True)
    matlab(out, "compare_jonswap_time('%s');" % run, 'run', 128)
    return str(out)


def main():
    # A second invocation must not alter status, results or send another email.
    assert not (ROOT / 'started.utc').exists() and not (ROOT / 'medium').exists()
    (ROOT / 'started.utc').write_text(utc() + '\n')
    run = ROOT / 'medium'
    try:
        status('preparing_source')
        prepare_runtime()
        status('preparing_order_consistent_initial_conditions')
        matlab(ROOT, "issues=checkcode('%s/prepare_jonswap_rescaled.m');for j=1:numel(issues),fprintf('Line %%d: %%s\\n',issues(j).line,issues(j).message);end;assert(isempty(issues));prepare_jonswap_rescaled('%s','%s',.06);" % (ROOT, ROOT, SOURCE), 'initialization', 32)
        prepare_family(run)
        status('HOS_running', active_family='medium', continuous_run=True, past_10s=False)
        with (run / 'controller.log').open('x') as log:
            process = subprocess.Popen(['python3', str(run / 'run_directional_hos.py')],
                                       cwd=run, stdout=log, stderr=subprocess.STDOUT)
            while process.poll() is None:
                if (run / 'status.json').exists() and not STATE['past_10s']:
                    current = json.loads((run / 'status.json').read_text())
                    cases = current.get('cases', {})
                    if len(cases) == 4 and all(v.get('last_logged_time_s', 0) >= 10 for v in cases.values()):
                        status('HOS_running', past_10s=True, passed_10s_utc=utc())
                time.sleep(10)
            assert process.returncode == 0, (run / 'controller.log').read_text()[-3000:]
        assert json.loads((run / 'status.json').read_text())['state'] == 'completed'
        STATE['outcomes']['medium'] = {'HOS': 'completed'}
        status('GL_processing')
        result = compare(run)
        STATE['outcomes']['medium'].update(GL='completed', GL_result=result)
        status('completed_HOS_and_GL', active_family=None)
    except Exception as error:
        if STATE['stage'] == 'GL_processing':
            STATE['outcomes']['medium'].update(GL='failed', GL_error=str(error))
        status('failed', error=str(error))
    finally:
        atomic(ROOT / 'summary.json', STATE)
        (ROOT / 'status.txt').write_text(STATE['stage'] + '\n')
        (ROOT / 'exit-code.txt').write_text(('0' if STATE['stage'] == 'completed_HOS_and_GL' else '1') + '\n')
        (ROOT / 'finished.utc').write_text(utc() + '\n')
        if (ROOT / 'completion_email.sh').is_file():
            subprocess.run(['sh', str(ROOT / 'completion_email.sh'), str(ROOT), ROOT.name], cwd=ROOT)


if __name__ == '__main__':
    main()
