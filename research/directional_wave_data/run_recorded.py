"""Run one bounded Linux job in its own process group; record resource evidence."""
import argparse
import datetime
import json
import os
import pathlib
import resource
import signal
import subprocess
import time

parser = argparse.ArgumentParser()
parser.add_argument('--log', required=True)
parser.add_argument('--report', required=True)
parser.add_argument('--timeout', required=True, type=float)
parser.add_argument('command', nargs=argparse.REMAINDER)
args = parser.parse_args()
command = args.command[1:] if args.command[:1] == ['--'] else args.command
assert command and args.timeout > 0
log = pathlib.Path(args.log)
report_path = pathlib.Path(args.report)
assert not log.exists() and not report_path.exists(), 'Refusing to overwrite a prior run.'
started = datetime.datetime.now(datetime.timezone.utc).isoformat()
start = time.monotonic()
peak_rss = peak_threads = 0
timed_out = False
with log.open('xb') as output:
    job = subprocess.Popen(command, stdout=output, stderr=subprocess.STDOUT,
                           stdin=subprocess.DEVNULL, start_new_session=True)
    while job.poll() is None:
        rss = threads = 0
        for entry in pathlib.Path('/proc').iterdir():
            if not entry.name.isdigit():
                continue
            try:
                if os.getpgid(int(entry.name)) != job.pid:
                    continue
                for line in (entry / 'status').read_text().splitlines():
                    if line.startswith('VmRSS:'):
                        rss += int(line.split()[1])
                    elif line.startswith('Threads:'):
                        threads += int(line.split()[1])
            except (FileNotFoundError, ProcessLookupError, PermissionError):
                continue
        peak_rss, peak_threads = max(peak_rss, rss), max(peak_threads, threads)
        if time.monotonic() - start > args.timeout:
            timed_out = True
            try:
                os.killpg(job.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
            try:
                job.wait(timeout=10)
            except subprocess.TimeoutExpired:
                try:
                    os.killpg(job.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
            break
        time.sleep(.25)
    code = job.wait()
usage = resource.getrusage(resource.RUSAGE_CHILDREN)
report = dict(command=command, started_utc=started, exit_code=code,
              timed_out=timed_out, elapsed_seconds=time.monotonic()-start,
              user_cpu_seconds=usage.ru_utime, system_cpu_seconds=usage.ru_stime,
              largest_child_peak_RSS_KiB=usage.ru_maxrss,
              sampled_process_group_peak_RSS_KiB=peak_rss,
              sampled_process_group_peak_threads=peak_threads,
              sampling_interval_seconds=.25,
              note='Group RSS sums can double-count shared pages; child peak is not a sum.')
with report_path.open('x') as output:
    json.dump(report, output, indent=2)
print(json.dumps(report, indent=2))
raise SystemExit(124 if timed_out else code)
