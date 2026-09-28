"""Prepare and launch the isolated Akp=.12 off-axis HOS replay."""
import hashlib
import json
import pathlib
import shlex
import subprocess

root=pathlib.Path(__file__).resolve().parents[2]
research=root/'research/unidirectional_time_series'
remote='/home/lxy/green-laplace-unidirectional-time-series-runs/hos-directional-high-xy-probes-20260928-v1'
source='/home/lxy/green-laplace-unidirectional-time-series-runs/hos-directional-fourphase-20260926-v2'
files=[research/'prepare_hos_high_xy_probe_run.m',research/'run_hos_high_xy_probe_replay.py']
hashes={path.name:hashlib.sha256(path.read_bytes()).hexdigest() for path in files}
local=root/'artifacts/unidirectional_time_series/hos-directional-high-xy-probes-v1'
local.mkdir(parents=True,exist_ok=False);(local/'source_hashes.json').write_text(json.dumps(hashes,indent=2))
key=str(pathlib.Path.home()/'.ssh/id_ed25519_cursor')
ssh=['ssh','-p','60093','-i',key,'-o','BatchMode=yes','lxy@60.188.112.99']
subprocess.run(ssh+[f'test ! -e {remote} && mkdir -p {remote}/bin'],check=True)
subprocess.run(['scp','-P','60093','-i',key,str(local/'source_hashes.json'),f'lxy@60.188.112.99:{remote}/'],check=True)
for path in files:
    subprocess.run(['scp','-P','60093','-i',key,str(path),f'lxy@60.188.112.99:{remote}/'],check=True)
subprocess.run(ssh+[f'cp {source}/bin/HOS-Ocean {source}/bin/time {remote}/bin/ && cp {source}/environment.json {remote}/'],check=True)
prepare=f"prepare_hos_high_xy_probe_run('{remote}','{source}');"
subprocess.run(ssh+[f'cd {remote} && /home/lxy/Desktop/matlabr2026a/bin/matlab -batch {shlex.quote(prepare)} > prepare.log 2>&1'],check=True)
result=subprocess.run(ssh+[f'cd {remote} && python3 run_hos_high_xy_probe_replay.py > controller.log 2>&1'])
subprocess.run(ssh+[f'tail -50 {remote}/controller.log; cat {remote}/status.json'])
if result.returncode:raise SystemExit(result.returncode)
for name in ['settings.json','status.json','prepare.log','controller.log']:
    subprocess.run(['scp','-P','60093','-i',key,f'lxy@60.188.112.99:{remote}/{name}',str(local/name)],check=True)
