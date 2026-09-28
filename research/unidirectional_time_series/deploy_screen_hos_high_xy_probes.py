"""Deploy the Akp=.12 off-axis probe amplitude screen."""
import hashlib
import json
import pathlib
import shlex
import subprocess

root=pathlib.Path(__file__).resolve().parents[2]
source=root/'research/unidirectional_time_series/screen_hos_high_xy_probes.m'
remote='/home/lxy/green-laplace-unidirectional-time-series-runs/hos-high-xy-probe-screen-20260928-v1'
local=root/'artifacts/unidirectional_time_series/hos-high-xy-probe-screen-v1'
local.mkdir(parents=True,exist_ok=False)
digest=hashlib.sha256(source.read_bytes()).hexdigest()
(local/'source_hashes.json').write_text(json.dumps({source.name:digest},indent=2))
key=str(pathlib.Path.home()/'.ssh/id_ed25519_cursor')
ssh=['ssh','-p','60093','-i',key,'-o','BatchMode=yes','lxy@60.188.112.99']
subprocess.run(ssh+[f'test ! -e {remote} && mkdir {remote}'],check=True)
subprocess.run(['scp','-P','60093','-i',key,str(source),f'lxy@60.188.112.99:{remote}/'],check=True)
call="maxNumCompThreads(4);screen_hos_high_xy_probes(pwd);"
result=subprocess.run(ssh+[f'cd {remote} && /home/lxy/Desktop/matlabr2026a/bin/matlab -batch {shlex.quote(call)} > run.log 2>&1'])
subprocess.run(ssh+[f'tail -40 {remote}/run.log'])
if result.returncode:raise SystemExit(result.returncode)
for name in ['report.json','screen.csv','screen.png','run.log']:
    subprocess.run(['scp','-P','60093','-i',key,f'lxy@60.188.112.99:{remote}/{name}',str(local/name)],check=True)
