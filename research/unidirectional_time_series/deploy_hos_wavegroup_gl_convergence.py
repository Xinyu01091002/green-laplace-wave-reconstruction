"""Deploy focused-wave-group band, grid and GL-rank convergence checks."""
import hashlib
import json
import pathlib
import shlex
import subprocess
import tarfile

root=pathlib.Path(__file__).resolve().parents[2]
remote='/home/lxy/green-laplace-unidirectional-time-series-runs/hos-wavegroup-gl-convergence-20260928-v1'
local=root/'artifacts/unidirectional_time_series/hos-wavegroup-gl-convergence-v1'
local.mkdir(parents=True,exist_ok=False)
names=[
    'gl_directional_time_modal_grid.m',
    'gl_directional_time_modal_grid_eta33.m',
    'prepare_directional_joint_wavegroup.m',
    'run_hos_wavegroup_gl_convergence.m',
]
files={name:root/'research/unidirectional_time_series'/name for name in names}
files['allocate_directional_record.m']=root/'research/directional_wave_data/allocate_directional_record.m'
hashes={name:hashlib.sha256(path.read_bytes()).hexdigest() for name,path in files.items()}
(local/'source_hashes.json').write_text(json.dumps(hashes,indent=2))
with tarfile.open(local/'source.tar.gz','w:gz') as archive:
    for name,path in files.items():
        archive.add(path,arcname=name)

key=str(pathlib.Path.home()/'.ssh/id_ed25519_cursor')
ssh=['ssh','-p','60093','-i',key,'-o','BatchMode=yes','lxy@60.188.112.99']
subprocess.run(ssh+[f'test ! -e {remote} && mkdir {remote}'],check=True)
subprocess.run(['scp','-P','60093','-i',key,str(local/'source.tar.gz'),f'lxy@60.188.112.99:{remote}/'],check=True)
quoted=' '.join(shlex.quote(name) for name in files)
subprocess.run(ssh+[f'cd {remote} && tar -xzf source.tar.gz && sha256sum {quoted}'],check=True)
call="maxNumCompThreads(4);run_hos_wavegroup_gl_convergence(pwd);"
result=subprocess.run(ssh+[f'cd {remote} && /home/lxy/Desktop/matlabr2026a/bin/matlab -batch {shlex.quote(call)} > run.log 2>&1'])
subprocess.run(ssh+[f'tail -60 {remote}/run.log'])
if result.returncode:
    raise SystemExit(result.returncode)
for name in ['report.json','band_metrics.csv','grid_metrics.csv','rank_metrics.csv','convergence.png','run.log']:
    subprocess.run(['scp','-P','60093','-i',key,f'lxy@60.188.112.99:{remote}/{name}',str(local/name)],check=True)
