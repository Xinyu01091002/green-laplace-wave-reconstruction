"""Deploy current joint-modal GL processing for focused XY and random HOS."""
import hashlib,json,pathlib,shlex,subprocess,tarfile
root=pathlib.Path(__file__).resolve().parents[2];src=root/'research/unidirectional_time_series'
remote='/home/lxy/green-laplace-unidirectional-time-series-runs/hos-xy-random-joint-modal-20260928-v1'
local=root/'artifacts/unidirectional_time_series/hos-xy-random-joint-modal-v1';local.mkdir(parents=True,exist_ok=False)
names=['postprocess_hos_xy_and_random.m','gl_directional_time_modal_grid.m','gl_directional_time_modal_grid_eta20.m','gl_directional_time_modal_grid_eta33.m','prepare_directional_joint_wavegroup.m']
files={name:src/name for name in names};files['allocate_directional_record.m']=root/'research/directional_wave_data/allocate_directional_record.m'
hashes={name:hashlib.sha256(path.read_bytes()).hexdigest() for name,path in files.items()};(local/'source_hashes.json').write_text(json.dumps(hashes,indent=2))
with tarfile.open(local/'source.tar.gz','w:gz') as archive:
    for name,path in files.items():archive.add(path,arcname=name)
key=str(pathlib.Path.home()/'.ssh/id_ed25519_cursor');ssh=['ssh','-p','60093','-i',key,'-o','BatchMode=yes','lxy@60.188.112.99']
subprocess.run(ssh+[f'test ! -e {remote} && mkdir -p {remote}/focused {remote}/random'],check=True)
subprocess.run(['scp','-P','60093','-i',key,str(local/'source.tar.gz'),f'lxy@60.188.112.99:{remote}/'],check=True)
subprocess.run(ssh+[f'cd {remote} && tar -xzf source.tar.gz && sha256sum '+ ' '.join(shlex.quote(name) for name in files)],check=True)
call="maxNumCompThreads(4);postprocess_hos_xy_and_random(fullfile(pwd,'focused'),fullfile(pwd,'random'));"
result=subprocess.run(ssh+[f'cd {remote} && /home/lxy/Desktop/matlabr2026a/bin/matlab -batch {shlex.quote(call)} > run.log 2>&1'])
subprocess.run(ssh+[f'tail -80 {remote}/run.log']);
if result.returncode:raise SystemExit(result.returncode)
for family in ['focused','random']:
    target=local/family;target.mkdir()
    for name in ['report.json','metrics.csv','eta20.png','eta22.png','eta33.png']:
        subprocess.run(['scp','-P','60093','-i',key,f'lxy@60.188.112.99:{remote}/{family}/{name}',str(target/name)],check=True)
subprocess.run(['scp','-P','60093','-i',key,f'lxy@60.188.112.99:{remote}/run.log',str(local/'run.log')],check=True)
