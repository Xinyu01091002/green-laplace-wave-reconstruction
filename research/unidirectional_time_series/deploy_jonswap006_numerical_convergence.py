"""Deploy fixed-input grid and GL-rank audit for JONSWAP kpHs/2=.06."""
import hashlib,json,pathlib,shlex,subprocess,tarfile
root=pathlib.Path(__file__).resolve().parents[2];src=root/'research/unidirectional_time_series';remote='/home/lxy/green-laplace-unidirectional-time-series-runs/hos-jonswap006-numerical-convergence-20260928-v1';local=root/'artifacts/unidirectional_time_series/hos-jonswap006-numerical-convergence-v1';local.mkdir(parents=True,exist_ok=False)
names=['run_jonswap006_numerical_convergence.m','gl_directional_time_modal_grid.m','gl_directional_time_modal_grid_eta33.m','prepare_directional_joint_wavegroup.m'];files={n:src/n for n in names};files['allocate_directional_record.m']=root/'research/directional_wave_data/allocate_directional_record.m';(local/'source_hashes.json').write_text(json.dumps({n:hashlib.sha256(p.read_bytes()).hexdigest() for n,p in files.items()},indent=2))
with tarfile.open(local/'source.tar.gz','w:gz') as a:
    for n,p in files.items():a.add(p,arcname=n)
key=str(pathlib.Path.home()/'.ssh/id_ed25519_cursor');ssh=['ssh','-p','60093','-i',key,'-o','BatchMode=yes','lxy@60.188.112.99'];subprocess.run(ssh+[f'test ! -e {remote} && mkdir {remote}'],check=True);subprocess.run(['scp','-P','60093','-i',key,str(local/'source.tar.gz'),f'lxy@60.188.112.99:{remote}/'],check=True);subprocess.run(ssh+[f'cd {remote} && tar -xzf source.tar.gz && sha256sum '+' '.join(shlex.quote(n) for n in files)],check=True)
call="maxNumCompThreads(4);run_jonswap006_numerical_convergence(pwd);";result=subprocess.run(ssh+[f'cd {remote} && /home/lxy/Desktop/matlabr2026a/bin/matlab -batch {shlex.quote(call)} > run.log 2>&1']);subprocess.run(ssh+[f'tail -100 {remote}/run.log']);
if result.returncode:raise SystemExit(result.returncode)
for n in ['report.json','grid_metrics.csv','rank_metrics.csv','run.log']:subprocess.run(['scp','-P','60093','-i',key,f'lxy@60.188.112.99:{remote}/{n}',str(local/n)],check=True)
