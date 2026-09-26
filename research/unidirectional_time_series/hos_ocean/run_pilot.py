import pathlib, subprocess, os, json, hashlib, time
root=pathlib.Path('/home/lxy/green-laplace-unidirectional-time-series/artifacts/hos-ocean-v2.1.0')
case=root/'wavegroup-phase000-pilot'
source=pathlib.Path('/home/lxy/green-laplace-unidirectional-time-series-runs/ow3d16-kpd1-akp012-20260925T212822Z/inputs/wavegroup/initial_fields.mat')
env=os.environ.copy()
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']: env[k]='1'
cmd="addpath('%s');prepare_hos_pilot('%s','%s');"%(root,source,case)
with (root/'prepare-pilot.log').open('w') as f:
    p=subprocess.run(['/home/lxy/Desktop/matlabr2026a/bin/matlab','-singleCompThread','-batch',cmd],stdout=f,stderr=subprocess.STDOUT,env=env)
assert p.returncode==0,(root/'prepare-pilot.log').read_text()[-4000:]
libs=root/'deps/usr/lib/x86_64-linux-gnu'
env['LD_LIBRARY_PATH']=':'.join(str(libs/p) for p in ['', 'lapack','blas'])
record={'state':'running','start_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'input_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'binary_sha256':hashlib.sha256((root/'install/bin/HOS-Ocean').read_bytes()).hexdigest(),'scripts':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [root/'prepare_hos_pilot.m',root/'run_pilot.py']},'files':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [case/'input_HOS-Ocean.yml',case/'Results/3d_ini.dat']}}
(case/'status.json').write_text(json.dumps(record,indent=2))
start=time.time()
with (case/'run.log').open('w') as f:
    p=subprocess.run([str(root/'deps/usr/bin/time'),'-v','-o','resources.txt',str(root/'install/bin/HOS-Ocean'),'input_HOS-Ocean.yml'],cwd=case,env=env,stdout=f,stderr=subprocess.STDOUT)
record.update(state='solver_exited_pending_validation',exit_code=p.returncode,wall_seconds=time.time()-start)
(case/'status.json').write_text(json.dumps(record,indent=2))
print(json.dumps(record,indent=2));print((case/'resources.txt').read_text());print((case/'run.log').read_text()[-4000:])