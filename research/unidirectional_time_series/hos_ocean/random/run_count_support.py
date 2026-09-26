import pathlib,os,subprocess,time,hashlib,json
r=pathlib.Path(__file__).resolve().parent
p=r/'count_energy_support.m';(r/'jonswap-energy-count-source.json').write_text(json.dumps({'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'scope':'Design and linear preview only, no solver launch'}))
env=os.environ.copy();env.update({k:'1' for k in ['OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS']})
start=time.time()
with (r/'jonswap-energy-count.log').open('w') as f:
 p=subprocess.run(['/usr/bin/taskset','-c','40',str(r/'akp012/bin/time'),'-v','-o',str(r/'jonswap-energy-count-resources.txt'),'/home/lxy/Desktop/matlabr2026a/bin/matlab','-singleCompThread','-batch',"addpath('%s');count_energy_support('%s');"%(r,r)],cwd=r,env=env,stdout=f,stderr=subprocess.STDOUT)
(r/'jonswap-energy-count-status.json').write_text(json.dumps({'exit_code':p.returncode,'wall_seconds':time.time()-start}))