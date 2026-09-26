import os,pathlib,subprocess,hashlib,json,time
r=pathlib.Path(__file__).resolve().parent
p=r/'plot_random_space.m'
(r/'space-plot-source.json').write_text(json.dumps({'source':p.name,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}))
env=os.environ.copy();env.update({k:'1' for k in ['OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS']})
start=time.time()
with (r/'space-plot.log').open('w') as f:
 p=subprocess.run(['/usr/bin/taskset','-c','40',str(r/'akp012/bin/time'),'-v','-o',str(r/'space-plot-resources.txt'),'/home/lxy/Desktop/matlabr2026a/bin/matlab','-singleCompThread','-batch',"addpath('%s');plot_random_space('%s');"%(r,r)],env=env,cwd=r,stdout=f,stderr=subprocess.STDOUT)
(r/'space-plot-status.json').write_text(json.dumps({'exit_code':p.returncode,'wall_seconds':time.time()-start}))