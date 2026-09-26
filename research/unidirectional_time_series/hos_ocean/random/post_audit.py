import os,pathlib,subprocess
r=pathlib.Path(__file__).resolve().parent
out=r/'record-direction-audit-v2'
assert (out/'audit.mat').exists()
env=os.environ.copy();env.update({k:'1' for k in ['OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS']})
cmd="addpath('%s');summarize_audit('%s');"%(r,out)
with (out/'summary.log').open('w') as f:
 subprocess.run(['/usr/bin/taskset','-c','40','/home/lxy/Desktop/matlabr2026a/bin/matlab','-singleCompThread','-batch',cmd],env=env,cwd=r,stdout=f,stderr=subprocess.STDOUT,check=True)
