import hashlib,json,os,pathlib,subprocess,time
r=pathlib.Path(__file__).resolve().parent
out=r/'record-direction-center-finer-v1'
assert not out.exists()
p=r/'audit_center_finer.m'
(r/'center-finer-source.json').write_text(json.dumps({'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'source':p.name,'scope':'existing records only; fixed refinements and time-window diagnostics'},indent=2))
env=os.environ.copy();env.update({k:'1' for k in ['OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS']})
cmd="addpath('%s');audit_center_finer('%s','%s');"%(r,r/'akp012',out)
start=time.time()
with (r/'center-finer.log').open('w') as f:
 p=subprocess.run(['/usr/bin/taskset','-c','40',str(r/'akp012/bin/time'),'-v','-o',str(r/'center-finer-resources.txt'),'/home/lxy/Desktop/matlabr2026a/bin/matlab','-singleCompThread','-batch',cmd],env=env,cwd=r,stdout=f,stderr=subprocess.STDOUT)
(r/'center-finer-status.json').write_text(json.dumps({'exit_code':p.returncode,'wall_seconds':time.time()-start,'output':str(out)},indent=2))
