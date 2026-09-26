import hashlib,json,os,pathlib,subprocess,time
r=pathlib.Path(__file__).resolve().parent
while True:
 state=json.loads((r/'status.json').read_text())
 if state['stage']=='completed':break
 if state['stage']=='failed':raise RuntimeError(state)
 time.sleep(15)
dest=r/'spark_diagnostic';dest.mkdir()
manifest={}
for kind in ['eta','psi']:
 p=r/'spark'/('spark_'+kind+'20_r_series.m');text=p.read_text();needle='if outsideDefect>1e-8'
 assert text.count(needle)==1
 instrumentation="capturePath=getenv('R4_CAPTURE_PATH');if ~isempty(capturePath),save(capturePath,'%s20','defect','outsideDefect','repairAudit','-v7.3');end\n"%kind
 new=text.replace(needle,instrumentation+needle);(dest/p.name).write_text(new)
 manifest[p.name]={'original_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'instrumented_sha256':hashlib.sha256((dest/p.name).read_bytes()).hexdigest(),'change':'Save intermediate field before unchanged rejection only'}
(r/'diagnostic-source.json').write_text(json.dumps(manifest,indent=2))
env=os.environ.copy();env.update({k:'1' for k in ['OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS']})
cmd="addpath('%s');capture_r4_rejected('%s');"%(r,r)
start=time.time()
with (r/'diagnostic.log').open('w') as f:
 p=subprocess.run(['/usr/bin/taskset','-c','40','/home/lxy/green-laplace-unidirectional-time-series-runs/hos-directional-akp002-20260926-v1/bin/time','-v','-o',str(r/'diagnostic-resources.txt'),'/home/lxy/Desktop/matlabr2026a/bin/matlab','-singleCompThread','-batch',cmd],env=env,cwd=r,stdout=f,stderr=subprocess.STDOUT)
(r/'diagnostic-status.json').write_text(json.dumps({'exit_code':p.returncode,'wall_seconds':time.time()-start},indent=2))
