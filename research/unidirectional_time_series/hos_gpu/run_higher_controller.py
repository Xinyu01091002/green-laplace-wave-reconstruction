import json,os,pathlib,re,signal,subprocess,time,traceback
r=pathlib.Path(__file__).resolve().parent;state={'state':'algebra_checks','started_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())}
def save():
 state['updated_utc']=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime());p=r/'status.tmp';p.write_text(json.dumps(state,indent=2));p.replace(r/'status.json')
save()
try:
 with (r/'algebra.log').open('x') as f:q=subprocess.run(['taskset','-c','40','wolframscript','-file',str(r/'verify_directional_transfer.wl')],cwd=r,stdout=f,stderr=subprocess.STDOUT,timeout=120)
 assert q.returncode==0,'Algebra check failed or unavailable; see algebra.log'
 env=dict(os.environ,OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1');start=time.monotonic();peak=0
 with (r/'matlab.log').open('x') as f:
  p=subprocess.Popen(['taskset','-c','40','/home/lxy/Desktop/matlabr2026a/bin/matlab','-singleCompThread','-batch',"addpath('%s');run_higher_harmonics('%s');"%(r,r)],cwd=r,env=env,stdout=f,stderr=subprocess.STDOUT,start_new_session=True);state.update(state='running',pid=p.pid);save()
  while p.poll() is None:
   try:
    status=pathlib.Path('/proc',str(p.pid),'status').read_text();m=re.search(r'^VmRSS:\s*(\d+)',status,re.M);rss=int(m[1]) if m else 0;peak=max(peak,rss)
    if rss>64*1024*1024:os.killpg(p.pid,signal.SIGTERM);raise RuntimeError('Own pilot exceeded 64 GiB memory budget')
   except FileNotFoundError:pass
   state.update(wall_seconds=time.monotonic()-start,peak_sampled_rss_kib=peak)
   try:state['progress']=json.loads((r/'progress.json').read_text())
   except (OSError,json.JSONDecodeError):pass
   save();time.sleep(10)
 assert p.returncode==0,'Scientific pilot failed; see matlab.log'
 report=json.loads((r/'report.json').read_text());assert report['state']=='completed';state.update(state='completed',wall_seconds=time.monotonic()-start);save()
except Exception:state.update(state='failed',error=traceback.format_exc());save();raise
