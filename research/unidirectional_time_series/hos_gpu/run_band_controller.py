import json,os,pathlib,re,signal,subprocess,time,traceback
r=pathlib.Path(__file__).resolve().parent;out=r/'first-band-check-v1';out.mkdir(exist_ok=False)
state={'state':'starting','started_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())}
def save():
 state['updated_utc']=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime());p=out/'status.tmp';p.write_text(json.dumps(state,indent=2));p.replace(out/'status.json')
save()
try:
 start=time.monotonic();peak=0;env=dict(os.environ,OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
 with (out/'matlab.log').open('x') as f:
  p=subprocess.Popen(['taskset','-c','41','/home/lxy/Desktop/matlabr2026a/bin/matlab','-singleCompThread','-batch',"addpath('%s');check_first_harmonic_band('%s');"%(r,r)],cwd=r,env=env,stdout=f,stderr=subprocess.STDOUT,start_new_session=True);state.update(state='running',pid=p.pid);save()
  while p.poll() is None:
   try:
    text=pathlib.Path('/proc',str(p.pid),'status').read_text();v=re.search(r'^VmRSS:\s*(\d+)',text,re.M);rss=int(v[1]) if v else 0;peak=max(peak,rss)
    if rss>32*1024*1024:os.killpg(p.pid,signal.SIGTERM);raise RuntimeError('Own band pilot memory guard')
   except FileNotFoundError:pass
   state.update(wall_seconds=time.monotonic()-start,peak_sampled_rss_kib=peak)
   try:state['progress']=json.loads((out/'progress.json').read_text())
   except (OSError,json.JSONDecodeError):pass
   save();time.sleep(10)
 assert p.returncode==0,'See matlab.log';assert json.loads((out/'report.json').read_text())['state']=='completed';state.update(state='completed',wall_seconds=time.monotonic()-start);save()
except Exception:state.update(state='failed',error=traceback.format_exc());save();raise
