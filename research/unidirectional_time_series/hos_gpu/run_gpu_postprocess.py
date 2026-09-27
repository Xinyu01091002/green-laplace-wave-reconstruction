"""Run the frozen scientific MATLAB postprocessor with isolated status/logs."""
import json,os,pathlib,re,signal,subprocess,time,traceback
r=pathlib.Path(__file__).resolve().parent
state={'state':'starting','started_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())}
def save():
 state['updated_utc']=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime());p=r/'status.tmp';p.write_text(json.dumps(state,indent=2));p.replace(r/'status.json')
try:
 env=dict(os.environ,OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
 cmd=['taskset','-c','40','/home/lxy/Desktop/matlabr2026a/bin/matlab','-singleCompThread','-batch',"addpath('%s');compare_jonswap_time('%s');"%(r,r)]
 start=time.monotonic();peak=0
 with (r/'matlab.log').open('x') as f:
  p=subprocess.Popen(cmd,cwd=r,env=env,stdout=f,stderr=subprocess.STDOUT,start_new_session=True);state.update(state='running',pid=p.pid);save()
  while p.poll() is None:
   try:
    text=pathlib.Path('/proc',str(p.pid),'status').read_text();m=re.search(r'^VmRSS:\s*(\d+)',text,re.M);rss=int(m[1]) if m else 0
    available=int(re.search(r'MemAvailable:\s*(\d+)',pathlib.Path('/proc/meminfo').read_text())[1]);peak=max(peak,rss)
    if rss>128*1024*1024 or available<64*1024*1024:os.killpg(p.pid,signal.SIGTERM);raise RuntimeError('Own postprocessor memory guard')
   except FileNotFoundError:pass
   state.update(wall_seconds=time.monotonic()-start,peak_sampled_rss_kib=peak);save();time.sleep(10)
 assert p.returncode==0,'MATLAB failed; see matlab.log'
 out=r/'gl-time-comparison-v1';report=json.loads((out/'report.json').read_text());assert report['status']=='FULL_RECORD_COMPARISON_COMPLETED'
 assert all((out/n).is_file() for n in ['metrics.csv','full_comparison.png','comparison.mat'])
 state.update(state='completed',wall_seconds=time.monotonic()-start,result_directory=str(out));save()
except Exception:state.update(state='failed',error=traceback.format_exc());save();raise
