import hashlib,json,os,pathlib,re,shutil,signal,subprocess,time
r=pathlib.Path(__file__).resolve().parent
env=os.environ.copy();env.update({k:'1' for k in ['OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS']})
env['MF12_PROGRESS']='1'
def utc():return time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())
state={'stage':'starting','started_utc':utc(),'stages':{},'memory_guard_GiB':160,'cpu':40}
def status():
 p=r/'status.json';q=p.with_suffix('.tmp');q.write_text(json.dumps(state,indent=2));q.replace(p)
def memory(pid):
 info={}
 for p in pathlib.Path('/proc').iterdir():
  if not p.name.isdigit():continue
  try:
   t=(p/'status').read_text();pp=int(re.search(r'^PPid:\s*(\d+)',t,re.M)[1]);rss=re.search(r'^VmRSS:\s*(\d+)',t,re.M);info[int(p.name)]=(pp,int(rss[1]) if rss else 0)
  except (OSError,TypeError,ValueError):continue
 ids={pid}
 while True:
  add={p for p,(pp,_) in info.items() if pp in ids}-ids
  if not add:break
  ids|=add
 return sum(info[p][1] for p in ids if p in info)
try:
 assert not (r/'input.mat').exists(),'Do not overwrite a previous run'
 mf=pathlib.Path('/home/lxy/green-laplace-unidirectional-time-series-runs/hos-directional-akp002-20260926-v1/mf12')
 shutil.copytree(mf,r/'mf12')
 files=[*r.glob('*.py'),*r.glob('*.m'),*r.glob('spark/*.m'),*r.glob('mf12/*.m')]
 (r/'manifest.json').write_text(json.dumps({str(p.relative_to(r)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},indent=2))
 for stage in ['prepare','r4_eta','r4_psi','mf12','summary']:
  available=int(re.search(r'MemAvailable:\s*(\d+)',pathlib.Path('/proc/meminfo').read_text())[1]);assert available>192*1024*1024,'Insufficient headroom'
  state['stage']=stage;state['updated_utc']=utc();status();started=time.monotonic()
  cmd="addpath('%s');benchmark_r4_mf12('%s','%s');"%(r,r,stage)
  with (r/(stage+'.log')).open('w') as log:
   p=subprocess.Popen(['/usr/bin/taskset','-c','40','/home/lxy/green-laplace-unidirectional-time-series-runs/hos-directional-akp002-20260926-v1/bin/time','-v','-o',str(r/(stage+'-resources.txt')),'/home/lxy/Desktop/matlabr2026a/bin/matlab','-singleCompThread','-batch',cmd],cwd=r,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
   peak=0
   while p.poll() is None:
    rss=memory(p.pid);peak=max(peak,rss);available=int(re.search(r'MemAvailable:\s*(\d+)',pathlib.Path('/proc/meminfo').read_text())[1])
    state['current']={'stage':stage,'pid':p.pid,'wall_seconds':time.monotonic()-started,'RSS_KiB':rss,'peak_RSS_KiB':peak};status()
    if rss>160*1024*1024 or available<64*1024*1024:
     os.killpg(p.pid,signal.SIGTERM);raise RuntimeError('Own-job memory guard triggered')
    time.sleep(2)
  state['stages'][stage]={'exit_code':p.returncode,'wall_seconds':time.monotonic()-started,'peak_sampled_RSS_KiB':peak}
  assert p.returncode==0,(r/(stage+'.log')).read_text()[-2500:]
 state.update(stage='completed',finished_utc=utc());state.pop('current',None);status()
except Exception as exc:
 state.update(stage='failed',error=str(exc),updated_utc=utc());status();raise
