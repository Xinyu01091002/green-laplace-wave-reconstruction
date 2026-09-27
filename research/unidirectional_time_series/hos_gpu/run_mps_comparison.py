"""Matched four-process A/B pilot, then an isolated MPS full-record replay."""
import csv,ctypes,hashlib,json,math,os,pathlib,subprocess,time,traceback
r=pathlib.Path(__file__).resolve().parent
source=pathlib.Path('/root/hos-low-80tp-20260927T052004Z');exe=pathlib.Path('/root/hos-long-20260927T040628Z/build/hos_run');compare=pathlib.Path('/root/hos-global-zero-audit-20260927T043227Z/hos_compare_runs')
phases=['phi000','phi090','phi180','phi270'];settings=json.loads((source/'settings.json').read_text())
env=dict(os.environ);env.pop('CUDA_MPS_PIPE_DIRECTORY',None);env.pop('CUDA_MPS_LOG_DIRECTORY',None)
pipe=r/'mps-pipe';logs=r/'mps-logs';pipe.mkdir();logs.mkdir()
mpsenv=dict(env,CUDA_MPS_PIPE_DIRECTORY=str(pipe),CUDA_MPS_LOG_DIRECTORY=str(logs))
state={'state':'starting','source_run':str(source),'precision':'fp64','tolerance':1e-12,'started_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'batches':{},'solver_sha256':hashlib.sha256(exe.read_bytes()).hexdigest()}
started=False;live=[]
nv=ctypes.CDLL('libnvidia-ml.so.1');assert nv.nvmlInit_v2()==0;device=ctypes.c_void_p();assert nv.nvmlDeviceGetHandleByIndex_v2(0,ctypes.byref(device))==0
class Memory(ctypes.Structure):_fields_=[('total',ctypes.c_ulonglong),('free',ctypes.c_ulonglong),('used',ctypes.c_ulonglong)]
def used():
 m=Memory();assert nv.nvmlDeviceGetMemoryInfo(device,ctypes.byref(m))==0;return m.used
def save():
 state['updated_utc']=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime());p=r/'status.tmp';p.write_text(json.dumps(state,indent=2));p.replace(r/'status.json')
def control(command):
 q=subprocess.run(['nvidia-cuda-mps-control'],input=command+'\n',text=True,capture_output=True,env=mpsenv,timeout=15);q.check_returncode();return q.stdout.strip()
def batch(name,duration,environment,mps=False):
 global live
 out=r/name;out.mkdir();handles=[];start=time.monotonic();proc={};peak=used();clients_verified=False;last=0
 for phase,cpu in zip(phases,[28,29,30,31]):
  f=(out/(phase+'.log')).open('x');handles.append(f)
  proc[phase]=subprocess.Popen(['taskset','-c',str(cpu),str(exe),'double',str(source/(phase+'.bin')),str(out/phase),str(duration),'1e-12'],stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT,env=environment)
 live=list(proc.values());record={'state':'running','duration_s':duration,'pids':{n:p.pid for n,p in proc.items()},'phases':{}};state['batches'][name]=record;state['state']=name;save()
 while any(p.poll() is None for p in proc.values()):
  peak=max(peak,used())
  if mps and not clients_verified and time.monotonic()-start>1:
   servers=control('get_server_list').split();connections={}
   for server in servers:
    if server.isdigit():connections[server]=control('get_client_list '+server).split()
   clients={int(v) for row in connections.values() for v in row if v.isdigit()}
   if set(record['pids'].values())<=clients:record['mps_connections']=connections;clients_verified=True
  if time.monotonic()-last>=5:
   for phase in phases:
    try:record['phases'][phase]=json.loads((out/phase/'status.json').read_text())
    except (OSError,json.JSONDecodeError):pass
   record.update(elapsed_wall_seconds=time.monotonic()-start,peak_device_bytes=peak);save();last=time.monotonic()
  if any(p.poll() not in [None,0] for p in proc.values()):raise RuntimeError(name+' solver failed')
  time.sleep(.1)
 wall=time.monotonic()-start
 for f in handles:f.close()
 if mps:assert clients_verified,'MPS clients were not verified; no MPS speed claim'
 summaries={}
 for phase,p in proc.items():
  assert p.returncode==0;summaries[phase]=json.loads((out/phase/'summary.json').read_text())
  with (out/phase/'diagnostics.csv').open() as f:rows=list(csv.DictReader(f))
  assert len(rows)==round(duration/.2)+1 and abs(float(rows[-1]['time'])-duration)<1e-7
  assert all(math.isfinite(float(v)) for row in rows for v in row.values())
  if duration==settings['duration_s']:
   dest=out/phase/'Results';dest.mkdir()
   with (dest/'probes.dat').open('w') as f:
    f.write('VARIABLES = t eta1 eta2 eta3 eta4 eta5\n')
    for row in rows:f.write(' '.join(row[k] for k in ['time','p1','p2','p3','p4','p5'])+'\n')
 record.update(state='completed',batch_wall_seconds=wall,peak_device_bytes=peak,summaries=summaries);save();live=[]
 return record
try:
 assert subprocess.run(['pgrep','-f','(^|/)nvidia-cuda-mps-(control|server)( |$)'],capture_output=True).returncode==1,'Existing MPS service; do not modify it'
 for n,h in json.loads((source/'initial-manifest.json').read_text()).items():assert hashlib.sha256((source/n).read_bytes()).hexdigest()==h
 baseline=batch('ordinary-10s',10,env)
 subprocess.run(['nvidia-cuda-mps-control','-d'],env=mpsenv,check=True);started=True
 candidate=batch('mps-10s',10,mpsenv,True)
 checks={}
 for phase in phases:
  q=subprocess.run([str(compare),'fields',str(r/'ordinary-10s'/phase),str(r/'mps-10s'/phase),'1e-12'],text=True,capture_output=True)
  checks[phase]={'exit':q.returncode,'output':q.stdout,'stderr':q.stderr};(r/'pilot-comparisons.json').write_text(json.dumps(checks,indent=2));assert q.returncode==0
  for count in ['accepted','rejected']:assert baseline['summaries'][phase][count]==candidate['summaries'][phase][count]
 state['pilot_batch_speedup']=baseline['batch_wall_seconds']/candidate['batch_wall_seconds'];save()
 batch('mps-full',settings['duration_s'],mpsenv,True)
 fullchecks={}
 for phase in phases:
  q=subprocess.run([str(compare),'fields',str(source/phase),str(r/'mps-full'/phase),'1e-10'],text=True,capture_output=True)
  fullchecks[phase]={'exit':q.returncode,'output':q.stdout,'stderr':q.stderr}
 (r/'full-comparisons.json').write_text(json.dumps(fullchecks,indent=2))
 state['state']='completed' if all(q['exit']==0 for q in fullchecks.values()) else 'completed_comparison_failed';save()
except Exception:
 state.update(state='failed',error=traceback.format_exc());save()
 for p in live:
  if p.poll() is None:p.terminate()
 for p in live:p.wait()
 raise
finally:
 if started:
  try:state['mps_shutdown']=control('quit')
  except Exception as e:state['mps_shutdown_error']=str(e)
  save()
 nv.nvmlShutdown()
