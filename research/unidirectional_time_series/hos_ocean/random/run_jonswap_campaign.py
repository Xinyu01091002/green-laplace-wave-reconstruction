"""Hybrid initialization, short high preflight, low/high HOS, then GL records."""
import csv,hashlib,json,math,os,pathlib,re,shutil,signal,subprocess,time
root=pathlib.Path(__file__).resolve().parent
previous=pathlib.Path('/home/lxy/green-laplace-unidirectional-time-series-runs/hos-directional-akp002-20260926-v1')
env=os.environ.copy();env.update({k:'1' for k in ['OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS']})
state={'stage':'starting','outcomes':{},'started_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())}
def utc():return time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())
def status(stage,**kw):
 state.update(stage=stage,updated_utc=utc(),**kw);p=root/'pipeline-status.json';q=p.with_suffix('.tmp');q.write_text(json.dumps(state,indent=2));q.replace(p)
def rss_tree(pid):
 info={}
 for p in pathlib.Path('/proc').iterdir():
  if not p.name.isdigit():continue
  try:
   t=(p/'status').read_text();pp=int(re.search(r'^PPid:\s*(\d+)',t,re.M)[1]);v=re.search(r'^VmRSS:\s*(\d+)',t,re.M);info[int(p.name)]=(pp,int(v[1]) if v else 0)
  except (OSError,TypeError,ValueError):pass
 selected={pid}
 while True:
  add={p for p,(pp,_) in info.items() if pp in selected}-selected
  if not add:break
  selected|=add
 return sum(info[p][1] for p in selected if p in info)
def matlab(run,command,stem,budget=32):
 start=time.monotonic();peak=0
 with (run/(stem+'.log')).open('w') as f:
  p=subprocess.Popen(['/usr/bin/taskset','-c','40',str(previous/'bin/time'),'-v','-o',str(run/(stem+'-resources.txt')),'/home/lxy/Desktop/matlabr2026a/bin/matlab','-singleCompThread','-batch',"addpath('%s');%s"%(root,command)],cwd=run,env=env,stdout=f,stderr=subprocess.STDOUT,start_new_session=True)
  while p.poll() is None:
   rss=rss_tree(p.pid);peak=max(peak,rss);available=int(re.search(r'MemAvailable:\s*(\d+)',pathlib.Path('/proc/meminfo').read_text())[1])
   (run/(stem+'-monitor.json')).write_text(json.dumps(dict(pid=p.pid,wall_seconds=time.monotonic()-start,peak_RSS_KiB=peak,current_RSS_KiB=rss)))
   if rss>budget*1024*1024 or available<64*1024*1024:os.killpg(p.pid,signal.SIGTERM);raise RuntimeError('MATLAB own-job memory guard')
   time.sleep(3)
 assert p.returncode==0,(run/(stem+'.log')).read_text()[-3000:]
def setup(run):
 (run/'bin').mkdir()
 for name in ['HOS-Ocean','time']:shutil.copy2(previous/'bin'/name,run/'bin'/name)
 for name in ['environment.json','build-provenance.json','benchmark-validation.json']:shutil.copy2(previous/name,run/name)
 for name in ['run_directional_hos.py','export_jonswap_hos.m']:shutil.copy2(root/name,run/name)
 (run/'selection.json').write_text(json.dumps({'chosen_ranks':8,'reason':'Prior short MPI validation; new broadband high startup verified separately'}))
def hos(run):
 with (run/'controller.log').open('w') as f:p=subprocess.run(['python3',str(run/'run_directional_hos.py')],cwd=run,stdout=f,stderr=subprocess.STDOUT)
 assert p.returncode==0,(run/'controller.log').read_text()[-3000:]
 assert json.loads((run/'status.json').read_text())['state']=='completed'
def gl(run):
 out=run/'gl-time-comparison-v1';out.mkdir();s=json.loads((run/'settings.json').read_text());files={}
 for phase in [0,90,180,270]:
  lines=(run/'cases'/('phi%03d'%phase)/'Results/probes.dat').read_text().splitlines();i=next(i for i,x in enumerate(lines) if x.startswith('VARIABLES'));rows=[[float(v) for v in x.split()] for x in lines[i+1:] if x.strip()]
  assert len(rows)==s['expected_samples'] and all(len(row)==6 and all(math.isfinite(v) for v in row) for row in rows)
  p=out/('phi%03d.csv'%phase)
  with p.open('w',newline='') as f:csv.writer(f).writerows(rows)
  files[p.name]=hashlib.sha256(p.read_bytes()).hexdigest()
 (out/'snapshot.json').write_text(json.dumps(dict(sample_count=s['expected_samples'],end_time_s=s['duration_s'],partial=False,files_sha256=files,source_run=str(run)),indent=2))
 (out/'source').mkdir();subprocess.run(['tar','-xzf',str(previous/'directional-preview-source.tar.gz'),'-C',str(out/'source')],check=True)
 matlab(out,"compare_jonswap_time('%s');"%run,'run',128)
 return str(out)
try:
 assert not (root/'low').exists() and not (root/'high').exists(),'Never overwrite existing campaign'
 status('preparing_R4_GL_initial_conditions')
 matlab(root,"prepare_jonswap_hybrid('%s');"%root,'initialization',32)
 for name in ['low','high']:setup(root/name)
 short=root/'short-high';short.mkdir();shutil.copytree(root/'high/inputs',short/'inputs');shutil.copytree(root/'high/cases',short/'cases');setup(short)
 s=json.loads((root/'high/settings.json').read_text());s.update(duration_s=2.0,duration_Tp=2/s['Tp'],expected_samples=11)
 (short/'settings.json').write_text(json.dumps(s,indent=2))
 for p in short.glob('cases/*/input.yml'):
  t=p.read_text();t,n=re.subn(r'(?m)^  duration: .*$', '  duration: 2.0',t);assert n==1;p.write_text(t)
 status('short_high_2_seconds');hos(short)
 status('short_high_passed',short_high_status=str(short/'status.json'))
 for name in ['low','high']:
  status('HOS_running',active_family=name,next_family='high' if name=='low' else None)
  hos(root/name);state['outcomes'][name]={'HOS':'completed'}
 for name in ['low','high']:
  status('GL_processing',active_family=name)
  try:state['outcomes'][name].update(GL='completed',GL_result=gl(root/name))
  except Exception as exc:state['outcomes'][name].update(GL='failed',GL_error=str(exc))
 final='completed_HOS_and_GL' if all(x.get('GL')=='completed' for x in state['outcomes'].values()) else 'HOS_completed_GL_needs_attention'
 status(final,active_family=None,next_family=None)
except Exception as exc:
 status('failed',error=str(exc))
finally:
 (root/'summary.json').write_text(json.dumps(state,indent=2));(root/'status.txt').write_text(state['stage']+'\n');(root/'exit-code.txt').write_text(('0' if state['stage']=='completed_HOS_and_GL' else '1')+'\n');(root/'finished.utc').write_text(utc()+'\n')
 subprocess.run(['sh',str(root/'completion_email.sh'),str(root),'hos-jonswap-r4gl-low-high-80tp-20260926'],cwd=root)
