"""Run an exact repeated-field doubled-domain HOS boundary check."""
import concurrent.futures,hashlib,json,os,pathlib,shlex,subprocess,time
repo=pathlib.Path(__file__).resolve().parents[2];research=repo/'research/unidirectional_time_series'
source='/home/lxy/green-laplace-unidirectional-time-series-runs/hos-eta33-random-steepness-20260929-v1';solver_source='/home/lxy/green-laplace-unidirectional-time-series-runs/hos-gl-timeseries-20260926-v1'
remote='/home/lxy/green-laplace-unidirectional-time-series-runs/hos-random-double-domain-20260929-v2';key=str(pathlib.Path.home()/'.ssh/id_ed25519_cursor');ssh=['ssh','-p','60093','-i',key,'-o','BatchMode=yes','lxy@60.188.112.99']
local=repo/'artifacts/unidirectional_time_series/hos-random-double-domain-20260929-v2';local.mkdir(parents=True,exist_ok=False)
files={'prepare_hos_random_double_domain.m':research/'prepare_hos_random_double_domain.m','analyze_hos_random_double_domain.m':research/'analyze_hos_random_double_domain.m'};manifest={n:hashlib.sha256(p.read_bytes()).hexdigest() for n,p in files.items()};(local/'manifest.json').write_text(json.dumps(manifest,indent=2))
controller=r'''import concurrent.futures,hashlib,json,os,pathlib,subprocess,time
r=pathlib.Path(__file__).parent;manifest=json.loads((r/'manifest.json').read_text());source=SOURCE;solver_source=SOLVER_SOURCE
for n,h in manifest.items():assert hashlib.sha256((r/n).read_bytes()).hexdigest()==h
solver=pathlib.Path(solver_source)/'build/sources/HOS-Ocean';timer=pathlib.Path('/home/lxy/green-laplace-unidirectional-time-series/artifacts/hos-ocean-v2.1.0/deps/usr/bin/time');matlab='/home/lxy/Desktop/matlabr2026a/bin/matlab';state={'state':'preparing','runs':[],'solver_sha256':hashlib.sha256(solver.read_bytes()).hexdigest()}
def save():(r/'status.json').write_text(json.dumps(state,indent=2))
env=os.environ.copy()
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:env[key]='1'
libs=pathlib.Path('/home/lxy/green-laplace-unidirectional-time-series/artifacts/hos-ocean-v2.1.0/deps/usr/lib/x86_64-linux-gnu');se=env.copy();se['LD_LIBRARY_PATH']=':'.join(str(libs/p) for p in ['', 'lapack','blas'])
def run(args,folder,name,environ):
 start=time.time()
 with (folder/(name+'.log')).open('x') as f:result=subprocess.run([str(timer),'-v','-o',str(folder/(name+'-resources.txt')),*map(str,args)],cwd=folder,env=environ,stdout=f,stderr=subprocess.STDOUT)
 item={'name':str(folder.relative_to(r))+'/'+name,'seconds':time.time()-start,'exit_code':result.returncode}
 if result.returncode:raise RuntimeError(str(item)+'\n'+(folder/(name+'.log')).read_text()[-3000:])
 return item
try:
 save();command="addpath('%s');prepare_hos_random_double_domain('%s','%s');"%(r,source,r/'cases');state['runs'].append(run([matlab,'-singleCompThread','-batch',command],r,'prepare',env));state['state']='running';save()
 def simulate(phase):return run([solver,'input.yml'],r/'cases'/f'phase{phase:03d}','hos',se)
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:state['runs']+=list(pool.map(simulate,[0,90,180,270]))
 state['state']='analyzing';save();command="addpath('%s');analyze_hos_random_double_domain('%s','%s');"%(r,source,r/'cases');state['runs'].append(run([matlab,'-singleCompThread','-batch',command],r,'analyze',env));state['state']='completed';save()
except Exception as exc:state.update(state='failed',error=str(exc));save();raise
'''.replace('SOLVER_SOURCE',repr(solver_source)).replace('SOURCE',repr(source));(local/'controller.py').write_text(controller)
def call(args):return subprocess.run(args,text=True,capture_output=True,check=True,timeout=180).stdout
assert 'ABSENT' in call([*ssh,f'if test -e {remote}; then echo EXISTS; else echo ABSENT; fi']);call([*ssh,f'mkdir -p {remote}'])
for n,p in files.items():call(['scp','-P','60093','-i',key,'-o','BatchMode=yes',str(p),f'lxy@60.188.112.99:{remote}/{n}'])
call(['scp','-P','60093','-i',key,'-o','BatchMode=yes',str(local/'manifest.json'),str(local/'controller.py'),f'lxy@60.188.112.99:{remote}/'])
bootstrap=f"import pathlib,subprocess;r=pathlib.Path({remote!r});s=(r/'launch.log').open('x');p=subprocess.Popen(['python3',str(r/'controller.py')],cwd=r,stdin=subprocess.DEVNULL,stdout=s,stderr=subprocess.STDOUT,start_new_session=True);print(p.pid)";pid=int(call([*ssh,'python3 -c '+shlex.quote(bootstrap)]).strip());(local/'location.json').write_text(json.dumps({'remote':remote,'pid':pid},indent=2));print(pid)
