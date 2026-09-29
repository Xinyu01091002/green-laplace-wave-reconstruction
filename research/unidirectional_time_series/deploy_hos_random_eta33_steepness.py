"""Deploy a three-seed unidirectional random-phase eta33 steepness campaign."""
import hashlib,json,pathlib,shlex,subprocess

repo=pathlib.Path(__file__).resolve().parents[2];research=repo/'research/unidirectional_time_series'
remote='/home/lxy/green-laplace-unidirectional-time-series-runs/hos-eta33-random-steepness-20260929-v1'
source='/home/lxy/green-laplace-unidirectional-time-series-runs/hos-gl-timeseries-20260926-v1'
key=str(pathlib.Path.home()/'.ssh/id_ed25519_cursor');ssh=['ssh','-p','60093','-i',key,'-o','BatchMode=yes','lxy@60.188.112.99']
local=repo/'artifacts/unidirectional_time_series/hos-eta33-random-steepness-20260929-v1';local.mkdir(parents=True,exist_ok=False)
files={
 'prepare_hos_random_eta33_steepness.m':research/'prepare_hos_random_eta33_steepness.m',
 'analyze_hos_random_eta33_steepness.m':research/'analyze_hos_random_eta33_steepness.m',
 'analyze_hos_random_eta33_amplitude_order.m':research/'analyze_hos_random_eta33_amplitude_order.m',
 'plot_hos_random_eta33_results.m':research/'plot_hos_random_eta33_results.m',
 'diagnose_random_hilbert_leakage.m':research/'diagnose_random_hilbert_leakage.m',
 'diagnose_random_tapered_harmonics.m':research/'diagnose_random_tapered_harmonics.m',
 'analyze_hos_random_eta33_tapered_gl.m':research/'analyze_hos_random_eta33_tapered_gl.m',
 'snapshot/repo/research/unidirectional_time_series/gl_unidirectional_time_series.m':research/'gl_unidirectional_time_series.m',
 'snapshot/repo/symbolic/generated/finite_depth_directional_order2_eta22_pure_gl8.json':repo/'symbolic/generated/finite_depth_directional_order2_eta22_pure_gl8.json',
 'snapshot/repo/symbolic/generated/finite_depth_directional_order3_nested_green_laplace.json':repo/'symbolic/generated/finite_depth_directional_order3_nested_green_laplace.json'}
manifest={name:hashlib.sha256(path.read_bytes()).hexdigest() for name,path in files.items()};(local/'manifest.json').write_text(json.dumps(manifest,indent=2))
controller=r'''import concurrent.futures,hashlib,json,os,pathlib,re,shutil,subprocess,time
r=pathlib.Path(__file__).parent;source=pathlib.Path(SOURCE);manifest=json.loads((r/'manifest.json').read_text())
for name,digest in manifest.items():assert hashlib.sha256((r/name).read_bytes()).hexdigest()==digest
solver=source/'build/sources/HOS-Ocean';timer=pathlib.Path('/home/lxy/green-laplace-unidirectional-time-series/artifacts/hos-ocean-v2.1.0/deps/usr/bin/time');matlab='/home/lxy/Desktop/matlabr2026a/bin/matlab'
state={'state':'started','started_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'runs':[],'solver_sha256':hashlib.sha256(solver.read_bytes()).hexdigest(),'max_concurrent_hos':24}
def save():(r/'status.json').write_text(json.dumps(state,indent=2))
env=os.environ.copy()
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:env[key]='1'
libs=pathlib.Path('/home/lxy/green-laplace-unidirectional-time-series/artifacts/hos-ocean-v2.1.0/deps/usr/lib/x86_64-linux-gnu');solver_env=env.copy();solver_env['LD_LIBRARY_PATH']=':'.join(str(libs/p) for p in ['', 'lapack','blas'])
def run(args,folder,name,environ):
 start=time.time()
 with (folder/(name+'.log')).open('x') as stream:result=subprocess.run([str(timer),'-v','-o',str(folder/(name+'-resources.txt')),*map(str,args)],cwd=folder,env=environ,stdout=stream,stderr=subprocess.STDOUT)
 item={'name':str(folder.relative_to(r))+'/'+name,'wall_seconds':time.time()-start,'exit_code':result.returncode}
 if result.returncode:raise RuntimeError(str(item)+'\n'+(folder/(name+'.log')).read_text()[-3000:])
 return item
try:
 state['state']='preparing';save();command="addpath('%s');prepare_hos_random_eta33_steepness('%s');"%(r,r)
 state['runs'].append(run([matlab,'-singleCompThread','-batch',command],r,'prepare',env));save()
 state['state']='smoke_akp018';save();case=r/'seed20260925/akp018/phase000';smoke=r/'smoke-akp018';shutil.copytree(case,smoke)
 text=(smoke/'input.yml').read_text();settings=json.loads((r/'seed20260925/akp018/initialization.json').read_text());text=re.sub(r'duration: [^\n]+','duration: '+str(settings['Tp']),text);(smoke/'input.yml').write_text(text)
 state['runs'].append(run([solver,'input.yml'],smoke,'hos',solver_env));save()
 tasks=[(seed,akp,phase) for seed in [20260925,20260926,20260927] for akp in range(2,19,2) for phase in [0,90,180,270]]
 state['state']='running_hos';save()
 def simulate(task):
  seed,akp,phase=task;folder=r/f'seed{seed}'/f'akp{akp:03d}'/f'phase{phase:03d}';return run([solver,'input.yml'],folder,'hos',solver_env)
 with concurrent.futures.ThreadPoolExecutor(max_workers=24) as pool:
  for item in pool.map(simulate,tasks):state['runs'].append(item);save()
 state['state']='analyzing';save();command="addpath('%s');analyze_hos_random_eta33_steepness('%s');"%(r/'snapshot/repo/research/unidirectional_time_series',r)
 state['runs'].append(run([matlab,'-singleCompThread','-batch',command],r,'analyze',env));save()
 state['state']='linear_hilbert_control';save();command="addpath('%s');diagnose_random_hilbert_leakage('%s');"%(r,r)
 state['runs'].append(run([matlab,'-singleCompThread','-batch',command],r,'linear-hilbert-control',env));save()
 state['state']='taper_diagnostic';save();command="addpath('%s');diagnose_random_tapered_harmonics('%s');"%(r,r)
 state['runs'].append(run([matlab,'-singleCompThread','-batch',command],r,'taper-diagnostic',env));save()
 state['state']='tapered_gl';save();command="addpath('%s');addpath('%s');analyze_hos_random_eta33_tapered_gl('%s');"%(r,r/'snapshot/repo/research/unidirectional_time_series',r)
 state['runs'].append(run([matlab,'-singleCompThread','-batch',command],r,'tapered-gl',env));state.update(state='completed',finished_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()));save()
except Exception as exc:state.update(state='failed',error=str(exc),finished_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()));save();raise
'''.replace('SOURCE',repr(source))
(local/'controller.py').write_text(controller)
def call(args,timeout=180):return subprocess.run(args,text=True,capture_output=True,check=True,timeout=timeout).stdout
assert 'TARGET_ABSENT' in call([*ssh,f"if test -e {shlex.quote(remote)}; then echo TARGET_EXISTS; else echo TARGET_ABSENT; fi"])
call([*ssh,f"mkdir -p {shlex.quote(remote)}/snapshot/repo/research/unidirectional_time_series {shlex.quote(remote)}/snapshot/repo/symbolic/generated && ln -s {shlex.quote(source)}/mf12 {shlex.quote(remote)}/mf12"])
for name,path in files.items():call(['scp','-P','60093','-i',key,'-o','BatchMode=yes',str(path),f'lxy@60.188.112.99:{remote}/{name}'])
call(['scp','-P','60093','-i',key,'-o','BatchMode=yes',str(local/'manifest.json'),str(local/'controller.py'),f'lxy@60.188.112.99:{remote}/'])
bootstrap=f"import pathlib,subprocess;r=pathlib.Path({remote!r});s=(r/'controller-launch.log').open('x');p=subprocess.Popen(['python3',str(r/'controller.py')],cwd=r,stdin=subprocess.DEVNULL,stdout=s,stderr=subprocess.STDOUT,start_new_session=True);print(p.pid)"
pid=int(call([*ssh,'python3 -c '+shlex.quote(bootstrap)]).strip());(local/'location.json').write_text(json.dumps({'remote':remote,'source':source,'pid':pid},indent=2));print(pid)
