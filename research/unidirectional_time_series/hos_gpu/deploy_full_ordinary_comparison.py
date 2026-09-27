"""Add the missing simultaneous-start ordinary full-length timing baseline."""
import datetime,hashlib,json,pathlib,shlex,subprocess
s=pathlib.Path(__file__).resolve().parent;stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
root='/root/hos-ordinary-full-'+stamp;local=s.parents[2]/'artifacts/hos_gpu'/('ordinary-full-'+stamp);local.mkdir(parents=True)
base=(s/'run_mps_comparison.py').read_text();marker="try:\n assert subprocess.run(['pgrep'";assert base.count(marker)==1
code=base.split(marker)[0]+r'''
try:
 assert subprocess.run(['pgrep','-f','(^|/)nvidia-cuda-mps-(control|server)( |$)'],capture_output=True).returncode==1,'MPS must be off for ordinary baseline'
 mpsroot=pathlib.Path('/root/hos-mps-20260927T090455Z');mps=json.loads((mpsroot/'status.json').read_text());assert mps['state']=='completed'
 assert mps['solver_sha256']==state['solver_sha256']
 for n,h in json.loads((source/'initial-manifest.json').read_text()).items():assert hashlib.sha256((source/n).read_bytes()).hexdigest()==h
 baseline=batch('ordinary-full',settings['duration_s'],env)
 checks={}
 for phase in phases:
  q=subprocess.run([str(compare),'fields',str(mpsroot/'mps-full'/phase),str(r/'ordinary-full'/phase),'1e-10'],text=True,capture_output=True)
  checks[phase]={'exit':q.returncode,'output':q.stdout,'stderr':q.stderr}
 (r/'full-comparisons.json').write_text(json.dumps(checks,indent=2))
 mps_seconds=mps['batches']['mps-full']['batch_wall_seconds']
 state.update(state='completed' if all(x['exit']==0 for x in checks.values()) else 'completed_comparison_failed',mps_full_batch_seconds=mps_seconds,matched_full_speedup=baseline['batch_wall_seconds']/mps_seconds,matched_wait_reduction=1-mps_seconds/baseline['batch_wall_seconds']);save()
except Exception:
 state.update(state='failed',error=traceback.format_exc());save()
 for p in live:
  if p.poll() is None:p.terminate()
 for p in live:p.wait()
 raise
finally:nv.nvmlShutdown()
'''
script=local/'run_ordinary_full.py';script.write_text(code);compile(code,str(script),'exec')
key=(pathlib.Path.home()/'.ssh/id_ed25519_cursor').as_posix();opts=['-i',key,'-o','BatchMode=yes','-o',f'ProxyCommand=ssh -p 60093 -i {key} -o BatchMode=yes -W %h:%p root@60.188.112.99'];ssh=['ssh',*opts,'root@192.168.2.92']
def run(cmd):
 p=subprocess.run(cmd,text=True,capture_output=True,timeout=120);print(p.stdout+p.stderr,flush=True);p.check_returncode();return p.stdout
run([*ssh,'test ! -e '+root+' && mkdir '+root]);run(['scp',*opts,str(script),'root@192.168.2.92:'+root+'/'])
digest=hashlib.sha256(script.read_bytes()).hexdigest()
launch="import hashlib,pathlib,subprocess; r=pathlib.Path("+repr(root)+"); assert hashlib.sha256((r/'run_ordinary_full.py').read_bytes()).hexdigest()=="+repr(digest)+"; f=(r/'controller.log').open('x'); p=subprocess.Popen(['python3',str(r/'run_ordinary_full.py')],cwd=r,stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT,start_new_session=True); print(p.pid)"
pid=int(run([*ssh,'python3 -c '+shlex.quote(launch)]).strip());(local/'location.json').write_text(json.dumps({'remote':root,'controller_pid':pid,'script_sha256':digest},indent=2));print(root,flush=True)
