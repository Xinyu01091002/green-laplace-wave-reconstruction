import datetime,hashlib,json,pathlib,shlex,subprocess
s=pathlib.Path(__file__).resolve().parent;stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ');r='/root/hos-mps-'+stamp
local=s.parents[2]/'artifacts/hos_gpu'/('mps-'+stamp);local.mkdir(parents=True)
key=(pathlib.Path.home()/'.ssh/id_ed25519_cursor').as_posix();opts=['-i',key,'-o','BatchMode=yes','-o',f'ProxyCommand=ssh -p 60093 -i {key} -o BatchMode=yes -W %h:%p root@60.188.112.99'];ssh=['ssh',*opts,'root@192.168.2.92']
def run(cmd):
 p=subprocess.run(cmd,text=True,capture_output=True,timeout=120);print(p.stdout+p.stderr,flush=True);p.check_returncode();return p.stdout
run([*ssh,'test ! -e '+r+' && mkdir '+r]);run(['scp',*opts,str(s/'run_mps_comparison.py'),'root@192.168.2.92:'+r+'/'])
digest=hashlib.sha256((s/'run_mps_comparison.py').read_bytes()).hexdigest()
code="import hashlib,pathlib,subprocess; r=pathlib.Path("+repr(r)+"); assert hashlib.sha256((r/'run_mps_comparison.py').read_bytes()).hexdigest()=="+repr(digest)+"; f=(r/'controller.log').open('x'); p=subprocess.Popen(['python3',str(r/'run_mps_comparison.py')],cwd=r,stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT,start_new_session=True); print(p.pid)"
pid=int(run([*ssh,'python3 -c '+shlex.quote(code)]).strip());(local/'location.json').write_text(json.dumps({'remote':r,'controller_pid':pid,'script_sha256':digest},indent=2));print(r,flush=True)
