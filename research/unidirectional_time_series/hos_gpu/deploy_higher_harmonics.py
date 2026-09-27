"""Deploy an isolated private scientific pilot; keep full scientific data remote."""
import datetime,hashlib,json,pathlib,shlex,subprocess,tarfile
s=pathlib.Path(__file__).resolve().parent;repo=s.parents[2];inputs=repo/'artifacts/hos_gpu/higher-harmonics-v1'
stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ');local=repo/'artifacts/hos_gpu'/('higher-'+stamp);local.mkdir(parents=True)
remote='/home/lxy/green-laplace-unidirectional-time-series-runs/gpu-higher-harmonics-'+stamp
files=list(inputs.glob('*.m'))+list(inputs.glob('*.wl'))+[s/'run_higher_controller.py',repo/'research/unidirectional_time_series/gl_laplace_terms.m']
manifest={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in files};(local/'manifest.json').write_text(json.dumps(manifest,indent=2))
with tarfile.open(local/'source.tar.gz','w:gz') as tar:
 for p in files+[local/'manifest.json']:tar.add(p,arcname=p.name)
key=str(pathlib.Path.home()/'.ssh/id_ed25519_cursor');opts=['-i',key,'-o','BatchMode=yes'];ssh=['ssh','-p','60093',*opts,'root@60.188.112.99']
def run(cmd):
 p=subprocess.run(cmd,text=True,capture_output=True,timeout=120);print(p.stdout+p.stderr,flush=True);p.check_returncode();return p.stdout
run([*ssh,'test ! -e '+remote+' && mkdir '+remote]);run(['scp','-P','60093',*opts,str(local/'source.tar.gz'),'root@60.188.112.99:'+remote+'/'])
code=r'''
import hashlib,json,os,pathlib,pwd,subprocess,tarfile
r=pathlib.Path(ROOT)
with tarfile.open(r/'source.tar.gz') as t:t.extractall(r)
for name,h in json.loads((r/'manifest.json').read_text()).items():assert hashlib.sha256((r/name).read_bytes()).hexdigest()==h
origins=[pathlib.Path('/home/lxy/green-laplace-unidirectional-time-series/src/internal/gl_no_stokes_eta33.m'),pathlib.Path('/home/lxy/green-laplace-unidirectional-time-series-runs/hos-jonswap-r4gl-80tp-20260926-v2/spark/spark_eta20_r_series.m')]
(r/'authority_hashes.json').write_text(json.dumps({str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in origins},indent=2))
user=pwd.getpwnam('lxy')
for p in [r,*r.iterdir()]:os.chown(p,user.pw_uid,user.pw_gid)
with (r/'controller.log').open('x') as f:p=subprocess.Popen(['runuser','-u','lxy','--','python3',str(r/'run_higher_controller.py')],cwd=r,stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT,start_new_session=True)
print(p.pid)
'''.replace('ROOT',repr(remote))
pid=int(run([*ssh,'python3 -c '+shlex.quote(code)]).strip());(local/'location.json').write_text(json.dumps({'remote':remote,'controller_pid':pid},indent=2));print(remote,flush=True)
