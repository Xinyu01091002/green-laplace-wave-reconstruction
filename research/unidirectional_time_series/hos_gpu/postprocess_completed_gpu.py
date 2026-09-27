"""Prepare the existing MATLAB GL workflow on 93 from completed GPU probes."""
import datetime,hashlib,json,pathlib,shlex,subprocess,sys
s=pathlib.Path(__file__).resolve().parent;repo=s.parents[2]
stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
local=repo/'artifacts/hos_gpu'/('postprocess-'+stamp);local.mkdir(parents=True)
gpu='/root/hos-low-80tp-20260927T052004Z';remote='/home/lxy/green-laplace-unidirectional-time-series-runs/gpu-low-postprocess-'+stamp
key=(pathlib.Path.home()/'.ssh/id_ed25519_cursor').as_posix();opts=['-i',key,'-o','BatchMode=yes'];proxy=f'ssh -p 60093 -i {key} -o BatchMode=yes -W %h:%p root@60.188.112.99'
ssh92=['ssh',*opts,'-o','ProxyCommand='+proxy,'root@192.168.2.92'];ssh93=['ssh','-p','60093',*opts,'root@60.188.112.99']
def run(cmd):
    p=subprocess.run(cmd,text=True,capture_output=True,timeout=180)
    with (local/'deployment.log').open('a') as f:f.write(p.stdout+p.stderr)
    print(p.stdout+p.stderr,flush=True);p.check_returncode();return p.stdout
code=r'''
import csv,hashlib,http.server,json,pathlib,secrets,tarfile,threading,io
r=pathlib.Path(ROOT);state=json.loads((r/'parallel-status.json').read_text());assert state['state']=='completed'
parts={};hashes={}
for phase in [0,90,180,270]:
 name='phi%03d'%phase;p=r/name/'Results/probes.dat';assert hashlib.sha256(p.read_bytes()).hexdigest()==state['cases'][name]['probes_sha256']
 rows=[s.split() for s in p.read_text().splitlines()[1:] if s.strip()];assert len(rows)==5505 and abs(float(rows[-1][0])-1100.8)<1e-8
 data=(chr(10).join(','.join(row) for row in rows)+chr(10)).encode()
 parts[name+'.csv']=data;hashes[name+'.csv']=hashlib.sha256(data).hexdigest()
parts['snapshot.json']=json.dumps(dict(sample_count=5505,end_time_s=1100.8,partial=False,files_sha256=hashes,source_run=str(r),backend='GPU FP64 global-zero convention',completion=state),indent=2).encode()
buf=io.BytesIO()
with tarfile.open(fileobj=buf,mode='w:gz') as tar:
 for name,data in parts.items():
  info=tarfile.TarInfo(name);info.size=len(data);tar.addfile(info,io.BytesIO(data))
data=buf.getvalue();token=secrets.token_hex(24)
class Handler(http.server.BaseHTTPRequestHandler):
 def do_GET(self):
  if self.client_address[0]!='192.168.2.93' or self.path!='/'+token:self.send_error(403);return
  self.send_response(200);self.send_header('Content-Length',str(len(data)));self.end_headers();self.wfile.write(data);self.wfile.flush();threading.Thread(target=server.shutdown,daemon=True).start()
 def log_message(self,*args):pass
server=http.server.ThreadingHTTPServer(('192.168.2.92',0),Handler)
print(json.dumps({'url':'http://192.168.2.92:'+str(server.server_port)+'/'+token,'sha256':hashlib.sha256(data).hexdigest()}),flush=True)
timer=threading.Timer(180,server.shutdown);timer.daemon=True;timer.start()
try:server.serve_forever()
finally:timer.cancel();server.server_close()
'''.replace('ROOT',repr(gpu))
server=subprocess.Popen([*ssh92,'python3 -c '+shlex.quote(code)],text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
try:
 info=json.loads(server.stdout.readline())
 prep=r'''
import hashlib,json,pathlib,shutil,tarfile,urllib.request
r=pathlib.Path(DEST);assert not r.exists();r.mkdir();out=r/'gl-time-comparison-v1';out.mkdir()
data=urllib.request.urlopen(URL,timeout=120).read();assert hashlib.sha256(data).hexdigest()==DIGEST;(r/'gpu-probes.tar.gz').write_bytes(data)
with tarfile.open(r/'gpu-probes.tar.gz') as t:t.extractall(out)
snapshot=json.loads((out/'snapshot.json').read_text())
for name,h in snapshot['files_sha256'].items():assert hashlib.sha256((out/name).read_bytes()).hexdigest()==h
base=pathlib.Path('/home/lxy/green-laplace-unidirectional-time-series-runs/hos-jonswap-r4gl-80tp-20260926-v2')
shutil.copy2(base/'low/settings.json',r/'settings.json');(r/'inputs').symlink_to(base/'low/inputs',target_is_directory=True)
shutil.copy2(base/'compare_jonswap_time.m',r/'compare_jonswap_time.m')
(out/'source').mkdir()
archive=base.parent/'hos-directional-akp002-20260926-v1/directional-preview-source.tar.gz'
with tarfile.open(archive) as t:t.extractall(out/'source')
meta={'source_campaign':str(base),'probe_archive_sha256':DIGEST,'matlab_script_sha256':hashlib.sha256((r/'compare_jonswap_time.m').read_bytes()).hexdigest(),'source_archive_sha256':hashlib.sha256(archive.read_bytes()).hexdigest()}
(r/'provenance.json').write_text(json.dumps(meta,indent=2));print(json.dumps(meta))
'''.replace('DEST',repr(remote)).replace('URL',repr(info['url'])).replace('DIGEST',repr(info['sha256']))
 run([*ssh93,'python3 -c '+shlex.quote(prep)]);server.wait(timeout=30)
finally:
 if server.poll() is None:server.terminate()
run(['scp','-P','60093',*opts,str(s/'run_gpu_postprocess.py'),'root@60.188.112.99:'+remote+'/'])
launch=r'''
import os,pathlib,pwd,subprocess
r=pathlib.Path(DEST);user=pwd.getpwnam('lxy')
for base,dirs,files in os.walk(r,followlinks=False):
 os.chown(base,user.pw_uid,user.pw_gid)
 for name in dirs+files:os.chown(pathlib.Path(base)/name,user.pw_uid,user.pw_gid,follow_symlinks=False)
with (r/'controller.log').open('x') as f:
 p=subprocess.Popen(['runuser','-u','lxy','--','python3',str(r/'run_gpu_postprocess.py')],cwd=r,stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT,start_new_session=True)
print(p.pid)
'''.replace('DEST',repr(remote))
pid=int(run([*ssh93,'python3 -c '+shlex.quote(launch)]).strip())
(local/'location.json').write_text(json.dumps({'remote':remote,'pid':pid},indent=2));print(remote,flush=True)
