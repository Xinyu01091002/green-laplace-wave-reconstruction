"""LAN-copy the frozen official CPU runtime and probe prefixes to 92."""
import hashlib,json,pathlib,re,subprocess,sys
source=pathlib.Path(__file__).resolve().parent;repo=source.parents[2];stamp=sys.argv[1];assert re.fullmatch(r'\d{8}T\d{6}Z',stamp)
cpu='/home/lxy/green-laplace-unidirectional-time-series-runs/hos-benchmark-package-'+stamp;dest='/root/hos-official-benchmark-'+stamp
local=repo/'artifacts/hos_gpu'/('package-'+stamp);key=pathlib.Path.home()/'.ssh/id_ed25519_cursor';common=['-i',str(key),'-o','BatchMode=yes']
proxy=f'ssh -p 60093 -i {key.as_posix()} -o BatchMode=yes -W %h:%p root@60.188.112.99'
ssh93=['ssh','-p','60093',*common,'root@60.188.112.99'];ssh92=['ssh',*common,'-o','ProxyCommand='+proxy,'root@192.168.2.92']
subprocess.run([*ssh92,'test ! -e '+dest+' && mkdir '+dest],check=True)
hashes={name:subprocess.check_output([*ssh93,'sha256sum '+cpu+'/'+name],text=True).split()[0] for name in ['cpu-runtime.tar.gz','reference-records.tar.gz']}
server=subprocess.Popen([*ssh93,'python3 '+cpu+'/serve_reference_once.py '+cpu+' benchmark'],stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
info=json.loads(server.stdout.readline())
config={key:{'url':f"http://192.168.2.93:{info['port']}/{info['token']}/{key}",'path':dest+'/'+name,'sha256':hashes[name]} for key,name in [('runtime','cpu-runtime.tar.gz'),('reference','reference-records.tar.gz')]}
script='''import urllib.request,pathlib,hashlib,json,subprocess
config=CONFIG
for key,spec in config.items():
 digest=hashlib.sha256();p=pathlib.Path(spec['path'])
 with urllib.request.urlopen(spec['url'],timeout=30) as src,p.open('xb') as dst:
  while chunk:=src.read(1024*1024):dst.write(chunk);digest.update(chunk)
 assert digest.hexdigest()==spec['sha256']
 subprocess.run(['tar','-xzf',str(p),'-C',str(p.parent)],check=True)
 print(json.dumps({'file':p.name,'bytes':p.stat().st_size,'sha256':digest.hexdigest()}),flush=True)
'''
try:
    r=subprocess.run([*ssh92,'python3 -'],input=script.replace('CONFIG',repr(config)),text=True,capture_output=True,timeout=180)
    print(r.stdout+r.stderr,flush=True);r.check_returncode();server.wait(timeout=30)
finally:
    if server.poll() is None:server.terminate()
(local/'deployment.json').write_text(json.dumps({'destination':dest,'archives':hashes},indent=2));print('Deployed:',dest,flush=True)
