"""Resume reference transfer/checks without rerunning the already completed probe."""
import hashlib,json,pathlib,re,subprocess,sys
source=pathlib.Path(__file__).resolve().parent;repo=source.parents[2]
stamp=sys.argv[1];assert re.fullmatch(r'\d{8}T\d{6}Z',stamp)
adaptive='--adaptive' in sys.argv[2:]
filename='adaptive_case5.bin' if adaptive else 'ocean_case5.bin'
suffix='-adaptive' if adaptive else '-lan'
cpu='/home/lxy/green-laplace-unidirectional-time-series-runs/gpu-production-probe-'+stamp;gpu='/root/hos-production-check-'+stamp
local=repo/'artifacts/hos_gpu'/('production-'+stamp);assert local.is_dir()
key=pathlib.Path.home()/'.ssh/id_ed25519_cursor';common=['-i',str(key),'-o','BatchMode=yes','-o','ConnectTimeout=10']
proxy=f'ssh -p 60093 -i {key.as_posix()} -o BatchMode=yes -W %h:%p root@60.188.112.99'
ssh93=['ssh','-p','60093',*common,'root@60.188.112.99'];ssh92=['ssh',*common,'-o','ProxyCommand='+proxy,'root@192.168.2.92']
helper=source/'serve_reference_once.py'
subprocess.run(['scp','-P','60093',*common,str(helper),'root@60.188.112.99:'+cpu+'/serve_reference_once.py'],check=True)
expected={}
for family in ['low','high']:
    expected[family]=subprocess.check_output([*ssh93,'sha256sum '+cpu+'/jonswap_'+family+'_phi000/'+filename],text=True).split()[0]
    if adaptive:
        trace=local/(family+suffix+'.trace')
        subprocess.run(['scp','-P','60093',*common,'root@60.188.112.99:'+cpu+'/jonswap_'+family+'_phi000/adaptive_case5.trace',str(trace)],check=True)
        subprocess.run(['scp',*common,'-o','ProxyCommand='+proxy,str(trace),'root@192.168.2.92:'+gpu+'/'],check=True)
server=subprocess.Popen([*ssh93,'python3 '+cpu+'/serve_reference_once.py '+cpu+(' adaptive' if adaptive else '')],stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
info=json.loads(server.stdout.readline())
fetch='''import urllib.request,pathlib,hashlib,json
targets=CONFIG
for family,target in targets.items():
 p=pathlib.Path(target['destination']); digest=hashlib.sha256()
 with urllib.request.urlopen(target['url'],timeout=30) as response,p.open('xb') as out:
  while chunk:=response.read(1024*1024):out.write(chunk);digest.update(chunk)
 actual=digest.hexdigest();assert actual==target['sha256'],(family,actual)
 print(json.dumps({'family':family,'destination':str(p),'sha256':actual,'bytes':p.stat().st_size}),flush=True)
'''
targets={family:{'url':f"http://192.168.2.93:{info['port']}/{info['token']}/{family}",'destination':gpu+'/'+family+suffix+'.bin','sha256':expected[family]} for family in ['low','high']}
try:
    r=subprocess.run([*ssh92,'python3 -'],input=fetch.replace('CONFIG',repr(targets)),text=True,capture_output=True,timeout=180)
    print(r.stdout+r.stderr,flush=True);r.check_returncode();server.wait(timeout=30)
finally:
    if server.poll() is None:server.terminate()
transfers=[json.loads(line) for line in r.stdout.splitlines() if line.startswith('{')]
(local/'transfers.json').write_text(json.dumps(transfers,indent=2))
records=[]
for family in ['low','high']:
    tests=[('double',None),('float','1e-7'),('float','1e-8'),('float','1e-9')] if adaptive else [('double',None),('float',None)]
    for precision,tol in tests:
        executable='hos_adaptive_check' if adaptive else 'hos_ocean_check'
        command='timeout 300 '+gpu+'/build/'+executable+' '+precision+' '+gpu+'/'+family+suffix+'.bin'+(' '+tol if tol else '')
        r=subprocess.run([*ssh92,command],capture_output=True,text=True,timeout=330)
        row={'family':family,'phase':0,'precision':precision,'tolerance_override':tol,'exit':r.returncode,'fields':[json.loads(line) for line in r.stdout.splitlines() if line.startswith('{')],'stderr':r.stderr}
        records.append(row);(local/'results.json').write_text(json.dumps(records,indent=2))
        display=dict(row)
        if adaptive:display['fields']=[f for f in row['fields'] if f.get('type')=='summary']
        print(json.dumps(display),flush=True)
(local/'locations.json').write_text(json.dumps({'cpu':cpu,'gpu':gpu,'transfer_helper_sha256':hashlib.sha256(helper.read_bytes()).hexdigest(),'method':'temporary IP-restricted one-time LAN transfer; server exits after two files'},indent=2))
print('Compact evidence:',local,flush=True)
if not all(r['exit']==0 for r in records):raise SystemExit('Some precision checks failed; all results retained')
