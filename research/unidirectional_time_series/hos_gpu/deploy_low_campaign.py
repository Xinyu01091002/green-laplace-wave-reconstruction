"""Prepare on 93, transfer through restricted LAN, and detach a fresh GPU queue."""
import datetime,hashlib,json,pathlib,shlex,subprocess,tarfile
s=pathlib.Path(__file__).resolve().parent;repo=s.parents[2]
stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
cpu='/home/lxy/green-laplace-unidirectional-time-series-runs/gpu-campaign-inputs-'+stamp;gpu='/root/hos-low-80tp-'+stamp
local=repo/'artifacts/hos_gpu'/('campaign-'+stamp);local.mkdir(parents=True)
key=(pathlib.Path.home()/'.ssh/id_ed25519_cursor').as_posix();opts=['-i',key,'-o','BatchMode=yes'];proxy=f'ssh -p 60093 -i {key} -o BatchMode=yes -W %h:%p root@60.188.112.99'
ssh93=['ssh','-p','60093',*opts,'root@60.188.112.99'];gpuopts=[*opts,'-o','ProxyCommand='+proxy];ssh92=['ssh',*gpuopts,'root@192.168.2.92']
def run(cmd,timeout=600):
    p=subprocess.run(cmd,text=True,capture_output=True,timeout=timeout)
    with (local/'deploy.log').open('a') as f:f.write(p.stdout+p.stderr)
    print(p.stdout+p.stderr,flush=True);p.check_returncode();return p.stdout
files=[s/n for n in ['build_ocean_probe.py','ocean_probe.f90','ocean_adaptive_probe.f90','prepare_campaign_initials.py','serve_reference_once.py']]
manifest={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in files+[s/'run_low_campaign.py']};(local/'source-manifest.json').write_text(json.dumps(manifest,indent=2))
bundle=local/'source.tar.gz'
with tarfile.open(bundle,'w:gz') as tar:
    for p in files:tar.add(p,arcname=p.name)
run([*ssh93,'test ! -e '+cpu+' && mkdir '+cpu]);run(['scp','-P','60093',*opts,str(bundle),'root@60.188.112.99:'+cpu+'/'])
run([*ssh93,'cd '+cpu+' && tar -xzf source.tar.gz && python3 prepare_campaign_initials.py'],1200)
expected=run([*ssh93,'sha256sum '+cpu+'/campaign-initials.tar.gz']).split()[0]
run([*ssh92,'test ! -e '+gpu+' && mkdir '+gpu])
server=subprocess.Popen([*ssh93,'python3 '+cpu+'/serve_reference_once.py '+cpu+' campaign'],text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
try:
    info=json.loads(server.stdout.readline())
    code="import urllib.request,hashlib,pathlib,tarfile; p=pathlib.Path("+repr(gpu)+"); data=urllib.request.urlopen("+repr('http://192.168.2.93:'+str(info['port'])+'/'+info['token']+'/initials')+",timeout=120).read(); assert hashlib.sha256(data).hexdigest()=="+repr(expected)+"; (p/'initials.tar.gz').write_bytes(data); tarfile.open(p/'initials.tar.gz').extractall(p); print('Initial archive verified')"
    run([*ssh92,'python3 -c '+shlex.quote(code)]);server.wait(timeout=30)
finally:
    if server.poll() is None:server.terminate()
run(['scp',*gpuopts,str(s/'run_low_campaign.py'),str(local/'source-manifest.json'),'root@192.168.2.92:'+gpu+'/'])
# Exact phase-zero modal bytes must match the previously verified official export.
code="import pathlib; p=pathlib.Path("+repr(gpu)+"); a=(p/'phi000.bin').read_bytes(); b=pathlib.Path('/root/hos-production-check-20260927T025257Z/low-adaptive.bin').read_bytes(); assert a[108:]==b[108:len(a)]; print('Phase-zero state bytes match previous reference')"
run([*ssh92,'python3 -c '+shlex.quote(code)])
code="import pathlib,subprocess,json; r=pathlib.Path("+repr(gpu)+"); f=(r/'controller.log').open('w'); p=subprocess.Popen(['python3',str(r/'run_low_campaign.py')],cwd=r,stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT,start_new_session=True); (r/'launcher.json').write_text(json.dumps({'pid':p.pid})); print(p.pid)"
pid=int(run([*ssh92,'python3 -c '+shlex.quote(code)]).strip())
location={'input_export_93':cpu,'gpu_run_92':gpu,'launcher_pid':pid,'initial_archive_sha256':expected};(local/'location.json').write_text(json.dumps(location,indent=2));print(json.dumps(location),flush=True)
