import datetime,hashlib,json,pathlib,subprocess
source=pathlib.Path(__file__).resolve().parent;repo=source.parents[2];stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
remote='/home/lxy/green-laplace-unidirectional-time-series-runs/hos-benchmark-package-'+stamp;local=repo/'artifacts/hos_gpu'/('package-'+stamp);local.mkdir(parents=True)
key=pathlib.Path.home()/'.ssh/id_ed25519_cursor';opts=['-i',str(key),'-o','BatchMode=yes'];ssh=['ssh','-p','60093',*opts,'root@60.188.112.99']
files=[source/'prepare_official_benchmark.py',source/'serve_reference_once.py'];manifest=local/'SHA256SUMS';manifest.write_text(''.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+p.name+'\n' for p in files))
def run(args):
    r=subprocess.run(args,text=True,capture_output=True,timeout=300)
    with (local/'run.log').open('a') as f:f.write(r.stdout+r.stderr)
    print(r.stdout+r.stderr,flush=True);r.check_returncode()
run([*ssh,'test ! -e '+remote+' && mkdir '+remote]);run(['scp','-P','60093',*opts,*map(str,files),str(manifest),'root@60.188.112.99:'+remote+'/'])
run([*ssh,f'cd {remote} && sha256sum -c SHA256SUMS && python3 prepare_official_benchmark.py'])
run(['scp','-P','60093',*opts,'root@60.188.112.99:'+remote+'/package-provenance.json',str(local)])
(local/'location.json').write_text(json.dumps({'remote':remote}));print('Package root:',remote,flush=True)
