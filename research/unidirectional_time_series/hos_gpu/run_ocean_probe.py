"""Deploy source-only official-module probe on 93; retrieve its small test files."""
import argparse,datetime,hashlib,pathlib,subprocess
parser=argparse.ArgumentParser();parser.add_argument('--cash-karp',action='store_true');parser.add_argument('--adaptive',action='store_true');options_cli=parser.parse_args()
source=pathlib.Path(__file__).resolve().parent;repo=source.parents[2]
stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
remote='/home/lxy/green-laplace-unidirectional-time-series-runs/gpu-rhs-probe-'+stamp
local=repo/'artifacts/hos_gpu'/('ocean-'+stamp);local.mkdir(parents=True)
key=pathlib.Path.home()/'.ssh/id_ed25519_cursor';options=['-i',str(key),'-o','BatchMode=yes']
ssh=['ssh','-p','60093',*options,'root@60.188.112.99']
files=[source/'ocean_probe.f90',source/'build_ocean_probe.py']
if options_cli.cash_karp:files.append(source/'ocean_ck_probe.f90')
if options_cli.adaptive:files.append(source/'ocean_adaptive_probe.f90')
manifest=local/'SHA256SUMS';manifest.write_text(''.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+p.name+'\n' for p in files))
def run(args):
    r=subprocess.run(args,text=True,capture_output=True,timeout=180)
    with (local/'run.log').open('a') as f:f.write(r.stdout+r.stderr)
    print(r.stdout+r.stderr,flush=True);r.check_returncode()
run([*ssh,f'test ! -e {remote} && mkdir {remote}'])
run(['scp','-P','60093',*options,*map(str,files),str(manifest),'root@60.188.112.99:'+remote+'/'])
run([*ssh,f'cd {remote} && sha256sum -c SHA256SUMS && python3 build_ocean_probe.py'+(' --cash-karp' if options_cli.cash_karp else '')+(' --adaptive' if options_cli.adaptive else '')])
for name in ['provenance.json','grid64x32','grid128x64','grid64x32-run.log','grid128x64-run.log']:
    run(['scp','-r','-P','60093',*options,'root@60.188.112.99:'+remote+'/'+name,str(local)])
print('Official reference snapshot:',remote);print('Local references:',local)
