"""Read existing initial fields, compare GPU precision, retain raw arrays remotely."""
import argparse,datetime,hashlib,json,pathlib,subprocess,sys
parser=argparse.ArgumentParser();parser.add_argument('--adaptive',action='store_true');args_cli=parser.parse_args()
source=pathlib.Path(__file__).resolve().parent;repo=source.parents[2]
stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
cpu='/home/lxy/green-laplace-unidirectional-time-series-runs/gpu-production-probe-'+stamp
gpu='/root/hos-production-check-'+stamp
local=repo/'artifacts/hos_gpu'/('production-'+stamp);local.mkdir(parents=True)
key=pathlib.Path.home()/'.ssh/id_ed25519_cursor'
common=['-i',str(key),'-o','BatchMode=yes','-o','ConnectTimeout=10']
proxy=f'ssh -p 60093 -i {key.as_posix()} -o BatchMode=yes -W %h:%p root@60.188.112.99'
gpu_opts=[*common,'-o','ProxyCommand='+proxy]
ssh93=['ssh','-p','60093',*common,'root@60.188.112.99'];ssh92=['ssh',*gpu_opts,'root@192.168.2.92']
def run(args,timeout=300):
    r=subprocess.run(args,text=True,capture_output=True,timeout=timeout)
    with (local/'run.log').open('a') as f:f.write(r.stdout+r.stderr)
    print(r.stdout+r.stderr,flush=True);r.check_returncode();return r.stdout
def manifest(files,name):
    p=local/name;p.write_text(''.join(hashlib.sha256(f.read_bytes()).hexdigest()+'  '+f.name+'\n' for f in files));return p
run([*ssh93,f'test ! -e {cpu} && mkdir {cpu}'])
files=[source/'ocean_probe.f90',source/'build_ocean_probe.py']
if args_cli.adaptive:files.append(source/'ocean_adaptive_probe.f90')
m=manifest(files,'cpu-SHA256SUMS')
run(['scp','-P','60093',*common,*map(str,files),str(m),'root@60.188.112.99:'+cpu+'/'])
print('Generating unchanged-CPU-module reference from low/high production initial fields',flush=True)
run([*ssh93,f'cd {cpu} && sha256sum -c cpu-SHA256SUMS && python3 build_ocean_probe.py --production'+(' --adaptive' if args_cli.adaptive else '')],timeout=1800 if args_cli.adaptive else 600)
for name in ['provenance.json','actual_input_provenance.json','jonswap_low_phi000-run.log','jonswap_high_phi000-run.log']:
    run(['scp','-P','60093',*common,'root@60.188.112.99:'+cpu+'/'+name,str(local/name)])
run([*ssh92,f'test ! -e {gpu} && mkdir {gpu}'])
files=sorted(p for p in source.iterdir() if p.suffix in {'.cu','.cuh'} or p.name=='CMakeLists.txt');m=manifest(files,'gpu-SHA256SUMS')
run(['scp',*gpu_opts,*map(str,files),str(m),'root@192.168.2.92:'+gpu+'/'])
target='hos_adaptive_check' if args_cli.adaptive else 'hos_ocean_check'
run([*ssh92,f'cd {gpu} && sha256sum -c gpu-SHA256SUMS && cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_CUDA_COMPILER=/usr/local/cuda/bin/nvcc -DCMAKE_CUDA_ARCHITECTURES=61 && cmake --build build --target {target} -j2'])
# Transfer over the bounded, IP-restricted LAN helper; no raw local copies.
subprocess.run([sys.executable,str(source/'resume_ocean_production.py'),stamp]+(['--adaptive'] if args_cli.adaptive else []),check=True)
