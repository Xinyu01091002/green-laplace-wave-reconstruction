"""Fresh GPU snapshot and a 0.2 s diagnostic pilot; longer matrix is a separate command."""
import datetime,hashlib,pathlib,subprocess,tarfile,json
source=pathlib.Path(__file__).resolve().parent;repo=source.parents[2]
stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ');remote='/root/hos-long-'+stamp
local=repo/'artifacts/hos_gpu'/('long-'+stamp);local.mkdir(parents=True)
key=pathlib.Path.home()/'.ssh/id_ed25519_cursor';proxy=f'ssh -p 60093 -i {key.as_posix()} -o BatchMode=yes -W %h:%p root@60.188.112.99'
opts=['-i',str(key),'-o','BatchMode=yes','-o','ProxyCommand='+proxy];ssh=['ssh',*opts,'root@192.168.2.92']
files=sorted(p for p in source.iterdir() if p.suffix in ['.cu','.cuh','.cpp'] or p.name in ['CMakeLists.txt','long_gpu_matrix.py'])
manifest=local/'SHA256SUMS';manifest.write_text(''.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+p.name+'\n' for p in files))
bundle=local/'source.tar.gz'
with tarfile.open(bundle,'w:gz') as tar:
    for p in [*files,manifest]:tar.add(p,arcname=p.name)
def run(args):
    r=subprocess.run(args,text=True,capture_output=True,timeout=300)
    with (local/'deploy.log').open('a') as f:f.write(r.stdout+r.stderr)
    print(r.stdout+r.stderr,flush=True);r.check_returncode()
run([*ssh,'test ! -e '+remote+' && mkdir '+remote]);run(['scp',*opts,str(bundle),'root@192.168.2.92:'+remote+'/'])
run([*ssh,f'cd {remote} && tar -xzf source.tar.gz && sha256sum -c SHA256SUMS && cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_CUDA_COMPILER=/usr/local/cuda/bin/nvcc -DCMAKE_CUDA_ARCHITECTURES=61 && cmake --build build --target hos_run hos_compare_runs -j2 && ./build/hos_run double /root/hos-production-check-20260927T025257Z/low-adaptive.bin {remote}/pilot 0.2 1e-12 && cat {remote}/pilot/diagnostics.csv'])
(local/'location.json').write_text(json.dumps({'remote':remote},indent=2));print('Long-run snapshot:',remote,flush=True)
