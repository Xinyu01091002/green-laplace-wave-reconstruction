"""Continue a still-running/completed reference after a client wait timeout.
Never reruns or stops the CPU scientific process. Prepare GPU while it finishes.
"""
import hashlib,pathlib,re,subprocess,sys,tarfile,time
source=pathlib.Path(__file__).resolve().parent;repo=source.parents[2];stamp=sys.argv[1]
assert re.fullmatch(r'\d{8}T\d{6}Z',stamp)
cpu='/home/lxy/green-laplace-unidirectional-time-series-runs/gpu-production-probe-'+stamp;gpu='/root/hos-production-check-'+stamp
local=repo/'artifacts/hos_gpu'/('production-'+stamp);assert local.is_dir()
key=pathlib.Path.home()/'.ssh/id_ed25519_cursor';common=['-i',str(key),'-o','BatchMode=yes','-o','ConnectTimeout=10']
proxy=f'ssh -p 60093 -i {key.as_posix()} -o BatchMode=yes -W %h:%p root@60.188.112.99'
options=[*common,'-o','ProxyCommand='+proxy]
ssh93=['ssh','-p','60093',*common,'root@60.188.112.99'];ssh92=['ssh',*options,'root@192.168.2.92']
def run(args,timeout=180):
    r=subprocess.run(args,capture_output=True,text=True,timeout=timeout)
    with (local/'continuation.log').open('a') as log:log.write(r.stdout+r.stderr)
    print(r.stdout+r.stderr,flush=True);r.check_returncode();return r.stdout
run([*ssh93,'test -f '+cpu+'/provenance.json'])
run([*ssh92,'test ! -e '+gpu+' && mkdir '+gpu])
files=sorted(p for p in source.iterdir() if p.suffix in {'.cu','.cuh'} or p.name=='CMakeLists.txt')
manifest=local/'gpu-SHA256SUMS';manifest.write_text(''.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+p.name+'\n' for p in files))
bundle=local/'gpu-source.tar.gz'
with tarfile.open(bundle,'w:gz') as archive:
    for p in [*files,manifest]:archive.add(p,arcname=p.name,recursive=False)
run(['scp',*options,str(bundle),'root@192.168.2.92:'+gpu+'/'])
run([*ssh92,f'cd {gpu} && tar -xzf gpu-source.tar.gz && sha256sum -c gpu-SHA256SUMS && cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_CUDA_COMPILER=/usr/local/cuda/bin/nvcc -DCMAKE_CUDA_ARCHITECTURES=61 && cmake --build build --target hos_adaptive_check -j2'])
start=time.monotonic();polls=0
while True:
    r=subprocess.run([*ssh93,'test -f '+cpu+'/actual_input_provenance.json'],capture_output=True,text=True,timeout=30)
    if r.returncode==0:break
    if time.monotonic()-start>1200:raise TimeoutError('Reference still pending; not restarted or killed')
    if polls%4==0:print('Waiting for original CPU reference to finish; no restart',flush=True)
    polls+=1;time.sleep(15)
for name in ['provenance.json','actual_input_provenance.json','jonswap_low_phi000-run.log','jonswap_high_phi000-run.log']:
    run(['scp','-P','60093',*common,'root@60.188.112.99:'+cpu+'/'+name,str(local/name)])
subprocess.run([sys.executable,str(source/'resume_ocean_production.py'),stamp,'--adaptive'],check=True)
