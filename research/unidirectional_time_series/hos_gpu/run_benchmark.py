"""Deploy fresh benchmark source snapshot and collect compact evidence only."""
import datetime
import hashlib
import pathlib
import subprocess

source=pathlib.Path(__file__).resolve().parent;repo=source.parents[2]
stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
remote='/root/hos-benchmark-'+stamp;local=repo/'artifacts/hos_gpu'/('benchmark-'+stamp)
local.mkdir(parents=True,exist_ok=False)
key=pathlib.Path.home()/'.ssh/id_ed25519_cursor'
proxy=f'ssh -p 60093 -i {key.as_posix()} -o BatchMode=yes -W %h:%p root@60.188.112.99'
options=['-i',str(key),'-o','BatchMode=yes','-o','ConnectTimeout=10','-o','ProxyCommand='+proxy]
ssh=['ssh',*options,'root@192.168.2.92']
files=sorted(p for p in source.iterdir() if p.suffix in {'.cu','.cuh','.cpp','.hpp'} or p.name in {'CMakeLists.txt','generate_cpu_backend.py','benchmark_matrix.py'})
manifest=local/'SHA256SUMS';manifest.write_text(''.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+p.name+'\n' for p in files))
def run(args):
    with (local/'run.log').open('a',encoding='utf-8') as log:
        p=subprocess.Popen(args,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
        for line in p.stdout:log.write(line);log.flush();print(line,end='',flush=True)
        if p.wait():raise RuntimeError('Command failed; see '+str(local/'run.log'))
run([*ssh,f'test ! -e {remote} && mkdir {remote}'])
run(['scp',*options,*map(str,files),str(manifest),'root@192.168.2.92:'+remote+'/'])
run([*ssh,f'cd {remote} && sha256sum -c SHA256SUMS && cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_CUDA_COMPILER=/usr/local/cuda/bin/nvcc -DCMAKE_CUDA_ARCHITECTURES=61 -DHOS_BUILD_BENCHMARK=ON && cmake --build build -j2 && ctest --test-dir build --output-on-failure && python3 benchmark_matrix.py'])
for name in ['timings.json','comparisons.json','hardware.json']:
    run(['scp',*options,'root@192.168.2.92:'+remote+'/measurements/'+name,str(local/name)])
print('Benchmark snapshot: '+remote,flush=True)
print('Local evidence: '+str(local),flush=True)
