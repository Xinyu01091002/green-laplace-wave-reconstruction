"""Deploy a hashed fresh source snapshot; compile and run bounded checks.

Uses the workstation key through SSH forwarding (no agent/key copying).
Does not alter GPU hardware/driver settings or existing run directories.
Python performs orchestration only, never scientific reference calculations.
"""
import datetime
import argparse
import hashlib
import pathlib
import shlex
import subprocess
import tarfile

parser = argparse.ArgumentParser()
parser.add_argument('--bands', action='store_true', help='Run the fixed-support precision diagnosis')
parser.add_argument('--orders', type=pathlib.Path, help='Directory containing MATLAB-generated order-five reference binaries')
parser.add_argument('--time-reference', type=pathlib.Path, help='MATLAB-generated hos_time.bin for short RK4 check')
parser.add_argument('--ocean', type=pathlib.Path, help='Directory containing official MPI reference grid folders')
parser.add_argument('--cash-karp', type=pathlib.Path, help='Directory containing original CK attempted-step references')
parser.add_argument('--adaptive', type=pathlib.Path, help='Directory containing original adaptive trajectory references')
args = parser.parse_args()

source = pathlib.Path(__file__).resolve().parent
repo = source.parents[2]
stamp = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
remote = '/root/hos-rhs-check-' + stamp
local = repo / 'artifacts' / 'hos_gpu' / stamp
local.mkdir(parents=True, exist_ok=False)
key = pathlib.Path.home() / '.ssh' / 'id_ed25519_cursor'
proxy = f'ssh -p 60093 -i {key.as_posix()} -o BatchMode=yes -o ConnectTimeout=10 -W %h:%p root@60.188.112.99'
options = ['-i', str(key), '-o', 'IdentitiesOnly=yes', '-o', 'BatchMode=yes',
           '-o', 'ConnectTimeout=10', '-o', 'ProxyCommand=' + proxy]
ssh = ['ssh', *options, 'root@192.168.2.92']
files = sorted(p for p in source.iterdir() if p.suffix in {'.cu', '.cuh'} or p.name == 'CMakeLists.txt')
if args.orders:
    references = sorted(args.orders.glob('hos_h*.bin'))
    if len(references) != 3:
        raise ValueError('Expected three independently generated MATLAB reference files')
    files += [source / 'make_hos_reference.m', *references]
if args.time_reference:
    if not args.time_reference.is_file():
        raise ValueError('Missing independent time reference')
    files += [source / 'make_time_reference.m', args.time_reference]
if args.ocean:
    ocean_files=sorted(args.ocean.glob('grid*/ocean_case*.bin'))
    if len(ocean_files)!=8:raise ValueError('Expected eight official-module reference files')
    for reference in ocean_files:
        staged=local/(reference.parent.name+'_'+reference.name)
        staged.write_bytes(reference.read_bytes());files.append(staged)
    files.append(source/'run_ocean_checks.py')
if args.cash_karp:
    ck_files=sorted(args.cash_karp.glob('grid*/ck_case*.bin'))
    if len(ck_files)!=24:raise ValueError('Expected 24 attempted-step reference files')
    for reference in ck_files:
        staged=local/(reference.parent.name+'_'+reference.name);staged.write_bytes(reference.read_bytes());files.append(staged)
    files.append(source/'run_ck_checks.py')
if args.adaptive:
    adaptive_files=sorted(args.adaptive.glob('grid*/adaptive_case*.bin'))
    if len(adaptive_files)!=6:raise ValueError('Expected six bounded adaptive trajectories')
    for reference in adaptive_files:
        for p in [reference,reference.with_suffix('.trace')]:
            staged=local/(p.parent.name+'_'+p.name);staged.write_bytes(p.read_bytes());files.append(staged)
    files.append(source/'run_adaptive_checks.py')
manifest = local / 'SHA256SUMS'
manifest.write_text(''.join(hashlib.sha256(p.read_bytes()).hexdigest() + '  ' + p.name + '\n' for p in files), encoding='ascii')

def execute(args,timeout=180):
    result = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=timeout)
    with (local / 'run.log').open('a', encoding='utf-8') as log:
        log.write(result.stdout + '\nexit=' + str(result.returncode) + '\n')
    print(result.stdout, flush=True)
    result.check_returncode()

execute([*ssh, f'test ! -e {shlex.quote(remote)} && mkdir {shlex.quote(remote)}'])
bundle=local/'snapshot.tar.gz'
with tarfile.open(bundle,'w:gz') as archive:
    for path in [*files,manifest]:archive.add(path,arcname=path.name,recursive=False)
print('Compressed snapshot bytes:',bundle.stat().st_size,flush=True)
execute(['scp', *options, str(bundle), 'root@192.168.2.92:' + remote + '/'],timeout=300)
commands = [f'cd {remote}', 'tar -xzf snapshot.tar.gz', 'sha256sum -c SHA256SUMS',
    'cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_CUDA_COMPILER=/usr/local/cuda/bin/nvcc -DCMAKE_CUDA_ARCHITECTURES=61',
    'cmake --build build -j2', 'ctest --test-dir build --output-on-failure -V',
    './build/hos_rhs_check float 256 128', './build/hos_rhs_check double 256 128']
if args.bands:
    commands += ['./build/hos_band_check float 2048 1024 0.15',
                 './build/hos_band_check double 2048 1024 0.15',
                 './build/hos_band_check float 2048 1024 1.3',
                 './build/hos_band_check float 2048 1024 20']
if args.orders:
    for precision in ['float', 'double']:
        for ref in references:
            commands.append('./build/hos_order_check ' + precision + ' ' + shlex.quote(ref.name))
        commands.append('./build/hos_order_check ' + precision + ' hos_h1.3.bin 2048 1024')
if args.time_reference:
    commands += ['./build/hos_time_check float ' + shlex.quote(args.time_reference.name),
                 './build/hos_time_check double ' + shlex.quote(args.time_reference.name)]
if args.ocean:
    commands.append('python3 run_ocean_checks.py')
if args.cash_karp:
    commands.append('python3 run_ck_checks.py')
if args.adaptive:
    commands.append('python3 run_adaptive_checks.py')
execute([*ssh, ' && '.join(commands)])
print('Source snapshot:', remote)
print('Local log:', local / 'run.log')
