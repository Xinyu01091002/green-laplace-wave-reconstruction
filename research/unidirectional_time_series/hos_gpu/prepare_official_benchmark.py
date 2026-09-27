"""On 93: build timed official MPI main, package dependencies and frozen prefixes.
Original module objects, production inputs and running jobs remain unchanged.
"""
import hashlib,json,os,pathlib,re,shlex,shutil,subprocess,tarfile
root=pathlib.Path(__file__).resolve().parent
runs=pathlib.Path('/home/lxy/green-laplace-unidirectional-time-series-runs')
base=runs/'hos-directional-fourphase-20260926-v2';dependency=runs/'hos-mpi-check-20260926-v1'
artifact=pathlib.Path('/home/lxy/green-laplace-unidirectional-time-series/artifacts/hos-ocean-v2.1.0')
build=base/'build';source=base/'source';deps=dependency/'deps/usr';lib=deps/'lib/x86_64-linux-gnu'
env=dict(os.environ);env.update(json.loads((dependency/'environment.json').read_text()))
runtime=root/'runtime';runtime.mkdir();(runtime/'bin').mkdir()
main=(source/'sources/HOS/HOS-ocean.f90').read_text()
needle='REAL(RP) :: kp_real';assert main.count(needle)==1
main=main.replace(needle,needle+'\nREAL(RP) :: bench_started,bench_elapsed')
needle='DO WHILE (T_stop_star-time_cur >= -tiny)';assert main.count(needle)==1
main=main.replace(needle,'CALL MPI_BARRIER(MPI_COMM_WORLD,Statinfo)\nbench_started=MPI_WTIME()\n'+needle)
main=main.replace('CALL MPI_BARRIER(Statinfo)','CALL MPI_BARRIER(MPI_COMM_WORLD,Statinfo)')
needle='CLOSE(1)';assert main.count(needle)==1
main=main.replace(needle,'''    bench_elapsed=MPI_WTIME()-bench_started
    IF(rank==0)THEN
        WRITE(*,'(A,ES25.16)')'BENCH_ADVANCE_SECONDS ',bench_elapsed
        WRITE(*,'(A,I12,1X,I12)')'BENCH_STEPS ',n_rk_tot,n_er_tot
    ENDIF
CLOSE(1)''')
driver=root/'HOS-benchmark.f90';driver.write_text(main)
compiler=artifact/'precision-v1/gfortran'
includes=[build/'sources/mod',build/'_deps/yaml_parser-build/mod',deps/'include',lib/'fortran/gfortran-mod-15/openmpi',lib/'openmpi/lib']
cmd=[str(compiler),'-cpp','-DMPI','-D_GNU_FORTRAN_COMPILE_RULE_','-D__GIT_BRANCH__="benchmark"','-D__GIT_COMMIT_HASH__="original-modules"',
 '-O3','-DNDEBUG','-fPIC','-msse2','-funroll-loops','-fno-protect-parens','-ffast-math','-pthread','-fallow-argument-mismatch','-ffree-line-length-none',
 *['-I'+str(p) for p in includes],'-c',str(driver),'-o',str(root/'main.o')]
with (root/'compile.log').open('w') as f:subprocess.run(cmd,env=env,stdout=f,stderr=subprocess.STDOUT,check=True)
args=shlex.split((build/'sources/CMakeFiles/HOS-Ocean.dir/link.txt').read_text());objects=[]
for i,a in enumerate(args):
    if a.endswith('.o'):
        path=(build/'sources'/a).resolve();objects.append({'path':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
        args[i]=str(root/'main.o') if a.endswith('/HOS/HOS-ocean.f90.o') else str(path)
exe=runtime/'bin/HOS-benchmark';args[args.index('-o')+1]=str(exe)
with (root/'link.log').open('w') as f:subprocess.run(args,env=env,stdout=f,stderr=subprocess.STDOUT,check=True)
shutil.copytree(deps,runtime/'mpi/usr',symlinks=True)
support=runtime/'support';support.mkdir()
# Resolve direct/transitive ELF dependencies of the executable and MPI launcher.
copied={}
for binary in [exe,deps/'bin/mpirun.openmpi',*list(lib.glob('libpmix.so*')),*list(lib.glob('libhwloc.so*'))]:
    result=subprocess.run(['ldd',str(binary)],env=env,text=True,capture_output=True)
    for line in result.stdout.splitlines():
        if 'not found' in line:raise RuntimeError(line)
        match=re.search(r'(\S+)\s+=>\s+(/\S+)',line)
        if not match:continue
        name,path=match.groups()
        if name in ['libc.so.6','libm.so.6','libmvec.so.1','libpthread.so.0','libdl.so.2','librt.so.1']:continue
        if name not in copied:
            shutil.copy2(path,support/name);copied[name]=hashlib.sha256((support/name).read_bytes()).hexdigest()
campaign=runs/'hos-jonswap-r4gl-80tp-20260926-v2';input_records=[]
for family in ['low','high']:
    origin=campaign/family/'cases/phi000';dest=runtime/'cases'/family;dest.mkdir(parents=True);(dest/'Results').mkdir()
    for p in sorted((origin/'Results').glob('3d_ini_*.dat')):
        shutil.copy2(p,dest/'Results'/p.name);input_records.append({'family':family,'file':p.name,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
    (dest/'input.yml').write_text((origin/'input.yml').read_text().replace('duration: 1100.8','duration: 2.0'))
    shutil.copy2(origin/'prob.inp',dest/'prob.inp')
(runtime/'rank.sh').write_text('#!/bin/sh\nexec taskset -c "$OMPI_COMM_WORLD_LOCAL_RANK" "$HOS_BENCH_EXE" input.yml\n');(runtime/'rank.sh').chmod(0o755)
refs=root/'references';refs.mkdir();provenance=[]
def records(path,end):
    lines=path.read_text().splitlines();index=next(i for i,line in enumerate(lines) if 'VARIABLES' in line)
    return [line.split() for line in lines[index+1:] if line.strip() and float(line.split()[0])<=end+1e-8]
for family,end,folder in [('low',27.6,campaign/'low/cases/phi000/Results'),('high',2.0,campaign/'short-high/cases/phi000/Results')]:
    probes=records(folder/'probes.dat',end);energy=records(folder/'vol_energy.dat',end);assert len(probes)==round(end/.2)+1
    assert probes==records(folder/'probes.dat',end) and energy==records(folder/'vol_energy.dat',end)
    edict={round(float(r[0]),6):r for r in energy}
    path=refs/(family+'-official.csv')
    with path.open('w') as f:
        f.write('time,p1,p2,p3,p4,p5,potential,kinetic,total,relative_energy_change\n')
        for row in probes:
            e=edict[round(float(row[0]),6)];f.write(','.join(row[:6]+[e[2],e[3],e[4],e[5]])+'\n')
    provenance.append({'family':family,'source':str(folder),'end_time_s':end,'rows':len(probes),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
meta={'production_binary_sha256':hashlib.sha256((build/'sources/HOS-Ocean').read_bytes()).hexdigest(),'benchmark_binary_sha256':hashlib.sha256(exe.read_bytes()).hexdigest(),
 'instrumented_main_sha256':hashlib.sha256(driver.read_bytes()).hexdigest(),'original_objects':objects,'support_libraries':copied,'initial_files':input_records,'frozen_references':provenance}
(runtime/'provenance.json').write_text(json.dumps(meta,indent=2));(refs/'provenance.json').write_text(json.dumps(provenance,indent=2));(root/'package-provenance.json').write_text(json.dumps(meta,indent=2))
for name,folder in [('cpu-runtime.tar.gz',runtime),('reference-records.tar.gz',refs)]:
    with tarfile.open(root/name,'w:gz') as tar:tar.add(folder,arcname=folder.name)
    print(name,(root/name).stat().st_size,hashlib.sha256((root/name).read_bytes()).hexdigest(),flush=True)
