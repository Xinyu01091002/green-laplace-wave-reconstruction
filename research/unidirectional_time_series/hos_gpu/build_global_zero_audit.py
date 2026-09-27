"""Isolated diagnostic correction; never edits the production source/build."""
import difflib,hashlib,json,os,pathlib,shlex,subprocess,tarfile
root=pathlib.Path(__file__).resolve().parent
runs=pathlib.Path('/home/lxy/green-laplace-unidirectional-time-series-runs');base=runs/'hos-directional-fourphase-20260926-v2';dependency=runs/'hos-mpi-check-20260926-v1'
build=base/'build';source=base/'source';deps=dependency/'deps/usr';lib=deps/'lib/x86_64-linux-gnu';env=dict(os.environ);env.update(json.loads((dependency/'environment.json').read_text()))
original=source/'sources/HOS/resol_HOS.f90';text=original.read_text();needle='da_phisrk(1,1) = da_phisrk(1,1) - g_star*a_etark(1,1)';assert text.count(needle)==1
replacement='#ifdef MPI\nIF (local_y_start == 0) THEN\n    '+needle+'\nENDIF\n#else\n'+needle+'\n#endif'
changed=root/'resol_HOS_global_zero.f90';changed.write_text(text.replace(needle,replacement))
(root/'mpi-global-zero.patch').write_text(''.join(difflib.unified_diff(text.splitlines(True),changed.read_text().splitlines(True),fromfile='a/sources/HOS/resol_HOS.f90',tofile='b/sources/HOS/resol_HOS.f90')))
compiler='/home/lxy/green-laplace-unidirectional-time-series/artifacts/hos-ocean-v2.1.0/precision-v1/gfortran'
inc=[build/'sources/mod',build/'_deps/yaml_parser-build/mod',deps/'include',lib/'fortran/gfortran-mod-15/openmpi',lib/'openmpi/lib']
cmd=[compiler,'-cpp','-DMPI','-D_GNU_FORTRAN_COMPILE_RULE_','-O3','-DNDEBUG','-fPIC','-msse2','-funroll-loops','-fno-protect-parens','-ffast-math','-pthread','-fallow-argument-mismatch','-ffree-line-length-none','-J'+str(root),*['-I'+str(p) for p in inc],'-c',str(changed),'-o',str(root/'resol-fixed.o')]
with (root/'compile.log').open('w') as log:subprocess.run(cmd,env=env,cwd=root,stdout=log,stderr=subprocess.STDOUT,check=True)
main=(source/'sources/HOS/HOS-ocean.f90').read_text();needle_main='REAL(RP) :: kp_real';assert main.count(needle_main)==1
main=main.replace(needle_main,needle_main+'\nREAL(RP) :: bench_started,bench_elapsed')
needle_main='DO WHILE (T_stop_star-time_cur >= -tiny)';assert main.count(needle_main)==1
main=main.replace(needle_main,'CALL MPI_BARRIER(MPI_COMM_WORLD,Statinfo)\nbench_started=MPI_WTIME()\n'+needle_main)
main=main.replace('CALL MPI_BARRIER(Statinfo)','CALL MPI_BARRIER(MPI_COMM_WORLD,Statinfo)')
assert main.count('CLOSE(1)')==1
main=main.replace('CLOSE(1)', '''bench_elapsed=MPI_WTIME()-bench_started
IF(rank==0)THEN
    WRITE(*,'(A,ES25.16)')'BENCH_ADVANCE_SECONDS ',bench_elapsed
    WRITE(*,'(A,I12,1X,I12)')'BENCH_STEPS ',n_rk_tot,n_er_tot
ENDIF
CLOSE(1)''')
main_path=root/'main_timed.f90';main_path.write_text(main)
main_cmd=cmd.copy();main_cmd[main_cmd.index(str(changed))]=str(main_path);main_cmd[main_cmd.index(str(root/'resol-fixed.o'))]=str(root/'main.o')
main_cmd+=['-D__GIT_BRANCH__="benchmark"','-D__GIT_COMMIT_HASH__="original-modules"']
with (root/'main-compile.log').open('w') as log:subprocess.run(main_cmd,env=env,cwd=root,stdout=log,stderr=subprocess.STDOUT,check=True)
args=shlex.split((build/'sources/CMakeFiles/HOS-Ocean.dir/link.txt').read_text())
for i,a in enumerate(args):
    if a.endswith('.o'):
        if a.endswith('/HOS/HOS-ocean.f90.o'):args[i]=str(root/'main.o')
        elif a.endswith('/HOS/resol_HOS.f90.o'):args[i]=str(root/'resol-fixed.o')
        else:args[i]=str((build/'sources'/a).resolve())
exe=root/'HOS-benchmark-globalzero-fixed';args[args.index('-o')+1]=str(exe)
with (root/'link.log').open('w') as log:subprocess.run(args,env=env,cwd=root,stdout=log,stderr=subprocess.STDOUT,check=True)
unmodified=root/'HOS-benchmark-mpi-sync';original_args=args.copy()
original_args[original_args.index(str(root/'resol-fixed.o'))]=str(build/'sources/CMakeFiles/HOS-Ocean.dir/HOS/resol_HOS.f90.o')
original_args[original_args.index('-o')+1]=str(unmodified)
with (root/'unmodified-link.log').open('w') as log:subprocess.run(original_args,env=env,cwd=root,stdout=log,stderr=subprocess.STDOUT,check=True)
meta={'scope':'Diagnostic copy only. Nonlinear gravity belongs only to the global spatial zero mode. Nonzero modes already have linear gravity in the integrating factor.',
 'original_source_sha256':hashlib.sha256(original.read_bytes()).hexdigest(),'changed_source_sha256':hashlib.sha256(changed.read_bytes()).hexdigest(),
 'binary_sha256':hashlib.sha256(exe.read_bytes()).hexdigest(),'unmodified_physics_binary_sha256':hashlib.sha256(unmodified.read_bytes()).hexdigest(),
 'timed_main_sha256':hashlib.sha256(main_path.read_bytes()).hexdigest(),'driver_note':'MPI_BARRIER calls use the required communicator and ierr arguments. Evolution/controller formulas unchanged.', 'replacement':replacement}
(root/'provenance.json').write_text(json.dumps(meta,indent=2))
with tarfile.open(root/'audit.tar.gz','w:gz') as tar:
    for p in [exe,unmodified,root/'provenance.json',changed,main_path,root/'mpi-global-zero.patch']:tar.add(p,arcname=p.name)
print('Global-zero audit binary ready',hashlib.sha256((root/'audit.tar.gz').read_bytes()).hexdigest(),flush=True)
