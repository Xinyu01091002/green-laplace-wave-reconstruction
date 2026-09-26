import pathlib,subprocess,os,shutil,hashlib,json
base=pathlib.Path('/home/lxy/green-laplace-unidirectional-time-series/artifacts/hos-ocean-v2.1.0'); root=base/'precision-v1'
def run(args,log):
    with (root/log).open('w') as f:
        p=subprocess.run(args,stdout=f,stderr=subprocess.STDOUT,env=env)
    if p.returncode: raise RuntimeError((root/log).read_text()[-5000:])
    print(log+' OK',flush=True)
env=os.environ.copy(); libs=base/'deps/usr/lib/x86_64-linux-gnu'; env['LD_LIBRARY_PATH']=':'.join(str(libs/p) for p in ['', 'lapack','blas'])
wrapper=root/'gfortran'
wrapper.write_text('#!/bin/sh\nexec '+str(root/'compiler/usr/bin/x86_64-linux-gnu-gfortran-13')+' -B'+str(root/'compiler/usr/libexec/gcc/x86_64-linux-gnu/13/')+'/ -B'+str(root/'compiler/usr/lib/gcc/x86_64-linux-gnu/13/')+'/ -B/usr/lib/gcc/x86_64-linux-gnu/13/ "$@"\n');wrapper.chmod(0o755)
# Private development-library links resolve against installed runtime libraries.
private=root/'compiler/usr/lib/gcc/x86_64-linux-gnu/13/libgfortran.so'
if private.is_symlink() and not private.exists():
    private.unlink();private.symlink_to('/lib/x86_64-linux-gnu/libgfortran.so.5')
for variant in ['baseline','precise']:
    src=root/('source-'+variant)
    subprocess.run(['git','clone','--no-hardlinks',str(base/'source'),str(src)],check=True)
    assert subprocess.check_output(['git','-C',str(src),'rev-parse','HEAD'],text=True).strip()=='4deb3b4913d993c4e6ea16f736e5fc5792e14f12'
    if variant=='precise':
        p=src/'sources/HOS/initial_condition.f90';s=p.read_text();old='READ(unit,1004) eta(i1,i2),phis(i1,i2)';assert s.count(old)==1;s=s.replace(old,'READ(unit,*) eta(i1,i2),phis(i1,i2)');p.write_text(s)
        p=src/'sources/HOS/output.f90';s=p.read_text();start=s.index('IF (i_3D == 1) THEN',s.index('SUBROUTINE output_time_step('));end=s.index('ELSEIF (i_3D == 2) THEN',start);block=s[start:end]
        # Only ordinary physical surface coordinates, eta/psi, and frame-time format.
        for old,new in [('102 FORMAT(3(ES12.5,X),ES12.5)','102 FORMAT(3(ES25.16E3,X),ES25.16E3)'),('103 FORMAT(A,ES12.5,A,I5,A,I5)','103 FORMAT(A,ES25.16E3,A,I5,A,I5)'),('104 FORMAT((ES12.5,X),ES12.5)','104 FORMAT((ES25.16E3,X),ES25.16E3)')]:
            assert block.count(old)==1;block=block.replace(old,new)
        p.write_text(s[:start]+block+s[end:])
        (root/'precision.patch').write_bytes(subprocess.check_output(['git','-C',str(src),'diff','--','sources/HOS/initial_condition.f90','sources/HOS/output.f90']))
    build=root/('build-'+variant)
    run(['cmake','-S',str(src),'-B',str(build),'-DCMAKE_Fortran_COMPILER='+str(wrapper),'-DUSE_MPI=OFF','-DUSE_HDF5=OFF','-DBUILD_TESTING=OFF','-DBLAS_LIBRARIES='+str(libs/'blas/libblas.so.3'),'-DLAPACK_LIBRARIES='+str(libs/'lapack/liblapack.so.3')+';'+str(libs/'blas/libblas.so.3')],variant+'-configure.log')
    run(['cmake','--build',str(build),'--target','HOS-Ocean','-j','2'],variant+'-build.log')
    print([str(p) for p in build.rglob('HOS-Ocean') if p.is_file()],flush=True)
(root/'build-provenance.json').write_text(json.dumps({'source_commit':'4deb3b4913d993c4e6ea16f736e5fc5792e14f12','patch_sha256':hashlib.sha256((root/'precision.patch').read_bytes()).hexdigest(),'binaries':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in root.glob('build-*/**/HOS-Ocean') if p.is_file()},'compiler_packages':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (root/'compiler-debs').glob('*.deb')}},indent=2))