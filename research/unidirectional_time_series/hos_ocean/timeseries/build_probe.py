import pathlib,subprocess,os,hashlib,json
root=pathlib.Path('/home/lxy/green-laplace-unidirectional-time-series-runs/hos-gl-timeseries-20260926-v1');base=pathlib.Path('/home/lxy/green-laplace-unidirectional-time-series/artifacts/hos-ocean-v2.1.0');src=root/'hos-source'
subprocess.run(['git','clone','--no-hardlinks',str(base/'source'),str(src)],check=True)
subprocess.run(['git','-C',str(src),'apply',str(base/'precision-v1/precision.patch')],check=True)
p=src/'sources/HOS/output.f90';s=p.read_text();old="WRITE(99,'(201(ES13.5,X))')";assert s.count(old)==1;s=s.replace(old,"WRITE(99,'(201(ES25.16E3,X))')");p.write_text(s)
(root/'hos-io.patch').write_bytes(subprocess.check_output(['git','-C',str(src),'diff']))
libs=base/'deps/usr/lib/x86_64-linux-gnu';env=os.environ.copy();env['LD_LIBRARY_PATH']=':'.join(str(libs/p) for p in ['', 'lapack','blas']);build=root/'build'
commands=[(['cmake','-S',str(src),'-B',str(build),'-DCMAKE_Fortran_COMPILER='+str(base/'precision-v1/gfortran'),'-DUSE_MPI=OFF','-DUSE_HDF5=OFF','-DBUILD_TESTING=OFF','-DBLAS_LIBRARIES='+str(libs/'blas/libblas.so.3'),'-DLAPACK_LIBRARIES='+str(libs/'lapack/liblapack.so.3')+';'+str(libs/'blas/libblas.so.3')],'configure.log'),(['cmake','--build',str(build),'--target','HOS-Ocean','-j','2'],'build.log')]
for args,name in commands:
    with (root/name).open('w') as f:p=subprocess.run(args,env=env,stdout=f,stderr=subprocess.STDOUT)
    assert p.returncode==0,(root/name).read_text()[-3000:]
(root/'build.json').write_text(json.dumps({'upstream':'4deb3b4913d993c4e6ea16f736e5fc5792e14f12','binary_sha256':hashlib.sha256((build/'sources/HOS-Ocean').read_bytes()).hexdigest(),'patch_sha256':hashlib.sha256((root/'hos-io.patch').read_bytes()).hexdigest()},indent=2))