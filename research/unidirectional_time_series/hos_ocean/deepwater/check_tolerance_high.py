import pathlib,subprocess,os,json,time,shutil,concurrent.futures
base=pathlib.Path('/home/lxy/green-laplace-unidirectional-time-series-runs/hos-deepwater-fourphase-20260926-v1');ins=pathlib.Path('/home/lxy/green-laplace-unidirectional-time-series/artifacts/hos-ocean-v2.1.0');timer=ins/'deps/usr/bin/time';solver=ins/'precision-v1/build-precise/sources/HOS-Ocean';env=os.environ.copy()
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:env[k]='1'
se=env.copy();lib=ins/'deps/usr/lib/x86_64-linux-gnu';se['LD_LIBRARY_PATH']=':'.join(str(lib/p) for p in ['', 'lapack','blas'])
records=[]
for tol in [12,14]:
    root=base/('tol%d'%tol);root.mkdir(exist_ok=True);folder=root/'akp012';folder.mkdir();shutil.copy2(base/'akp012/initial.mat',folder/'initial.mat')
    for phase in [0,90,180,270]:
        name='phase%03d'%phase;d=folder/name;d.mkdir();(d/'Results').mkdir();shutil.copy2(base/'akp012'/name/'Results/3d_ini.dat',d/'Results/3d_ini.dat');y=(base/'akp012'/name/'input.yml').read_text().replace('1.e-10','1.e-%d'%tol);(d/'input.yml').write_text(y)
def sim(args):
    tol,phase=args;d=base/('tol%d'%tol)/'akp012'/('phase%03d'%phase);start=time.time()
    with (d/'hos.log').open('w') as f:p=subprocess.run([str(timer),'-v','-o','hos-resources.txt',str(solver),'input.yml'],cwd=d,env=se,stdout=f,stderr=subprocess.STDOUT)
    assert p.returncode==0;return {'case':str(d.relative_to(base)),'wall_seconds':time.time()-start,'exit_code':p.returncode}
with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:records=list(pool.map(sim,[(t,p) for t in [12,14] for p in [0,90,180,270]]))
for tol in [12,14]:
    root=base/('tol%d'%tol);mat=root/'akp012/initial.mat'
    command="addpath('%s');d=load('%s');d.r.tolerance=1e-%d;save('%s','-struct','d','-v7.3');hos_deepwater('analyze','%s',.12);"%(base,mat,tol,mat,root)
    start=time.time()
    with (root/'analysis.log').open('w') as f:p=subprocess.run([str(timer),'-v','-o',str(root/'analysis-resources.txt'),'/home/lxy/Desktop/matlabr2026a/bin/matlab','-singleCompThread','-batch',command],env=env,cwd=root,stdout=f,stderr=subprocess.STDOUT)
    assert p.returncode==0,(root/'analysis.log').read_text()[-2000:];records.append({'stage':'analyze_tol%d'%tol,'wall_seconds':time.time()-start,'exit_code':p.returncode})
(base/'tolerance-high-status.json').write_text(json.dumps({'state':'completed','runs':records},indent=2))