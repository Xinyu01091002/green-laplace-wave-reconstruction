import pathlib,subprocess,os,shutil,json,time,hashlib,concurrent.futures
base=pathlib.Path('/home/lxy/green-laplace-unidirectional-time-series/artifacts/hos-ocean-v2.1.0');r=base/'precision-v1'
env=os.environ.copy()
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:env[k]='1'
src='/home/lxy/green-laplace-unidirectional-time-series-runs/ow3d16-kpd1-akp012-20260925T212822Z/inputs/wavegroup/initial_fields.mat'
with (r/'prepare.log').open('w') as f:
    p=subprocess.run(['/home/lxy/Desktop/matlabr2026a/bin/matlab','-singleCompThread','-batch',"addpath('%s');prepare_hos_precision('%s','%s');"%(r,src,r/'case-precise-full')],env=env,stdout=f,stderr=subprocess.STDOUT)
assert p.returncode==0,(r/'prepare.log').read_text()[-3000:]
for name in ['case-baseline-low','case-precise-low']:
    case=r/name;case.mkdir();(case/'Results').mkdir()
    for f in ['input_HOS-Ocean.yml','import_reference.mat','input-audit.json','Results/3d_ini.dat']:
        shutil.copy2(base/'wavegroup-phase000-pilot'/f,case/f)
libs=base/'deps/usr/lib/x86_64-linux-gnu';env['LD_LIBRARY_PATH']=':'.join(str(libs/p) for p in ['', 'lapack','blas'])
def execute(name,variant):
    case=r/name;binary=next((r/('build-'+variant)).rglob('HOS-Ocean'))
    report={'state':'running','source_commit':'4deb3b4913d993c4e6ea16f736e5fc5792e14f12','binary':str(binary),'binary_sha256':hashlib.sha256(binary.read_bytes()).hexdigest(),'input_sha256':hashlib.sha256((case/'Results/3d_ini.dat').read_bytes()).hexdigest(),'yaml_sha256':hashlib.sha256((case/'input_HOS-Ocean.yml').read_bytes()).hexdigest(),'started_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())}
    (case/'status.json').write_text(json.dumps(report,indent=2));start=time.time()
    with (case/'run.log').open('w') as f:
        p=subprocess.run([str(base/'deps/usr/bin/time'),'-v','-o','resources.txt',str(binary),'input_HOS-Ocean.yml'],cwd=case,env=env,stdout=f,stderr=subprocess.STDOUT)
    report.update(state='exited_pending_validation',exit_code=p.returncode,wall_seconds=time.time()-start)
    (case/'status.json').write_text(json.dumps(report,indent=2)); print(name,json.dumps(report),flush=True)
    assert p.returncode==0
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
    list(pool.map(lambda args:execute(*args),[('case-baseline-low','baseline'),('case-precise-low','precise'),('case-precise-full','precise')]))