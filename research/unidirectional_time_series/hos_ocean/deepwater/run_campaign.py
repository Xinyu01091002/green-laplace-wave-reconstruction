import pathlib,subprocess,os,json,time,hashlib,concurrent.futures
root=pathlib.Path('/home/lxy/green-laplace-unidirectional-time-series-runs/hos-deepwater-fourphase-20260926-v1')
install=pathlib.Path('/home/lxy/green-laplace-unidirectional-time-series/artifacts/hos-ocean-v2.1.0')
solver=install/'precision-v1/build-precise/sources/HOS-Ocean';matlab='/home/lxy/Desktop/matlabr2026a/bin/matlab';timer=install/'deps/usr/bin/time'
env=os.environ.copy()
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:env[k]='1'
state={'state':'started','started_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'solver_sha256':hashlib.sha256(solver.read_bytes()).hexdigest(),'source_hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [root/'hos_deepwater.m',root/'run_campaign.py',root/'mf12-no31-source.tar.gz']},'runs':[]}
def save(): (root/'status.json').write_text(json.dumps(state,indent=2))
def run(cmd,folder,name,runenv):
    start=time.time()
    with (folder/(name+'.log')).open('w') as f:p=subprocess.run([str(timer),'-v','-o',str(folder/(name+'-resources.txt')),*cmd],cwd=folder,env=runenv,stdout=f,stderr=subprocess.STDOUT)
    record={'name':str(folder.relative_to(root))+'/'+name,'exit_code':p.returncode,'wall_seconds':time.time()-start}
    if p.returncode: raise RuntimeError(json.dumps(record)+'\n'+(folder/(name+'.log')).read_text()[-3000:])
    return record
try:
    save()
    for akp in [.02,.12]:
        name='akp%03d'%round(100*akp);state['state']='preparing_'+name;save()
        command="addpath('%s');hos_deepwater('prepare','%s',%.17g);"%(root,root,akp)
        state['runs'].append(run([matlab,'-singleCompThread','-batch',command],root,'prepare-'+name,env));save()
        state['state']='running_'+name;save();se=env.copy();libs=install/'deps/usr/lib/x86_64-linux-gnu';se['LD_LIBRARY_PATH']=':'.join(str(libs/p) for p in ['', 'lapack','blas'])
        def simulate(phase):return run([str(solver),'input.yml'],root/name/('phase%03d'%phase),'hos',se)
        with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:state['runs']+=list(pool.map(simulate,[0,90,180,270]))
        state['state']='analyzing_'+name;save()
        command="addpath('%s');hos_deepwater('analyze','%s',%.17g);"%(root,root,akp)
        state['runs'].append(run([matlab,'-singleCompThread','-batch',command],root,'analyze-'+name,env));save()
    state['state']='completed';state['finished_utc']=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime());save()
except Exception as e:
    state['state']='failed';state['error']=str(e);save();raise