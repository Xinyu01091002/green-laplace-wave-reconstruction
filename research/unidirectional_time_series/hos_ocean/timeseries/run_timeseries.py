import pathlib,subprocess,os,json,time,hashlib,concurrent.futures,shutil,re
root=pathlib.Path('/home/lxy/green-laplace-unidirectional-time-series-runs/hos-gl-timeseries-20260926-v1');old=pathlib.Path('/home/lxy/green-laplace-unidirectional-time-series-runs/hos-deepwater-fourphase-20260926-v1');ins=pathlib.Path('/home/lxy/green-laplace-unidirectional-time-series/artifacts/hos-ocean-v2.1.0')
(root/'mf12').symlink_to(old/'mf12',target_is_directory=True);(root/'ow-reference').mkdir();subprocess.run(['tar','-xzf',str(root/'ow-timeseries-reference.tar.gz'),'-C',str(root/'ow-reference')],check=True)
solver=root/'build/sources/HOS-Ocean';timer=ins/'deps/usr/bin/time';matlab='/home/lxy/Desktop/matlabr2026a/bin/matlab';env=os.environ.copy()
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:env[k]='1'
se=env.copy();libs=ins/'deps/usr/lib/x86_64-linux-gnu';se['LD_LIBRARY_PATH']=':'.join(str(libs/p) for p in ['', 'lapack','blas'])
state={'state':'started','started_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'runs':[],'hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [root/'prepare_hos_timeseries.m',root/'analyze_hos_timeseries.m',root/'run_timeseries.py',root/'gl-timeseries-source.tar.gz',root/'ow-timeseries-reference.tar.gz',solver]}}
def save(): (root/'status.json').write_text(json.dumps(state,indent=2))
def run(args,folder,name,environ):
    start=time.time()
    with (folder/(name+'.log')).open('w') as f:p=subprocess.run([str(timer),'-v','-o',str(folder/(name+'-resources.txt')),*args],cwd=folder,env=environ,stdout=f,stderr=subprocess.STDOUT)
    entry={'name':str(folder.relative_to(root))+'/'+name,'wall_seconds':time.time()-start,'exit_code':p.returncode}
    assert p.returncode==0,str(entry)+'\n'+(folder/(name+'.log')).read_text()[-2500:]
    return entry
try:
    for akp in [.02,.12]:
        name='akp%03d'%round(akp*100);state['state']='preparing_'+name;save()
        cmd="addpath('%s');prepare_hos_timeseries('%s',%.17g);"%(root,root,akp)
        state['runs'].append(run([matlab,'-singleCompThread','-batch',cmd],root,'prepare-'+name,env));save()
        if akp==.02:
            case=root/name/'phase000';smoke=root/'probe-smoke';smoke.mkdir();(smoke/'Results').mkdir();shutil.copy2(case/'Results/3d_ini.dat',smoke/'Results/3d_ini.dat');shutil.copy2(case/'prob.inp',smoke/'prob.inp')
            settings=json.loads((root/name/'initialization.json').read_text());text=(case/'input.yml').read_text();text=re.sub(r'duration: [^\n]+','duration: '+str(settings['output_dt']),text);text=text.replace('physical space: false','physical space: true');(smoke/'input.yml').write_text(text)
            state['runs'].append(run([str(solver),'input.yml'],smoke,'hos',se))
            lines=(smoke/'Results/probes.dat').read_text().splitlines();i=next(i for i,x in enumerate(lines) if x.startswith('ZONE'));pv=[float(x) for x in lines[i+1].split()]
            lines=(smoke/'Results/3d.dat').read_text().splitlines();i=next(i for i,x in enumerate(lines) if x.startswith('ZONE'));field=[float(x) for x in lines[i+settings['probe_index']].split()]
            input_lines=(smoke/'Results/3d_ini.dat').read_text().splitlines();expected=float(input_lines[67+settings['probe_index']].split()[0]);assert abs(pv[1]-field[2])<1e-11 and abs(pv[1]-expected)<1e-11
            (smoke/'io-check.json').write_text(json.dumps({'probe_eta':pv[1],'native_eta':field[2],'input_eta':expected,'passed':True},indent=2));save()
        state['state']='running_'+name;save()
        def sim(phase):return run([str(solver),'input.yml'],root/name/('phase%03d'%phase),'hos',se)
        with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:state['runs']+=list(pool.map(sim,[0,90,180,270]))
        state['state']='analyzing_'+name;save();cmd="addpath('%s');analyze_hos_timeseries('%s',%.17g);"%(root,root,akp)
        state['runs'].append(run([matlab,'-singleCompThread','-batch',cmd],root,'analyze-'+name,env));save()
    state['state']='completed';state['finished_utc']=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime());save()
except Exception as e:state['state']='failed';state['error']=str(e);save();raise