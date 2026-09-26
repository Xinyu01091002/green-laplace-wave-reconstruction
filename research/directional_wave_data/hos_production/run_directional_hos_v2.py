import pathlib,subprocess,os,json,time,hashlib,shutil,re,math,signal
r=pathlib.Path('/home/lxy/green-laplace-unidirectional-time-series-runs/hos-directional-fourphase-20260926-v2');mpi=pathlib.Path('/home/lxy/green-laplace-unidirectional-time-series-runs/hos-mpi-check-20260926-v1');base=pathlib.Path('/home/lxy/green-laplace-unidirectional-time-series/artifacts/hos-ocean-v2.1.0')
source=pathlib.Path('/home/lxy/green-laplace-unidirectional-time-series-runs/ow3d16-kpd1-akp012-20260925T212822Z/inputs/wavegroup/initial_fields.mat')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def utc():return time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())
def atomic(p,data):
    q=p.with_suffix(p.suffix+'.tmp');q.write_text(json.dumps(data,indent=2));q.replace(p)
state={'state':'preparing','created_utc':utc(),'cases':{},'MPI_ranks_per_case':4,'concurrent_cases':4,'total_compute_processes':16,'cpu_sets':{'phi000':'16-19','phi090':'20-23','phi180':'24-27','phi270':'28-31'},'aggregate_RSS_limit_GiB':8}
atomic(r/'status.json',state);jobs=[];logs=[]
def stop_own():
    for job in jobs:
        if job['proc'].poll() is None:
            try:os.killpg(job['proc'].pid,signal.SIGTERM)
            except ProcessLookupError:pass
try:
    available=int(re.search(r'MemAvailable:\s*(\d+)',pathlib.Path('/proc/meminfo').read_text())[1]);assert available>64*1024*1024,'Insufficient available memory for launch'
    assert json.loads((r/'probe-check.json').read_text())['status']=='passed'
    assert sha(r/'bin/HOS-Ocean')==json.loads((r/'build-provenance.json').read_text())['binary_sha256']
    env=os.environ.copy();env.update(json.loads((r/'environment.json').read_text()));env['GFORTRAN_UNBUFFERED_PRECONNECTED']='y'
    settings=json.loads((r/'settings.json').read_text())
    wrapper=r/'rank.sh';wrapper.write_text('#!/bin/sh\ncpu=$((HOS_CPU_BASE + OMPI_COMM_WORLD_RANK))\nprintf "HOS rank %s assigned CPU %s\\n" "$OMPI_COMM_WORLD_RANK" "$cpu"\nexec /usr/bin/taskset -c "$cpu" '+str(r/'bin/time')+' -v -o "rank_${OMPI_COMM_WORLD_RANK}.resources.txt" '+str(r/'bin/HOS-Ocean')+' input.yml\n');wrapper.chmod(0o755)
    manifest={'created_utc':utc(),'binary_sha256':sha(r/'bin/HOS-Ocean'),'initial_fields_sha256':sha(source),'initial_source':str(source),'settings':settings,'files_sha256':{str(p.relative_to(r)):sha(p) for p in [r/'run_directional_hos.py',r/'prepare_directional_hos.m',wrapper,r/'environment.json',r/'settings.json',*r.glob('cases/*/input.yml'),*r.glob('cases/*/prob.inp'),*r.glob('cases/*/Results/3d_ini_*.dat')]}}
    atomic(r/'manifest.json',manifest);state['state']='running';state['started_utc']=utc();(r/'started.utc').write_text(state['started_utc']+'\n')
    for index,phase in enumerate([0,90,180,270]):
        name='phi%03d'%phase;case=r/'cases'/name;en=env.copy();en['HOS_CPU_BASE']=str(16+4*index);log=(case/'run.log').open('w');logs.append(log)
        args=[str(mpi/'deps/usr/bin/mpirun.openmpi'),'--mca','btl','self,vader,tcp','--bind-to','none','-x','HOS_CPU_BASE','-np','4',str(wrapper)]
        proc=subprocess.Popen(args,cwd=case,env=en,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
        jobs.append({'name':name,'case':case,'proc':proc,'index':index,'start':time.monotonic(),'done':False});state['cases'][name]={'state':'running','launcher_pid':proc.pid,'initial_probe_check':'pending','started_utc':utc(),'peak_sampled_job_RSS_KiB':0,'last_logged_time_s':0.0}
    atomic(r/'status.json',state);peak=0
    while not all(j['done'] for j in jobs):
        info={}
        for p in pathlib.Path('/proc').iterdir():
            if not p.name.isdigit():continue
            try:
                st=(p/'status').read_text();pp=int(re.search(r'^PPid:\s*(\d+)',st,re.M)[1]);rss=re.search(r'^VmRSS:\s*(\d+)',st,re.M);name=re.search(r'^Name:\s*(.+)',st,re.M)[1];cpus=re.search(r'^Cpus_allowed_list:\s*(.+)',st,re.M)[1];info[int(p.name)]=(pp,int(rss[1]) if rss else 0,name,cpus)
            except (OSError,TypeError,ValueError):pass
        total=0
        for j in jobs:
            caseState=state['cases'][j['name']];selected={j['proc'].pid};changed=True
            while changed:
                add={pid for pid,(pp,_,_,_) in info.items() if pp in selected}-selected;changed=bool(add);selected|=add
            rss=sum(info[p][1] for p in selected if p in info);total+=rss;caseState['current_job_RSS_KiB']=rss;caseState['peak_sampled_job_RSS_KiB']=max(rss,caseState['peak_sampled_job_RSS_KiB']);caseState['solver_processes']=[{'pid':p,'cpus_allowed':info[p][3]} for p in selected if p in info and info[p][2]=='HOS-Ocean']
            try:
                with (j['case']/'run.log').open('rb') as f:f.seek(max(0,(j['case']/'run.log').stat().st_size-32768));tail=f.read().decode(errors='replace')
                matches=re.findall(r'CPU time for Output time step\s+\S+\s+\d+\s+\d+\s+(\S+)',tail)
                if matches:caseState['last_logged_time_s']=float(matches[-1]);caseState['progress_percent']=100*float(matches[-1])/220
            except (OSError,ValueError):pass
            probe=j['case']/'Results/probes.dat'
            if probe.exists() and caseState['initial_probe_check']=='pending':
                lines=probe.read_text(errors='replace').splitlines();idx=next((i for i,x in enumerate(lines) if x.startswith('VARIABLES')),None)
                if idx is not None and len(lines)>idx+1:
                    try:first=[float(x) for x in lines[idx+1].split()]
                    except ValueError:first=[]
                    if len(first)==6:
                        error=max(abs(a-b) for a,b in zip(first[1:],settings['expected_initial_probe_eta'][j['index']]));assert abs(first[0])<1e-10 and error<1e-10,'MPI probe input mismatch '+j['name'];caseState['initial_probe_check']='passed';caseState['initial_probe_max_abs_error_m']=error
            code=j['proc'].poll()
            if code is not None and not j['done']:
                assert code==0,'Solver failed '+j['name']+' exit '+str(code)
                lines=probe.read_text().splitlines();idx=next(i for i,x in enumerate(lines) if x.startswith('VARIABLES'));values=[[float(x) for x in line.split()] for line in lines[idx+1:] if line.strip()]
                assert len(values)==1101 and all(len(row)==6 and all(math.isfinite(x) for x in row) for row in values),'Incomplete/nonfinite probe output '+j['name']
                assert max(abs(row[0]-.2*i) for i,row in enumerate(values))<1e-8,'Unexpected probe times '+j['name']
                caseState.update(state='completed_and_probe_validated',exit_code=code,wall_seconds=time.monotonic()-j['start'],finished_utc=utc(),sample_count=len(values),last_logged_time_s=220.0,progress_percent=100.0);j['done']=True;atomic(j['case']/'completion.json',caseState)
        peak=max(peak,total);state.update(updated_utc=utc(),current_aggregate_job_RSS_KiB=total,peak_sampled_aggregate_job_RSS_KiB=peak)
        assert total<8*1024*1024,'Aggregate HOS RSS exceeded 8 GiB budget'
        atomic(r/'status.json',state)
        if not all(j['done'] for j in jobs):time.sleep(5)
    state.update(state='completed',finished_utc=utc());atomic(r/'status.json',state);(r/'finished.utc').write_text(state['finished_utc']+'\n')
except Exception as e:
    stop_own();state.update(state='failed',error=str(e),updated_utc=utc());atomic(r/'status.json',state);raise
finally:
    for log in logs:log.close()