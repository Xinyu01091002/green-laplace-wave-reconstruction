"""Run the frozen Akp=.12 four-phase HOS replay at off-axis probes."""
import hashlib
import json
import math
import os
import pathlib
import re
import signal
import subprocess
import time

ROOT=pathlib.Path('/home/lxy/green-laplace-unidirectional-time-series-runs/hos-directional-high-xy-probes-20260928-v1')
SOURCE=pathlib.Path('/home/lxy/green-laplace-unidirectional-time-series-runs/hos-directional-fourphase-20260926-v2')
MPI=pathlib.Path('/home/lxy/green-laplace-unidirectional-time-series-runs/hos-mpi-check-20260926-v1')

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def utc():return time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())
def atomic(path,value):
    temporary=path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value,indent=2));temporary.replace(path)

settings=json.loads((ROOT/'settings.json').read_text())
binary=ROOT/'bin/HOS-Ocean'
assert sha(binary)==json.loads((SOURCE/'build-provenance.json').read_text())['binary_sha256']
environment=os.environ.copy();environment.update(json.loads((SOURCE/'environment.json').read_text()))
environment['GFORTRAN_UNBUFFERED_PRECONNECTED']='y'
wrapper=ROOT/'rank.sh'
wrapper.write_text('#!/bin/sh\nrank=$((OMPI_COMM_WORLD_RANK + 1))\ncpu=$(printf "%s" "$HOS_CPU_LIST" | cut -d, -f "$rank")\nprintf "HOS rank %s assigned CPU %s\\n" "$OMPI_COMM_WORLD_RANK" "$cpu"\nexec /usr/bin/taskset -c "$cpu" '+str(ROOT/'bin/time')+' -v -o "rank_${OMPI_COMM_WORLD_RANK}.resources.txt" '+str(binary)+' input.yml\n')
wrapper.chmod(0o755)
state={'state':'running','created_utc':utc(),'cases':{},'cpu_sets':settings['cpu_sets'],'peak_sampled_aggregate_job_RSS_KiB':0}
jobs=[];logs=[]
try:
    for index,phase in enumerate([0,90,180,270]):
        name=f'phi{phase:03d}';case=ROOT/'cases'/name
        cpu_range=settings['cpu_sets'][index]
        first,last=(int(x) for x in cpu_range.split('-'))
        env=environment.copy();env['HOS_CPU_LIST']=','.join(str(x) for x in range(first,last+1))
        log=(case/'run.log').open('w');logs.append(log)
        args=[str(MPI/'deps/usr/bin/mpirun.openmpi'),'--mca','btl','self,vader,tcp','--bind-to','none','-x','HOS_CPU_LIST','-np','4',str(wrapper)]
        process=subprocess.Popen(args,cwd=case,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
        jobs.append({'name':name,'case':case,'process':process,'start':time.monotonic(),'done':False,'index':index})
        state['cases'][name]={'state':'running','launcher_pid':process.pid,'cpu_set':cpu_range,'initial_probe_check':'pending','peak_sampled_job_RSS_KiB':0,'last_logged_time_s':0.0}
    atomic(ROOT/'status.json',state)
    while not all(job['done'] for job in jobs):
        process_info={}
        for path in pathlib.Path('/proc').iterdir():
            if not path.name.isdigit():continue
            try:
                status=(path/'status').read_text();ppid=int(re.search(r'^PPid:\s*(\d+)',status,re.M)[1]);rss=re.search(r'^VmRSS:\s*(\d+)',status,re.M)
                process_info[int(path.name)]=(ppid,int(rss[1]) if rss else 0)
            except (OSError,TypeError,ValueError):pass
        aggregate=0
        for job in jobs:
            selected={job['process'].pid};changed=True
            while changed:
                additions={pid for pid,(parent,_) in process_info.items() if parent in selected}-selected
                changed=bool(additions);selected|=additions
            rss=sum(process_info[pid][1] for pid in selected if pid in process_info);aggregate+=rss
            case_state=state['cases'][job['name']]
            case_state['peak_sampled_job_RSS_KiB']=max(rss,case_state['peak_sampled_job_RSS_KiB'])
            probe=job['case']/'Results/probes.dat'
            if probe.exists() and case_state['initial_probe_check']=='pending':
                lines=probe.read_text(errors='replace').splitlines();header=next((i for i,line in enumerate(lines) if line.startswith('VARIABLES')),None)
                if header is not None and len(lines)>header+1:
                    try:first=[float(x) for x in lines[header+1].split()]
                    except ValueError:first=[]
                    expected=settings['expected_initial_probe_eta'][job['index']]
                    if len(first)==1+len(expected):
                        error=max(abs(a-b) for a,b in zip(first[1:],expected))
                        assert abs(first[0])<1e-10 and error<1e-10
                        case_state['initial_probe_check']='passed';case_state['initial_probe_max_abs_error_m']=error
            try:
                with (job['case']/'run.log').open('rb') as stream:
                    stream.seek(max(0,(job['case']/'run.log').stat().st_size-32768));tail=stream.read().decode(errors='replace')
                matches=re.findall(r'CPU time for Output time step\s+\S+\s+\d+\s+\d+\s+(\S+)',tail)
                if matches:case_state['last_logged_time_s']=float(matches[-1]);case_state['progress_percent']=100*float(matches[-1])/220
            except (OSError,ValueError):pass
            code=job['process'].poll()
            if code is not None and not job['done']:
                assert code==0,f"Solver failed {job['name']} exit {code}"
                lines=probe.read_text().splitlines();header=next(i for i,line in enumerate(lines) if line.startswith('VARIABLES'))
                values=[[float(x) for x in line.split()] for line in lines[header+1:] if line.strip()]
                assert len(values)==1101 and all(len(row)==13 and all(math.isfinite(x) for x in row) for row in values)
                assert max(abs(row[0]-.2*i) for i,row in enumerate(values))<1e-8
                case_state.update(state='completed_and_probe_validated',exit_code=code,wall_seconds=time.monotonic()-job['start'],sample_count=len(values),finished_utc=utc(),progress_percent=100.0,last_logged_time_s=220.0)
                atomic(job['case']/'completion.json',case_state);job['done']=True
        state['peak_sampled_aggregate_job_RSS_KiB']=max(aggregate,state['peak_sampled_aggregate_job_RSS_KiB'])
        state['updated_utc']=utc();atomic(ROOT/'status.json',state)
        assert aggregate<8*1024*1024
        if not all(job['done'] for job in jobs):time.sleep(5)
    state['state']='completed';state['finished_utc']=utc();atomic(ROOT/'status.json',state)
except Exception as error:
    for job in jobs:
        if job['process'].poll() is None:
            try:os.killpg(job['process'].pid,signal.SIGTERM)
            except ProcessLookupError:pass
    state['state']='failed';state['error']=str(error);state['updated_utc']=utc();atomic(ROOT/'status.json',state)
    raise
finally:
    for log in logs:log.close()
