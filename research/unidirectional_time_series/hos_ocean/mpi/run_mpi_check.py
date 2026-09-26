import pathlib,subprocess,os,json,time,hashlib,re,signal
r=pathlib.Path('/home/lxy/green-laplace-unidirectional-time-series-runs/hos-mpi-check-20260926-v1');b=pathlib.Path('/home/lxy/green-laplace-unidirectional-time-series/artifacts/hos-ocean-v2.1.0');reference=b/'precision-v1/case-precise-full';env=os.environ.copy();env.update(json.loads((r/'environment.json').read_text()))
source=reference/'Results/3d_ini.dat';lines=source.read_text().splitlines();assert len(lines)==68+1024*256
for np in [1,4]:
    case=r/('np%d'%np);case.mkdir();(case/'Results').mkdir();(case/'input.yml').write_text((reference/'input_HOS-Ocean.yml').read_text());ny=256//np
    for rank in range(np):
        header=lines[:67]+['%-20s%12.5E%4s%5d%4s%5d'%('ZONE SOLUTIONTIME = ',0,', I=',1024,', J=',ny)]
        data=lines[68+rank*ny*1024:68+(rank+1)*ny*1024]
        (case/'Results'/('3d_ini_%03d.dat'%rank)).write_text('\n'.join(header+data)+'\n')
wrapper=r/'rank.sh';wrapper.write_text('#!/bin/sh\nexec '+str(b/'deps/usr/bin/time')+' -v -o "rank_${OMPI_COMM_WORLD_RANK}.resources.txt" '+str(r/'build/sources/HOS-Ocean')+' input.yml\n');wrapper.chmod(0o755)
state={'state':'prepared','max_compute_processes':4,'source_input_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'binary_sha256':hashlib.sha256((r/'build/sources/HOS-Ocean').read_bytes()).hexdigest(),'runs':[]}
def save(): (r/'status.json').write_text(json.dumps(state,indent=2))
def snapshot(leader):
    info={}
    for p in pathlib.Path('/proc').iterdir():
        if not p.name.isdigit():continue
        try:
            s=(p/'status').read_text();ppid=int(re.search(r'^PPid:\s*(\d+)',s,re.M)[1]);rss=re.search(r'^VmRSS:\s*(\d+)',s,re.M);comm=re.search(r'^Name:\s*(.+)',s,re.M)[1];info[int(p.name)]=(ppid,int(rss[1]) if rss else 0,comm)
        except (OSError,ValueError,TypeError):pass
    selected={leader};changed=True
    while changed:
        add={pid for pid,(pp,_,_) in info.items() if pp in selected}-selected;changed=bool(add);selected|=add
    total=sum(info[p][1] for p in selected if p in info);hos=sum(info[p][1] for p in selected if p in info and info[p][2]=='HOS-Ocean');pss=0
    for pid in selected:
        try:pss+=int(re.search(r'^Pss:\s*(\d+)',pathlib.Path('/proc',str(pid),'smaps_rollup').read_text(),re.M)[1])
        except (OSError,TypeError):pass
    return total,hos,pss
try:
    save()
    for np in [1,4]:
        case=r/('np%d'%np);state['state']='running_np%d'%np;save();start=time.time();samples=[]
        args=[str(r/'deps/usr/bin/mpirun.openmpi'),'--mca','btl','self,vader,tcp','--bind-to','core','--map-by','core','--report-bindings','-np',str(np),str(wrapper)]
        with (case/'run.log').open('w') as log:
            proc=subprocess.Popen(args,cwd=case,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
            while proc.poll() is None:
                a,c,p=snapshot(proc.pid);samples.append([time.time()-start,a,c,p]);
                if a>8*1024*1024:os.killpg(proc.pid,signal.SIGTERM);raise RuntimeError('Own MPI job exceeded 8 GiB sampled aggregate RSS budget')
                time.sleep(.25)
        wall=time.time()-start
        report={'mpi_processes':np,'exit_code':proc.returncode,'wall_seconds':wall,'peak_sampled_job_RSS_KiB':max(x[1] for x in samples),'peak_sampled_HOS_RSS_sum_KiB':max(x[2] for x in samples),'peak_sampled_job_PSS_KiB':max(x[3] for x in samples)}
        (case/'memory_samples.json').write_text(json.dumps(samples));(case/'runtime.json').write_text(json.dumps(report,indent=2));state['runs'].append(report);save();print(json.dumps(report),flush=True)
        assert proc.returncode==0,(case/'run.log').read_text()[-3000:]
    state['state']='solver_runs_completed_pending_validation';save()
except Exception as e:state['state']='failed';state['error']=str(e);save();raise