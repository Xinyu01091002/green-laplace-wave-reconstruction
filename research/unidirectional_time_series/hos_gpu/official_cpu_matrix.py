"""Run isolated official MPI benchmarks on 92, with measured loop wall time/RSS."""
import argparse,json,os,pathlib,re,shutil,subprocess,time
parser=argparse.ArgumentParser();parser.add_argument('--runtime',required=True,type=pathlib.Path);parser.add_argument('--out',required=True,type=pathlib.Path)
parser.add_argument('--exe',type=pathlib.Path);parser.add_argument('--duration',default=2.,type=float);parser.add_argument('--family',choices=['low','high','both'],default='both');parser.add_argument('--ranks',type=int,choices=[8,16,32],default=8);args=parser.parse_args()
r=args.runtime.resolve();out=args.out.resolve();out.mkdir(exist_ok=False);u=r/'mpi/usr';lib=u/'lib/x86_64-linux-gnu'
env=dict(os.environ,OPAL_PREFIX=str(u),PMIX_MCA_mca_base_component_path=str(lib/'pmix2/lib/pmix'),LD_LIBRARY_PATH=':'.join(map(str,[r/'support',lib,lib/'openmpi/lib'])),
         OMPI_ALLOW_RUN_AS_ROOT='1',OMPI_ALLOW_RUN_AS_ROOT_CONFIRM='1',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMPI_MCA_pml='ob1',OMPI_MCA_btl='self,vader,tcp')
exe=(args.exe or r/'bin/HOS-benchmark').resolve();env['HOS_BENCH_EXE']=str(exe)
def rss_tree(pid):
    data={}
    for p in pathlib.Path('/proc').iterdir():
        if not p.name.isdigit():continue
        try:
            s=(p/'status').read_text();parent=re.search(r'^PPid:\s*(\d+)',s,re.M);rss=re.search(r'^VmRSS:\s*(\d+)',s,re.M)
            data[int(p.name)]=(int(parent[1]),int(rss[1]) if rss else 0)
        except (OSError,TypeError):pass
    selected={pid}
    while True:
        extra={p for p,(pp,_) in data.items() if pp in selected}-selected
        if not extra:break
        selected|=extra
    return sum(data[p][1] for p in selected if p in data)
def rows(path):
    text=path.read_text().splitlines();i=next(i for i,s in enumerate(text) if 'VARIABLES' in s)
    return [s.split() for s in text[i+1:] if s.strip()]
results=[]
for family in (['low','high'] if args.family=='both' else [args.family]):
    case=out/family
    if args.ranks==8:shutil.copytree(r/'cases'/family,case)
    else:
        case.mkdir();(case/'Results').mkdir()
        for name in ['input.yml','prob.inp']:shutil.copy2(r/'cases'/family/name,case/name)
        data=[];header=None
        for path in sorted((r/'cases'/family/'Results').glob('3d_ini_*.dat')):
            lines=path.read_text().splitlines(keepends=True)
            if header is None:header=lines[:67]
            data.extend(lines[68:])
        assert len(data)==1024*512;ny=512//args.ranks;count=1024*ny
        for rank in range(args.ranks):
            with (case/'Results'/('3d_ini_%03d.dat'%rank)).open('w') as f:
                f.writelines(header);f.write(f'{"ZONE SOLUTIONTIME = ":20s}{0:12.5E}{", I=":4s}{1024:5d}{", J=":4s}{ny:5d}\n');f.writelines(data[rank*count:(rank+1)*count])
    p=case/'input.yml';s=p.read_text();s,n=re.subn(r'(?m)^  duration: .*$',f'  duration: {args.duration:.17g}',s);assert n==1;p.write_text(s)
    command=[str(u/'bin/mpirun.openmpi'),'--allow-run-as-root','--oversubscribe','--bind-to','none','-np',str(args.ranks),str(r/'rank.sh')]
    start=time.monotonic();last=start;peak=0
    print('Starting official MPI'+str(args.ranks)+' '+family+' duration='+str(args.duration),flush=True)
    with (case/'run.log').open('w') as log:
        process=subprocess.Popen(command,cwd=case,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
        while process.poll() is None:
            peak=max(peak,rss_tree(process.pid))
            if time.monotonic()-last>=30:
                tail=(case/'run.log').read_text(errors='replace').splitlines()[-2:]
                print(json.dumps({'family':family,'elapsed':time.monotonic()-start,'peak_rss_kib':peak,'tail':tail}),flush=True);last=time.monotonic()
            if time.monotonic()-start>3600:
                os.killpg(process.pid,15);process.wait();raise TimeoutError('Own benchmark exceeded budget')
            time.sleep(.1)
    log=(case/'run.log').read_text();match=re.search(r'BENCH_ADVANCE_SECONDS\s+([\d.E+-]+)',log);counts=re.search(r'BENCH_STEPS\s+(\d+)\s+(\d+)',log)
    record={'family':family,'duration_s':args.duration,'ranks':args.ranks,'exit':process.returncode,'whole_process_seconds':time.monotonic()-start,'aggregate_peak_rss_kib':peak,
            'advance_seconds':float(match[1]) if match else None,'accepted':int(counts[1]) if counts else None,'rejected':int(counts[2]) if counts else None}
    if (case/'Results/probes.dat').exists() and 'VARIABLES' in (case/'Results/probes.dat').read_text():
        probes=rows(case/'Results/probes.dat');energy={round(float(a[0]),6):a for a in rows(case/'Results/vol_energy.dat')}
        with (case/'diagnostics.csv').open('w') as f:
            f.write('time,p1,p2,p3,p4,p5,potential,kinetic,total,relative_energy_change\n')
            for a in probes:
                e=energy.get(round(float(a[0]),6))
                if e is not None and len(a)>=6 and len(e)>=6:f.write(','.join(a[:6]+[e[2],e[3],e[4],e[5]])+'\n')
        record['last_output_time_s']=float(probes[-1][0]);record['output_rows']=len(probes)
    results.append(record);(out/'summary.json').write_text(json.dumps(results,indent=2));print(json.dumps(record),flush=True)
