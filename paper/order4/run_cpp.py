"""Process orchestration only; all numerical reconstruction is in C++/MATLAB."""
from pathlib import Path
import argparse,concurrent.futures,hashlib,json,os,platform,subprocess,time

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('method',choices=['gl','wit'])
    parser.add_argument('input',type=Path)
    parser.add_argument('output',type=Path)
    parser.add_argument('--build',type=Path,default=Path('artifacts/order4-build'))
    parser.add_argument('--rank',type=int,choices=[6,8,10],default=6)
    parser.add_argument('--threads',type=int,default=1)
    parser.add_argument('--split-workers',type=int,default=0)
    args=parser.parse_args()
    if args.threads<1 or args.split_workers<0:parser.error('Invalid thread/worker count')
    if args.split_workers and (args.method!='gl' or args.threads!=1):parser.error('Split GL uses one thread per process')
    args.output.parent.mkdir(parents=True,exist_ok=True)
    binary=(args.build/('gl_eta44' if args.method=='gl' else 'wit_eta44')).resolve()
    env=dict(os.environ,OMP_NUM_THREADS=str(args.threads))
    started=time.perf_counter()
    def run(cmd):subprocess.run(list(map(str,cmd)),env=env,check=True)
    if args.split_workers:
        parts=[];commands=[]
        folder=args.output.parent/(args.output.stem+'_parts');folder.mkdir(exist_ok=True)
        for outer in range(args.rank):
            for part in range(-1,args.rank):
                p=folder/f'outer{outer}_part{part}.bin';parts.append(p)
                commands.append([binary,args.input,p,args.rank,1,outer,part])
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.split_workers) as pool:
            list(pool.map(run,commands))
        run([(args.build/'merge_eta44_outer').resolve(),args.output,*parts])
    elif args.method=='gl':
        run([binary,args.input,args.output,args.rank,args.threads])
    else:
        run([binary,args.input,args.output,'--sector','44'])
    elapsed=time.perf_counter()-started
    report={}
    report.update(method=args.method,rank=args.rank if args.method=='gl' else None,
                  threads_per_process=args.threads,split_workers=args.split_workers,
                  process_wall_seconds=elapsed,timing_statistic='single measured run',
                  warmup_runs=0,precision='double',host=platform.platform(),
                  cpu_count=os.cpu_count(),input_sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(),
                  executable_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
                  timing_boundary='process launch through completion; includes serialization and split merge',
                  stokes_correction_used=False)
    args.output.with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report))

if __name__=='__main__':main()
