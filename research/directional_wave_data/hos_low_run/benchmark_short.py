import pathlib,subprocess,os,json,time,hashlib,re,signal
r=pathlib.Path('/home/lxy/green-laplace-unidirectional-time-series-runs/hos-directional-akp002-20260926-v1');mpi=pathlib.Path('/home/lxy/green-laplace-unidirectional-time-series-runs/hos-mpi-check-20260926-v1');env=os.environ.copy();env.update(json.loads((r/'environment.json').read_text()));env['HOS_CPU_BASE']='8';env['GFORTRAN_UNBUFFERED_PRECONNECTED']='y'
stage=r/'staging4/cases/phi000';data=[]
for rank in range(4):data+=(stage/'Results'/('3d_ini_%03d.dat'%rank)).read_text().splitlines()[68:]
assert len(data)==1024*256
wrapper=r/'rank.sh';wrapper.write_text('#!/bin/sh\ncpu=$((HOS_CPU_BASE + OMPI_COMM_WORLD_RANK))\nexec /usr/bin/taskset -c "$cpu" '+str(r/'bin/time')+' -v -o "rank_${OMPI_COMM_WORLD_RANK}.resources.txt" '+str(r/'bin/HOS-Ocean')+' input.yml\n');wrapper.chmod(0o755)
summary={'physical_duration_s':2,'production_duration_s':220,'runs':[],'state':'running_short_benchmarks'}
def save(): (r/'benchmark.json').write_text(json.dumps(summary,indent=2))
def memory(pid):
 info={}
 for p in pathlib.Path('/proc').iterdir():
  if not p.name.isdigit():continue
  try:
   s=(p/'status').read_text();pp=int(re.search(r'^PPid:\s*(\d+)',s,re.M)[1]);m=re.search(r'^VmRSS:\s*(\d+)',s,re.M);info[int(p.name)]=(pp,int(m[1]) if m else 0)
  except (OSError,TypeError,ValueError):pass
 selected={pid}
 while True:
  add={p for p,(pp,_) in info.items() if pp in selected}-selected
  if not add:break
  selected|=add
 return sum(info[p][1] for p in selected if p in info)
for np in [4,8]:
 case=r/('bench_np%d'%np);(case/'Results').mkdir(parents=True);ny=256//np
 for rank in range(np):
  header=['# frozen Akp=.02 MF12 order-2 input']*67+['%-20s%12.5E%4s%5d%4s%5d'%('ZONE SOLUTIONTIME = ',0,', I=',1024,', J=',ny)]
  (case/'Results'/('3d_ini_%03d.dat'%rank)).write_text('\n'.join(header+data[rank*ny*1024:(rank+1)*ny*1024])+'\n')
 text=(stage/'input.yml').read_text();assert 'duration: 220.0' in text;text=text.replace('duration: 220.0','duration: 2.0').replace('physical space: false','physical space: true');assert 'duration: 2.0' in text and 'duration: 220' not in text;(case/'input.yml').write_text(text);(case/'prob.inp').write_bytes((stage/'prob.inp').read_bytes())
 save();start=time.monotonic();peak=0
 with (case/'run.log').open('w') as f:
  p=subprocess.Popen([str(mpi/'deps/usr/bin/mpirun.openmpi'),'--mca','btl','self,vader,tcp','--bind-to','none','-x','HOS_CPU_BASE','-np',str(np),str(wrapper)],cwd=case,env=env,stdout=f,stderr=subprocess.STDOUT,start_new_session=True)
  while p.poll() is None:
   rss=memory(p.pid);peak=max(peak,rss)
   if rss>8*1024*1024:os.killpg(p.pid,signal.SIGTERM);raise RuntimeError('Benchmark exceeded 8 GiB RSS')
   time.sleep(.5)
 report={'ranks':np,'wall_seconds':time.monotonic()-start,'exit_code':p.returncode,'peak_sampled_job_RSS_KiB':peak,'physical_duration_s':2,'input_sha256':hashlib.sha256((case/'input.yml').read_bytes()).hexdigest()};summary['runs'].append(report);save();print(json.dumps(report),flush=True);assert p.returncode==0,(case/'run.log').read_text()[-2500:]
summary['state']='awaiting_result_parity';save()