import pathlib,json,subprocess,os,shutil,time,csv,hashlib
r=pathlib.Path('/home/lxy/green-laplace-unidirectional-time-series-runs/hos-directional-akp002-20260926-v1');high=pathlib.Path('/home/lxy/green-laplace-unidirectional-time-series-runs/hos-directional-fourphase-20260926-v2');env=os.environ.copy();env.update({k:'1' for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']})
def utc():return time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())
def status(stage,**extra):
 p=r/'pipeline-status.json';d={'stage':stage,'updated_utc':utc(),**extra};tmp=p.with_suffix('.tmp');tmp.write_text(json.dumps(d,indent=2));tmp.replace(p)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def full_gl(run):
 out=run/'full-gl-comparison-20260926-v1';assert not out.exists();out.mkdir();files={}
 for phase in [0,90,180,270]:
  p=run/'cases'/('phi%03d'%phase)/'Results/probes.dat';lines=p.read_text().splitlines();i=next(i for i,x in enumerate(lines) if x.startswith('VARIABLES'));rows=[[float(v) for v in x.split()] for x in lines[i+1:] if x.strip()];assert len(rows)==1101 and all(len(x)==6 for x in rows);assert max(abs(x[0]-.2*i) for i,x in enumerate(rows))<1e-8
  target=out/('phi%03d.csv'%phase)
  with target.open('w',newline='') as f:csv.writer(f).writerows(rows)
  files[target.name]=sha(target)
 (out/'snapshot.json').write_text(json.dumps({'sample_count':1101,'end_time_s':220,'partial':False,'source_run':str(run),'frozen_utc':utc(),'files_sha256':files},indent=2));(out/'source').mkdir();subprocess.run(['tar','-xzf',str(r/'directional-preview-source.tar.gz'),'-C',str(out/'source')],check=True)
 command="addpath('%s');compare_directional_full('%s');"%(r,run)
 with (out/'run.log').open('w') as f:p=subprocess.run(['/usr/bin/taskset','-c','40',str(r/'bin/time'),'-v','-o',str(out/'resources.txt'),'/home/lxy/Desktop/matlabr2026a/bin/matlab','-singleCompThread','-batch',command],cwd=run,env=env,stdout=f,stderr=subprocess.STDOUT)
 assert p.returncode==0,(out/'run.log').read_text()[-2000:]
 return str(out)
try:
 b=json.loads((r/'benchmark.json').read_text());v=json.loads((r/'benchmark-validation.json').read_text());assert v['status']=='PASSED_2_SECOND_MPI4_MPI8_AND_PROBE_CHECK';assert b['physical_duration_s']==2 and len(b['runs'])==2
 times={x['ranks']:x for x in b['runs']};np=8 if times[8]['wall_seconds']<times[4]['wall_seconds'] and 4*times[8]['peak_sampled_job_RSS_KiB']<8*1024*1024 else 4
 selection={'chosen_ranks':np,'reason':'Lower measured wall-clock time under memory/CPU budget; total CPU hours not selection criterion','benchmark_physical_duration_s':2,'benchmark':b,'max_HOS_compute_slots':4*np,'expected_retained_OW3D':8,'automatic_GL_after_HOS':True};(r/'selection.json').write_text(json.dumps(selection,indent=2));status('preparing_selected_production',selection=selection)
 settings=json.loads((r/'staging4/settings.json').read_text());assert abs(settings['Akp']-.02)<1e-12;settings['MPI_ranks_per_phase']=np;settings['short_benchmark_duration_s']=2;(r/'settings.json').write_text(json.dumps(settings,indent=2));ny=256//np
 for phase in [0,90,180,270]:
  src=r/'staging4/cases'/('phi%03d'%phase);case=r/'cases'/('phi%03d'%phase);(case/'Results').mkdir(parents=True);data=[]
  for rank in range(4):data+=(src/'Results'/('3d_ini_%03d.dat'%rank)).read_text().splitlines()[68:]
  assert len(data)==1024*256
  for rank in range(np):
   header=['# Fresh MF12 Akp=.02, same geometry; production 220 seconds']*67+['%-20s%12.5E%4s%5d%4s%5d'%('ZONE SOLUTIONTIME = ',0,', I=',1024,', J=',ny)];(case/'Results'/('3d_ini_%03d.dat'%rank)).write_text('\n'.join(header+data[rank*ny*1024:(rank+1)*ny*1024])+'\n')
  shutil.copy2(src/'input.yml',case/'input.yml');shutil.copy2(src/'prob.inp',case/'prob.inp');assert 'duration: 220.0' in (case/'input.yml').read_text()
 status('HOS_production_running',selected_MPI_ranks=np,GL_status='queued_after_all_four_phases')
 with (r/'production-controller.log').open('w') as f:p=subprocess.run(['python3',str(r/'run_directional_hos.py')],cwd=r,stdout=f,stderr=subprocess.STDOUT)
 assert p.returncode==0 and json.loads((r/'status.json').read_text())['state']=='completed',(r/'production-controller.log').read_text()[-2500:]
 status('GL_low_processing');lowout=full_gl(r);status('GL_high_full_processing',low_result=lowout);highout=full_gl(high)
 status('completed_HOS_and_GL',low_result=lowout,high_result=highout,selected_MPI_ranks=np)
except Exception as e:status('failed',error=str(e));raise