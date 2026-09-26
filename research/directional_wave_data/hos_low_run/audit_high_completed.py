import pathlib,json,re,datetime
r=pathlib.Path('/home/lxy/green-laplace-unidirectional-time-series-runs/hos-directional-fourphase-20260926-v2');s=json.loads((r/'status.json').read_text());assert s['state']=='completed';rows=[]
for phase in [0,90,180,270]:
 d=r/'cases'/('phi%03d'%phase);lines=(d/'Results/probes.dat').read_text().splitlines();i=next(i for i,x in enumerate(lines) if x.startswith('VARIABLES'));v=[[float(y) for y in x.split()] for x in lines[i+1:] if x.strip()];assert len(v)==1101 and abs(v[-1][0]-220)<1e-8
 wall=[]
 for p in d.glob('rank_*.resources.txt'):
  text=p.read_text();value=re.search(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*(\S+)',text)[1];sec=0
  for q in value.split(':'):sec=60*sec+float(q)
  wall.append(sec)
 rows.append({'phase':phase,'sample_count':len(v),'last_physical_time_s':v[-1][0],'max_rank_wall_seconds':max(wall),'rank_count':len(wall)})
start=datetime.datetime.fromisoformat(s['started_utc'].replace('Z','+00:00'));end=datetime.datetime.fromisoformat(s['finished_utc'].replace('Z','+00:00'))
a={'all_four_completed':True,'batch_wait_seconds':(end-start).total_seconds(),'started_utc':s['started_utc'],'finished_utc':s['finished_utc'],'cases':rows,'peak_sampled_aggregate_RSS_KiB':s['peak_sampled_aggregate_job_RSS_KiB'],'note':'Native records and completion states are authoritative; historical last_logged_time fields may lag after a case exits.'}
out=pathlib.Path('/home/lxy/green-laplace-unidirectional-time-series-runs/hos-directional-akp002-20260926-v1/high-completion-audit.json');out.write_text(json.dumps(a,indent=2));print(json.dumps(a,indent=2))