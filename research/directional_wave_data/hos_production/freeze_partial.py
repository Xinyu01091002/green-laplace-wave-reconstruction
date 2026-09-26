import pathlib,subprocess,json,hashlib,math,time,csv
run=pathlib.Path('/home/lxy/green-laplace-unidirectional-time-series-runs/hos-directional-fourphase-20260926-v2');out=run/'partial-comparison-20260926-v1';assert not out.exists();out.mkdir();records=[]
for phase in [0,90,180,270]:
 p=run/'cases'/('phi%03d'%phase)/'Results/probes.dat';lines=p.read_text().splitlines();i=next(i for i,x in enumerate(lines) if x.startswith('VARIABLES'));rows=[]
 for line in lines[i+1:]:
  try:v=[float(x) for x in line.split()]
  except ValueError:continue
  if len(v)==6:assert all(math.isfinite(x) for x in v);rows.append(v)
 records.append(rows)
n=min(map(len,records));n-=1-n%2
assert (n-1)*.2>150
files={}
for phase,rows in zip([0,90,180,270],records):
 assert max(abs(row[0]-.2*i) for i,row in enumerate(rows[:n]))<1e-8
 p=out/('phi%03d.csv'%phase)
 with p.open('w',newline='') as f:csv.writer(f).writerows(rows[:n])
 files[p.name]=hashlib.sha256(p.read_bytes()).hexdigest()
manifest={'frozen_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'sample_count':n,'end_time_s':.2*(n-1),'source_counts':list(map(len,records)),'source_run':str(run),'files_sha256':files,'partial':True,'solver_not_paused':True}
(out/'snapshot.json').write_text(json.dumps(manifest,indent=2));print(json.dumps(manifest,indent=2))