import pathlib,re,json,csv,tarfile,hashlib
r=pathlib.Path('/home/lxy/green-laplace-unidirectional-time-series-runs/hos-gl-timeseries-20260926-v1');state=json.loads((r/'status.json').read_text());assert state['state']=='completed'
rows=[]
for p in sorted(r.rglob('*-resources.txt')):
    s=p.read_text();get=lambda pat:re.search(pat,s).group(1)
    wall=0
    for v in get(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*(\S+)').split(':'):wall=60*wall+float(v)
    rows.append({'file':str(p.relative_to(r)),'wall_seconds':wall,'cpu_seconds':float(get(r'User time \(seconds\):\s*(\S+)'))+float(get(r'System time \(seconds\):\s*(\S+)')),'peak_RSS_KiB':int(get(r'Maximum resident set size \(kbytes\):\s*(\d+)')),'exit_code':int(get(r'Exit status:\s*(\d+)'))})
with (r/'performance.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
hos=[x for x in rows if x['file'].endswith('/hos-resources.txt')];summary={'state':'completed','HOS_processes':len(hos),'sum_HOS_process_wall_seconds':sum(x['wall_seconds'] for x in hos),'sum_HOS_cpu_seconds':sum(x['cpu_seconds'] for x in hos),'performance':rows,'metrics':{},'source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in r.glob('*.m')}}
for name in ['akp002','akp012']:
    metrics=list(csv.DictReader((r/name/'metrics.csv').open()));summary['metrics'][name]=metrics
    print(name)
    for row in metrics:
        if row['window']=='main_group' and row['method'] in ['GL6','GL12','MF12']:print(row)
summary['cumulative_HOS_process_wall_seconds_including_previous_33_runs']=2168.31+summary['sum_HOS_process_wall_seconds']
(r/'summary.json').write_text(json.dumps(summary,indent=2));print('PERFORMANCE',json.dumps(rows));print('HOS TOTAL',summary['sum_HOS_process_wall_seconds'])
files=[r/'summary.json',r/'performance.csv',r/'status.json',r/'build.json',r/'hos-io.patch',r/'probe-smoke/io-check.json']
for name in ['akp002','akp012']:
    files.extend([r/name/x for x in ['comparison.png','comparison.pdf','metrics.csv','report.json','initialization.json']])
with tarfile.open(r/'compact-results.tar.gz','w:gz') as tar:
    for p in files:tar.add(p,arcname=str(p.relative_to(r)))