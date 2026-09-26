import pathlib,re,json,csv
root=pathlib.Path('/home/lxy/green-laplace-unidirectional-time-series-runs/hos-deepwater-fourphase-20260926-v1');base=pathlib.Path('/home/lxy/green-laplace-unidirectional-time-series/artifacts/hos-ocean-v2.1.0')
rows=[]
for scope,files in [('deepwater',sorted(root.rglob('hos-resources.txt'))),('earlier_pilots',sorted(base.rglob('resources.txt')))]:
    for p in files:
        s=p.read_text()
        if 'HOS-Ocean' not in s:continue
        def field(pattern):return re.search(pattern,s).group(1)
        clock=field(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*(\S+)');wall=0
        for part in clock.split(':'):wall=wall*60+float(part)
        user=float(field(r'User time \(seconds\):\s*(\S+)'));system=float(field(r'System time \(seconds\):\s*(\S+)'))
        rows.append(dict(scope=scope,file=str(p),wall_seconds=wall,cpu_seconds=user+system,peak_RSS_KiB=int(field(r'Maximum resident set size \(kbytes\):\s*(\d+)'))))
summary={}
for name in ['deepwater','earlier_pilots','all']:
    selected=rows if name=='all' else [x for x in rows if x['scope']==name]
    summary[name]={'HOS_processes':len(selected),'sum_process_wall_seconds':sum(x['wall_seconds'] for x in selected),'sum_user_plus_system_cpu_seconds':sum(x['cpu_seconds'] for x in selected)}
batches=[['akp002'],['akp012'],['tol12/akp002','tol14/akp002'],['tol12/akp012','tol14/akp012'],['tol16/akp002']]
summary['deepwater']['approx_active_batch_wait_seconds']=sum(max(x['wall_seconds'] for x in rows if x['scope']=='deepwater' and any(str(pathlib.Path(x['file']).relative_to(root)).startswith(prefix+'/') for prefix in batch)) for batch in batches)
(root/'hos-total-runtime.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
with (root/'hos-runtime-inventory.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)