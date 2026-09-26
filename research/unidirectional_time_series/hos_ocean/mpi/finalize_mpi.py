import pathlib,json,re,hashlib
r=pathlib.Path('/home/lxy/green-laplace-unidirectional-time-series-runs/hos-mpi-check-20260926-v1');v=json.loads((r/'validation.json').read_text());assert v['status']=='MPI_1_AND_4_AGREE_WITH_SERIAL';s=json.loads((r/'status.json').read_text());assert all(x['exit_code']==0 for x in s['runs']);s['state']='completed_and_validated';(r/'status.json').write_text(json.dumps(s,indent=2))
for a in s['runs']:
    np=a['mpi_processes'];timers=[]
    for p in sorted((r/('np%d'%np)).glob('rank_*.resources.txt')):
        text=p.read_text();get=lambda pat:float(re.search(pat,text)[1]);timers.append({'rank_file':p.name,'peak_RSS_KiB':int(get(r'Maximum resident set size \(kbytes\):\s*(\d+)')),'cpu_seconds':get(r'User time \(seconds\):\s*(\S+)')+get(r'System time \(seconds\):\s*(\S+)')})
    a['ranks']=timers;a['sum_rank_cpu_seconds']=sum(x['cpu_seconds'] for x in timers);a['sum_rank_peak_RSS_KiB_upper_bound']=sum(x['peak_RSS_KiB'] for x in timers)
summary={'runs':s['runs'],'speedup_np4_vs_np1':s['runs'][0]['wall_seconds']/s['runs'][1]['wall_seconds'],'validation':v,'hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [*r.glob('*.py'),*r.glob('*.m'),r/'build-provenance.json']}}
(r/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps({'runs':summary['runs'],'speedup':summary['speedup_np4_vs_np1'],'np1_serial_max':[max(x[j] for x in v['np1_vs_serial']['relative_L2_by_frame']) for j in [0,1]],'np4_serial_max':[max(x[j] for x in v['np4_vs_serial']['relative_L2_by_frame']) for j in [0,1]]},indent=2))