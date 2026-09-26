import pathlib,json,re,csv,hashlib,tarfile
r=pathlib.Path('/home/lxy/green-laplace-unidirectional-time-series-runs/hos-deepwater-fourphase-20260926-v1')
assert json.loads((r/'status.json').read_text())['state']=='completed'
assert json.loads((r/'tolerance16-status.json').read_text())['state']=='completed'
assert json.loads((r/'tolerance-status.json').read_text())['state']=='completed'
assert json.loads((r/'tolerance-high-status.json').read_text())['state']=='completed'
records=[]
for p in sorted(r.rglob('*resources.txt')):
    text=p.read_text();rss=re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)',text);wall=re.search(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*(\S+)',text);user=re.search(r'User time \(seconds\):\s*(\S+)',text);exitcode=re.search(r'Exit status:\s*(\d+)',text)
    if not (rss and wall and user and exitcode):continue
    sec=0
    for n in wall[1].split(':'):sec=60*sec+float(n)
    records.append({'file':str(p.relative_to(r)),'wall_seconds':sec,'user_seconds':float(user[1]),'peak_RSS_KiB':int(rss[1]),'exit_code':int(exitcode[1])})
with (r/'performance.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=records[0].keys());w.writeheader();w.writerows(records)
summary={'state':'completed_with_tolerance_audit','final_tolerance':{'akp002':1e-16,'akp012':1e-14},'metrics':{},'performance':records,'hashes':{str(p.relative_to(r)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [*r.glob('*.m'),*r.glob('*.py'),r/'mf12-no31-source.tar.gz']}}
for name in ['akp002','akp012']:
    summary['metrics'][name]=list(csv.DictReader((r/('tol16' if name=='akp002' else 'tol14')/name/'metrics.csv').open()))
(r/'final-summary.json').write_text(json.dumps(summary,indent=2))
files=[r/'final-summary.json',r/'performance.csv',r/'support_and_tolerance_audit_final.csv',r/'status.json',r/'tolerance-status.json',r/'tolerance-high-status.json']
for name in ['akp002','akp012']:
    files.extend([r/name/'initialization.json',r/('tol16' if name=='akp002' else 'tol14')/name/'metrics.csv',r/('tol16' if name=='akp002' else 'tol14')/name/'analysis.json',r/('tol16' if name=='akp002' else 'tol14')/name/'comparison.png',r/('tol16' if name=='akp002' else 'tol14')/name/'comparison_zoom.png'])
with tarfile.open(r/'compact-results.tar.gz','w:gz') as tar:
    for p in files:tar.add(p,arcname=str(p.relative_to(r)))
for name in summary['metrics']:print(name,summary['metrics'][name])
for d in records:
    if '/tol' in d['file'] or d['file'].startswith('tol14/') or d['file'].startswith('prepare-') or d['file'].startswith('analyze-'):print(d)