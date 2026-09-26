import pathlib,time,json,subprocess
r=pathlib.Path('/home/lxy/green-laplace-unidirectional-time-series-runs/hos-gl-timeseries-20260926-v1')
while True:
    s=json.loads((r/'status.json').read_text())
    if s['state']=='completed':break
    if s['state']=='failed':raise RuntimeError(s.get('error','failed'))
    time.sleep(5)
cmd="addpath('%s');plot_timeseries_summary('%s');"%(r,r)
with (r/'summary-plot.log').open('w') as f:
    p=subprocess.run(['/home/lxy/Desktop/matlabr2026a/bin/matlab','-singleCompThread','-batch',cmd],cwd=r,stdout=f,stderr=subprocess.STDOUT)
assert p.returncode==0,(r/'summary-plot.log').read_text()[-2000:]
with (r/'collect.log').open('w') as f:p=subprocess.run(['python3',str(r/'collect_timeseries.py')],stdout=f,stderr=subprocess.STDOUT)
assert p.returncode==0,(r/'collect.log').read_text()[-2000:]
print((r/'collect.log').read_text())