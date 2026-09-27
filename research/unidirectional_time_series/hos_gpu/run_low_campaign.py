"""Detached, sequential FP64 replay of the frozen four-phase 80Tp low family."""
import csv,hashlib,json,math,pathlib,subprocess,time,traceback
r=pathlib.Path(__file__).resolve().parent
exe=pathlib.Path('/root/hos-long-20260927T040628Z/build/hos_run')
compare=pathlib.Path('/root/hos-global-zero-audit-20260927T043227Z/hos_compare_runs')
s=json.loads((r/'settings.json').read_text());duration=float(s['duration_s']);expected=int(s['expected_samples'])
state={'state':'preflight','precision':'fp64','tolerance':1e-12,'duration_s':duration,'expected_samples':expected,'cases':{},'executable_sha256':hashlib.sha256(exe.read_bytes()).hexdigest(),'started_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())}
def save():
    state['updated_utc']=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime());p=r/'status.tmp';p.write_text(json.dumps(state,indent=2));p.replace(r/'status.json')
def run(name,out,seconds):
    with (r/(out+'.log')).open('w') as f:
        p=subprocess.Popen(['taskset','-c','31',str(exe),'double',str(r/(name+'.bin')),str(r/out),str(seconds),'1e-12'],stdout=f,stderr=subprocess.STDOUT)
        state.update(active_phase=name,pid=p.pid);save()
        while p.poll() is None:
            try:state['progress']=json.loads((r/out/'status.json').read_text())
            except (OSError,json.JSONDecodeError):pass
            save();time.sleep(10)
    if p.returncode:raise RuntimeError(out+' failed; see preserved log')
    return json.loads((r/out/'summary.json').read_text())
try:
    for name,digest in json.loads((r/'initial-manifest.json').read_text()).items():assert hashlib.sha256((r/name).read_bytes()).hexdigest()==digest,name
    # All four initial probe gates must pass before any full record starts.
    for phase in [0,90,180,270]:
        name='phi%03d'%phase;run(name,name+'-preflight',.2)
        q=subprocess.run([str(compare),'probes',str(r/(name+'-initial.csv')),str(r/(name+'-preflight')),'1e-8'],capture_output=True,text=True)
        (r/(name+'-initial-check.jsonl')).write_text(q.stdout+q.stderr);assert q.returncode==0,name
    state['state']='running';save()
    for phase in [0,90,180,270]:
        name='phi%03d'%phase;record=run(name,name,duration)
        with (r/name/'diagnostics.csv').open() as f:rows=list(csv.DictReader(f))
        assert len(rows)==expected and abs(float(rows[-1]['time'])-duration)<1e-7
        assert all(math.isfinite(float(v)) for row in rows for v in row.values())
        dest=r/name/'Results';dest.mkdir()
        with (dest/'probes.dat').open('w') as f:
            f.write('VARIABLES = t eta1 eta2 eta3 eta4 eta5\n')
            for row in rows:f.write(' '.join(row[k] for k in ['time','p1','p2','p3','p4','p5'])+'\n')
        state['cases'][name]={'state':'completed','summary':record,'samples':len(rows),'probes_sha256':hashlib.sha256((dest/'probes.dat').read_bytes()).hexdigest()};save()
    state.update(state='completed',active_phase=None,postprocessing='HOS probe records ready; GL comparison not run');save()
except Exception:
    state.update(state='failed',error=traceback.format_exc());save();raise
