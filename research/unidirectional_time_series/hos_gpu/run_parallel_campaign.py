"""Adopt the running first phase and launch the other three on the same GPU."""
import csv,ctypes,hashlib,json,math,os,pathlib,subprocess,time,traceback
r=pathlib.Path(__file__).resolve().parent
exe=pathlib.Path('/root/hos-long-20260927T040628Z/build/hos_run')
s=json.loads((r/'settings.json').read_text());adopt=json.loads((r/'parallel-adoption.json').read_text())
state={'state':'starting','precision':'fp64','tolerance':1e-12,'duration_s':s['duration_s'],'cases':{},'started_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'solver_sha256':hashlib.sha256(exe.read_bytes()).hexdigest()}
children={};handles=[]
def identity(pid):
    try:return pathlib.Path('/proc',str(pid),'stat').read_text().split()[21]
    except FileNotFoundError:return None
def save():
    state['updated_utc']=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime());p=r/'parallel-status.tmp';p.write_text(json.dumps(state,indent=2));p.replace(r/'parallel-status.json')
nv=ctypes.CDLL('libnvidia-ml.so.1');assert nv.nvmlInit_v2()==0
device=ctypes.c_void_p();assert nv.nvmlDeviceGetHandleByIndex_v2(0,ctypes.byref(device))==0
class Memory(ctypes.Structure):_fields_=[('total',ctypes.c_ulonglong),('free',ctypes.c_ulonglong),('used',ctypes.c_ulonglong)]
def memory():
    v=Memory();assert nv.nvmlDeviceGetMemoryInfo(device,ctypes.byref(v))==0;return v.used
try:
    for name,digest in json.loads((r/'initial-manifest.json').read_text()).items():assert hashlib.sha256((r/name).read_bytes()).hexdigest()==digest
    for phase in [0,90,180,270]:
        name='phi%03d'%phase
        checks=[json.loads(x) for x in (r/(name+'-initial-check.jsonl')).read_text().splitlines()];assert len(checks)==6 and all(x['pass'] for x in checks)
    pid=adopt['solver_pid'];assert identity(pid)==adopt['solver_start_ticks']
    state['cases']['phi000']={'pid':pid,'start_ticks':identity(pid),'cpu':31,'state':'running','adopted':True}
    for phase,cpu in [(90,30),(180,29),(270,28)]:
        name='phi%03d'%phase;assert not (r/name).exists()
        f=(r/(name+'.log')).open('x');handles.append(f)
        p=subprocess.Popen(['taskset','-c',str(cpu),str(exe),'double',str(r/(name+'.bin')),str(r/name),str(s['duration_s']),'1e-12'],stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT)
        children[name]=p;state['cases'][name]={'pid':p.pid,'start_ticks':identity(p.pid),'cpu':cpu,'state':'running','adopted':False}
    state['state']='running';save()
    with (r/'parallel-samples.jsonl').open('x') as samples:
        while True:
            for name,case in state['cases'].items():
                if case['state']!='running':continue
                try:case['progress']=json.loads((r/name/'status.json').read_text())
                except (OSError,json.JSONDecodeError):pass
                alive=children[name].poll() is None if name in children else identity(case['pid'])==case['start_ticks']
                if alive:continue
                try:
                    summary=json.loads((r/name/'summary.json').read_text())
                    with (r/name/'diagnostics.csv').open() as f:rows=list(csv.DictReader(f))
                    assert len(rows)==s['expected_samples'] and abs(float(rows[-1]['time'])-s['duration_s'])<1e-7
                    assert all(math.isfinite(float(v)) for row in rows for v in row.values())
                    if name in children:assert children[name].returncode==0
                    dest=r/name/'Results';dest.mkdir()
                    with (dest/'probes.dat').open('w') as f:
                        f.write('VARIABLES = t eta1 eta2 eta3 eta4 eta5\n')
                        for row in rows:f.write(' '.join(row[k] for k in ['time','p1','p2','p3','p4','p5'])+'\n')
                    case.update(state='completed',summary=summary,samples=len(rows),probes_sha256=hashlib.sha256((dest/'probes.dat').read_bytes()).hexdigest())
                except Exception:case.update(state='failed',error=traceback.format_exc())
            state['device_used_bytes']=memory();state['peak_sampled_device_bytes']=max(state.get('peak_sampled_device_bytes',0),state['device_used_bytes'])
            samples.write(json.dumps({'monotonic_seconds':time.monotonic(),'device_used_bytes':state['device_used_bytes'],'phases':{n:c.get('progress',{}) for n,c in state['cases'].items()}})+'\n');samples.flush();save()
            if all(c['state']!='running' for c in state['cases'].values()):break
            time.sleep(10)
    state['state']='completed' if all(c['state']=='completed' for c in state['cases'].values()) else 'completed_with_failures';save()
except Exception:
    state.update(state='controller_failed',error=traceback.format_exc());save();raise
finally:
    nv.nvmlShutdown()
