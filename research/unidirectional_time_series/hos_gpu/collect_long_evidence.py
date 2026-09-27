"""Fetch compact measurements/comparisons only; leave all modal fields remote."""
import argparse,datetime,json,pathlib,shlex,subprocess
p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=pathlib.Path);a=p.parse_args()
key=(pathlib.Path.home()/'.ssh/id_ed25519_cursor').as_posix()
ssh=['ssh','-i',key,'-o','BatchMode=yes','-o',f'ProxyCommand=ssh -p 60093 -i {key} -o BatchMode=yes -W %h:%p root@60.188.112.99','root@192.168.2.92']
code=r'''
import hashlib,json,pathlib,subprocess
b=pathlib.Path('/root/hos-official-benchmark-20260927T041158Z')
g=pathlib.Path('/root/hos-long-20260927T040628Z')
audit=pathlib.Path('/root/hos-global-zero-audit-20260927T043227Z')
result={'remote_roots':{'benchmark':str(b),'gpu_long':str(g),'audit':str(audit)},'measurements':{},'comparisons':{},'hashes':{}}
for name,path in [('cpu8_corrected_2s',b/'corrected-2s/summary.json'),('cpu8_original_2s',b/'unmodified-valid-2s/summary.json'),('gpu_2s',b/'gpu-2s/summary.json'),('gpu_long',g/'matrix.json'),('cpu32_corrected_long',b/'corrected-long32/summary.json'),('audit_provenance',audit/'provenance.json')]:
    if path.exists():result['measurements'][name]=json.loads(path.read_text())
for path in [g/'build/hos_run',audit/'HOS-benchmark-globalzero-fixed',audit/'HOS-benchmark-mpi-sync',audit/'hos_compare_runs',audit/'official_cpu_matrix.py',audit/'short_gpu_benchmark.py']:
    result['hashes'][str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
for precision in ['double','float']:
    path=g/('low-'+precision+'-138s')/'summary.json'
    if path.exists():result['measurements']['gpu_low_138s_'+precision]=json.loads(path.read_text())
if (g/'low-float-138s/summary.json').exists():
    proc=subprocess.run([str(audit/'hos_compare_runs'),'fields',str(g/'low-double-138s'),str(g/'low-float-138s'),'1e-3'],capture_output=True,text=True)
    result['comparisons']['precision_low_138s']={'exit':proc.returncode,'rows':[json.loads(s) for s in proc.stdout.splitlines()],'stderr':proc.stderr}
for family in ['low','high']:
    pairs=[('cpu8_2s_'+family,'probes',b/'corrected-2s'/family/'diagnostics.csv',b/'gpu-2s'/(family+'-double'),'1e-8'),('precision_long_'+family,'fields',g/(family+'-double'),g/(family+'-float'),'1e-3')]
    longref=b/'corrected-long32'/family/'diagnostics.csv'
    if longref.exists():pairs.append(('cpu32_long_'+family,'probes',longref,g/(family+'-double'),'1e-8'))
    for name,kind,ref,cand,tol in pairs:
        proc=subprocess.run([str(audit/'hos_compare_runs'),kind,str(ref),str(cand),tol],capture_output=True,text=True)
        result['comparisons'][name]={'exit':proc.returncode,'rows':[json.loads(s) for s in proc.stdout.splitlines()],'stderr':proc.stderr}
print(json.dumps(result))
'''
r=subprocess.run([*ssh,'python3 -c '+shlex.quote(code)],text=True,capture_output=True,timeout=180);r.check_returncode()
result=json.loads(r.stdout);result['collected_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
result['speed_ratios']=[]
for cpu in result['measurements']['cpu8_corrected_2s']:
    for gpu in result['measurements']['gpu_2s']:
        if cpu['family']==gpu['family']:
            result['speed_ratios'].append({'family':cpu['family'],'gpu_precision':gpu['precision'],'cpu8_over_gpu_advance':cpu['advance_seconds']/gpu['advance_with_diagnostics_seconds'],'same_precision_and_tolerance':gpu['precision']=='double'})
a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result['speed_ratios'],indent=2))
