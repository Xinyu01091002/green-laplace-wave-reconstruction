"""Additional tolerance sensitivity after the initial short-run plateau.
Reuses frozen references and executable; no CPU rerun or hardware changes.
"""
import json,pathlib,re,subprocess,sys
source=pathlib.Path(__file__).resolve().parent;repo=source.parents[2];stamp=sys.argv[1]
assert re.fullmatch(r'\d{8}T\d{6}Z',stamp)
root='/root/hos-production-check-'+stamp;local=repo/'artifacts/hos_gpu'/('production-'+stamp)
key=pathlib.Path.home()/'.ssh/id_ed25519_cursor'
proxy=f'ssh -p 60093 -i {key.as_posix()} -o BatchMode=yes -W %h:%p root@60.188.112.99'
ssh=['ssh','-i',str(key),'-o','BatchMode=yes','-o','ProxyCommand='+proxy,'root@192.168.2.92']
rows=[]
for family in ['low','high']:
    for tol in ['1e-11','1e-12']:
        r=subprocess.run([*ssh,'timeout 300 '+root+'/build/hos_adaptive_check float '+root+'/'+family+'-adaptive.bin '+tol],capture_output=True,text=True,timeout=330)
        row={'family':family,'phase':0,'precision':'float','tolerance_override':tol,'exit':r.returncode,
             'fields':[json.loads(line) for line in r.stdout.splitlines() if line.startswith('{')],'stderr':r.stderr}
        rows.append(row);(local/'strict_fp32.json').write_text(json.dumps(rows,indent=2))
        view=dict(row);view['fields']=[f for f in row['fields'] if f.get('type')=='summary'];print(json.dumps(view),flush=True)
if not all(r['exit']==0 for r in rows):raise SystemExit('Some strict-tolerance checks failed; evidence retained')
