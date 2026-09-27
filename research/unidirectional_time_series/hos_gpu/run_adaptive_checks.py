import json,pathlib,subprocess,sys
root=pathlib.Path(__file__).resolve().parent;failed=False
for path in sorted(root.glob('grid*_adaptive_case*.bin')):
    for precision,tol in [('double',None),('float','1e-8')]:
        cmd=[str(root/'build/hos_adaptive_check'),precision,str(path)]+([tol] if tol else [])
        r=subprocess.run(cmd,capture_output=True,text=True,timeout=180)
        name=path.stem+'-'+precision
        (root/(name+'.log')).write_text(r.stdout+r.stderr)
        print(name,flush=True)
        for line in r.stdout.splitlines():
            if line.startswith('{'):
                row=json.loads(line)
                if row.get('type')!='trace':print(line,flush=True)
        if r.stderr:print(r.stderr,flush=True)
        print('exit='+str(r.returncode),flush=True);failed=failed or r.returncode!=0
sys.exit(1 if failed else 0)
