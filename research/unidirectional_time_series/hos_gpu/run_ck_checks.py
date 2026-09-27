import pathlib,subprocess,sys
root=pathlib.Path(__file__).resolve().parent;failed=False
for p in sorted(root.glob('grid*_ck_case*.bin')):
    r=subprocess.run([str(root/'build/hos_ocean_ck_check'),str(p)],capture_output=True,text=True)
    print(p.name,flush=True);print(r.stdout+r.stderr,flush=True);print('exit='+str(r.returncode),flush=True)
    failed=failed or r.returncode!=0
sys.exit(1 if failed else 0)
