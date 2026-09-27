"""Report every official-reference comparison, including failures."""
import pathlib,subprocess,sys
root=pathlib.Path(__file__).resolve().parent
failed=False
for path in sorted(root.glob('grid*_ocean_case*.bin')):
    r=subprocess.run([str(root/'build/hos_ocean_check'),str(path)],capture_output=True,text=True)
    print(path.name,flush=True);print(r.stdout+r.stderr,flush=True);print('exit='+str(r.returncode),flush=True)
    failed=failed or r.returncode!=0
sys.exit(1 if failed else 0)
