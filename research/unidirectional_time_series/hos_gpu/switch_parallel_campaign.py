"""One-time handoff of this queue only; keep its existing solver alive."""
import hashlib,json,os,pathlib,signal,subprocess,time
r=pathlib.Path(__file__).resolve().parent
old=json.loads((r/'status.json').read_text());controller=json.loads((r/'launcher.json').read_text())['pid'];pid=old['pid']
assert old['state']=='running' and old['active_phase']=='phi000'
assert str(r/'run_low_campaign.py').encode() in pathlib.Path('/proc',str(controller),'cmdline').read_bytes().split(b'\0')
args=pathlib.Path('/proc',str(pid),'cmdline').read_bytes().split(b'\0');assert str(r/'phi000').encode() in args and str(r/'phi000.bin').encode() in args
assert all(not (r/('phi%03d'%p)).exists() for p in [90,180,270])
assert not (r/'parallel-adoption.json').exists()
adoption={'old_controller':controller,'solver_pid':pid,'solver_start_ticks':pathlib.Path('/proc',str(pid),'stat').read_text().split()[21],'old_status':old,'solver_progress':json.loads((r/'phi000/status.json').read_text()),'parallel_script_sha256':hashlib.sha256((r/'run_parallel_campaign.py').read_bytes()).hexdigest()}
with (r/'parallel-adoption.json').open('x') as f:json.dump(adoption,f,indent=2)
os.kill(controller,signal.SIGTERM) # Deliberately NOT killpg: the solver is retained.
for _ in range(30):
    if not pathlib.Path('/proc',str(controller)).exists():break
    time.sleep(.1)
else:raise RuntimeError('Old controller did not exit; do not launch duplicate queue')
assert pathlib.Path('/proc',str(pid),'stat').read_text().split()[21]==adoption['solver_start_ticks']
old.update(state='superseded_by_parallel_controller',status_file='parallel-status.json');(r/'status.json').write_text(json.dumps(old,indent=2))
with (r/'parallel-controller.log').open('x') as log:
    p=subprocess.Popen(['python3',str(r/'run_parallel_campaign.py')],cwd=r,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
(r/'parallel-launcher.json').write_text(json.dumps({'pid':p.pid},indent=2));print(json.dumps({'parallel_controller':p.pid,'retained_solver':pid}))
