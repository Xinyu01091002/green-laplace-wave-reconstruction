"""Preserve the failed root-account launch and use the original MATLAB owner."""
import pathlib,shlex,subprocess
root='/home/lxy/green-laplace-unidirectional-time-series-runs/gpu-low-postprocess-20260927T085331Z'
code=r'''
import json,os,pathlib,pwd,subprocess
r=pathlib.Path(ROOT);assert json.loads((r/'status.json').read_text())['state']=='failed'
for name in ['status.json','controller.log','matlab.log']:(r/name).rename(r/(name+'.root-account-failed'))
user=pwd.getpwnam('lxy')
for base,dirs,files in os.walk(r,followlinks=False):
 os.chown(base,user.pw_uid,user.pw_gid)
 for name in dirs+files:os.chown(pathlib.Path(base)/name,user.pw_uid,user.pw_gid,follow_symlinks=False)
with (r/'controller.log').open('x') as f:
 p=subprocess.Popen(['runuser','-u','lxy','--','python3',str(r/'run_gpu_postprocess.py')],cwd=r,stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT,start_new_session=True)
print(json.dumps({'pid':p.pid,'account':'lxy','root':str(r)}))
'''.replace('ROOT',repr(root))
subprocess.run(['ssh','-p','60093','-i',str(pathlib.Path.home()/'.ssh/id_ed25519_cursor'),'-o','BatchMode=yes','root@60.188.112.99','python3 -c '+shlex.quote(code)],check=True)
