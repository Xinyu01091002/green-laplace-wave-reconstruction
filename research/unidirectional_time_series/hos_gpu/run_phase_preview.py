import pathlib,shlex,subprocess
r='/home/lxy/green-laplace-unidirectional-time-series-runs/gpu-low-postprocess-20260927T085331Z'
cmd=['runuser','-u','lxy','--','taskset','-c','41','/home/lxy/Desktop/matlabr2026a/bin/matlab','-singleCompThread','-batch',"addpath('%s');preview_gpu_phase_records('%s');"%(r,r)]
subprocess.run(['ssh','-p','60093','-i',str(pathlib.Path.home()/'.ssh/id_ed25519_cursor'),'-o','BatchMode=yes','root@60.188.112.99',shlex.join(cmd)],check=True)
