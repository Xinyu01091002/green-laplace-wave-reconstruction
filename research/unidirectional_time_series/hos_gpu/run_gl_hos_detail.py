import pathlib,shlex,subprocess
r='/home/lxy/green-laplace-unidirectional-time-series-runs/gpu-low-postprocess-20260927T085331Z'
s=pathlib.Path(__file__).resolve().parent;local=s.parents[2]/'artifacts/hos_gpu/postprocess-20260927T085331Z';key=str(pathlib.Path.home()/'.ssh/id_ed25519_cursor');opts=['-i',key,'-o','BatchMode=yes']
subprocess.run(['scp','-P','60093',*opts,str(s/'plot_gl_hos_detail.m'),'root@60.188.112.99:'+r+'/'],check=True)
cmd=['runuser','-u','lxy','--','taskset','-c','40','/home/lxy/Desktop/matlabr2026a/bin/matlab','-singleCompThread','-batch',"addpath('%s');plot_gl_hos_detail('%s');"%(r,r)]
subprocess.run(['ssh','-p','60093',*opts,'root@60.188.112.99','cd '+shlex.quote(r)+' && '+shlex.join(cmd)],check=True)
files=['gl_hos_center_detail.png','gl_hos_five_probe_detail.png','gl_hos_detail.json','gl_hos_center_detail.csv','gl_hos_center_detail.pdf']
subprocess.run(['scp','-P','60093',*opts,*['root@60.188.112.99:'+r+'/gl-time-comparison-v1/'+n for n in files],str(local)],check=True)
