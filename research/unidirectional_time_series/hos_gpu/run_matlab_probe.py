"""Copy and run one named diagnostic/plot function in an existing owned run."""
import argparse,pathlib,shlex,subprocess
p=argparse.ArgumentParser();p.add_argument('source',type=pathlib.Path);p.add_argument('remote');p.add_argument('--cpu',default='41');a=p.parse_args();assert a.source.suffix=='.m' and a.source.stem.isidentifier()
opts=['-i',str(pathlib.Path.home()/'.ssh/id_ed25519_cursor'),'-o','BatchMode=yes'];ssh=['ssh','-p','60093',*opts,'root@60.188.112.99']
subprocess.run(['scp','-P','60093',*opts,str(a.source),'root@60.188.112.99:'+a.remote+'/'],check=True)
cmd=['runuser','-u','lxy','--','taskset','-c',a.cpu,'/home/lxy/Desktop/matlabr2026a/bin/matlab','-singleCompThread','-batch',"addpath('%s');%s('%s');"%(a.remote,a.source.stem,a.remote)]
subprocess.run([*ssh,'cd '+shlex.quote(a.remote)+' && '+shlex.join(cmd)],check=True)
