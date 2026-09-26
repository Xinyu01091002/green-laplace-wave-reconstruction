import pathlib,subprocess
r=pathlib.Path(__file__).resolve().parent
with (r/'waveform-review.log').open('w') as f:
 subprocess.run(['/usr/bin/taskset','-c','40','/home/lxy/Desktop/matlabr2026a/bin/matlab','-singleCompThread','-batch',"addpath('%s');plot_r4_waveforms('%s');"%(r,r)],cwd=r,stdout=f,stderr=subprocess.STDOUT,check=True)