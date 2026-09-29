"""Collect compact outputs from the isolated eta33 steepness supplement."""
import json
import pathlib
import subprocess

repo=pathlib.Path(__file__).resolve().parents[2]
local=repo/'artifacts/unidirectional_time_series/hos-eta33-steepness-20260929-v1'
location=json.loads((local/'location.json').read_text())
remote=location['remote'];key=str(pathlib.Path.home()/'.ssh/id_ed25519_cursor')
ssh=['ssh','-p','60093','-i',key,'-o','BatchMode=yes','lxy@60.188.112.99']
status=json.loads(subprocess.check_output([*ssh,'cat '+remote+'/status.json'],text=True))
assert status['state']=='completed' and status.get('exit_code',0)==0,status
(local/'status.json').write_text(json.dumps(status,indent=2))
out=local/'analysis';out.mkdir(exist_ok=False)
names=['eta33_steepness_metrics.csv','eta33_steepness_report.json','eta33_steepness_ladder.png','eta33_steepness_ladder.pdf',
       'eta33_steepness_experiment_final.png','eta33_steepness_experiment_final.pdf','eta33_steepness_fit.json']
subprocess.run(['scp','-P','60093','-i',key,'-o','BatchMode=yes',*[f'lxy@60.188.112.99:{remote}/analysis/{n}' for n in names],str(out)],check=True)
rank=local/'rank_analysis';rank.mkdir(exist_ok=False)
rank_names=['eta33_rank_metrics.csv','eta33_rank_ladder.png','eta33_rank_ladder.pdf']
subprocess.run(['scp','-P','60093','-i',key,'-o','BatchMode=yes',*[f'lxy@60.188.112.99:{remote}/rank_analysis/{n}' for n in rank_names],str(rank)],check=True)
for amplitude in ['004','006','008','010']:
    source=f'{remote}/akp{amplitude}/initialization.json'
    subprocess.run(['scp','-P','60093','-i',key,'-o','BatchMode=yes',f'lxy@60.188.112.99:{source}',str(local/f'akp{amplitude}-initialization.json')],check=True)
print((out/'eta33_steepness_metrics.csv').read_text())
print((out/'eta33_steepness_report.json').read_text())
