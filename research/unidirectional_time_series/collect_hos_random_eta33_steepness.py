"""Collect compact metrics and plots from the random-phase steepness campaign."""
import json,pathlib,subprocess
repo=pathlib.Path(__file__).resolve().parents[2];local=repo/'artifacts/unidirectional_time_series/hos-eta33-random-steepness-20260929-v1'
location=json.loads((local/'location.json').read_text());remote=location['remote'];key=str(pathlib.Path.home()/'.ssh/id_ed25519_cursor');ssh=['ssh','-p','60093','-i',key,'-o','BatchMode=yes','lxy@60.188.112.99']
status=json.loads(subprocess.check_output([*ssh,'cat '+remote+'/status.json'],text=True));assert status['state']=='completed' and status.get('exit_code',0)==0,status;(local/'status.json').write_text(json.dumps(status,indent=2))
out=local/'analysis';out.mkdir(exist_ok=False);names=['random_eta33_metrics.csv','random_eta33_aggregate.csv','random_eta33_seed_fits.csv','random_eta33_report.json','random_eta33_steepness.png','random_eta33_steepness.pdf']
subprocess.run(['scp','-P','60093','-i',key,'-o','BatchMode=yes',*[f'lxy@60.188.112.99:{remote}/analysis/{n}' for n in names],str(out)],check=True)
for seed in [20260925,20260926,20260927]:
 for akp in range(2,19,2):
  name=f'seed{seed}-akp{akp:03d}-initialization.json';subprocess.run(['scp','-P','60093','-i',key,'-o','BatchMode=yes',f'lxy@60.188.112.99:{remote}/seed{seed}/akp{akp:03d}/initialization.json',str(local/name)],check=True)
print((out/'random_eta33_aggregate.csv').read_text());print((out/'random_eta33_seed_fits.csv').read_text());print((out/'random_eta33_report.json').read_text())
for folder,names in {
 'hilbert_leakage':['hilbert_leakage_metrics.csv','hilbert_leakage_diagnostic.png','hilbert_leakage_diagnostic.pdf'],
 'taper_diagnostic':['taper_metrics.csv','taper_fits.csv','taper_diagnostic.png','taper_diagnostic.pdf'],
 'tapered_gl':['tapered_gl_metrics.csv','tapered_gl_aggregate.csv','tapered_gl_result.png','tapered_gl_result.pdf']}.items():
 target=local/folder;target.mkdir(exist_ok=False)
 subprocess.run(['scp','-P','60093','-i',key,'-o','BatchMode=yes',*[f'lxy@60.188.112.99:{remote}/{folder}/{n}' for n in names],str(target)],check=True)
print((local/'tapered_gl/tapered_gl_aggregate.csv').read_text())
for folder,names in {
 'boundary_plot':['random_eta33_boundary_effect.png','random_eta33_boundary_effect.pdf'],
 'tapered_rank':['tapered_rank_metrics.csv','tapered_rank.png','tapered_rank.pdf']}.items():
 target=local/folder;target.mkdir(exist_ok=False)
 subprocess.run(['scp','-P','60093','-i',key,'-o','BatchMode=yes',*[f'lxy@60.188.112.99:{remote}/{folder}/{n}' for n in names],str(target)],check=True)
