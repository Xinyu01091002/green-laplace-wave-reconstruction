"""Freeze the existing numerical routines and adapt only orchestration/scoring."""
import ast
import hashlib
import json
import pathlib
import subprocess
import tarfile

repo=pathlib.Path(__file__).resolve().parents[4]
source=pathlib.Path(__file__).resolve().parent
out=repo/'artifacts/hos_ocean/random-launch-source-v2'
assert not out.exists(), 'Use a new snapshot directory rather than overwriting'
out.mkdir(parents=True)

def change(text,old,new):
    assert text.count(old)==1,old
    return text.replace(old,new)

controller=(repo/'research/directional_wave_data/hos_low_run/run_low_production.py').read_text(encoding='utf-8-sig')
controller=change(controller,"r=pathlib.Path('/home/lxy/green-laplace-unidirectional-time-series-runs/hos-directional-akp002-20260926-v1')","r=pathlib.Path(__file__).resolve().parent")
(out/'run_directional_hos.py').write_text(controller,encoding='utf-8')
comparison=(repo/'research/directional_wave_data/hos_low_run/compare_directional_full.m').read_text(encoding='utf-8-sig')
comparison=change(comparison,'function compare_directional_full(run)','function compare_random_full(run)')
comparison=change(comparison,"'full-gl-comparison-20260926-v1'","'random-gl-comparison-v1'")
comparison=change(comparison,"[~,peak]=max(abs(hilbert(input)));limits=t(peak)+[-2 2]*s.Tp;main=t>=limits(1)&t<=limits(2);assert(limits(1)>=t(1)&&limits(2)<=t(end),'Main group incomplete');","limits=s.scoring_window_s(:).';main=t>=limits(1)&t<=limits(2);assert(limits(1)>=t(1)&&limits(2)<=t(end),'Fixed interior window incomplete');")
comparison=comparison.replace('main_group','fixed_interior').replace('main_window','interior_window')
comparison=change(comparison,'Directional GL time-series reconstruction — full HOS record, common output filtering','Phase-only random waves — fixed interior window, common output filtering')
(out/'compare_random_full.m').write_text(comparison,encoding='utf-8')
for src in [source/'pipeline.py',source/'prepare_random_run.m',
            repo/'research/directional_wave_data/prepare_compact_wavegroup.m',
            repo/'research/directional_wave_data/hos_production/prepare_directional_hos.m']:
    (out/src.name).write_bytes(src.read_bytes())
for p in out.glob('*.py'):
    ast.parse(p.read_text(encoding='utf-8-sig'),filename=str(p))
manifest={'base_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip(),
          'files_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in out.iterdir() if p.is_file()},
          'scope':'User approved sequential phase-only random high/low labels .12/.02, same seed 20260925'}
(out/'source-manifest.json').write_text(json.dumps(manifest,indent=2))
with tarfile.open(out.parent/'random-launch-source-v2.tar.gz','w:gz') as tar:
    for p in out.iterdir():
        tar.add(p,arcname=p.name)
print(json.dumps(manifest,indent=2))
