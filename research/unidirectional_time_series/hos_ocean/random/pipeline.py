"""Sequential high/low phase-only random HOS, with per-family GL processing."""
import csv
import hashlib
import json
import math
import os
import pathlib
import shutil
import subprocess
import time

root = pathlib.Path(__file__).resolve().parent
previous = pathlib.Path('/home/lxy/green-laplace-unidirectional-time-series-runs/hos-directional-akp002-20260926-v1')
env = os.environ.copy()
env.update({k: '1' for k in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS']})
outcomes = {}

def utc():
    return time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())

def status(stage, **extra):
    p = root / 'pipeline-status.json'
    q = p.with_suffix('.tmp')
    q.write_text(json.dumps(dict(stage=stage, updated_utc=utc(), outcomes=outcomes, **extra), indent=2))
    q.replace(p)

def matlab(run, command, stem):
    with (run / (stem + '.log')).open('w') as f:
        p = subprocess.run(['/usr/bin/taskset', '-c', '40', str(previous/'bin/time'), '-v', '-o',
                            str(run/(stem+'-resources.txt')), '/home/lxy/Desktop/matlabr2026a/bin/matlab',
                            '-singleCompThread', '-batch', "addpath('%s');%s" % (root, command)],
                           cwd=run, env=env, stdout=f, stderr=subprocess.STDOUT)
    if p.returncode:
        raise RuntimeError((run/(stem+'.log')).read_text()[-3000:])

def prepare(run, label):
    run.mkdir()
    (run/'bin').mkdir()
    for f in ['HOS-Ocean', 'time']:
        shutil.copy2(previous/'bin'/f, run/'bin'/f)
    for f in ['environment.json', 'build-provenance.json', 'benchmark-validation.json']:
        shutil.copy2(previous/f, run/f)
    shutil.copy2(root/'prepare_directional_hos.m', run/'prepare_directional_hos.m')
    shutil.copy2(root/'run_directional_hos.py', run/'run_directional_hos.py')
    (run/'selection.json').write_text(json.dumps({'chosen_ranks':8, 'reason':'Reuse completed 2 s MPI4/8 equivalence and measured wait comparison; no new benchmark', 'prior_benchmark':str(previous/'benchmark.json')}))
    matlab(run, "prepare_random_run('%s',%.17g,'%s');" % (run,label,previous/'mf12'), 'prepare')
    settings=json.loads((run/'staging4/settings.json').read_text())
    settings['MPI_ranks_per_phase']=8
    (run/'settings.json').write_text(json.dumps(settings,indent=2))
    for phase in [0,90,180,270]:
        src=run/'staging4/cases'/('phi%03d'%phase)
        case=run/'cases'/('phi%03d'%phase)
        (case/'Results').mkdir(parents=True)
        data=[]
        for rank in range(4):
            data += (src/'Results'/('3d_ini_%03d.dat'%rank)).read_text().splitlines()[68:]
        assert len(data)==1024*256
        for rank in range(8):
            header=['# Phase-only random MF12 order2; eight MPI ranks']*67+['%-20s%12.5E%4s%5d%4s%5d'%('ZONE SOLUTIONTIME = ',0,', I=',1024,', J=',32)]
            (case/'Results'/('3d_ini_%03d.dat'%rank)).write_text('\n'.join(header+data[rank*32768:(rank+1)*32768])+'\n')
        for f in ['input.yml','prob.inp']:
            shutil.copy2(src/f,case/f)

def compare(run):
    out=run/'random-gl-comparison-v1'
    out.mkdir()
    hashes={}
    for phase in [0,90,180,270]:
        lines=(run/'cases'/('phi%03d'%phase)/'Results/probes.dat').read_text().splitlines()
        i=next(i for i,x in enumerate(lines) if x.startswith('VARIABLES'))
        rows=[[float(v) for v in x.split()] for x in lines[i+1:] if x.strip()]
        assert len(rows)==1101 and all(len(x)==6 and all(math.isfinite(v) for v in x) for x in rows)
        assert max(abs(row[0]-.2*i) for i,row in enumerate(rows))<1e-8
        p=out/('phi%03d.csv'%phase)
        with p.open('w',newline='') as f:
            csv.writer(f).writerows(rows)
        hashes[p.name]=hashlib.sha256(p.read_bytes()).hexdigest()
    (out/'snapshot.json').write_text(json.dumps(dict(sample_count=1101,end_time_s=220,partial=False,source_run=str(run),files_sha256=hashes,frozen_utc=utc()),indent=2))
    (out/'source').mkdir()
    subprocess.run(['tar','-xzf',str(previous/'directional-preview-source.tar.gz'),'-C',str(out/'source')],check=True)
    matlab(out,"compare_random_full('%s');"%run,'run')
    return str(out)

try:
    assert not (root/'akp012').exists() and not (root/'akp002').exists(), 'Do not overwrite earlier runs'
    for name,label in [('akp012',.12),('akp002',.02)]:
        run=root/name
        status('preparing',active_family=name)
        prepare(run,label)
        status('HOS_running',active_family=name,next_family='akp002' if name=='akp012' else None)
        with (run/'controller.log').open('w') as f:
            p=subprocess.run(['python3',str(run/'run_directional_hos.py')],cwd=run,stdout=f,stderr=subprocess.STDOUT)
        assert p.returncode==0 and json.loads((run/'status.json').read_text())['state']=='completed', (run/'controller.log').read_text()[-3000:]
        outcomes[name]={'HOS':'completed'}
        status('GL_processing',active_family=name)
        try:
            outcomes[name]['GL_result']=compare(run)
            outcomes[name]['GL']='completed'
        except Exception as exc:
            outcomes[name]['GL']='failed'
            outcomes[name]['GL_error']=str(exc)
    status('completed_HOS_and_GL' if all(v.get('GL')=='completed' for v in outcomes.values()) else 'HOS_complete_GL_needs_attention')
except Exception as exc:
    status('failed',error=str(exc))
    raise
