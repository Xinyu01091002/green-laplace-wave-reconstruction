import pathlib,os,json,shutil,subprocess,re,time
r=pathlib.Path('/home/lxy/green-laplace-unidirectional-time-series-runs/hos-directional-fourphase-20260926-v2');mpi=pathlib.Path('/home/lxy/green-laplace-unidirectional-time-series-runs/hos-mpi-check-20260926-v1');case=r/'probe-check';assert not case.exists();(case/'Results').mkdir(parents=True)
source=r/'cases/phi090';shutil.copy2(source/'prob.inp',case/'prob.inp')
for p in (source/'Results').glob('3d_ini_*.dat'):shutil.copy2(p,case/'Results'/p.name)
s=(source/'input.yml').read_text().replace('duration: 220.0','duration: 0.2').replace('physical space: false','physical space: true');(case/'input.yml').write_text(s)
env=os.environ.copy();env.update(json.loads((r/'environment.json').read_text()));start=time.time()
with (case/'run.log').open('w') as f:p=subprocess.run([str(mpi/'deps/usr/bin/mpirun.openmpi'),'--mca','btl','self,vader,tcp','--bind-to','core','--map-by','core','-np','4',str(r/'bin/HOS-Ocean'),'input.yml'],cwd=case,env=env,stdout=f,stderr=subprocess.STDOUT)
assert p.returncode==0,(case/'run.log').read_text()[-2500:]
lines=(case/'Results/probes.dat').read_text().splitlines();i=next(i for i,x in enumerate(lines) if x.startswith('VARIABLES'));probes=[[float(x) for x in line.split()] for line in lines[i+1:] if line.strip()];assert len(probes)==2
settings=json.loads((r/'settings.json').read_text());checks=[]
for rank in range(4):
    lines=(case/'Results'/('3d_%03d.dat'%rank)).read_text().splitlines();zones=[i for i,x in enumerate(lines) if x.startswith('ZONE')];assert len(zones)==2
    for j,(ix,iy) in enumerate(settings['probe_indices']):
        if (iy-1)//64!=rank:continue
        index=((iy-1)%64)*1024+ix-1
        for it,z in enumerate(zones):
            values=[float(x) for x in lines[z+1+index].split()];value=values[2] if it==0 else values[0];error=abs(value-probes[it][j+1]);checks.append({'time':probes[it][0],'probe':j+1,'native_field':value,'probe_value':probes[it][j+1],'abs_error':error});assert error<1e-10
assert len(checks)==10
initial_error=max(abs(x-y) for x,y in zip(probes[0][1:],settings['expected_initial_probe_eta'][1]));assert initial_error<1e-10
report={'status':'passed','checks':checks,'initial_max_abs_error_m':initial_error,'wall_seconds':time.time()-start,'MPI_ranks':4,'only_IO_broadcast_and_flush_changed':True};(r/'probe-check.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))