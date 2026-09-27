"""Run on 93 in a fresh directory: re-link unchanged actual MPI solver objects."""
import argparse,hashlib,json,os,pathlib,shlex,subprocess
parser=argparse.ArgumentParser();parser.add_argument('--production',action='store_true');parser.add_argument('--cash-karp',action='store_true');parser.add_argument('--adaptive',action='store_true');options=parser.parse_args()
root=pathlib.Path(__file__).resolve().parent
base=pathlib.Path('/home/lxy/green-laplace-unidirectional-time-series-runs/hos-directional-fourphase-20260926-v2')
dependency=base.parent/'hos-mpi-check-20260926-v1'
build=base/'build';src=base/'source';deps=dependency/'deps/usr';lib=deps/'lib/x86_64-linux-gnu'
compiler='/home/lxy/green-laplace-unidirectional-time-series/artifacts/hos-ocean-v2.1.0/precision-v1/gfortran'
env=dict(os.environ);env.update(json.loads((dependency/'environment.json').read_text()))
env.update(OMPI_ALLOW_RUN_AS_ROOT='1',OMPI_ALLOW_RUN_AS_ROOT_CONFIRM='1')
main=src/'sources/HOS/HOS-ocean.f90';original=main.read_text()
hook='! guessing deta_dt at t=0 for energy output at t=0' if options.production else 'CALL build_derivatives()'
assert original.count(hook)==1
insertion='\nCALL export_gpu_reference()\nCALL MPI_FINALIZE(Statinfo)\nSTOP\n'
patched=original.replace(hook,insertion+hook if options.production else hook+insertion)
marker='END PROGRAM HOS_ocean';assert patched.count(marker)==1
patched=patched.replace(marker,(root/'ocean_probe.f90').read_text()+'\n'+marker)
if (root/'ocean_ck_probe.f90').exists():
    patched=patched.replace(marker,(root/'ocean_ck_probe.f90').read_text()+'\n'+marker)
else:
    patched=patched.replace('IF(cash_karp_check)CALL export_ck_reference(test_id)','IF(cash_karp_check)STOP "Missing CK probe source"')
controller_provenance={}
if (root/'ocean_adaptive_probe.f90').exists():
    begin=original.index('    ! Going to next time step')
    end=original.index('\n',original.index('    dt = h_rk ! saving the step size for next start',begin))+1
    controller=original[begin:end];instrumented=controller
    needle='        CALL CPU_TIME(t_i_indiv)';assert instrumented.count(needle)==1
    instrumented=instrumented.replace(needle,needle+'\n        attempts=attempts+1\n        IF(attempts>10000)STOP "Adaptive reference budget exceeded"\n        trace_before=time_cur')
    needle='        IF (time_next - time_cur < h_rk) h_loc = time_next - time_cur';assert instrumented.count(needle)==1
    instrumented=instrumented.replace(needle,needle+'\n        trace_proposal=h_rk;trace_actual=h_loc')
    needle='            h_rk      = MIN(h_rk, dt_out, dt_lin)';assert instrumented.count(needle)==1
    instrumented=instrumented.replace(needle,needle+'''\n            trace_accept=0
            IF(time_cur>trace_before)trace_accept=1
            WRITE(trace_unit,'(I8,1X,5(ES25.16,1X),I1,1X,ES25.16)') &
                attempts,trace_before,trace_proposal,trace_actual,error,time_cur,trace_accept,h_rk''')
    helper=(root/'ocean_adaptive_probe.f90').read_text().replace('!__ORIGINAL_CONTROLLER__',instrumented)
    patched=patched.replace(marker,helper+'\n'+marker)
    (root/'original_controller.f90.txt').write_text(controller)
    controller_provenance={'original_controller_sha256':hashlib.sha256(controller.encode()).hexdigest(),'instrumented_controller_sha256':hashlib.sha256(instrumented.encode()).hexdigest()}
else:
    patched=patched.replace('IF(adaptive_check .AND. test_id/=3)CALL export_adaptive_reference(test_id)','IF(adaptive_check)STOP "Missing adaptive probe source"')
driver=root/'HOS-ocean-probe.f90';driver.write_text(patched)
includes=[build/'sources/mod',build/'_deps/yaml_parser-build/mod',deps/'include',lib/'fortran/gfortran-mod-15/openmpi',lib/'openmpi/lib']
compile_cmd=[compiler,'-cpp','-DMPI','-D_GNU_FORTRAN_COMPILE_RULE_', '-D__GIT_BRANCH__="rhs-probe"','-D__GIT_COMMIT_HASH__="unchanged-modules"','-O2','-ffree-line-length-none','-fallow-argument-mismatch',*[('-I'+str(p)) for p in includes],'-c',str(driver),'-o',str(root/'probe.o')]
def run(cmd,cwd,log):
    with (root/log).open('w') as f:r=subprocess.run(cmd,cwd=cwd,env=env,stdout=f,stderr=subprocess.STDOUT)
    if r.returncode:raise RuntimeError((root/log).read_text()[-6000:])
run(compile_cmd,root,'compile.log')
args=shlex.split((build/'sources/CMakeFiles/HOS-Ocean.dir/link.txt').read_text())
objects=[]
for i,a in enumerate(args):
    if a.endswith('.o'):
        old=(build/'sources'/a).resolve();objects.append({'path':str(old),'sha256':hashlib.sha256(old.read_bytes()).hexdigest()})
        args[i]=str(root/'probe.o') if a.endswith('/HOS/HOS-ocean.f90.o') else str(old)
args[args.index('-o')+1]=str(root/'HOS-probe')
run(args,root,'link.log')
provenance={'upstream_commit':subprocess.check_output(['git','-C',str(src),'rev-parse','HEAD'],text=True).strip(),
 'main_before_sha256':hashlib.sha256(main.read_bytes()).hexdigest(),'patched_main_sha256':hashlib.sha256(driver.read_bytes()).hexdigest(),
 'module_objects':objects,'source_hashes':{str(p.relative_to(src)):hashlib.sha256(p.read_bytes()).hexdigest() for p in src.glob('sources/**/*.f90')},
 'probe_binary_sha256':hashlib.sha256((root/'HOS-probe').read_bytes()).hexdigest()}
provenance.update(controller_provenance)
(root/'provenance.json').write_text(json.dumps(provenance,indent=2))
template=(base.parent/'hos-jonswap-r4gl-80tp-20260926-v2/low/cases/phi000/input.yml').read_text()
if options.production:
    campaign=base.parent/'hos-jonswap-r4gl-80tp-20260926-v2'
    source_records=[]
    for family in ['low','high']:
        origin=campaign/family/'cases/phi000';case=root/('jonswap_'+family+'_phi000');case.mkdir(exist_ok=False)
        (case/'Results').mkdir();(case/'actual_initial.flag').write_text('Use official post-initialization state\n')
        if options.cash_karp:(case/'cash_karp.flag').write_text('Fixed attempted steps, no campaign\n')
        if options.adaptive:(case/'adaptive.flag').write_text('Bounded 0.3 second trajectory\n')
        cfg=(origin/'input.yml').read_text().replace('duration: 1100.8','duration: 0.01').replace('activate: true','activate: false')
        (case/'input.yml').write_text(cfg)
        initial=case/'Results/3d_ini_000.dat';entries=[]
        with initial.open('w') as dest:
            for rank in range(8):
                p=origin/'Results'/('3d_ini_%03d.dat'%rank);before=hashlib.sha256(p.read_bytes()).hexdigest()
                with p.open() as inp:
                    header=[next(inp) for _ in range(68)]
                    if rank==0:
                        dest.writelines(header[:67]);dest.write(f'{"ZONE SOLUTIONTIME = ":20s}{0:12.5E}{", I=":4s}{1024:5d}{", J=":4s}{512:5d}\n')
                    count=0
                    for line in inp:dest.write(line);count+=1
                    assert count==1024*64,(p,count)
                assert hashlib.sha256(p.read_bytes()).hexdigest()==before
                entries.append({'path':str(p),'sha256':before,'records':count})
        source_records.append({'family':family,'phase':0,'source_parts':entries,'assembled_sha256':hashlib.sha256(initial.read_bytes()).hexdigest(),
          'settings':json.loads((campaign/family/'settings.json').read_text())})
        run([str(deps/'bin/mpirun.openmpi'),'--allow-run-as-root','--bind-to','none','-np','1','taskset','-c','40',str(root/'HOS-probe'),'input.yml'],case,case.name+'-run.log')
    (root/'actual_input_provenance.json').write_text(json.dumps(source_records,indent=2))
    print('Production-size initial-state probe complete:',root,flush=True)
    raise SystemExit(0)
for nx,ny in [(64,32),(128,64)]:
    case=root/f'grid{nx}x{ny}';case.mkdir(exist_ok=False)
    if options.cash_karp:(case/'cash_karp.flag').write_text('Fixed attempted steps, no campaign\n')
    if options.adaptive:(case/'adaptive.flag').write_text('Bounded synthetic trajectory\n')
    cfg=template.replace('    x: 1024','    x: '+str(nx)).replace('    y: 512','    y: '+str(ny))
    cfg=cfg.replace('11260.188722544061','9.0').replace('4504.075489017624','7.0').replace('35.842293906810035','1.3')
    cfg=cfg.replace('duration: 1100.8','duration: 0.01').replace('activate: true','activate: false')
    (case/'input.yml').write_text(cfg)
    run([str(deps/'bin/mpirun.openmpi'),'--allow-run-as-root','-np','1',str(root/'HOS-probe'),'input.yml'],case,case.name+'-run.log')
print('Official module probe complete:',root,flush=True)
