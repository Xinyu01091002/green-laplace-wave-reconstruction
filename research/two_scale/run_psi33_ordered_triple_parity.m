function result = run_psi33_ordered_triple_parity()
%RUN_PSI33_ORDERED_TRIPLE_PARITY FFT Psi33 versus direct ordered kernel.
here=fileparts(mfilename('fullpath')); addpath(here);
nx=32; ny=16; kxv=[0:nx/2-1,-nx/2:-1]; kyv=[0:ny/2-1,-ny/2:-1];
[kx,ky]=meshgrid(kxv,kyv);
modes.kx=[2;2;3]; modes.ky=[0;1;1];
modes.amplitude=[1.1e-3;-0.7e-3+0.2e-3i;0.9e-3-0.1e-3i];
eta11_hat=zeros(ny,nx);
for ii=1:numel(modes.kx)
    ix=find(kxv==modes.kx(ii),1); iy=find(kyv==modes.ky(ii),1);
    eta11_hat(iy,ix)=modes.amplitude(ii)*nx*ny;
end
base=struct('q_min',1.9,'q_max',3.2,'forward_axis',[1,0], ...
    'cone_half_angle',pi/6,'support_tolerance',1e-13);
cases={'two_scale','inner_shared_outer_two_scale','shared_scale'};
case_results=struct(); all_pass=true;
for icase=1:numel(cases)
    cfg=base; cfg.mode=cases{icase};
    [~,field_hat,meta]=finite_depth_directional_psi33_two_scale_gl(eta11_hat,kx,ky,cfg);
    direct=direct_ordered_psi33_two_scale_gl(modes,cfg);
    targets=unique([[direct.records.kx].',[direct.records.ky].'],'rows');
    max_abs=0; max_rel=0;
    for it=1:size(targets,1)
        tx=targets(it,1); ty=targets(it,2);
        selected=[direct.records.kx]==tx & [direct.records.ky]==ty;
        expected=sum([direct.records(selected).coefficient]);
        actual=field_hat(find(kyv==ty,1),find(kxv==tx,1))/(nx*ny);
        discrepancy=abs(actual-expected);
        max_abs=max(max_abs,discrepancy);
        max_rel=max(max_rel,discrepancy/max(abs(expected),1e-30));
    end
    passed=max_abs<2e-18 && max_rel<2e-10; all_pass=all_pass&&passed;
    case_results.(cases{icase})=struct('max_absolute',max_abs, ...
        'max_relative',max_rel,'pass',passed,'counts',meta.counts);
    fprintf('%s Psi33 ordered parity max abs/rel = %.6e / %.6e\n', ...
        cases{icase},max_abs,max_rel);
end
result=struct('schema_version',1,'matlab_release',version('-release'), ...
    'grid',[ny,nx],'cases',case_results,'pass',all_pass);
artifact_dir=fullfile(here,'..','..','artifacts','finite_depth_directional_psi33_two_scale_gl');
if ~exist(artifact_dir,'dir'), mkdir(artifact_dir); end
artifact=fullfile(artifact_dir,'ordered_triple_parity.json');
fid=fopen(artifact,'w'); cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'%s',jsonencode(result,'PrettyPrint',true));
fprintf('artifact = %s\n',artifact);
fprintf('OVERALL_PSI33_TWO_SCALE_GL_ORDERED_PARITY = %s\n',pass_text(all_pass));
if ~all_pass, error('psi33_two_scale_gl:parity','Psi33 parity failed.'); end
end

function out=pass_text(value)
if value, out='PASS'; else, out='FAIL'; end
end
