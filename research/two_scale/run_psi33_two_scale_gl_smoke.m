function result = run_psi33_two_scale_gl_smoke()
%RUN_PSI33_TWO_SCALE_GL_SMOKE Diagonal, support, zero and cost smoke.
here=fileparts(mfilename('fullpath')); addpath(here);
nx=32; ny=8; kxv=[0:nx/2-1,-nx/2:-1]; kyv=[0:ny/2-1,-ny/2:-1];
[kx,ky]=meshgrid(kxv,kyv); amplitude=1e-3;
eta11_hat=zeros(ny,nx); eta11_hat(1,2)=amplitude*nx*ny;
modes=struct('kx',1,'ky',0,'amplitude',amplitude);
base=struct('q_min',0.9,'q_max',1.1,'forward_axis',[1,0], ...
    'cone_half_angle',pi/6,'support_tolerance',1e-13);
cases={'two_scale','inner_shared_outer_two_scale','shared_scale'};
outputs=struct(); all_pass=true;
for ii=1:numel(cases)
    cfg=base; cfg.mode=cases{ii};
    started=tic;
    [~,hat,meta]=finite_depth_directional_psi33_two_scale_gl(eta11_hat,kx,ky,cfg);
    runtime=toc(started);
    direct=direct_ordered_psi33_two_scale_gl(modes,cfg);
    expected=direct.records(1).coefficient;
    actual=hat(1,4)/(nx*ny);
    rel=abs(actual-expected)/max(abs(expected),1e-30);
    mask=stage3_mask(kx,ky,base);
    passed=rel<5e-11 && hat(1,1)==0 && norm(hat(~mask),2)==0;
    all_pass=all_pass&&passed;
    outputs.(cases{ii})=struct('actual_real',real(actual),'actual_imag',imag(actual), ...
        'expected_real',real(expected),'expected_imag',imag(expected), ...
        'relative_error',rel,'zero_output',abs(hat(1,1)), ...
        'off_support_norm',norm(hat(~mask),2),'runtime_seconds',runtime, ...
        'counts',meta.counts,'pass',passed);
    fprintf('%s Psi33 diagonal rel = %.6e; FFT/IFFT/products = %d/%d/%d\n', ...
        cases{ii},rel,meta.counts.fft,meta.counts.ifft,meta.counts.products);
end
result=struct('schema_version',1,'matlab_release',version('-release'), ...
    'grid',[ny,nx],'cases',outputs,'pass',all_pass);
artifact_dir=fullfile(here,'..','..','artifacts','finite_depth_directional_psi33_two_scale_gl');
if ~exist(artifact_dir,'dir'), mkdir(artifact_dir); end
artifact=fullfile(artifact_dir,'smoke_result.json');
fid=fopen(artifact,'w'); cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'%s',jsonencode(result,'PrettyPrint',true));
fprintf('artifact = %s\n',artifact);
fprintf('OVERALL_PSI33_TWO_SCALE_GL_SMOKE = %s\n',pass_text(all_pass));
if ~all_pass, error('psi33_two_scale_gl:smoke','Psi33 smoke failed.'); end
end

function mask=stage3_mask(kx,ky,cfg)
q=hypot(kx,ky); axis=cfg.forward_axis./norm(cfg.forward_axis);
projection=kx.*axis(1)+ky.*axis(2);
mask=q>=3*cfg.q_min*cos(cfg.cone_half_angle) & q<=3*cfg.q_max ...
    & projection>0 & q>0;
end

function out=pass_text(value)
if value, out='PASS'; else, out='FAIL'; end
end
