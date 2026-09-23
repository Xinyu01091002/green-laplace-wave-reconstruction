function result = run_two_scale_gl_smoke()
%RUN_TWO_SCALE_GL_SMOKE Small-support diagonal and scope smoke gate.

here = fileparts(mfilename('fullpath'));
addpath(here);
nx = 32;
ny = 8;
kxv = [0:nx/2-1,-nx/2:-1];
kyv = [0:ny/2-1,-ny/2:-1];
[kx,ky] = meshgrid(kxv,kyv);

amplitude = 1e-3;
eta11_hat = zeros(ny,nx);
eta11_hat(1,2) = amplitude*nx*ny; % q=1, positive analytic coefficient A.

base = struct( ...
    'q_min',0.9, ...
    'q_max',1.1, ...
    'forward_axis',[1,0], ...
    'cone_half_angle',pi/6, ...
    'support_tolerance',1e-13);

cfg = base;
cfg.mode = 'two_scale';
tic;
[~,hat_ts,meta_ts] = finite_depth_directional_eta33_two_scale_gl( ...
    eta11_hat,kx,ky,cfg);
time_ts = toc;

cfg.mode = 'shared_scale';
tic;
[~,hat_ss,meta_ss] = finite_depth_directional_eta33_two_scale_gl( ...
    eta11_hat,kx,ky,cfg);
time_ss = toc;

q0 = 1;
d3 = 3*q0^2*(14+15*cosh(2*q0)+6*cosh(4*q0)+cosh(6*q0)) ...
    /(256*sinh(q0)^6);
expected_single_mode = d3*amplitude^3;
actual_ts = hat_ts(1,4)/(nx*ny);
actual_ss = hat_ss(1,4)/(nx*ny);
rel_ts = abs(actual_ts-expected_single_mode)/abs(expected_single_mode);
rel_ss = abs(actual_ss-expected_single_mode)/abs(expected_single_mode);

zero_ts = abs(hat_ts(1,1));
zero_ss = abs(hat_ss(1,1));
off_support_ts = norm(hat_ts(~stage3_mask(kx,ky,base)),2);
off_support_ss = norm(hat_ss(~stage3_mask(kx,ky,base)),2);

result = struct();
result.schema_version = 1;
result.matlab_release = version('-release');
result.grid = [ny,nx];
result.precision = class(eta11_hat);
result.diagonal_expected = expected_single_mode;
result.two_scale = struct('actual',actual_ts,'relative_error',rel_ts, ...
    'zero_output',zero_ts,'off_support_norm',off_support_ts, ...
    'runtime_seconds',time_ts,'counts',meta_ts.counts);
result.shared_scale = struct('actual',actual_ss,'relative_error',rel_ss, ...
    'zero_output',zero_ss,'off_support_norm',off_support_ss, ...
    'runtime_seconds',time_ss,'counts',meta_ss.counts);
result.pass = rel_ts < 5e-11 && rel_ss < 5e-2 ...
    && zero_ts==0 && zero_ss==0 ...
    && off_support_ts==0 && off_support_ss==0;

artifact_dir = fullfile(here,'..','..','artifacts', ...
    'finite_depth_directional_eta33_two_scale_gl');
if ~exist(artifact_dir,'dir')
    mkdir(artifact_dir);
end
artifact = fullfile(artifact_dir,'smoke_result.json');
fid = fopen(artifact,'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid,'%s',jsonencode(result,'PrettyPrint',true));

fprintf('two-scale diagonal relative error = %.6e\n',rel_ts);
fprintf('shared-scale diagonal relative error = %.6e\n',rel_ss);
fprintf('two-scale FFT/IFFT/products = %d/%d/%d\n', ...
    meta_ts.counts.fft,meta_ts.counts.ifft,meta_ts.counts.products);
fprintf('shared-scale FFT/IFFT/products = %d/%d/%d\n', ...
    meta_ss.counts.fft,meta_ss.counts.ifft,meta_ss.counts.products);
fprintf('artifact = %s\n',artifact);
fprintf('OVERALL_ETA33_TWO_SCALE_GL_SMOKE = %s\n',pass_text(result.pass));
if ~result.pass
    error('eta33_two_scale_gl:smoke','Two-Scale GL smoke gate failed.');
end
end

function mask = stage3_mask(kx,ky,cfg)
q = hypot(kx,ky);
axis = cfg.forward_axis./norm(cfg.forward_axis);
projection = kx.*axis(1)+ky.*axis(2);
mask = q>=3*cfg.q_min*cos(cfg.cone_half_angle) ...
    & q<=3*cfg.q_max & projection>0 & q>0;
end

function out = pass_text(value)
if value
    out = 'PASS';
else
    out = 'FAIL';
end
end
