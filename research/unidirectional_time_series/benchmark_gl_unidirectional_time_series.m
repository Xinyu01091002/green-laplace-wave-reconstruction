function benchmark_gl_unidirectional_time_series
% Fresh-process timing of the complete non-enumerating API. No references.
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));addpath(fileparts(mfilename('fullpath')));
out=fullfile(root,'artifacts','unidirectional_time_series','unidirectional-gl-demo-refined');
assert(isfile(fullfile(out,'demo.mat')) && ~isfile(fullfile(out,'timing.json')));
source=fullfile(root,'results','unidirectional_time_series','ow3d_boundary_kh1_alpha1_akp002','pilot.mat');
d=load(source,'t','raw','report');N=numel(d.t);mask=zeros(N,1);mask(1)=1;mask(2:(N+1)/2)=2;
ht=imag(ifft(fft(d.raw).*mask));eta1=(d.raw(:,1)-d.raw(:,3)-ht(:,2)+ht(:,4))/4;
g=9.81;h=d.report.depth_m;kp=d.report.kp_rad_m;wp=sqrt(g*kp*tanh(kp*h));
opts=struct('omega_max',4*wp,'peak_wavenumber',kp,'quadrature_rank',8);
timer=tic;[first,audit]=gl_unidirectional_time_series(eta1,d.t,g,h,opts);firstSeconds=toc(timer);
warm=zeros(3,1);
for repeat=1:3
    timer=tic;value=gl_unidirectional_time_series(eta1,d.t,g,h,opts);warm(repeat)=toc(timer);
    assert(norm(value.eta22-first.eta22)/norm(first.eta22)<1e-13);
end
% A declared support change tests complete re-preparation, not a cached plan.
changed=opts;changed.omega_max=3*wp;
timer=tic;[~,changedAudit]=gl_unidirectional_time_series(eta1,d.t,g,h,changed);changedSeconds=toc(timer);
report=struct('first_call_seconds',firstSeconds,'warm_seconds',warm,'median_warm_seconds',median(warm), ...
    'changed_support_seconds',changedSeconds,'original_parent_bins',numel(audit.parent_bins), ...
    'changed_parent_bins',numel(changedAudit.parent_bins),'all_preparation_and_refinement_included',true, ...
    'MATLAB_startup_excluded',true,'profile_enabled',false,'changed_support_is_timing_only',true);
opts.order=3;opts.memory_budget_MiB=2048;
timer=tic;[cubic,cubicAudit]=gl_unidirectional_time_series(eta1,d.t,g,h,opts);report.first_cubic_call_seconds=toc(timer);
cubicTimes=zeros(3,1);
for repeat=1:3
    timer=tic;value=gl_unidirectional_time_series(eta1,d.t,g,h,opts);cubicTimes(repeat)=toc(timer);
    assert(norm(value.eta33-cubic.eta33)/norm(cubic.eta33)<1e-13);
end
report.cubic_warm_seconds=cubicTimes;report.cubic_median_seconds=median(cubicTimes);
report.cubic_converged=cubicAudit.converged;
fid=fopen(fullfile(out,'timing.json'),'w');fprintf(fid,'%s\n',jsonencode(report));fclose(fid);disp(report);
end
