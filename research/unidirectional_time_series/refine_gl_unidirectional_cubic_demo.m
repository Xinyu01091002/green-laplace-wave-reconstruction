function refine_gl_unidirectional_cubic_demo
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));addpath(fileparts(mfilename('fullpath')));
out=fullfile(root,'artifacts','unidirectional_time_series','unidirectional-gl-cubic-demo');
assert(~isfile(fullfile(out,'high_refinement.mat')));old=load(fullfile(out,'demo.mat'));
d=load(fullfile(root,'results','unidirectional_time_series','ow3d_boundary_kh1_alpha1_akp012','pilot.mat'),'t','raw','report');
N=numel(d.t);mask=zeros(N,1);mask(1)=1;mask(2:(N+1)/2)=2;ht=imag(ifft(fft(d.raw).*mask));
eta1=(d.raw(:,1)-d.raw(:,3)-ht(:,2)+ht(:,4))/4;
g=9.81;h=d.report.depth_m;kp=d.report.kp_rad_m;
opts=struct('order',3,'omega_max',4*sqrt(g*kp*tanh(kp*h)),'peak_wavenumber',kp, ...
    'quadrature_rank',8,'domain_lengths',[800 1600],'memory_budget_MiB',4096);
profile clear;profile on;
[result,audit]=gl_unidirectional_time_series(eta1,d.t,g,h,opts);
profile off;p=profile('info');names=string({p.FunctionTable.FunctionName});
assert(~any(contains(names,["time_pairs","time_triples","ordered_pair","pair_reference","integral_reference"])));
z=old.traces{2};score=z.main;
report=struct('audit',audit,'eta22_main_relative',norm(result.eta22(score)-z.reference2(score))/norm(z.reference2(score)), ...
    'eta33_main_relative',norm(result.eta33(score)-z.reference3(score))/norm(z.reference3(score)), ...
    'eta33_full_relative',norm(result.eta33-z.reference3)/norm(z.reference3));
save(fullfile(out,'high_refinement.mat'),'result','report');fid=fopen(fullfile(out,'high_refinement.json'),'w');fprintf(fid,'%s\n',jsonencode(report));fclose(fid);
disp(report);
end
