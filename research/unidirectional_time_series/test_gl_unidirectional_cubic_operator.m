function test_gl_unidirectional_cubic_operator(outputDirectory)
% The ONLY numerical GL comparator is the original spatial FFT graph.
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));addpath(root);setup_green_laplace;
addpath(fileparts(mfilename('fullpath')));
if nargin<1,out=fullfile(root,'artifacts','unidirectional_time_series','unidirectional-gl-cubic-operator');
else,out=char(outputDirectory);end
assert(~isfolder(out));mkdir(out);profile clear;profile on;
g=9.81;h=1;q0=atanh(sqrt(7)/3);k=q0*[1,2];omega=sqrt(g*k.*tanh(k));
A=[.009+.002i,.006-.001i];N=129;t=(0:N-1)'*(4*pi/omega(1))/N;
eta1=real(exp(-1i*t*omega)*A.');
opts=struct('order',3,'omega_max',omega(2)*(1+1e-12),'peak_wavenumber',k(1), ...
    'quadrature_rank',8,'domain_lengths',2*pi/q0,'spatial_points',64);
[prediction,audit]=gl_unidirectional_time_series(eta1,t,g,h,opts);
[qx,qy]=meshgrid(q0*[0:31,-32:-1],[0,1,-2,-1]);
spatial=complex(zeros(N,1));
for it=1:N
    spectrum=complex(zeros(4,64));spectrum(1,[2,3])=256*A/h.*exp(-1i*omega*t(it));
    [field,nativeAudit]=gl_no_stokes_eta33(spectrum,qx,qy,q0,8,string(root));
    assert(nativeAudit.pair_loops==0 && nativeAudit.triple_loops==0);
    spatial(it)=h*field(1,1);
end
nativeError=norm(prediction.eta33_analytic-spatial)/norm(spatial);assert(nativeError<1e-8);
scaled=gl_unidirectional_time_series(2*eta1,t,g,h,opts);
scaleError=norm(scaled.eta33-8*prediction.eta33)/norm(prediction.eta33);assert(scaleError<1e-10);
profile off;p=profile('info');names=string({p.FunctionTable.FunctionName});
assert(~any(contains(names,["time_pairs","time_triples","ordered_pair","pair_reference","integral_reference"])));
report=struct('native_spatial_GL8_cubic_relative',nativeError,'cubic_homogeneity_relative',scaleError, ...
    'no_interaction_enumerator_called',true,'audit',audit);
save(fullfile(out,'checks.mat'),'prediction','spatial','eta1','t','report');
fid=fopen(fullfile(out,'checks.json'),'w');fprintf(fid,'%s\n',jsonencode(report));fclose(fid);
fid=fopen(fullfile(out,'called_functions.json'),'w');fprintf(fid,'%s\n',jsonencode(names));fclose(fid);
disp(report);
end
