function test_gl_unidirectional_time_series(outputDirectory)
% Only original FFT-GL spatial operators and input/convergence checks.
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));addpath(root);setup_green_laplace;
addpath(fileparts(mfilename('fullpath')));
if nargin<1,out=fullfile(root,'artifacts','unidirectional_time_series','unidirectional-gl-operator-v3');
else,out=char(outputDirectory);end
assert(~isfolder(out));mkdir(out);profile clear;profile on;
g=9.81;h=1;q0=atanh(sqrt(7)/3);k=q0*[1,2];omega=sqrt(g*k.*tanh(k));
assert(abs(omega(2)/omega(1)-1.5)<1e-13);
A=[.009+.002i,.006-.001i];N=129;period=4*pi/omega(1);t=(0:N-1)'*period/N;
eta1=real(exp(-1i*t*omega)*A.');
opts=struct('omega_max',omega(2)*(1+1e-12),'peak_wavenumber',k(1), ...
    'quadrature_rank',8,'domain_lengths',2*pi/q0,'spatial_points',64);
[prediction,audit]=gl_unidirectional_time_series(eta1,t,g,h,opts);
[qx,qy]=meshgrid(q0*[0:31,-32:-1],[0,1,-2,-1]);
spatial=complex(zeros(N,1));
for it=1:N
    spectrum=complex(zeros(4,64));spectrum(1,[2,3])=256*A/h.*exp(-1i*omega*t(it));
    field=green_laplace_eta22(spectrum,qx,qy,q0,8,string(root));
    spatial(it)=h*field(1,1);
end
nativeError=norm(prediction.eta22_analytic-spatial)/norm(spatial);assert(nativeError<1e-9);
inputError=norm(prediction.eta1_used-eta1)/norm(eta1);assert(inputError<1e-12);
shifted=gl_unidirectional_time_series(eta1,t+7.3,g,h,opts);
assert(norm(shifted.eta22-prediction.eta22)/norm(prediction.eta22)<1e-11);
scaled=gl_unidirectional_time_series(2*eta1,t,g,h,opts);
assert(norm(scaled.eta22-4*prediction.eta22)/norm(prediction.eta22)<1e-11);
% Native-band projection keeps parents and omits only unrepresentable sums.
smallN=33;ts=(0:smallN-1)'*.2;ws=2*pi*[7,10]/(smallN*.2);
es=real(exp(-1i*ts*ws)*[.001;.0007]);
limited=struct('omega_max',max(ws)*(1+1e-12),'domain_lengths',25,'quadrature_rank',8);
[bandResult,bandAudit]=gl_unidirectional_time_series(es,ts,g,10,limited);
assert(max(bandAudit.parent_bins)==10 && max(bandAudit.output_bins)==16 && bandAudit.omitted_output_bins==4);
% The second parent's self/cross sums are all outside the native output band.
onlyFirst=real(.001*exp(-1i*ws(1)*ts));
firstOnlyResult=gl_unidirectional_time_series(onlyFirst,ts,g,10,limited);
bandParity=norm(bandResult.eta22-firstOnlyResult.eta22)/norm(firstOnlyResult.eta22);
assert(bandParity<1e-10);
profile off;p=profile('info');names=string({p.FunctionTable.FunctionName});
forbidden=["time_pairs","time_triples","ordered_pair","pair_reference","integral_reference"];
assert(~any(contains(names,forbidden)),'Forbidden interaction enumerator was called.');
report=struct('native_spatial_GL8_relative',nativeError,'first_input_relative',inputError, ...
    'native_band_projection_passed',true,'out_of_band_parent_effect',bandParity, ...
    'no_interaction_enumerator_called',true,'audit',audit);
save(fullfile(out,'native_checks.mat'),'report','prediction','spatial','eta1','t');
fid=fopen(fullfile(out,'native_checks.json'),'w');fprintf(fid,'%s\n',jsonencode(report));fclose(fid);
fid=fopen(fullfile(out,'called_functions.json'),'w');fprintf(fid,'%s\n',jsonencode(names));fclose(fid);
disp(report);
end
