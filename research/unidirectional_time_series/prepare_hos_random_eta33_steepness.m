function prepare_hos_random_eta33_steepness(root)
% Fixed-amplitude, random-phase unidirectional families for eta33 steepness.
addpath(fullfile(root,'mf12'));
amplitudes=.02:.02:.18;seeds=[20260925 20260926 20260927];
for seed=seeds
    for akp=amplitudes,prepare_case(root,seed,akp);end
end
end

function prepare_case(root,seed,akp)
g=9.81;kp=.0279;h=1/kp;L=68*2*pi/kp;N=4096;dk=2*pi/L;
ka=(1:N/2-1)*dk;w=kp/sqrt(2*log(10))*ones(size(ka));w(ka<kp)=.004606;
shape=exp(-(ka-kp).^2./(2*w.^2));[v,index]=sort(shape,'descend');
n=find(cumsum(v)>=.99999*sum(v),1);bins=sort(index(1:n));k=ka(bins);
stream=RandStream('mt19937ar','Seed',seed);randomPhase=2*pi*rand(stream,size(k));
A=shape(bins)/sum(shape(bins))*(akp/kp).*exp(1i*randomPhase);
Tp=2*pi/sqrt(g*kp*tanh(kp*h));E=zeros(N,4);P=E;times=zeros(4,1);
out=fullfile(root,sprintf('seed%d',seed),sprintf('akp%03d',round(100*akp)));
assert(~isfolder(out));mkdir(out);
for phase=1:4
    timer=tic;c=coeff(A*exp(1i*(phase-1)*pi/2),k,g,h);
    [eta,psi]=mf12_spectral_surface(c,L,1,N,1,0);times(phase)=toc(timer);
    assert(all(isfinite(eta))&&all(isfinite(psi)));E(:,phase)=eta(:);P(:,phase)=psi(:);
    folder=fullfile(out,sprintf('phase%03d',(phase-1)*90));mkdir(folder);mkdir(fullfile(folder,'Results'));
    fid=fopen(fullfile(folder,'Results','3d_ini.dat'),'w');
    for line=1:67,fprintf(fid,'# MF12 11+20+22+33; unidirectional random phase; line %d\n',line);end
    fprintf(fid,'%-20s%12.5E%4s%5d%4s%5d\n','ZONE SOLUTIONTIME = ',0,', I=',N,', J=',1);
    fprintf(fid,'%.17e %.17e\n',[eta(:),psi(:)].');fclose(fid);
    fid=fopen(fullfile(folder,'input.yml'),'w');
    fprintf(fid,'restart: disable\ngravity: %.17g\ndomain size:\n  x: %.17g\ncomputation case:\n  3D: false\n  duration: %.17g\n  type: Initial surface quantities\n',g,L,50*Tp);
    fprintf(fid,'numerical parameters:\n  time:\n    integration tolerance: 1.e-10\n    relative tolerance: true\n    Dommermuth initialisation:\n      n: 4\n      Ta: 0.0\n  discretization:\n    x: %d\n  surface nonlinearity order: 5\n  dealiasing:\n    x: 5\nbathymetry:\n  depth: %.17g\noutput:\n  directory: Results\n  dimensional: true\n  frequency: %.17g\n  free surface:\n    physical space: false\n  probes:\n    activate: true\n',N,h,40/Tp);fclose(fid);
    fid=fopen(fullfile(folder,'prob.inp'),'w');fprintf(fid,'%.17g\n',round(.67*N)*L/N);fclose(fid);
end
[first,~,third]=separate(E);F=fft(first)/N;
assert(norm(F(bins+1).'-A)/norm(A)<1e-11);
linearHs=4*sqrt(sum(abs(A).^2)/2);actualSteepness=kp*linearHs/2;
report=struct('nominal_Akp',akp,'random_seed',seed,'phase_model','independent uniform positive-mode phases', ...
    'g',g,'kp',kp,'h',h,'kph',1,'Tp',Tp,'L',L,'N',N,'parents',n, ...
    'retained_shape_L1_fraction',sum(shape(bins))/sum(shape),'linear_Hs_4sigma_m',linearHs, ...
    'actual_kp_Hs_over_2',actualSteepness,'amplitude_rule','same modal magnitudes as focused family; phases randomized; no Hs renormalization', ...
    'initial_eta_max',max(abs(E),[],'all'),'initial_third_max',max(abs(third)), ...
    'initialization','11+20+22+33; muStar=0; linear frequencies; no mixed-sign order3', ...
    'M',5,'qx',5,'Ta',0,'tolerance',1e-10,'relative_tolerance',true, ...
    'probe_index',round(.67*N)+1,'probe_x',round(.67*N)*L/N,'duration_Tp',50,'output_dt',Tp/40, ...
    'scoring_window_Tp',[10 40],'phase_generation_seconds',times);
save(fullfile(out,'initial.mat'),'A','k','bins','E','P','report','-v7.3');writejson(fullfile(out,'initialization.json'),report);
disp(report);
end

function c=coeff(A,k,g,h)
c=mf12_spectral_coefficients(3,g,h,real(A),imag(A),k,zeros(size(k)),0,0, ...
    struct('enable_subharmonic',false,'disable_third_order_correction',true));
c.muStar(:)=0;c.omega=sqrt(g*k.*tanh(k*h));c.third_order_subharmonic_mode='skip';
assert(c.superharmonic_only&&all(c.muStar==0));
end

function [first,second,third]=separate(E)
N=size(E,1);c1=(E(:,1)-1i*E(:,2)-E(:,3)+1i*E(:,4))/4;
c3=(E(:,1)+1i*E(:,2)-E(:,3)-1i*E(:,4))/4;
f=fft(c1);q=zeros(N,1);q(2:N/2)=2*f(2:N/2);first=ifft(q);
f=fft(c3);q(:)=0;q(2:N/2)=2*f(2:N/2);third=real(ifft(q));
second=(E(:,1)-E(:,2)+E(:,3)-E(:,4))/4;
end

function writejson(file,value)
fid=fopen(file,'w');fprintf(fid,'%s',jsonencode(value,PrettyPrint=true));fclose(fid);
end
