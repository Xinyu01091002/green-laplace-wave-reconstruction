function prepare_hos_timeseries(root,akp)
addpath(fullfile(root,'mf12')); out=fullfile(root,sprintf('akp%03d',round(100*akp)));
prepare(out,akp);
end
function c=coeff(A,k,g,h)
c=mf12_spectral_coefficients(3,g,h,real(A),imag(A),k,zeros(size(k)),0,0,struct('enable_subharmonic',false,'disable_third_order_correction',true));
% Explicitly exclude every primary-harmonic potential correction.
c.muStar(:)=0;c.omega=sqrt(g*k.*tanh(k*h));c.third_order_subharmonic_mode='skip';
assert(c.superharmonic_only && all(c.muStar==0));
end
function [e22,p22,e33,p33]=higher(c,L,N)
[e,p]=mf12_spectral_surface(c,L,1,N,1,0);
c2=rmfield(c,{'G_3','G_np2m','G_2npm','G_npmpp'});
[e2,p2]=mf12_spectral_surface(c2,L,1,N,1,0);e33=e-e2;p33=p-p2;
c2.a(:)=0;c2.b(:)=0;c2.G_npm(2:2:end)=0;c2.mu_npm(2:2:end)=0;
[e22,p22]=mf12_spectral_surface(c2,L,1,N,1,0);
end
function prepare(out,akp)
assert(~isfolder(out));mkdir(out);g=9.81;kp=.0279;h=1/kp;L=68*2*pi/kp;N=4096;dk=2*pi/L;
ka=(1:N/2-1)*dk;w=kp/sqrt(2*log(10))*ones(size(ka));w(ka<kp)=.004606;
shape=exp(-(ka-kp).^2./(2*w.^2));[v,idx]=sort(shape,'descend');n=find(cumsum(v)>=.99999*sum(v),1);assert(n<600);bins=sort(idx(1:n));k=ka(bins);
Tp=2*pi/sqrt(g*kp*tanh(kp*h));tf=40*Tp;om=sqrt(g*k.*tanh(k*h));
A=shape(bins)/sum(shape(bins))*(akp/kp).*exp(-1i*k*(.67*L)+1i*om*tf);
assert(3*max(bins)<N/2);E=zeros(N,4);P=E;times=zeros(4,1);
for j=1:4
    timer=tic;c=coeff(A*exp(1i*(j-1)*pi/2),k,g,h);[e,p]=mf12_spectral_surface(c,L,1,N,1,0);times(j)=toc(timer);
    assert(all(isfinite(e)) && all(isfinite(p)));E(:,j)=e(:);P(:,j)=p(:);
    if j==1,[expected22,~,expected33,~]=higher(c,L,N);end
    folder=fullfile(out,sprintf('phase%03d',(j-1)*90));mkdir(folder);mkdir(fullfile(folder,'Results'));
    f=fopen(fullfile(folder,'Results','3d_ini.dat'),'w');for z=1:67,fprintf(f,'# independent MF12 11+20+22+33, no31, line %d\n',z);end
    fprintf(f,'%-20s%12.5E%4s%5d%4s%5d\n','ZONE SOLUTIONTIME = ',0,', I=',N,', J=',1);
    fprintf(f,'%.17e %.17e\n',[e(:),p(:)].');fclose(f);
    f=fopen(fullfile(folder,'input.yml'),'w');
    fprintf(f,'restart: disable\ngravity: %.17g\ndomain size:\n  x: %.17g\ncomputation case:\n  3D: false\n  duration: %.17g\n  type: Initial surface quantities\n',g,L,50*Tp);
    fprintf(f,'numerical parameters:\n  time:\n    integration tolerance: 1.e-10\n    relative tolerance: true\n    Dommermuth initialisation:\n      n: 4\n      Ta: 0.0\n  discretization:\n    x: %d\n  surface nonlinearity order: 5\n  dealiasing:\n    x: 5\nbathymetry:\n  depth: %.17g\noutput:\n  directory: Results\n  dimensional: true\n  frequency: %.17g\n  free surface:\n    physical space: false\n  probes:\n    activate: true\n',N,h,40/Tp);fclose(f);
    f=fopen(fullfile(folder,'prob.inp'),'w');fprintf(f,'%.17g\n',round(.67*N)*L/N);fclose(f);
end
[z1,e2,e3]=separate(E);F=fft(z1)/N;
assert(norm(F(bins+1).'-A)/norm(A)<1e-11);assert(norm(e2-expected22(:))/norm(e2)<1e-10);assert(norm(e3-expected33(:))/norm(e3)<1e-10);
r=struct('akp',akp,'g',g,'kp',kp,'h',h,'kph',1,'min_parent_kh',min(k*h),'Tp',Tp,'L',L,'N',N,'alpha',1,'nominal_focus_time_Tp',40,'retained_shape_L1_fraction',sum(shape(bins))/sum(shape),'parents',n,'phase_generation_seconds',times,'initial_eta_max',max(abs(E),[],'all'),'initialization','11+20+22+33; muStar=0; linear frequencies; no mixed-sign order3','M',5,'qx',5,'Ta',0,'tolerance',1e-10,'relative_tolerance',true,'probe_index',round(.67*N)+1,'probe_x',round(.67*N)*L/N,'duration_Tp',50,'output_dt',Tp/40);
save(fullfile(out,'initial.mat'),'A','k','bins','E','P','r','-v7.3');writejson(fullfile(out,'initialization.json'),r);disp(r);
end
function [z1,e2,e3]=separate(E)
N=size(E,1);c1=(E(:,1)-1i*E(:,2)-E(:,3)+1i*E(:,4))/4;c3=(E(:,1)+1i*E(:,2)-E(:,3)-1i*E(:,4))/4;
f=fft(c1);q=zeros(N,1);q(2:N/2)=2*f(2:N/2);z1=ifft(q);
f=fft(c3);q(:)=0;q(2:N/2)=2*f(2:N/2);e3=real(ifft(q));
e2=(E(:,1)-E(:,2)+E(:,3)-E(:,4))/4;
end

function writejson(file,r)
f=fopen(file,'w');fprintf(f,'%s',jsonencode(r,PrettyPrint=true));fclose(f);
end
