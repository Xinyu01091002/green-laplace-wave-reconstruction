function hos_deepwater(action,root,akp)
addpath(fullfile(root,'mf12')); out=fullfile(root,sprintf('akp%03d',round(100*akp)));
if strcmp(action,'prepare'), prepare(out,akp); else, analyze(out); end
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
assert(~isfolder(out));mkdir(out);g=9.81;kp=.0279;h=20/kp;L=68*2*pi/kp;N=4096;dk=2*pi/L;
ka=(1:N/2-1)*dk;w=kp/sqrt(2*log(10))*ones(size(ka));w(ka<kp)=.004606;
shape=exp(-(ka-kp).^2./(2*w.^2));[v,idx]=sort(shape,'descend');n=find(cumsum(v)>=.99999*sum(v),1);assert(n<600);bins=sort(idx(1:n));k=ka(bins);
Tp=2*pi/sqrt(g*kp*tanh(kp*h));tf=40*Tp;om=sqrt(g*k.*tanh(k*h));
A=shape(bins)/sum(shape(bins))*(akp/kp).*exp(-1i*k*L/2+1i*om*tf);
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
    fprintf(f,'restart: disable\ngravity: %.17g\ndomain size:\n  x: %.17g\ncomputation case:\n  3D: false\n  duration: %.17g\n  type: Initial surface quantities\n',g,L,20*Tp);
    fprintf(f,'numerical parameters:\n  time:\n    integration tolerance: 1.e-10\n    Dommermuth initialisation:\n      n: 4\n      Ta: 0.0\n  discretization:\n    x: %d\n  surface nonlinearity order: 5\n  dealiasing:\n    x: 5\nbathymetry:\n  depth: %.17g\noutput:\n  directory: Results\n  dimensional: true\n  frequency: %.17g\n  free surface:\n    physical space: true\n',N,h,4/Tp);fclose(f);
end
[z1,e2,e3]=separate(E);F=fft(z1)/N;
assert(norm(F(bins+1).'-A)/norm(A)<1e-11);assert(norm(e2-expected22(:))/norm(e2)<1e-10);assert(norm(e3-expected33(:))/norm(e3)<1e-10);
r=struct('akp',akp,'g',g,'kp',kp,'h',h,'kph',20,'min_parent_kh',min(k*h),'Tp',Tp,'L',L,'N',N,'alpha',1,'nominal_focus_time_Tp',40,'retained_shape_L1_fraction',sum(shape(bins))/sum(shape),'parents',n,'phase_generation_seconds',times,'initial_eta_max',max(abs(E),[],'all'),'initialization','11+20+22+33; muStar=0; linear frequencies; no mixed-sign order3','M',5,'qx',5,'Ta',0,'tolerance',1e-10);
save(fullfile(out,'initial.mat'),'A','k','bins','E','P','r','-v7.3');writejson(fullfile(out,'initialization.json'),r);disp(r);
end
function [z1,e2,e3]=separate(E)
N=size(E,1);c1=(E(:,1)-1i*E(:,2)-E(:,3)+1i*E(:,4))/4;c3=(E(:,1)+1i*E(:,2)-E(:,3)-1i*E(:,4))/4;
f=fft(c1);q=zeros(N,1);q(2:N/2)=2*f(2:N/2);z1=ifft(q);
f=fft(c3);q(:)=0;q(2:N/2)=2*f(2:N/2);e3=real(ifft(q));
e2=(E(:,1)-E(:,2)+E(:,3)-E(:,4))/4;
end
function analyze(out)
s=load(fullfile(out,'initial.mat'));r=s.r;N=r.N;wanted=[0,3,20]*r.Tp;E=zeros(N,4,3);P=E;
for j=1:4
    file=fullfile(out,sprintf('phase%03d',(j-1)*90),'Results','3d.dat');f=fopen(file);assert(f>=0);hits=false(1,3);frames=0;
    while ~feof(f)
        line=fgetl(f);if ~ischar(line),break;end;if ~startsWith(strtrim(line),'ZONE'),continue;end
        tok=regexp(line,'SOLUTIONTIME\s*=\s*([+\-\d.Ee]+)','tokens','once');t=str2double(tok{1});frames=frames+1;
        if frames==1,nc=4;else,nc=2;end;v=fscanf(f,'%f',[nc,N]);assert(isequal(size(v),[nc,N])&&all(isfinite(v),'all'));
        ix=find(abs(t-wanted)<1e-8*r.Tp);if ~isempty(ix),E(:,j,ix)=v(end-1,:).';P(:,j,ix)=v(end,:).';hits(ix)=true;end
    end
    fclose(f);assert(all(hits));
end
assert(norm(E(:,:,1)-s.E,'fro')/norm(s.E,'fro')<1e-12);assert(norm(P(:,:,1)-s.P,'fro')/norm(s.P,'fro')<1e-12);
rows=[];pred2=zeros(N,3);pred3=pred2;obs2=pred2;obs3=pred2;first=pred2;timings=zeros(3,1);x=(0:N-1)'*r.L/N;
for it=1:3
    [z,e2,e3]=separate(E(:,:,it));F=fft(z)/N;A=F(s.bins+1).';proj=zeros(N,1);proj(s.bins+1)=F(s.bins+1);zp=ifft(proj*N);
    omitted=norm(z-zp)/norm(z);timer=tic;c=coeff(A,s.k,r.g,r.h);[p2,~,p3,~]=higher(c,r.L,N);timings(it)=toc(timer);
    first(:,it)=real(z);pred2(:,it)=p2(:);pred3(:,it)=p3(:);obs2(:,it)=e2;obs3(:,it)=e3;
    for harmonic=[2,3]
        if harmonic==2,p=p2(:);o=e2;else,p=p3(:);o=e3;end
        % Main group window selected exclusively from the first-harmonic envelope.
        [~,center]=max(abs(z));dist=mod(x-x(center)+r.L/2,r.L)-r.L/2;mask=abs(dist)<=2*(2*pi/r.kp);
        rows(end+1,:)=[wanted(it)/r.Tp,harmonic,norm(p-o)/norm(o),max(abs(p-o))/max(abs(o)),norm(p(mask)-o(mask))/norm(o(mask)),omitted,max(abs(o)),max(abs(p))];
    end
end
tab=array2table(rows,'VariableNames',{'time_Tp','harmonic','raw_relative_L2','raw_relative_Linf','main_group_relative_L2','first_projection_relative_L2','observed_peak_m','predicted_peak_m'});writetable(tab,fullfile(out,'metrics.csv'));
report=struct('settings',r,'metrics',table2struct(tab),'reconstruction_seconds',timings,'phase_sector_caveat','four-phase + spatial positive-k projection; not perturbation-order isolation','comparison','same fixed parent support; raw units; no fit/shift','status','completed');writejson(fullfile(out,'analysis.json'),report);
save(fullfile(out,'comparison.mat'),'x','first','pred2','pred3','obs2','obs3','report','-v7.3');
f=figure('Visible','off','Position',[80,80,1500,800]);tiledlayout(2,3);
for h=1:2,for it=1:3,nexttile;if h==1,o=obs2(:,it);p=pred2(:,it);else,o=obs3(:,it);p=pred3(:,it);end
plot(x/(2*pi/r.kp),o,'k-',x/(2*pi/r.kp),p,'r--');xlabel('x / lambda_p');ylabel(sprintf('harmonic %d, eta (m)',h+1));title(sprintf('Akp=%.2f, t=%g Tp',r.akp,wanted(it)/r.Tp));grid on;if h==1&&it==1,legend('HOS four-phase','MF12 from current first harmonic');end
end,end;exportgraphics(f,fullfile(out,'comparison.png'),'Resolution',150);exportgraphics(f,fullfile(out,'comparison.pdf'),'ContentType','vector');close(f);disp(tab);
end
function writejson(file,r)
f=fopen(file,'w');fprintf(f,'%s',jsonencode(r,PrettyPrint=true));fclose(f);
end