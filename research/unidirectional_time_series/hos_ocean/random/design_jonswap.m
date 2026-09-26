function design_jonswap(root)
out=fullfile(root,'jonswap-design-v1');assert(~isfolder(out));mkdir(out);
g=9.81;kp=.0279;h=1/kp;lambda=2*pi/kp;wp=sqrt(g*kp*tanh(kp*h));fp=wp/(2*pi);Tp=1/fp;
q=linspace(.05,30,300001)';f=fp*q;rows={};shapes=zeros(numel(q),2);
for j=1:2
 gamma=[3.3 1];gamma=gamma(j);sig=.07*ones(size(q));sig(q>1)=.09;
 S=q.^(-5).*exp(-1.25*q.^(-4)).*gamma.^exp(-.5*((q-1)./sig).^2);S=S/trapz(f,S);shapes(:,j)=S;
 for cut=[2 2.5 3]
  ix=q>=.5&q<=cut;ff=f(ix);ss=S(ix);m0=trapz(ff,ss);mu=trapz(ff,ff.*ss)/m0;beta=sqrt(trapz(ff,(ff-mu).^2.*ss)/m0)/mu;
  rows(end+1,:)={gamma,cut,1-m0,beta};
 end
end
T=cell2table(rows,'VariableNames',{'gamma','fmax_over_fp','omitted_variance_fraction','relative_bandwidth'});writetable(T,fullfile(out,'spectrum_options.csv'));
nx=1024;ny=512;Lx=50*lambda;Ly=20*lambda;
[mx,my]=meshgrid([0:nx/2-1,-nx/2:-1],[0:ny/2-1,-ny/2:-1]);kx=mx*2*pi/Lx;ky=my*2*pi/Ly;k=hypot(kx,ky);theta=atan2(ky,kx);om=sqrt(g*k.*tanh(k*h));freq=om/(2*pi);
use=kx>0&freq>=.5*fp&freq<=2.5*fp&abs(theta)<=pi/3;
kk=k(use);ww=om(use);ff=freq(use);th=theta(use);qr=ff/fp;sig=.07*ones(size(qr));sig(qr>1)=.09;
S=qr.^(-5).*exp(-1.25*qr.^(-4)).*3.3.^exp(-.5*((qr-1)./sig).^2);
sigmaTheta=deg2rad(25/sqrt(2));D=exp(-.5*(th/sigmaTheta).^2)/(sqrt(2*pi)*sigmaTheta*erf((pi/3)/(sqrt(2)*sigmaTheta)));
cg=g*(tanh(kk*h)+kk*h.*sech(kk*h).^2)./(2*ww);
variance=S.*D.*cg./(2*pi*kk)*(2*pi/Lx)*(2*pi/Ly);variance=variance/sum(variance);
Hs=.12/kp;C=sqrt(2*(Hs/4)^2*variance);rng(20260925,'twister');C=C.*exp(2i*pi*rand(size(C)));
F=complex(zeros(ny,nx));F(use)=nx*ny*C;eta=real(ifft2(F));
discreteHs=4*std(eta(:),1);mu=sum(variance.*ff);beta=sqrt(sum(variance.*(ff-mu).^2))/mu;
N=numel(C);report=struct('status','DESIGN_ONLY_LINEAR_PREVIEW_NO_HOS_NO_MF12','kp',kp,'h',h,'kph',kp*h,'Tp',Tp,'fp',fp,'lambda',lambda, ...
 'kpHs',[.02 .12],'Hs',[.02 .12]/kp,'gamma_recommended',3.3,'frequency_cutoffs_over_fp',[.5 2.5], ...
 'energy_direction_sigma_deg',25/sqrt(2),'theta_support_deg',[-60 60],'domain_lambda',[50 20],'domain_m',[Lx Ly],'grid',[nx ny], ...
 'parents',N,'ordered_cross_pairs',N*(N-1),'one_real_array_per_ordered_pair_GiB',8*N*(N-1)/2^30, ...
 'discrete_Hs_high',discreteHs,'discrete_relative_bandwidth',beta,'max_k_over_kp',max(kk)/kp,'nyquist_over_kp',[nx/100 ny/40], ...
 'second_harmonic_axis_bounds_over_kp',[2*max(abs(kx(use))) 2*max(abs(ky(use)))]/kp, ...
 'third_harmonic_axis_bounds_over_kp',[3*max(abs(kx(use))) 3*max(abs(ky(use)))]/kp, ...
 'duration_Tp',80,'duration_s',80*Tp,'scoring_Tp',[10 70],'sample_dt_s',.2,'seed',20260925, ...
 'HOS_RAM_geometric_estimate_GiB',4.25*2,'wall_geometric_scale_from_random_baseline_hours',(4934/3600)*2*(80*Tp/220), ...
 'resource_note','Grid/time scaling is not a benchmark or guarantee; broad high frequencies and larger Hs can reduce adaptive dt. MF12 pair memory is a separate preparation constraint.');
assert(abs(discreteHs-Hs)/Hs<1e-12);
save(fullfile(out,'linear_design.mat'),'C','kx','ky','use','report','-v7.3');fout=fopen(fullfile(out,'report.json'),'w');fprintf(fout,'%s',jsonencode(report,PrettyPrint=true));fclose(fout);
figh=figure('Visible','off','Position',[40 40 1400 950]);tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
nexttile;hold on;for j=1:2,plot(q,fp*shapes(:,j),'LineWidth',1.3);end;xline(2.5,'k:');xlim([.4 3.2]);grid on;xlabel('f / f_p');ylabel('Normalized variance density');legend('JONSWAP gamma=3.3','gamma=1 (PM limit)','Proposed cutoff','Location','northeast');title('Spectrum design before retained-band Hs normalization');
nexttile;imagesc((0:nx-1)*50/nx,(0:ny-1)*20/ny,eta);set(gca,'YDir','normal');axis image;colorbar;xlabel('x / lambda_p');ylabel('y / lambda_p');title('Linear design preview only, kp Hs=0.12; not a HOS result');
nexttile;tt=(0:.2:80*Tp)';A=C.*exp(1i*(kx(use)*Lx/2+ky(use)*Ly/2));z=zeros(size(tt));for start=1:256:N,ii=start:min(start+255,N);z=z+real(exp(-1i*tt*ww(ii).')*A(ii));end;plot(tt/Tp,z,'k');grid on;xlim([0 80]);xlabel('t / T_p');ylabel('Linear elevation (m)');title('Linear center-probe preview over 80 Tp, same proposed random realization');
exportgraphics(figh,fullfile(out,'jonswap_design.png'),'Resolution',160);close(figh);disp(T);disp(report);
end
