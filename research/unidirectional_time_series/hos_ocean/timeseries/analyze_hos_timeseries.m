function analyze_hos_timeseries(root,akp)
addpath(fullfile(root,'gl-source','src'));addpath(fullfile(root,'gl-source','research','unidirectional_time_series'));addpath(fullfile(root,'mf12'));
name=sprintf('akp%03d',round(100*akp));out=fullfile(root,name);s=load(fullfile(out,'initial.mat'));r=s.r;
raw=[];
for j=1:4
    file=fullfile(out,sprintf('phase%03d',(j-1)*90),'Results','probes.dat');f=fopen(file);assert(f>=0);
    while ~feof(f),line=fgetl(f);if startsWith(strtrim(line),'VARIABLES'),break;end;end
    v=fscanf(f,'%f',[2,Inf]).';fclose(f);assert(size(v,1)>=2001&&all(isfinite(v),'all'));
    t=v(:,1);assert(max(abs(t-(0:numel(t)-1)'*r.output_dt))<1e-8);
    assert(abs(v(1,2)-s.E(r.probe_index,j))<1e-11);
    mask=t>=32*r.Tp-1e-8&t<=48*r.Tp+1e-8;assert(nnz(mask)==641);raw(:,j)=v(mask,2);time=t(mask);
end
old=load(fullfile(root,'ow-reference',sprintf('ow3d_boundary_kh1_alpha1_akp%03d',round(100*akp)),'pilot.mat'));
[H,ht]=process(raw,time,r.h,r.kp,r.g,r.Tp,'HOS');
[O,ot]=process(old.raw,old.t,old.report.depth_m,old.report.kp_rad_m,9.81,old.report.Tp_s,'OW3D');
parity=norm(O.pred-old.pred,'fro')/norm(old.pred,'fro');assert(parity<1e-10,'Legacy processing does not reproduce saved predictions');
metrics=[ht;ot];writetable(metrics,fullfile(out,'metrics.csv'));
report=struct('settings',r,'legacy_prediction_relative_L2',parity,'HOS_method_seconds',H.seconds,'OW3D_method_seconds',O.seconds,'HOS_parent_count',numel(H.A),'OW3D_parent_count',numel(O.A),'HOS_input_projection_relative_L2',H.projection,'OW3D_input_projection_relative_L2',O.projection,'common_filter','Identity for raw metrics; same native FFT sum-support mask on target and every prediction for common_sum_band','window','32--48 Tp after startup, 641 samples; main group = input envelope peak +/-2Tp','initialization_difference','HOS no31 versus legacy OW3D enabled third primary corrections','boundary_difference','periodic HOS versus wall-bounded OW3D; nominal focus at .67L with wide clearance','comparison','near-matched conditions, not solver parity; no fitted gain/offset/time shift');
writejson(fullfile(out,'report.json'),report);save(fullfile(out,'comparison.mat'),'H','O','r','report','metrics','raw','time','-v7.3');
f=figure('Visible','off','Position',[50 50 1350 900]);tiledlayout(3,2,'TileSpacing','compact');
for col=1:2
    if col==1,d=H;label='HOS';else,d=O;label='OW3D';end
    tx=(d.t-d.t(1))/r.Tp-8;[~,idx]=max(abs(hilbert(d.eta1)));limits=tx(idx)+[-2 2];
    nexttile(col);plot(tx,d.eta1,'k');grid on;xlabel('(t - 40Tp)/Tp, nominal');ylabel('first harmonic (m)');title(sprintf('%s, Akp=%.2f',label,akp));
    nexttile(col+2);plot(tx,d.reference,'k-',tx,d.pred(:,1),'b--',tx,d.pred(:,3),'g-.',tx,d.pred(:,4),'r:','LineWidth',1.1);xlim(limits);grid on;ylabel('second harmonic (m)');legend('Reference','GL6','GL12','MF12','Location','best');
    nexttile(col+4);plot(tx,d.pred(:,[1 3 4])-d.reference,'LineWidth',1);xlim(limits);grid on;xlabel('(t - 40Tp)/Tp, nominal');ylabel('prediction - reference (m)');
end
exportgraphics(f,fullfile(out,'comparison.png'),'Resolution',160);exportgraphics(f,fullfile(out,'comparison.pdf'),'ContentType','vector');close(f);
writetable(table(H.t,H.eta1,H.reference,H.pred(:,1),H.pred(:,3),H.pred(:,4),'VariableNames',{'time_s','first_m','HOS_second_m','GL6_m','GL12_m','MF12_m'}),fullfile(out,'hos_timeseries.csv'));disp(metrics);disp(report);
end
function [d,metrics]=process(raw,t,h,kp,g,Tp,label)
N=numel(t);dt=mean(diff(t));tr=t-t(1);hh=imag(hilbert(raw));eta1phase=(raw(:,1)-raw(:,3)-hh(:,2)+hh(:,4))/4;reference=(raw(:,1)-raw(:,2)+raw(:,3)-raw(:,4))/4;
F=fft(eta1phase)/N;bins=(1:floor((N-1)/2))';w=2*pi*bins/(N*dt);k=zeros(size(w));
for i=1:numel(w),k(i)=fzero(@(z)g*z*tanh(z*h)-w(i)^2,[0,max(1,2*w(i)^2/g+2/h)]);end
keep=k*h>=.3 & 2*w<pi/dt;A=2*conj(F(bins(keep)+1));omega=w(keep);k=k(keep);bb=bins(keep);E=exp(-1i*tr*omega.');z=E*A;eta1=real(z);pred=zeros(N,6);seconds=zeros(1,6);
for j=1:3,rank=[6 8 12];timer=tic;pred(:,j)=real(gl_eta22_time_pairs(A,omega,k,h,kp,tr,rank(j)));seconds(j)=toc(timer);end
timer=tic;[pred(:,4),linear]=mf12_eta22_time(A,k,g,h,tr);seconds(4)=toc(timer);assert(norm(linear-eta1)/norm(eta1)<1e-11);
timer=tic;B=(3-tanh(k*h).^2)./(4*tanh(k*h).^3);pred(:,5)=real(z.*(E*(A.*k.*B)));seconds(5)=toc(timer);
timer=tic;bp=(3-tanh(kp*h)^2)/(4*tanh(kp*h)^3);pred(:,6)=kp*bp*real(z.^2);seconds(6)=toc(timer);
names=["GL6","GL8","GL12","MF12","VWA","Walker"];[~,peak]=max(abs(z));main=abs(t-t(peak))<=2*Tp;
signed=[0:floor(N/2),-floor(N/2):-1]';filter=abs(signed)>=2*min(bb)&abs(signed)<=2*max(bb);all=[reference,pred];filtered=real(ifft(fft(all).*filter));rows={};
for mode=1:2
    if mode==1,data=all;kind='raw';else,data=filtered;kind='common_sum_band';end
    for region=1:2
        if region==1,mask=true(N,1);window='full';else,mask=main;window='main_group';end
        ref=data(mask,1);
        for j=1:6,v=data(mask,j+1);rows(end+1,:)={label,kind,window,names(j),norm(v-ref)/norm(ref),max(abs(v-ref))/max(abs(ref)),norm(v)/norm(ref)};end
    end
end
metrics=cell2table(rows,'VariableNames',{'solver','filter','window','method','relative_L2','relative_Linf','norm_ratio'});
d=struct('t',t,'eta1',eta1,'eta1phase',eta1phase,'reference',reference,'pred',pred,'seconds',seconds,'A',A,'omega',omega,'k',k,'projection',norm(eta1-eta1phase)/norm(eta1phase));
end
function writejson(file,r)
f=fopen(file,'w');fprintf(f,'%s',jsonencode(r,PrettyPrint=true));fclose(f);
end