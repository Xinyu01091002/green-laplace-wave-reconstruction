function plot_r4_waveforms(root)
s=load(fullfile(root,'input.mat'));ref=load(fullfile(root,'mf12.mat'));a=load(fullfile(root,'r4_eta_rejected_capture.mat'));b=load(fullfile(root,'r4_psi_rejected_capture.mat'));
E=2*a.eta20;P=2*sqrt(s.g)*b.psi20;x=(0:s.nx-1)*s.Lx/s.nx/(2*pi/s.kp);iy=s.ny/2+1;
out=fullfile(root,'waveform-review-v1');assert(~isfolder(out));mkdir(out);
for kind=1:2
 if kind==1,truth=ref.eta20;pred=E;unit='Elevation (m)';name='eta20';else,truth=ref.psi20;pred=P;unit='Surface potential (m^2/s)';name='psi20';end
 f=figure('Visible','off','Position',[40 40 1450 1050]);tiledlayout(4,1,'TileSpacing','compact','Padding','compact');
 nexttile;plot(x,truth(iy,:),'k-',x,pred(iy,:),'r--','LineWidth',1.05);grid on;xlim([0 50]);ylabel(unit);title([name ': centerline, full domain, same physical units']);legend('MF12','R4','Location','northwest');
 nexttile;plot(x,truth(iy,:),'k-',x,pred(iy,:),'r--','LineWidth',1.3);grid on;xlim([20 30]);ylabel(unit);title('Fixed central interval: 20-30 peak wavelengths');
 [~,peak]=max(abs(truth(iy,:)));xp=x(peak);
 nexttile;plot(x,truth(iy,:),'k-',x,pred(iy,:),'r--','LineWidth',1.3);grid on;xlim([max(0,xp-2) min(50,xp+2)]);ylabel(unit);title('Around largest absolute MF12 centerline excursion (selected from reference only)');
 nexttile;plot(x,pred(iy,:)-truth(iy,:),'Color',[.1 .35 .8],'LineWidth',1);grid on;xlim([0 50]);ylabel(['R4 - MF12; ' unit]);xlabel('x / lambda_p');title('Raw difference on its own vertical scale; no fitting, shifting or filtering');
 exportgraphics(f,fullfile(out,[name '_waveforms.png']),'Resolution',160);close(f);
end
end
