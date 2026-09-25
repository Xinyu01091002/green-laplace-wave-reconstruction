function plot_ow3d_eta22_pilot(out)
d=load(fullfile(out,'pilot.mat'));
t=d.t; rawinput=d.eta1phase; input=d.eta1; reference=d.reference;
pred=d.pred; names=d.names; dt=d.report.sample_dt_s;
set(groot,'defaultFigureVisible','off');
colors=[0 .38 .65;.2 .65 .75;.15 .5 .3;.85 .25 .1;.55 .25 .65;.7 .55 .05];
styles={'-','--',':','--','-.',':'};
fig=figure('Color','w','Position',[100,100,1400,1050]);
tiledlayout(3,2,'TileSpacing','compact','Padding','compact');
nexttile([1,2]); plot(t,rawinput,'Color',[.6 .6 .6]); hold on; plot(t,input,'k');
ylabel('\eta_1 (m)'); xlabel('Simulation elapsed time (s)');
legend('OW3D first phase sector','Common first-order input','Location','best'); grid on;
nexttile([1,2]); plot(t,reference,'k','LineWidth',1.7); hold on;
for j=1:6,plot(t,pred(:,j),styles{j},'Color',colors(j,:),'LineWidth',1.2);end
xlim([t(1),t(end)]);
ylabel('\eta_{22} (m)'); xlabel('Simulation elapsed time (s)'); grid on;
legend(["OW3D second phase sector",names],'Location','eastoutside');
nexttile; hold on;
for j=1:6,plot(t,pred(:,j)-reference,styles{j},'Color',colors(j,:),'LineWidth',1);end
ylabel('Prediction - OW3D (m)'); xlabel('Simulation elapsed time (s)'); grid on;
nexttile; N=numel(t); f=(0:floor(N/2))'/(N*dt);
spec=2*abs(fft([reference,pred]))/N;
semilogy(f,spec(1:numel(f),1),'k','LineWidth',1.4); hold on;
for j=1:6,semilogy(f,spec(1:numel(f),j+1),styles{j},'Color',colors(j,:));end
xlabel('Frequency (Hz)'); ylabel('Fourier amplitude (m)'); grid on; xlim([0,.5]); ylim([1e-10,1e-1]);
akp=0.02; if isfield(d.report,'steepness_akp'),akp=d.report.steepness_akp;end
sgtitle(sprintf('OW3D: k_p h = %.3g, Alpha = 1, Ak_p = %.3g | same input, no alignment', ...
    d.report.kp_rad_m*d.report.depth_m,akp));
exportgraphics(fig,fullfile(out,'eta22_comparison.png'),'Resolution',180);
exportgraphics(fig,fullfile(out,'eta22_comparison.pdf'),'ContentType','vector'); close(fig);
end
