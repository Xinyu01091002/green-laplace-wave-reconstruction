function plot_ow3d_main_group(out)
% Display-only crop fixed from the shared first-order envelope, not errors.
d=load(fullfile(out,'pilot.mat'));
t=d.t; z=exp(-1i*(t-t(1))*d.omega.')*d.A;
[~,peak]=max(abs(z)); center=t(peak); limits=center+[-2,2]*d.report.Tp_s;
mask=t>=limits(1) & t<=limits(2);
assert(nnz(mask)>5);
set(groot,'defaultFigureVisible','off');
f=figure('Color','w','Position',[100,100,1200,850]);
tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
nexttile; plot(t,d.eta1,'k','LineWidth',1.3); xlim(limits); grid on;
ylabel('\eta_1 (m)'); title('Common first-order input');
colors=[0 .38 .65;.2 .65 .75;.15 .5 .3;.85 .25 .1;.55 .25 .65;.7 .55 .05];
styles={'-','--',':','--','-.',':'};
nexttile; plot(t,d.reference,'k','LineWidth',1.7); hold on;
for j=1:6,plot(t,d.pred(:,j),styles{j},'Color',colors(j,:),'LineWidth',1.2);end
xlim(limits); range=[d.reference(mask),d.pred(mask,:)];
ylim(1.1*[min(range,[],'all'),max(range,[],'all')]); grid on; ylabel('\eta_{22} (m)');
legend(["OW3D second phase sector",d.names],'Location','eastoutside');
nexttile; hold on;
for j=1:6,plot(t,d.pred(:,j)-d.reference,styles{j},'Color',colors(j,:),'LineWidth',1.2);end
xlim(limits); range=d.pred(mask,:)-d.reference(mask);
ylim(1.1*[min(range,[],'all'),max(range,[],'all')]); grid on;
ylabel('Prediction - OW3D (m)'); xlabel('Simulation elapsed time (s)');
akp=0.02; if isfield(d.report,'steepness_akp'),akp=d.report.steepness_akp;end
sgtitle(sprintf('Main group: k_p h = %.3g, Alpha = 1, Ak_p = %.3g | envelope peak +/- 2T_p', ...
    d.report.kp_rad_m*d.report.depth_m,akp));
exportgraphics(f,fullfile(out,'eta22_main_group.png'),'Resolution',180);
exportgraphics(f,fullfile(out,'eta22_main_group.pdf'),'ContentType','vector'); close(f);
fid=fopen(fullfile(out,'main_group_display.json'),'w'); assert(fid>=0);
cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(struct('center_elapsed_s',center,'limits_elapsed_s',limits, ...
    'definition','shared first-order envelope peak +/- 2Tp; display only', ...
    'alignment_or_refit',false)));
diff=d.pred(mask,:)-d.reference(mask);
metric=table(d.names',vecnorm(diff)'/norm(d.reference(mask)), ...
    'VariableNames',{'method','main_group_relative_L2'});
writetable(metric,fullfile(out,'main_group_metrics.csv'));
end
