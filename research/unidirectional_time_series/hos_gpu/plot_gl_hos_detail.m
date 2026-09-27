function plot_gl_hos_detail(run)
% Plot the already computed reconstruction, without rerunning or fitting it.
out=fullfile(run,'gl-time-comparison-v1');z=load(fullfile(out,'comparison.mat'),'t','results','report');
t=z.t;s=z.report.source_settings;window=[10 15]*s.Tp;mask=t>=window(1)&t<=window(2);
d=z.results{1};hos=d.reference;gl=d.pred(:,2);score=t>=d.limits(1)&t<=d.limits(2);
err=norm(gl(score)-hos(score))/norm(hos(score));zoomerr=norm(gl(mask)-hos(mask))/norm(hos(mask));
f=figure('Visible','off','Position',[50 50 1250 760]);tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
nexttile;plot(t(mask),hos(mask),'k-','LineWidth',1.7);hold on;plot(t(mask),gl(mask),'r--','LineWidth',1.5);grid on;xlim(window);ylabel('Second-harmonic elevation (m)');
title(sprintf('Center probe: HOS vs GL time-series reconstruction | raw 10--70 Tp L2 = %.2f%%',100*err));legend('HOS four-phase reference','GL reconstruction, 1.875 deg','Location','northwest','Box','off');
nexttile;plot(t(mask),gl(mask)-hos(mask),'Color',[.1 .35 .75],'LineWidth',1.3);yline(0,'k:');grid on;xlim(window);xlabel('Time (s)');ylabel('GL - HOS (m)');title('Residual on the same unfiltered records; no amplitude or time adjustment');
exportgraphics(f,fullfile(out,'gl_hos_center_detail.png'),'Resolution',160);exportgraphics(f,fullfile(out,'gl_hos_center_detail.pdf'),'ContentType','vector');close(f);
f=figure('Visible','off','Position',[50 50 1350 1150]);tiledlayout(5,1,'TileSpacing','compact','Padding','compact');
rows=zeros(5,4);
for p=1:5
 d=z.results{p};hos=d.reference;gl=d.pred(:,2);err=norm(gl(score)-hos(score))/norm(hos(score));
 rows(p,:)=[p,s.probes_xy_m(p,2)-s.probes_xy_m(1,2),err,max(abs(gl(score)-hos(score)))];
 nexttile;plot(t(mask),hos(mask),'k-','LineWidth',1.3);hold on;plot(t(mask),gl(mask),'r--','LineWidth',1.2);grid on;xlim(window);ylabel('Elevation (m)');title(sprintf('y offset %.2f m | raw 10--70 Tp L2 = %.2f%%',rows(p,2),100*err));
 if p==1,legend('HOS four-phase reference','GL time-series reconstruction','Location','northwest','Box','off');end
end
xlabel('Time (s)');sgtitle('Second harmonic: fixed 10--15 Tp detail, same units and time origin');exportgraphics(f,fullfile(out,'gl_hos_five_probe_detail.png'),'Resolution',150);close(f);
report=struct('plotted_window_s',window,'scoring_window_s',z.results{1}.limits,'center_raw_scoring_L2',rows(1,3),'center_zoom_L2',zoomerr,'probe_metrics',rows,'definition','Second-harmonic phase sector. GL input is the separated first harmonic plus the initial directional prior; this is not a total-elevation reconstruction.');
fid=fopen(fullfile(out,'gl_hos_detail.json'),'w');fprintf(fid,'%s',jsonencode(report,PrettyPrint=true));fclose(fid);
writematrix([t(mask),z.results{1}.reference(mask),z.results{1}.pred(mask,2),z.results{1}.pred(mask,2)-z.results{1}.reference(mask)],fullfile(out,'gl_hos_center_detail.csv'));
end
