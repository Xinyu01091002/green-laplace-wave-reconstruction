function plot_directional_eta20(out)
d=load(fullfile(out,'eta20_fields.mat'));set(groot,'defaultFigureVisible','off');
f=figure('Color','w','Position',[100,100,1200,850]);tiledlayout(3,1,'Padding','compact','TileSpacing','compact');
nexttile;plot(d.t,d.d.eta0,'ko','MarkerSize',4);hold on;plot(d.t,d.eta20raw(:,3),'-','Color',[0,.4,.7]);
grid on;xlim(d.report.display_limits_s);ylabel('\eta (m)');title('Raw phase-zero record and nonzero spatial-difference prediction');
legend('OW3D phase zero','GL16 including stationary nonzero-K terms','Location','southoutside','Orientation','horizontal');
for j=1:2
    nexttile;plot(d.t,d.reference20{j},'ko','MarkerSize',4);hold on;
    plot(d.t,d.varying{j},'LineWidth',1.2);grid on;xlim(d.report.display_limits_s);ylabel('\eta_{20} (m)');
    title(sprintf('Oscillatory part: 0 < |omega| < %.1f omega_p (limited by sampled Nyquist)',d.ratios(j)));
    legend('OW3D','GL6','GL12','GL16','Location','southoutside','Orientation','horizontal');
end
xlabel('Elapsed time (s); original samples only');
sgtitle(sprintf('Directional eta20: kh=%.3g, spread=%g deg, Akp=%.3g', ...
    d.report.metadata.kph,d.report.metadata.spread_label_degrees,d.report.metadata.Akp));
exportgraphics(f,fullfile(out,'eta20_comparison.png'),'Resolution',160);exportgraphics(f,fullfile(out,'eta20_comparison.pdf'),'ContentType','vector');close(f);
end
