function plot_gl_unidirectional_delivery
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
out=fullfile(root,'artifacts','unidirectional_time_series','unidirectional-gl-cubic-demo');
d=load(fullfile(out,'demo.mat'));refined=load(fullfile(out,'high_refinement.mat'));
d.traces{2}.result=refined.result;
fig=figure('Visible','off','Color','w','Position',[100,100,1200,730]);tiledlayout(2,2,'TileSpacing','compact');
for j=1:2
    z=d.traces{j};tx=(z.result.t-z.result.t(1))/z.Tp;
    nexttile(j);plot(tx,z.reference2,'k-','DisplayName','OW3D second phase sector');hold on;
    plot(tx,z.result.eta22,'--','DisplayName','GL8 operator');grid on;xlim([min(tx(z.main)),max(tx(z.main))]);
    title(sprintf('Akp %.2f, eta22',d.reports{j}.akp));ylabel('Elevation (m)');legend('Location','best');
    nexttile(j+2);plot(tx,z.reference3,'k-','DisplayName','OW3D third phase sector');hold on;
    plot(tx,z.result.eta33,'--','DisplayName','GL8 operator');grid on;xlim([min(tx(z.main)),max(tx(z.main))]);
    title('eta33');xlabel('Time from record start / Tp');ylabel('Elevation (m)');legend('Location','best');
end
exportgraphics(fig,fullfile(out,'final_demo.png'),'Resolution',150);
exportgraphics(fig,fullfile(out,'final_demo.pdf'),'ContentType','vector');close(fig);
end
