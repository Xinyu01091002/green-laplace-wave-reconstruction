function plot_directional_multiprobe(base)
s=load(fullfile(base,'summary.mat'));set(groot,'defaultFigureVisible','off');
fig=figure('Color','w','Position',[50,50,1500,1500]);
layout=tiledlayout(5,2,'Padding','compact','TileSpacing','compact');
for ip=1:5
    d=load(fullfile(s.folders(ip),'joint_pilot.mat'));
    for variable=1:2
        ref=d.d.eta2;pred=d.eta;units='\eta_{22} (m)';
        if variable==2,ref=d.d.psi2;pred=d.psi;units='\psi_{22} (m^2/s)';end
        nexttile;plot(d.t,ref,'ko','MarkerSize',4,'LineWidth',1);hold on;
        plot(d.t,pred(:,2),'-','Color',[0,.4,.7],'LineWidth',1.4);
        plot(d.t,pred(:,1),'--','Color',[.1,.6,.4],'LineWidth',1);
        plot(d.t,pred(:,4),':','Color',[.65,.25,.65],'LineWidth',1.1);
        xlim(d.report.display_limits_s);grid on;ylabel(units);
        title(sprintf('%s: x=%.1f m, y=%.1f m',strrep(char(s.labels(ip)),'_',' '),d.report.metadata.probe));
        if ip==5,xlabel('Elapsed time (s)');end
    end
end
lg=legend('OW3D saved samples','Joint 7.5 deg','Joint 15 deg','Single direction','Orientation','horizontal');
lg.Layout.Tile='south';
title(layout,'Five fixed probes | k_p h=1, spreading label 25 deg, Ak_p=0.02 | 4 s samples');
exportgraphics(fig,fullfile(base,'multiprobe_comparison.png'),'Resolution',180);
exportgraphics(fig,fullfile(base,'multiprobe_comparison.pdf'),'ContentType','vector');close(fig);
end
