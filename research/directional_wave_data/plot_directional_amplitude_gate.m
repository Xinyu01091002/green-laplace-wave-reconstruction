function plot_directional_amplitude_gate(base)
s=load(fullfile(base,'summary.mat'));ids=find(s.ran);set(groot,'defaultFigureVisible','off');
f=figure('Color','w','Position',[50,50,1450,280*numel(ids)]);tl=tiledlayout(numel(ids),2,'Padding','compact','TileSpacing','compact');
for j=ids(:).'
    d=load(fullfile(s.newFolders(j),'joint_pilot.mat'));ratio=max(abs(d.d.eta2))/s.criteria.center_peak_m;
    for v=1:2
        ref=d.d.eta2;pred=d.eta(:,2);yl='eta22 (m)';
        if v==2,ref=d.d.psi2;pred=d.psi(:,2);yl='psi22 (m^2/s)';end
        nexttile;plot(d.t,ref,'ko','MarkerSize',4);hold on;plot(d.t,pred,'-','Color',[0,.4,.7],'LineWidth',1.3);
        xlim(d.report.display_limits_s);grid on;ylabel(yl);xlabel('Elapsed time (s)');
        title(sprintf('%s: x=%.1f, y=%.1f m; eta22 peak ratio=%.3f', ...
            strrep(char(s.newLabels(j)),'_',' '),d.report.metadata.probe,ratio));
    end
end
lg=legend('OW3D saved samples','Joint GL, 7.5 deg');lg.Layout.Tile='south';
title(tl,sprintf('Amplitude-qualified lateral probes | kh=1, spread=25 deg, Akp=%.2f | sampled eta22 peak >= centerline / 3',s.criteria.Akp));
exportgraphics(f,fullfile(base,'qualified_probes.png'),'Resolution',170);exportgraphics(f,fullfile(base,'qualified_probes.pdf'),'ContentType','vector');close(f);
end
