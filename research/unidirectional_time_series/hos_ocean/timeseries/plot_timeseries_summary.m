function plot_timeseries_summary(root)
f=figure('Visible','off','Position',[50 50 1450 820]);tl=tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
for akp=[.02 .12]
    name=sprintf('akp%03d',round(100*akp));s=load(fullfile(root,name,'comparison.mat'));
    for col=1:2
        if col==1,d=s.H;label='HOS';else,d=s.O;label='OW3D';end
        tx=(d.t-d.t(1))/s.r.Tp-8;[~,i]=max(abs(hilbert(d.eta1)));limits=tx(i)+[-2 2];
        nexttile;hold on;plot(tx,d.reference,'k-','LineWidth',1.8);plot(tx,d.pred(:,1),'--','Color',[.1 .4 .8],'LineWidth',1.25);plot(tx,d.pred(:,3),'-.','Color',[.1 .65 .3],'LineWidth',1.25);plot(tx,d.pred(:,4),':','Color',[.85 .15 .1],'LineWidth',1.25);xlim(limits);grid on;
        xlabel('(t - nominal focus time) / T_p');ylabel('Second harmonic, eta (m)');
        rows=strcmp(s.metrics.solver,label)&strcmp(s.metrics.filter,'raw')&strcmp(s.metrics.window,'main_group');
        err=s.metrics.relative_L2(rows&strcmp(s.metrics.method,'GL6'));title(sprintf('%s, Akp=%.2f | GL6 main-group L2 = %.3f%%',label,akp,100*err));
        if akp==.02&&col==1,legend('Reference','GL6','GL12','MF12','Location','northwest');end
    end
end
title(tl,'Single-point time-series reconstruction: identical processing, near-matched wavegroups');exportgraphics(f,fullfile(root,'timeseries_summary.png'),'Resolution',170);exportgraphics(f,fullfile(root,'timeseries_summary.pdf'),'ContentType','vector');close(f);
end