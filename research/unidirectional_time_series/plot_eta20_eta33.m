function plot_eta20_eta33(out,orders)
if nargin<2,orders=[20,33];end
d=load(fullfile(out,'fields.mat'));set(groot,'defaultFigureVisible','off');
colors=[0 .4 .7;.15 .65 .65;.15 .45 .2;.7 .2 .5;.8 .55 .05];styles={'-','--',':','-.',':'};
for order=orders
    pred=d.(sprintf('pred%d',order));obs=d.(sprintf('obs%d',order));names=d.(sprintf('names%d',order));
    f=figure('Color','w','Position',[100,100,1250,850]);tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
    nexttile;plot(d.t,d.eta1,'k','LineWidth',1.3);ylabel('\eta_1 (m)');grid on;xlim(d.report.display_limits_s);
    nexttile;plot(d.t,obs,'k','LineWidth',1.8);hold on;
    for j=1:5,plot(d.t,pred(:,j),styles{j},'Color',colors(j,:),'LineWidth',1.3);end
    ylabel(sprintf('\\eta_{%d} (m)',order));grid on;xlim(d.report.display_limits_s);
    legend(["OW3D phase sector",names],'Location','eastoutside');
    nexttile;hold on;
    for j=1:5,plot(d.t,pred(:,j)-obs,styles{j},'Color',colors(j,:),'LineWidth',1.2);end
    ylabel('Prediction - OW3D (m)');xlabel('Simulation elapsed time (s)');grid on;xlim(d.report.display_limits_s);
    suffix='positive triple sums; no MF12 third order';
    if order==20
        ratio=.5;
        if isfield(d.report,'eta20_cutoff_ratio'),ratio=d.report.eta20_cutoff_ratio;end
        suffix=sprintf('nonzero frequencies below %g omega_p; GL diagnostics',ratio);
    end
    sgtitle(sprintf('k_p h = 1, Alpha = 1, Ak_p = %.2f | eta%d | %s',d.report.Akp,order,suffix));
    exportgraphics(f,fullfile(out,sprintf('eta%d_main_group.png',order)),'Resolution',180);
    exportgraphics(f,fullfile(out,sprintf('eta%d_main_group.pdf',order)),'ContentType','vector');close(f);
end
end
