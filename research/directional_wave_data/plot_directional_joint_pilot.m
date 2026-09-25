function plot_directional_joint_pilot(out)
d=load(fullfile(out,'joint_pilot.mat'));set(groot,'defaultFigureVisible','off');
f=figure('Color','w','Position',[100,100,1350,1050]);tiledlayout(3,1,'Padding','compact','TileSpacing','compact');
nexttile;plot(d.t,d.realInput,'ko','MarkerSize',4);hold on;
plot(d.t,d.input,'-','Color',[0,.4,.7],'LineWidth',1.3);plot(d.t,d.predictedInput(:,1),'--','Color',[.7,.3,.2]);
xlim(d.report.display_limits_s);ylabel('\eta_1 (m)');grid on;
legend('Observed first phase sector','Input projection','Initial linear spectrum only', ...
    'Location','southoutside','Orientation','horizontal');
colors=[0 .4 .7;.1 .6 .4;.7 .3 .2;.55 .25 .65];styles={'-','--',':','-.'};
for variable=1:2
    ref=d.d.eta2;pred=d.eta;label='\eta_{22} (m)';
    if variable==2,ref=d.d.psi2;pred=d.psi;label='\psi_{22} (m^2/s)';end
    nexttile;plot(d.t,ref,'ko','MarkerSize',5,'LineWidth',1.2);hold on;
    for j=1:4,plot(d.t,pred(:,j),styles{j},'Color',colors(j,:),'LineWidth',1.3);end
    xlim(d.report.display_limits_s);grid on;ylabel(label);
    legend('OW3D samples','Joint input: 15 deg','Joint input: 7.5 deg','Initial only','Single direction', ...
        'Location','southoutside','Orientation','horizontal','NumColumns',3);
end
xlabel('Simulation elapsed time (s); lines connect the saved sample times');
sgtitle('Directional joint-input trial: k_p h = 1, spread label 25 deg, Ak_p = 0.02 | 4 s samples');
exportgraphics(f,fullfile(out,'joint_input_comparison.png'),'Resolution',180);
exportgraphics(f,fullfile(out,'joint_input_comparison.pdf'),'ContentType','vector');close(f);
end
