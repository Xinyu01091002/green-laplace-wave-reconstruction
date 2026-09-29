function plot_hos_eta33_steepness_experiment(resultRoot)
% Final compact plot from collected six-point steepness and rank metrics.
arguments
    resultRoot (1,:) char
end
main=readtable(fullfile(resultRoot,'analysis','eta33_steepness_metrics.csv'));
rank=readtable(fullfile(resultRoot,'rank_analysis','eta33_rank_metrics.csv'));
x=main.Akp.^2;y=main.eta33_main_relative;X=[ones(size(x)),x];b=X\y;fit=X*b;
R2=1-sum((y-fit).^2)/sum((y-mean(y)).^2);
pHos=polyfit(log(main.Akp),log(main.HOS_main_norm),1);
pGl=polyfit(log(main.Akp),log(main.GL_main_norm),1);

fig=figure('Visible','off','Color','w','Position',[100 100 1450 900]);
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
nexttile;
plot(main.Akp,100*y,'o','LineWidth',2,'MarkerSize',8,'DisplayName','Measured GL8--HOS');hold on;
dense=linspace(min(main.Akp),max(main.Akp),200)';plot(dense,100*(b(1)+b(2)*dense.^2),'-','LineWidth',2,'DisplayName','E_0+c(Ak_p)^2 fit');
grid on;box on;xlabel('Steepness A k_p');ylabel('Main-window relative L_2 error (%)');
title(sprintf('Monotone error: E_0=%.3f%%, R^2=%.6f',100*b(1),R2));legend('Location','northwest');

nexttile;
for J=[4 6 8]
    use=rank.GL_rank==J;semilogy(rank.Akp(use),100*rank.main_relative(use),'-o','LineWidth',1.8,'MarkerSize',7,'DisplayName',sprintf('GL%d',J));hold on;
end
grid on;box on;xlabel('Steepness A k_p');ylabel('Main-window relative L_2 error (%)');title('Fixed-input GL rank ladder');legend('Location','northwest');

nexttile;
loglog(main.Akp,main.HOS_main_norm,'-o','LineWidth',1.8,'MarkerSize',7,'DisplayName',sprintf('HOS third sector, p=%.2f',pHos(1)));hold on;
loglog(main.Akp,main.GL_main_norm,'-s','LineWidth',1.8,'MarkerSize',7,'DisplayName',sprintf('GL eta_{33}, p=%.2f',pGl(1)));
grid on;box on;xlabel('Steepness A k_p');ylabel('Main-window L_2 norm (m)');title('Both amplitudes remain approximately cubic');legend('Location','northwest');

nexttile;
plot(main.Akp,main.main_norm_ratio,'-o','LineWidth',2,'MarkerSize',8);grid on;box on;
xlabel('Steepness A k_p');ylabel('||GL eta_{33}||_2 / ||HOS third sector||_2');title('GL increasingly underpredicts the HOS-sector norm');
ylim([.88 1]);
sgtitle('Single-direction eta_{33} steepness experiment, k_p h=1');
out=fullfile(resultRoot,'analysis');
exportgraphics(fig,fullfile(out,'eta33_steepness_experiment_final.png'),'Resolution',180);
exportgraphics(fig,fullfile(out,'eta33_steepness_experiment_final.pdf'),'ContentType','vector');close(fig);
report=struct('error_model','E=E0+c*(Akp)^2','E0_fraction',b(1),'c',b(2),'R2',R2, ...
    'HOS_third_norm_power',pHos(1),'GL_eta33_norm_power',pGl(1));
fid=fopen(fullfile(out,'eta33_steepness_fit.json'),'w');fprintf(fid,'%s\n',jsonencode(report));fclose(fid);
disp(report);
end
