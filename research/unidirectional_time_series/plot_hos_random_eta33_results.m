function plot_hos_random_eta33_results(resultRoot)
% Historical untapered leakage plot; not a valid eta33 comparison.
arguments
    resultRoot (1,:) char
end
T=readtable(fullfile(resultRoot,'analysis','random_eta33_metrics.csv'));
A=readtable(fullfile(resultRoot,'analysis','random_eta33_aggregate.csv'));
T.cosine=(T.norm_ratio.^2+1-T.eta33_relative.^2)./(2*T.norm_ratio);
seeds=unique(T.seed);cosineMedian=zeros(height(A),1);cosineMin=cosineMedian;cosineMax=cosineMedian;
for i=1:height(A)
    use=T.nominal_Akp==A.nominal_Akp(i);cosineMedian(i)=median(T.cosine(use));
    cosineMin(i)=min(T.cosine(use));cosineMax(i)=max(T.cosine(use));
end
fig=figure('Visible','off','Color','w','Position',[100 100 1450 860]);tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
nexttile;
for seed=seeds.',use=T.seed==seed;plot(T.nominal_Akp(use),100*T.eta33_relative(use),'-o','LineWidth',1.3,'DisplayName',sprintf('seed %.0f',seed));hold on;end
plot(A.nominal_Akp,100*A.median_error,'k-o','LineWidth',2.2,'MarkerSize',7,'DisplayName','median');grid on;box on;
xlabel('Nominal A k_p scale');ylabel('Relative L_2 error (%)');title('GL8 versus HOS third phase sector');legend('Location','northwest');
nexttile;
errorbar(A.actual_kp_Hs_over_2,cosineMedian,cosineMedian-cosineMin,cosineMax-cosineMedian,'o-','LineWidth',1.8);grid on;box on;
xlabel('Actual k_p H_s / 2');ylabel('Normalized inner product');title('GL--HOS waveform correlation');yline(0,'k:','HandleVisibility','off');
nexttile;
for seed=seeds.'
    use=T.seed==seed;pH=polyfit(log(T.actual_kp_Hs_over_2(use)),log(T.HOS_third_norm(use)),1);
    pG=polyfit(log(T.actual_kp_Hs_over_2(use)),log(T.GL_eta33_norm(use)),1);
    loglog(T.actual_kp_Hs_over_2(use),T.HOS_third_norm(use),'-o','LineWidth',1.3,'DisplayName',sprintf('HOS %.0f, p=%.2f',seed,pH(1)));hold on;
    loglog(T.actual_kp_Hs_over_2(use),T.GL_eta33_norm(use),'--','LineWidth',1.2,'DisplayName',sprintf('GL %.0f, p=%.2f',seed,pG(1)));
end
grid on;box on;xlabel('Actual k_p H_s / 2');ylabel('L_2 norm (m)');title('HOS sector is not cubic; GL remains cubic');legend('Location','eastoutside');
nexttile;plot(A.actual_kp_Hs_over_2,A.median_norm_ratio,'k-o','LineWidth',2.1);hold on;
fill([A.actual_kp_Hs_over_2;flipud(A.actual_kp_Hs_over_2)],[A.min_norm_ratio;flipud(A.max_norm_ratio)],[.8 .85 1],'FaceAlpha',.3,'EdgeColor','none');
grid on;box on;xlabel('Actual k_p H_s / 2');ylabel('||GL eta_{33}||_2 / ||HOS third sector||_2');title('Median and three-seed range');
sgtitle('Unidirectional random-phase eta_{33}, k_p h=1');
out=fullfile(resultRoot,'analysis');exportgraphics(fig,fullfile(out,'random_eta33_final.png'),'Resolution',180);
exportgraphics(fig,fullfile(out,'random_eta33_final.pdf'),'ContentType','vector');close(fig);
summary=table(A.nominal_Akp,A.actual_kp_Hs_over_2,A.median_error,A.min_error,A.max_error,cosineMedian,cosineMin,cosineMax, ...
    'VariableNames',{'nominal_Akp','actual_kp_Hs_over_2','median_error','min_error','max_error','median_cosine','min_cosine','max_cosine'});
writetable(summary,fullfile(out,'random_eta33_final_summary.csv'));disp(summary);
end
