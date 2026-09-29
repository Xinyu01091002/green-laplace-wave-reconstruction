function analyze_hos_random_eta33_amplitude_order(runRoot)
% Extract the cubic-in-amplitude coefficient from HOS and GL random-wave records.
arguments
    runRoot (1,:) char
end
d=load(fullfile(runRoot,'analysis','random_eta33_fields.mat'),'fields');fields=d.fields;
amplitudes=(.02:.02:.18)';scaled=amplitudes/max(amplitudes);seeds=[20260925 20260926 20260927];
X3=[scaled,scaled.^3,scaled.^5];X4=[scaled,scaled.^3,scaled.^5,scaled.^7];rows=cell(numel(seeds),13);
out=fullfile(runRoot,'amplitude_order');assert(~isfolder(out));mkdir(out);
for iseed=1:numel(seeds)
    index=(iseed-1)*numel(amplitudes)+(1:numel(amplitudes));sample=fields{index(1)};N=numel(sample.t);
    hos=zeros(numel(amplitudes),N);gl=hos;first=hos;
    for i=1:numel(index)
        hos(i,:)=fields{index(i)}.third(:).';gl(i,:)=fields{index(i)}.eta33(:).';first(i,:)=fields{index(i)}.first(:).';
    end
    cH3=X3\hos;cG3=X3\gl;cF3=X3\first;cH4=X4\hos;cG4=X4\gl;cF4=X4\first;
    cubicH=cH4(2,:);cubicG=cG4(2,:);cosine=dot(cubicH,cubicG)/(norm(cubicH)*norm(cubicG));
    relative=norm(cubicG-cubicH)/norm(cubicH);ratio=norm(cubicG)/norm(cubicH);
    fitH=norm(X4*cH4-hos,'fro')/norm(hos,'fro');fitG=norm(X4*cG4-gl,'fro')/norm(gl,'fro');
    sensitivityH=norm(cH4(2,:)-cH3(2,:))/norm(cH4(2,:));sensitivityG=norm(cG4(2,:)-cG3(2,:))/norm(cG4(2,:));
    linearH=norm(cH4(1,:))/norm(cH4(2,:));linearG=norm(cG4(1,:))/norm(cG4(2,:));
    firstCubic=norm(cF4(2,:))/norm(cF4(1,:));
    rows(iseed,:)={seeds(iseed),relative,cosine,ratio,fitH,fitG,sensitivityH,sensitivityG,linearH,linearG,firstCubic,cond(X3),cond(X4)};
    save(fullfile(out,sprintf('seed%d_coefficients.mat',seeds(iseed))),'amplitudes','scaled','X3','X4','cH3','cG3','cF3','cH4','cG4','cF4','-v7.3');
    fig=figure('Visible','off','Color','w','Position',[100 100 1150 650]);tiledlayout(2,1,'Padding','compact');
    nexttile;plot(sample.t,cubicH,'k-','DisplayName','HOS cubic coefficient');hold on;plot(sample.t,cubicG,'--','DisplayName','GL cubic coefficient');grid on;
    ylabel('Coefficient (m at scaled amplitude 1)');title(sprintf('Seed %d: cubic coefficient, relative %.2f%%, cosine %.3f',seeds(iseed),100*relative,cosine));legend('Location','best');
    nexttile;plot(sample.t,cH4(1,:),'DisplayName','HOS linear leakage');hold on;plot(sample.t,cH4(3,:),'DisplayName','HOS fifth coefficient');plot(sample.t,cH4(4,:),'DisplayName','HOS seventh coefficient');grid on;
    xlabel('Time from 10 T_p (s)');ylabel('Coefficient (m)');legend('Location','best');
    exportgraphics(fig,fullfile(out,sprintf('seed%d_amplitude_order.png',seeds(iseed))),'Resolution',160);close(fig);
end
metrics=cell2table(rows,'VariableNames',{'seed','cubic_relative','cubic_cosine','cubic_norm_ratio','HOS_fit_residual','GL_fit_residual', ...
    'HOS_cubic_basis_sensitivity','GL_cubic_basis_sensitivity','HOS_linear_to_cubic_norm','GL_linear_to_cubic_norm', ...
    'first_cubic_to_linear_norm','condition_X3','condition_X4'});
writetable(metrics,fullfile(out,'amplitude_order_metrics.csv'));
report=struct('amplitude_coordinate','nominal Akp divided by 0.18','basis3','a,a^3,a^5','basis4','a,a^3,a^5,a^7', ...
    'primary','four-term all-nine-point least squares','sensitivity','relative change of cubic coefficient between three- and four-term bases', ...
    'no_interaction_enumerator',true);
fid=fopen(fullfile(out,'amplitude_order_report.json'),'w');fprintf(fid,'%s\n',jsonencode(report));fclose(fid);disp(metrics);disp(report);
end
