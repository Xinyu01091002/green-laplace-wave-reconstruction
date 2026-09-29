function analyze_hos_random_eta33_steepness(runRoot)
% Non-enumerating GL8 analysis of fixed-amplitude random-phase HOS families.
arguments
    runRoot (1,:) char
end
here=fileparts(mfilename('fullpath'));addpath(here);
amplitudes=(.02:.02:.18);seeds=[20260925 20260926 20260927];rows=cell(0,17);fields=cell(0,1);
for seed=seeds
    for akp=amplitudes
        name=sprintf('akp%03d',round(100*akp));folder=fullfile(runRoot,sprintf('seed%d',seed),name);
        d=load(fullfile(folder,'initial.mat'),'report');settings=d.report;raw=[];fulltime=[];
        for phase=1:4
            file=fullfile(folder,sprintf('phase%03d',(phase-1)*90),'Results','probes.dat');
            fid=fopen(file);assert(fid>=0);
            while ~feof(fid),line=fgetl(fid);if startsWith(strtrim(line),'VARIABLES'),break;end,end
            values=fscanf(fid,'%f',[2,Inf]).';fclose(fid);assert(size(values,1)>=2001&&all(isfinite(values),'all'));
            if phase==1,fulltime=values(:,1);raw=zeros(numel(fulltime),4);
            else,assert(max(abs(fulltime-values(:,1)))<1e-10);end
            raw(:,phase)=values(:,2); %#ok<AGROW>
        end
        ht=imag(analytic_signal(raw));firstFull=(raw(:,1)-raw(:,3)-ht(:,2)+ht(:,4))/4;
        thirdFull=(raw(:,1)-raw(:,3)+ht(:,2)-ht(:,4))/4;
        crop=fulltime>=10*settings.Tp-1e-8 & fulltime<=40*settings.Tp+1e-8;
        time=fulltime(crop);first=firstFull(crop);third=thirdFull(crop);t=time-time(1);assert(numel(t)==1201);
        wp=sqrt(settings.g*settings.kp*tanh(settings.kp*settings.h));
        options=struct('order',3,'omega_max',4*wp,'peak_wavenumber',settings.kp, ...
            'quadrature_rank',8,'domain_lengths',[25 50 100 200 400 800 1600], ...
            'relative_tolerance',5e-4,'memory_budget_MiB',8192);
        timer=tic;[prediction,audit]=gl_unidirectional_time_series(first,t,settings.g,settings.h,options);seconds=toc(timer);
        last=audit.levels(end);error=norm(prediction.eta33-third)/norm(third);ratio=norm(prediction.eta33)/norm(third);
        rows(end+1,:)={seed,akp,settings.actual_kp_Hs_over_2,settings.linear_Hs_4sigma_m,numel(t), ...
            numel(audit.parent_bins),audit.input_projection_relative,error,ratio,norm(third),norm(prediction.eta33), ...
            last.relative_change,last.L_over_h,last.Nx,seconds,audit.converged,settings.initial_eta_max}; %#ok<AGROW>
        fields{end+1}=struct('seed',seed,'Akp',akp,'actual_steepness',settings.actual_kp_Hs_over_2, ...
            't',t,'first',first,'third',third,'eta1_used',prediction.eta1_used,'eta33',prediction.eta33,'audit',audit); %#ok<AGROW>
    end
end
metrics=cell2table(rows,'VariableNames',{'seed','nominal_Akp','actual_kp_Hs_over_2','linear_Hs_m','samples', ...
    'parent_bins','input_projection_relative','eta33_relative','norm_ratio','HOS_third_norm','GL_eta33_norm', ...
    'last_domain_change','final_L_over_h','final_Nx','GL_seconds','converged','initial_eta_max'});

summary=cell(numel(amplitudes),8);
for i=1:numel(amplitudes)
    use=metrics.nominal_Akp==amplitudes(i);values=metrics.eta33_relative(use);ratios=metrics.norm_ratio(use);
    summary(i,:)={amplitudes(i),metrics.actual_kp_Hs_over_2(find(use,1)),median(values),min(values),max(values), ...
        median(ratios),min(ratios),max(ratios)};
end
aggregate=cell2table(summary,'VariableNames',{'nominal_Akp','actual_kp_Hs_over_2','median_error','min_error','max_error', ...
    'median_norm_ratio','min_norm_ratio','max_norm_ratio'});
x=aggregate.actual_kp_Hs_over_2.^2;y=aggregate.median_error;X=[ones(size(x)),x];b=X\y;fitted=X*b;
R2=1-sum((y-fitted).^2)/sum((y-mean(y)).^2);
seedFits=cell(numel(seeds),4);
for i=1:numel(seeds)
    use=metrics.seed==seeds(i);xs=metrics.actual_kp_Hs_over_2(use).^2;ys=metrics.eta33_relative(use);bs=[ones(size(xs)),xs]\ys;
    yh=[ones(size(xs)),xs]*bs;seedFits(i,:)={seeds(i),bs(1),bs(2),1-sum((ys-yh).^2)/sum((ys-mean(ys)).^2)};
end
fits=cell2table(seedFits,'VariableNames',{'seed','E0','quadratic_coefficient','R2'});
report=struct('seeds',seeds,'nominal_Akp',amplitudes,'fit_coordinate','actual kp*Hs/2', ...
    'median_fit_E0',b(1),'median_fit_quadratic_coefficient',b(2),'median_fit_R2',R2, ...
    'scoring_window_Tp',[10 40],'harmonic_extraction','complete 50Tp Hilbert-four-phase extraction before crop', ...
    'GL','non-enumerating gl_unidirectional_time_series GL8','no_fitted_adjustments',true);
out=fullfile(runRoot,'analysis');assert(~isfolder(out));mkdir(out);
writetable(metrics,fullfile(out,'random_eta33_metrics.csv'));writetable(aggregate,fullfile(out,'random_eta33_aggregate.csv'));
writetable(fits,fullfile(out,'random_eta33_seed_fits.csv'));
fid=fopen(fullfile(out,'random_eta33_report.json'),'w');fprintf(fid,'%s\n',jsonencode(report));fclose(fid);
save(fullfile(out,'random_eta33_fields.mat'),'metrics','aggregate','fits','fields','report','-v7.3');

fig=figure('Visible','off','Color','w','Position',[100 100 1450 850]);tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
nexttile;
for seed=seeds,use=metrics.seed==seed;plot(metrics.nominal_Akp(use),100*metrics.eta33_relative(use),'-o','LineWidth',1.3,'DisplayName',sprintf('seed %d',seed));hold on;end
plot(aggregate.nominal_Akp,100*aggregate.median_error,'k-o','LineWidth',2.3,'MarkerSize',7,'DisplayName','median');grid on;box on;
xlabel('Nominal A k_p scale');ylabel('Relative L_2 error (%)');title('Random-phase eta_{33} error');legend('Location','northwest');
nexttile;
errorbar(aggregate.actual_kp_Hs_over_2,100*aggregate.median_error, ...
    100*(aggregate.median_error-aggregate.min_error),100*(aggregate.max_error-aggregate.median_error),'o','LineWidth',1.6);hold on;
dense=linspace(min(aggregate.actual_kp_Hs_over_2),max(aggregate.actual_kp_Hs_over_2),200)';
plot(dense,100*(b(1)+b(2)*dense.^2),'-','LineWidth',2);grid on;box on;
xlabel('Actual k_p H_s / 2');ylabel('Median relative L_2 error (%)');title(sprintf('Median fit: E_0+c s^2, R^2=%.4f',R2));
nexttile;
for seed=seeds,use=metrics.seed==seed;loglog(metrics.actual_kp_Hs_over_2(use),metrics.HOS_third_norm(use),'-o','LineWidth',1.3,'DisplayName',sprintf('HOS seed %d',seed));hold on;end
grid on;box on;xlabel('Actual k_p H_s / 2');ylabel('HOS third-sector L_2 norm (m)');title('Third-sector amplitude scaling');legend('Location','northwest');
nexttile;plot(aggregate.actual_kp_Hs_over_2,aggregate.median_norm_ratio,'k-o','LineWidth',2.2);hold on;
fill([aggregate.actual_kp_Hs_over_2;flipud(aggregate.actual_kp_Hs_over_2)], ...
    [aggregate.min_norm_ratio;flipud(aggregate.max_norm_ratio)],[.8 .85 1],'FaceAlpha',.3,'EdgeColor','none');
grid on;box on;xlabel('Actual k_p H_s / 2');ylabel('||GL eta_{33}||_2 / ||HOS third sector||_2');title('Median and seed range');
sgtitle('Unidirectional random-phase eta_{33} steepness experiment, k_p h=1');
exportgraphics(fig,fullfile(out,'random_eta33_steepness.png'),'Resolution',180);
exportgraphics(fig,fullfile(out,'random_eta33_steepness.pdf'),'ContentType','vector');close(fig);
disp(aggregate);disp(fits);disp(report);
end

function z=analytic_signal(y)
N=size(y,1);mask=zeros(N,1);mask(1)=1;
if mod(N,2),mask(2:(N+1)/2)=2;else,mask(2:N/2)=2;mask(N/2+1)=1;end
z=ifft(fft(y).*mask);
end
