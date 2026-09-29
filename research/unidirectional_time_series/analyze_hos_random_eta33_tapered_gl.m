function analyze_hos_random_eta33_tapered_gl(runRoot)
% Recompute GL from a 5Tp-tapered full record and score the central 15--35Tp.
arguments
    runRoot (1,:) char
end
amplitudes=.02:.02:.18;seeds=[20260925 20260926 20260927];rows=cell(0,15);fields=cell(0,1);
for seed=seeds
    for akp=amplitudes
        folder=fullfile(runRoot,sprintf('seed%d',seed),sprintf('akp%03d',round(100*akp)));
        d=load(fullfile(folder,'initial.mat'),'report');r=d.report;raw=[];fulltime=[];
        for phase=1:4
            file=fullfile(folder,sprintf('phase%03d',(phase-1)*90),'Results','probes.dat');fid=fopen(file);assert(fid>=0);
            while ~feof(fid),line=fgetl(fid);if startsWith(strtrim(line),'VARIABLES'),break;end,end
            values=fscanf(fid,'%f',[2,Inf]).';fclose(fid);assert(all(isfinite(values),'all'));
            if phase==1,fulltime=values(:,1);raw=zeros(numel(fulltime),4);else,assert(max(abs(fulltime-values(:,1)))<1e-10);end
            raw(:,phase)=values(:,2); %#ok<AGROW>
        end
        t=fulltime-fulltime(1);window=cosine_taper(t,50*r.Tp,5*r.Tp);raw=raw.*window;ht=imag(analytic_signal(raw));
        first=(raw(:,1)-raw(:,3)-ht(:,2)+ht(:,4))/4;third=(raw(:,1)-raw(:,3)+ht(:,2)-ht(:,4))/4;
        wp=sqrt(r.g*r.kp*tanh(r.kp*r.h));options=struct('order',3,'omega_max',4*wp,'peak_wavenumber',r.kp, ...
            'quadrature_rank',8,'domain_lengths',[25 50 100 200 400 800 1600], ...
            'relative_tolerance',5e-4,'memory_budget_MiB',16384);
        timer=tic;[prediction,audit]=gl_unidirectional_time_series(first,t,r.g,r.h,options);seconds=toc(timer);
        score=t>=15*r.Tp-1e-8&t<=35*r.Tp+1e-8;gl=prediction.eta33(score);hos=third(score);
        cosine=dot(gl,hos)/(norm(gl)*norm(hos));error=norm(gl-hos)/norm(hos);ratio=norm(gl)/norm(hos);last=audit.levels(end);
        rows(end+1,:)={seed,akp,r.actual_kp_Hs_over_2,error,cosine,ratio,norm(hos),norm(gl), ...
            audit.input_projection_relative,last.relative_change,last.L_over_h,last.Nx,seconds,audit.converged,5}; %#ok<AGROW>
        fields{end+1}=struct('seed',seed,'Akp',akp,'actual_steepness',r.actual_kp_Hs_over_2,'t',t(score), ...
            'HOS_third',hos,'GL_eta33',gl,'first',first,'taper',window,'audit',audit); %#ok<AGROW>
    end
end
metrics=cell2table(rows,'VariableNames',{'seed','nominal_Akp','actual_kp_Hs_over_2','eta33_relative','cosine','norm_ratio', ...
    'HOS_third_norm','GL_eta33_norm','input_projection_relative','last_domain_change','final_L_over_h','final_Nx','GL_seconds','converged','taper_width_Tp'});
summary=cell(numel(amplitudes),11);
for i=1:numel(amplitudes)
    use=metrics.nominal_Akp==amplitudes(i);summary(i,:)={amplitudes(i),metrics.actual_kp_Hs_over_2(find(use,1)), ...
        median(metrics.eta33_relative(use)),min(metrics.eta33_relative(use)),max(metrics.eta33_relative(use)), ...
        median(metrics.cosine(use)),min(metrics.cosine(use)),max(metrics.cosine(use)), ...
        median(metrics.norm_ratio(use)),min(metrics.norm_ratio(use)),max(metrics.norm_ratio(use))};
end
aggregate=cell2table(summary,'VariableNames',{'nominal_Akp','actual_kp_Hs_over_2','median_error','min_error','max_error', ...
    'median_cosine','min_cosine','max_cosine','median_norm_ratio','min_norm_ratio','max_norm_ratio'});
out=fullfile(runRoot,'tapered_gl');assert(~isfolder(out));mkdir(out);writetable(metrics,fullfile(out,'tapered_gl_metrics.csv'));writetable(aggregate,fullfile(out,'tapered_gl_aggregate.csv'));
save(fullfile(out,'tapered_gl_fields.mat'),'metrics','aggregate','fields','-v7.3');
fig=figure('Visible','off','Color','w','Position',[100 100 1350 760]);tiledlayout(2,2,'Padding','compact');
nexttile;
for seed=seeds,use=metrics.seed==seed;plot(metrics.actual_kp_Hs_over_2(use),100*metrics.eta33_relative(use),'-o','DisplayName',sprintf('seed %d',seed));hold on;end
grid on;xlabel('Actual k_p H_s/2');ylabel('Relative L_2 error (%)');title('Tapered GL8 versus tapered HOS eta_{33} sector');legend('Location','best');
nexttile;errorbar(aggregate.actual_kp_Hs_over_2,aggregate.median_cosine,aggregate.median_cosine-aggregate.min_cosine,aggregate.max_cosine-aggregate.median_cosine,'o-','LineWidth',1.6);grid on;
xlabel('Actual k_p H_s/2');ylabel('Normalized inner product');title('Median and three-seed correlation range');
nexttile;errorbar(aggregate.actual_kp_Hs_over_2,aggregate.median_norm_ratio,aggregate.median_norm_ratio-aggregate.min_norm_ratio,aggregate.max_norm_ratio-aggregate.median_norm_ratio,'o-','LineWidth',1.6);grid on;
xlabel('Actual k_p H_s/2');ylabel('||GL||/||HOS||');title('Norm ratio');
nexttile;plot(metrics.actual_kp_Hs_over_2,100*metrics.input_projection_relative,'o');hold on;plot(metrics.actual_kp_Hs_over_2,100*metrics.last_domain_change,'x');grid on;
xlabel('Actual k_p H_s/2');ylabel('Percent');legend('First-input projection','Last domain change','Location','best');title('Numerical diagnostics');
sgtitle('Random-phase eta_{33} after 5T_p edge taper; score 15--35T_p');
exportgraphics(fig,fullfile(out,'tapered_gl_result.png'),'Resolution',180);exportgraphics(fig,fullfile(out,'tapered_gl_result.pdf'),'ContentType','vector');close(fig);
disp(aggregate);
end

function window=cosine_taper(t,duration,width)
window=ones(size(t));left=t<width;window(left)=.5*(1-cos(pi*t(left)/width));right=t>duration-width;
window(right)=.5*(1-cos(pi*(duration-t(right))/width));window=max(0,min(1,window));
end

function z=analytic_signal(y)
N=size(y,1);mask=zeros(N,1);mask(1)=1;
if mod(N,2),mask(2:(N+1)/2)=2;else,mask(2:N/2)=2;mask(N/2+1)=1;end
z=ifft(fft(y).*mask);
end
