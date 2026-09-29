function analyze_hos_eta33_steepness(runRoot)
% Analyze a fixed-spectrum six-point steepness ladder without interaction enumeration.
arguments
    runRoot (1,:) char
end
here=fileparts(mfilename('fullpath'));addpath(here);
amplitudes=[.02 .04 .06 .08 .10 .12];
rows=cell(numel(amplitudes),15);fields=cell(numel(amplitudes),1);
for icase=1:numel(amplitudes)
    akp=amplitudes(icase);name=sprintf('akp%03d',round(100*akp));
    d=load(fullfile(runRoot,name,'initial.mat'),'r');r=d.r;raw=[];fulltime=[];
    for phase=1:4
        file=fullfile(runRoot,name,sprintf('phase%03d',(phase-1)*90),'Results','probes.dat');
        fid=fopen(file);assert(fid>=0,'Missing %s',file);
        while ~feof(fid),line=fgetl(fid);if startsWith(strtrim(line),'VARIABLES'),break;end,end
        values=fscanf(fid,'%f',[2,Inf]).';fclose(fid);assert(size(values,1)>=2001&&all(isfinite(values),'all'));
        if phase==1,fulltime=values(:,1);raw=zeros(numel(fulltime),4);
        else,assert(max(abs(fulltime-values(:,1)))<1e-10);end
        raw(:,phase)=values(:,2); %#ok<AGROW> Preallocated to four columns when phase==1.
    end
    ht=imag(analytic_signal(raw));
    firstFull=(raw(:,1)-raw(:,3)-ht(:,2)+ht(:,4))/4;
    thirdFull=(raw(:,1)-raw(:,3)+ht(:,2)-ht(:,4))/4;
    crop=fulltime>=32*r.Tp-1e-8 & fulltime<=48*r.Tp+1e-8;
    time=fulltime(crop);first=firstFull(crop);third=thirdFull(crop);t=time-time(1);
    assert(numel(t)==641);wp=sqrt(r.g*r.kp*tanh(r.kp*r.h));
    options=struct('order',3,'omega_max',4*wp,'peak_wavenumber',r.kp, ...
        'quadrature_rank',8,'domain_lengths',[25 50 100 200 400 800 1600], ...
        'relative_tolerance',5e-4,'memory_budget_MiB',4096);
    timer=tic;[prediction,audit]=gl_unidirectional_time_series(first,t,r.g,r.h,options);seconds=toc(timer);
    [~,peak]=max(abs(analytic_signal(prediction.eta1_used)));main=abs(t-t(peak))<=2*r.Tp;
    mainError=norm(prediction.eta33(main)-third(main))/norm(third(main));
    fullError=norm(prediction.eta33-third)/norm(third);
    normRatio=norm(prediction.eta33(main))/norm(third(main));
    levels=audit.levels;last=levels(end);
    rows(icase,:)={akp,numel(t),numel(audit.parent_bins),audit.input_projection_relative, ...
        mainError,fullError,normRatio,norm(third(main)),norm(prediction.eta33(main)), ...
        last.relative_change,last.L_over_h,last.Nx,seconds,audit.converged,r.initial_eta_max};
    fields{icase}=struct('akp',akp,'t',t,'time',time,'first',first,'third',third, ...
        'eta1_used',prediction.eta1_used,'eta33',prediction.eta33,'main',main,'audit',audit);
end
metrics=cell2table(rows,'VariableNames',{'Akp','samples','parent_bins','input_projection_relative', ...
    'eta33_main_relative','eta33_full_relative','main_norm_ratio','HOS_main_norm','GL_main_norm', ...
    'last_domain_change','final_L_over_h','final_Nx','GL_seconds','converged','initial_eta_max'});

base=fields{1};firstChange=zeros(height(metrics),1);thirdChange=firstChange;glChange=firstChange;
for i=1:height(metrics)
    scale=metrics.Akp(i)/metrics.Akp(1);z=fields{i};
    firstChange(i)=norm(z.first/scale-base.first)/norm(base.first);
    thirdChange(i)=norm(z.third/scale^3-base.third)/norm(base.third);
    glChange(i)=norm(z.eta33/scale^3-base.eta33)/norm(base.eta33);
end
metrics.first_over_A_change=firstChange;metrics.HOS_third_over_A3_change=thirdChange;
metrics.GL_eta33_over_A3_change=glChange;

pError=polyfit(log(metrics.Akp),log(metrics.eta33_main_relative),1);
pHos=polyfit(log(metrics.Akp),log(metrics.HOS_main_norm),1);
pGl=polyfit(log(metrics.Akp),log(metrics.GL_main_norm),1);
guide=metrics.eta33_main_relative(1)*(metrics.Akp/metrics.Akp(1)).^2;
report=struct('error_power_all_six',pError(1),'HOS_third_norm_power',pHos(1), ...
    'GL_eta33_norm_power',pGl(1),'expected_relative_contamination_power',2, ...
    'harmonic_extraction','complete 50Tp Hilbert-four-phase extraction, then fixed 32--48Tp crop', ...
    'GL','non-enumerating gl_unidirectional_time_series, GL8','no_fitted_adjustments',true);

out=fullfile(runRoot,'analysis');assert(~isfolder(out));mkdir(out);
writetable(metrics,fullfile(out,'eta33_steepness_metrics.csv'));
fid=fopen(fullfile(out,'eta33_steepness_report.json'),'w');fprintf(fid,'%s\n',jsonencode(report));fclose(fid);
save(fullfile(out,'eta33_steepness_fields.mat'),'metrics','fields','report','-v7.3');

fig=figure('Visible','off','Color','w','Position',[100 100 1550 520]);
tiledlayout(1,3,'Padding','compact','TileSpacing','compact');
nexttile;loglog(metrics.Akp,100*metrics.eta33_main_relative,'-o','LineWidth',2,'MarkerSize',7);hold on;
loglog(metrics.Akp,100*guide,'--','LineWidth',1.7);grid on;box on;
xlabel('Steepness A k_p');ylabel('Main-window relative L_2 error (%)');
title(sprintf('GL--HOS eta_{33}, slope %.2f',pError(1)));legend('Measured','O((Ak_p)^2) guide','Location','northwest');
nexttile;plot(metrics.Akp,metrics.main_norm_ratio,'-o','LineWidth',2,'MarkerSize',7);grid on;box on;
xlabel('Steepness A k_p');ylabel('||GL||_2 / ||HOS third sector||_2');title('Main-window norm ratio');
nexttile;semilogy(metrics.Akp,max(firstChange,eps),'-o','LineWidth',1.8);hold on;
semilogy(metrics.Akp,max(thirdChange,eps),'-s','LineWidth',1.8);
semilogy(metrics.Akp,max(glChange,eps),'-d','LineWidth',1.8);grid on;box on;
xlabel('Steepness A k_p');ylabel('Change relative to Ak_p=0.02');
title('Departure after amplitude normalization');
legend('eta_1/A','HOS third/A^3','GL eta_{33}/A^3','Location','northwest');
sgtitle('Single-direction eta_{33} steepness ladder, k_p h=1');
exportgraphics(fig,fullfile(out,'eta33_steepness_ladder.png'),'Resolution',180);
exportgraphics(fig,fullfile(out,'eta33_steepness_ladder.pdf'),'ContentType','vector');close(fig);
disp(metrics);disp(report);
end

function z=analytic_signal(y)
N=size(y,1);mask=zeros(N,1);mask(1)=1;
if mod(N,2),mask(2:(N+1)/2)=2;else,mask(2:N/2)=2;mask(N/2+1)=1;end
z=ifft(fft(y).*mask);
end
