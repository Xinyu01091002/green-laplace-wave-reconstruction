function analyze_hos_eta33_rank_steepness(runRoot)
% Fixed-input GL4/6/8 rank ladder for the six-point HOS steepness family.
arguments
    runRoot (1,:) char
end
here=fileparts(mfilename('fullpath'));addpath(here);
amplitudes=[.02 .04 .06 .08 .10 .12];ranks=[4 6 8];rows=cell(0,10);
for icase=1:numel(amplitudes)
    akp=amplitudes(icase);name=sprintf('akp%03d',round(100*akp));
    d=load(fullfile(runRoot,name,'initial.mat'),'r');settings=d.r;raw=[];fulltime=[];
    for phase=1:4
        file=fullfile(runRoot,name,sprintf('phase%03d',(phase-1)*90),'Results','probes.dat');
        fid=fopen(file);assert(fid>=0);
        while ~feof(fid),line=fgetl(fid);if startsWith(strtrim(line),'VARIABLES'),break;end,end
        values=fscanf(fid,'%f',[2,Inf]).';fclose(fid);assert(all(isfinite(values),'all'));
        if phase==1,fulltime=values(:,1);raw=zeros(numel(fulltime),4);
        else,assert(max(abs(fulltime-values(:,1)))<1e-10);end
        raw(:,phase)=values(:,2); %#ok<AGROW>
    end
    ht=imag(analytic_signal(raw));firstFull=(raw(:,1)-raw(:,3)-ht(:,2)+ht(:,4))/4;
    thirdFull=(raw(:,1)-raw(:,3)+ht(:,2)-ht(:,4))/4;
    crop=fulltime>=32*settings.Tp-1e-8 & fulltime<=48*settings.Tp+1e-8;
    time=fulltime(crop);first=firstFull(crop);third=thirdFull(crop);t=time-time(1);
    wp=sqrt(settings.g*settings.kp*tanh(settings.kp*settings.h));results=cell(1,numel(ranks));
    for irank=1:numel(ranks)
        J=ranks(irank);options=struct('order',3,'omega_max',4*wp,'peak_wavenumber',settings.kp, ...
            'quadrature_rank',J,'domain_lengths',[25 50 100 200 400 800 1600], ...
            'relative_tolerance',5e-4,'memory_budget_MiB',4096);
        timer=tic;[prediction,audit]=gl_unidirectional_time_series(first,t,settings.g,settings.h,options);seconds=toc(timer);
        [~,peak]=max(abs(analytic_signal(prediction.eta1_used)));main=abs(t-t(peak))<=2*settings.Tp;
        last=audit.levels(end);results{irank}=prediction.eta33;
        rows(end+1,:)={akp,J,norm(prediction.eta33(main)-third(main))/norm(third(main)), ...
            norm(prediction.eta33-third)/norm(third),norm(prediction.eta33(main))/norm(third(main)), ...
            last.relative_change,last.L_over_h,last.Nx,seconds,audit.converged}; %#ok<AGROW>
    end
    caseOut=fullfile(runRoot,'rank_analysis',name);mkdir(caseOut);
    save(fullfile(caseOut,'rank_fields.mat'),'akp','time','first','third','main','ranks','results');
end
metrics=cell2table(rows,'VariableNames',{'Akp','GL_rank','main_relative','full_relative','main_norm_ratio', ...
    'last_domain_change','final_L_over_h','final_Nx','seconds','converged'});
for i=1:height(metrics)
    reference=metrics.main_relative(metrics.Akp==metrics.Akp(i)&metrics.GL_rank==8);
    metrics.error_change_from_GL8(i)=metrics.main_relative(i)-reference;
end
out=fullfile(runRoot,'rank_analysis');writetable(metrics,fullfile(out,'eta33_rank_metrics.csv'));
fig=figure('Visible','off','Color','w','Position',[100 100 1100 480]);tiledlayout(1,2,'Padding','compact');
nexttile;
for J=ranks,mask=metrics.GL_rank==J;semilogy(metrics.Akp(mask),100*metrics.main_relative(mask),'-o','LineWidth',1.8,'DisplayName',sprintf('GL%d',J));hold on;end
grid on;box on;xlabel('Steepness A k_p');ylabel('Main-window relative L_2 error (%)');title('Rank ladder against HOS');legend('Location','northwest');
nexttile;
for J=ranks(1:2),mask=metrics.GL_rank==J;plot(metrics.Akp(mask),100*metrics.error_change_from_GL8(mask),'-o','LineWidth',1.8,'DisplayName',sprintf('GL%d error - GL8 error',J));hold on;end
yline(0,'k:','HandleVisibility','off');grid on;box on;xlabel('Steepness A k_p');ylabel('Percentage-point difference');title('Effect of GL rank');legend('Location','best');
sgtitle('Single-direction eta_{33}: fixed-input GL rank diagnostic');
exportgraphics(fig,fullfile(out,'eta33_rank_ladder.png'),'Resolution',180);
exportgraphics(fig,fullfile(out,'eta33_rank_ladder.pdf'),'ContentType','vector');close(fig);
disp(metrics);
end

function z=analytic_signal(y)
N=size(y,1);mask=zeros(N,1);mask(1)=1;
if mod(N,2),mask(2:(N+1)/2)=2;else,mask(2:N/2)=2;mask(N/2+1)=1;end
z=ifft(fft(y).*mask);
end
