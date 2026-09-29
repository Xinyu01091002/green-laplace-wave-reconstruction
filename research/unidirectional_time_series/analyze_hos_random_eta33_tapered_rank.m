function analyze_hos_random_eta33_tapered_rank(runRoot,akp)
% GL4/6/8 rank diagnostic on the corrected tapered random-wave comparison.
arguments
    runRoot (1,:) char
    akp (1,1) double = .12
end
seeds=[20260925 20260926 20260927];ranks=[4 6];saved=load(fullfile(runRoot,'tapered_gl','tapered_gl_fields.mat'),'fields');
allAmplitudes=.02:.02:.18;ia=find(abs(allAmplitudes-akp)<1e-12,1);rows=cell(0,10);
baseMetrics=readtable(fullfile(runRoot,'tapered_gl','tapered_gl_metrics.csv'));
for iseed=1:numel(seeds)
    field=saved.fields{(iseed-1)*numel(allAmplitudes)+ia};folder=fullfile(runRoot,sprintf('seed%d',seeds(iseed)),sprintf('akp%03d',round(100*akp)));
    d=load(fullfile(folder,'initial.mat'),'report');r=d.report;raw=[];fulltime=[];
    for phase=1:4
        file=fullfile(folder,sprintf('phase%03d',(phase-1)*90),'Results','probes.dat');fid=fopen(file);assert(fid>=0);
        while ~feof(fid),line=fgetl(fid);if startsWith(strtrim(line),'VARIABLES'),break;end,end
        values=fscanf(fid,'%f',[2,Inf]).';fclose(fid);
        if phase==1,fulltime=values(:,1);raw=zeros(numel(fulltime),4);else,assert(max(abs(fulltime-values(:,1)))<1e-10);end
        raw(:,phase)=values(:,2); %#ok<AGROW>
    end
    t=fulltime-fulltime(1);window=cosine_taper(t,50*r.Tp,5*r.Tp);raw=raw.*window;ht=imag(analytic_signal(raw));
    first=(raw(:,1)-raw(:,3)-ht(:,2)+ht(:,4))/4;score=t>=15*r.Tp-1e-8&t<=35*r.Tp+1e-8;hos=field.HOS_third;
    wp=sqrt(r.g*r.kp*tanh(r.kp*r.h));
    for J=ranks
        options=struct('order',3,'omega_max',4*wp,'peak_wavenumber',r.kp,'quadrature_rank',J, ...
            'domain_lengths',[25 50 100 200 400 800 1600],'relative_tolerance',5e-4,'memory_budget_MiB',16384);
        timer=tic;[prediction,audit]=gl_unidirectional_time_series(first,t,r.g,r.h,options);seconds=toc(timer);gl=prediction.eta33(score);last=audit.levels(end);
        rows(end+1,:)={seeds(iseed),akp,J,norm(gl-hos)/norm(hos),dot(gl,hos)/(norm(gl)*norm(hos)),norm(gl)/norm(hos), ...
            last.relative_change,last.L_over_h,last.Nx,seconds}; %#ok<AGROW>
    end
    use=baseMetrics.seed==seeds(iseed)&abs(baseMetrics.nominal_Akp-akp)<1e-12;b=baseMetrics(use,:);
    rows(end+1,:)={seeds(iseed),akp,8,b.eta33_relative,b.cosine,b.norm_ratio,b.last_domain_change,b.final_L_over_h,b.final_Nx,b.GL_seconds}; %#ok<AGROW>
end
metrics=cell2table(rows,'VariableNames',{'seed','nominal_Akp','GL_rank','relative_error','cosine','norm_ratio','last_domain_change','final_L_over_h','final_Nx','seconds'});
out=fullfile(runRoot,'tapered_rank');assert(~isfolder(out));mkdir(out);writetable(metrics,fullfile(out,'tapered_rank_metrics.csv'));
fig=figure('Visible','off','Color','w','Position',[100 100 1100 460]);tiledlayout(1,2,'Padding','compact');
nexttile;
for seed=seeds,use=metrics.seed==seed;plot(metrics.GL_rank(use),100*metrics.relative_error(use),'-o','LineWidth',1.6,'DisplayName',sprintf('seed %d',seed));hold on;end
grid on;xticks([4 6 8]);xlabel('Shared inner/outer GL rank');ylabel('Relative L_2 error (%)');title('Tapered random-wave eta_{33} rank ladder');legend('Location','best');
nexttile;
for seed=seeds,use=metrics.seed==seed;plot(metrics.GL_rank(use),metrics.cosine(use),'-o','LineWidth',1.6,'DisplayName',sprintf('seed %d',seed));hold on;end
grid on;xticks([4 6 8]);xlabel('Shared inner/outer GL rank');ylabel('Normalized inner product');title('Waveform correlation');legend('Location','best');
sgtitle(sprintf('Fixed tapered input, nominal A k_p=%.2f',akp));exportgraphics(fig,fullfile(out,'tapered_rank.png'),'Resolution',180);exportgraphics(fig,fullfile(out,'tapered_rank.pdf'),'ContentType','vector');close(fig);
disp(metrics);
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
