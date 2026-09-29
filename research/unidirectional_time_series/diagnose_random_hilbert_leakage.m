function diagnose_random_hilbert_leakage(runRoot)
% Test finite-record Hilbert leakage with an exact linear four-phase control.
arguments
    runRoot (1,:) char
end
amplitudes=.02:.02:.18;seeds=[20260925 20260926 20260927];rows=cell(0,12);
saved=load(fullfile(runRoot,'analysis','random_eta33_fields.mat'),'fields');fields=saved.fields;
for iseed=1:numel(seeds)
    for ia=1:numel(amplitudes)
        akp=amplitudes(ia);folder=fullfile(runRoot,sprintf('seed%d',seeds(iseed)),sprintf('akp%03d',round(100*akp)));
        d=load(fullfile(folder,'initial.mat'),'A','k','report');r=d.report;N=2001;t=(0:N-1)'*r.output_dt;
        omega=sqrt(r.g*d.k.*tanh(d.k*r.h));base=exp(1i*(r.probe_x*d.k-t*omega))*d.A.';linear=zeros(N,4);
        for phase=1:4,linear(:,phase)=real(base*exp(1i*(phase-1)*pi/2));end
        linearHt=imag(analytic_signal(linear));linearFirst=(linear(:,1)-linear(:,3)-linearHt(:,2)+linearHt(:,4))/4;
        linearThird=(linear(:,1)-linear(:,3)+linearHt(:,2)-linearHt(:,4))/4;
        raw=[];fulltime=[];
        for phase=1:4
            file=fullfile(folder,sprintf('phase%03d',(phase-1)*90),'Results','probes.dat');fid=fopen(file);assert(fid>=0);
            while ~feof(fid),line=fgetl(fid);if startsWith(strtrim(line),'VARIABLES'),break;end,end
            values=fscanf(fid,'%f',[2,Inf]).';fclose(fid);
            if phase==1,fulltime=values(:,1);raw=zeros(numel(fulltime),4);else,assert(max(abs(fulltime-values(:,1)))<1e-10);end
            raw(:,phase)=values(:,2); %#ok<AGROW>
        end
        ht=imag(analytic_signal(raw));hosThird=(raw(:,1)-raw(:,3)+ht(:,2)-ht(:,4))/4;
        crop=t>=10*r.Tp-1e-8&t<=40*r.Tp+1e-8;leak=linearThird(crop);hos=hosThird(crop);
        field=fields{(iseed-1)*numel(amplitudes)+ia};gl=field.eta33(:);assert(numel(gl)==nnz(crop));corrected=hos-leak;
        leakCosine=dot(leak,hos)/(norm(leak)*norm(hos));rawError=norm(gl-hos)/norm(hos);
        correctedError=norm(gl-corrected)/norm(corrected);correctedCosine=dot(gl,corrected)/(norm(gl)*norm(corrected));
        rows(end+1,:)={seeds(iseed),akp,r.actual_kp_Hs_over_2,norm(leak)/norm(linearFirst(crop)), ...
            norm(leak)/norm(hos),leakCosine,rawError,correctedError,correctedCosine,norm(gl)/norm(corrected),norm(leak),norm(hos)}; %#ok<AGROW>
    end
end
metrics=cell2table(rows,'VariableNames',{'seed','nominal_Akp','actual_kp_Hs_over_2','linear_control_third_to_first', ...
    'linear_leak_to_HOS_third','linear_leak_HOS_cosine','raw_GL_error','leak_subtracted_GL_error', ...
    'leak_subtracted_GL_cosine','leak_subtracted_GL_norm_ratio','linear_leak_norm','HOS_third_norm'});
out=fullfile(runRoot,'hilbert_leakage');assert(~isfolder(out));mkdir(out);writetable(metrics,fullfile(out,'hilbert_leakage_metrics.csv'));
fig=figure('Visible','off','Color','w','Position',[100 100 1350 760]);tiledlayout(2,2,'Padding','compact');
nexttile;
for seed=seeds,use=metrics.seed==seed;loglog(metrics.actual_kp_Hs_over_2(use),metrics.linear_leak_norm(use),'-o','DisplayName',sprintf('linear control %d',seed));hold on;end
grid on;xlabel('Actual k_p H_s/2');ylabel('Spurious third-sector norm');title('Pure linear control should be zero');legend('Location','northwest');
nexttile;
for seed=seeds,use=metrics.seed==seed;plot(metrics.actual_kp_Hs_over_2(use),metrics.linear_leak_to_HOS_third(use),'-o','DisplayName',sprintf('seed %d',seed));hold on;end
grid on;xlabel('Actual k_p H_s/2');ylabel('||linear leakage|| / ||HOS third sector||');title('Leakage magnitude relative to HOS sector');legend('Location','best');
nexttile;
for seed=seeds,use=metrics.seed==seed;plot(metrics.actual_kp_Hs_over_2(use),100*metrics.raw_GL_error(use),'-o','DisplayName',sprintf('raw %d',seed));hold on;end
grid on;xlabel('Actual k_p H_s/2');ylabel('GL relative error (%)');title('Before leakage subtraction');legend('Location','best');
nexttile;
for seed=seeds,use=metrics.seed==seed;plot(metrics.actual_kp_Hs_over_2(use),100*metrics.leak_subtracted_GL_error(use),'-o','DisplayName',sprintf('corrected %d',seed));hold on;end
grid on;xlabel('Actual k_p H_s/2');ylabel('GL relative error (%)');title('Diagnostic subtraction of linear-control leakage');legend('Location','best');
sgtitle('Finite-record Hilbert leakage diagnostic');
exportgraphics(fig,fullfile(out,'hilbert_leakage_diagnostic.png'),'Resolution',180);exportgraphics(fig,fullfile(out,'hilbert_leakage_diagnostic.pdf'),'ContentType','vector');close(fig);
disp(metrics);
end

function z=analytic_signal(y)
N=size(y,1);mask=zeros(N,1);mask(1)=1;
if mod(N,2),mask(2:(N+1)/2)=2;else,mask(2:N/2)=2;mask(N/2+1)=1;end
z=ifft(fft(y).*mask);
end
