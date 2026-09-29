function diagnose_random_tapered_harmonics(runRoot)
% Test raised-cosine edge tapers before Hilbert--four-phase extraction.
arguments
    runRoot (1,:) char
end
amplitudes=.02:.02:.18;seeds=[20260925 20260926 20260927];widths=[0 2 5 10];rows=cell(0,11);
saved=load(fullfile(runRoot,'analysis','random_eta33_fields.mat'),'fields');fields=saved.fields;
for width=widths
    for iseed=1:numel(seeds)
        for ia=1:numel(amplitudes)
            akp=amplitudes(ia);folder=fullfile(runRoot,sprintf('seed%d',seeds(iseed)),sprintf('akp%03d',round(100*akp)));
            d=load(fullfile(folder,'initial.mat'),'A','k','report');r=d.report;N=2001;t=(0:N-1)'*r.output_dt;
            omega=sqrt(r.g*d.k.*tanh(d.k*r.h));base=exp(1i*(r.probe_x*d.k-t*omega))*d.A.';linear=zeros(N,4);
            for phase=1:4,linear(:,phase)=real(base*exp(1i*(phase-1)*pi/2));end
            raw=[];fulltime=[];
            for phase=1:4
                file=fullfile(folder,sprintf('phase%03d',(phase-1)*90),'Results','probes.dat');fid=fopen(file);assert(fid>=0);
                while ~feof(fid),line=fgetl(fid);if startsWith(strtrim(line),'VARIABLES'),break;end,end
                values=fscanf(fid,'%f',[2,Inf]).';fclose(fid);
                if phase==1,fulltime=values(:,1);raw=zeros(numel(fulltime),4);else,assert(max(abs(fulltime-values(:,1)))<1e-10);end
                raw(:,phase)=values(:,2); %#ok<AGROW>
            end
            window=cosine_taper(t,50*r.Tp,width*r.Tp);linear=linear.*window;raw=raw.*window;
            linearHt=imag(analytic_signal(linear));rawHt=imag(analytic_signal(raw));
            linearFirst=(linear(:,1)-linear(:,3)-linearHt(:,2)+linearHt(:,4))/4;
            linearThird=(linear(:,1)-linear(:,3)+linearHt(:,2)-linearHt(:,4))/4;
            hosThird=(raw(:,1)-raw(:,3)+rawHt(:,2)-rawHt(:,4))/4;
            score=t>=15*r.Tp-1e-8&t<=35*r.Tp+1e-8;leak=linearThird(score);hos=hosThird(score);
            field=fields{(iseed-1)*numel(amplitudes)+ia};scoreGl=field.t>=5*r.Tp-1e-8&field.t<=25*r.Tp+1e-8;gl=field.eta33(scoreGl);
            assert(numel(gl)==nnz(score));
            cosine=dot(leak,hos)/(norm(leak)*norm(hos));glError=norm(gl-hos)/norm(hos);glCosine=dot(gl,hos)/(norm(gl)*norm(hos));
            rows(end+1,:)={width,seeds(iseed),akp,r.actual_kp_Hs_over_2,norm(leak)/norm(linearFirst(score)), ...
                norm(leak)/norm(hos),cosine,norm(hos),glError,glCosine,norm(gl)/norm(hos)}; %#ok<AGROW>
        end
    end
end
metrics=cell2table(rows,'VariableNames',{'taper_width_Tp','seed','nominal_Akp','actual_kp_Hs_over_2', ...
    'linear_control_third_to_first','linear_leak_to_HOS_third','linear_leak_HOS_cosine','HOS_third_norm', ...
    'existing_GL_error','existing_GL_cosine','existing_GL_norm_ratio'});
fitRows=cell(0,6);
for width=widths
    for seed=seeds
        use=metrics.taper_width_Tp==width&metrics.seed==seed;
        power=polyfit(log(metrics.actual_kp_Hs_over_2(use)),log(metrics.HOS_third_norm(use)),1);
        fitRows(end+1,:)={width,seed,power(1),median(metrics.linear_leak_to_HOS_third(use)), ...
            median(metrics.existing_GL_error(use)),median(metrics.existing_GL_cosine(use))}; %#ok<AGROW>
    end
end
fits=cell2table(fitRows,'VariableNames',{'taper_width_Tp','seed','HOS_third_power', ...
    'median_linear_leak_fraction','median_existing_GL_error','median_existing_GL_cosine'});
out=fullfile(runRoot,'taper_diagnostic');assert(~isfolder(out));mkdir(out);writetable(metrics,fullfile(out,'taper_metrics.csv'));writetable(fits,fullfile(out,'taper_fits.csv'));
fig=figure('Visible','off','Color','w','Position',[100 100 1350 760]);tiledlayout(2,2,'Padding','compact');
nexttile;
for seed=seeds,use=fits.seed==seed;plot(fits.taper_width_Tp(use),fits.HOS_third_power(use),'-o','DisplayName',sprintf('seed %d',seed));hold on;end
yline(3,'k--','cubic','HandleVisibility','off');grid on;xlabel('Taper width per edge / T_p');ylabel('HOS third-sector amplitude power');title('Does taper restore cubic scaling?');legend('Location','best');
nexttile;
for seed=seeds,use=fits.seed==seed;semilogy(fits.taper_width_Tp(use),fits.median_linear_leak_fraction(use),'-o','DisplayName',sprintf('seed %d',seed));hold on;end
grid on;xlabel('Taper width per edge / T_p');ylabel('Median ||linear leakage|| / ||HOS third||');title('Pure-linear leakage');legend('Location','best');
nexttile;
for seed=seeds,use=fits.seed==seed;plot(fits.taper_width_Tp(use),100*fits.median_existing_GL_error(use),'-o','DisplayName',sprintf('seed %d',seed));hold on;end
grid on;xlabel('Taper width per edge / T_p');ylabel('Median existing-GL error (%)');title('Central 15--35 T_p diagnostic');legend('Location','best');
nexttile;
for seed=seeds,use=fits.seed==seed;plot(fits.taper_width_Tp(use),fits.median_existing_GL_cosine(use),'-o','DisplayName',sprintf('seed %d',seed));hold on;end
grid on;xlabel('Taper width per edge / T_p');ylabel('Median GL--HOS cosine');title('Central waveform correlation');legend('Location','best');
sgtitle('Raised-cosine taper diagnostic for random-wave Hilbert separation');
exportgraphics(fig,fullfile(out,'taper_diagnostic.png'),'Resolution',180);exportgraphics(fig,fullfile(out,'taper_diagnostic.pdf'),'ContentType','vector');close(fig);
disp(fits);
end

function window=cosine_taper(t,duration,width)
window=ones(size(t));if width==0,return,end
left=t<width;window(left)=.5*(1-cos(pi*t(left)/width));right=t>duration-width;
window(right)=.5*(1-cos(pi*(duration-t(right))/width));window=max(0,min(1,window));
end

function z=analytic_signal(y)
N=size(y,1);mask=zeros(N,1);mask(1)=1;
if mod(N,2),mask(2:(N+1)/2)=2;else,mask(2:N/2)=2;mask(N/2+1)=1;end
z=ifft(fft(y).*mask);
end
