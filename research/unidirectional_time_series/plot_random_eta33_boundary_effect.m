function plot_random_eta33_boundary_effect(runRoot,akp)
% Plot direct time-series evidence of finite-record Hilbert boundary leakage.
arguments
    runRoot (1,:) char
    akp (1,1) double = .12
end
seeds=[20260925 20260926 20260927];amplitudes=.02:.02:.18;
saved=load(fullfile(runRoot,'tapered_gl','tapered_gl_fields.mat'),'fields');fields=saved.fields;
fig=figure('Visible','off','Color','w','Position',[50 50 1800 1100]);tiledlayout(3,3,'Padding','compact','TileSpacing','compact');
for iseed=1:numel(seeds)
    folder=fullfile(runRoot,sprintf('seed%d',seeds(iseed)),sprintf('akp%03d',round(100*akp)));
    d=load(fullfile(folder,'initial.mat'),'A','k','report');r=d.report;N=2001;t=(0:N-1)'*r.output_dt;raw=[];
    for phase=1:4
        file=fullfile(folder,sprintf('phase%03d',(phase-1)*90),'Results','probes.dat');fid=fopen(file);assert(fid>=0);
        while ~feof(fid),line=fgetl(fid);if startsWith(strtrim(line),'VARIABLES'),break;end,end
        values=fscanf(fid,'%f',[2,Inf]).';fclose(fid);
        if phase==1,raw=zeros(size(values,1),4);time=values(:,1);else,assert(max(abs(time-values(:,1)))<1e-10);end
        raw(:,phase)=values(:,2); %#ok<AGROW>
    end
    window=cosine_taper(t,50*r.Tp,5*r.Tp);rawTaper=raw.*window;
    ht=imag(analytic_signal(raw));htTaper=imag(analytic_signal(rawTaper));
    third=(raw(:,1)-raw(:,3)+ht(:,2)-ht(:,4))/4;
    thirdTaper=(rawTaper(:,1)-rawTaper(:,3)+htTaper(:,2)-htTaper(:,4))/4;

    omega=sqrt(r.g*d.k.*tanh(d.k*r.h));base=exp(1i*(r.probe_x*d.k-t*omega))*d.A.';linear=zeros(N,4);
    for phase=1:4,linear(:,phase)=real(base*exp(1i*(phase-1)*pi/2));end
    linearTaper=linear.*window;lh=imag(analytic_signal(linear));lht=imag(analytic_signal(linearTaper));
    linearThird=(linear(:,1)-linear(:,3)+lh(:,2)-lh(:,4))/4;
    linearThirdTaper=(linearTaper(:,1)-linearTaper(:,3)+lht(:,2)-lht(:,4))/4;

    ia=find(abs(amplitudes-akp)<1e-12,1);field=fields{(iseed-1)*numel(amplitudes)+ia};tx=t/r.Tp;
    nexttile((iseed-1)*3+1);plot(tx,third,'Color',[.65 .65 .65],'DisplayName','untapered');hold on;
    plot(tx,thirdTaper,'k-','DisplayName','5T_p taper');xline(5,'k:','HandleVisibility','off');xline(45,'k:','HandleVisibility','off');
    xlim([0 50]);grid on;ylabel('Third sector (m)');title(sprintf('Seed %d: HOS extraction',seeds(iseed)));if iseed==1,legend('Location','best');end
    nexttile((iseed-1)*3+2);plot(tx,linearThird,'Color',[.8 .35 .2],'DisplayName','untapered linear leakage');hold on;
    plot(tx,linearThirdTaper,'b-','DisplayName','tapered linear leakage');xlim([0 50]);grid on;title('Pure-linear control (truth = 0)');if iseed==1,legend('Location','best');end
    nexttile((iseed-1)*3+3);plot(field.t/r.Tp,field.HOS_third,'k-','DisplayName','tapered HOS third');hold on;
    plot(field.t/r.Tp,field.GL_eta33,'--','DisplayName','tapered-input GL8 eta_{33}');xlim([15 35]);grid on;title('Central comparison');if iseed==1,legend('Location','best');end
    if iseed==3,xlabel(nexttile(7),'t/T_p');xlabel(nexttile(8),'t/T_p');xlabel(nexttile(9),'t/T_p');end
end
sgtitle(sprintf('Random-phase third-order boundary diagnostic, nominal A k_p=%.2f, actual k_pH_s/2=%.5f',akp,r.actual_kp_Hs_over_2));
out=fullfile(runRoot,'boundary_plot');assert(~isfolder(out));mkdir(out);
exportgraphics(fig,fullfile(out,'random_eta33_boundary_effect.png'),'Resolution',180);
exportgraphics(fig,fullfile(out,'random_eta33_boundary_effect.pdf'),'ContentType','vector');close(fig);
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
