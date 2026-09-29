function plot_hos_random_harmonics_akp012(runRoot)
% HOS second/third harmonics and spectra for nominal Akp=.12.
arguments
    runRoot (1,:) char
end
akp=.12;seeds=[20260925 20260926 20260927];fig=figure('Visible','off','Color','w','Position',[40 40 1900 1100]);
tiledlayout(3,4,'Padding','compact','TileSpacing','compact');
for iseed=1:numel(seeds)
    folder=fullfile(runRoot,sprintf('seed%d',seeds(iseed)),'akp012');d=load(fullfile(folder,'initial.mat'),'report');r=d.report;raw=[];
    for phase=1:4
        file=fullfile(folder,sprintf('phase%03d',(phase-1)*90),'Results','probes.dat');fid=fopen(file);assert(fid>=0);
        while ~feof(fid),line=fgetl(fid);if startsWith(strtrim(line),'VARIABLES'),break;end,end
        values=fscanf(fid,'%f',[2,Inf]).';fclose(fid);
        if phase==1,time=values(:,1);raw=zeros(size(values,1),4);else,assert(max(abs(time-values(:,1)))<1e-10);end
        raw(:,phase)=values(:,2); %#ok<AGROW>
    end
    t=time-time(1);window=cosine_taper(t,50*r.Tp,5*r.Tp);tapered=raw.*window;ht=imag(analytic_signal(tapered));
    second=(tapered(:,1)-tapered(:,2)+tapered(:,3)-tapered(:,4))/4;
    third=(tapered(:,1)-tapered(:,3)+ht(:,2)-ht(:,4))/4;
    score=t>=15*r.Tp-1e-8&t<=35*r.Tp+1e-8;tx=t/r.Tp;N=numel(t);dt=mean(diff(t));wp=2*pi/r.Tp;
    omega=(0:floor(N/2))'*2*pi/(N*dt);frequency=omega/wp;
    F2=fft(second)/N;F3=fft(third)/N;amp2=abs(F2(1:numel(frequency)));amp3=abs(F3(1:numel(frequency)));
    amp2(2:end)=2*amp2(2:end);amp3(2:end)=2*amp3(2:end);

    nexttile((iseed-1)*4+1);plot(tx(score),1e3*second(score),'k-');grid on;xlim([15 35]);
    ylabel('Elevation (mm)');title(sprintf('Seed %d: second harmonic',seeds(iseed)));
    nexttile((iseed-1)*4+2);plot(tx(score),1e3*third(score),'k-');grid on;xlim([15 35]);
    ylabel('Elevation (mm)');title('Third harmonic');
    nexttile((iseed-1)*4+3);semilogy(frequency(2:end),amp2(2:end),'b-');grid on;xlim([0 6]);
    ylabel('One-sided amplitude (m)');title('Second-harmonic spectrum');
    nexttile((iseed-1)*4+4);semilogy(frequency(2:end),amp3(2:end),'r-');grid on;xlim([0 6]);
    ylabel('One-sided amplitude (m)');title('Third-harmonic spectrum');
    if iseed==3
        xlabel(nexttile(9),'t/T_p');xlabel(nexttile(10),'t/T_p');
        xlabel(nexttile(11),'\omega/\omega_p');xlabel(nexttile(12),'\omega/\omega_p');
    end
end
sgtitle(sprintf('Unidirectional random-phase HOS harmonics, nominal A k_p=%.2f, actual k_pH_s/2=%.5f',akp,r.actual_kp_Hs_over_2));
out=fullfile(runRoot,'harmonics_akp012');assert(~isfolder(out));mkdir(out);
exportgraphics(fig,fullfile(out,'hos_second_third_harmonics_akp012.png'),'Resolution',180);
exportgraphics(fig,fullfile(out,'hos_second_third_harmonics_akp012.pdf'),'ContentType','vector');close(fig);
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
