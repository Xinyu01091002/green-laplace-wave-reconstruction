function diagnose_hos_wavegroup_third_spectrum(runRoot)
% Quantify the low-frequency peak in the focused HOS third harmonic.
arguments
    runRoot (1,:) char
end
folder=fullfile(runRoot,'akp012');d=load(fullfile(folder,'initial.mat'),'A','k','r');r=d.r;raw=[];
for phase=1:4
    file=fullfile(folder,sprintf('phase%03d',(phase-1)*90),'Results','probes.dat');fid=fopen(file);assert(fid>=0);
    while ~feof(fid),line=fgetl(fid);if startsWith(strtrim(line),'VARIABLES'),break;end,end
    values=fscanf(fid,'%f',[2,Inf]).';fclose(fid);
    if phase==1,time=values(:,1);raw=zeros(size(values,1),4);else,assert(max(abs(time-values(:,1)))<1e-10);end
    raw(:,phase)=values(:,2); %#ok<AGROW>
end
t=time-time(1);window=cosine_taper(t,50*r.Tp,5*r.Tp);tapered=raw.*window;ht=imag(analytic_signal(tapered));
third=(tapered(:,1)-tapered(:,3)+ht(:,2)-ht(:,4))/4;
N=numel(t);dt=mean(diff(t));wp=2*pi/r.Tp;omega=(0:floor(N/2))'*2*pi/(N*dt);ratio=omega/wp;
F=fft(third)/N;positive=F(1:numel(ratio));energy=abs(positive).^2;energy(2:end)=2*energy(2:end);
parentOmega=sqrt(r.g*d.k.*tanh(d.k*r.h));support=[3*min(parentOmega)/wp,3*max(parentOmega)/wp];
low=ratio>0&ratio<2;main=ratio>=2&ratio<=5;outside=ratio>0&ratio<support(1);
[lowAmplitude,lowIndex]=max(abs(positive(low)));lowFrequencies=ratio(low);lowPeak=lowFrequencies(lowIndex);
[mainAmplitude,mainIndex]=max(abs(positive(main)));mainFrequencies=ratio(main);mainPeak=mainFrequencies(mainIndex);

base=exp(1i*(r.probe_x*d.k-t*parentOmega))*d.A.';linear=zeros(N,4);
for phase=1:4,linear(:,phase)=real(base*exp(1i*(phase-1)*pi/2));end
linear=linear.*window;lh=imag(analytic_signal(linear));linearThird=(linear(:,1)-linear(:,3)+lh(:,2)-lh(:,4))/4;
LF=fft(linearThird)/N;linearPositive=LF(1:numel(ratio));linearLow=max(abs(linearPositive(low)));
report=struct('parent_frequency_range_over_wp',[min(parentOmega)/wp,max(parentOmega)/wp], ...
    'positive_triple_sum_support_over_wp',support,'low_peak_over_wp',lowPeak,'main_peak_over_wp',mainPeak, ...
    'low_to_main_peak_amplitude_ratio',lowAmplitude/mainAmplitude, ...
    'energy_below_2wp_fraction',sum(energy(low))/sum(energy(2:end)), ...
    'energy_below_triple_sum_support_fraction',sum(energy(outside))/sum(energy(2:end)), ...
    'linear_control_low_peak_to_HOS_low_peak',linearLow/lowAmplitude,'taper_width_Tp',5);
out=fullfile(runRoot,'wavegroup_harmonics_akp012');fid=fopen(fullfile(out,'third_spectrum_diagnostic.json'),'w');fprintf(fid,'%s\n',jsonencode(report));fclose(fid);disp(report);
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
