function analyze_hos_random_double_domain(sourceRoot,targetRoot)
% Compare original and exactly repeated doubled-domain HOS records.
arguments
    sourceRoot (1,:) char
    targetRoot (1,:) char
end
source=fullfile(sourceRoot,'seed20260927','akp012');d=load(fullfile(source,'initial.mat'),'report');r=d.report;
original=zeros(2001,4);double1=original;double2=original;time=[];rows=cell(4,7);
for phase=1:4
    a=read_probes(fullfile(source,sprintf('phase%03d',(phase-1)*90),'Results','probes.dat'),2);
    b=read_probes(fullfile(targetRoot,sprintf('phase%03d',(phase-1)*90),'Results','probes.dat'),3);
    assert(size(a,2)==2&&size(b,2)==3&&size(a,1)==size(b,1));
    if phase==1,time=a(:,1);else,assert(max(abs(time-a(:,1)))<1e-10);end
    assert(max(abs(time-b(:,1)))<1e-10);original(:,phase)=a(:,2);double1(:,phase)=b(:,2);double2(:,phase)=b(:,3);
    central=time>=15*r.Tp-1e-8&time<=35*r.Tp+1e-8;
    rows(phase,:)={(phase-1)*90,norm(double1(:,phase)-original(:,phase))/norm(original(:,phase)), ...
        norm(double2(:,phase)-double1(:,phase))/norm(double1(:,phase)), ...
        norm(double1(central,phase)-original(central,phase))/norm(original(central,phase)), ...
        norm(double2(central,phase)-double1(central,phase))/norm(double1(central,phase)), ...
        max(abs(double1(:,phase)-original(:,phase))),max(abs(double2(:,phase)-double1(:,phase)))};
end
t=time-time(1);window=cosine_taper(t,50*r.Tp,5*r.Tp);original=original.*window;double1=double1.*window;double2=double2.*window;
[second0,third0]=harmonics(original);[second1,third1]=harmonics(double1);[second2,third2]=harmonics(double2);
central=t>=15*r.Tp-1e-8&t<=35*r.Tp+1e-8;
harmonic=table( ...
    norm(second1(central)-second0(central))/norm(second0(central)), ...
    norm(second2(central)-second1(central))/norm(second1(central)), ...
    norm(third1(central)-third0(central))/norm(third0(central)), ...
    norm(third2(central)-third1(central))/norm(third1(central)), ...
    'VariableNames',{'second_double_vs_original','second_repeat_probe_difference','third_double_vs_original','third_repeat_probe_difference'});
metrics=cell2table(rows,'VariableNames',{'global_phase_deg','full_double_vs_original','full_repeat_probe_difference','central_double_vs_original','central_repeat_probe_difference','max_abs_double_vs_original_m','max_abs_repeat_probe_difference_m'});
writetable(metrics,fullfile(targetRoot,'double_domain_metrics.csv'));writetable(harmonic,fullfile(targetRoot,'double_domain_harmonics.csv'));
fig=figure('Visible','off','Color','w','Position',[100 100 1300 620]);tiledlayout(2,2,'Padding','compact');tx=t/r.Tp;
nexttile;plot(tx(central),1e3*second0(central),'k-',tx(central),1e3*second1(central),'--');grid on;title('Second harmonic: original and doubled domain');ylabel('Elevation (mm)');legend('68\lambda_p','136\lambda_p');
nexttile;plot(tx(central),1e3*third0(central),'k-',tx(central),1e3*third1(central),'--');grid on;title('Third harmonic: original and doubled domain');ylabel('Elevation (mm)');legend('68\lambda_p','136\lambda_p');
nexttile;plot(tx(central),1e6*(second1(central)-second0(central)));grid on;xlabel('t/T_p');ylabel('Difference (\mum)');title('Second-harmonic difference');
nexttile;plot(tx(central),1e6*(third1(central)-third0(central)));grid on;xlabel('t/T_p');ylabel('Difference (\mum)');title('Third-harmonic difference');
exportgraphics(fig,fullfile(targetRoot,'double_domain_comparison.png'),'Resolution',180);exportgraphics(fig,fullfile(targetRoot,'double_domain_comparison.pdf'),'ContentType','vector');close(fig);
disp(metrics);disp(harmonic);
end

function values=read_probes(file,ncols)
fid=fopen(file);assert(fid>=0);while ~feof(fid),line=fgetl(fid);if startsWith(strtrim(line),'VARIABLES'),break;end,end
values=fscanf(fid,'%f',[ncols,Inf]).';fclose(fid);
end

function [second,third]=harmonics(raw)
ht=imag(analytic_signal(raw));second=(raw(:,1)-raw(:,2)+raw(:,3)-raw(:,4))/4;
third=(raw(:,1)-raw(:,3)+ht(:,2)-ht(:,4))/4;
end

function window=cosine_taper(t,duration,width)
window=ones(size(t));left=t<width;window(left)=.5*(1-cos(pi*t(left)/width));right=t>duration-width;
window(right)=.5*(1-cos(pi*(duration-t(right))/width));window=max(0,min(1,window));
end

function z=analytic_signal(y)
N=size(y,1);mask=zeros(N,1);mask(1)=1;if mod(N,2),mask(2:(N+1)/2)=2;else,mask(2:N/2)=2;mask(N/2+1)=1;end
z=ifft(fft(y).*mask);
end
