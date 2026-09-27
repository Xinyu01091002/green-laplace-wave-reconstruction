function preview_gpu_phase_records(run)
% Same phase-sector definitions as the unchanged full GL postprocessor.
out=fullfile(run,'gl-time-comparison-v1');s=jsondecode(fileread(fullfile(run,'settings.json')));
raw=zeros(s.expected_samples,5,4);
for j=1:4
 a=readmatrix(fullfile(out,sprintf('phi%03d.csv',(j-1)*90)));
 assert(isequal(size(a),[s.expected_samples,6])&&all(isfinite(a),'all'));
 if j==1,t=a(:,1);else,assert(max(abs(a(:,1)-t))<1e-9);end
 raw(:,:,j)=a(:,2:end);
end
HH=imag(hilbert(raw));first=(raw(:,:,1)-raw(:,:,3)-HH(:,:,2)+HH(:,:,4))/4;
second=(raw(:,:,1)-raw(:,:,2)+raw(:,:,3)-raw(:,:,4))/4;
amplitudes=max(abs(second),[],1);ratios=amplitudes/amplitudes(1);
report=struct('samples',numel(t),'end_time_s',t(end),'second_sector_peak_m',amplitudes,'off_axis_amplitude_ratios',ratios,'eligible',ratios>=1/3,'note','Four-phase sectors, not exact perturbation orders. No alignment or amplitude correction.');
f=fopen(fullfile(out,'phase_preview.json'),'w');fprintf(f,'%s',jsonencode(report,PrettyPrint=true));fclose(f);
fig=figure('Visible','off','Position',[40 40 1400 900]);tiledlayout(3,1,'TileSpacing','compact');
nexttile;plot(t,raw(:,1,1),'Color',[.65 .65 .65]);hold on;plot(t,first(:,1),'b');grid on;xlim([0 t(end)]);ylabel('Elevation (m)');title('Center probe: complete GPU HOS record');legend('Phase 0 total','Separated first harmonic','Box','off');
nexttile;plot(t,second,'LineWidth',.8);grid on;xlim([0 t(end)]);ylabel('Elevation (m)');title('Second-harmonic phase sector: five probes');
nexttile;plot(t,second(:,1),'k','LineWidth',1);grid on;xlim([10 15]*s.Tp);xlabel('Time (s)');ylabel('Elevation (m)');title('Center probe: fixed 10--15 Tp detail (no GL prediction yet)');
exportgraphics(fig,fullfile(out,'phase_preview.png'),'Resolution',160);close(fig);
end
