function run_jonswap006_native_eta22(out)
run='/home/lxy/green-laplace-unidirectional-time-series-runs/hos-jonswap-kphs006-20tp-20260927-v1/medium';
s=jsondecode(fileread(fullfile(run,'settings.json')));initial=load(fullfile(run,'inputs','initial_fields.mat'),'C','kx','ky','om');
phases=[0,90,180,270];v=readprobe(fullfile(run,'cases','phi000','Results','probes.dat'));t=v(:,1);raw=zeros(numel(t),5,4);raw(:,:,1)=v(:,2:end);
for p=2:4,v=readprobe(fullfile(run,'cases',sprintf('phi%03d',phases(p)),'Results','probes.dat'));raw(:,:,p)=v(:,2:end);end
H=zeros(size(raw));for p=1:4,H(:,:,p)=imag(hilbert(raw(:,:,p)));end
first=(raw(:,:,1)-raw(:,:,3)-H(:,:,2)+H(:,:,4))/4;reference=(raw(:,:,1)-raw(:,:,2)+raw(:,:,3)-raw(:,:,4))/4;score=t>=s.scoring_window_s(1)&t<=s.scoring_window_s(2);
prediction=zeros(size(first));rows=zeros(5,9);total=tic;
for probe=1:5,clock=tic;[eta,audit]=gl_native_joint_eta22(initial,first(:,probe),t,s,probe,8);seconds=toc(clock);prediction(:,probe)=real(eta);rows(probe,:)=[probe,seconds,audit.first_projection_relative,audit.retained_first_energy,audit.input_bins(1),audit.input_bins(end),audit.energy_weighted_condition,norm(prediction(score,probe)-reference(score,probe))/norm(reference(score,probe)),norm(prediction(score,probe))/norm(reference(score,probe))];end
metrics=array2table(rows,'VariableNames',{'probe','seconds','first_projection_relative','retained_first_energy','first_bin','last_bin','energy_weighted_condition','HOS_relative','norm_ratio'});writetable(metrics,fullfile(out,'metrics.csv'));
report=struct('source',run,'method','native HOS k grid plus observed temporal-bin complex correction','total_seconds',toc(total),'median_seconds',median(metrics.seconds),'metrics',table2struct(metrics));fid=fopen(fullfile(out,'report.json'),'w');fprintf(fid,'%s\n',jsonencode(report));fclose(fid);save(fullfile(out,'comparison.mat'),'report','metrics','t','first','reference','prediction');disp(metrics);disp(report)
end
function v=readprobe(file),lines=readlines(file);header=find(startsWith(strtrim(lines),'VARIABLES'),1,'last');v=readmatrix(file,'NumHeaderLines',header);v=v(all(isfinite(v),2),:);end
