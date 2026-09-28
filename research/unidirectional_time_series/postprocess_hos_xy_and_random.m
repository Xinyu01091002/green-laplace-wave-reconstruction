function postprocess_hos_xy_and_random(focused_out,random_out)
%POSTPROCESS_HOS_XY_AND_RANDOM Current GL workflow on focused XY and random HOS.
focused_root='/home/lxy/green-laplace-unidirectional-time-series-runs/hos-directional-high-xy-probes-20260928-v1';
focused_initial='/home/lxy/green-laplace-unidirectional-time-series-runs/hos-directional-fourphase-20260926-v2/inputs/initial_fields.mat';
random_root='/home/lxy/green-laplace-unidirectional-time-series-runs/hos-random-fourphase-20260926-v1/akp012';
process_family(focused_root,focused_initial,focused_out,"focused_xy");
process_family(random_root,fullfile(random_root,'inputs','initial_fields.mat'),random_out,"random_phase");
end

function process_family(run_root,initial_file,out,family)
s=jsondecode(fileread(fullfile(run_root,'settings.json')));
initial=load(initial_file,'C','kx','ky','om');
phases=[0,90,180,270];
value=read_probe_file(fullfile(run_root,'cases','phi000','Results','probes.dat'));
t=value(:,1);raw=zeros(numel(t),size(value,2)-1,4);raw(:,:,1)=value(:,2:end);
for phase_index=2:4
    file=fullfile(run_root,'cases',sprintf('phi%03d',phases(phase_index)),'Results','probes.dat');
    value=read_probe_file(file);assert(isequal(t,value(:,1)));
    raw(:,:,phase_index)=value(:,2:end);
end
H=zeros(size(raw));
for phase_index=1:4,H(:,:,phase_index)=imag(hilbert(raw(:,:,phase_index)));end
first=(raw(:,:,1)-raw(:,:,3)-H(:,:,2)+H(:,:,4))/4;
eta22_ref=(raw(:,:,1)-raw(:,:,2)+raw(:,:,3)-raw(:,:,4))/4;
eta33_ref=(raw(:,:,1)-raw(:,:,3)+H(:,:,2)-H(:,:,4))/4;
eta20_raw=mean(raw,3);N=numel(t);dt=mean(diff(t));
native_frequency=2*pi*[0:floor(N/2),-floor(N/2):-1]'/(N*dt);
omega_p=sqrt(s.g*s.kp*tanh(s.kp*s.h));
low_mask=abs(native_frequency)>0&abs(native_frequency)<.5*omega_p;
eta20_ref=real(ifft(fft(eta20_raw).*low_mask));

if family=="focused_xy"
    rows=(5:12)';centres=s.centre_for_offaxis_rows(:);
    first_max=max(abs(first),[],1);total_max=max(abs(raw(:,:,1)),[],1);
    first_ratio=first_max(rows).'./first_max(centres).';
    total_ratio=total_max(rows).'./total_max(centres).';
    assert(all(first_ratio>=1/3 | total_ratio>=1/3));
    limits=zeros(numel(rows),2);
    for index=1:numel(rows)
        [~,peak]=max(abs(first(:,rows(index))));limits(index,:)=t(peak)+[-2,2]*s.Tp;
    end
else
    rows=(1:size(first,2))';centres=rows;first_ratio=ones(size(rows));total_ratio=ones(size(rows));
    limits=repmat(s.scoring_window_s(:).',numel(rows),1);
end

pred20=zeros(N,numel(rows));pred22=pred20;pred33=pred20;
metrics_rows=zeros(numel(rows),24);input_info=cell(numel(rows),1);
for index=1:numel(rows)
    probe=rows(index);[model,input_info{index}]=prepare_directional_joint_wavegroup( ...
        initial,first(:,probe),t,s,probe,7.5,.999,20);
    clock=tic;eta22=gl_directional_time_modal_grid( ...
        model.A,model.omega,model.kx,model.ky,s.g,s.h,s.kp,t-t(1), ...
        struct('Nx',256,'Ny',256,'J',8,'Baseband',true));seconds22=toc(clock);
    clock=tic;eta33=gl_directional_time_modal_grid_eta33( ...
        model.A,model.omega,model.kx,model.ky,s.g,s.h,s.kp,t-t(1), ...
        struct('Nx',256,'Ny',256,'J',10,'Baseband',true));seconds33=toc(clock);
    clock=tic;eta20=gl_directional_time_modal_grid_eta20( ...
        model.A,model.omega,model.kx,model.ky,s.g,s.h,t-t(1), ...
        struct('Nx',512,'Ny',512,'J',16));seconds20=toc(clock);
    pred20(:,index)=real(ifft(fft(eta20).*low_mask));
    pred22(:,index)=real(eta22);pred33(:,index)=real(eta33);
    score=t>=limits(index,1)&t<=limits(index,2);
    xy=s.probes_xy_m(probe,:);
    if family=="random_phase"
        offset=[NaN,NaN];
    else
        offset=(xy-s.nominal_focus_xy(:).')/s.lambda_p;
    end
    metrics_rows(index,:)=[probe,centres(index),xy,offset,first_ratio(index),total_ratio(index), ...
        input_info{index}.retained_bins(1),input_info{index}.retained_bins(end), ...
        input_info{index}.frequency_count,input_info{index}.retained_energy_fraction, ...
        input_info{index}.projection_relative, ...
        relative(pred20(:,index),eta20_ref(:,probe),score), ...
        relative(pred22(:,index),eta22_ref(:,probe),score), ...
        relative(pred33(:,index),eta33_ref(:,probe),score), ...
        norm(pred20(score,index))/norm(eta20_ref(score,probe)), ...
        norm(pred22(score,index))/norm(eta22_ref(score,probe)), ...
        norm(pred33(score,index))/norm(eta33_ref(score,probe)), ...
        seconds20,seconds22,seconds33,limits(index,:)];
end
names={'probe','same_x_centre_probe','x_m','y_m','dx_lambda','dy_lambda', ...
    'first_max_ratio','total_phi000_max_ratio','first_bin','last_bin', ...
    'frequency_count','retained_first_energy','first_projection_relative', ...
    'eta20_HOS_relative','eta22_HOS_relative','eta33_HOS_relative', ...
    'eta20_norm_ratio','eta22_norm_ratio','eta33_norm_ratio', ...
    'eta20_seconds','eta22_seconds','eta33_seconds','window_start_s','window_end_s'};
metrics=array2table(metrics_rows,'VariableNames',names);
writetable(metrics,fullfile(out,'metrics.csv'));
audits=cellfun(@(x)rmfield(x,'input'),input_info,'UniformOutput',false);
report=struct('family',family,'source_run',run_root,'Akp',s.Akp, ...
    'random_seed',field_or(s,'seed',NaN),'linear_Hs_m',field_or(s,'linear_Hs_4sigma_m',NaN), ...
    'kp_Hs_over_2',field_or(s,'kp_Hs_over_2',NaN), ...
    'eta20_scope','nonzero temporal frequencies below 0.5 omega_p; strict spatial Q=0 excluded', ...
    'eta22_scope','positive pure-sum phase sector','eta33_scope','Hilbert third phase sector', ...
    'no_fitting',true,'pair_loops',0,'triple_loops',0, ...
    'amplitude_gate_pass',all(first_ratio>=1/3 | total_ratio>=1/3), ...
    'input_audits',{audits},'metrics',table2struct(metrics));
writejson(fullfile(out,'report.json'),report);
save(fullfile(out,'comparison.mat'),'report','metrics','t','first','eta20_ref','eta22_ref','eta33_ref', ...
    'pred20','pred22','pred33','raw','limits');
plot_component(t,rows,limits,eta20_ref,pred20,'eta20',out);
plot_component(t,rows,limits,eta22_ref,pred22,'eta22',out);
plot_component(t,rows,limits,eta33_ref,pred33,'eta33',out);
disp(metrics)
end

function value=read_probe_file(file)
lines=readlines(file);header=find(startsWith(strtrim(lines),'VARIABLES'),1,'last');
assert(~isempty(header));value=readmatrix(file,'NumHeaderLines',header);
value=value(all(isfinite(value),2),:);
end

function plot_component(t,rows,limits,reference,candidate,name,out)
count=numel(rows);figure_handle=figure('Visible','off','Color','w','Position',[100,100,1300,900]);
tiledlayout(ceil(count/2),2);
for index=1:count
    nexttile;plot(t,reference(:,rows(index)),'k-',t,candidate(:,index),'r--');grid on;
    xlim(limits(index,:));xlabel('Time (s)');ylabel('Elevation (m)');title(sprintf('Probe %d %s',rows(index),name));
    if index==1,legend('HOS phase sector','Joint-modal GL','Location','best');end
end
exportgraphics(figure_handle,fullfile(out,name+".png"),'Resolution',150);close(figure_handle);
end

function value=relative(candidate,reference,mask)
value=norm(candidate(mask)-reference(mask))/norm(reference(mask));
end

function value=field_or(s,name,default)
if isfield(s,name),value=s.(name);else,value=default;end
end

function writejson(file,value)
fid=fopen(file,'w');assert(fid>=0);fprintf(fid,'%s\n',jsonencode(value));fclose(fid);
end
