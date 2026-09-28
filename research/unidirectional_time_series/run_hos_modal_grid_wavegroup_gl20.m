function run_hos_modal_grid_wavegroup_gl20(out)
%RUN_HOS_MODAL_GRID_WAVEGROUP_GL20 Focused-wave-group eta20 reconstruction.
source='/home/lxy/green-laplace-unidirectional-time-series-runs/hos-directional-fourphase-20260926-v2';
s=jsondecode(fileread(fullfile(source,'settings.json')));
base=fullfile(source,'full-gl-comparison-20260926-v1');
d=load(fullfile(base,'comparison.mat'),'first');first=d.first;
initial=load(fullfile(source,'inputs','initial_fields.mat'),'C','kx','ky','om');
raw=zeros(s.expected_samples,5,4);
for phase=1:4
    values=readmatrix(fullfile(base,sprintf('phi%03d.csv',(phase-1)*90)));
    if phase==1,t=values(:,1);else,assert(isequal(t,values(:,1)));end
    raw(:,:,phase)=values(:,2:end);
end
observed_raw=mean(raw,3);N=numel(t);dt=mean(diff(t));
native_frequency=2*pi*[0:floor(N/2),-floor(N/2):-1]'/(N*dt);
omega_p=sqrt(s.g*s.kp*tanh(s.kp*s.h));
subharmonic_mask=abs(native_frequency)>0 & abs(native_frequency)<.5*omega_p;
reference=real(ifft(fft(observed_raw).*subharmonic_mask));

models20=cell(5,1);info20=cell(5,1);limits=zeros(5,2);
for probe=1:5
    [models20{probe},info20{probe}]=prepare_directional_joint_wavegroup( ...
        initial,first(:,probe),t,s,probe,7.5,.999,20);
    [~,peak]=max(abs(hilbert(first(:,probe))));
    limits(probe,:)=t(peak)+[-2,2]*s.Tp;
end

grid_sizes=[128,256,512];grid_prediction=cell(5,3);grid_rows=cell(0,15);
for probe=1:5
    model=models20{probe};score=window_mask(t,limits(probe,:));
    [qx_limit,qy_limit]=fixed_limits(model,s.h,2,256);
    temporary=cell(3,15);
    for grid_index=1:3
        grid_size=grid_sizes(grid_index);fprintf('GRID p%d size=%d\n',probe,grid_size);
        options=struct('Nx',grid_size,'Ny',grid_size,'J',16, ...
            'QxLimit',qx_limit,'QyLimit',qy_limit);
        clock=tic;[eta,audit]=gl_directional_time_modal_grid_eta20( ...
            model.A,model.omega,model.kx,model.ky,s.g,s.h,t-t(1),options);
        seconds=toc(clock);value=project_subharmonic(eta,subharmonic_mask);
        grid_prediction{probe,grid_index}=value;
        temporary(grid_index,:)={probe,grid_size,seconds,NaN, ...
            relative_error(value,reference(:,probe),score), ...
            relative_error(value,reference(:,probe),true(N,1)), ...
            norm(value(score))/norm(reference(score,probe)), ...
            audit.wavevector_projection_rms_relative, ...
            audit.dispersion_residual_rms_relative, ...
            audit.strict_zero_source_energy_fraction, ...
            audit.nonzero_invalid_laplace_source_energy_fraction, ...
            audit.delta_q_rms_pair,audit.work_frequency_count, ...
            mean(observed_raw(:,probe)),max(abs(value(score)))};
    end
    for grid_index=1:3
        temporary{grid_index,4}=relative_change( ...
            grid_prediction{probe,grid_index},grid_prediction{probe,3},score);
    end
    grid_rows=[grid_rows;temporary]; %#ok<AGROW>
    grid_metrics=cell2table(grid_rows,'VariableNames',{ ...
        'probe','grid_size','GL_seconds','main_window_change_to_512', ...
        'main_window_HOS_relative','full_record_HOS_relative','norm_ratio', ...
        'wavevector_projection_rms','dispersion_residual_rms', ...
        'strict_zero_source_energy','nonzero_invalid_source_energy', ...
        'delta_q_rms_pair','work_frequency_count', ...
        'excluded_observed_mean_m','main_window_max_abs_m'});
    writetable(grid_metrics,fullfile(out,'grid_metrics.csv'));
end

ranks=[6,12,16];rank_prediction=cell(5,3);rank_rows=cell(0,8);
for probe=1:5
    model=models20{probe};score=window_mask(t,limits(probe,:));
    [qx_limit,qy_limit]=fixed_limits(model,s.h,2,256);
    temporary=cell(3,8);
    for rank_index=1:3
        rank=ranks(rank_index);fprintf('RANK p%d J=%d\n',probe,rank);
        options=struct('Nx',256,'Ny',256,'J',rank, ...
            'QxLimit',qx_limit,'QyLimit',qy_limit);
        clock=tic;eta=gl_directional_time_modal_grid_eta20( ...
            model.A,model.omega,model.kx,model.ky,s.g,s.h,t-t(1),options);
        seconds=toc(clock);value=project_subharmonic(eta,subharmonic_mask);
        rank_prediction{probe,rank_index}=value;
        temporary(rank_index,:)={probe,rank,seconds,NaN, ...
            relative_error(value,reference(:,probe),score), ...
            relative_error(value,reference(:,probe),true(N,1)), ...
            norm(value(score))/norm(reference(score,probe)),max(abs(value(score)))};
    end
    for rank_index=1:3
        temporary{rank_index,4}=relative_change( ...
            rank_prediction{probe,rank_index},rank_prediction{probe,3},score);
    end
    rank_rows=[rank_rows;temporary]; %#ok<AGROW>
    rank_metrics=cell2table(rank_rows,'VariableNames',{ ...
        'probe','GL_rank','GL_seconds','main_window_change_to_GL16', ...
        'main_window_HOS_relative','full_record_HOS_relative','norm_ratio', ...
        'main_window_max_abs_m'});
    writetable(rank_metrics,fullfile(out,'rank_metrics.csv'));
end

probes=[1,4];counts=[20,36];band_prediction=cell(2,2);band_rows=cell(0,10);
for probe_index=1:2
    probe=probes(probe_index);score=window_mask(t,limits(probe,:));
    model=cell(1,2);input_info=cell(1,2);
    for count_index=1:2
        [model{count_index},input_info{count_index}]= ...
            prepare_directional_joint_wavegroup( ...
            initial,first(:,probe),t,s,probe,7.5,.999,counts(count_index));
    end
    [qx_limit,qy_limit]=fixed_limits(model{2},s.h,2,256);
    temporary=cell(2,10);
    for count_index=1:2
        options=struct('Nx',256,'Ny',256,'J',16, ...
            'QxLimit',qx_limit,'QyLimit',qy_limit);
        fprintf('BAND p%d count=%d\n',probe,counts(count_index));
        clock=tic;eta=gl_directional_time_modal_grid_eta20( ...
            model{count_index}.A,model{count_index}.omega, ...
            model{count_index}.kx,model{count_index}.ky, ...
            s.g,s.h,t-t(1),options);
        seconds=toc(clock);value=project_subharmonic(eta,subharmonic_mask);
        band_prediction{probe_index,count_index}=value;
        temporary(count_index,:)={probe,counts(count_index), ...
            input_info{count_index}.retained_bins(1), ...
            input_info{count_index}.retained_bins(end),seconds,NaN, ...
            relative_error(value,reference(:,probe),score), ...
            norm(value(score))/norm(reference(score,probe)), ...
            input_info{count_index}.retained_energy_fraction, ...
            2*counts(count_index)-1};
    end
    for count_index=1:2
        temporary{count_index,6}=relative_change( ...
            band_prediction{probe_index,count_index}, ...
            band_prediction{probe_index,2},score);
    end
    band_rows=[band_rows;temporary]; %#ok<AGROW>
    band_metrics=cell2table(band_rows,'VariableNames',{ ...
        'probe','minimum_count','first_bin','last_bin','GL_seconds', ...
        'main_window_change_to_36','main_window_HOS_relative','norm_ratio', ...
        'retained_first_energy','work_frequency_count'});
    writetable(band_metrics,fullfile(out,'band_metrics.csv'));
end

report=struct('sector','eta20 nonzero spatial difference', ...
    'reference','four-phase HOS average projected identically to nonzero abs(omega)<0.5 omega_p', ...
    'strict_zero_spatial_mode','excluded','temporal_DC','excluded by comparison mask', ...
    'output_cutoff_rad_s',.5*omega_p,'grid_sizes',grid_sizes,'GL_ranks',ranks, ...
    'band_counts',counts,'main_window_seconds',limits, ...
    'grid_metrics',table2struct(grid_metrics), ...
    'rank_metrics',table2struct(rank_metrics), ...
    'band_metrics',table2struct(band_metrics), ...
    'pair_loops',0,'fitted_coefficients',0);
writejson(fullfile(out,'report.json'),report);
save(fullfile(out,'comparison.mat'),'report','grid_metrics','rank_metrics', ...
    'band_metrics','grid_prediction','rank_prediction','band_prediction', ...
    't','first','reference','observed_raw','limits','subharmonic_mask');

figure_handle=figure('Visible','off','Color','w','Position',[100,100,1300,900]);
tiledlayout(3,2);
for probe=1:5
    nexttile;plot(t,reference(:,probe),'k-',t,grid_prediction{probe,3},'r--');
    grid on;xlim(limits(probe,:));xlabel('Time (s)');ylabel('Elevation (m)');
    title(sprintf('Probe %d: eta20, GL16, 512^2',probe));
    if probe==1,legend('HOS even-sector subharmonic','Joint-modal GL16 eta20','Location','best');end
end
exportgraphics(figure_handle,fullfile(out,'comparison.png'),'Resolution',160);
close(figure_handle);
end

function value=project_subharmonic(eta,mask)
value=real(ifft(fft(eta).*mask));
end

function [qx_limit,qy_limit]=fixed_limits(model,h,order,reference_grid)
qx=h*model.kx;qy=h*model.ky;
qx_limit=(order+.05)*max(abs(qx))/(1-2/reference_grid);
qy_limit=(order+.05)*max(abs(qy))/(1-2/reference_grid);
end

function mask=window_mask(t,limits)
mask=t>=limits(1)&t<=limits(2);
end

function value=relative_error(candidate,reference,mask)
value=norm(candidate(mask)-reference(mask))/norm(reference(mask));
end

function value=relative_change(candidate,reference,mask)
value=norm(candidate(mask)-reference(mask))/norm(reference(mask));
end

function writejson(file,value)
fid=fopen(file,'w');assert(fid>=0);fprintf(fid,'%s\n',jsonencode(value));fclose(fid);
end
