function run_hos_modal_grid_wavegroup_gl22(out)
%RUN_HOS_MODAL_GRID_WAVEGROUP_GL22 Five-probe GL8 modal-grid trial.
run='/home/lxy/green-laplace-unidirectional-time-series-runs/hos-directional-fourphase-20260926-v2';
s=jsondecode(fileread(fullfile(run,'settings.json')));
base=fullfile(run,'full-gl-comparison-20260926-v1');
d=load(fullfile(base,'comparison.mat'),'first');
initial=load(fullfile(run,'inputs','initial_fields.mat'),'C','kx','ky','om');
raw=zeros(s.expected_samples,5,4);
for phase=1:4
    values=readmatrix(fullfile(base,sprintf('phi%03d.csv',(phase-1)*90)));
    if phase==1,t=values(:,1);else,assert(isequal(t,values(:,1)));end
    raw(:,:,phase)=values(:,2:end);
end
reference=(raw(:,:,1)-raw(:,:,2)+raw(:,:,3)-raw(:,:,4))/4;
first=d.first;

models=cell(5,1);info=cell(5,1);limits=zeros(5,2);
prepare_seconds=zeros(5,1);
for probe=1:5
    clock=tic;
    [models{probe},info{probe}]=prepare_directional_joint_wavegroup( ...
        initial,first(:,probe),t,s,probe,7.5,.999,20);
    prepare_seconds(probe)=toc(clock);
    [~,peak]=max(abs(hilbert(first(:,probe))));
    limits(probe,:)=t(peak)+[-2,2]*s.Tp;
    fprintf('INPUT p%d bins=%d:%d count=%d energy=%.10f%% projection=%.6g\n', ...
        probe,info{probe}.retained_bins(1),info{probe}.retained_bins(end), ...
        info{probe}.frequency_count,100*info{probe}.retained_energy_fraction, ...
        info{probe}.projection_relative);
end
input_audits=cellfun(@(x)rmfield(x,'input'),info,'UniformOutput',false);
writejson(fullfile(out,'input_audit.json'),input_audits);

grid_sizes=[64,128,256];
predictions=cell(5,numel(grid_sizes));
rows=cell(0,15);
for probe=1:5
    model=models{probe};
    score=t>=limits(probe,1)&t<=limits(probe,2);
    previous=[];
    for level=1:numel(grid_sizes)
        grid_size=grid_sizes(level);
        options=struct('Nx',grid_size,'Ny',grid_size,'J',8);
        fprintf('START p%d grid=%dx%d\n',probe,grid_size,grid_size);
        clock=tic;
        [eta,audit]=gl_directional_time_modal_grid( ...
            model.A,model.omega,model.kx,model.ky, ...
            s.g,s.h,s.kp,t-t(1),options);
        seconds=toc(clock);
        value=real(eta);predictions{probe,level}=value;
        grid_change=NaN;
        if ~isempty(previous)
            grid_change=norm(value(score)-previous(score))/norm(value(score));
        end
        previous=value;
        hos_relative=norm(value(score)-reference(score,probe)) ...
            /norm(reference(score,probe));
        rows(end+1,:)={probe,grid_size,seconds,grid_change,hos_relative, ...
            info{probe}.retained_bins(1),info{probe}.retained_bins(end), ...
            info{probe}.frequency_count,info{probe}.retained_energy_fraction, ...
            info{probe}.projection_relative, ...
            audit.wavevector_projection_rms_relative, ...
            audit.dispersion_residual_rms_relative, ...
            audit.energy_weighted_invalid_laplace_source_fraction, ...
            audit.significant_bin_invalid_laplace_source_fraction, ...
            audit.station_input_relative}; %#ok<AGROW>
        metrics=cell2table(rows,'VariableNames',{ ...
            'probe','grid_size','GL_seconds','main_window_grid_change', ...
            'main_window_HOS_relative','first_bin','last_bin', ...
            'frequency_count','retained_energy','first_projection_relative', ...
            'wavevector_projection_rms','dispersion_residual_rms', ...
            'invalid_source_energy','worst_significant_bin_invalid', ...
            'station_input_relative'});
        writetable(metrics,fullfile(out,'metrics.csv'));
        save(fullfile(out,'progress.mat'),'predictions','metrics','t', ...
            'reference','first','limits','audit');
        disp(metrics(end,:));
    end
end

report=struct( ...
    'settings',s,'input_audits',{input_audits}, ...
    'input_preparation_seconds',prepare_seconds,'metrics',table2struct(metrics), ...
    'main_window_seconds',limits,'whole_record_samples',numel(t), ...
    'selection_rule','peak-centred contiguous observed wave-group band; >=20 bins and >=99.9 percent full-record first-harmonic energy', ...
    'acceptance_focus','final eta22 in the main wave-group window; lowest-frequency worst-bin behaviour is diagnostic only', ...
    'finite_GL_rank',8,'pair_loops',0,'time_snapshot_loops',0, ...
    'polynomial_resolvent',false,'fitted_coefficients',0);
writejson(fullfile(out,'report.json'),report);
save(fullfile(out,'comparison.mat'),'report','predictions','metrics', ...
    't','first','reference','limits');

figure_handle=figure('Visible','off','Color','w','Position',[100,100,1300,900]);
tiledlayout(3,2);
for probe=1:5
    nexttile;
    plot(t,reference(:,probe),'k-',t,predictions{probe,end},'r--');
    grid on;xlim(limits(probe,:));
    xlabel('Time (s)');ylabel('Elevation (m)');
    title(sprintf('Probe %d: bins %d--%d',probe, ...
        info{probe}.retained_bins(1),info{probe}.retained_bins(end)));
    if probe==1
        legend('HOS second sector','Joint-modal GL8','Location','best');
    end
end
exportgraphics(figure_handle,fullfile(out,'comparison.png'),'Resolution',150);
close(figure_handle);
end

function writejson(file,value)
fid=fopen(file,'w');assert(fid>=0);
fprintf(fid,'%s\n',jsonencode(value));fclose(fid);
end
