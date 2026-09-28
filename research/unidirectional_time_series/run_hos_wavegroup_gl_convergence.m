function run_hos_wavegroup_gl_convergence(out)
%RUN_HOS_WAVEGROUP_GL_CONVERGENCE Band, grid and rank checks at two probes.
source='/home/lxy/green-laplace-unidirectional-time-series-runs/hos-directional-fourphase-20260926-v2';
third_root='/home/lxy/green-laplace-unidirectional-time-series-runs/hos-third-extraction-focused-high-20260928-v2';
s=jsondecode(fileread(fullfile(source,'settings.json')));
base=fullfile(source,'full-gl-comparison-20260926-v1');
d=load(fullfile(third_root,'third_sector.mat'),'first','third','t');
initial=load(fullfile(source,'inputs','initial_fields.mat'),'C','kx','ky','om');
t=d.t;first=d.first;reference3=d.third;
raw=zeros(s.expected_samples,5,4);
for phase=1:4
    values=readmatrix(fullfile(base,sprintf('phi%03d.csv',(phase-1)*90)));
    assert(isequal(t,values(:,1)));raw(:,:,phase)=values(:,2:end);
end
reference2=(raw(:,:,1)-raw(:,:,2)+raw(:,:,3)-raw(:,:,4))/4;

probes=[1,4];counts=[20,24,28,36];models=cell(2,4);info=cell(2,4);
limits=zeros(2,2);
for probe_index=1:2
    probe=probes(probe_index);
    [~,peak]=max(abs(hilbert(first(:,probe))));
    limits(probe_index,:)=t(peak)+[-2,2]*s.Tp;
    for count_index=1:4
        [models{probe_index,count_index},info{probe_index,count_index}]= ...
            prepare_directional_joint_wavegroup( ...
            initial,first(:,probe),t,s,probe,7.5,.999,counts(count_index));
    end
end

band2=cell(2,4);band3=cell(2,4);band_rows=cell(0,16);
for probe_index=1:2
    probe=probes(probe_index);score=window_mask(t,limits(probe_index,:));
    broad=models{probe_index,4};
    [qx2,qy2]=fixed_limits(broad,s.h,2,256);
    [qx3,qy3]=fixed_limits(broad,s.h,3,256);
    temporary=cell(4,16);
    for count_index=1:4
        model=models{probe_index,count_index};input_info=info{probe_index,count_index};
        options2=struct('Nx',256,'Ny',256,'J',8, ...
            'QxLimit',qx2,'QyLimit',qy2,'Baseband',true);
        options3=struct('Nx',256,'Ny',256,'J',8, ...
            'QxLimit',qx3,'QyLimit',qy3,'Baseband',true);
        fprintf('BAND p%d count=%d bins=%d:%d\n',probe,counts(count_index), ...
            input_info.retained_bins(1),input_info.retained_bins(end));
        clock=tic;[eta2,audit2]=gl_directional_time_modal_grid( ...
            model.A,model.omega,model.kx,model.ky,s.g,s.h,s.kp,t-t(1),options2);
        seconds2=toc(clock);
        clock=tic;[eta3,audit3]=gl_directional_time_modal_grid_eta33( ...
            model.A,model.omega,model.kx,model.ky,s.g,s.h,s.kp,t-t(1),options3);
        seconds3=toc(clock);
        band2{probe_index,count_index}=real(eta2);
        band3{probe_index,count_index}=real(eta3);
        temporary(count_index,:)={probe,counts(count_index), ...
            input_info.retained_bins(1),input_info.retained_bins(end), ...
            input_info.retained_energy_fraction,input_info.projection_relative, ...
            seconds2,seconds3,NaN,NaN, ...
            relative_error(real(eta2),reference2(:,probe),score), ...
            relative_error(real(eta3),reference3(:,probe),score), ...
            audit2.wavevector_projection_rms_relative, ...
            audit3.wavevector_projection_rms_relative, ...
            audit2.work_frequency_count,audit3.work_frequency_count};
    end
    reference_band2=band2{probe_index,4};reference_band3=band3{probe_index,4};
    for count_index=1:4
        temporary{count_index,9}=relative_change( ...
            band2{probe_index,count_index},reference_band2,score);
        temporary{count_index,10}=relative_change( ...
            band3{probe_index,count_index},reference_band3,score);
    end
    band_rows=[band_rows;temporary]; %#ok<AGROW>
    band_metrics=cell2table(band_rows,'VariableNames',{ ...
        'probe','minimum_count','first_bin','last_bin','retained_energy', ...
        'first_projection_relative','eta22_seconds','eta33_seconds', ...
        'eta22_change_to_36','eta33_change_to_36', ...
        'eta22_HOS_relative','eta33_HOS_relative', ...
        'eta22_q_projection_rms','eta33_q_projection_rms', ...
        'eta22_work_frequency_count','eta33_work_frequency_count'});
    writetable(band_metrics,fullfile(out,'band_metrics.csv'));
end

grid2=cell(2,2);grid3=cell(2,2);grid_rows=cell(0,12);grid_sizes=[256,512];
for probe_index=1:2
    probe=probes(probe_index);score=window_mask(t,limits(probe_index,:));
    model=models{probe_index,1};
    [qx2,qy2]=fixed_limits(model,s.h,2,256);
    [qx3,qy3]=fixed_limits(model,s.h,3,256);
    temporary=cell(2,12);
    for grid_index=1:2
        grid_size=grid_sizes(grid_index);
        fprintf('GRID p%d size=%d\n',probe,grid_size);
        options2=struct('Nx',grid_size,'Ny',grid_size,'J',8, ...
            'QxLimit',qx2,'QyLimit',qy2,'Baseband',true);
        options3=struct('Nx',grid_size,'Ny',grid_size,'J',8, ...
            'QxLimit',qx3,'QyLimit',qy3,'Baseband',true);
        clock=tic;[eta2,audit2]=gl_directional_time_modal_grid( ...
            model.A,model.omega,model.kx,model.ky,s.g,s.h,s.kp,t-t(1),options2);
        seconds2=toc(clock);
        clock=tic;[eta3,audit3]=gl_directional_time_modal_grid_eta33( ...
            model.A,model.omega,model.kx,model.ky,s.g,s.h,s.kp,t-t(1),options3);
        seconds3=toc(clock);
        grid2{probe_index,grid_index}=real(eta2);
        grid3{probe_index,grid_index}=real(eta3);
        temporary(grid_index,:)={probe,grid_size,seconds2,seconds3,NaN,NaN, ...
            relative_error(real(eta2),reference2(:,probe),score), ...
            relative_error(real(eta3),reference3(:,probe),score), ...
            audit2.wavevector_projection_rms_relative, ...
            audit3.wavevector_projection_rms_relative, ...
            audit2.work_frequency_count,audit3.work_frequency_count};
    end
    for grid_index=1:2
        temporary{grid_index,5}=relative_change( ...
            grid2{probe_index,grid_index},grid2{probe_index,2},score);
        temporary{grid_index,6}=relative_change( ...
            grid3{probe_index,grid_index},grid3{probe_index,2},score);
    end
    grid_rows=[grid_rows;temporary]; %#ok<AGROW>
    grid_metrics=cell2table(grid_rows,'VariableNames',{ ...
        'probe','grid_size','eta22_seconds','eta33_seconds', ...
        'eta22_change_to_512','eta33_change_to_512', ...
        'eta22_HOS_relative','eta33_HOS_relative', ...
        'eta22_q_projection_rms','eta33_q_projection_rms', ...
        'eta22_work_frequency_count','eta33_work_frequency_count'});
    writetable(grid_metrics,fullfile(out,'grid_metrics.csv'));
end

rank2=cell(2,3);rank3=cell(2,3);rank_rows=cell(0,10);ranks=[6,8,10];
for probe_index=1:2
    probe=probes(probe_index);score=window_mask(t,limits(probe_index,:));
    model=models{probe_index,1};
    [qx2,qy2]=fixed_limits(model,s.h,2,256);
    [qx3,qy3]=fixed_limits(model,s.h,3,256);
    temporary=cell(3,10);
    for rank_index=1:3
        rank=ranks(rank_index);fprintf('RANK p%d J=%d\n',probe,rank);
        options2=struct('Nx',256,'Ny',256,'J',rank, ...
            'QxLimit',qx2,'QyLimit',qy2,'Baseband',true);
        options3=struct('Nx',256,'Ny',256,'J',rank, ...
            'QxLimit',qx3,'QyLimit',qy3,'Baseband',true);
        clock=tic;eta2=gl_directional_time_modal_grid( ...
            model.A,model.omega,model.kx,model.ky,s.g,s.h,s.kp,t-t(1),options2);
        seconds2=toc(clock);
        clock=tic;eta3=gl_directional_time_modal_grid_eta33( ...
            model.A,model.omega,model.kx,model.ky,s.g,s.h,s.kp,t-t(1),options3);
        seconds3=toc(clock);
        rank2{probe_index,rank_index}=real(eta2);
        rank3{probe_index,rank_index}=real(eta3);
        temporary(rank_index,:)={probe,rank,seconds2,seconds3,NaN,NaN, ...
            relative_error(real(eta2),reference2(:,probe),score), ...
            relative_error(real(eta3),reference3(:,probe),score),39,58};
    end
    for rank_index=1:3
        temporary{rank_index,5}=relative_change( ...
            rank2{probe_index,rank_index},rank2{probe_index,3},score);
        temporary{rank_index,6}=relative_change( ...
            rank3{probe_index,rank_index},rank3{probe_index,3},score);
    end
    rank_rows=[rank_rows;temporary]; %#ok<AGROW>
    rank_metrics=cell2table(rank_rows,'VariableNames',{ ...
        'probe','GL_rank','eta22_seconds','eta33_seconds', ...
        'eta22_change_to_GL10','eta33_change_to_GL10', ...
        'eta22_HOS_relative','eta33_HOS_relative', ...
        'eta22_work_frequency_count','eta33_work_frequency_count'});
    writetable(rank_metrics,fullfile(out,'rank_metrics.csv'));
end

report=struct('probes',probes,'minimum_frequency_counts',counts, ...
    'grid_sizes',grid_sizes,'GL_ranks',ranks,'main_window_seconds',limits, ...
    'band_metrics',table2struct(band_metrics), ...
    'grid_metrics',table2struct(grid_metrics), ...
    'rank_metrics',table2struct(rank_metrics), ...
    'baseband_enabled',true,'pair_loops',0,'triple_loops',0, ...
    'selection_independent_of_HOS_error',true);
writejson(fullfile(out,'report.json'),report);
save(fullfile(out,'comparison.mat'),'report','band_metrics','grid_metrics', ...
    'rank_metrics','band2','band3','grid2','grid3','rank2','rank3', ...
    't','reference2','reference3','first','limits');

figure_handle=figure('Visible','off','Color','w','Position',[100,100,1200,750]);
tiledlayout(2,2);
nexttile;plot_band(band_metrics,'minimum_count','eta22_change_to_36','eta22 band change');
nexttile;plot_band(band_metrics,'minimum_count','eta33_change_to_36','eta33 band change');
nexttile;plot_band(grid_metrics,'grid_size','eta22_change_to_512','eta22 grid change');
nexttile;plot_band(grid_metrics,'grid_size','eta33_change_to_512','eta33 grid change');
exportgraphics(figure_handle,fullfile(out,'convergence.png'),'Resolution',160);
close(figure_handle);
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

function plot_band(table_data,x_column,y_column,label)
probes=unique(table_data.probe);
x=table_data.(x_column);y=table_data.(y_column);
hold on
for index=1:numel(probes)
    use=table_data.probe==probes(index);
    plot(x(use),100*y(use),'-o', ...
        'DisplayName',sprintf('Probe %d',probes(index)));
end
grid on;xlabel(strrep(x_column,'_',' '));ylabel([label,' (%)']);legend('Location','best');
end

function writejson(file,value)
fid=fopen(file,'w');assert(fid>=0);fprintf(fid,'%s\n',jsonencode(value));fclose(fid);
end
