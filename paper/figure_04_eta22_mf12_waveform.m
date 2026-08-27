function output_path=figure_04_eta22_mf12_waveform( ...
        input_root,output_root,figure_root)
%FIGURE_04_ETA22_MF12_WAVEFORM GL6/MF12 120-degree crossing field.
arguments
    input_root (1,1) string = string(fullfile('artifacts', ...
        'eta22_playground_rank_inputs'))
    output_root (1,1) string = string(fullfile('artifacts', ...
        'eta22_gl6_cross120_s30_waveform'))
    figure_root (1,1) string = string(fullfile('results','paper'))
end
if ~isfolder(output_root),mkdir(output_root);end
if ~isfolder(figure_root),mkdir(figure_root);end
project_root=string(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(project_root,'src','internal'));

case_id="kph1_crossing_a120_s30";
state_path=fullfile(input_root,case_id,'state_and_reference.mat');
if ~isfile(state_path)
    error(['Missing matched 120-degree rank-ladder input. Run ' ...
        'generate_eta22_rank_data with no output cutoff first.']);
end
saved=load(state_path,'state','mf12_reference','mf12_audit');
state=saved.state; depth=state.depth;
if state.config.crossing_angle_deg~=120 || state.config.spread_deg~=30
    error('The saved input is not the declared 120-degree/30-degree case.');
end
[ny,nx]=size(state.eta11_analytic_spectrum);
qx=depth*state.kx; qy=depth*state.ky;
[candidate_dimensionless,audit]=green_laplace_eta22( ...
    state.eta11_analytic_spectrum,qx,qy,depth,6,project_root);
reference=real(saved.mf12_reference);
candidate=real(depth*candidate_dimensionless);
difference=candidate-reference;
reference_norm=norm(reference(:));
candidate_norm=norm(candidate(:));
difference_norm=norm(difference(:));
raw_l2=difference_norm/reference_norm;
raw_linf=max(abs(difference(:)))/max(abs(reference(:)));
norm_ratio=candidate_norm/reference_norm;
surface_similarity_Q=difference_norm/(reference_norm+candidate_norm);

x=((0:nx-1)-floor(nx/2))*state.config.domain_wavelengths/nx;
y=((0:ny-1)-floor(ny/2))*state.config.domain_wavelengths/ny;
peak_amplitude=0.10;peak_wavenumber=1;
eta_scale=peak_amplitude^2*peak_wavenumber;
reference_normalized=fftshift(reference)/eta_scale;
candidate_normalized=fftshift(candidate)/eta_scale;
difference_normalized=fftshift(difference)/eta_scale;
reference_peak_normalized=max(abs(reference_normalized),[],'all');
difference_percent=100*difference_normalized/reference_peak_normalized;
field_limit=max(abs([reference_normalized(:);candidate_normalized(:)]));
difference_percent_limit=max(abs(difference_percent(:)));
orange=[0.8500 0.3250 0.0980]; blue=[0.0000 0.4470 0.7410];
display_window=1.0;

fig=figure('Color','w','Units','centimeters','Position',[2 2 18 10]);
ax1=axes(fig,'Position',[0.035 0.59 0.255 0.34]);
field_image(ax1,x,y,reference_normalized,field_limit,display_window, ...
    '(a) MF12 reference',false);
ax2=axes(fig,'Position',[0.300 0.59 0.255 0.34]);
cb2=field_image(ax2,x,y,candidate_normalized,field_limit,display_window, ...
    '(b) GL6',true);
ylabel(ax2,'');
ax3=axes(fig,'Position',[0.610 0.59 0.255 0.34]);
cb3=field_image(ax3,x,y,difference_percent,difference_percent_limit, ...
    display_window,'(c) Difference',true); ylabel(ax3,'');
set(ax1,'Position',[0.035 0.59 0.255 0.34]);
set(ax2,'Position',[0.300 0.59 0.255 0.34]);
set(cb2,'Position',[0.563 0.59 0.013 0.34]);
set(ax3,'Position',[0.610 0.59 0.255 0.34]);
set(cb3,'Position',[0.873 0.59 0.013 0.34]);
set_colorbar_label(fig,[0.592 0.59 0.014 0.34], ...
    '$\eta_{22}/(A^2k_p)$');
set_colorbar_label(fig,[0.892 0.59 0.075 0.34], ...
    '$\mathrm{Difference}\ (\%\ \mathrm{of\ MF12\ peak})$');

row=floor(ny/2)+1;
ax4=axes(fig,'Position',[0.06 0.08 0.60 0.38]);
h_reference=plot(ax4,x,reference_normalized(row,:),'-','Color',[0.15 0.15 0.15], ...
    'LineWidth',1.25,'DisplayName','MF12'); hold(ax4,'on');
h_gl6=plot(ax4,x,candidate_normalized(row,:),'--','Color',blue, ...
    'LineWidth',1.25,'DisplayName','GL6');
line_style(ax4,display_window*[-1 1],[-field_limit field_limit], ...
    '(d) Centerline comparison', ...
    '$\eta_{22}/(A^2k_p)$');
legend(ax4,[h_reference h_gl6],{'MF12','GL6'},'Location','best', ...
    'Box','off','Interpreter','latex','FontName','CMU Serif','FontSize',8);
ax5=axes(fig,'Position',[0.735 0.08 0.24 0.38]);
plot(ax5,x,difference_percent(row,:),'-','Color',orange,'LineWidth',1.15);
line_style(ax5,display_window*[-1 1], ...
    [-difference_percent_limit difference_percent_limit], ...
    sprintf('(e) Centerline error; $Q=%.4f$',surface_similarity_Q), ...
    '$\mathrm{Difference}\ (\%\ \mathrm{of\ MF12\ peak})$');
for ax=[ax1 ax2 ax3 ax4 ax5]
    set(ax,'FontName','CMU Serif','FontSize',8,'LineWidth',0.75, ...
        'TickDir','out','Layer','top');
    latex_axes(ax);
    disableDefaultInteractivity(ax);
    ax.Toolbar.Visible='off';
end
ax5.YAxis.Exponent=0;ytickformat(ax5,'%.2f');
ax4.YLabel.FontSize=7;ax5.YLabel.FontSize=7;

output_path=fullfile(figure_root,'fig_eta22_gl6_waveform.png');
exportgraphics(fig,output_path,'Resolution',600);
centerlines=table(x(:),reference_normalized(row,:).', ...
    candidate_normalized(row,:).',difference_percent(row,:).', ...
    'VariableNames',{'x_over_lambda_p','mf12_eta22_hat', ...
    'gl6_eta22_hat','difference_pct_of_mf12_peak'});
writetable(centerlines,fullfile(output_root,'centerlines.csv'));
writetable(centerlines,fullfile(figure_root, ...
    'fig_eta22_gl6_waveform_centerline.csv'));
metrics=table(raw_l2,100*raw_l2,raw_linf,100*raw_linf,norm_ratio, ...
    surface_similarity_Q,saved.mf12_audit.active_positive_parent_count, ...
    'VariableNames',{'physical_raw_relative_l2', ...
    'physical_raw_relative_l2_pct','physical_raw_relative_linf', ...
    'physical_raw_relative_linf_pct','gl_over_mf12_norm_ratio', ...
    'surface_similarity_Q','mf12_active_positive_parent_count'});
writetable(metrics,fullfile(figure_root, ...
    'fig_eta22_gl6_waveform_metrics.csv'));
summary=struct('schema_version',1,'status', ...
    'post-freeze-nonselecting-representative-gl6-waveform', ...
    'case_id',case_id,'quadrature_rank',6, ...
    'peak_depth_kph',depth,'crossing_angle_deg',120, ...
    'nominal_spread_deg',30, ...
    'mf12_active_positive_parent_count', ...
        saved.mf12_audit.active_positive_parent_count, ...
    'display_normalization','eta22_hat=eta22/(A^2*k_p)', ...
    'difference_display','100*(GL6-MF12)/max(abs(MF12))', ...
    'peak_amplitude',peak_amplitude,'peak_wavenumber',peak_wavenumber, ...
    'raw_relative_l2',raw_l2,'raw_relative_linf',raw_linf, ...
    'candidate_over_reference_norm_ratio',norm_ratio, ...
    'surface_similarity_Q',surface_similarity_Q, ...
    'fft_ifft_count',audit.fft_ifft_count, ...
    'pointwise_product_count',audit.pointwise_product_count, ...
    'low_output_cutoff_applied',false, ...
    'rescaling_alignment_or_gain',false,'matlab_release',version('-release'));
write_json(fullfile(output_root,'summary.json'),summary);
save(fullfile(output_root,'fields.mat'),'reference','candidate', ...
    'difference','x','y','summary','-v7.3');
close(fig);
fprintf(['GL6_WAVEFORM_COMPLETE cross=120 spread=30 ' ...
    'Q=%.9g rawL2=%.9g%% %s\n'],surface_similarity_Q,100*raw_l2,output_path);
end

function cb=field_image(ax,x,y,z,limit,display_window,title_text,show_colorbar)
imagesc(ax,x,y,z); axis(ax,'image'); set(ax,'YDir','normal');
clim(ax,[-limit limit]); colormap(ax,redblue_map(256));
xlim(ax,display_window*[-1 1]);ylim(ax,display_window*[-1 1]);
if show_colorbar
    cb=colorbar(ax); cb.FontName='CMU Serif';
    cb.TickLabelInterpreter='latex';
else
    cb=gobjects(0);
end
xlabel(ax,'$x/\lambda_p$'); ylabel(ax,'$y/\lambda_p$');
title(ax,title_text,'FontWeight','normal');
end

function line_style(ax,x_limits,y_limits,title_text,y_text)
grid(ax,'on'); box(ax,'on'); xlim(ax,x_limits); ylim(ax,y_limits);
xlabel(ax,'$x/\lambda_p$ at $y/\lambda_p=0$'); ylabel(ax,y_text);
title(ax,title_text,'FontWeight','normal');
end

function set_colorbar_label(fig,position,label_text)
label_axes=axes(fig,'Position',position,'Visible','off','HitTest','off');
text(label_axes,0.5,0.5,label_text,'Units','normalized', ...
    'Interpreter','latex','FontName','CMU Serif','FontSize',7, ...
    'Rotation',90,'HorizontalAlignment','center','VerticalAlignment','middle');
end

function latex_axes(ax)
ax.TickLabelInterpreter='latex';
ax.XLabel.Interpreter='latex'; ax.YLabel.Interpreter='latex';
ax.Title.Interpreter='latex';
ax.XLabel.FontName='CMU Serif'; ax.YLabel.FontName='CMU Serif';
ax.Title.FontName='CMU Serif';
end

function map=redblue_map(n)
half=floor(n/2);
map=[linspace(0.15,1,half).',linspace(0.35,1,half).',ones(half,1); ...
    ones(n-half,1),linspace(1,0.25,n-half).',linspace(1,0.15,n-half).'];
end

function write_json(path,value)
fid=fopen(path,'w'); if fid<0,error('Cannot write %s.',path);end
cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(value,PrettyPrint=true));
end
