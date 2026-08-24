function output_path=figure_05_psi33_mf12_waveform( ...
        mf12_root,figure_root,quadrature_rank)
%FIGURE_05_PSI33_MF12_WAVEFORM Figure from matched MF12/GL fields.
%   Recomputes the ranked direct graph from the eta11 spectrum saved with MF12.
%   Both the displayed fields and error metrics use the physical completion
%   2*real(Psi33_plus). No gain, phase, time, or spatial alignment is applied.
arguments
    mf12_root (1,1) string
    figure_root (1,1) string = string(fullfile('results','paper'))
    quadrature_rank (1,1) double {mustBeInteger,mustBePositive} = 6
end
if ~isfolder(mf12_root),error('MF12 field root not found: %s',mf12_root);end
if ~isfolder(figure_root),mkdir(figure_root);end
project_root=string(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(project_root,'paper','internal'));

case_ids=["kph-1p0-spread-30-cross-000", ...
    "kph-1p0-spread-30-cross-120"];
crossings=[0,120]; mf12_parent_count=3000; domain_wavelengths=NaN;
metrics=repmat(struct('crossing_angle_deg',NaN, ...
    'physical_raw_relative_l2_pct',NaN, ...
    'physical_raw_relative_linf_pct',NaN,'physical_peak_change_pct',NaN, ...
    'gl_over_mf12_norm_ratio',NaN),numel(case_ids),1);
fields=cell(numel(case_ids),2);
for case_index=1:numel(case_ids)
    mf12_path=fullfile(mf12_root,case_ids(case_index),'fields.mat');
    if ~isfile(mf12_path)
        error('A required fixed-domain field is missing for %s.',case_ids(case_index));
    end
    loaded=load(mf12_path,'fields');source=loaded.fields;
    n=size(source.eta11_spectrum,1);
    if n~=512||size(source.eta11_spectrum,2)~=n|| ...
            source.input_stats.executed_component_count~=mf12_parent_count
        error('The matched-input grid/parent-count contract is not satisfied.');
    end
    if isnan(domain_wavelengths)
        domain_wavelengths=source.config.domain_wavelengths;
    elseif source.config.domain_wavelengths~=domain_wavelengths
        error('The source cases do not share one physical domain.');
    end
    depth=source.config.peak_kh/source.config.peak_wavenumber;
    dk=source.config.delta_k_over_kp*source.config.peak_wavenumber;
    fft_axis=[0:(n/2-1),-n/2:-1];
    [mode_x,mode_y]=meshgrid(fft_axis,fft_axis);
    qx=depth*dk*mode_x;qy=depth*dk*mode_y;
    [~,~,~,candidate_analytic]=green_laplace_eta33_ranked( ...
        source.eta11_spectrum,qx,qy,source.config.peak_kh, ...
        quadrature_rank,project_root);
    reference_physical=2*real(source.reference_Psi33_plus);
    candidate_physical=2*real(candidate_analytic);
    difference=candidate_physical-reference_physical;
    reference_peak=max(abs(reference_physical),[],'all');
    candidate_peak=max(abs(candidate_physical),[],'all');
    metrics(case_index)=struct('crossing_angle_deg',crossings(case_index), ...
        'physical_raw_relative_l2_pct', ...
        100*norm(difference(:))/norm(reference_physical(:)), ...
        'physical_raw_relative_linf_pct', ...
        100*max(abs(difference),[],'all')/reference_peak, ...
        'physical_peak_change_pct', ...
        100*abs(candidate_peak-reference_peak)/reference_peak, ...
        'gl_over_mf12_norm_ratio', ...
        norm(candidate_physical(:))/norm(reference_physical(:)));
    fields{case_index,1}=reference_physical;
    fields{case_index,2}=candidate_physical;
end
metric_table=struct2table(metrics);
writetable(metric_table,fullfile(figure_root, ...
    'fig_psi33_fixed_domain16_waveform_metrics.csv'));

% The 120-degree crossing is the more demanding of the two fixed-domain cases.
reference=fields{2,1};candidate=fields{2,2};
n=size(reference,1);x=((0:n-1)-n/2)*domain_wavelengths/n;
display_window=1.5;
peak_amplitude=0.02;peak_wavenumber=1;gravity=1;depth=1;
peak_omega=sqrt(gravity*peak_wavenumber*tanh(peak_wavenumber*depth));
psi_scale=peak_omega*peak_amplitude^3*peak_wavenumber;
reference_normalized=fftshift(reference)/psi_scale;
candidate_normalized=fftshift(candidate)/psi_scale;
difference_normalized=candidate_normalized-reference_normalized;
reference_peak_normalized=max(abs(reference_normalized),[],'all');
difference_percent=100*difference_normalized/reference_peak_normalized;
center=n/2+1;
centerlines=table(x(:),reference_normalized(center,:).', ...
    candidate_normalized(center,:).',difference_percent(center,:).', ...
    'VariableNames',{'x_over_lambda_p','mf12_m3000_psi33_hat', ...
    sprintf('gl%d_psi33_hat',quadrature_rank), ...
    'difference_pct_of_mf12_peak'});
writetable(centerlines,fullfile(figure_root, ...
    'fig_psi33_fixed_domain16_waveform_centerline.csv'));

field_limit=max(abs([reference_normalized(:);candidate_normalized(:)]));
difference_percent_limit=max(abs(difference_percent),[],'all');
blue=[0.0000 0.4470 0.7410];orange=[0.8500 0.3250 0.0980];
fig=figure('Color','w','Units','centimeters','Position',[2 2 18 10]);
ax1=axes(fig,'Position',[0.035 0.59 0.255 0.34]);
field_image(ax1,x,reference_normalized,field_limit,display_window, ...
    '(a) MF12 ($M=3000$)',false);
ax2=axes(fig,'Position',[0.300 0.59 0.255 0.34]);
cb2=field_image(ax2,x,candidate_normalized,field_limit,display_window, ...
    sprintf('(b) GL%d',quadrature_rank),true);ylabel(ax2,'');
ax3=axes(fig,'Position',[0.610 0.59 0.255 0.34]);
cb3=field_image(ax3,x,difference_percent,difference_percent_limit,display_window, ...
    '(c) Difference (\%)',true);ylabel(ax3,'');
set(ax1,'Position',[0.035 0.59 0.255 0.34]);
set(ax2,'Position',[0.300 0.59 0.255 0.34]);
set(cb2,'Position',[0.563 0.59 0.013 0.34]);
set(ax3,'Position',[0.610 0.59 0.255 0.34]);
set(cb3,'Position',[0.873 0.59 0.013 0.34]);

ax4=axes(fig,'Position',[0.06 0.08 0.60 0.38]);
h_mf12=plot(ax4,x,reference_normalized(center,:),'-','Color',[0.12 0.12 0.12], ...
    'LineWidth',1.2,'DisplayName','MF12 ($M=3000$)');hold(ax4,'on');
h_gl=plot(ax4,x,candidate_normalized(center,:),'--','Color',blue, ...
    'LineWidth',1.25,'DisplayName',sprintf('GL%d',quadrature_rank));
line_style(ax4,display_window*[-1,1],[-field_limit field_limit], ...
    '(d) Centerline comparison', ...
    '$\widehat{\Psi}_{33}\equiv\Psi_{33}/(\omega_p A^3k_p)$');
legend(ax4,[h_mf12 h_gl],{'MF12 ($M=3000$)',sprintf('GL%d',quadrature_rank)}, ...
    'Location','best','Box','off','Interpreter','latex', ...
    'FontName','CMU Serif','FontSize',8);

ax5=axes(fig,'Position',[0.735 0.08 0.24 0.38]);
plot(ax5,x,difference_percent(center,:),'-','Color',orange,'LineWidth',1.15);
line_style(ax5,display_window*[-1,1], ...
    [-difference_percent_limit difference_percent_limit], ...
    sprintf('(e) Difference; $L_2=%.3f\\%%$', ...
    metric_table.physical_raw_relative_l2_pct(2)), ...
    '$100\,\Delta\Psi_{33}/\max|\Psi_{33}^{\mathrm{MF12}}|\;(\%)$');
for ax=[ax1 ax2 ax3 ax4 ax5]
    axis_style(ax);
end
ax5.YAxis.Exponent=0;ytickformat(ax5,'%.2f');
ax4.YLabel.FontSize=7;ax5.YLabel.FontSize=7;
cb3.Ruler.Exponent=0;
cb3.TickLabels=compose('%.2f',cb3.Ticks);
output_path=fullfile(figure_root,'fig_psi33_fixed_domain16_waveform.png');
exportgraphics(fig,output_path,'Resolution',600);close(fig);
fprintf('PSI33_FIXED_DOMAIN16_WAVEFORM_COMPLETE %s\n',output_path);
end

function cb=field_image(ax,x,field,limit,display_window,title_text,show_colorbar)
imagesc(ax,x,x,field);axis(ax,'xy','equal','tight');
clim(ax,[-limit limit]);colormap(ax,redblue_map(256));
xlim(ax,display_window*[-1,1]);ylim(ax,display_window*[-1,1]);
xlabel(ax,'$x/\lambda_p$');
if show_colorbar
    cb=colorbar(ax);cb.FontName='CMU Serif';
    cb.TickLabelInterpreter='latex';
else
    cb=gobjects(0);
end
ylabel(ax,'$y/\lambda_p$');title(ax,title_text,'FontWeight','normal');
end

function line_style(ax,x_limits,y_limits,title_text,y_text)
grid(ax,'on');box(ax,'on');xlim(ax,x_limits);ylim(ax,y_limits);
xlabel(ax,'$x/\lambda_p$ at $y/\lambda_p=0$');ylabel(ax,y_text);
title(ax,title_text,'FontWeight','normal');
end

function axis_style(ax)
set(ax,'FontName','CMU Serif','FontSize',8,'LineWidth',0.75, ...
    'TickDir','out','Layer','top','TickLabelInterpreter','latex');
ax.XLabel.Interpreter='latex';ax.YLabel.Interpreter='latex';
ax.Title.Interpreter='latex';
ax.XLabel.FontName='CMU Serif';ax.YLabel.FontName='CMU Serif';
ax.Title.FontName='CMU Serif';
end

function map=redblue_map(n)
half=floor(n/2);
map=[linspace(0.15,1,half).',linspace(0.35,1,half).',ones(half,1); ...
    ones(n-half,1),linspace(1,0.25,n-half).',linspace(1,0.15,n-half).'];
end
