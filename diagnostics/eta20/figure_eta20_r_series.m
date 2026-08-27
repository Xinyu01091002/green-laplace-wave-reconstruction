function output_path=figure_eta20_r_series( ...
        figure_root,mf12_root,peak_kh)
%FIGURE_ETA20_R_SERIES GL12 versus Neumann R2/R4/R6 formulas.
% R2/R4/R6 use independent ordered-pair reconstruction because the high-power
% fixed-FFT expansions are conditioning diagnostics, not production evidence.
arguments
    figure_root (1,1) string = string(fullfile('results','paper'))
    mf12_root (1,1) string = string(getenv('MF12_ROOT'))
    peak_kh (1,1) double {mustBePositive,mustBeFinite} = 1
end
if ~isfolder(figure_root),mkdir(figure_root);end
project_root=string(fileparts(fileparts(fileparts(mfilename('fullpath')))));
addpath(fullfile(project_root,'diagnostics','eta20'));
addpath(fullfile(project_root,'diagnostics','eta20','generated'));
addpath(fullfile(project_root,'paper','internal'));
setup_green_laplace("MF12Root",mf12_root);
mf12_source=string(fileparts(which('mf12_spectral_coefficients')));

config=struct('retained_mass_percent',99,'jonswap_gamma',3.3, ...
    'spread_deg',30,'crossing_angle_deg',120,'peak_kh',peak_kh, ...
    'domain_wavelengths',10,'grid_size',256,'auto_expand_grid',true, ...
    'selection_mode','nested-energy','maximum_components',2500, ...
    'radial_bins',48,'angular_bins',96,'minimum_k_over_kp',0.3, ...
    'maximum_k_over_kp',3.5,'crest_amplitude',0.1, ...
    'time_peak_periods',0,'minimum_output_midpoint_q',0.3);
state=generate_eta22_input(config);
eta11_spectrum=state.eta11_spectrum;kx=state.kx;ky=state.ky;depth=state.depth;
n=size(eta11_spectrum,1);domain_wavelengths=config.domain_wavelengths;
[reference,reference_audit]=mf12_eta20_reference( ...
    eta11_spectrum,kx,ky,depth,mf12_source);
[gl12_half,gl12_audit]=eta20_green_laplace_shared( ...
    eta11_spectrum,kx,ky,depth,"shared12",project_root);
gl12=2*gl12_half;

powers=[2,4,6];r_fields=cell(1,3);r_audits=cell(1,3);
for index=1:numel(powers)
    [r_half,r_audits{index}]= ...
        eta20_neumann_r_series_ordered_pair( ...
        eta11_spectrum,kx,ky,depth,powers(index),project_root);
    r_fields{index}=2*r_half;
end
if ~strcmp(gl12_audit.information_boundary,'eta11-only') ...
        || any(cellfun(@(a)~strcmp(a.information_boundary,'eta11-only'),r_audits))
    error('The GL12/R2/R4/R6 eta11-only information boundary failed.');
end

methods=["GL12";"R2";"R4";"R6"];
candidates=[{gl12},r_fields];method_order=[12;2;4;6];
records=struct([]);reference_norm=norm(reference(:));
reference_peak=max(abs(reference),[],'all');
for index=1:numel(candidates)
    candidate=candidates{index};difference=candidate-reference;
    difference_norm=norm(difference(:));candidate_norm=norm(candidate(:));
    if index==1
        evaluation_graph="fixed-FFT";
    else
        evaluation_graph="ordered-pair-validation";
    end
    record=struct('method',methods(index),'peak_kpd',depth, ...
        'crossing_angle_deg',config.crossing_angle_deg, ...
        'spread_deg',config.spread_deg,'grid_size',n, ...
        'selected_component_count',reference_audit.active_positive_parent_count, ...
        'method_order_or_rank',method_order(index), ...
        'physical_raw_relative_l2_pct',100*difference_norm/reference_norm, ...
        'physical_raw_relative_linf_pct', ...
            100*max(abs(difference),[],'all')/reference_peak, ...
        'candidate_over_mf12_norm_ratio',candidate_norm/reference_norm, ...
        'surface_similarity_Q',difference_norm/(reference_norm+candidate_norm), ...
        'evaluation_graph',evaluation_graph, ...
        'strict_zero_mode','set-to-zero-separate-sector');
    if index==1,records=record;else,records(index)=record;end
end
metrics=struct2table(records);
base=fullfile(figure_root,'fig_eta20_gl12_r2_r4_r6_appendix');
writetable(metrics,base+"_metrics.csv");

carrier_wavenumber=1;field_scale=config.crest_amplitude^2*carrier_wavenumber;
x=(-n/2:n/2-1)*domain_wavelengths/n;center=n/2+1;
reference_n=fftshift(reference)/field_scale;
candidate_n=cellfun(@(f)fftshift(f)/field_scale,candidates,'UniformOutput',false);
reference_peak_n=max(abs(reference_n),[],'all');
difference_pct=cellfun(@(f)100*(f-reference_n)/reference_peak_n, ...
    candidate_n,'UniformOutput',false);
centerlines=table(x(:),reference_n(center,:).',candidate_n{1}(center,:).', ...
    candidate_n{2}(center,:).',candidate_n{3}(center,:).', ...
    candidate_n{4}(center,:).',difference_pct{1}(center,:).', ...
    difference_pct{2}(center,:).',difference_pct{3}(center,:).', ...
    difference_pct{4}(center,:).', ...
    'VariableNames',{'x_over_lambda_p','mf12_eta20_scaled','gl12_eta20_scaled', ...
    'r2_eta20_scaled','r4_eta20_scaled','r6_eta20_scaled', ...
    'gl12_difference_pct','r2_difference_pct','r4_difference_pct', ...
    'r6_difference_pct'});
writetable(centerlines,base+"_centerline.csv");

mode_axis=0:n/2;kx_over_k0=mode_axis/domain_wavelengths;
kx_depth=kx_over_k0*depth;
spectra=cell(1,5);spectra{1}=abs(fft(reference_n(center,:)))/n;
for index=1:4,spectra{index+1}=abs(fft(candidate_n{index}(center,:)))/n;end
for index=1:5,spectra{index}=spectra{index}(1:n/2+1);end
spectrum_table=table(kx_over_k0(:),kx_depth(:),spectra{1}(:),spectra{2}(:), ...
    spectra{3}(:),spectra{4}(:),spectra{5}(:), ...
    'VariableNames',{'abs_kx_over_k0','abs_kx_depth', ...
    'mf12_eta20_spectrum_scaled', ...
    'gl12_eta20_spectrum_scaled','r2_eta20_spectrum_scaled', ...
    'r4_eta20_spectrum_scaled','r6_eta20_spectrum_scaled'});
writetable(spectrum_table,base+"_centerline_spectrum.csv");

display_window=1.5;
field_limit=max(abs(reference_n),[],'all');
difference_limit=max(cellfun(@(f)max(abs(f),[],'all'),difference_pct));
line_limit=max(abs([reference_n(center,:),candidate_n{1}(center,:), ...
    candidate_n{2}(center,:),candidate_n{3}(center,:),candidate_n{4}(center,:)]));
difference_line_limit=max(abs([difference_pct{1}(center,:), ...
    difference_pct{2}(center,:),difference_pct{3}(center,:), ...
    difference_pct{4}(center,:)]));
colors=[0.0000 0.4470 0.7410;0.8500 0.3250 0.0980; ...
    0.4940 0.1840 0.5560;0.4660 0.6740 0.1880];
styles={'--',':','-.','--'};
fig=figure('Color','w','Units','centimeters','Position',[2 2 32 12.5]);
ax1=axes(fig,'Position',[0.025 0.60 0.145 0.32]);
field_image(ax1,x,reference_n,field_limit,display_window,'(a) MF12 reference',false);
cb1=colorbar(ax1);set(cb1,'Position',[0.174 0.60 0.007 0.32]);
set_colorbar_label(fig,[0.181 0.60 0.014 0.32],'$\eta_{20}/(A^2k_0)$');

contour_positions=[0.215 0.400 0.585 0.770];
contour_axes=gobjects(1,4);cb_error=gobjects(0);
for index=1:4
    contour_axes(index)=axes(fig,'Position', ...
        [contour_positions(index) 0.60 0.145 0.32]);
    title_text=sprintf('(%c) %s error, $Q=%.4f$', ...
        'a'+index,methods(index),records(index).surface_similarity_Q);
    show_colorbar=index==4;
    cb=field_image(contour_axes(index),x,difference_pct{index}, ...
        difference_limit,display_window,title_text,show_colorbar);
    if show_colorbar,cb_error=cb;end
end
set(cb_error,'Position',[0.921 0.60 0.007 0.32]);
set_colorbar_label(fig,[0.932 0.60 0.050 0.32], ...
    '$\mathrm{Error}\ (\%\ \mathrm{of\ MF12\ peak})$');

ax6=axes(fig,'Position',[0.045 0.085 0.270 0.37]);
h=plot(ax6,x,reference_n(center,:),'-','Color',[.12 .12 .12], ...
    'LineWidth',1.2,'DisplayName','MF12');hold(ax6,'on');
handles=h;
for index=1:4
    handles(end+1)=plot(ax6,x,candidate_n{index}(center,:),styles{index}, ...
        'Color',colors(index,:),'LineWidth',1.15, ...
        'DisplayName',methods(index)); %#ok<AGROW>
end
line_style(ax6,display_window*[-1,1],[-line_limit line_limit], ...
    '(f) Centerline comparison','$\eta_{20}/(A^2k_0)$');
legend(ax6,handles,'Location','best','Box','off','Interpreter','latex', ...
    'FontName','CMU Serif','FontSize',6.5);

ax7=axes(fig,'Position',[0.375 0.085 0.270 0.37]);hold(ax7,'on');
for index=1:4
    plot(ax7,x,difference_pct{index}(center,:),styles{index}, ...
        'Color',colors(index,:),'LineWidth',1.25, ...
        'DisplayName',sprintf('%s; full-field Q=%.4f', ...
        methods(index),records(index).surface_similarity_Q));
end
line_style(ax7,display_window*[-1,1], ...
    [-difference_line_limit difference_line_limit],'(g) Centerline error', ...
    '$\mathrm{Error}\ (\%\ \mathrm{of\ MF12\ peak})$');
legend(ax7,'Location','best','Box','off','Interpreter','none', ...
    'FontName','CMU Serif','FontSize',6.25);

ax8=axes(fig,'Position',[0.700 0.085 0.270 0.37]);hold(ax8,'on');
spectrum_x_limit=1.5;
visible_band=kx_depth<=spectrum_x_limit;
spectrum_ceiling=1.08*max(cellfun(@(s)max(s(visible_band)),spectra));
low_band=patch(ax8,[0 .3 .3 0],[0 0 spectrum_ceiling spectrum_ceiling], ...
    [.88 .88 .88],'FaceAlpha',.35,'EdgeColor','none', ...
    'HandleVisibility','off');
uistack(low_band,'bottom');
h_spec=plot(ax8,kx_depth,spectra{1},'-', ...
    'Color',[.12 .12 .12],'LineWidth',1.3,'DisplayName','MF12');
spec_handles=h_spec;
for index=1:4
    spec_handles(end+1)=plot(ax8,kx_depth,spectra{index+1},styles{index}, ...
        'Color',colors(index,:),'LineWidth',1.15,'DisplayName',methods(index)); %#ok<AGROW>
end
grid(ax8,'on');box(ax8,'on');xlim(ax8,[0 spectrum_x_limit]);
ylim(ax8,[0 spectrum_ceiling]);set(ax8,'YAxisLocation','right');
xline(ax8,.3,':','Color',[.35 .35 .35],'LineWidth',.9, ...
    'HandleVisibility','off');
text(ax8,.15,.94*spectrum_ceiling,'$|k_x|d<0.3$', ...
    'Interpreter','latex','HorizontalAlignment','center', ...
    'VerticalAlignment','top','FontSize',6.5,'Color',[.3 .3 .3]);
xlabel(ax8,'$|k_x|d$');ylabel(ax8,'$|\widehat{\eta}_{20}|/(A^2k_0)$');
title(ax8,'(h) Centerline spectra','FontWeight','normal');
legend(ax8,spec_handles,'Location','northeast','Box','off', ...
    'Interpreter','latex','FontName','CMU Serif','FontSize',6);

for ax=[ax1 contour_axes ax6 ax7 ax8],axis_style(ax);end
ax7.YAxis.Exponent=0;ytickformat(ax7,'%.2f');
cb_error.Ruler.Exponent=0;
output_path=base+".png";exportgraphics(fig,output_path,'Resolution',600);close(fig);
fprintf(['ETA20_GL12_R246_APPENDIX_COMPLETE Q_GL12=%.9g Q_R2=%.9g ' ...
    'Q_R4=%.9g Q_R6=%.9g %s\n'],records(1).surface_similarity_Q, ...
    records(2).surface_similarity_Q,records(3).surface_similarity_Q, ...
    records(4).surface_similarity_Q,output_path);
end

function [field,audit]=mf12_eta20_reference( ...
        eta11_spectrum,kx,ky,depth,mf12_root)
source_file=fullfile(mf12_root,'mf12_spectral_coefficients.m');
if ~isfile(source_file),error('Official MF12 source not found.');end
addpath(mf12_root,'-begin');resolved=string(which('mf12_spectral_coefficients'));
if ~strcmpi(char(resolved),char(source_file)),error('Unexpected MF12 source.');end
positive=kx>0;threshold=1e-13*max(1,max(abs(eta11_spectrum),[],'all'));
parents=find(positive&abs(eta11_spectrum)>threshold);count=numel(eta11_spectrum);
coefficients_eta=eta11_spectrum(parents)/count;
coefficients=mf12_spectral_coefficients(2,1,depth,2*real(coefficients_eta), ...
    2*imag(coefficients_eta),kx(parents),ky(parents),0,0);
indices=2:2:numel(coefficients.G_npm);
values=(coefficients.A_npm(indices)+1i*coefficients.B_npm(indices)) ...
    .*coefficients.G_npm(indices);
dkx=kx(1,2)-kx(1,1);dky=ky(2,1)-ky(1,1);n=size(eta11_spectrum,1);
columns=mod(round(coefficients.kx_npm(indices)/dkx),n)+1;
rows=mod(round(coefficients.ky_npm(indices)/dky),n)+1;
one_sided=accumarray([rows(:),columns(:)],values(:),[n n],@sum,complex(0));
one_sided(1,1)=0;field=real(ifft2(one_sided))*count;
audit=struct('active_positive_parent_count',numel(parents),'external_source',resolved);
end

function cb=field_image(ax,x,field,limit,window,title_text,show_colorbar)
imagesc(ax,x,x,field);axis(ax,'xy','equal','tight');clim(ax,[-limit limit]);
colormap(ax,redblue_map(256));xlim(ax,window*[-1,1]);ylim(ax,window*[-1,1]);
xlabel(ax,'$x/\lambda_p$');ylabel(ax,'$y/\lambda_p$');
title(ax,title_text,'FontWeight','normal');
if show_colorbar,cb=colorbar(ax);else,cb=gobjects(0);end
end
function line_style(ax,xlimits,ylimits,title_text,ytext)
grid(ax,'on');box(ax,'on');xlim(ax,xlimits);ylim(ax,ylimits);
xlabel(ax,'$x/\lambda_p$ at $y/\lambda_p=0$');ylabel(ax,ytext);
title(ax,title_text,'FontWeight','normal');
end
function set_colorbar_label(fig,position,label_text)
label_axes=axes(fig,'Position',position,'Visible','off','HitTest','off');
text(label_axes,.5,.5,label_text,'Units','normalized','Interpreter','latex', ...
    'FontName','CMU Serif','FontSize',7,'Rotation',90, ...
    'HorizontalAlignment','center','VerticalAlignment','middle');
end
function axis_style(ax)
set(ax,'FontName','CMU Serif','FontSize',8,'LineWidth',.75, ...
    'TickDir','out','Layer','top','TickLabelInterpreter','latex');
ax.XLabel.Interpreter='latex';ax.YLabel.Interpreter='latex';
ax.Title.Interpreter='latex';
end
function map=redblue_map(n)
half=floor(n/2);map=[linspace(.15,1,half).',linspace(.35,1,half).', ...
    ones(half,1);ones(n-half,1),linspace(1,.25,n-half).', ...
    linspace(1,.15,n-half).'];
end
