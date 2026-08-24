function output_path=figure_02_input_cases(output_root)
%FIGURE_02_INPUT_CASES Paper rendering of the eta22 input cases.
arguments
    output_root (1,1) string = string(fullfile( ...
        fileparts(fileparts(mfilename('fullpath'))),'results','paper'))
end
if ~isfolder(output_root), mkdir(output_root); end
project_root=string(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(project_root,'paper','internal'));

definitions={ ...
    struct('title','$0^\circ$ crossing sea', ...
        'subtitle','$30^\circ$ spreading', ...
        'spread',30,'crossing',0), ...
    struct('title','Crossing sea', ...
        'subtitle','$90^\circ$ crossing; $30^\circ$ spreading', ...
        'spread',30,'crossing',90), ...
    struct('title','Crossing sea', ...
        'subtitle','$120^\circ$ crossing; $30^\circ$ spreading', ...
        'spread',30,'crossing',120)};

states=cell(1,3);
for case_index=1:3
    config=playground_config(definitions{case_index});
    states{case_index}=generate_eta22_input(config);
end

fig=figure('Color','w','Units','centimeters','Position',[2 2 20 13.5]);
layout=tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
layout.OuterPosition=[0 0 0.92 1];
spectrum_axes=gobjects(1,3);
for case_index=1:3
    state=states{case_index};
    spectral=continuous_directional_spectrum(definitions{case_index});
    ax=nexttile(case_index); spectrum_axes(case_index)=ax;
    imagesc(ax,spectral.kx,spectral.ky,spectral.log_density);
    set(ax,'YDir','normal'); hold(ax,'on');
    contour(ax,spectral.kx,spectral.ky,spectral.log_density, ...
        [-3 -2 -1],'LineColor',[0.90 0.93 0.96],'LineWidth',0.45);
    hold(ax,'off');
    xlim(ax,[0 3.5]); ylim(ax,[-3.5 3.5]); axis(ax,'equal');
    xticks(ax,[0 1 2 3]); yticks(ax,[-3 -2 -1 0 1 2 3]);
    grid(ax,'on'); box(ax,'on');
    xlabel(ax,'$k_x/k_p$');
    if case_index==1, ylabel(ax,'$k_y/k_p$'); end
    title(ax,{sprintf('(%c) %s',char('a'+case_index-1), ...
        definitions{case_index}.title),definitions{case_index}.subtitle}, ...
        'FontWeight','normal');
    colormap(ax,parula(256)); clim(ax,[-4 0]); format_axes(ax);

    n=size(state.eta11,1);
    coordinate=(-n/2:n/2-1)*state.config.domain_wavelengths/n;
    display_field=fftshift(real(state.eta11))*state.depth ...
        /state.config.crest_amplitude;
    center=n/2+1;
    ax=nexttile(case_index+3);
    surf(ax,coordinate,coordinate,display_field,'EdgeColor','none');
    hold(ax,'on');
    plot3(ax,coordinate,zeros(size(coordinate)),display_field(center,:), ...
        '-','Color',[0.82 0.20 0.16],'LineWidth',1.15);
    hold(ax,'off');
    view(ax,38,27); box(ax,'on'); grid(ax,'on');
    xlim(ax,[-5 5]); ylim(ax,[-5 5]); zlim(ax,[-0.5 1]);
    xlabel(ax,'$x/\lambda_p$'); ylabel(ax,'$y/\lambda_p$');
    if case_index==1, zlabel(ax,'$\eta^{(11)}/A$'); end
    title(ax,sprintf('(%c) Focused linear field',char('d'+case_index-1)), ...
        'FontWeight','normal');
    colormap(ax,parula(256)); format_axes(ax);
end
colourbar_axes=axes(fig,'Position',[0.925 0.59 0.01 0.30], ...
    'Visible','off');
colormap(colourbar_axes,parula(256)); clim(colourbar_axes,[-4 0]);
cb=colorbar(colourbar_axes,'Location','eastoutside');
cb.Position=[0.945 0.59 0.014 0.30];
cb.Ticks=[-4 -3 -2 -1 0];
cb.TickLabels={'10^{-4}','10^{-3}','10^{-2}','10^{-1}','1'};

base=fullfile(output_root,'fig_eta22_linear_input_cases_ce');
exportgraphics(fig,base+".png",'Resolution',900);
close(fig);
output_path=base+".png";
end

function value=continuous_directional_spectrum(definition)
kx=linspace(0,3.5,900); ky=linspace(-3.5,3.5,1200);
[kx_grid,ky_grid]=meshgrid(kx,ky);
radial=hypot(kx_grid,ky_grid); angle=rad2deg(atan2(ky_grid,kx_grid));
radial_density=jonswap_wavenumber_spectrum(radial,1,1,3.3);
if definition.crossing==0
    centers=0;
else
    centers=[-0.5 0.5]*definition.crossing;
end
directional_density=zeros(size(angle));
for center=centers
    directional_density=directional_density+exp( ...
        -0.5*((angle-center)/definition.spread).^2) ...
        /(sqrt(2*pi)*definition.spread*numel(centers));
end
density=radial_density./max(radial,realmin).*directional_density;
density(radial<0.3 | radial>3.5 | kx_grid<=0)=0;
relative_density=density/max(density,[],'all');
value=struct('kx',kx,'ky',ky, ...
    'log_density',log10(max(relative_density,1e-4)));
end

function spectrum=jonswap_wavenumber_spectrum(k,kp,gravity,gamma)
k_safe=max(k,realmin);
frequency=sqrt(gravity*k_safe)/(2*pi);
peak_frequency=sqrt(gravity*kp)/(2*pi);
sigma=0.09*ones(size(k_safe)); sigma(frequency<=peak_frequency)=0.07;
peak_shape=exp(-0.5*((frequency-peak_frequency) ...
    ./(sigma*peak_frequency)).^2);
spectrum_f=frequency.^(-5).*exp(-(5/4)*(peak_frequency./frequency).^4) ...
    .*gamma.^peak_shape;
spectrum=spectrum_f.*sqrt(gravity./k_safe)/(4*pi);
end

function config=playground_config(definition)
config=struct('retained_mass_percent',99,'jonswap_gamma',3.3, ...
    'spread_deg',definition.spread, ...
    'crossing_angle_deg',definition.crossing,'peak_kh',1, ...
    'domain_wavelengths',10,'grid_size',256,'auto_expand_grid',true, ...
    'selection_mode','nested-energy','maximum_components',2500, ...
    'radial_bins',48,'angular_bins',96,'minimum_k_over_kp',0.3, ...
    'maximum_k_over_kp',3.5,'crest_amplitude',0.1, ...
    'time_peak_periods',0,'minimum_output_midpoint_q',0.3);
end

function format_axes(ax)
set(ax,'FontName','CMU Serif','FontSize',7.6,'LineWidth',0.7, ...
    'TickDir','out','Layer','top','TickLabelInterpreter','latex');
disableDefaultInteractivity(ax);
ax.Toolbar.Visible='off';
ax.XLabel.Interpreter='latex'; ax.YLabel.Interpreter='latex';
ax.ZLabel.Interpreter='latex'; ax.Title.Interpreter='latex';
end
