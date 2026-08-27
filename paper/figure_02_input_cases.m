function output_path=figure_02_input_cases(output_root,domain_wavelengths, ...
        maximum_components,show_decomposition)
%FIGURE_02_INPUT_CASES Paper rendering of the eta22 input cases.
arguments
    output_root (1,1) string = fullfile('results','paper')
    domain_wavelengths (1,1) double {mustBePositive} = 10
    maximum_components (1,1) double {mustBeInteger,mustBePositive} = 2500
    show_decomposition (1,1) logical = false
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
    config=playground_config(definitions{case_index},domain_wavelengths, ...
        maximum_components);
    states{case_index}=generate_eta22_input(config);
end

if show_decomposition
    figure_height=15.2;
    layout_rows=7;
else
    figure_height=12.4;
    layout_rows=2;
end
fig=figure('Color','w','Units','centimeters', ...
    'Position',[2 2 17.2 figure_height]);
layout=tiledlayout(fig,layout_rows,3,'TileSpacing','tight','Padding','compact');
layout.OuterPosition=[0.06 0.015 0.88 0.97];
spectrum_colours=directional_spectrum_colormap(256);
for case_index=1:3
    state=states{case_index};
    spectral=continuous_directional_spectrum(definitions{case_index});
    if show_decomposition
        ax=nexttile(layout,case_index,[3 1]);
    else
        ax=nexttile(layout,case_index);
    end
    imagesc(ax,spectral.kx,spectral.ky,spectral.log_density);
    set(ax,'YDir','normal'); hold(ax,'on');
    contour(ax,spectral.kx,spectral.ky,spectral.log_density, ...
        [-3 -2 -1],'LineColor',[0.16 0.16 0.18],'LineWidth',0.55);
    hold(ax,'off');
    if show_decomposition
        xlim(ax,[0 2.8]);
    else
        xlim(ax,[0 2.35]);
    end
    ylim(ax,[-2.1 2.1]); axis(ax,'equal');
    xticks(ax,[0 1 2]); yticks(ax,[-2 -1 0 1 2]);
    grid(ax,'on'); box(ax,'on');
    ax.GridColor=[0.12 0.15 0.18]; ax.GridAlpha=0.13;
    xlabel(ax,'$k_x/k_p$');
    if case_index==1, ylabel(ax,'$k_y/k_p$'); end
    title(ax,{sprintf('(%c) %s',char('a'+case_index-1), ...
        definitions{case_index}.title),definitions{case_index}.subtitle}, ...
        'FontWeight','normal');
    colormap(ax,spectrum_colours); clim(ax,[-4 0]); format_axes(ax);

    n=size(state.eta11,1);
    coordinate=(-n/2:n/2-1)*state.config.domain_wavelengths/n;
    display_field=fftshift(real(state.eta11))*state.depth ...
        /state.config.crest_amplitude;
    center=n/2+1;
    if show_decomposition
        ax=nexttile(layout,case_index+12,[3 1]);
    else
        ax=nexttile(layout,case_index+3);
    end
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
if show_decomposition
    for case_index=1:3
        add_directional_decomposition(layout,case_index+9, ...
            definitions{case_index},case_index==1);
    end
end
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

function add_directional_decomposition(layout,tile_number,definition,show_ylabel)
if definition.crossing==0
    centers=0;
else
    centers=[-0.5 0.5]*definition.crossing;
end
component_colours=[0.00 0.50 0.75;0.92 0.42 0.08];
if numel(centers)==1
    component_colours=component_colours(1,:);
end
ax=nexttile(layout,tile_number);
theta=linspace(-120,120,961);
component_values=zeros(numel(centers),numel(theta));
for component_index=1:numel(centers)
    component_values(component_index,:)=exp(-0.5* ...
        ((theta-centers(component_index))/definition.spread).^2) ...
        /numel(centers);
end
total=sum(component_values,1);
normalization=max(total);
hold(ax,'on');
if numel(centers)>1
    for component_index=1:numel(centers)
        plot(ax,theta,component_values(component_index,:)/normalization, ...
            '--','Color',component_colours(component_index,:),'LineWidth',0.85);
    end
end
plot(ax,theta,total/normalization,'k-','LineWidth',1.05);
hold(ax,'off');
xlim(ax,[-120 120]); ylim(ax,[0 1.08]);
xticks(ax,[-90 0 90]); yticks(ax,[0 1]);
xlabel(ax,'$\theta$');
if show_ylabel, ylabel(ax,'$D(\theta)$'); else, yticklabels(ax,{}); end
grid(ax,'on'); box(ax,'on');
ax.GridColor=[0.25 0.25 0.28]; ax.GridAlpha=0.10;
format_axes(ax);
end

function colours=directional_spectrum_colormap(number_of_colours)
% Portable blue-white-red approximation to matplotlib's coolwarm map.
anchors=[0.230 0.299 0.754;0.554 0.690 0.996;0.865 0.865 0.865; ...
    0.957 0.598 0.477;0.706 0.016 0.150];
colours=interp1(linspace(0,1,size(anchors,1)),anchors, ...
    linspace(0,1,number_of_colours),'pchip');
colours=max(0,min(1,colours));
end

function config=playground_config(definition,domain_wavelengths, ...
        maximum_components)
config=struct('retained_mass_percent',99,'jonswap_gamma',3.3, ...
    'spread_deg',definition.spread, ...
    'crossing_angle_deg',definition.crossing,'peak_kh',1, ...
    'domain_wavelengths',domain_wavelengths, ...
    'grid_size',256,'auto_expand_grid',true, ...
    'selection_mode','nested-energy', ...
    'maximum_components',maximum_components, ...
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
