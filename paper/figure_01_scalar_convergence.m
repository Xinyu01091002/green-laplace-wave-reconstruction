function outputs = figure_01_scalar_convergence(output_directory)
%FIGURE_01_SCALAR_CONVERGENCE Modal-resolvent convergence.
% Plots relative modal-resolvent errors at the peak-centered determinant-root
% scale.  The second-order subharmonic curve uses the regular difference-
% frequency pair km=1.5*kp, kn=0.5*kp, beta=0.  The p=2:5 curves are
% equal-parent superharmonics.  These are not full-kernel or field errors.

arguments
    output_directory (1,1) string = string(fullfile( ...
        fileparts(fileparts(mfilename('fullpath'))), ...
        'results','paper'))
end
if ~isfolder(output_directory), mkdir(output_directory); end
project_root = string(fileparts(fileparts(mfilename('fullpath'))));
paper_figure_directory = fullfile(project_root,'results','paper');
if ~isfolder(paper_figure_directory), mkdir(paper_figure_directory); end

finite_peak_depth = 1;
subharmonic_parent_ratios = [1.5,0.5];
sector_orders = [0,2,3,4,5];
quadrature_ranks = 1:8;
display_ranks = [2,4,6,8];
display_indices = ismember(quadrature_ranks,display_ranks);
colors = [0.1500,0.1500,0.1500; ...
          0.0000,0.4470,0.7410; ...
          0.8500,0.3250,0.0980; ...
          0.4660,0.6740,0.1880; ...
          0.4940,0.1840,0.5560];
markers = {'o','s','^','d','v'};

depth_cases = [finite_peak_depth,Inf];
detroot_errors = zeros(numel(depth_cases),numel(sector_orders), ...
    numel(quadrature_ranks));
standard_errors = zeros(size(detroot_errors));
r_squared_values = zeros(numel(depth_cases),numel(sector_orders));

for depth_index = 1:numel(depth_cases)
    peak_depth = depth_cases(depth_index);
    for order_index = 1:numel(sector_orders)
        n = sector_orders(order_index);
        if n == 0
            high_ratio = subharmonic_parent_ratios(1);
            low_ratio = subharmonic_parent_ratios(2);
            if isfinite(peak_depth)
                parent_high = sqrt(high_ratio*peak_depth* ...
                    tanh(high_ratio*peak_depth));
                parent_low = sqrt(low_ratio*peak_depth* ...
                    tanh(low_ratio*peak_depth));
                output_frequency = sqrt((high_ratio-low_ratio)*peak_depth* ...
                    tanh((high_ratio-low_ratio)*peak_depth));
            else
                parent_high = sqrt(high_ratio);
                parent_low = sqrt(low_ratio);
                output_frequency = sqrt(high_ratio-low_ratio);
            end
            r = abs(parent_high-parent_low)/output_frequency;
        elseif isfinite(peak_depth)
            r = sqrt(tanh(n*peak_depth)/(n*tanh(peak_depth)));
        else
            r = 1/sqrt(n);
        end
        r_squared = r^2;
        kappa = sqrt(1-r_squared);
        exact = 1/(1-r_squared);
        r_squared_values(depth_index,order_index) = r_squared;
        for rank_index = 1:numel(quadrature_ranks)
            J = quadrature_ranks(rank_index);
            [nodes,weights] = gauss_laguerre_rule(J);
            standard = sum(weights.*sinh(r*nodes)/r);
            slow = exp((1-(1-r)/kappa)*nodes);
            fast = exp((1-(1+r)/kappa)*nodes);
            determinant_root = sum(weights.*(slow-fast))/(2*r*kappa);
            standard_errors(depth_index,order_index,rank_index) = ...
                abs(standard-exact)/exact;
            detroot_errors(depth_index,order_index,rank_index) = ...
                abs(determinant_root-exact)/exact;
        end
    end
end

figure_handle = figure('Color','w','Units','centimeters', ...
    'Position',[2,2,18.0,6.6],'PaperPositionMode','auto');
layout = tiledlayout(figure_handle,1,2,'Padding','compact', ...
    'TileSpacing','compact');
axes_handles = gobjects(1,2);
for depth_index = 1:numel(depth_cases)
    axes_handles(depth_index) = nexttile(layout,depth_index);
    current_axis = axes_handles(depth_index);
    hold(current_axis,'on');
    set(current_axis,'YScale','log');
    for order_index = 1:numel(sector_orders)
        plotted_errors = squeeze(detroot_errors(depth_index,order_index, ...
            display_indices));
        if sector_orders(order_index) == 0
            display_name = 'Second-order subharmonic';
        else
            display_name = sprintf('$p=%d$',sector_orders(order_index));
        end
        semilogy(current_axis,display_ranks,100*plotted_errors, ...
            'LineWidth',1.15,'Marker',markers{order_index}, ...
            'MarkerFaceColor','w','MarkerSize',4.5, ...
            'Color',colors(order_index,:), ...
            'DisplayName',display_name);
    end
    xlabel(current_axis,'Gauss--Laguerre order','Interpreter','latex');
    if depth_index == 1
        ylabel(current_axis,'Relative modal-resolvent error (\%)', ...
            'Interpreter','latex');
    end
    xlim(current_axis,[1.7,8.3]);
    ylim(current_axis,[1e-8,20]);
    xticks(current_axis,display_ranks);
    xticklabels(current_axis,compose('GL%d',display_ranks));
    set(current_axis,'FontName','CMU Serif','FontSize',8, ...
        'LineWidth',0.75,'TickDir','out','Layer','top','Box','on', ...
        'XMinorGrid','off','YMinorGrid','off');
    grid(current_axis,'on');
    set(current_axis,'GridAlpha',0.14,'TickLabelInterpreter','latex');
    current_axis.Toolbar.Visible = 'off';
    if depth_index == 1
        panel_label = '(a) $k_p d=1$';
    else
        panel_label = '(b) Deep-water limit';
    end
    text(current_axis,0.96,0.95,panel_label,'Units','normalized', ...
        'HorizontalAlignment','right','VerticalAlignment','top', ...
        'Interpreter','latex','FontName','CMU Serif','FontSize',8);
end
legend_handle = legend(axes_handles(1),'Orientation','horizontal', ...
    'NumColumns',5,'Interpreter','latex','Box','off','FontSize',8, ...
    'FontName','CMU Serif');
legend_handle.Layout.Tile = 'north';

png_path = fullfile(output_directory, ...
    'bound_harmonic_actual_convergence_ce.png');
exportgraphics(figure_handle,png_path,'Resolution',1000);
paper_png_path = fullfile(paper_figure_directory, ...
    'fig_gl_scalar_convergence.png');
exportgraphics(figure_handle,paper_png_path,'Resolution',1000);
close(figure_handle);

rows = zeros(numel(depth_cases)*numel(sector_orders)* ...
    numel(quadrature_ranks),7);
row_index = 0;
for depth_index = 1:numel(depth_cases)
    for order_index = 1:numel(sector_orders)
        for rank_index = 1:numel(quadrature_ranks)
            row_index = row_index+1;
            rows(row_index,:) = [depth_index,depth_cases(depth_index), ...
                sector_orders(order_index),quadrature_ranks(rank_index), ...
                r_squared_values(depth_index,order_index), ...
                100*standard_errors(depth_index,order_index,rank_index), ...
                100*detroot_errors(depth_index,order_index,rank_index)];
        end
    end
end
data_table = array2table(rows,'VariableNames', ...
    {'depth_case','peak_kph','n','J','r_squared', ...
    'standard_actual_percent','detroot_actual_percent'});
csv_path = fullfile(output_directory, ...
    'bound_harmonic_actual_convergence_ce_data.csv');
writetable(data_table,csv_path);
paper_csv_path = fullfile(paper_figure_directory, ...
    'fig_gl_scalar_convergence_data.csv');
writetable(data_table,paper_csv_path);

outputs = struct('status','complete-ce-main-modal-resolvent-figure', ...
    'scope',['relative modal-resolvent errors; representative second-order ' ...
    'subharmonic pair; p=2:5 equal-parent superharmonics'], ...
    'peak_kph',depth_cases,'sector_orders',sector_orders, ...
    'subharmonic_parent_ratios',subharmonic_parent_ratios, ...
    'quadrature_ranks',quadrature_ranks,'double_column_width_mm',180, ...
    'png_dpi',1000,'png',png_path,'csv',csv_path, ...
    'paper_png',paper_png_path, ...
    'paper_csv',paper_csv_path);
fprintf('BOUND_HARMONIC_ACTUAL_CE_FIGURE=COMPLETE\n');
fprintf('PNG=%s\n',png_path);
end

function [nodes,weights] = gauss_laguerre_rule(J)
diagonal = (2*(1:J)-1).';
off_diagonal = (1:J-1).';
jacobi = diag(diagonal)+diag(off_diagonal,1)+diag(off_diagonal,-1);
[vectors,values] = eig(jacobi,'vector');
[nodes,order] = sort(values.');
vectors = vectors(:,order);
weights = vectors(1,:).^2;
end
