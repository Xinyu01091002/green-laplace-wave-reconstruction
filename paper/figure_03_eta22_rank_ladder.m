function output_path=figure_03_eta22_rank_ladder( ...
        input_root,output_root,output_name,apply_low_output_cutoff)
%FIGURE_03_ETA22_RANK_LADDER Field-level convergence figure.
arguments
    input_root (1,1) string = fullfile('artifacts', ...
        'eta22_field_rank_ladder_ce')
    output_root (1,1) string = fullfile('results','paper')
    output_name (1,1) string = "fig_eta22_field_rank_ladder_ce_no_output_cutoff"
    apply_low_output_cutoff (1,1) logical = false
end
if ~isfolder(output_root),mkdir(output_root);end
T=readtable(fullfile(input_root,'raw_records.csv'),'TextType','string');
run_summary=jsondecode(fileread(fullfile(input_root,'summary.json')));
if logical(run_summary.apply_low_output_cutoff)~=apply_low_output_cutoff
    error('Plot cutoff declaration does not match the stored MATLAB run.');
end
display_ranks=[2 4 6 8];
display_table=T(ismember(T.J,display_ranks),:);
depths=[1 2 4]; cases=["single_s30","balanced_spreading_a090_s30", ...
    "crossing_a120_s30"];
labels={'$0^\circ$ crossing sea, $\sigma_\theta=30^\circ$', ...
    '$90^\circ$ crossing sea, $\sigma_\theta=30^\circ$', ...
    '$120^\circ$ crossing sea, $\sigma_\theta=30^\circ$'};
colors=[0.15 0.15 0.15;0.0000 0.4470 0.7410;0.8500 0.3250 0.0980];
markers={'o','s','^'};
positive=display_table.raw_relative_l2_pct( ...
    display_table.raw_relative_l2_pct>0);
y_limits=10.^([floor(log10(min(positive)))-0.15, ...
    ceil(log10(max(positive)))+0.15]);

fig=figure('Color','w','Units','centimeters','Position',[2 2 18 7.6]);
positions=[0.07 0.17 0.275 0.67;0.385 0.17 0.275 0.67; ...
    0.70 0.17 0.275 0.67];
legend_handles=gobjects(1,numel(cases));
for depth_index=1:numel(depths)
    ax=axes(fig,'Position',positions(depth_index,:));
    set(ax,'YScale','log'); hold(ax,'on');
    for case_index=1:numel(cases)
        selected=display_table.kph==depths(depth_index) & ...
            display_table.case_id==cases(case_index);
        subset=sortrows(display_table(selected,:),"J");
        handle=semilogy(ax,subset.J,subset.raw_relative_l2_pct, ...
            ['-' markers{case_index}],'Color',colors(case_index,:), ...
            'MarkerFaceColor','w','MarkerSize',4.5,'LineWidth',1.15);
        if depth_index==1,legend_handles(case_index)=handle;end
    end
    grid(ax,'on'); box(ax,'on'); xlim(ax,[1.6 8.4]); ylim(ax,y_limits);
    ax.XMinorGrid='off'; ax.YMinorGrid='off';
    xticks(ax,display_ranks); xticklabels(ax,compose('GL%d',display_ranks));
    xlabel(ax,'Gauss--Laguerre order');
    if depth_index==1
        ylabel(ax,'Raw relative $L_2$ error (\%)');
    else
        ylabel(ax,'');
    end
    title(ax,sprintf('(%c) $k_p h=%g$',char('a'+depth_index-1), ...
        depths(depth_index)),'FontWeight','normal');
    set(ax,'FontName','CMU Serif','FontSize',8,'LineWidth',0.75, ...
        'TickDir','out','Layer','top','TickLabelInterpreter','latex');
    ax.XLabel.Interpreter='latex'; ax.YLabel.Interpreter='latex';
    ax.Title.Interpreter='latex';
end
legend(legend_handles,labels,'Orientation','horizontal','NumColumns',3, ...
    'Position',[0.18 0.91 0.64 0.05],'Box','off','Interpreter','latex', ...
    'FontName','CMU Serif','FontSize',8);

base=fullfile(output_root,output_name);
exportgraphics(fig,base+".png",'Resolution',600);
writetable(T,base+"_data.csv");
close(fig);
output_path=base+".png";
end
