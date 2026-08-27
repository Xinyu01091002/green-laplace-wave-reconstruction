function result=figure_psi33_mixed_rank( ...
        field_path,output_root,figure_root)
%FIGURE_PSI33_MIXED_RANK Independent inner/outer pure-GL ranks.
arguments
    field_path (1,1) string
    output_root (1,1) string = string(fullfile('artifacts', ...
        'phi_s33_mixed_rank_diagnostic'))
    figure_root (1,1) string = string(fullfile('results','paper'))
end
if ~isfile(field_path),error('Matched surface-potential field is missing.');end
if ~isfolder(output_root),mkdir(output_root);end
if ~isfolder(figure_root),mkdir(figure_root);end
project_root=string(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(project_root,'paper','internal'));
loaded=load(field_path,'fields');source=loaded.fields;
n=size(source.eta11_spectrum,1);
depth=source.config.peak_kh/source.config.peak_wavenumber;
dk=source.config.delta_k_over_kp*source.config.peak_wavenumber;
axis_fft=[0:(n/2-1),-n/2:-1];
[mode_x,mode_y]=meshgrid(axis_fft,axis_fft);
qx=depth*dk*mode_x;qy=depth*dk*mode_y;
reference=2*real(source.reference_Psi33_plus);
reference_norm=norm(reference(:));reference_peak=max(abs(reference),[],'all');

% First slice: hold J3=4 and vary J2. Second slice: hold J2=10 and vary J3.
% The uniform (6,6) baseline and the similar-product (12,3) allocation are
% included explicitly. Formula, scales, and input remain unchanged.
rank_pairs=[4 4;6 4;8 4;10 4;12 4;6 6;8 6;10 6;12 6; ...
    8 8;10 8;12 8;10 10;12 3];
records=struct([]);
for index=1:size(rank_pairs,1)
    inner_rank=rank_pairs(index,1);outer_rank=rank_pairs(index,2);
    started=tic;
    [~,audit,components,~]=green_laplace_eta33_ranked( ...
        source.eta11_spectrum,qx,qy,source.config.peak_kh, ...
        inner_rank,project_root,outer_rank);
    elapsed=toc(started);
    candidate=2*real(ifft2(components.Psi33_raw_resolvent_spectrum));
    difference=candidate-reference;
    candidate_norm=norm(candidate(:));difference_norm=norm(difference(:));
    record=struct( ...
        'inner_rank_J2',inner_rank,'outer_rank_J3',outer_rank, ...
        'surface_similarity_Q',difference_norm/(reference_norm+candidate_norm), ...
        'physical_raw_relative_l2_pct',100*difference_norm/reference_norm, ...
        'physical_raw_relative_linf_pct', ...
            100*max(abs(difference),[],'all')/reference_peak, ...
        'gl_over_mf12_norm_ratio',candidate_norm/reference_norm, ...
        'diagonal_correction_applied',false, ...
        'fft_ifft_count',audit.fft_ifft_count, ...
        'pointwise_product_count',audit.pointwise_product_count, ...
        'elapsed_seconds',elapsed,'pareto_nondominated',false);
    if index==1,records=record;else,records(index)=record;end
    fprintf('MIXED_RANK J2=%d J3=%d Q=%.9g FFT=%d elapsed=%.3fs\n', ...
        inner_rank,outer_rank,record.surface_similarity_Q, ...
        record.fft_ifft_count,elapsed);
end

for index=1:numel(records)
    dominated=false;
    for other=1:numel(records)
        if other==index,continue,end
        no_worse=records(other).fft_ifft_count<=records(index).fft_ifft_count ...
            && records(other).surface_similarity_Q<=records(index).surface_similarity_Q;
        strict=records(other).fft_ifft_count<records(index).fft_ifft_count ...
            || records(other).surface_similarity_Q<records(index).surface_similarity_Q;
        dominated=dominated||(no_worse&&strict);
    end
    records(index).pareto_nondominated=~dominated;
end
result_table=struct2table(records);
writetable(result_table,fullfile(output_root,'mixed_rank_metrics.csv'));
writetable(result_table,fullfile(figure_root, ...
    'fig_phi_s33_mixed_rank_cost_accuracy_data.csv'));

fig=figure('Color','w','Units','centimeters','Position',[2 2 12 8]);
ax=axes(fig);hold(ax,'on');grid(ax,'on');box(ax,'on');
is_pareto=result_table.pareto_nondominated;
scatter(ax,result_table.fft_ifft_count(~is_pareto), ...
    result_table.surface_similarity_Q(~is_pareto),48,[0.55 0.55 0.55],'filled');
scatter(ax,result_table.fft_ifft_count(is_pareto), ...
    result_table.surface_similarity_Q(is_pareto),58,[0.00 0.45 0.74],'filled');
for index=1:height(result_table)
    text(ax,result_table.fft_ifft_count(index)+12, ...
        result_table.surface_similarity_Q(index), ...
        sprintf('(%d,%d)',result_table.inner_rank_J2(index), ...
        result_table.outer_rank_J3(index)), ...
        'FontName','CMU Serif','FontSize',8,'Interpreter','latex');
end
set(ax,'YScale','log','FontName','CMU Serif','FontSize',9, ...
    'TickLabelInterpreter','latex','TickDir','out');
xlabel(ax,'FFT/IFFT count','Interpreter','latex');
ylabel(ax,'Full-field $Q$','Interpreter','latex');
title(ax,'Pure GL inner--outer rank allocation $(J_2,J_3)$', ...
    'Interpreter','latex','FontWeight','normal');
legend(ax,{'Dominated','Pareto nondominated'},'Interpreter','latex', ...
    'Location','northeast','Box','off');
figure_path=fullfile(figure_root,'fig_phi_s33_mixed_rank_cost_accuracy.png');
exportgraphics(fig,figure_path,'Resolution',600);close(fig);

result=struct('status','post-freeze-nonselecting-pure-gl-mixed-rank', ...
    'input_case','kph-1p0-spread-30-cross-120', ...
    'metric','full-field symmetric Q','diagonal_correction_applied',false, ...
    'records',records,'figure_path',figure_path,'overall_pass',true);
file_id=fopen(fullfile(output_root,'summary.json'),'w');
if file_id<0,error('Could not open mixed-rank summary.');end
cleanup=onCleanup(@()fclose(file_id));
fwrite(file_id,jsonencode(result,'PrettyPrint',true),'char');
fprintf('PHI_S33_MIXED_RANK_DIAGNOSTIC_COMPLETE %s\n',figure_path);
end
