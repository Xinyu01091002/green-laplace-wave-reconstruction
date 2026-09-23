function data=reproduce_order4_runtime()
%REPRODUCE_ORDER4_RUNTIME Plot archived measured times, without refitting.
root=setup_green_laplace();here=fileparts(mfilename('fullpath'));
out=fullfile(root,'results','order4');if ~isfolder(out),mkdir(out);end
data=readtable(fullfile(here,'data','runtime','eta44_runtime_figure_data.csv'), ...
    'TextType','string');
assert(all(data.cores==32));
fig=figure('Color','w','Position',[100 100 780 470]);hold on;
methods=["WIT","GL6","GL8","GL10"];colors=[0 0 0;lines(3)];
for k=1:numel(methods)
    rows=data(data.method==methods(k),:);rows=sortrows(rows,'component_count');
    plot(rows.component_count,rows.wall_seconds,'-o','Color',colors(k,:), ...
        'LineWidth',1.3,'DisplayName',methods(k));
end
set(gca,'YScale','log');grid on;xlabel('Number of retained linear components');
ylabel('Wall-clock time (s)');legend('Location','best');
title('Archived 32-core measurements: GL median of 3; WIT single runs');
exportgraphics(fig,fullfile(out,'eta44_gl_wit_runtime.png'),'Resolution',200);close(fig);
writetable(data,fullfile(out,'archived_runtime_rows.csv'));
end
