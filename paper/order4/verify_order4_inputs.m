function verify_order4_inputs()
root=setup_green_laplace();here=fileparts(mfilename('fullpath'));addpath(here);
cases=dir(fullfile(here,'inputs','runtime','n*'));
for k=1:numel(cases)
    prepare_order4_inputs(fullfile(cases(k).folder,cases(k).name), ...
        fullfile(root,'artifacts','order4-input-check',cases(k).name));
end
prepare_order4_inputs(fullfile(here,'inputs','field_n0382'), ...
    fullfile(root,'artifacts','order4-input-check','field_n0382'));
rows=readtable(fullfile(here,'data','runtime','eta44_runtime_figure_data.csv'),'TextType','string');
raw=[readtable(fullfile(here,'data','runtime','timings.csv')); ...
    readtable(fullfile(here,'data','runtime','timings_n1500.csv'))];
for k=1:height(rows)
    if rows.method(k)=="WIT",continue,end
    selected=raw.component_count==rows.component_count(k) & raw.rank==rows.quadrature_rank(k);
    assert(sum(selected)==3,'Each GL point must retain three measured runs.');
    assert(abs(median(raw.wall_seconds(selected))-rows.wall_seconds(k))<1e-8);
end
fprintf('PASS: %d frozen input cases and all GL timing medians\n',numel(cases)+1);
end
