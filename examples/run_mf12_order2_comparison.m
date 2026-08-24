%RUN_MF12_ORDER2_COMPARISON Same-input pure-sum eta22 and psi22 comparison.
project_root = string(fileparts(fileparts(mfilename('fullpath'))));
addpath(project_root);
mf12_root = string(getenv('MF12_ROOT'));
if strlength(mf12_root)==0
    error(['Set MF12_ROOT to the external MF12 repository, for example: ' ...
        '$env:MF12_ROOT="C:\path\to\spectral-domain-implementation...' ...
        '"']);
end
project_root = setup_green_laplace("MF12Root",mf12_root);

g = 9.81;
h = 1.0;
Lx = 2*pi;
Ly = 2*pi;
Nx = 64;
Ny = 64;
t = 0.0;
coeffs = gl_spectral_coefficients(2,g,h, ...
    [0.010,0.006],[0.000,0.002],[2,3],[0,1],0,0, ...
    struct('eta22_rank',6));
comparison = gl_mf12_order2_comparison( ...
    coeffs,Lx,Ly,Nx,Ny,t);

fprintf('eta22 raw relative L2: %.6g\n', ...
    comparison.eta22_metrics.raw_relative_l2);
fprintf('psi22 raw relative L2: %.6g\n', ...
    comparison.psi22_metrics.raw_relative_l2);

figure('Color','w','Position',[100,100,1200,720]);
tiledlayout(2,3,'Padding','compact','TileSpacing','compact');
fields = {comparison.mf12_eta22,comparison.gl_eta22, ...
    comparison.gl_eta22-comparison.mf12_eta22, ...
    comparison.mf12_psi22,comparison.gl_psi22, ...
    comparison.gl_psi22-comparison.mf12_psi22};
titles = {'MF12 \eta_{22}','GL \eta_{22}','GL-MF12 \eta_{22}', ...
    'MF12 \psi_{22}','GL \psi_{22}','GL-MF12 \psi_{22}'};
for index = 1:6
    nexttile;
    imagesc(comparison.X(1,:),comparison.Y(:,1),fields{index});
    axis image xy; colorbar; title(titles{index});
end
output_root = fullfile(project_root,'results','mf12_order2_comparison');
if ~isfolder(output_root), mkdir(output_root); end
exportgraphics(gcf,fullfile(output_root,'eta22_psi22_comparison.png'), ...
    'Resolution',200);
save(fullfile(output_root,'comparison.mat'),'comparison');
