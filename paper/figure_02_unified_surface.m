function output = figure_02_unified_surface(output_directory)
%FIGURE_02_UNIFIED_SURFACE Paired eta/psi outputs from one API call.
arguments
    output_directory (1,1) string = string(fullfile( ...
        fileparts(fileparts(mfilename('fullpath'))),'results','paper'))
end
if ~isfolder(output_directory), mkdir(output_directory); end
setup_green_laplace;
c = gl_spectral_coefficients(3,9.81,1,[0.01,0.006],[0,0.002], ...
    [2,3],[0,1],0,0,struct('eta22_rank',6));
[~,~,X,Y,parts,audit] = gl_spectral_surface(c,2*pi,2*pi,64,64,0);

figure_handle = figure('Color','w','Position',[100,100,1100,760]);
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
fields = {parts.eta22,parts.psi22,parts.eta33,parts.psi33};
titles = {'(a) \eta_{22}','(b) \psi_{22}', ...
    '(c) \eta_{33}','(d) \psi_{33}'};
for index = 1:4
    nexttile;
    imagesc(X(1,:),Y(:,1),fields{index});
    axis image xy; colorbar; title(titles{index});
    xlabel('x'); ylabel('y');
end
png = fullfile(output_directory,'fig_unified_eta_psi_surface.png');
exportgraphics(figure_handle,png,'Resolution',300);
close(figure_handle);
mat = fullfile(output_directory,'fig_unified_eta_psi_surface.mat');
save(mat,'X','Y','parts','audit');
output = struct('png',png,'mat',mat,'audit',audit);
end
