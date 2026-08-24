%RUN_MINIMAL_EXAMPLE Unified eta/psi reconstruction through order three.
project_root = string(fileparts(fileparts(mfilename('fullpath'))));
addpath(project_root);
project_root = setup_green_laplace();

g = 9.81;
h = 1.0;
Lx = 2*pi;
Ly = 2*pi;
Nx = 64;
Ny = 64;
t = 0.0;
a = [0.010,0.006];
b = [0.000,0.002];
kx = [2,3];
ky = [0,1];

options = struct('eta22_rank',6);
coeffs = gl_spectral_coefficients( ...
    3,g,h,a,b,kx,ky,0,0,options);
[eta,psi,X,Y,components,audit] = ...
    gl_spectral_surface(coeffs,Lx,Ly,Nx,Ny,t);

fprintf('Green--Laplace order-three example\n');
fprintf('  eta L2 norm: %.8e\n',norm(eta(:)));
fprintf('  psi L2 norm: %.8e\n',norm(psi(:)));
fprintf('  components: %s\n',strjoin(audit.included_components,', '));

figure('Color','w','Position',[100,100,1100,420]);
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
nexttile;
imagesc(X(1,:),Y(:,1),components.eta33);
axis image xy; colorbar; title('\eta_{33}');
xlabel('x'); ylabel('y');
nexttile;
imagesc(X(1,:),Y(:,1),components.psi33);
axis image xy; colorbar; title('\psi_{33}');
xlabel('x'); ylabel('y');

output_root = fullfile(project_root,'results','minimal_example');
if ~isfolder(output_root), mkdir(output_root); end
exportgraphics(gcf,fullfile(output_root,'eta33_psi33.png'),'Resolution',200);
save(fullfile(output_root,'fields.mat'), ...
    'eta','psi','X','Y','components','audit');
