function [eta,psi,X,Y,components,audit] = ...
    gl_spectral_surface(coeffs,Lx,Ly,Nx,Ny,t)
%GL_SPECTRAL_SURFACE Reconstruct surface elevation and surface potential.
%   The first four outputs match mf12_spectral_surface. Additional outputs
%   expose same-variable Green--Laplace components and execution metadata.

validate_coefficients(coeffs);
validateattributes(Lx,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(Ly,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(Nx,{'numeric'},{'scalar','integer','positive','even'});
validateattributes(Ny,{'numeric'},{'scalar','integer','positive','even'});
validateattributes(t,{'numeric'},{'scalar','real','finite'});

project_root = string(fileparts(fileparts(mfilename('fullpath'))));
dx = Lx/Nx;
dy = Ly/Ny;
[X,Y] = meshgrid((0:Nx-1)*dx,(0:Ny-1)*dy);
[kx_grid,ky_grid] = fft_wavenumber_grid(Lx,Ly,Nx,Ny);
eta11_spectrum = deposit_analytic_input( ...
    coeffs,kx_grid,ky_grid,Nx,Ny,t);
qx = coeffs.h*kx_grid;
qy = coeffs.h*ky_grid;
q = hypot(qx,qy);
nu = sqrt(q.*tanh(q));

eta11_plus = coeffs.h*ifft2(eta11_spectrum);
psi11_spectrum = complex(zeros(Ny,Nx));
nonzero = nu>0;
psi11_spectrum(nonzero) = -1i*eta11_spectrum(nonzero)./nu(nonzero);
potential_scale = coeffs.h*sqrt(coeffs.g*coeffs.h);
psi11_plus = potential_scale*ifft2(psi11_spectrum);

components = struct();
components.eta11_plus = eta11_plus;
components.psi11_plus = psi11_plus;
components.eta11 = real(eta11_plus);
components.psi11 = real(psi11_plus);
component_audits = struct();
eta_plus = eta11_plus;
psi_plus = psi11_plus;
included = ["eta11","psi11"];

if coeffs.order>=2
    [eta22_dimensionless,eta22_audit] = green_laplace_eta22( ...
        eta11_spectrum,qx,qy,coeffs.peak_depth, ...
        coeffs.options.eta22_rank,project_root);
    [psi22_dimensionless,psi22_audit] = ...
        finite_depth_directional_dual_branch_green_laplace_psi22( ...
        eta11_spectrum,qx,qy,coeffs.peak_depth,project_root, ...
        coeffs.options.allow_below_model_parent_domain,"radial2","none");
    eta22_plus = coeffs.h*eta22_dimensionless;
    psi22_plus = potential_scale*psi22_dimensionless;
    components.eta22_plus = eta22_plus;
    components.psi22_plus = psi22_plus;
    components.eta22 = real(eta22_plus);
    components.psi22 = real(psi22_plus);
    eta_plus = eta_plus+eta22_plus;
    psi_plus = psi_plus+psi22_plus;
    eta22_audit.candidate_id = 'gl-eta22-pure-prescribed-rank-v1';
    psi22_audit.candidate_id = 'gl-psi22-dual-branch-radial2-v1';
    component_audits.eta22 = eta22_audit;
    component_audits.psi22 = psi22_audit;
    included = [included,"eta22","psi22"];
end

if coeffs.order>=3
    [eta33_dimensionless,order3_audit,order3_parts, ...
        psi33_dimensionless] = finite_depth_directional_nested_gl_eta33_stable( ...
        eta11_spectrum,qx,qy,coeffs.peak_depth,project_root);
    eta33_plus = coeffs.h*eta33_dimensionless;
    psi33_plus = potential_scale*psi33_dimensionless;
    components.eta33_plus = eta33_plus;
    components.psi33_plus = psi33_plus;
    components.eta33 = real(eta33_plus);
    components.psi33 = real(psi33_plus);
    components.order3_internal = order3_parts;
    eta_plus = eta_plus+eta33_plus;
    psi_plus = psi_plus+psi33_plus;
    order3_audit.candidate_id = 'gl-eta33-psi33-stable-v1';
    component_audits.order3 = order3_audit;
    included = [included,"eta33","psi33"];
end

eta = real(eta_plus);
psi = real(psi_plus);
components.eta_plus = eta_plus;
components.psi_plus = psi_plus;
components.eta = eta;
components.psi = psi;

audit = struct( ...
    'api',coeffs.api, ...
    'order',coeffs.order, ...
    'sector',char(coeffs.options.sector), ...
    'included_components',{cellstr(included)}, ...
    'surface_potential_definition','psi=phi(x,z=eta,t)', ...
    'flat_potential_exposed',false, ...
    'rescaling_alignment_or_gain',false, ...
    'grid',[Ny,Nx], ...
    'domain',[Ly,Lx], ...
    'time',t, ...
    'component_audits',component_audits);
end

function validate_coefficients(coeffs)
required = {'api','order','g','h','a','b','kx','ky','omega', ...
    'peak_depth','options'};
if ~isstruct(coeffs) || ~all(isfield(coeffs,required)) ...
        || ~strcmp(coeffs.api,'green-laplace-spectral-v1')
    error('green_laplace:Coefficients', ...
        'Input must come from gl_spectral_coefficients.');
end
end

function [kx,ky] = fft_wavenumber_grid(Lx,Ly,Nx,Ny)
mode_x = [0:(Nx/2-1),-Nx/2:-1];
mode_y = [0:(Ny/2-1),-Ny/2:-1];
[mx,my] = meshgrid(mode_x,mode_y);
kx = (2*pi/Lx)*mx;
ky = (2*pi/Ly)*my;
end

function spectrum = deposit_analytic_input(coeffs,kx_grid,ky_grid,Nx,Ny,t)
delta_kx = kx_grid(1,2)-kx_grid(1,1);
delta_ky = ky_grid(2,1)-ky_grid(1,1);
spectrum = complex(zeros(Ny,Nx));
scale = Nx*Ny;
for index = 1:numel(coeffs.a)
    mode_x = round(coeffs.kx(index)/delta_kx);
    mode_y = round(coeffs.ky(index)/delta_ky);
    error_x = abs(coeffs.kx(index)-mode_x*delta_kx);
    error_y = abs(coeffs.ky(index)-mode_y*delta_ky);
    tolerance = coeffs.options.grid_tolerance*max(1,hypot( ...
        coeffs.kx(index),coeffs.ky(index)));
    if error_x>tolerance || error_y>tolerance
        error('green_laplace:OffGridComponent', ...
            'Every component must lie on the requested FFT grid.');
    end
    if abs(mode_x)>=Nx/2 || abs(mode_y)>=Ny/2
        error('green_laplace:Nyquist', ...
            'An input component reaches or exceeds a Nyquist boundary.');
    end
    column = mod(mode_x,Nx)+1;
    row = mod(mode_y,Ny)+1;
    amplitude = complex(coeffs.a(index),coeffs.b(index))/coeffs.h;
    value = scale*amplitude*exp(-1i*coeffs.omega(index)*t);
    spectrum(row,column) = spectrum(row,column)+value;
end
end
