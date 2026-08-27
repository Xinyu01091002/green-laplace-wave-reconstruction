function prepared = eta20_prepare_real_input( ...
    eta11_spectrum, kx, ky, depth)
%ETA20_PREPARE_REAL_INPUT Eta11-only real-spectrum input audit.

arguments
    eta11_spectrum (:,:) {mustBeNumeric}
    kx (:,:) double
    ky (:,:) double
    depth (1,1) double {mustBePositive, mustBeFinite}
end
if ~isequal(size(eta11_spectrum),size(kx),size(ky))
    error('eta11_spectrum, kx, and ky must have identical sizes.');
end
if any(~isfinite(eta11_spectrum),'all') || ...
        any(~isfinite(kx),'all') || any(~isfinite(ky),'all')
    error('Neumann eta20 inputs must be finite.');
end
[ny,nx] = size(eta11_spectrum);
if min(nx,ny) < 16 || mod(nx,2) ~= 0 || mod(ny,2) ~= 0
    error('The frozen graph requires an even 2D grid of at least 16 by 16.');
end
partner_x = mod(-(0:nx-1),nx)+1;
partner_y = mod(-(0:ny-1),ny)+1;
mirror = conj(eta11_spectrum(partner_y,partner_x));
input_scale = max(1,norm(eta11_spectrum(:)));
hermitian_error = norm(eta11_spectrum(:)-mirror(:))/input_scale;
if hermitian_error > 1.0e-12
    error('eta11_spectrum is not the spectrum of a real field.');
end
if abs(eta11_spectrum(1,1)) > 1.0e-13*input_scale
    error('The first-order strict vector-zero mode must be zero.');
end
support_threshold = 1.0e-13*max(1,max(abs(eta11_spectrum),[],'all'));
support = abs(eta11_spectrum) > support_threshold;
if any(support & kx == 0,'all')
    error('Active parents on the kx=0 analytic boundary are excluded.');
end
positive = kx > 0;
if ~any(support & positive,'all')
    error('No strict-forward positive-frequency parent is active.');
end
u_spectrum = zeros(size(eta11_spectrum),'like',eta11_spectrum);
u_spectrum(positive) = eta11_spectrum(positive);

mode_x = [0:nx/2-1,-nx/2:-1];
mode_y = [0:ny/2-1,-ny/2:-1].';
[mode_x_grid,mode_y_grid] = meshgrid(mode_x,mode_y);
active = abs(u_spectrum) > support_threshold;
active_x = mode_x_grid(active);
active_y = mode_y_grid(active);
if max(active_x)-min(active_x) >= nx/2 || ...
        max(active_y)-min(active_y) >= ny/2
    error('Active parent difference support is not alias-safe.');
end
radial_k = hypot(kx,ky);
q = depth*radial_k;
energy = abs(u_spectrum).^2;
energy_sum = sum(energy,'all');
q0 = sum(q.*energy,'all')/energy_sum;
[~,peak_linear] = max(energy(:));
peak_q = q(peak_linear);
tail = positive & q < 0.3;
prepared = struct( ...
    'nx',nx,'ny',ny, ...
    'q',q,'q_output',q,'q0',q0,'peak_q',peak_q, ...
    'u_spectrum',u_spectrum, ...
    'support_threshold',support_threshold, ...
    'active_positive_count',nnz(active), ...
    'tail_q_below_point3_active_count',nnz(tail & active), ...
    'tail_q_below_point3_energy_fraction', ...
        sum(energy(tail),'all')/energy_sum, ...
    'hermitian_input_error',hermitian_error);
end
