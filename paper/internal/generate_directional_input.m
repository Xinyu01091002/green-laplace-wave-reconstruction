function state = generate_directional_input(config)
%GENERATE_DIRECTIONAL_INPUT Configurable dual-JONSWAP eta11 field.
%
% This is an input-spectrum generator, not a candidate or reference kernel.
% It follows the frozen spectral-cell-mass convention used by the official
% dual-JONSWAP validation, while exposing the discretization parameters.

arguments
    config (1,1) struct
end
required = {'retained_mass_percent','jonswap_gamma','spread_deg', ...
    'crossing_angle_deg','peak_kh','domain_wavelengths','grid_size', ...
    'selection_mode','maximum_components','radial_bins','angular_bins', ...
    'minimum_k_over_kp','maximum_k_over_kp','crest_amplitude', ...
    'time_peak_periods'};
for index = 1:numel(required)
    if ~isfield(config,required{index})
        error('Missing playground configuration field: %s',required{index});
    end
end
validateattributes(config.retained_mass_percent,{'numeric'}, ...
    {'scalar','>',0,'<=',100});
validateattributes(config.jonswap_gamma,{'numeric'},{'scalar','>=',1});
validateattributes(config.spread_deg,{'numeric'},{'scalar','positive'});
validateattributes(config.crossing_angle_deg,{'numeric'}, ...
    {'scalar','>=',0,'<',180});
validateattributes(config.peak_kh,{'numeric'},{'scalar','positive'});
validateattributes(config.domain_wavelengths,{'numeric'}, ...
    {'scalar','positive'});
validateattributes(config.grid_size,{'numeric'}, ...
    {'scalar','integer','>=',64});
if ~isfield(config,'auto_expand_grid')
    config.auto_expand_grid = false;
end
validateattributes(config.auto_expand_grid,{'logical','numeric'}, ...
    {'scalar'});
validateattributes(config.maximum_components,{'numeric'}, ...
    {'scalar','integer','positive'});
validateattributes(config.minimum_k_over_kp,{'numeric'}, ...
    {'scalar','positive'});
validateattributes(config.maximum_k_over_kp,{'numeric'}, ...
    {'scalar','>',config.minimum_k_over_kp});
validateattributes(config.crest_amplitude,{'numeric'},{'scalar','positive'});
validateattributes(config.time_peak_periods,{'numeric'},{'scalar','>=',0});

kp = 1.0;
gravity = 1.0;
depth = config.peak_kh/kp;
delta_k = kp/config.domain_wavelengths;
maximum_mode = ceil(config.maximum_k_over_kp*kp/delta_k);
requested_grid_size = config.grid_size;
minimum_alias_safe_grid = 4*maximum_mode+2;
if logical(config.auto_expand_grid)
    if mod(config.grid_size,2) ~= 0
        config.grid_size = config.grid_size+1;
    end
    if config.grid_size < minimum_alias_safe_grid
        config.grid_size = 2^nextpow2(minimum_alias_safe_grid);
    end
else
    if mod(config.grid_size,2) ~= 0
        error('FFT grid size must be even; received N=%d.',config.grid_size);
    end
    if config.grid_size < minimum_alias_safe_grid
        error([ ...
            'Difference support needs N >= %d for this domain and k/kp ' ...
            'tail, but N=%d. Enable automatic grid expansion or enter a ' ...
            'larger even grid.'],minimum_alias_safe_grid,config.grid_size);
    end
end

modes = -maximum_mode:maximum_mode;
[mode_x,mode_y] = meshgrid(modes,modes);
kx = mode_x*delta_k;
ky = mode_y*delta_k;
radial = hypot(kx,ky);
angle = atan2(ky,kx);
radial_density = jonswap_wavenumber_spectrum( ...
    radial,kp,gravity,config.jonswap_gamma);
directional_density = directional_mixture_density( ...
    angle,deg2rad(config.spread_deg),config.crossing_angle_deg);
cell_mass = radial_density./max(radial,realmin) ...
    .*directional_density*delta_k^2;
declared_envelope = radial/kp >= config.minimum_k_over_kp ...
    &radial/kp <= config.maximum_k_over_kp ...
    &isfinite(cell_mass) & cell_mass > 0;
admissible = kx > 0 & declared_envelope;
declared_envelope_mass = sum(cell_mass(declared_envelope),'all');
total_cell_mass = sum(cell_mass(admissible),'all');
if total_cell_mass <= 0
    error('No admissible positive-frequency JONSWAP cells were generated.');
end

if string(config.selection_mode) == "radial-angular"
    pool = aggregate_cells(mode_x,mode_y,kx,ky,radial,angle,cell_mass, ...
        admissible,config.radial_bins,config.angular_bins, ...
        config.minimum_k_over_kp*kp,config.maximum_k_over_kp*kp);
else
    indices = find(admissible);
    pool = struct( ...
        'mode_x',mode_x(indices), ...
        'mode_y',mode_y(indices), ...
        'kx',kx(indices), ...
        'ky',ky(indices), ...
        'radial',radial(indices), ...
        'angle',angle(indices), ...
        'cell_mass',cell_mass(indices));
end
[sorted_mass,order] = sort(pool.cell_mass(:),'descend');
target_mass = config.retained_mass_percent/100*total_cell_mass;
target_count = find(cumsum(sorted_mass) >= target_mass,1);
if isempty(target_count)
    target_count = numel(sorted_mass);
end
selected_count = min([target_count,config.maximum_components, ...
    numel(sorted_mass)]);
selected_order = order(1:selected_count);
selection = select_pool(pool,selected_order);
selected_mass = selection.cell_mass(:);

eta11_spectrum = complex(zeros(config.grid_size));
fft_scale = config.grid_size^2;
amplitude = 0.5*config.crest_amplitude ...
    *selected_mass/sum(selected_mass);
for component = 1:selected_count
    row = mod(selection.mode_y(component),config.grid_size)+1;
    column = mod(selection.mode_x(component),config.grid_size)+1;
    eta11_spectrum(row,column) = eta11_spectrum(row,column) ...
        +fft_scale*amplitude(component)/depth;
    mirror_row = mod(-selection.mode_y(component),config.grid_size)+1;
    mirror_column = mod(-selection.mode_x(component),config.grid_size)+1;
    eta11_spectrum(mirror_row,mirror_column) = ...
        conj(eta11_spectrum(row,column));
end
[kx_grid,ky_grid] = wave_number_grid(config.grid_size,delta_k);
peak_omega = sqrt(gravity*kp*tanh(depth*kp));
peak_period = 2*pi/peak_omega;
physical_time = config.time_peak_periods*peak_period;
omega = sqrt(gravity*hypot(kx_grid,ky_grid) ...
    .*tanh(depth*hypot(kx_grid,ky_grid)));
phase = ones(size(eta11_spectrum));
phase(kx_grid > 0) = exp(-1i*omega(kx_grid > 0)*physical_time);
phase(kx_grid < 0) = exp(1i*omega(kx_grid < 0)*physical_time);
eta11_spectrum = eta11_spectrum.*phase;
eta11 = real(ifft2(eta11_spectrum));

normalized_mass = selected_mass/sum(selected_mass);
mean_angle = sum(normalized_mass.*selection.angle(:));
combined_spread = sqrt(sum(normalized_mass ...
    .*(selection.angle(:)-mean_angle).^2));
state = struct( ...
    'config',config, ...
    'eta11_spectrum',eta11_spectrum, ...
    'eta11',eta11, ...
    'kx',kx_grid,'ky',ky_grid, ...
    'selection',selection, ...
    'delta_k',delta_k,'depth',depth, ...
    'requested_grid_size',requested_grid_size, ...
    'minimum_alias_safe_grid',minimum_alias_safe_grid, ...
    'grid_was_expanded',config.grid_size ~= requested_grid_size, ...
    'peak_period',peak_period,'physical_time',physical_time, ...
    'available_component_count',nnz(admissible), ...
    'admissible_cell_mass',total_cell_mass, ...
    'selected_component_count',selected_count, ...
    'declared_envelope_cell_mass',declared_envelope_mass, ...
    'forward_admissible_cell_mass',total_cell_mass, ...
    'forward_support_mass_fraction', ...
        total_cell_mass/declared_envelope_mass, ...
    'represented_declared_envelope_mass_fraction', ...
        sum(selected_mass)/declared_envelope_mass, ...
    'represented_cell_mass_fraction',sum(selected_mass)/total_cell_mass, ...
    'target_cell_mass_fraction',config.retained_mass_percent/100, ...
    'target_reached',sum(selected_mass) >= target_mass, ...
    'combined_realized_rms_spread_deg',rad2deg(combined_spread), ...
    'r2',[],'r2_seconds',NaN,'r2_audit',struct(), ...
    'shared4',[],'shared4_seconds',NaN,'shared4_audit',struct(), ...
    'endpoint4',[],'endpoint4_seconds',NaN,'endpoint4_audit',struct(), ...
    'endpoint6',[],'endpoint6_seconds',NaN,'endpoint6_audit',struct(), ...
    'mf12',[],'mf12_audit',struct(),'active_output',"");
end

function density = directional_mixture_density(angle,sigma,crossing_angle_deg)
if crossing_angle_deg == 0
    centers = 0.0;
else
    centers = deg2rad([-0.5,0.5]*crossing_angle_deg);
end
density = zeros(size(angle));
for center = centers
    density = density+exp(-0.5*((angle-center)/sigma).^2) ...
        /(sqrt(2*pi)*sigma*numel(centers));
end
end

function pool = aggregate_cells(mode_x,mode_y,kx,ky,radial,angle, ...
        cell_mass,admissible,radial_count,angular_count,minimum_k,maximum_k)
indices = find(admissible);
radial_edges = logspace(log10(minimum_k),log10(maximum_k),radial_count+1);
angular_edges = linspace(-pi/2,pi/2,angular_count+1);
radial_bin = discretize(radial(indices),radial_edges);
angular_bin = discretize(angle(indices),angular_edges);
maximum_bins = radial_count*angular_count;
values = zeros(maximum_bins,7);
count = 0;
for radial_index = 1:radial_count
    for angular_index = 1:angular_count
        members = find(radial_bin == radial_index ...
            &angular_bin == angular_index);
        if isempty(members), continue; end
        member_indices = indices(members);
        member_mass = cell_mass(member_indices);
        bin_mass = sum(member_mass);
        centroid_x = sum(member_mass.*kx(member_indices))/bin_mass;
        centroid_y = sum(member_mass.*ky(member_indices))/bin_mass;
        [~,offset] = min((kx(member_indices)-centroid_x).^2 ...
            +(ky(member_indices)-centroid_y).^2);
        representative = member_indices(offset);
        count = count+1;
        values(count,:) = [mode_x(representative),mode_y(representative), ...
            kx(representative),ky(representative),radial(representative), ...
            angle(representative),bin_mass];
    end
end
values = values(1:count,:);
pool = struct('mode_x',values(:,1),'mode_y',values(:,2), ...
    'kx',values(:,3),'ky',values(:,4),'radial',values(:,5), ...
    'angle',values(:,6),'cell_mass',values(:,7));
end

function selected = select_pool(pool,order)
names = fieldnames(pool);
selected = struct();
for index = 1:numel(names)
    selected.(names{index}) = pool.(names{index})(order);
end
end

function [kx,ky] = wave_number_grid(grid_size,delta_k)
indices = [0:grid_size/2-1,-grid_size/2:-1];
[kx,ky] = meshgrid(indices*delta_k,indices*delta_k);
end

function spectrum = jonswap_wavenumber_spectrum(k,kp,gravity,gamma)
k_safe = max(k,realmin);
frequency = sqrt(gravity*k_safe)/(2*pi);
peak_frequency = sqrt(gravity*kp)/(2*pi);
sigma = 0.09*ones(size(k_safe));
sigma(frequency <= peak_frequency) = 0.07;
peak_shape = exp(-0.5*((frequency-peak_frequency) ...
    ./(sigma*peak_frequency)).^2);
spectrum_f = frequency.^(-5) ...
    .*exp(-(5/4)*(peak_frequency./frequency).^4) ...
    .*gamma.^peak_shape;
df_dk = sqrt(gravity./k_safe)/(4*pi);
spectrum = spectrum_f.*df_dk;
spectrum(k <= 0) = 0;
end
