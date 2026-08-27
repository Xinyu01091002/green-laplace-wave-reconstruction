function [eta20, audit, state] = ...
    eta20_green_laplace_shared( ...
        eta11_spectrum, kx, ky, depth, representation, project_root, ...
        quadrature_scale_override)
%ETA20_GREEN_LAPLACE_SHARED Frozen eta11-only GL diagnostic graphs.
%
% representation is shared GL2/4/6/8/12/16, endpoint-subtracted GL4/6/8, or the
% exploratory "two-scale" branch graph.

arguments
    eta11_spectrum (:,:) {mustBeNumeric}
    kx (:,:) double
    ky (:,:) double
    depth (1,1) double {mustBePositive, mustBeFinite}
    representation (1,1) string {mustBeMember( ...
        representation,["shared","shared2","shared4","shared6", ...
        "shared8","shared12","shared16","endpoint4","endpoint6","endpoint8", ...
        "endpoint8r2","endpoint8r4","endpoint8h4","endpoint8h42", ...
        "endpoint8m2h42","endpointshell864m2", ...
        "endpointls10r8","endpointc44","two-scale"])} = "shared"
    project_root (1,1) string = ...
        string(fileparts(fileparts(fileparts(mfilename('fullpath')))))
    quadrature_scale_override (1,1) double = NaN
end

prepared = eta20_prepare_real_input( ...
    eta11_spectrum,kx,ky,depth);
config = finite_depth_eta20_two_scale_gl_configuration();
[shared_nodes,shared_weights,shared_cost,shared_freeze,endpoint_coefficient, ...
    endpoint_copies,endpoint_terms] = ...
    shared_rule(representation,config);
q = prepared.q;
nu = sqrt(q.*tanh(q));
dno = q.*tanh(q);
parent_nonzero = q > 0;
unit_x = zeros(size(q));
unit_y = zeros(size(q));
unit_x(parent_nonzero) = depth*kx(parent_nonzero)./q(parent_nonzero);
unit_y(parent_nonzero) = depth*ky(parent_nonzero)./q(parent_nonzero);
coth_q = zeros(size(q));
coth_q(parent_nonzero) = 1./tanh(q(parent_nonzero));
base = struct( ...
    'one',ones(size(q)), ...
    'dno',dno, ...
    'nu',nu, ...
    'bx',unit_x.*coth_q.*nu, ...
    'by',unit_y.*coth_q.*nu, ...
    'qx',q.*unit_x, ...
    'qy',q.*unit_y, ...
    'a',q.*coth_q.*nu);

energy = abs(prepared.u_spectrum).^2;
energy_sum = sum(energy,'all');
qx_dimensionless = depth*kx;
qy_dimensionless = depth*ky;
mean_qx = sum(qx_dimensionless.*energy,'all')/energy_sum;
mean_qy = sum(qy_dimensionless.*energy,'all')/energy_sum;
vector_variance = sum(( ...
    (qx_dimensionless-mean_qx).^2 ...
    +(qy_dimensionless-mean_qy).^2).*energy,'all')/energy_sum;
delta_q = sqrt(2*vector_variance);
quadrature_scale = delta_q;
scale_override_active = isfinite(quadrature_scale_override);
if scale_override_active
    if quadrature_scale_override <= 0
        error('Diagnostic quadrature scale override must be positive.');
    end
    if representation == "two-scale" || representation == "endpointshell864m2"
        error('Diagnostic scale override is unavailable for this representation.');
    end
    quadrature_scale = quadrature_scale_override;
end
active = abs(prepared.u_spectrum) > prepared.support_threshold;
projection = qx_dimensionless;
active_nu = nu(active);
active_projection = projection(active);
if any(active_projection <= 0)
    error('The exponential balance requires strict-forward parents.');
end
% Exact difference-frequency projection balance.  For every contributing
% pair, p1-p2=P, so shifting nu_j by v*p_j in the parent exponentials and
% restoring exp(sign*v*P*tau) at the output leaves the frozen GL formula
% unchanged.  The half-minimum phase velocity also keeps both output
% restoration exponents non-growing on the declared strict-forward support.
balance_velocity = 0.5*min(active_nu./active_projection);
radial_balance_channel_count = 1;
if representation == "endpoint8r2"
    radial_balance_channel_count = 2;
elseif representation == "endpoint8r4"
    radial_balance_channel_count = 4;
elseif representation == "endpointls10r8"
    radial_balance_channel_count = 8;
end
[balance_masks,balance_centers,balance_pair_supports] = ...
    difference_balance_channels( ...
        active,nu-balance_velocity*projection,radial_balance_channel_count);
balance_center = mean(balance_centers);
support_mask_transforms = radial_balance_channel_count ...
    +radial_balance_channel_count^2;
terminal_balance_channel_count = radial_balance_channel_count;
terminal_balance_masks = balance_masks;
terminal_balance_centers = balance_centers;
terminal_balance_pair_supports = balance_pair_supports;
terminal_radial4_node_count = 0;
if representation == "endpoint8h4"
    terminal_radial4_node_count = 1;
elseif representation == "endpoint8h42" || representation == "endpoint8m2h42"
    terminal_radial4_node_count = 2;
end
terminal_radial4_active = terminal_radial4_node_count > 0;
if terminal_radial4_active
    terminal_balance_channel_count = 4;
    [terminal_balance_masks,terminal_balance_centers, ...
        terminal_balance_pair_supports] = difference_balance_channels( ...
            active,nu-balance_velocity*projection, ...
            terminal_balance_channel_count);
    support_mask_transforms = support_mask_transforms ...
        +terminal_balance_channel_count+terminal_balance_channel_count^2;
end
q0 = prepared.q0;
nu0 = sqrt(q0*tanh(q0));
group_velocity = (tanh(q0)+q0*sech(q0)^2)/(2*nu0);

eta20_spectrum = complex(zeros(size(eta11_spectrum)));
degenerate_threshold = 64*eps(max(1,q0));
shell_active = representation == "endpointshell864m2";
if shell_active && delta_q > degenerate_threshold
    [eta20_spectrum,shell_state] = shell_accumulation( ...
        base,nu,q,projection,prepared.u_spectrum,active,balance_velocity, ...
        delta_q,group_velocity,endpoint_terms,endpoint_copies);
    expected_transforms = shared_cost.transform_count;
    expected_products = shared_cost.product_count;
    expected_parent_iffts = shared_cost.parent_iffts;
    expected_product_ffts = shared_cost.product_ffts;
    support_mask_transforms = shell_state.support_mask_transforms;
    candidate_id = ...
        'gl-eta20-multiterm-endpoint-shell864-balanced-v1';
    hermitian_projection = true;
elseif delta_q > degenerate_threshold
    if representation ~= "two-scale"
        scale = quadrature_scale;
        for node_index = 1:numel(shared_nodes)
            x = shared_nodes(node_index);
            measure_weight = shared_weights(node_index);
            tau = x/scale;
            node_masks = balance_masks;
            node_centers = balance_centers;
            node_pair_supports = balance_pair_supports;
            if terminal_radial4_active ...
                    && node_index > numel(shared_nodes)-terminal_radial4_node_count
                node_masks = terminal_balance_masks;
                node_centers = terminal_balance_centers;
                node_pair_supports = terminal_balance_pair_supports;
            end
            suppress_node_roundoff = (representation=="shared16" ...
                && node_index>15) || scale_override_active;
            [fd_plus,fk_plus] = forcing_at_time( ...
                base,nu,projection,prepared.u_spectrum,node_masks, ...
                node_centers,node_pair_supports, ...
                balance_velocity,tau,+1,suppress_node_roundoff);
            [fd_minus,fk_minus] = forcing_at_time( ...
                base,nu,projection,prepared.u_spectrum,node_masks, ...
                node_centers,node_pair_supports, ...
                balance_velocity,tau,-1,suppress_node_roundoff);
            eta20_spectrum = eta20_spectrum ...
                +(measure_weight/(4*scale))*( ...
                    1i*fk_plus-nu.*fd_plus ...
                    -1i*fk_minus-nu.*fd_minus);
        end
        expected_transforms = shared_cost.transform_count;
        expected_products = shared_cost.product_count;
        expected_parent_iffts = shared_cost.parent_iffts;
        expected_product_ffts = shared_cost.product_ffts;
        candidate_id = sprintf( ...
            'gl-eta20-shared-scale-n%d-v1',numel(shared_nodes));
        hermitian_projection = false;
    else
        slow_scale = delta_q*(1-group_velocity);
        fast_scale = delta_q*(1+group_velocity);
        if slow_scale <= 0 || fast_scale <= 0
            error('Frozen Two-Scale eta20 branch scales are not positive.');
        end
        for node_index = 1:numel(config.slow_nodes)
            x = config.slow_nodes(node_index);
            weight = config.slow_weights(node_index);
            tau = x/slow_scale;
            [fd_plus,fk_plus] = forcing_at_time( ...
                base,nu,projection,prepared.u_spectrum,balance_masks, ...
                balance_centers,balance_pair_supports, ...
                balance_velocity,tau,+1);
            eta20_spectrum = eta20_spectrum ...
                +(weight*exp(x)/(4*slow_scale))*( ...
                    1i*fk_plus-nu.*fd_plus);
        end
        for node_index = 1:numel(config.fast_nodes)
            x = config.fast_nodes(node_index);
            weight = config.fast_weights(node_index);
            tau = x/fast_scale;
            [fd_minus,fk_minus] = forcing_at_time( ...
                base,nu,projection,prepared.u_spectrum,balance_masks, ...
                balance_centers,balance_pair_supports, ...
                balance_velocity,tau,-1);
            eta20_spectrum = eta20_spectrum ...
                +(weight*exp(x)/(4*fast_scale))*( ...
                    -1i*fk_minus-nu.*fd_minus);
        end
        eta20_spectrum = hermitian_project(eta20_spectrum);
        expected_transforms = 57;
        expected_products = 33;
        expected_parent_iffts = 48;
        expected_product_ffts = 6;
        candidate_id = 'gl-eta20-two-scale-s2-f1-hp-v1';
        hermitian_projection = true;
    end
else
    if representation ~= "two-scale"
        expected_transforms = shared_cost.transform_count;
        expected_products = shared_cost.product_count;
        expected_parent_iffts = shared_cost.parent_iffts;
        expected_product_ffts = shared_cost.product_ffts;
        candidate_id = sprintf( ...
            'gl-eta20-shared-scale-n%d-v1',numel(shared_nodes));
        hermitian_projection = false;
    else
        expected_transforms = 57;
        expected_products = 33;
        expected_parent_iffts = 48;
        expected_product_ffts = 6;
        candidate_id = 'gl-eta20-two-scale-s2-f1-hp-v1';
        hermitian_projection = true;
    end
end

if endpoint_coefficient ~= 0 && ~shell_active ...
        && delta_q > degenerate_threshold
    eta20_spectrum = eta20_spectrum+endpoint_correction( ...
        nu,q,projection,prepared.u_spectrum,balance_masks, ...
        balance_centers,balance_pair_supports,terminal_balance_masks, ...
        terminal_balance_centers,terminal_balance_pair_supports, ...
        terminal_radial4_node_count,balance_velocity,quadrature_scale, ...
        shared_nodes,shared_weights,endpoint_terms,endpoint_copies, ...
        [-1,1],true);
end

% Exact shared GL fields are Hermitian.  Higher nodes amplify roundoff, so
% enforce that exact symmetry as an O(N) numerical guard (no transform).
if (startsWith(representation,"shared") ...
        || startsWith(representation,"endpoint")) && numel(shared_nodes) > 2
    eta20_spectrum = hermitian_project(eta20_spectrum);
    hermitian_projection = true;
    if endpoint_coefficient ~= 0
        if representation == "endpointc44"
            candidate_id = ...
                'gl-eta20-endpoint-composite-leg4-lag4-hp-v1';
        elseif representation == "endpoint8r2" || representation == "endpoint8r4"
            candidate_id = sprintf( ...
                'gl-eta20-endpoint-subtracted-n8-radial%d-hp-v1', ...
                radial_balance_channel_count);
        elseif representation == "endpoint8h4" || representation == "endpoint8h42"
            candidate_id = sprintf( ...
                'gl-eta20-endpoint-gl8-terminal%d-radial4-hp-v1', ...
                terminal_radial4_node_count);
        elseif representation == "endpoint8m2h42"
            candidate_id = ...
                'gl-eta20-multiterm-endpoint-gl8-terminal2-radial4-v1';
        elseif representation == "endpointls10r8"
            candidate_id = ...
                'gl-eta20-endpoint-log-sinc10-radial8-hp-v1';
        else
            candidate_id = sprintf( ...
                'gl-eta20-endpoint-subtracted-n%d-hp-v1', ...
                numel(shared_nodes));
        end
    else
        candidate_id = sprintf( ...
            'gl-eta20-shared-scale-n%d-hp-v1', ...
            numel(shared_nodes));
    end
end
if scale_override_active
    candidate_id = [candidate_id '-scale-override-diagnostic'];
end

eta20_spectrum(q == 0) = 0;
eta20_complex = ifft2(eta20_spectrum/depth);
eta20 = real(eta20_complex);
imaginary_leakage = norm(imag(eta20_complex(:))) ...
    /max(1,norm(eta20(:)));
if imaginary_leakage > 1.0e-8
    error('Two-Scale GL eta20 output is not Hermitian (leakage %.17g).', ...
        imaginary_leakage);
end
if any(~isfinite(eta20),'all') || any(~isfinite(eta20_spectrum),'all')
    error('Two-Scale GL eta20 output is not finite.');
end

audit = struct( ...
    'candidate_id',candidate_id, ...
    'project_root',project_root, ...
    'formula_provenance','wolfram-exact-outer-resolvent-laplace-freeze', ...
    'freeze_sha256',shared_freeze, ...
    'feasibility_sha256',config.feasibility_sha256, ...
    'information_boundary','eta11-only', ...
    'strict_forward_parent_rule','kx-positive-analytic-spectrum', ...
    'strict_zero_mode','set-to-zero-separate-sector', ...
    'opposition_scope','strict-forward parents only; exact opposition excluded', ...
    'carrier_q0',q0, ...
    'peak_kh',prepared.peak_q, ...
    'delta_q_rms_pair',delta_q, ...
    'quadrature_scale',quadrature_scale, ...
    'quadrature_scale_override_active',scale_override_active, ...
    'group_velocity_q0',group_velocity, ...
    'slow_endpoint_subtraction',endpoint_coefficient ~= 0, ...
    'slow_endpoint_coefficient',endpoint_coefficient, ...
    'endpoint_separable_term_count',size(endpoint_terms,1), ...
    'endpoint_copies_added_exactly',endpoint_copies, ...
    'exponential_balance_mode','exact-half-forward-projection-plus-center', ...
    'balance_velocity',balance_velocity, ...
    'balance_center_nu',balance_center, ...
    'radial_balance_channel_count',radial_balance_channel_count, ...
    'radial_balance_centers',balance_centers, ...
    'terminal_balance_channel_count',terminal_balance_channel_count, ...
    'terminal_balance_centers',terminal_balance_centers, ...
    'terminal_radial4_active',terminal_radial4_active, ...
    'terminal_radial4_node_count',terminal_radial4_node_count, ...
    'degenerate_zero_output',delta_q <= degenerate_threshold, ...
    'tail_kh_below_point3_active_count', ...
        prepared.tail_q_below_point3_active_count, ...
    'tail_kh_below_point3_energy_fraction', ...
        prepared.tail_q_below_point3_energy_fraction, ...
    'active_positive_parent_count',prepared.active_positive_count, ...
    'spectrum_input_fft_or_ifft',expected_transforms, ...
    'unique_parent_filtered_ifft',expected_parent_iffts, ...
    'product_accumulator_fft',expected_product_ffts, ...
    'final_output_ifft',1, ...
    'complex_physical_products',expected_products, ...
    'support_mask_transforms',support_mask_transforms, ...
    'support_mask_rule','exact channel-pair S_i-S_j support; explicit q_output=0 removal', ...
    'scalar_spectral_reductions',9, ...
    'pair_loops',0, ...
    'hermitian_projection',hermitian_projection, ...
    'asymptotic_work','constant-times N log N', ...
    'hermitian_input_error',prepared.hermitian_input_error, ...
    'imaginary_output_leakage',imaginary_leakage);
state = struct( ...
    'u_spectrum',prepared.u_spectrum, ...
    'eta20_spectrum',eta20_spectrum/depth, ...
    'delta_q_rms_pair',delta_q, ...
    'quadrature_scale',quadrature_scale, ...
    'quadrature_scale_override_active',scale_override_active, ...
    'group_velocity_q0',group_velocity, ...
    'balance_velocity',balance_velocity, ...
    'balance_center_nu',balance_center);
end

function [nodes,weights,cost,freeze_sha256,endpoint_coefficient, ...
        endpoint_copies,endpoint_terms] = ...
        shared_rule(representation,base)
endpoint_coefficient = 0;
endpoint_copies = 0;
endpoint_terms = zeros(0,3);
if representation == "two-scale" || representation == "shared" ...
        || representation == "shared2"
    nodes = base.shared_nodes;
    weights = base.shared_weights.*exp(base.shared_nodes);
    node_count = numel(nodes);
    cost = struct('transform_count',36*node_count+3, ...
        'product_count',22*node_count, ...
        'parent_iffts',32*node_count,'product_ffts',4*node_count);
    freeze_sha256 = base.freeze_sha256;
    return
end
if representation == "endpointc44"
    levels = finite_depth_eta20_endpoint_composite44_configuration();
    endpoint_coefficient = levels.endpoint_coefficient;
    endpoint_copies = levels.endpoint_copies;
    endpoint_terms = [-1,-1,endpoint_coefficient];
    nodes = levels.scaled_times;
    weights = levels.scaled_measure_weights;
    node_count = numel(nodes);
    cost = struct('transform_count',48*node_count+5, ...
        'product_count',28*node_count+1, ...
        'parent_iffts',40*node_count+1,'product_ffts',8*node_count+1);
    freeze_sha256 = levels.freeze_sha256;
    return
end
if representation == "endpointls10r8"
    alternative = finite_depth_eta20_alternative_nlogn_configuration();
    endpoint = finite_depth_eta20_endpoint_subtracted_gl_configuration();
    endpoint_coefficient = endpoint.endpoint_coefficient;
    endpoint_copies = endpoint.endpoint_copies;
    endpoint_terms = [-1,-1,endpoint_coefficient];
    nodes = alternative.log_sinc_nodes;
    weights = alternative.log_sinc_weights;
    channel_count = alternative.log_sinc_radial_balance_channels;
    node_count = numel(nodes);
    parent_iffts = 40*channel_count*node_count+1;
    product_ffts = 8*channel_count^2*node_count+1;
    support_transforms = channel_count+channel_count^2;
    cost = struct( ...
        'transform_count',parent_iffts+product_ffts ...
            +support_transforms+1, ...
        'product_count',28*channel_count^2*node_count+1, ...
        'parent_iffts',parent_iffts,'product_ffts',product_ffts);
    freeze_sha256 = alternative.freeze_sha256;
    return
end
if representation == "endpoint8m2h42" ...
        || representation == "endpointshell864m2"
    multi = finite_depth_eta20_multiterm_shell_configuration();
    endpoint_terms = [multi.endpoint_first_exponents(:), ...
        multi.endpoint_second_exponents(:),multi.endpoint_coefficients(:)];
    endpoint_coefficient = endpoint_terms(1,3);
    endpoint_copies = 2;
    freeze_sha256 = multi.freeze_sha256;
    if representation == "endpoint8m2h42"
        endpoint = finite_depth_eta20_endpoint_subtracted_gl_configuration();
        nodes = endpoint.gl8_nodes;
        weights = endpoint.gl8_weights.*exp(endpoint.gl8_nodes);
        cost = struct('transform_count',multi.multiterm_h42_transform_count, ...
            'product_count',multi.multiterm_h42_product_count, ...
            'parent_iffts',1013,'product_ffts',305);
    else
        nodes = zeros(1,0);
        weights = zeros(1,0);
        cost = struct('transform_count',multi.shell864_transform_count, ...
            'product_count',multi.shell864_product_count, ...
            'parent_iffts',5981,'product_ffts',6233);
    end
    return
end
if representation == "endpoint8r2" || representation == "endpoint8r4" ...
        || representation == "endpoint8h4" || representation == "endpoint8h42"
    order = 8;
else
    order = str2double(extract(representation,digitsPattern));
end
if startsWith(representation,"endpoint")
    levels = finite_depth_eta20_endpoint_subtracted_gl_configuration();
    endpoint_coefficient = levels.endpoint_coefficient;
    endpoint_copies = levels.endpoint_copies;
    endpoint_terms = [-1,-1,endpoint_coefficient];
    nodes = levels.(sprintf('gl%d_nodes',order));
    weights = levels.(sprintf('gl%d_weights',order)) ...
        .*exp(levels.(sprintf('gl%d_nodes',order)));
    if representation == "endpoint8h4" || representation == "endpoint8h42"
        if representation == "endpoint8h4"
            parent_iffts = 441;
            product_ffts = 185;
            transform_count = 649;
            product_count = 645;
        else
            parent_iffts = 561;
            product_ffts = 305;
            transform_count = 889;
            product_count = 1065;
        end
        cost = struct( ...
            'transform_count',transform_count, ...
            'product_count',product_count, ...
            'parent_iffts',parent_iffts,'product_ffts',product_ffts);
    elseif representation == "endpoint8r2" || representation == "endpoint8r4"
        if representation == "endpoint8r2"
            channel_count = 2;
        else
            channel_count = 4;
        end
        parent_iffts = 40*channel_count*order+1;
        product_ffts = 8*channel_count^2*order+1;
        support_transforms = channel_count+channel_count^2;
        cost = struct( ...
            'transform_count',parent_iffts+product_ffts ...
                +support_transforms+1, ...
            'product_count',28*channel_count^2*order+1, ...
            'parent_iffts',parent_iffts,'product_ffts',product_ffts);
    else
        cost = struct( ...
            'transform_count',48*order+5, ...
            'product_count',28*order+1, ...
            'parent_iffts',40*order+1,'product_ffts',8*order+1);
    end
    freeze_sha256 = levels.freeze_sha256;
    return
end
if order==12
    levels = finite_depth_eta20_shared_gl12_diagnostic_configuration();
elseif order==16
    levels = finite_depth_eta20_shared_gl16_diagnostic_configuration();
else
    levels = finite_depth_eta20_shared_gl_levels_configuration();
end
nodes = levels.(sprintf('gl%d_nodes',order));
weights = levels.(sprintf('gl%d_weights',order)) ...
    .*exp(levels.(sprintf('gl%d_nodes',order)));
cost = struct( ...
    'transform_count',36*order+3, ...
    'product_count',22*order, ...
    'parent_iffts',32*order,'product_ffts',4*order);
freeze_sha256 = levels.freeze_sha256;
end

function correction = endpoint_correction( ...
        nu,q,projection,u_spectrum,masks,centers,pair_supports, ...
        terminal_masks,terminal_centers,terminal_pair_supports, ...
        terminal_radial4_node_count,balance_velocity,scale,nodes,weights, ...
        terms,endpoint_copies,branch_signs,include_explicit)
if size(terms,1)==1 && terms(1,1)==-1 && terms(1,2)==-1
    correction = endpoint_correction_single(nu,q,projection,u_spectrum, ...
        masks,centers,pair_supports,terminal_masks,terminal_centers, ...
        terminal_pair_supports,terminal_radial4_node_count,balance_velocity, ...
        scale,nodes,weights,terms(1,3),endpoint_copies,branch_signs, ...
        include_explicit);
    return
end
powers = unique([terms(:,1);terms(:,2)]).';
power_multipliers = cell(1,numel(powers));
for power_index=1:numel(powers)
    power_multipliers{power_index} = radial_power(q,powers(power_index));
end
if include_explicit
    endpoint_banks = cell(1,numel(powers));
    for power_index=1:numel(powers)
        endpoint_banks{power_index} = ifft2( ...
            power_multipliers{power_index}.*u_spectrum);
    end
    explicit_field = complex(zeros(size(q)));
    for term_index=1:size(terms,1)
        first_power = find(powers==terms(term_index,1),1);
        second_power = find(powers==terms(term_index,2),1);
        explicit_field = explicit_field+terms(term_index,3) ...
            *endpoint_banks{first_power}.*conj(endpoint_banks{second_power});
    end
    correction = endpoint_copies*fft2(explicit_field);
else
    correction = complex(zeros(size(q)));
end
for node_index = 1:numel(nodes)
    x = nodes(node_index);
    tau = x/scale;
    node_masks = masks;
    node_centers = centers;
    node_pair_supports = pair_supports;
    if terminal_radial4_node_count > 0 ...
            && node_index > numel(nodes)-terminal_radial4_node_count
        node_masks = terminal_masks;
        node_centers = terminal_centers;
        node_pair_supports = terminal_pair_supports;
    end
    channel_count = numel(node_masks);
    for branch_sign = branch_signs
        left = cell(channel_count,numel(powers));
        right = cell(channel_count,numel(powers));
        left_nu = cell(channel_count,numel(powers));
        right_nu = cell(channel_count,numel(powers));
        for channel = 1:channel_count
            mask = node_masks{channel};
            residual = nu-balance_velocity*projection-node_centers(channel);
            left_exponential = zeros(size(nu));
            right_exponential = zeros(size(nu));
            left_exponential(mask) = exp( ...
                branch_sign*tau*residual(mask));
            right_exponential(mask) = exp( ...
                -branch_sign*tau*residual(mask));
            for power_index=1:numel(powers)
                multiplier = power_multipliers{power_index};
                left{channel,power_index} = ifft2( ...
                    multiplier.*left_exponential.*u_spectrum);
                right{channel,power_index} = ifft2( ...
                    multiplier.*right_exponential.*u_spectrum);
                left_nu{channel,power_index} = ifft2( ...
                    nu.*multiplier.*left_exponential.*u_spectrum);
                right_nu{channel,power_index} = ifft2( ...
                    nu.*multiplier.*right_exponential.*u_spectrum);
            end
        end
        for first = 1:channel_count
            for second = 1:channel_count
                pair_support = node_pair_supports{first,second};
                pair_field = complex(zeros(size(q)));
                pair_sigma_field = complex(zeros(size(q)));
                for term_index=1:size(terms,1)
                    first_power = find(powers==terms(term_index,1),1);
                    second_power = find(powers==terms(term_index,2),1);
                    coefficient = terms(term_index,3);
                    pair_field = pair_field+coefficient ...
                        *left{first,first_power} ...
                        .*conj(right{second,second_power});
                    pair_sigma_field = pair_sigma_field+coefficient*( ...
                        left_nu{first,first_power} ...
                        .*conj(right{second,second_power}) ...
                        -left{first,first_power} ...
                        .*conj(right_nu{second,second_power}));
                end
                pair = fft2(pair_field);
                pair_sigma = fft2(pair_sigma_field);
                pair(~pair_support) = 0;
                pair_sigma(~pair_support) = 0;
                shift_difference = node_centers(first)-node_centers(second);
                factor = balanced_output_factor( ...
                    nu,projection,pair_support,balance_velocity, ...
                    shift_difference,tau,branch_sign);
                correction = correction ...
                    -(weights(node_index)/scale) ...
                    *factor.*(nu.*pair-branch_sign*pair_sigma);
            end
        end
    end
end
end

function correction = endpoint_correction_single( ...
        nu,q,projection,u_spectrum,masks,centers,pair_supports, ...
        terminal_masks,terminal_centers,terminal_pair_supports, ...
        terminal_radial4_node_count,balance_velocity,scale,nodes,weights, ...
        coefficient,endpoint_copies,branch_signs,include_explicit)
q_inverse = zeros(size(q));
q_inverse(q > 0) = 1./q(q > 0);
if include_explicit
    endpoint_bank = ifft2(q_inverse.*u_spectrum);
    correction = endpoint_copies*coefficient ...
        *fft2(endpoint_bank.*conj(endpoint_bank));
else
    correction = complex(zeros(size(q)));
end
for node_index = 1:numel(nodes)
    tau = nodes(node_index)/scale;
    node_masks = masks;
    node_centers = centers;
    node_pair_supports = pair_supports;
    if terminal_radial4_node_count > 0 ...
            && node_index > numel(nodes)-terminal_radial4_node_count
        node_masks = terminal_masks;
        node_centers = terminal_centers;
        node_pair_supports = terminal_pair_supports;
    end
    channel_count = numel(node_masks);
    for branch_sign = branch_signs
        left_zero = cell(1,channel_count);
        right_zero = cell(1,channel_count);
        left_nu = cell(1,channel_count);
        right_nu = cell(1,channel_count);
        for channel = 1:channel_count
            mask = node_masks{channel};
            residual = nu-balance_velocity*projection-node_centers(channel);
            left_exponential = zeros(size(nu));
            right_exponential = zeros(size(nu));
            left_exponential(mask) = exp( ...
                branch_sign*tau*residual(mask));
            right_exponential(mask) = exp( ...
                -branch_sign*tau*residual(mask));
            left_zero{channel} = ifft2( ...
                q_inverse.*left_exponential.*u_spectrum);
            right_zero{channel} = ifft2( ...
                q_inverse.*right_exponential.*u_spectrum);
            left_nu{channel} = ifft2( ...
                nu.*q_inverse.*left_exponential.*u_spectrum);
            right_nu{channel} = ifft2( ...
                nu.*q_inverse.*right_exponential.*u_spectrum);
        end
        for first = 1:channel_count
            for second = 1:channel_count
                pair_support = node_pair_supports{first,second};
                pair = fft2(left_zero{first}.*conj(right_zero{second}));
                pair_sigma = fft2( ...
                    left_nu{first}.*conj(right_zero{second}) ...
                    -left_zero{first}.*conj(right_nu{second}));
                pair(~pair_support) = 0;
                pair_sigma(~pair_support) = 0;
                shift_difference = node_centers(first)-node_centers(second);
                factor = balanced_output_factor( ...
                    nu,projection,pair_support,balance_velocity, ...
                    shift_difference,tau,branch_sign);
                correction = correction ...
                    -coefficient*(weights(node_index)/scale) ...
                    *factor.*(nu.*pair-branch_sign*pair_sigma);
            end
        end
    end
end
end

function multiplier = radial_power(q,power)
multiplier = zeros(size(q));
nonzero = q > 0;
multiplier(nonzero) = q(nonzero).^power;
end

function [spectrum,state] = shell_accumulation( ...
        base,nu,q,projection,u_spectrum,active,balance_velocity, ...
        delta_q,group_velocity,terms,endpoint_copies)
configuration = finite_depth_eta20_multiterm_shell_configuration();
levels = finite_depth_eta20_endpoint_subtracted_gl_configuration();
[masks1,centers1,supports1] = difference_balance_channels( ...
    active,nu-balance_velocity*projection,1);
[masks8,centers8,supports8] = difference_balance_channels( ...
    active,nu-balance_velocity*projection,8);
[masks12,centers12,supports12] = difference_balance_channels( ...
    active,nu-balance_velocity*projection,12);
spectrum = endpoint_correction(nu,q,projection,u_spectrum, ...
    masks1,centers1,supports1,masks1,centers1,supports1,0, ...
    balance_velocity,1,zeros(1,0),zeros(1,0),terms,endpoint_copies, ...
    zeros(1,0),true);
for shell_index=1:3
    order = configuration.shell_orders(shell_index);
    nodes = levels.(sprintf('gl%d_nodes',order));
    weights = levels.(sprintf('gl%d_weights',order)).*exp(nodes);
    if shell_index < 3
        anchor = configuration.shell_anchor_q(shell_index);
        slow_scale = anchor*(1-group_velocity);
        fast_scale = anchor*(1+group_velocity);
    else
        slow_scale = delta_q;
        fast_scale = delta_q;
    end
    if slow_scale <= 0 || fast_scale <= 0
        error('Frozen output-shell branch scale is not positive.');
    end
    shell_spectrum = complex(zeros(size(q)));
    for branch_sign=[1,-1]
        if branch_sign==1
            channel_count = configuration.slow_balance_channels(shell_index);
            scale = slow_scale;
        else
            channel_count = configuration.fast_balance_channels(shell_index);
            scale = fast_scale;
        end
        [masks,centers,pair_supports] = select_balance_partition( ...
            channel_count,masks1,centers1,supports1, ...
            masks8,centers8,supports8,masks12,centers12,supports12);
        shell_spectrum = shell_spectrum+branch_forcing_integral( ...
            base,nu,projection,u_spectrum,masks,centers,pair_supports, ...
            balance_velocity,scale,nodes,weights,branch_sign);
        shell_spectrum = shell_spectrum+endpoint_correction( ...
            nu,q,projection,u_spectrum,masks,centers,pair_supports, ...
            masks,centers,pair_supports,0,balance_velocity,scale, ...
            nodes,weights,terms,0,branch_sign,false);
    end
    if shell_index==1
        output_mask = q>0 & q<=configuration.shell_edges(2);
    elseif shell_index==2
        output_mask = q>configuration.shell_edges(2) ...
            & q<=configuration.shell_edges(3);
    else
        output_mask = q>configuration.shell_edges(3);
    end
    shell_spectrum(~output_mask) = 0;
    spectrum = spectrum+shell_spectrum;
end
spectrum = hermitian_project(spectrum);
state = struct('support_mask_transforms',230);
end

function result = branch_forcing_integral( ...
        base,nu,projection,u_spectrum,masks,centers,pair_supports, ...
        balance_velocity,scale,nodes,weights,branch_sign)
result = complex(zeros(size(nu)));
for node_index=1:numel(nodes)
    tau = nodes(node_index)/scale;
    [fd,fk] = forcing_at_time(base,nu,projection,u_spectrum,masks, ...
        centers,pair_supports,balance_velocity,tau,branch_sign);
    result = result+(weights(node_index)/(4*scale)) ...
        *(branch_sign*1i*fk-nu.*fd);
end
end

function [masks,centers,supports] = select_balance_partition( ...
        count,masks1,centers1,supports1,masks8,centers8,supports8, ...
        masks12,centers12,supports12)
if count==1
    masks=masks1; centers=centers1; supports=supports1;
elseif count==8
    masks=masks8; centers=centers8; supports=supports8;
elseif count==12
    masks=masks12; centers=centers12; supports=supports12;
else
    error('Unsupported frozen shell balance-channel count.');
end
end

function [fd_spectrum,fk_spectrum] = forcing_at_time( ...
    base,nu,projection,u_spectrum,masks,centers,pair_supports, ...
    balance_velocity,tau,branch_sign,suppress_roundoff)
if nargin<11,suppress_roundoff=false;end
channel_count = numel(masks);
left = cell(1,channel_count);
right = cell(1,channel_count);
for channel = 1:channel_count
    mask = masks{channel};
    residual = nu-balance_velocity*projection-centers(channel);
    left_exponential = zeros(size(nu));
    right_exponential = zeros(size(nu));
    left_exponential(mask) = exp(branch_sign*tau*residual(mask));
    right_exponential(mask) = exp(-branch_sign*tau*residual(mask));
    left{channel} = filtered_bank(base,u_spectrum,left_exponential);
    right{channel} = filtered_bank(base,u_spectrum,right_exponential);
end
fd_spectrum = complex(zeros(size(nu)));
fk_spectrum = complex(zeros(size(nu)));
for first = 1:channel_count
    for second = 1:channel_count
        forcing_d = ...
            -left{first}.dno.*conj(right{second}.one) ...
            -left{first}.one.*conj(right{second}.dno) ...
            +left{first}.nu.*conj(right{second}.nu) ...
            +left{first}.bx.*conj(right{second}.bx) ...
            +left{first}.by.*conj(right{second}.by);
        forcing_k = ...
             1i*left{first}.a.*conj(right{second}.one) ...
            -1i*left{first}.bx.*conj(right{second}.qx) ...
            -1i*left{first}.by.*conj(right{second}.qy) ...
            +1i*left{first}.qx.*conj(right{second}.bx) ...
            +1i*left{first}.qy.*conj(right{second}.by) ...
            -1i*left{first}.one.*conj(right{second}.a);
        pair_d = fft2(forcing_d);
        pair_k = fft2(forcing_k);
        pair_support = pair_supports{first,second};
        pair_d(~pair_support) = 0;
        pair_k(~pair_support) = 0;
        if suppress_roundoff
            pair_d = suppress_fft_roundoff(pair_d,pair_support);
            pair_k = suppress_fft_roundoff(pair_k,pair_support);
        end
        shift_difference = centers(first)-centers(second);
        factor = balanced_output_factor( ...
            nu,projection,pair_support,balance_velocity, ...
            shift_difference,tau,branch_sign);
        fd_spectrum = fd_spectrum+factor.*pair_d;
        fk_spectrum = fk_spectrum+factor.*pair_k;
    end
end
end

function spectrum = suppress_fft_roundoff(spectrum,support)
% Remove coefficients below a conservative double-precision FFT backward
% error before the separated output exponential can amplify them.
active_values=abs(spectrum(support));
if isempty(active_values),return,end
scale=max(active_values);
if scale==0,return,end
relative_floor=16*eps*max(1,log2(numel(spectrum)));
spectrum(support & abs(spectrum)<relative_floor*scale)=0;
end

function factor = balanced_output_factor( ...
        nu,projection,output_support,balance_velocity,shift_difference, ...
        tau,branch_sign)
factor = zeros(size(nu));
factor(output_support) = exp(-tau*(nu(output_support) ...
    -branch_sign*balance_velocity*projection(output_support) ...
    -branch_sign*shift_difference));
end

function [masks,centers,pair_supports] = ...
        difference_balance_channels(active,residual,count)
values = residual(active);
edges = linspace(min(values),max(values),count+1);
masks = cell(1,count);
centers = zeros(1,count);
use_equal_count_partition = false;
for channel = 1:count
    if channel < count
        masks{channel} = active & residual >= edges(channel) ...
            & residual < edges(channel+1);
    else
        masks{channel} = active & residual >= edges(channel) ...
            & residual <= edges(channel+1);
    end
    if ~any(masks{channel},'all')
        use_equal_count_partition = true;
        break
    end
    channel_values = residual(masks{channel});
    centers(channel) = 0.5*(min(channel_values)+max(channel_values));
end
if use_equal_count_partition
    active_indices = find(active);
    if numel(active_indices) < count
        error('The active support has fewer parents than balance channels.');
    end
    [~,order] = sort(residual(active_indices),'ascend');
    ordered_indices = active_indices(order);
    cuts = floor((0:count)*numel(ordered_indices)/count);
    for channel = 1:count
        masks{channel} = false(size(active));
        members = ordered_indices(cuts(channel)+1:cuts(channel+1));
        masks{channel}(members) = true;
        channel_values = residual(members);
        centers(channel) = 0.5*(min(channel_values)+max(channel_values));
    end
end
mask_transforms = cell(1,count);
for channel = 1:count
    mask_transforms{channel} = fft2(double(masks{channel}));
end
pair_supports = cell(count,count);
for first = 1:count
    for second = 1:count
        pair_supports{first,second} = real(ifft2( ...
            mask_transforms{first}.*conj(mask_transforms{second}))) > 0.5;
    end
end
end

function bank = filtered_bank(base,u_spectrum,exponential)
families = {'one','dno','nu','bx','by','qx','qy','a'};
bank = struct();
for family_index = 1:numel(families)
    family = families{family_index};
    bank.(family) = ifft2(base.(family).*exponential.*u_spectrum);
end
end

function projected = hermitian_project(spectrum)
[ny,nx] = size(spectrum);
partner_x = mod(-(0:nx-1),nx)+1;
partner_y = mod(-(0:ny-1),ny)+1;
projected = 0.5*(spectrum+conj(spectrum(partner_y,partner_x)));
end
