function [eta33_plus,audit,components,psi33_plus] = ...
    gl_no_stokes_eta33( ...
    eta11_spectrum,qx,qy,peak_depth_qp,quadrature_rank,project_root, ...
    outer_quadrature_rank)
%GL_NO_STOKES_ETA33 Pure Green--Laplace order-three elevation and surface potential.
% Port of the frozen graph with all diagonal-repair branches removed.

arguments
    eta11_spectrum (:,:) {mustBeNumeric}
    qx (:,:) double
    qy (:,:) double
    peak_depth_qp (1,1) double {mustBePositive}
    quadrature_rank (1,1) double {mustBeInteger,mustBePositive} = 6
    project_root (1,1) string = string( ...
        fileparts(fileparts(fileparts(mfilename('fullpath')))))
    outer_quadrature_rank (1,1) double = NaN
end
if isnan(outer_quadrature_rank),outer_quadrature_rank=quadrature_rank;end
if ~ismember(quadrature_rank,[4,6,8,10,12]) ...
        || ~ismember(outer_quadrature_rank,[3,4,6,8,10])
    error('Unsupported bounded inner/outer Green--Laplace rank pair.');
end
if ~isequal(size(eta11_spectrum),size(qx),size(qy))
    error('Spectrum, qx, and qy must have identical sizes.');
end
[plan,plan_cache_hit] = cached_execution_plan( ...
    qx,qy,peak_depth_qp,quadrature_rank,outer_quadrature_rank,project_root);
freeze = plan.freeze;
q = plan.q;
nu = plan.nu;
A = plan.A;
sqrt_A = plan.sqrt_A;
nonzero = plan.nonzero;
multipliers = plan.multipliers;
compute_psi = nargout>=4;
surface_only = false;

threshold = 1e-13*max(1,max(abs(eta11_spectrum(:))));
support = abs(eta11_spectrum)>threshold;
if ~any(support,'all') || any(qx(support)<=0) || ...
        any(q(support)<=0)
    error('Strict-forward eta11 support with every parent q>0 is required.');
end
assert_cubic_support_is_alias_safe(support);
pair_support = real(ifft2(fft2(double(support)).^2))>0.5;
triple_support = real(ifft2(fft2(double(support)).^3))>0.5;
support_fft_count = 4;

% Balance the separated Laplace factors without changing their algebra.
% The unbalanced graph forms an exponentially small convolution and then
% multiplies it by cosh(sqrt(A)*t).  Far from the peak-centred scale that
% exposes FFT roundoff to enormous output multipliers.  A common shift c
% on every parent leaf is exactly cancelled by exp(-n*c*t) after the
% order-n product.  The midpoint below balances the worst low-parent and
% high-output exponents on the declared support.
minimum_parent_nu = min(nu(support));
maximum_pair_a = max(sqrt_A(pair_support));
maximum_triple_a = max(sqrt_A(triple_support));
inner_balance_shift = 0.5*(minimum_parent_nu+maximum_pair_a/2);
outer_balance_shift = 0.5*(minimum_parent_nu+maximum_triple_a/3);

inner_nodes = plan.inner_nodes;
inner_weights = plan.inner_weights;
outer_nodes = plan.outer_nodes;
outer_weights = plan.outer_weights;
lambda2 = plan.lambda2;
lambda3 = plan.lambda3;
if ~isfinite(lambda2) || ~isfinite(lambda3) ...
        || lambda2<=0 || lambda3<=0
    error('Frozen Green--Laplace tail scales are invalid.');
end

candidate_spectrum = complex(zeros(size(q)));
if compute_psi
    psi_spectrum = complex(zeros(size(q)));
else
    psi_spectrum = [];
end
if use_cached_node_factors()
    inner_node_factors=node_factor_bank( ...
        nu,sqrt_A,inner_nodes,inner_weights,lambda2,inner_balance_shift,2);
    outer_node_factors=node_factor_bank( ...
        nu,sqrt_A,outer_nodes,outer_weights,lambda3,outer_balance_shift,3);
else
    inner_node_factors=[];
    outer_node_factors=[];
end
fft_count = support_fft_count;
product_count = 0;
outer_node_diagnostics = repmat(struct( ...
    'node',0,'time',0,'pair_eta_max',0,'pair_phi_max',0, ...
    'pair_eta_core_spectrum_max',0,'pair_eta_correction_spectrum_max',0, ...
    'forcing_k_spectrum_max',0,'forcing_d_spectrum_max',0, ...
    'increment_spectrum_max',0),1,outer_quadrature_rank);
for outer_index = 1:outer_quadrature_rank
    node = outer_nodes(outer_index);
    time = node/lambda3;
    if isempty(outer_node_factors)
        damped_spectrum = eta11_spectrum.*exp( ...
            -time*(nu-outer_balance_shift));
    else
        damped_spectrum=eta11_spectrum.* ...
            outer_node_factors.damping(:,:,outer_index);
    end
    [first,first_fft] = first_order_state( ...
        damped_spectrum,multipliers);
    [pair,pair_fft,pair_products,pair_diagnostics] = inner_pair_fields( ...
        damped_spectrum,first,q,A,sqrt_A,multipliers, ...
        pair_support,inner_nodes,inner_weights,lambda2, ...
        inner_balance_shift,compute_psi, ...
        inner_node_factors);
    fft_count = fft_count+first_fft+pair_fft;
    product_count = product_count+pair_products;

    forcing_k_pair = first.phix.*pair.etax+first.phiy.*pair.etay ...
        +pair.phix.*first.etax+pair.phiy.*first.etay ...
        -first.eta.*pair.phizz-pair.eta.*first.phizz;
    forcing_k_direct = first.eta.*(first.phixz.*first.etax ...
        +first.phiyz.*first.etay) ...
        -0.5*first.eta.^2.*first.phizzz;
    forcing_d_pair = first.eta.*pair.phitz ...
        +pair.eta.*first.phitz ...
        +first.phix.*pair.phix+first.phiy.*pair.phiy ...
        +first.phiz.*pair.phiz;
    forcing_d_direct = 0.5*first.eta.^2.*first.phitzz ...
        +first.eta.*(first.phix.*first.phixz ...
        +first.phiy.*first.phiyz+first.phiz.*first.phizz);
    forcing_k = forcing_k_direct/4+forcing_k_pair/2;
    forcing_d = forcing_d_direct/4+forcing_d_pair/2;
    product_count = product_count+19;
    forcing_k_spectrum = fft2(forcing_k);
    forcing_d_spectrum = fft2(forcing_d);
    fft_count = fft_count+2;
    forcing_k_spectrum(~triple_support) = 0;
    forcing_d_spectrum(~triple_support) = 0;
    forcing_k_spectrum = suppress_fft_roundoff( ...
        forcing_k_spectrum,triple_support);
    forcing_d_spectrum = suppress_fft_roundoff( ...
        forcing_d_spectrum,triple_support);
    if compute_psi
        [contact,contact_t,contact_products] = ...
            surface_taylor_contact(first,pair);
        contact_spectrum = fft2(contact);
        contact_t_spectrum = fft2(contact_t);
        contact_spectrum(~triple_support) = 0;
        contact_t_spectrum(~triple_support) = 0;
        contact_spectrum = suppress_fft_roundoff( ...
            contact_spectrum,triple_support);
        contact_t_spectrum = suppress_fft_roundoff( ...
            contact_t_spectrum,triple_support);
        forcing_k_psi_spectrum = forcing_k_spectrum+A.*contact_spectrum;
        forcing_d_psi_spectrum = forcing_d_spectrum-contact_t_spectrum;
        forcing_k_psi_spectrum = suppress_fft_roundoff( ...
            forcing_k_psi_spectrum,triple_support);
        forcing_d_psi_spectrum = suppress_fft_roundoff( ...
            forcing_d_psi_spectrum,triple_support);
        fft_count = fft_count+2;
        product_count = product_count+contact_products;
    end
    if isempty(outer_node_factors)
        sinh_over_a = time*ones(size(q));
        sinh_over_a(nonzero) = ...
            sinh(sqrt_A(nonzero)*time)./sqrt_A(nonzero);
        cosh_factor=cosh(sqrt_A*time);
        coefficient = outer_weights(outer_index)*exp(node)/lambda3 ...
            *exp(-3*outer_balance_shift*time);
    else
        sinh_over_a=outer_node_factors.sinh_over_a(:,:,outer_index);
        cosh_factor=outer_node_factors.cosh_factor(:,:,outer_index);
        coefficient=outer_node_factors.coefficient(outer_index);
    end
    if surface_only
        increment_spectrum_max=0;
    else
        increment_spectrum = coefficient*( ...
            A.*sinh_over_a.*forcing_d_spectrum ...
            -1i*cosh_factor.*forcing_k_spectrum);
        candidate_spectrum = candidate_spectrum+increment_spectrum;
        increment_spectrum_max=max(abs(increment_spectrum),[],'all');
    end
    if compute_psi
        psi_spectrum = psi_spectrum+coefficient*( ...
            -1i*cosh_factor.*forcing_d_psi_spectrum ...
            -sinh_over_a.*forcing_k_psi_spectrum);
    end
    outer_node_diagnostics(outer_index) = struct( ...
        'node',node,'time',time, ...
        'pair_eta_max',max(abs(pair.eta),[],'all'), ...
        'pair_phi_max',max(abs(pair.phix),[],'all'), ...
        'pair_eta_core_spectrum_max', ...
            pair_diagnostics.eta_core_spectrum_max, ...
        'pair_eta_correction_spectrum_max', ...
            pair_diagnostics.eta_correction_spectrum_max, ...
        'forcing_k_spectrum_max',max(abs(forcing_k_spectrum),[],'all'), ...
        'forcing_d_spectrum_max',max(abs(forcing_d_spectrum),[],'all'), ...
        'increment_spectrum_max',increment_spectrum_max);
end

if surface_only
    raw_resolvent_spectrum=[];
else
    raw_resolvent_spectrum = candidate_spectrum;
end
crossing_correction_spectrum = complex(zeros(size(q)));
if compute_psi,psi_raw_resolvent_spectrum = psi_spectrum;end
if compute_psi
    psi_crossing_correction_spectrum = complex(zeros(size(q)));
end
if ~surface_only,candidate_spectrum(~triple_support) = 0;end
psi_spectrum(~triple_support) = 0;
if surface_only
    eta33_plus=[];
else
    eta33_plus = ifft2(candidate_spectrum);
    fft_count = fft_count+1;
end
if compute_psi
    psi33_plus = ifft2(psi_spectrum);
    fft_count = fft_count+1;
else
    psi33_plus = [];
end
if ~surface_only && any(~isfinite(eta33_plus),'all')
    error('Nested Green--Laplace eta33 output is not finite.');
end
if compute_psi && any(~isfinite(psi33_plus),'all')
    error('Direct Green--Laplace Psi33 output is not finite.');
end

audit = struct( ...
    'candidate_id',string(freeze.candidate_id), ...
    'executor_revision','stable-fft-v2', ...
    'numerical_stability_policy',[ ...
        'balanced separated exponentials plus source-local ' ...
        'double-precision FFT backward-error suppression'], ...
    'status','frozen-symbolic-pre-mf12', ...
    'information_boundary','eta11-only', ...
    'internal_generated_fields',{{'eta22_plus','Phi22_plus','Phi22_t_plus'}}, ...
    'input_chi_route_used',false, ...
    'stokes_correction_used',false, ...
    'oracle_or_mf12_used',false, ...
    'inner_quadrature_rank',quadrature_rank, ...
    'outer_quadrature_rank',outer_quadrature_rank, ...
    'inner_scale',lambda2, ...
    'outer_scale',lambda3, ...
    'inner_balance_shift',inner_balance_shift, ...
    'outer_balance_shift',outer_balance_shift, ...
    'outer_node_diagnostics',outer_node_diagnostics, ...
    'static_plan_cache_hit',plan_cache_hit, ...
    'crossing_gate','none', ...
    'core_fft_ifft_count',fft_count-support_fft_count, ...
    'support_mask_fft_ifft_count',support_fft_count, ...
    'fft_ifft_count',fft_count, ...
    'pointwise_product_count',product_count, ...
    'pair_loops',0, ...
    'triple_loops',0, ...
    'cost_class','fixed-constant N log N');
if compute_psi
    audit.output_variable = 'surface Psi33_plus';
    audit.direct_surface_source_transform = true;
end
components = struct( ...
    'raw_resolvent_spectrum',raw_resolvent_spectrum, ...
    'crossing_correction_spectrum',crossing_correction_spectrum, ...
    'total_spectrum',candidate_spectrum);
if compute_psi
    psi_raw_resolvent_spectrum(~triple_support) = 0;
    psi_crossing_correction_spectrum(~triple_support) = 0;
    components.Psi33_direct_spectrum = psi_spectrum;
    components.Psi33_raw_resolvent_spectrum = psi_raw_resolvent_spectrum;
    components.Psi33_crossing_correction_spectrum = ...
        psi_crossing_correction_spectrum;
end
end

function [state,count] = first_order_state(spectrum,multipliers)
state = struct();
if use_batched_fft()
    fields=ifft2(cat(3,spectrum, ...
        multipliers.dx.*spectrum,multipliers.dy.*spectrum, ...
        multipliers.qx_over_nu.*spectrum, ...
        multipliers.qy_over_nu.*spectrum, ...
        multipliers.minus_i_nu.*spectrum, ...
        multipliers.nu_qx.*spectrum,multipliers.nu_qy.*spectrum, ...
        multipliers.minus_i_q2_over_nu.*spectrum, ...
        multipliers.minus_i_q2_nu.*spectrum, ...
        multipliers.minus_A.*spectrum,multipliers.minus_q2.*spectrum));
    names={'eta','etax','etay','phix','phiy','phiz','phixz','phiyz', ...
        'phizz','phizzz','phitz','phitzz'};
    for field_index=1:numel(names)
        state.(names{field_index})=fields(:,:,field_index);
    end
else
    state.eta = ifft2(spectrum);
    state.etax = ifft2(multipliers.dx.*spectrum);
    state.etay = ifft2(multipliers.dy.*spectrum);
    state.phix = ifft2(multipliers.qx_over_nu.*spectrum);
    state.phiy = ifft2(multipliers.qy_over_nu.*spectrum);
    state.phiz = ifft2(multipliers.minus_i_nu.*spectrum);
    state.phixz = ifft2(multipliers.nu_qx.*spectrum);
    state.phiyz = ifft2(multipliers.nu_qy.*spectrum);
    state.phizz = ifft2(multipliers.minus_i_q2_over_nu.*spectrum);
    state.phizzz = ifft2(multipliers.minus_i_q2_nu.*spectrum);
    state.phitz = ifft2(multipliers.minus_A.*spectrum);
    state.phitzz = ifft2(multipliers.minus_q2.*spectrum);
end
state.eta_t = state.phiz;
state.phiz_t = state.phitz;
state.phizz_t = state.phitzz;
count = 12;
end

function [pair,count,products,diagnostics] = inner_pair_fields( ...
    spectrum,first,q,A,sqrt_A,multipliers,output_support, ...
    nodes,weights,lambda,balance_shift, ...
    compute_psi,node_factors)
eta_spectrum = complex(zeros(size(q)));
etat_spectrum = complex(zeros(size(q)));
phi_spectrum = complex(zeros(size(q)));
phit_spectrum = complex(zeros(size(q)));
count = 0;
products = 0;
nonzero = q>0;
for index = 1:numel(nodes)
    node = nodes(index);
    time = node/lambda;
    if isempty(node_factors)
        damped = spectrum.*exp(-time*(multipliers.nu-balance_shift));
    else
        damped=spectrum.*node_factors.damping(:,:,index);
    end
    if use_batched_fft()
        filtered=ifft2(cat(3,damped, ...
            damped.*multipliers.nu,damped.*multipliers.nu2, ...
            damped.*multipliers.nu3, ...
            damped.*multipliers.qx_over_nu, ...
            damped.*multipliers.qy_over_nu, ...
            damped.*multipliers.qx,damped.*multipliers.qy, ...
            damped.*multipliers.q2_over_nu,damped.*multipliers.q2, ...
            damped.*multipliers.nu_qx,damped.*multipliers.nu_qy));
        v=filtered(:,:,1);v_nu=filtered(:,:,2);v_nu2=filtered(:,:,3);
        v_nu3=filtered(:,:,4);hx=filtered(:,:,5);hy=filtered(:,:,6);
        jx=filtered(:,:,7);jy=filtered(:,:,8);
        radial_over_nu=filtered(:,:,9);q2_field=filtered(:,:,10);
        nu_jx=filtered(:,:,11);nu_jy=filtered(:,:,12);
    else
        v = ifft2(damped);
        v_nu = ifft2(damped.*multipliers.nu);
        v_nu2 = ifft2(damped.*multipliers.nu2);
        v_nu3 = ifft2(damped.*multipliers.nu3);
        hx = ifft2(damped.*multipliers.qx_over_nu);
        hy = ifft2(damped.*multipliers.qy_over_nu);
        jx = ifft2(damped.*multipliers.qx);
        jy = ifft2(damped.*multipliers.qy);
        radial_over_nu = ifft2(damped.*multipliers.q2_over_nu);
        q2_field = ifft2(damped.*multipliers.q2);
        nu_jx = ifft2(damped.*multipliers.nu_qx);
        nu_jy = ifft2(damped.*multipliers.nu_qy);
    end
    count = count+12;

    source_d = 2*v.*v_nu2+v_nu.^2-hx.^2-hy.^2;
    source_k = 2*v.*radial_over_nu+2*(hx.*jx+hy.*jy);
    source_d_s = 2*v.*v_nu3+4*v_nu.*v_nu2 ...
        -2*hx.*jx-2*hy.*jy;
    source_k_s = 2*(v_nu.*radial_over_nu+v.*q2_field) ...
        +2*(jx.^2+hx.*nu_jx+jy.^2+hy.*nu_jy);
    products = products+18;
    if use_batched_fft()
        source_spectra=fft2(cat(3,source_d,source_k,source_d_s,source_k_s));
        sd=source_spectra(:,:,1);sk=source_spectra(:,:,2);
        sd_s=source_spectra(:,:,3);sk_s=source_spectra(:,:,4);
    else
        sd = fft2(source_d);
        sk = fft2(source_k);
        sd_s = fft2(source_d_s);
        sk_s = fft2(source_k_s);
    end
    count = count+4;
    sd(~output_support) = 0;
    sk(~output_support) = 0;
    sd_s(~output_support) = 0;
    sk_s(~output_support) = 0;
    sd = suppress_fft_roundoff(sd,output_support);
    sk = suppress_fft_roundoff(sk,output_support);
    sd_s = suppress_fft_roundoff(sd_s,output_support);
    sk_s = suppress_fft_roundoff(sk_s,output_support);
    if isempty(node_factors)
        sinh_over_a = time*ones(size(q));
        sinh_over_a(nonzero) = ...
            sinh(sqrt_A(nonzero)*time)./sqrt_A(nonzero);
        cosh_factor = cosh(sqrt_A*time);
        coefficient = weights(index)*exp(node)/lambda ...
            *exp(-2*balance_shift*time)/4;
    else
        sinh_over_a=node_factors.sinh_over_a(:,:,index);
        cosh_factor=node_factors.cosh_factor(:,:,index);
        coefficient=node_factors.coefficient(index)/4;
    end
    eta_spectrum = eta_spectrum+coefficient*( ...
        -A.*sinh_over_a.*sd+cosh_factor.*sk);
    if compute_psi
        etat_spectrum = etat_spectrum-1i*coefficient*( ...
            -A.*sinh_over_a.*sd_s+cosh_factor.*sk_s);
    end
    phi_spectrum = phi_spectrum+1i*coefficient*( ...
        cosh_factor.*sd-sinh_over_a.*sk);
    phit_spectrum = phit_spectrum+coefficient*( ...
        cosh_factor.*sd_s-sinh_over_a.*sk_s);
end

eta_core_spectrum = eta_spectrum;

eta_correction_spectrum = complex(zeros(size(q)));
eta_spectrum(~output_support) = 0;
etat_spectrum(~output_support) = 0;
phi_spectrum(~output_support) = 0;
phit_spectrum(~output_support) = 0;

if compute_psi
    etat_spectrum = suppress_fft_roundoff(etat_spectrum,output_support);
end
if use_batched_fft()
    output_bank=cat(3,eta_spectrum, ...
        multipliers.dx.*eta_spectrum,multipliers.dy.*eta_spectrum, ...
        multipliers.dx.*phi_spectrum,multipliers.dy.*phi_spectrum, ...
        A.*phi_spectrum,multipliers.q2.*phi_spectrum,A.*phit_spectrum);
    if compute_psi,output_bank=cat(3,output_bank,etat_spectrum);end
    output_fields=ifft2(output_bank);
    pair.eta=output_fields(:,:,1);pair.etax=output_fields(:,:,2);
    pair.etay=output_fields(:,:,3);pair.phix=output_fields(:,:,4);
    pair.phiy=output_fields(:,:,5);pair.phiz=output_fields(:,:,6);
    pair.phizz=output_fields(:,:,7);pair.phitz=output_fields(:,:,8);
    count=count+8;
    if compute_psi
        pair.eta_t=output_fields(:,:,9);pair.phiz_t=pair.phitz;
        count=count+1;
    end
else
    pair.eta = ifft2(eta_spectrum);
    pair.etax = ifft2(multipliers.dx.*eta_spectrum);
    pair.etay = ifft2(multipliers.dy.*eta_spectrum);
    pair.phix = ifft2(multipliers.dx.*phi_spectrum);
    pair.phiy = ifft2(multipliers.dy.*phi_spectrum);
    pair.phiz = ifft2(A.*phi_spectrum);
    pair.phizz = ifft2(multipliers.q2.*phi_spectrum);
    pair.phitz = ifft2(A.*phit_spectrum);
    count = count+8;
    if compute_psi
        pair.eta_t = ifft2(etat_spectrum);
        pair.phiz_t = pair.phitz;
        count = count+1;
    end
end
diagnostics = struct( ...
    'eta_core_spectrum_max',max(abs(eta_core_spectrum),[],'all'), ...
    'eta_correction_spectrum_max', ...
        max(abs(eta_correction_spectrum),[],'all'));
end

function enabled=use_batched_fft()
persistent cached_raw cached_enabled
raw=string(getenv('GL_MATLAB_GL_BATCHED_FFT'));
if isempty(cached_raw) || raw~=cached_raw
    cached_raw=raw;
    cached_enabled=ismember(lower(strtrim(raw)),["1","true","yes","on"]);
end
enabled=cached_enabled;
end

function enabled=use_cached_node_factors()
persistent cached_raw cached_enabled
raw=string(getenv('GL_MATLAB_GL_CACHE_FACTORS'));
if isempty(cached_raw) || raw~=cached_raw
    cached_raw=raw;
    cached_enabled=ismember(lower(strtrim(raw)),["1","true","yes","on"]);
end
enabled=cached_enabled;
end

function enabled=use_surface_only_runtime()
persistent cached_raw cached_enabled
raw=string(getenv('GL_MATLAB_GL_SURFACE_ONLY'));
if isempty(cached_raw) || raw~=cached_raw
    cached_raw=raw;
    cached_enabled=ismember(lower(strtrim(raw)),["1","true","yes","on"]);
end
enabled=cached_enabled;
end

function factors=node_factor_bank(nu,sqrt_A,nodes,weights,lambda,shift,order)
page_count=numel(nodes);
shape=[size(nu),page_count];
factors=struct('damping',zeros(shape),'sinh_over_a',zeros(shape), ...
    'cosh_factor',zeros(shape),'coefficient',zeros(1,page_count));
nonzero=sqrt_A>0;
for index=1:page_count
    node=nodes(index);time=node/lambda;
    factors.damping(:,:,index)=exp(-time*(nu-shift));
    page=time*ones(size(nu));
    page(nonzero)=sinh(sqrt_A(nonzero)*time)./sqrt_A(nonzero);
    factors.sinh_over_a(:,:,index)=page;
    factors.cosh_factor(:,:,index)=cosh(sqrt_A*time);
    factors.coefficient(index)=weights(index)*exp(node)/lambda ...
        *exp(-order*shift*time);
end
end

function [contact,contact_t,products] = surface_taylor_contact(first,pair)
% Analytic-plus order-three surface Taylor contact frozen by Wolfram.
contact = 0.5*first.eta.*pair.phiz ...
    +0.5*pair.eta.*first.phiz ...
    +0.125*first.eta.^2.*first.phizz;
contact_t = 0.5*(first.eta_t.*pair.phiz ...
    +first.eta.*pair.phiz_t) ...
    +0.5*(pair.eta_t.*first.phiz ...
    +pair.eta.*first.phiz_t) ...
    +0.25*first.eta.*first.eta_t.*first.phizz ...
    +0.125*first.eta.^2.*first.phizz_t;
products = 12;
end

function spectrum = suppress_fft_roundoff(spectrum,support)
% Remove coefficients below a conservative double-precision FFT backward
% error.  Such bins cannot carry resolved convolution information, while a
% separated Green--Laplace output multiplier can amplify them enormously.
active_values = abs(spectrum(support));
if isempty(active_values), return, end
scale = max(active_values);
if scale==0, return, end
relative_floor = 16*eps*max(1,log2(numel(spectrum)));
spectrum(support & abs(spectrum)<relative_floor*scale) = 0;
end

function [plan,cache_hit] = cached_execution_plan( ...
    qx,qy,peak_depth_qp,inner_rank,outer_rank,project_root)
persistent cached_plan cached_qx cached_qy cached_peak cached_inner_rank ...
    cached_outer_rank cached_root
cache_hit = ~isempty(cached_plan) ...
    && isequal(cached_qx,qx) ...
    && isequal(cached_qy,qy) ...
    && isequal(cached_peak,peak_depth_qp) ...
    && isequal(cached_inner_rank,inner_rank) ...
    && isequal(cached_outer_rank,outer_rank) ...
    && isequal(cached_root,project_root);
if cache_hit
    plan = cached_plan;
    return
end
freeze = load_and_verify_interfaces(project_root);
q = hypot(qx,qy);
nu = sqrt(q.*tanh(q));
A = q.*tanh(q);
sqrt_A = sqrt(A);
safe_nu = max(nu,realmin);
nonzero = q>0;
unit_x = zeros(size(q));
unit_y = zeros(size(q));
unit_x(nonzero) = qx(nonzero)./q(nonzero);
unit_y(nonzero) = qy(nonzero)./q(nonzero);
[inner_nodes,inner_weights] = gauss_laguerre_rule(inner_rank);
[outer_nodes,outer_weights] = gauss_laguerre_rule(outer_rank);
lambda2 = 2*sqrt(peak_depth_qp*tanh(peak_depth_qp)) ...
    -sqrt(2*peak_depth_qp*tanh(2*peak_depth_qp));
lambda3 = 3*sqrt(peak_depth_qp*tanh(peak_depth_qp)) ...
    -sqrt(3*peak_depth_qp*tanh(3*peak_depth_qp));
q2 = q.^2;
unit_x2 = unit_x.^2;
unit_xy = unit_x.*unit_y;
multipliers = struct( ...
    'qx',qx, ...
    'qy',qy, ...
    'q2',q2, ...
    'nu',nu, ...
    'nu2',nu.^2, ...
    'nu3',nu.^3, ...
    'dx',1i*qx, ...
    'dy',1i*qy, ...
    'qx_over_nu',qx./safe_nu, ...
    'qy_over_nu',qy./safe_nu, ...
    'q2_over_nu',q2./safe_nu, ...
    'minus_i_nu',-1i*nu, ...
    'nu_qx',nu.*qx, ...
    'nu_qy',nu.*qy, ...
    'minus_i_q2_over_nu',-1i*q2./safe_nu, ...
    'minus_i_q2_nu',-1i*q2.*nu, ...
    'minus_A',-A, ...
    'minus_q2',-q2, ...
    'unit_x',unit_x, ...
    'unit_y',unit_y, ...
    'unit_x2',unit_x2, ...
    'unit_xy',unit_xy, ...
    'nu_unit_x',nu.*unit_x, ...
    'nu_unit_y',nu.*unit_y, ...
    'nu_unit_x2',nu.*unit_x2, ...
    'nu_unit_xy',nu.*unit_xy);
plan = struct( ...
    'freeze',freeze, ...
    'q',q, ...
    'nu',nu, ...
    'A',A, ...
    'sqrt_A',sqrt_A, ...
    'safe_nu',safe_nu, ...
    'nonzero',nonzero, ...
    'unit_x',unit_x, ...
    'unit_y',unit_y, ...
    'inner_nodes',inner_nodes, ...
    'inner_weights',inner_weights, ...
    'outer_nodes',outer_nodes, ...
    'outer_weights',outer_weights, ...
    'lambda2',lambda2, ...
    'lambda3',lambda3);
plan.multipliers = multipliers;
cached_plan = plan;
cached_qx = qx;
cached_qy = qy;
cached_peak = peak_depth_qp;
cached_inner_rank = inner_rank;
cached_outer_rank = outer_rank;
cached_root = project_root;
end

function freeze = load_and_verify_interfaces(project_root)
freeze_path = fullfile(project_root,'symbolic','generated', ...
    'finite_depth_directional_order3_nested_green_laplace.json');
if ~isfile(freeze_path)
    error('Frozen Wolfram nested eta33 interface is missing.');
end
freeze = jsondecode(fileread(freeze_path));
if ~freeze.overall_exact_gate_pass || freeze.oracle_or_mf12_used ...
        || freeze.input_chi_route_used ...
        || freeze.inner_quadrature_rank~=4 ...
        || freeze.outer_quadrature_rank~=4
    error('Nested eta33 freeze crossed its information boundary.');
end
end

function [nodes,weights] = gauss_laguerre_rule(rank)
% Golub--Welsch rule for the standard weight exp(-x) on [0,infinity).
indices = 1:rank;
jacobi = diag(2*indices-1)+diag(1:rank-1,1)+diag(1:rank-1,-1);
[vectors,values] = eig(jacobi,'vector');
[nodes,order] = sort(real(values));
nodes = nodes.';
weights = vectors(1,order).^2;
end

function assert_cubic_support_is_alias_safe(support)
[ny,nx] = size(support);
mode_x = [0:(ceil(nx/2)-1),-floor(nx/2):-1];
mode_y = [0:(ceil(ny/2)-1),-floor(ny/2):-1];
[column,row] = meshgrid(mode_x,mode_y);
if max(3*abs(column(support)))>=nx/2 ...
        || max(3*abs(row(support)))>=ny/2
    error('Cubic output support reaches a Nyquist boundary.');
end
end
