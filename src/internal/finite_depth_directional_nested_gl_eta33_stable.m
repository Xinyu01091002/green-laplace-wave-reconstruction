function [eta33_plus,audit,components,psi33_plus] = ...
    finite_depth_directional_nested_gl_eta33_stable( ...
    eta11_spectrum,qx,qy,peak_depth_qp,project_root)
%FINITE_DEPTH_DIRECTIONAL_NESTED_GREEN_LAPLACE_ETA33 Frozen eta11-only graph.
% The inner GL4 response generates eta22, flat Phi22, and d_t Phi22 from
% damped eta11 leaves. The original order-three kinematic/dynamic forcing is
% then passed through the outer GL4 response. A permutation-symmetric G3
% endpoint repair restores the exact finite-depth eta33 Stokes diagonal.
% Requesting the optional fourth output activates the separately frozen
% direct surface-Psi33 K1 graph on the same stable-fft-v2 backend.

arguments
    eta11_spectrum (:,:) {mustBeNumeric}
    qx (:,:) double
    qy (:,:) double
    peak_depth_qp (1,1) double {mustBePositive}
    project_root (1,1) string = string( ...
        fileparts(fileparts(fileparts(mfilename('fullpath')))))
end
if ~isequal(size(eta11_spectrum),size(qx),size(qy))
    error('Spectrum, qx, and qy must have identical sizes.');
end
[plan,plan_cache_hit] = cached_execution_plan( ...
    qx,qy,peak_depth_qp,project_root);
freeze = plan.freeze;
gate = plan.gate;
q = plan.q;
nu = plan.nu;
A = plan.A;
sqrt_A = plan.sqrt_A;
nonzero = plan.nonzero;
multipliers = plan.multipliers;
compute_psi = nargout>=4;

threshold = 1e-13*max(1,max(abs(eta11_spectrum(:))));
support = abs(eta11_spectrum)>threshold;
if ~any(support,'all') || any(qx(support)<=0) || ...
        any(q(support)<=0.5)
    error('Strict-forward eta11 support with every parent q>0.5 is required.');
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

nodes = plan.nodes;
weights = plan.weights;
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
fft_count = support_fft_count;
product_count = 0;
outer_node_diagnostics = repmat(struct( ...
    'node',0,'time',0,'pair_eta_max',0,'pair_phi_max',0, ...
    'pair_eta_core_spectrum_max',0,'pair_eta_correction_spectrum_max',0, ...
    'forcing_k_spectrum_max',0,'forcing_d_spectrum_max',0, ...
    'increment_spectrum_max',0),1,4);
for outer_index = 1:4
    node = nodes(outer_index);
    time = node/lambda3;
    damped_spectrum = eta11_spectrum.*exp( ...
        -time*(nu-outer_balance_shift));
    [first,first_fft] = first_order_state( ...
        damped_spectrum,multipliers);
    [pair,pair_fft,pair_products,pair_diagnostics] = inner_pair_fields( ...
        damped_spectrum,first,q,A,sqrt_A,multipliers, ...
        pair_support,nodes,weights,lambda2, ...
        plan.pair_eta_correction,plan.pair_phi_correction, ...
        inner_balance_shift,compute_psi);
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
    sinh_over_a = time*ones(size(q));
    sinh_over_a(nonzero) = ...
        sinh(sqrt_A(nonzero)*time)./sqrt_A(nonzero);
    coefficient = weights(outer_index)*exp(node)/lambda3 ...
        *exp(-3*outer_balance_shift*time);
    increment_spectrum = coefficient*( ...
        A.*sinh_over_a.*forcing_d_spectrum ...
        -1i*cosh(sqrt_A*time).*forcing_k_spectrum);
    candidate_spectrum = candidate_spectrum+increment_spectrum;
    if compute_psi
        psi_spectrum = psi_spectrum+coefficient*( ...
            -1i*cosh(sqrt_A*time).*forcing_d_psi_spectrum ...
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
        'increment_spectrum_max',max(abs(increment_spectrum),[],'all'));
end

[g3_source,g3_fft,g3_products] = triple_gate_source( ...
    eta11_spectrum,multipliers,gate);
correction_spectrum = fft2(g3_source);
fft_count = fft_count+g3_fft+1;
product_count = product_count+g3_products;
correction = plan.third_order_correction;
correction(q==0) = 0;
raw_resolvent_spectrum = candidate_spectrum;
crossing_correction_spectrum = correction_spectrum.*correction;
candidate_spectrum = candidate_spectrum ...
    +crossing_correction_spectrum;
if compute_psi
    psi_correction = cached_third_order_surface_stokes_correction( ...
        q/3,nodes,weights,lambda3);
    psi_correction(q==0) = 0;
    psi_spectrum = psi_spectrum+correction_spectrum.*psi_correction;
end
candidate_spectrum(~triple_support) = 0;
psi_spectrum(~triple_support) = 0;
eta33_plus = ifft2(candidate_spectrum);
fft_count = fft_count+1;
if compute_psi
    psi33_plus = ifft2(psi_spectrum);
    fft_count = fft_count+1;
else
    psi33_plus = [];
end
if any(~isfinite(eta33_plus),'all')
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
    'oracle_or_mf12_used',false, ...
    'inner_quadrature_rank',4, ...
    'outer_quadrature_rank',4, ...
    'inner_scale',lambda2, ...
    'outer_scale',lambda3, ...
    'inner_balance_shift',inner_balance_shift, ...
    'outer_balance_shift',outer_balance_shift, ...
    'outer_node_diagnostics',outer_node_diagnostics, ...
    'static_plan_cache_hit',plan_cache_hit, ...
    'crossing_gate','G3=product_pairs((1+cosine_ij)/2)', ...
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
    components.Psi33_direct_spectrum = psi_spectrum;
end
end

function [state,count] = first_order_state(spectrum,multipliers)
state = struct();
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
state.eta_t = state.phiz;
state.phiz_t = state.phitz;
state.phizz_t = state.phitzz;
count = 12;
end

function [pair,count,products,diagnostics] = inner_pair_fields( ...
    spectrum,first,q,A,sqrt_A,multipliers,output_support, ...
    nodes,weights,lambda,eta_correction,phi_correction,balance_shift, ...
    compute_psi)
eta_spectrum = complex(zeros(size(q)));
etat_spectrum = complex(zeros(size(q)));
phi_spectrum = complex(zeros(size(q)));
phit_spectrum = complex(zeros(size(q)));
count = 0;
products = 0;
nonzero = q>0;
for index = 1:4
    node = nodes(index);
    time = node/lambda;
    damped = spectrum.*exp(-time*(multipliers.nu-balance_shift));
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
    count = count+12;

    source_d = 2*v.*v_nu2+v_nu.^2-hx.^2-hy.^2;
    source_k = 2*v.*radial_over_nu+2*(hx.*jx+hy.*jy);
    source_d_s = 2*v.*v_nu3+4*v_nu.*v_nu2 ...
        -2*hx.*jx-2*hy.*jy;
    source_k_s = 2*(v_nu.*radial_over_nu+v.*q2_field) ...
        +2*(jx.^2+hx.*nu_jx+jy.^2+hy.*nu_jy);
    products = products+18;
    sd = fft2(source_d);
    sk = fft2(source_k);
    sd_s = fft2(source_d_s);
    sk_s = fft2(source_k_s);
    count = count+4;
    sd(~output_support) = 0;
    sk(~output_support) = 0;
    sd_s(~output_support) = 0;
    sk_s(~output_support) = 0;
    sd = suppress_fft_roundoff(sd,output_support);
    sk = suppress_fft_roundoff(sk,output_support);
    sd_s = suppress_fft_roundoff(sd_s,output_support);
    sk_s = suppress_fft_roundoff(sk_s,output_support);
    sinh_over_a = time*ones(size(q));
    sinh_over_a(nonzero) = ...
        sinh(sqrt_A(nonzero)*time)./sqrt_A(nonzero);
    cosh_factor = cosh(sqrt_A*time);
    coefficient = weights(index)*exp(node)/lambda ...
        *exp(-2*balance_shift*time)/4;
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

ux = ifft2(spectrum.*multipliers.unit_x);
uy = ifft2(spectrum.*multipliers.unit_y);
uxx = ifft2(spectrum.*multipliers.unit_x2);
uxy = ifft2(spectrum.*multipliers.unit_xy);
uyy = first.eta-uxx;
nu_u = 1i*first.phiz;
nu_ux = ifft2(spectrum.*multipliers.nu_unit_x);
nu_uy = ifft2(spectrum.*multipliers.nu_unit_y);
nu_uxx = ifft2(spectrum.*multipliers.nu_unit_x2);
nu_uxy = ifft2(spectrum.*multipliers.nu_unit_xy);
nu_uyy = nu_u-nu_uxx;
count = count+8;
g2_source = 0.25*first.eta.^2+0.5*(ux.^2+uy.^2) ...
    +0.25*(uxx.^2+2*uxy.^2+uyy.^2);
g2_source_s = 0.5*first.eta.*nu_u ...
    +ux.*nu_ux+uy.*nu_uy ...
    +0.5*uxx.*nu_uxx+uxy.*nu_uxy+0.5*uyy.*nu_uyy;
products = products+18;
g2_spectrum = fft2(g2_source);
g2_s_spectrum = fft2(g2_source_s);
count = count+2;
eta_spectrum = eta_spectrum+g2_spectrum.*eta_correction;
if compute_psi
    etat_spectrum = etat_spectrum-1i*g2_s_spectrum.*eta_correction;
end
phi_spectrum = phi_spectrum+g2_spectrum.*phi_correction;
phit_spectrum = phit_spectrum ...
    -1i*g2_s_spectrum.*phi_correction;
eta_correction_spectrum = g2_spectrum.*eta_correction;
eta_spectrum(~output_support) = 0;
etat_spectrum(~output_support) = 0;
phi_spectrum(~output_support) = 0;
phit_spectrum(~output_support) = 0;

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
    etat_spectrum = suppress_fft_roundoff(etat_spectrum,output_support);
    pair.eta_t = ifft2(etat_spectrum);
    pair.phiz_t = pair.phitz;
    count = count+1;
end
diagnostics = struct( ...
    'eta_core_spectrum_max',max(abs(eta_core_spectrum),[],'all'), ...
    'eta_correction_spectrum_max', ...
        max(abs(eta_correction_spectrum),[],'all'));
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

function [eta_correction,phi_correction] = ...
    pair_stokes_corrections(q,nodes,weights,lambda)
eta_correction = complex(zeros(size(q)));
phi_correction = complex(zeros(size(q)));
active = q>0;
if ~any(active,'all'), return, end
qa = q(active);
nu = sqrt(qa.*tanh(qa));
Q = 2*qa;
A = Q.*tanh(Q);
a = sqrt(A);
s = 2*nu;
sd = 3*nu.^2-qa.^2./nu.^2;
sk = 4*qa.^2./nu;
eta_core = zeros(size(qa));
phi_core = complex(zeros(size(qa)));
for index = 1:4
    node = nodes(index);
    time = node/lambda;
    exp_cosh = 0.5*(exp(-(s-a)*time)+exp(-(s+a)*time));
    exp_sinh_over_a = 0.5*( ...
        exp(-(s-a)*time)-exp(-(s+a)*time))./a;
    coefficient = weights(index)*exp(node)/lambda/4;
    eta_core = eta_core+coefficient*( ...
        -A.*exp_sinh_over_a.*sd+exp_cosh.*sk);
    phi_core = phi_core+1i*coefficient*( ...
        exp_cosh.*sd-exp_sinh_over_a.*sk);
end
t = tanh(qa);
eta_exact = qa.*(3-t.^2)./(4*t.^3);
detuning = A-s.^2;
phi_exact = -1i*(s.*sd-sk)./detuning/4;
eta_correction(active) = eta_exact-eta_core;
phi_correction(active) = phi_exact-phi_core;
end

function correction = third_order_stokes_correction( ...
    q,nodes,weights,lambda)
correction = zeros(size(q));
active = q>0;
if ~any(active,'all'), return, end
qa = q(active);
nu = sqrt(qa.*tanh(qa));
s = 3*nu;
A = 3*qa.*tanh(3*qa);
a = sqrt(A);
csch_q = 1./sinh(qa);
rk = 9i/8*qa.^2.*(13+32*cosh(2*qa)+3*cosh(4*qa)) ...
    .*csch_q.^4.*nu;
rd = -9/8*qa.^2.*(-25+16*cosh(2*qa)+cosh(4*qa)) ...
    .*csch_q.^4;
core = complex(zeros(size(qa)));
for index = 1:4
    node = nodes(index);
    time = node/lambda;
    exp_cosh = 0.5*(exp(-(s-a)*time)+exp(-(s+a)*time));
    exp_sinh = 0.5*(exp(-(s-a)*time)-exp(-(s+a)*time));
    coefficient = weights(index)*exp(node)/lambda/24;
    core = core+coefficient*(a.*exp_sinh.*rd-1i*exp_cosh.*rk);
end
stokes = 3/256*qa.^2.*(14+15*cosh(2*qa) ...
    +6*cosh(4*qa)+cosh(6*qa)).*csch_q.^6;
correction(active) = real(stokes-core);
end

function psi_correction = third_order_surface_stokes_correction( ...
    q,nodes,weights,lambda)
% Exact surface-Psi33 diagonal minus the unchanged GL4 K1 response.
psi_correction = complex(zeros(size(q)));
active = q>0;
if ~any(active,'all'), return, end
qa = q(active);
nu = sqrt(qa.*tanh(qa));
s = 3*nu;
A2 = 2*qa.*tanh(2*qa);
A3 = 3*qa.*tanh(3*qa);
a3 = sqrt(A3);
csch_q = 1./sinh(qa);
eta22 = qa.*coth(qa).^3.*(3-tanh(qa).^2)/4;
phi22 = -3i/8*cosh(2*qa).*csch_q.^4.*nu;
phi33 = 1i/64*qa.*(-11+2*cosh(2*qa)) ...
    .*(-1+2*cosh(2*qa)).*coth(qa).*csch_q.^6.*nu;
contact = 0.5*A2.*phi22-0.5i*eta22.*nu ...
    -0.125i*qa.^2./nu;
psi_exact = phi33+contact;
rk = 9i/8*qa.^2.*(13+32*cosh(2*qa)+3*cosh(4*qa)) ...
    .*csch_q.^4.*nu;
rd = -9/8*qa.^2.*(-25+16*cosh(2*qa)+cosh(4*qa)) ...
    .*csch_q.^4;
fk_psi = rk/24+A3.*contact;
fd_psi = rd/24+1i*s.*contact;
psi_core = complex(zeros(size(qa)));
for index = 1:4
    node = nodes(index);
    time = node/lambda;
    exp_cosh = 0.5*(exp(-(s-a3)*time)+exp(-(s+a3)*time));
    exp_sinh_over_a = 0.5*( ...
        exp(-(s-a3)*time)-exp(-(s+a3)*time))./a3;
    coefficient = weights(index)*exp(node)/lambda;
    psi_core = psi_core+coefficient*( ...
        -1i*exp_cosh.*fd_psi-exp_sinh_over_a.*fk_psi);
end
psi_correction(active) = psi_exact-psi_core;
end

function value = cached_third_order_surface_stokes_correction( ...
    q,nodes,weights,lambda)
% Keep the eta33-only stable path free of Psi33 initialization work.
persistent cached_q cached_nodes cached_weights cached_lambda cached_value
if isempty(cached_value) || ~isequal(cached_q,q) ...
        || ~isequal(cached_nodes,nodes) || ~isequal(cached_weights,weights) ...
        || ~isequal(cached_lambda,lambda)
    cached_value = third_order_surface_stokes_correction( ...
        q,nodes,weights,lambda);
    cached_q = q;
    cached_nodes = nodes;
    cached_weights = weights;
    cached_lambda = lambda;
end
value = cached_value;
end

function [source,count,products] = ...
    triple_gate_source(spectrum,multipliers,gate)
moments = cell(3,3);
moments{1,1} = ifft2(spectrum);
moments{2,1} = ifft2(spectrum.*multipliers.unit_x);
moments{1,2} = ifft2(spectrum.*multipliers.unit_y);
moments{3,1} = ifft2(spectrum.*multipliers.unit_x2);
moments{2,2} = ifft2(spectrum.*multipliers.unit_xy);
moments{1,3} = moments{1,1}-moments{3,1};
count = 5;
source = complex(zeros(size(spectrum)));
terms = gate.separated_terms;
products = 0;
for index = 1:numel(terms)
    p1 = terms(index).leaf1_exponents;
    p2 = terms(index).leaf2_exponents;
    p3 = terms(index).leaf3_exponents;
    coefficient = terms(index).coefficient_numerator ...
        /terms(index).coefficient_denominator;
    source = source+coefficient ...
        *moments{p1(1)+1,p1(2)+1} ...
        .*moments{p2(1)+1,p2(2)+1} ...
        .*moments{p3(1)+1,p3(2)+1};
    products = products+2;
end
end

function [plan,cache_hit] = cached_execution_plan( ...
    qx,qy,peak_depth_qp,project_root)
persistent cached_plan cached_qx cached_qy cached_peak cached_root
cache_hit = ~isempty(cached_plan) ...
    && isequal(cached_qx,qx) ...
    && isequal(cached_qy,qy) ...
    && isequal(cached_peak,peak_depth_qp) ...
    && isequal(cached_root,project_root);
if cache_hit
    plan = cached_plan;
    return
end
[freeze,gate] = load_and_verify_interfaces(project_root);
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
[nodes,weights] = gauss_laguerre_four();
lambda2 = 2*sqrt(peak_depth_qp*tanh(peak_depth_qp)) ...
    -sqrt(2*peak_depth_qp*tanh(2*peak_depth_qp));
lambda3 = 3*sqrt(peak_depth_qp*tanh(peak_depth_qp)) ...
    -sqrt(3*peak_depth_qp*tanh(3*peak_depth_qp));
[pair_eta_correction,pair_phi_correction] = pair_stokes_corrections( ...
    q/2,nodes,weights,lambda2);
third_order_correction = third_order_stokes_correction( ...
    q/3,nodes,weights,lambda3);
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
    'gate',gate, ...
    'q',q, ...
    'nu',nu, ...
    'A',A, ...
    'sqrt_A',sqrt_A, ...
    'safe_nu',safe_nu, ...
    'nonzero',nonzero, ...
    'unit_x',unit_x, ...
    'unit_y',unit_y, ...
    'nodes',nodes, ...
    'weights',weights, ...
    'lambda2',lambda2, ...
    'lambda3',lambda3, ...
    'pair_eta_correction',pair_eta_correction, ...
    'pair_phi_correction',pair_phi_correction, ...
    'third_order_correction',third_order_correction);
plan.multipliers = multipliers;
cached_plan = plan;
cached_qx = qx;
cached_qy = qy;
cached_peak = peak_depth_qp;
cached_root = project_root;
end

function [freeze,gate] = load_and_verify_interfaces(project_root)
freeze_path = fullfile(project_root,'symbolic','generated', ...
    'finite_depth_directional_order3_nested_green_laplace.json');
gate_path = fullfile(project_root,'symbolic','generated', ...
    'finite_depth_directional_order3_crossing_stokes_gate.json');
if ~isfile(freeze_path) || ~isfile(gate_path)
    error('Frozen Wolfram nested eta33 interfaces are missing.');
end
freeze = jsondecode(fileread(freeze_path));
gate = jsondecode(fileread(gate_path));
if ~freeze.overall_exact_gate_pass || freeze.oracle_or_mf12_used ...
        || freeze.input_chi_route_used ...
        || freeze.inner_quadrature_rank~=4 ...
        || freeze.outer_quadrature_rank~=4
    error('Nested eta33 freeze crossed its information boundary.');
end
if ~gate.overall_exact_gate_pass || gate.oracle_or_mf12_used ...
        || gate.input_chi_route_used || gate.separated_term_count~=27
    error('Crossing Stokes gate interface is invalid.');
end
end

function [nodes,weights] = gauss_laguerre_four()
persistent cached_nodes cached_weights
if isempty(cached_nodes)
    cached_nodes = sort(real(roots([1/24,-2/3,3,-4,1]))).';
    laguerre5 = [-1/120,5/24,-5/3,5,-5,1];
    cached_weights = cached_nodes ...
        ./(25*polyval(laguerre5,cached_nodes).^2);
end
nodes = cached_nodes;
weights = cached_weights;
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
