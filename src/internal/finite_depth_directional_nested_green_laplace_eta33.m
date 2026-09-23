function [eta33_plus,audit,components,state33,state12] = ...
    finite_depth_directional_nested_green_laplace_eta33( ...
    eta11_spectrum,qx,qy,peak_depth_qp,project_root,options)
%FINITE_DEPTH_DIRECTIONAL_NESTED_GREEN_LAPLACE_ETA33 Frozen eta11-only graph.
% The no-Stokes inner GL response generates eta22, flat Phi22, and d_t Phi22 from
% damped eta11 leaves. The original order-three kinematic/dynamic forcing is
% then passed through the no-Stokes outer GL response at the same rank.
% The optional fourth output is an internal lower-state package for order
% four: flat Phi33, d_t Phi33, and the required horizontal/vertical banks.
% Requesting the optional fifth output returns the internally generated
% first-/second-order state needed by the retained GL--Partition-G4 eta44
% executor.  These flat-potential fields are internal dependencies, not
% public Phi candidate outputs.

arguments
    eta11_spectrum (:,:) {mustBeNumeric}
    qx (:,:) double
    qy (:,:) double
    peak_depth_qp (1,1) double {mustBePositive}
    project_root (1,1) string = string( ...
        fileparts(fileparts(fileparts(mfilename('fullpath')))))
    options.include_eta_time (1,1) logical = false
    options.quadrature_rank (1,1) double {mustBeInteger,mustBePositive} = 4
end
if ~ismember(options.quadrature_rank,[4,6,8,10,12])
    error('Unsupported no-Stokes nested Green--Laplace rank.');
end
if ~isequal(size(eta11_spectrum),size(qx),size(qy))
    error('Spectrum, qx, and qy must have identical sizes.');
end
[plan,plan_cache_hit] = cached_execution_plan( ...
    qx,qy,peak_depth_qp,options.quadrature_rank,project_root);
freeze = plan.freeze;
q = plan.q;
nu = plan.nu;
A = plan.A;
sqrt_A = plan.sqrt_A;
nonzero = plan.nonzero;
multipliers = plan.multipliers;

threshold = 1e-13*max(1,max(abs(eta11_spectrum(:))));
support = abs(eta11_spectrum)>threshold;
if ~any(support,'all') || any(qx(support)<=0) || ...
        any(q(support)<=0)
    error('Strict-forward eta11 support with every parent q>0 is required.');
end
assert_cubic_support_is_alias_safe(support);
pair_support = real(ifft2(fft2(double(support)).^2))>0.5;
triple_support = real(ifft2(fft2(double(support)).^3))>0.5;
minimum_parent_nu=min(nu(support));
maximum_pair_a=max(sqrt_A(pair_support));
maximum_triple_a=max(sqrt_A(triple_support));
inner_balance_shift=0.5*(minimum_parent_nu+maximum_pair_a/2);
outer_balance_shift=0.5*(minimum_parent_nu+maximum_triple_a/3);

nodes = plan.nodes;
weights = plan.weights;
lambda2 = plan.lambda2;
lambda3 = plan.lambda3;
if ~isfinite(lambda2) || ~isfinite(lambda3) ...
        || lambda2<=0 || lambda3<=0
    error('Frozen Green--Laplace tail scales are invalid.');
end
if use_cached_node_factors()
    inner_node_factors=cached_node_factor_bank( ...
        nu,sqrt_A,nodes,weights,lambda2,inner_balance_shift,2);
    outer_node_factors=cached_node_factor_bank( ...
        nu,sqrt_A,nodes,weights,lambda3,outer_balance_shift,3);
else
    inner_node_factors=[];
    outer_node_factors=[];
end

candidate_spectrum = complex(zeros(size(q)));
if options.include_eta_time
    candidate_t_spectrum = complex(zeros(size(q)));
else
    candidate_t_spectrum = [];
end
phi_spectrum = complex(zeros(size(q)));
phit_spectrum = complex(zeros(size(q)));
fft_count = 0;
product_count = 0;

for outer_index = 1:options.quadrature_rank
    node = nodes(outer_index);
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
    [pair,pair_fft,pair_products] = inner_pair_fields( ...
        damped_spectrum,q,A,sqrt_A,multipliers, ...
        pair_support,nodes,weights,lambda2,inner_balance_shift, ...
        inner_node_factors,"full");
    fft_count = fft_count+first_fft+pair_fft;
    product_count = product_count+pair_products;

    [forcing_k,forcing_d,forcing_products] = ...
        order3_forcing(first,pair);
    [forcing_k_t,forcing_d_t,forcing_t_products] = ...
        order3_forcing_time(first,pair);
    product_count = product_count+forcing_products+forcing_t_products;
    forcing_k_spectrum = fft2(forcing_k);
    forcing_d_spectrum = fft2(forcing_d);
    forcing_k_t_spectrum = fft2(forcing_k_t);
    forcing_d_t_spectrum = fft2(forcing_d_t);
    fft_count = fft_count+4;
    forcing_k_spectrum(~triple_support) = 0;
    forcing_d_spectrum(~triple_support) = 0;
    forcing_k_t_spectrum(~triple_support) = 0;
    forcing_d_t_spectrum(~triple_support) = 0;
    forcing_k_spectrum=suppress_fft_roundoff( ...
        forcing_k_spectrum,triple_support);
    forcing_d_spectrum=suppress_fft_roundoff( ...
        forcing_d_spectrum,triple_support);
    forcing_k_t_spectrum=suppress_fft_roundoff( ...
        forcing_k_t_spectrum,triple_support);
    forcing_d_t_spectrum=suppress_fft_roundoff( ...
        forcing_d_t_spectrum,triple_support);
    if isempty(outer_node_factors)
        sinh_over_a = time*ones(size(q));
        sinh_over_a(nonzero) = ...
            sinh(sqrt_A(nonzero)*time)./sqrt_A(nonzero);
        cosh_factor=cosh(sqrt_A*time);
        coefficient = weights(outer_index)*exp(node)/lambda3 ...
            *exp(-3*outer_balance_shift*time);
    else
        sinh_over_a=outer_node_factors.sinh_over_a(:,:,outer_index);
        cosh_factor=outer_node_factors.cosh_factor(:,:,outer_index);
        coefficient=outer_node_factors.coefficient(outer_index);
    end
    candidate_spectrum = candidate_spectrum+coefficient*( ...
        A.*sinh_over_a.*forcing_d_spectrum ...
        -1i*cosh_factor.*forcing_k_spectrum);
    if options.include_eta_time
        candidate_t_spectrum = candidate_t_spectrum+coefficient*( ...
            A.*sinh_over_a.*forcing_d_t_spectrum ...
            -1i*cosh_factor.*forcing_k_t_spectrum);
    end
    phi_spectrum = phi_spectrum+coefficient*( ...
        -1i*cosh_factor.*forcing_d_spectrum ...
        -sinh_over_a.*forcing_k_spectrum);
    phit_spectrum = phit_spectrum+coefficient*( ...
        -1i*cosh_factor.*forcing_d_t_spectrum ...
        -sinh_over_a.*forcing_k_t_spectrum);
end

raw_resolvent_spectrum = candidate_spectrum;
crossing_correction_spectrum = complex(zeros(size(q)));
candidate_spectrum(~triple_support) = 0;
if options.include_eta_time,candidate_t_spectrum(~triple_support) = 0;end
phi_spectrum(~triple_support) = 0;
phit_spectrum(~triple_support) = 0;
eta33_plus = ifft2(candidate_spectrum);
state33 = third_order_state_bank(eta33_plus,candidate_spectrum, ...
    phi_spectrum,phit_spectrum,multipliers,~use_eta44_root_pruning());
fft_count = fft_count+1+state33.fft_ifft_count;
state33 = rmfield(state33,'fft_ifft_count');
if options.include_eta_time
    state33.eta_t = ifft2(candidate_t_spectrum);
    fft_count = fft_count+1;
end
state12 = struct();
if nargout>=5
    if use_eta44_root_pruning()
        [first12,first12_fft] = first_order_state_eta44_root( ...
            eta11_spectrum,multipliers,plan.safe_nu,nonzero);
        [pair12,pair12_fft,pair12_products] = inner_pair_fields( ...
            eta11_spectrum,q,A,sqrt_A,multipliers, ...
            pair_support,nodes,weights,lambda2,inner_balance_shift, ...
            inner_node_factors,"eta44_root");
        extra_fft_count=3;
    else
        [first12,first12_fft] = first_order_state( ...
            eta11_spectrum,multipliers);
        [pair12,pair12_fft,pair12_products] = inner_pair_fields( ...
            eta11_spectrum,q,A,sqrt_A,multipliers, ...
            pair_support,nodes,weights,lambda2,inner_balance_shift, ...
            inner_node_factors,"full");
        phi1_spectrum=-1i*eta11_spectrum./plan.safe_nu;
        phi1_spectrum(~nonzero)=0;
        first12.phixzz=ifft2(multipliers.dx.*multipliers.q2.*phi1_spectrum);
        first12.phiyzz=ifft2(multipliers.dy.*multipliers.q2.*phi1_spectrum);
        first12.phizzzz=ifft2(multipliers.q2.^2.*phi1_spectrum);
        extra_fft_count=6;
    end
    pair12.phixz = ifft2(multipliers.dx.*A.*pair12.Phi_spectrum);
    pair12.phiyz = ifft2(multipliers.dy.*A.*pair12.Phi_spectrum);
    pair12.phizzz = ifft2(multipliers.q2.*A.*pair12.Phi_spectrum);
    state12 = struct('first',first12,'pair',pair12);
    fft_count = fft_count+first12_fft+pair12_fft+extra_fft_count;
    product_count = product_count+pair12_products;
end
if any(~isfinite(eta33_plus),'all') ...
        || (isfield(state33,'Phi') && any(~isfinite(state33.Phi),'all')) ...
        || (isfield(state33,'Phi_t') && any(~isfinite(state33.Phi_t),'all'))
    error('Nested Green--Laplace eta33/Phi33 state is not finite.');
end

audit = struct( ...
    'candidate_id',string(freeze.candidate_id), ...
    'status','internal-lower-state-backend-for-eta44', ...
    'information_boundary','eta11-only', ...
    'internal_generated_fields',{{'eta22_plus','Phi22_plus', ...
        'Phi22_t_plus','eta33_plus','Phi33_plus','Phi33_t_plus'}}, ...
    'input_chi_route_used',false, ...
    'oracle_or_mf12_used',false, ...
    'inner_quadrature_rank',options.quadrature_rank, ...
    'outer_quadrature_rank',options.quadrature_rank, ...
    'inner_scale',lambda2, ...
    'outer_scale',lambda3, ...
    'inner_balance_shift',inner_balance_shift, ...
    'outer_balance_shift',outer_balance_shift, ...
    'numerical_stability_policy', ...
        'balanced separated exponentials plus FFT backward-error suppression', ...
    'static_plan_cache_hit',plan_cache_hit, ...
    'stokes_correction_used',false, ...
    'crossing_gate','none', ...
    'fft_ifft_count',fft_count, ...
    'pointwise_product_count',product_count, ...
    'pair_loops',0, ...
    'triple_loops',0, ...
    'cost_class','fixed-constant N log N');
components = struct( ...
    'raw_resolvent_spectrum',raw_resolvent_spectrum, ...
    'crossing_correction_spectrum',crossing_correction_spectrum, ...
    'total_spectrum',candidate_spectrum, ...
    'Phi33_spectrum',phi_spectrum, ...
    'Phi33_t_spectrum',phit_spectrum);
end

function [forcing_k,forcing_d,products] = order3_forcing(first,pair)
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
products = 19;
end

function [forcing_k_t,forcing_d_t,products] = ...
    order3_forcing_time(first,pair)
forcing_k_pair_t = ...
    first.phix_t.*pair.etax+first.phix.*pair.etax_t ...
    +first.phiy_t.*pair.etay+first.phiy.*pair.etay_t ...
    +pair.phix_t.*first.etax+pair.phix.*first.etax_t ...
    +pair.phiy_t.*first.etay+pair.phiy.*first.etay_t ...
    -first.eta_t.*pair.phizz-first.eta.*pair.phizz_t ...
    -pair.eta_t.*first.phizz-pair.eta.*first.phizz_t;
forcing_k_direct_t = first.eta_t.*( ...
    first.phixz.*first.etax+first.phiyz.*first.etay) ...
    +first.eta.*(first.phixz_t.*first.etax ...
    +first.phixz.*first.etax_t ...
    +first.phiyz_t.*first.etay+first.phiyz.*first.etay_t) ...
    -first.eta.*first.eta_t.*first.phizzz ...
    -0.5*first.eta.^2.*first.phizzz_t;
forcing_d_pair_t = ...
    first.eta_t.*pair.phitz+first.eta.*pair.phitz_t ...
    +pair.eta_t.*first.phitz+pair.eta.*first.phitz_t ...
    +first.phix_t.*pair.phix+first.phix.*pair.phix_t ...
    +first.phiy_t.*pair.phiy+first.phiy.*pair.phiy_t ...
    +first.phiz_t.*pair.phiz+first.phiz.*pair.phiz_t;
forcing_d_direct_t = first.eta.*first.eta_t.*first.phitzz ...
    +0.5*first.eta.^2.*first.phitzz_t ...
    +first.eta_t.*(first.phix.*first.phixz ...
    +first.phiy.*first.phiyz+first.phiz.*first.phizz) ...
    +first.eta.*(first.phix_t.*first.phixz ...
    +first.phix.*first.phixz_t ...
    +first.phiy_t.*first.phiyz+first.phiy.*first.phiyz_t ...
    +first.phiz_t.*first.phizz+first.phiz.*first.phizz_t);
forcing_k_t = forcing_k_direct_t/4+forcing_k_pair_t/2;
forcing_d_t = forcing_d_direct_t/4+forcing_d_pair_t/2;
products = 38;
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
state.etax_t = state.phixz;
state.etay_t = state.phiyz;
state.phix_t = -state.etax;
state.phiy_t = -state.etay;
state.phiz_t = state.phitz;
state.phizz_t = state.phitzz;
state.phitzz_t = -state.phizzz;
state.phixz_t = ifft2(multipliers.minus_i_nu2_qx.*spectrum);
state.phiyz_t = ifft2(multipliers.minus_i_nu2_qy.*spectrum);
state.phizzz_t = ifft2(multipliers.minus_q2_A.*spectrum);
state.phitz_t = ifft2(multipliers.i_A_nu.*spectrum);
count = 16;
end

function [state,count] = first_order_state_eta44_root( ...
        spectrum,multipliers,safe_nu,nonzero)
state=struct();
state.eta=ifft2(spectrum);
state.etax=ifft2(multipliers.dx.*spectrum);
state.etay=ifft2(multipliers.dy.*spectrum);
state.phix=ifft2(multipliers.qx_over_nu.*spectrum);
state.phiy=ifft2(multipliers.qy_over_nu.*spectrum);
state.phiz=ifft2(multipliers.minus_i_nu.*spectrum);
state.phixz=ifft2(multipliers.nu_qx.*spectrum);
state.phiyz=ifft2(multipliers.nu_qy.*spectrum);
state.phizz=ifft2(multipliers.minus_i_q2_over_nu.*spectrum);
state.phizzz=ifft2(multipliers.minus_i_q2_nu.*spectrum);
state.phitz=ifft2(multipliers.minus_A.*spectrum);
state.phitzz=ifft2(multipliers.minus_q2.*spectrum);
state.phizzz_t=ifft2(multipliers.minus_q2_A.*spectrum);
phi1_spectrum=-1i*spectrum./safe_nu;
phi1_spectrum(~nonzero)=0;
state.phixzz=ifft2(multipliers.dx.*multipliers.q2.*phi1_spectrum);
state.phiyzz=ifft2(multipliers.dy.*multipliers.q2.*phi1_spectrum);
state.phizzzz=ifft2(multipliers.q2.^2.*phi1_spectrum);
count=16;
end

function [pair,count,products] = inner_pair_fields( ...
    spectrum,q,A,sqrt_A,multipliers,output_support, ...
    nodes,weights,lambda,balance_shift,node_factors,output_mode)
eta_spectrum = complex(zeros(size(q)));
phi_spectrum = complex(zeros(size(q)));
phit_spectrum = complex(zeros(size(q)));
eta44_root=output_mode=="eta44_root";
if eta44_root
    etat_spectrum=[];
    phitt_spectrum=[];
else
    etat_spectrum = complex(zeros(size(q)));
    phitt_spectrum = complex(zeros(size(q)));
end
count = 0;
products = 0;
nonzero = q>0;
for index = 1:numel(nodes)
    node = nodes(index);
    time = node/lambda;
    if isempty(node_factors)
        damped = spectrum.*exp( ...
            -time*(multipliers.nu-balance_shift));
    else
        damped=spectrum.*node_factors.damping(:,:,index);
    end
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
    if ~eta44_root
        v_nu4 = ifft2(damped.*multipliers.nu4);
        q2_nu = ifft2(damped.*multipliers.q2_nu);
        nu2_jx = ifft2(damped.*multipliers.nu2_qx);
        nu2_jy = ifft2(damped.*multipliers.nu2_qy);
        count = count+16;
    else
        count = count+12;
    end

    source_d = 2*v.*v_nu2+v_nu.^2-hx.^2-hy.^2;
    source_k = 2*v.*radial_over_nu+2*(hx.*jx+hy.*jy);
    source_d_s = 2*v.*v_nu3+4*v_nu.*v_nu2 ...
        -2*hx.*jx-2*hy.*jy;
    source_k_s = 2*(v_nu.*radial_over_nu+v.*q2_field) ...
        +2*(jx.^2+hx.*nu_jx+jy.^2+hy.*nu_jy);
    if ~eta44_root
        source_d_ss = 2*v.*v_nu4+6*v_nu.*v_nu3 ...
            +4*v_nu2.^2-2*(jx.^2+hx.*nu_jx ...
            +jy.^2+hy.*nu_jy);
        source_k_ss = 2*(v_nu2.*radial_over_nu ...
            +2*v_nu.*q2_field+v.*q2_nu) ...
            +6*(jx.*nu_jx+jy.*nu_jy) ...
            +2*(hx.*nu2_jx+hy.*nu2_jy);
        products = products+34;
    else
        products = products+18;
    end
    sd = fft2(source_d);
    sk = fft2(source_k);
    sd_s = fft2(source_d_s);
    sk_s = fft2(source_k_s);
    if ~eta44_root
        sd_ss = fft2(source_d_ss);
        sk_ss = fft2(source_k_ss);
        count = count+6;
    else
        count = count+4;
    end
    sd(~output_support) = 0;
    sk(~output_support) = 0;
    sd_s(~output_support) = 0;
    sk_s(~output_support) = 0;
    sd=suppress_fft_roundoff(sd,output_support);
    sk=suppress_fft_roundoff(sk,output_support);
    sd_s=suppress_fft_roundoff(sd_s,output_support);
    sk_s=suppress_fft_roundoff(sk_s,output_support);
    if ~eta44_root
        sd_ss(~output_support) = 0;
        sk_ss(~output_support) = 0;
        sd_ss=suppress_fft_roundoff(sd_ss,output_support);
        sk_ss=suppress_fft_roundoff(sk_ss,output_support);
    end
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
    if ~eta44_root
        etat_spectrum = etat_spectrum-1i*coefficient*( ...
            -A.*sinh_over_a.*sd_s+cosh_factor.*sk_s);
    end
    phi_spectrum = phi_spectrum+1i*coefficient*( ...
        cosh_factor.*sd-sinh_over_a.*sk);
    phit_spectrum = phit_spectrum+coefficient*( ...
        cosh_factor.*sd_s-sinh_over_a.*sk_s);
    if ~eta44_root
        phitt_spectrum = phitt_spectrum-1i*coefficient*( ...
            cosh_factor.*sd_ss-sinh_over_a.*sk_ss);
    end
end

eta_spectrum(~output_support) = 0;
phi_spectrum(~output_support) = 0;
phit_spectrum(~output_support) = 0;
if ~eta44_root
    etat_spectrum(~output_support) = 0;
    phitt_spectrum(~output_support) = 0;
end

pair.eta = ifft2(eta_spectrum);
pair.etax = ifft2(multipliers.dx.*eta_spectrum);
pair.etay = ifft2(multipliers.dy.*eta_spectrum);
pair.phix = ifft2(multipliers.dx.*phi_spectrum);
pair.phiy = ifft2(multipliers.dy.*phi_spectrum);
pair.phiz = ifft2(A.*phi_spectrum);
pair.phizz = ifft2(multipliers.q2.*phi_spectrum);
pair.phitz = ifft2(A.*phit_spectrum);
pair.phizz_t = ifft2(multipliers.q2.*phit_spectrum);
if output_mode=="full"
    pair.eta_t = ifft2(etat_spectrum);
    pair.etax_t = ifft2(multipliers.dx.*etat_spectrum);
    pair.etay_t = ifft2(multipliers.dy.*etat_spectrum);
    pair.phix_t = ifft2(multipliers.dx.*phit_spectrum);
    pair.phiy_t = ifft2(multipliers.dy.*phit_spectrum);
    pair.phiz_t = pair.phitz;
    pair.phitz_t = ifft2(A.*phitt_spectrum);
    count=count+16;
else
    count=count+9;
end
pair.Phi_spectrum = phi_spectrum;
end

function spectrum = suppress_fft_roundoff(spectrum,support)
active_values=abs(spectrum(support));
if isempty(active_values),return,end
scale=max(active_values);
if scale==0,return,end
relative_floor=16*eps*max(1,log2(numel(spectrum)));
spectrum(support & abs(spectrum)<relative_floor*scale)=0;
end

function state = third_order_state_bank( ...
    eta,eta_spectrum,phi_spectrum,phit_spectrum,multipliers,include_traces)
state = struct();
state.eta = eta;
state.etax = ifft2(multipliers.dx.*eta_spectrum);
state.etay = ifft2(multipliers.dy.*eta_spectrum);
if include_traces
    state.Phi = ifft2(phi_spectrum);
    state.Phi_t = ifft2(phit_spectrum);
end
state.phix = ifft2(multipliers.dx.*phi_spectrum);
state.phiy = ifft2(multipliers.dy.*phi_spectrum);
state.phiz = ifft2(multipliers.A.*phi_spectrum);
state.phizz = ifft2(multipliers.q2.*phi_spectrum);
state.phizt = ifft2(multipliers.A.*phit_spectrum);
state.fft_ifft_count = 7+2*include_traces;
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

function enabled=use_eta44_root_pruning()
persistent cached_raw cached_enabled
raw=string(getenv('GL_MATLAB_GL_PRUNE_ETA44_ROOT'));
if isempty(cached_raw) || raw~=cached_raw
    cached_raw=raw;
    cached_enabled=ismember(lower(strtrim(raw)),["1","true","yes","on"]);
end
enabled=cached_enabled;
end

function factors=cached_node_factor_bank( ...
        nu,sqrt_A,nodes,weights,lambda,shift,order)
persistent cache
if isempty(cache),cache=cell(1,4);end
entry=cache{order};
hit=~isempty(entry) && isequal(entry.nu,nu) ...
    && isequal(entry.nodes,nodes) && isequal(entry.weights,weights) ...
    && isequal(entry.lambda,lambda) && isequal(entry.shift,shift);
if hit
    factors=entry.factors;
    return
end
page_count=numel(nodes);shape=[size(nu),page_count];
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
cache{order}=struct('nu',nu,'nodes',nodes,'weights',weights, ...
    'lambda',lambda,'shift',shift,'factors',factors);
end

function [plan,cache_hit] = cached_execution_plan( ...
    qx,qy,peak_depth_qp,quadrature_rank,project_root)
persistent cached_plan cached_qx cached_qy cached_peak cached_rank cached_root
cache_hit = ~isempty(cached_plan) ...
    && isequal(cached_qx,qx) ...
    && isequal(cached_qy,qy) ...
    && isequal(cached_peak,peak_depth_qp) ...
    && isequal(cached_rank,quadrature_rank) ...
    && isequal(cached_root,project_root);
if cache_hit
    plan = cached_plan;
    return
end
freeze = load_and_verify_interface(project_root);
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
[nodes,weights] = gauss_laguerre_rule(quadrature_rank);
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
    'A',A, ...
    'q2',q2, ...
    'nu',nu, ...
    'nu2',nu.^2, ...
    'nu3',nu.^3, ...
    'nu4',nu.^4, ...
    'dx',1i*qx, ...
    'dy',1i*qy, ...
    'qx_over_nu',qx./safe_nu, ...
    'qy_over_nu',qy./safe_nu, ...
    'q2_over_nu',q2./safe_nu, ...
    'minus_i_nu',-1i*nu, ...
    'nu_qx',nu.*qx, ...
    'nu_qy',nu.*qy, ...
    'nu2_qx',nu.^2.*qx, ...
    'nu2_qy',nu.^2.*qy, ...
    'q2_nu',q2.*nu, ...
    'minus_i_q2_over_nu',-1i*q2./safe_nu, ...
    'minus_i_q2_nu',-1i*q2.*nu, ...
    'minus_A',-A, ...
    'minus_q2',-q2, ...
    'minus_i_nu2_qx',-1i*nu.^2.*qx, ...
    'minus_i_nu2_qy',-1i*nu.^2.*qy, ...
    'minus_q2_A',-q2.*A, ...
    'i_A_nu',1i*A.*nu, ...
    'unit_x',unit_x, ...
    'unit_y',unit_y, ...
    'unit_x2',unit_x2, ...
    'unit_xy',unit_xy, ...
    'nu_unit_x',nu.*unit_x, ...
    'nu_unit_y',nu.*unit_y, ...
    'nu_unit_x2',nu.*unit_x2, ...
    'nu_unit_xy',nu.*unit_xy, ...
    'nu2_unit_x',nu.^2.*unit_x, ...
    'nu2_unit_y',nu.^2.*unit_y, ...
    'nu2_unit_x2',nu.^2.*unit_x2, ...
    'nu2_unit_xy',nu.^2.*unit_xy);
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
    'nodes',nodes, ...
    'weights',weights, ...
    'lambda2',lambda2, ...
    'lambda3',lambda3);
plan.multipliers = multipliers;
cached_plan = plan;
cached_qx = qx;
cached_qy = qy;
cached_peak = peak_depth_qp;
cached_rank = quadrature_rank;
cached_root = project_root;
end

function freeze = load_and_verify_interface(project_root)
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
required_states = ["Phi33_plus","Phi33_t_plus"];
if ~all(ismember(required_states,string(freeze.internal_generated_fields)))
    error('Frozen eta33 interface does not declare the order-four state bank.');
end
end

function [nodes,weights] = gauss_laguerre_rule(rank)
indices=1:rank;
jacobi=diag(2*indices-1)+diag(1:rank-1,1)+diag(1:rank-1,-1);
[vectors,values]=eig(jacobi,'vector');
[nodes,order]=sort(real(values));
vectors=vectors(:,order);
weights=real(vectors(1,:)).^2;
nodes=nodes(:).';weights=weights(:).';
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
