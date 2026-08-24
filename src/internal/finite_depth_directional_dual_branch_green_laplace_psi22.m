function [Psi22_plus,audit,components] = ...
    finite_depth_directional_dual_branch_green_laplace_psi22( ...
    eta11_spectrum,qx,qy,peak_depth_qp,project_root, ...
    allow_below_model_parent_domain,exponential_balance_mode, ...
    correction_mode)
%FINITE_DEPTH_DIRECTIONAL_DUAL_BRANCH_GREEN_LAPLACE_PSI22 Direct Psi22.
% The eta11-only stationary resolvent is split into its slow and fast
% detuning branches.  Each branch uses a frozen two-node Gauss--Laguerre
% rule and its own finite-depth peak scale; the scales tend to the exact
% deep-water endpoints 2-sqrt(2) and 2+sqrt(2), respectively.
arguments
    eta11_spectrum (:,:) {mustBeNumeric}
    qx (:,:) double
    qy (:,:) double
    peak_depth_qp (1,1) double {mustBePositive}
    project_root (1,1) string = string( ...
        fileparts(fileparts(fileparts(mfilename('fullpath')))))
    allow_below_model_parent_domain (1,1) logical = false
    exponential_balance_mode (1,1) string = "radial2"
    correction_mode (1,1) string = "none"
end
if ~isequal(size(eta11_spectrum),size(qx),size(qy))
    error('Spectrum, qx, and qy must have identical sizes.');
end
if exponential_balance_mode~="radial2" || correction_mode~="none"
    error(['Only the retained dual-branch radial2/no-angular-correction ', ...
        'Psi22 result is available in the active executor. Historical ', ...
        'ablations are recoverable from the pre-cleanup archive tag.']);
end
freeze = load_and_verify_freeze(project_root);
q = hypot(qx,qy);
nu = sqrt(q.*tanh(q));
sqrt_A = sqrt(q.*tanh(q));
safe_nu = max(nu,realmin);
threshold = 1e-13*max(1,max(abs(eta11_spectrum(:))));
support = abs(eta11_spectrum)>threshold;
if ~any(support,'all') || any(qx(support)<=0) || any(q(support)<=0)
    error('Strict-forward eta11 support with every parent q>0 is required.');
end
below_model_domain = any(q(support)<0.3);
below_model_spectral_mass = sum(abs(eta11_spectrum(q<0.3)).^2,'all') ...
    /max(sum(abs(eta11_spectrum).^2,'all'),realmin);
if below_model_domain && ~allow_below_model_parent_domain
    error(['A parent lies below the declared q=kd>=0.3 model domain. Set ' ...
        'allow_below_model_parent_domain=true only for a labelled ' ...
        'wave-group tail extrapolation.']);
end
assert_quadratic_support_is_alias_safe(support);
support_count = real(ifft2(fft2(double(support)).^2));
output_support = support_count>0.5;

[channel_masks,channel_shifts,channel_pair_supports] = ...
    radial_balance_channels(support,nu,2);
channel_count = numel(channel_masks);

[nodes,weights] = gauss_laguerre_two();
peak_nu = sqrt(peak_depth_qp*tanh(peak_depth_qp));
peak_output_nu = sqrt(2*peak_depth_qp*tanh(2*peak_depth_qp));
lambda_minus = 2*peak_nu-peak_output_nu;
lambda_plus = 2*peak_nu+peak_output_nu;
if lambda_minus<=0 || ~isfinite(lambda_minus) ...
        || lambda_plus<=0 || ~isfinite(lambda_plus)
    error('Frozen dual-branch Green--Laplace scales are invalid.');
end

safe_a = max(sqrt_A,realmin);
psi_resolvent_spectrum = complex(zeros(size(q)));
psi_slow_spectrum = complex(zeros(size(q)));
psi_fast_spectrum = complex(zeros(size(q)));
branch_names = ["slow","fast"];
branch_scales = [lambda_minus,lambda_plus];
branch_signs = [-1,1];
for branch_index = 1:2
    branch_name = branch_names(branch_index);
    branch_scale = branch_scales(branch_index);
    source_k_sign = branch_signs(branch_index);
    for node_index = 1:2
        node = nodes(node_index);
        time = node/branch_scale;
        coefficient = weights(node_index)*exp(node)/branch_scale/8;
        response = radial_balanced_branch_response( ...
            eta11_spectrum,qx,qy,q,nu,safe_nu,sqrt_A,safe_a, ...
            time,branch_name,source_k_sign,channel_masks, ...
            channel_shifts,channel_pair_supports);
        branch_contribution = 1i*coefficient*response;
        psi_resolvent_spectrum = psi_resolvent_spectrum ...
            +branch_contribution;
        if branch_name == "slow"
            psi_slow_spectrum = psi_slow_spectrum+branch_contribution;
        else
            psi_fast_spectrum = psi_fast_spectrum+branch_contribution;
        end
    end
end

% Denominator-free Taylor contact, evaluated directly from eta11.
eta11 = ifft2(eta11_spectrum);
phi1z = -1i*ifft2(eta11_spectrum.*nu);
taylor_contact = 0.5*eta11.*phi1z;
taylor_spectrum = fft2(taylor_contact);

% The selected engineering result uses no angular correction. Historical
% correction ablations are preserved only by the pre-cleanup archive tag.
psi_angular_correction_spectrum = complex(zeros(size(q)));
correction_transform_adjustment = -4;

psi_spectrum = psi_resolvent_spectrum+taylor_spectrum ...
    +psi_angular_correction_spectrum;
psi_spectrum(~output_support) = 0;
Psi22_plus = ifft2(psi_spectrum);
if any(~isfinite(Psi22_plus),'all')
    error('Dual-branch Green--Laplace Psi22 is not finite.');
end

fft_ifft_count = 10+33*channel_count+9*channel_count^2 ...
    +correction_transform_adjustment;
components = struct( ...
    'Psi22_spectrum',psi_spectrum, ...
    'Psi22_resolvent_spectrum',psi_resolvent_spectrum, ...
    'Psi22_slow_resolvent_spectrum',psi_slow_spectrum, ...
    'Psi22_fast_resolvent_spectrum',psi_fast_spectrum, ...
    'Psi22_taylor_contact_spectrum',taylor_spectrum, ...
    'Psi22_angular_correction_spectrum',psi_angular_correction_spectrum, ...
    'taylor_contact_plus',taylor_contact);
audit = struct( ...
    'candidate_id',char(freeze.candidate_id), ...
    'status','retained-engineering-selection', ...
    'output_variable','Psi22=phi(x,eta,t)|order2', ...
    'candidate_input_fields',{{'eta11'}}, ...
    'information_boundary','eta11-only', ...
    'explicit_eta22_dependency',false, ...
    'explicit_Phi22_dependency',false, ...
    'oracle_or_mf12_used',false, ...
    'model_parent_domain','q=kd>=0.3', ...
    'below_model_parent_domain',below_model_domain, ...
    'below_model_spectral_mass',below_model_spectral_mass, ...
    'diagnostic_tail_extrapolation_authorized', ...
        allow_below_model_parent_domain, ...
    'full_crossing_angle_dependence_retained',true, ...
    'exponential_balance_mode',char(exponential_balance_mode), ...
    'angular_correction_mode',char(correction_mode), ...
    'historical_angular_ablations_active',false, ...
    'radial_balance_channel_count',channel_count, ...
    'quadrature_rank_total',4, ...
    'quadrature_rank_per_branch',2, ...
    'lambda_minus',lambda_minus, ...
    'lambda_plus',lambda_plus, ...
    'fft_ifft_count',fft_ifft_count, ...
    'pair_loops',0, ...
    'cost_class','fixed-constant N log N');
end

function response = radial_balanced_branch_response( ...
        spectrum,qx,qy,q,nu,safe_nu,sqrt_A,safe_a,time,branch_name, ...
        source_k_sign,masks,shifts,pair_supports)
count = numel(masks);
fields = cell(1,count);
for channel = 1:count
    channel_spectrum = spectrum.*masks{channel} ...
        .*exp(-time*(nu-shifts(channel)));
    fields{channel} = source_fields( ...
        channel_spectrum,qx,qy,q,nu,safe_nu);
end
response = complex(zeros(size(q)));
for first = 1:count
    for second = 1:count
        [source_d_spectrum,source_k_spectrum] = ...
            cross_source_spectra(fields{first},fields{second});
        pair_support = pair_supports{first,second};
        source_d_spectrum(~pair_support) = 0;
        source_k_spectrum(~pair_support) = 0;
        shift_sum = shifts(first)+shifts(second);
        if branch_name == "slow"
            output_factor = exp((sqrt_A-shift_sum)*time);
        else
            output_factor = exp(-(sqrt_A+shift_sum)*time);
        end
        pair_response = output_factor.*(source_d_spectrum ...
            +source_k_sign*source_k_spectrum./safe_a);
        pair_response(~pair_support) = 0;
        response = response+pair_response;
    end
end
end

function fields = source_fields(spectrum,qx,qy,q,nu,safe_nu)
fields = struct( ...
    'v',ifft2(spectrum), ...
    'v_nu',ifft2(spectrum.*nu), ...
    'v_nu2',ifft2(spectrum.*nu.^2), ...
    'hx',ifft2(spectrum.*qx./safe_nu), ...
    'hy',ifft2(spectrum.*qy./safe_nu), ...
    'radial_over_nu',ifft2(spectrum.*q.^2./safe_nu), ...
    'jx',ifft2(spectrum.*qx), ...
    'jy',ifft2(spectrum.*qy));
end

function [source_d_spectrum,source_k_spectrum] = ...
        cross_source_spectra(left,right)
source_d = 2*left.v.*right.v_nu2+left.v_nu.*right.v_nu ...
    -left.hx.*right.hx-left.hy.*right.hy;
source_k = 2*left.v.*right.radial_over_nu ...
    +2*(left.hx.*right.jx+left.hy.*right.jy);
source_d_spectrum = fft2(source_d);
source_k_spectrum = fft2(source_k);
end

function [masks,shifts,pair_supports] = ...
        radial_balance_channels(support,nu,count)
active_nu = nu(support);
edges = linspace(min(active_nu),max(active_nu),count+1);
masks = cell(1,count);
shifts = zeros(1,count);
retained = false(1,count);
for channel = 1:count
    if channel<count
        masks{channel} = support & nu>=edges(channel) ...
            &nu<edges(channel+1);
    else
        masks{channel} = support & nu>=edges(channel) ...
            &nu<=edges(channel+1);
    end
    if ~any(masks{channel},'all'), continue, end
    retained(channel) = true;
    values = nu(masks{channel});
    shifts(channel) = 0.5*(min(values)+max(values));
end
masks = masks(retained);
shifts = shifts(retained);
count = numel(masks);
pair_supports = cell(count,count);
mask_transforms = cell(1,count);
for channel = 1:count
    mask_transforms{channel} = fft2(double(masks{channel}));
end
for first = 1:count
    for second = 1:count
        pair_supports{first,second} = real(ifft2( ...
            mask_transforms{first}.*mask_transforms{second}))>0.5;
    end
end
end

function freeze = load_and_verify_freeze(project_root)
path = fullfile(project_root,'symbolic','generated', ...
    'finite_depth_directional_order2_psi22_dual_branch_gl.json');
freeze = jsondecode(fileread(path));
if ~freeze.overall_exact_gate_pass || freeze.oracle_or_mf12_used ...
        || freeze.explicit_eta22_dependency ...
        || freeze.explicit_Phi22_dependency ...
        || freeze.time_evolution_used ...
        || ~isequal(string(freeze.candidate_input_fields),"eta11") ...
        || ~isequal(string(freeze.candidate_output_fields),"Psi22_plus") ...
        || freeze.quadrature.rank_total~=4 ...
        || freeze.quadrature.rank_per_branch~=2
    error('Frozen dual-branch Psi22 Green--Laplace interface is invalid.');
end
end


function [nodes,weights] = gauss_laguerre_two()
nodes = [2-sqrt(2),2+sqrt(2)];
weights = [nodes(1)/(4*(sqrt(2)-1)^2), ...
    nodes(2)/(4*(sqrt(2)+1)^2)];
end

function assert_quadratic_support_is_alias_safe(support)
[ny,nx] = size(support);
mode_x = [0:(ceil(nx/2)-1),-floor(nx/2):-1];
mode_y = [0:(ceil(ny/2)-1),-floor(ny/2):-1];
[column,row] = meshgrid(mode_x,mode_y);
if max(2*abs(column(support)))>=nx/2 ...
        || max(2*abs(row(support)))>=ny/2
    error('Quadratic output support reaches a Nyquist boundary.');
end
end
