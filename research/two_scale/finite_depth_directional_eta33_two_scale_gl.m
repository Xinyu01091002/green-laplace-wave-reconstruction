function [eta33, eta33_hat, meta] = finite_depth_directional_eta33_two_scale_gl(eta11_hat, kx, ky, cfg)
%FINITE_DEPTH_DIRECTIONAL_ETA33_TWO_SCALE_GL Frozen eta11-only eta33 graph.
%
% eta11_hat is the positive-frequency analytic eta11 spectrum.  No external
% eta22, Phi22, Psi22, eta33, MF12, or oracle field is accepted.  The
% internally generated eta22/Phi22 state uses the same frozen Green--Laplace
% solve and is included in the operation count.

if nargin<4
    cfg = struct();
end
cfg = default_field(cfg,'mode','two_scale');
cfg = default_field(cfg,'q_min',0.5000001);
cfg = default_field(cfg,'q_max',inf);
cfg = default_field(cfg,'forward_axis',[1,0]);
cfg = default_field(cfg,'cone_half_angle',pi/4);
cfg = default_field(cfg,'support_tolerance',1e-12);
cfg = default_field(cfg,'freeze_json','');
if ~ismember(cfg.mode,{'two_scale','shared_scale'})
    error('eta33_two_scale_gl:mode','mode must be two_scale or shared_scale.');
end

if ~isequal(size(eta11_hat),size(kx),size(ky))
    error('eta33_two_scale_gl:shape','eta11_hat, kx, and ky must have equal size.');
end
if cfg.cone_half_angle >= pi/2
    error('eta33_two_scale_gl:cone','cone_half_angle must be strictly below pi/2.');
end
if cfg.q_min<=0 || cfg.q_max<=0 || cfg.cone_half_angle<0 ...
        || cfg.support_tolerance<0 || numel(cfg.forward_axis)~=2 ...
        || norm(cfg.forward_axis)==0
    error('eta33_two_scale_gl:config','Invalid support or forward-axis configuration.');
end
if cfg.q_max <= cfg.q_min
    error('eta33_two_scale_gl:support','q_max must exceed q_min.');
end

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
if isempty(cfg.freeze_json)
    freeze_path = fullfile(root,'symbolic','generated', ...
        'finite_depth_directional_order3_two_scale_gl.json');
else
    freeze_path = cfg.freeze_json;
end
freeze = jsondecode(fileread(freeze_path));

q = hypot(kx,ky);
G = q.*tanh(q);
omega = sqrt(G);
axis = cfg.forward_axis./norm(cfg.forward_axis);
projection = kx.*axis(1)+ky.*axis(2);
amplitude_scale = max(abs(eta11_hat),[],'all');
active = abs(eta11_hat) > cfg.support_tolerance*max(amplitude_scale,1);
active_nonzero = active & q>0;
if any(active & q==0,'all')
    error('eta33_two_scale_gl:zeroInput','eta11 input contains a zero mode.');
end
if any(active_nonzero & (q < cfg.q_min | q > cfg.q_max),'all')
    error('eta33_two_scale_gl:parentSupport','eta11 parent support is outside [q_min,q_max].');
end
if any(active_nonzero & projection <= q*cos(cfg.cone_half_angle),'all')
    error('eta33_two_scale_gl:forwardCone','eta11 parent is outside the strict forward cone.');
end
if any(active_nonzero,'all')
    kx_limit = max(abs(kx),[],'all');
    ky_limit = max(abs(ky),[],'all');
    if 3*max(abs(kx(active_nonzero))) >= kx_limit ...
            || 3*max(abs(ky(active_nonzero))) >= ky_limit
        error('eta33_two_scale_gl:aliasing', ...
            'Cubic pure-sum support reaches an FFT Nyquist boundary.');
    end
end

mask2 = stage_mask(2);
mask3 = stage_mask(3);
nonzero_q = q>0;
invomega = zeros(size(q));
invomega(nonzero_q) = 1./omega(nonzero_q);

counts = struct('fft',0,'ifft',0,'products',0, ...
    'pair_solves',0,'outer_nodes',0,'pair_nodes',0);
actual_times = struct('order2_slow',[],'order2_fast',[], ...
    'order2_shared',[],'order3_slow',[],'order3_fast',[], ...
    'order3_shared',[]);

eta33_hat = zeros(size(eta11_hat));
outer_rules = stage_rules(3);
for ib = 1:numel(outer_rules)
    outer_rule = outer_rules{ib};
    for j = 1:numel(outer_rule.nodes)
        tau = outer_rule.nodes(j)/outer_rule.scale;
        counts.outer_nodes = counts.outer_nodes+1;
        actual_times.(outer_rule.time_field)(end+1) = tau;
        numerator_hat = order3_eta_numerator_hat(tau);
        eta33_hat = eta33_hat + resolvent_filter( ...
            numerator_hat,tau,outer_rule.weights(j),outer_rule.nodes(j), ...
            outer_rule.scale,outer_rule.kind);
    end
end
eta33_hat = 2*eta33_hat.*mask3;
eta33 = ifft_counted(eta33_hat);

meta = struct();
meta.schema_version = 1;
meta.candidate = 'finite_depth_directional_eta33_two_scale_gl';
meta.mode = cfg.mode;
meta.status = freeze.status;
meta.eta11_only = true;
meta.external_lower_order_fields = false;
meta.execution_graph = 'exact_channel_fusion_v1';
meta.runtime_kernel_dependency = false;
meta.graph_equivalence_reference = 'direct_ordered_eta33_two_scale_gl';
meta.scope = freeze.scope;
meta.q_parent = [cfg.q_min,cfg.q_max];
meta.cone_half_angle = cfg.cone_half_angle;
meta.output_support_order2 = [2*cfg.q_min*cos(cfg.cone_half_angle),2*cfg.q_max];
meta.output_support_order3 = [3*cfg.q_min*cos(cfg.cone_half_angle),3*cfg.q_max];
meta.zero_output_value = 0;
meta.cubic_support_unaliased = true;
meta.counts = counts;
meta.actual_auxiliary_times = actual_times;
meta.freeze_json = freeze_path;
meta.source_sha256 = freeze.source_sha256;
meta.claim_boundary = freeze.claim_boundary;

    function mask = stage_mask(stage)
        lower = stage*cfg.q_min*cos(cfg.cone_half_angle);
        upper = stage*cfg.q_max;
        mask = q>=lower & q<=upper & projection>0 & q>0;
    end

    function rules = stage_rules(stage)
        if stage==2
            prefix = 'order2';
        else
            prefix = 'order3';
        end
        if strcmp(cfg.mode,'two_scale')
            slow_n = freeze.primary_node_counts.(prefix).slow;
            fast_n = freeze.primary_node_counts.(prefix).fast;
            slow = gl_rule(slow_n);
            fast = gl_rule(fast_n);
            rules = {
                make_rule(slow,number(freeze.scales.([prefix '_slow'])),'slow',[prefix '_slow']), ...
                make_rule(fast,number(freeze.scales.([prefix '_fast'])),'fast',[prefix '_fast'])};
        else
            n = freeze.shared_control_node_counts.(prefix);
            shared = gl_rule(n);
            rules = {make_rule(shared,number(freeze.scales.([prefix '_shared'])), ...
                'shared',[prefix '_shared'])};
        end
    end

    function out = make_rule(gl,scale,kind,time_field)
        out = struct('nodes',gl.nodes,'weights',gl.weights, ...
            'scale',scale,'kind',kind,'time_field',time_field);
    end

    function out = gl_rule(n)
        records = freeze.gauss_laguerre.(['n' num2str(n)]);
        out.nodes = arrayfun(@(x) number(x.node),records);
        out.weights = arrayfun(@(x) number(x.weight),records);
    end

    function value = number(x)
        if isnumeric(x)
            value = double(x);
        else
            value = str2double(regexprep(x,'`.*$',''));
        end
    end

    function output_hat = resolvent_filter(input_hat,tau,wgl,xgl,scale,kind)
        actual_weight = wgl*exp(xgl)/scale;
        multiplier = zeros(size(q));
        switch kind
            case 'slow'
                multiplier(nonzero_q) = actual_weight*exp(tau*omega(nonzero_q)) ...
                    .*invomega(nonzero_q)/2;
            case 'fast'
                multiplier(nonzero_q) = -actual_weight*exp(-tau*omega(nonzero_q)) ...
                    .*invomega(nonzero_q)/2;
            case 'shared'
                multiplier(nonzero_q) = actual_weight*sinh(tau*omega(nonzero_q)) ...
                    .*invomega(nonzero_q);
        end
        output_hat = multiplier.*input_hat;
        output_hat(~nonzero_q) = 0;
    end

    function numerator_hat = order3_eta_numerator_hat(tau)
        b1 = first_bank_order1(tau);
        pair = solve_pair(tau);

        e2 = spectral_pair_jet(pair.eta0,pair.eta1,1,1);
        e2x = spectral_pair_jet(pair.eta0,pair.eta1,1i*kx,1);
        e2y = spectral_pair_jet(pair.eta0,pair.eta1,1i*ky,1);
        p2x = spectral_pair_jet(pair.phi0,pair.phi1,1i*kx,1);
        p2y = spectral_pair_jet(pair.phi0,pair.phi1,1i*ky,1);
        p2z = spectral_pair_jet(pair.phi0,pair.phi1,G,1);
        p2zz = spectral_pair_jet(pair.phi0,pair.phi1,q.^2,1);
        p2tz = 1i*p2z{2};

        fk = jet_add( ...
            jet_product(b1.phix,e2x,1), ...
            jet_product(b1.phiy,e2y,1), ...
            jet_product(p2x,b1.etax,1), ...
            jet_product(p2y,b1.etay,1), ...
            jet_product(b1.eta,jet_add( ...
                jet_product(b1.phizx,b1.etax,1), ...
                jet_product(b1.phizy,b1.etay,1)),1), ...
            jet_scale(jet_product(b1.eta,p2zz,1),-1), ...
            jet_scale(jet_product(e2,b1.phizz,1),-1), ...
            jet_scale(jet_product(jet_product(b1.eta,b1.eta,1),b1.phizzz,1),-0.5));

        fd = add_values( ...
            product_value(b1.eta{1},p2tz), ...
            product_value(e2{1},b1.phitz{1}), ...
            0.5*product_value(product_value(b1.eta{1},b1.eta{1}),b1.phitzz{1}), ...
            product_value(b1.phix{1},p2x{1}), ...
            product_value(b1.phiy{1},p2y{1}), ...
            product_value(b1.phiz{1},p2z{1}), ...
            product_value(b1.eta{1},add_values( ...
                product_value(b1.phix{1},b1.phizx{1}), ...
                product_value(b1.phiy{1},b1.phizy{1}), ...
                product_value(b1.phiz{1},b1.phizz{1}))));

        numerator_hat = G.*fft_counted(fd) + 1i*fft_counted(fk{2});
    end

    function pair = solve_pair(outer_tau)
        counts.pair_solves = counts.pair_solves+1;
        pair.eta0 = zeros(size(q));
        pair.eta1 = zeros(size(q));
        pair.phi0 = zeros(size(q));
        pair.phi1 = zeros(size(q));
        rules = stage_rules(2);
        for ir = 1:numel(rules)
            pair_rule = rules{ir};
            for jj = 1:numel(pair_rule.nodes)
                inner_tau = pair_rule.nodes(jj)/pair_rule.scale;
                counts.pair_nodes = counts.pair_nodes+1;
                actual_times.(pair_rule.time_field)(end+1) = inner_tau;
                b = first_bank_inner_order2(outer_tau+inner_tau);
                fk = jet_add( ...
                    jet_product(b.phix,b.etax,2), ...
                    jet_product(b.phiy,b.etay,2), ...
                    jet_scale(jet_product(b.eta,b.phizz,2),-1));
                fd = jet_add( ...
                    jet_product(b.eta,b.phitz,2), ...
                    jet_scale(jet_add( ...
                        jet_product(b.phix,b.phix,2), ...
                        jet_product(b.phiy,b.phiy,2), ...
                        jet_product(b.phiz,b.phiz,2)),0.5));

                fk_hat = {fft_counted(fk{1}),fft_counted(fk{2}), ...
                    fft_counted(fk{3})};
                fd_hat = {fft_counted(fd{1}),fft_counted(fd{2}), ...
                    fft_counted(fd{3})};
                neta0_hat = G.*fd_hat{1} + 1i*fk_hat{2};
                neta1_hat = G.*fd_hat{2} + 1i*fk_hat{3};
                nphi0_hat = 1i*fd_hat{2} - fk_hat{1};
                nphi1_hat = 1i*fd_hat{3} - fk_hat{2};

                pair.eta0 = pair.eta0 + resolvent_filter(neta0_hat, ...
                    inner_tau,pair_rule.weights(jj),pair_rule.nodes(jj),pair_rule.scale,pair_rule.kind);
                pair.eta1 = pair.eta1 + resolvent_filter(neta1_hat, ...
                    inner_tau,pair_rule.weights(jj),pair_rule.nodes(jj),pair_rule.scale,pair_rule.kind);
                pair.phi0 = pair.phi0 + resolvent_filter(nphi0_hat, ...
                    inner_tau,pair_rule.weights(jj),pair_rule.nodes(jj),pair_rule.scale,pair_rule.kind);
                pair.phi1 = pair.phi1 + resolvent_filter(nphi1_hat, ...
                    inner_tau,pair_rule.weights(jj),pair_rule.nodes(jj),pair_rule.scale,pair_rule.kind);
            end
        end
        pair.eta0 = pair.eta0.*mask2;
        pair.eta1 = pair.eta1.*mask2;
        pair.phi0 = pair.phi0.*mask2;
        pair.phi1 = pair.phi1.*mask2;
    end

    function b = first_bank_order1(tau)
        % Sixteen unique inverse transforms generate the original 24 jet
        % entries.  Scalar phase/sign relations are exact consequences of
        % the frozen linear eta11-to-Phi11 convention.
        damped = (eta11_hat/2).*exp(-tau*omega);
        a0 = ifft_counted(damped);
        aw1 = ifft_counted(omega.*damped);
        aw2 = ifft_counted(omega.^2.*damped);
        aw3 = ifft_counted(omega.^3.*damped);
        xm1 = ifft_counted(kx.*invomega.*damped);
        x0 = ifft_counted(kx.*damped);
        x1 = ifft_counted(kx.*omega.*damped);
        x2 = ifft_counted(kx.*omega.^2.*damped);
        ym1 = ifft_counted(ky.*invomega.*damped);
        y0 = ifft_counted(ky.*damped);
        y1 = ifft_counted(ky.*omega.*damped);
        y2 = ifft_counted(ky.*omega.^2.*damped);
        qm1 = ifft_counted(q.^2.*invomega.*damped);
        q0 = ifft_counted(q.^2.*damped);
        q1 = ifft_counted(q.^2.*omega.*damped);
        q2 = ifft_counted(q.^2.*omega.^2.*damped);

        b.eta = {a0,-aw1};
        b.etax = {1i*x0,-1i*x1};
        b.etay = {1i*y0,-1i*y1};
        b.phix = {xm1,-x0};
        b.phiy = {ym1,-y0};
        b.phiz = {-1i*aw1,1i*aw2};
        b.phizz = {-1i*qm1,1i*q0};
        b.phizzz = {-1i*q1,1i*q2};
        b.phizx = {x1,-x2};
        b.phizy = {y1,-y2};
        b.phitz = {-aw2,aw3};
        b.phitzz = {-q0,q1};
    end

    function b = first_bank_inner_order2(tau)
        % The order-two forcing uses only these eight fields.  Their 24 jet
        % entries collapse to sixteen unique spectral multipliers.
        damped = (eta11_hat/2).*exp(-tau*omega);
        a0 = ifft_counted(damped);
        aw1 = ifft_counted(omega.*damped);
        aw2 = ifft_counted(omega.^2.*damped);
        aw3 = ifft_counted(omega.^3.*damped);
        aw4 = ifft_counted(omega.^4.*damped);
        xm1 = ifft_counted(kx.*invomega.*damped);
        x0 = ifft_counted(kx.*damped);
        x1 = ifft_counted(kx.*omega.*damped);
        x2 = ifft_counted(kx.*omega.^2.*damped);
        ym1 = ifft_counted(ky.*invomega.*damped);
        y0 = ifft_counted(ky.*damped);
        y1 = ifft_counted(ky.*omega.*damped);
        y2 = ifft_counted(ky.*omega.^2.*damped);
        qm1 = ifft_counted(q.^2.*invomega.*damped);
        q0 = ifft_counted(q.^2.*damped);
        q1 = ifft_counted(q.^2.*omega.*damped);

        b.eta = {a0,-aw1,aw2};
        b.etax = {1i*x0,-1i*x1,1i*x2};
        b.etay = {1i*y0,-1i*y1,1i*y2};
        b.phix = {xm1,-x0,x1};
        b.phiy = {ym1,-y0,y1};
        b.phiz = {-1i*aw1,1i*aw2,-1i*aw3};
        b.phizz = {-1i*qm1,1i*q0,-1i*q1};
        b.phitz = {-aw2,aw3,-aw4};
    end

    function jet = spectral_pair_jet(value_hat,derivative_hat,multiplier,order)
        jet = cell(1,order+1);
        jet{1} = ifft_counted(multiplier.*value_hat);
        if order>=1
            jet{2} = ifft_counted(multiplier.*derivative_hat);
        end
    end

    function out = jet_product(a,b,order)
        out = cell(1,order+1);
        for rr = 0:order
            value = zeros(size(q));
            for ss = 0:rr
                counts.products = counts.products+1;
                value = value + nchoosek(rr,ss)*(a{ss+1}.*b{rr-ss+1});
            end
            out{rr+1} = value;
        end
    end

    function out = jet_scale(a,factor)
        out = cellfun(@(x) factor*x,a,'UniformOutput',false);
    end

    function out = jet_add(varargin)
        order = numel(varargin{1});
        out = cell(1,order);
        for rr = 1:order
            value = zeros(size(q));
            for aa = 1:nargin
                value = value+varargin{aa}{rr};
            end
            out{rr} = value;
        end
    end

    function out = product_value(a,b)
        counts.products = counts.products+1;
        out = a.*b;
    end

    function out = add_values(varargin)
        out = zeros(size(q));
        for aa = 1:nargin
            out = out+varargin{aa};
        end
    end

    function out = fft_counted(in)
        counts.fft = counts.fft+1;
        out = fft2(in);
    end

    function out = ifft_counted(in)
        counts.ifft = counts.ifft+1;
        out = ifft2(in);
    end
end

function out = default_field(in,name,value)
out = in;
if ~isfield(out,name)
    out.(name) = value;
end
end
