function [psi33, psi33_hat, meta] = finite_depth_directional_psi33_two_scale_gl(eta11_hat,kx,ky,cfg)
%FINITE_DEPTH_DIRECTIONAL_PSI33_TWO_SCALE_GL Frozen eta11-only surface Psi33.

if nargin<4, cfg = struct(); end
cfg = default_field(cfg,'mode','two_scale');
cfg = default_field(cfg,'q_min',0.5000001);
cfg = default_field(cfg,'q_max',inf);
cfg = default_field(cfg,'forward_axis',[1,0]);
cfg = default_field(cfg,'cone_half_angle',pi/4);
cfg = default_field(cfg,'support_tolerance',1e-12);
cfg = default_field(cfg,'freeze_json','');
if ~ismember(cfg.mode,{'two_scale','inner_shared_outer_two_scale','shared_scale'})
    error('psi33_two_scale_gl:mode', ...
        'mode must be two_scale, inner_shared_outer_two_scale, or shared_scale.');
end
if ~isequal(size(eta11_hat),size(kx),size(ky))
    error('psi33_two_scale_gl:shape','eta11_hat, kx, and ky must have equal size.');
end
if cfg.cone_half_angle>=pi/2 || cfg.cone_half_angle<0 || cfg.q_min<=0 ...
        || cfg.q_max<=cfg.q_min || cfg.support_tolerance<0 ...
        || numel(cfg.forward_axis)~=2 || norm(cfg.forward_axis)==0
    error('psi33_two_scale_gl:config','Invalid support or forward-axis configuration.');
end

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
eta_freeze_path = fullfile(root,'symbolic','generated', ...
    'finite_depth_directional_order3_two_scale_gl.json');
if isempty(cfg.freeze_json)
    psi_freeze_path = fullfile(root,'symbolic','generated', ...
        'finite_depth_directional_order3_psi33_two_scale_gl.json');
else
    psi_freeze_path = cfg.freeze_json;
end
freeze = jsondecode(fileread(eta_freeze_path));
psi_freeze = jsondecode(fileread(psi_freeze_path));

q = hypot(kx,ky);
q2 = q.^2;
G = q.*tanh(q);
omega = sqrt(G);
invomega = zeros(size(q));
nonzero_q = q>0;
invomega(nonzero_q) = 1./omega(nonzero_q);
axis = cfg.forward_axis./norm(cfg.forward_axis);
projection = kx.*axis(1)+ky.*axis(2);
amplitude_scale = max(abs(eta11_hat),[],'all');
active = abs(eta11_hat)>cfg.support_tolerance*max(amplitude_scale,1);
active_nonzero = active & nonzero_q;
if any(active & ~nonzero_q,'all')
    error('psi33_two_scale_gl:zeroInput','eta11 input contains a zero mode.');
end
if any(active_nonzero & (q<cfg.q_min | q>cfg.q_max),'all')
    error('psi33_two_scale_gl:parentSupport','eta11 parent support is outside [q_min,q_max].');
end
if any(active_nonzero & projection<=q*cos(cfg.cone_half_angle),'all')
    error('psi33_two_scale_gl:forwardCone','eta11 parent is outside the strict forward cone.');
end
if any(active_nonzero,'all')
    if 3*max(abs(kx(active_nonzero)))>=max(abs(kx),[],'all') ...
            || 3*max(abs(ky(active_nonzero)))>=max(abs(ky),[],'all')
        error('psi33_two_scale_gl:aliasing','Cubic support reaches an FFT Nyquist boundary.');
    end
    balance_velocity = 0.5*min(omega(active_nonzero)./projection(active_nonzero));
else
    balance_velocity = 0;
end

support_spectrum = fft2(double(eta11_hat~=0));
pair_support = round(real(ifft2(support_spectrum.^2)))>0;
triple_support = round(real(ifft2(support_spectrum.^3)))>0;
mask2 = stage_mask(2) & pair_support;
mask3 = stage_mask(3) & triple_support;
counts = struct('fft',1,'ifft',2,'support_fft',1,'support_ifft',2, ...
    'products',0,'pair_solves_full',0, ...
    'pair_solves_value',0,'outer_nodes',0,'pair_nodes_full',0, ...
    'pair_nodes_value',0);
actual_times = struct('order2_slow',[],'order2_fast',[], ...
    'order2_shared',[],'order3_slow',[],'order3_fast',[], ...
    'order3_shared',[]);

phi33_hat = zeros(size(q));
outer_rules = stage_rules(3);
for ib = 1:numel(outer_rules)
    outer_rule = outer_rules{ib};
    for jj = 1:numel(outer_rule.nodes)
        outer_tau_node = outer_rule.nodes(jj)/outer_rule.scale;
        counts.outer_nodes = counts.outer_nodes+1;
        actual_times.(outer_rule.time_field)(end+1) = outer_tau_node;
        numerator_hat = order3_phi_numerator_hat(outer_tau_node);
        phi33_hat = phi33_hat+resolvent_filter(numerator_hat.*mask3,outer_tau_node, ...
            outer_rule.weights(jj),outer_rule.nodes(jj), ...
            outer_rule.scale,outer_rule.kind);
    end
end

pair0 = solve_pair_value(0);
b0 = first_surface_bank();
eta20 = ifft_counted(pair0.eta0);
phi2z0 = ifft_counted(G.*pair0.phi0);
surface_correction = add_values( ...
    product_value(b0.eta,phi2z0), ...
    product_value(eta20,b0.phiz), ...
    0.5*product_value(product_value(b0.eta,b0.eta),b0.phizz));
surface_correction_hat = fft_counted(surface_correction);
psi33_hat = 2*(phi33_hat+surface_correction_hat).*mask3;
psi33_hat(~nonzero_q) = 0;
psi33 = ifft_counted(psi33_hat);

meta = struct();
meta.schema_version = 1;
meta.candidate = 'finite_depth_directional_psi33_two_scale_gl';
meta.target = 'surface_Psi33';
meta.flat_phi33_is_internal_only = true;
meta.mode = cfg.mode;
meta.status = psi_freeze.status;
meta.eta11_only = true;
meta.external_lower_order_fields = false;
meta.execution_graph = 'exact_channel_fusion_v1';
meta.exponential_balance = struct('type','additive_forward_projection', ...
    'velocity',balance_velocity,'identity','sum(k_parallel)=K_parallel');
meta.runtime_kernel_dependency = false;
meta.q_parent = [cfg.q_min,cfg.q_max];
meta.output_support_order3 = [3*cfg.q_min*cos(cfg.cone_half_angle),3*cfg.q_max];
meta.zero_output_value = 0;
meta.gauge = 'C3(t)=0; strict K=0 coefficient is zero';
meta.counts = counts;
meta.actual_auxiliary_times = actual_times;
meta.freeze_json = psi_freeze_path;
meta.reused_resolvent_json = eta_freeze_path;
meta.claim_boundary = psi_freeze.claim_boundary;

    function mask = stage_mask(stage)
        mask = q>=stage*cfg.q_min*cos(cfg.cone_half_angle) ...
            & q<=stage*cfg.q_max & projection>0 & nonzero_q;
    end

    function rules = stage_rules(stage)
        prefix = ['order' num2str(stage)];
        use_two_scale = strcmp(cfg.mode,'two_scale') ...
            || (strcmp(cfg.mode,'inner_shared_outer_two_scale') && stage==3);
        if use_two_scale
            slow = gl_rule(freeze.primary_node_counts.(prefix).slow);
            fast = gl_rule(freeze.primary_node_counts.(prefix).fast);
            rules = {make_rule(slow,number(freeze.scales.([prefix '_slow'])), ...
                'slow',[prefix '_slow']), ...
                make_rule(fast,number(freeze.scales.([prefix '_fast'])), ...
                'fast',[prefix '_fast'])};
        else
            shared = gl_rule(freeze.shared_control_node_counts.(prefix));
            rules = {make_rule(shared,number(freeze.scales.([prefix '_shared'])), ...
                'shared',[prefix '_shared'])};
        end
    end

    function out = make_rule(gl,scale,kind,time_field)
        out = struct('nodes',gl.nodes,'weights',gl.weights,'scale',scale, ...
            'kind',kind,'time_field',time_field);
    end

    function out = gl_rule(n)
        records = freeze.gauss_laguerre.(['n' num2str(n)]);
        out.nodes = arrayfun(@(x) number(x.node),records);
        out.weights = arrayfun(@(x) number(x.weight),records);
    end

    function value = number(x)
        if isnumeric(x), value = double(x);
        else, value = str2double(regexprep(x,'`.*$',''));
        end
    end

    function output_hat = resolvent_filter(input_hat,tau,wgl,xgl,scale,kind)
        aw = wgl*exp(xgl)/scale;
        multiplier = zeros(size(q));
        switch kind
            case 'slow'
                multiplier(nonzero_q) = aw*exp(tau*(omega(nonzero_q) ...
                    -balance_velocity*projection(nonzero_q))) ...
                    .*invomega(nonzero_q)/2;
            case 'fast'
                multiplier(nonzero_q) = -aw*exp(-tau*(omega(nonzero_q) ...
                    +balance_velocity*projection(nonzero_q))) ...
                    .*invomega(nonzero_q)/2;
            case 'shared'
                multiplier(nonzero_q) = aw*exp(-tau*balance_velocity ...
                    *projection(nonzero_q)).*sinh(tau*omega(nonzero_q)) ...
                    .*invomega(nonzero_q);
        end
        output_hat = multiplier.*input_hat;
        output_hat(~nonzero_q) = 0;
    end

    function numerator_hat = order3_phi_numerator_hat(outer_tau)
        b1 = first_bank_order1(outer_tau);
        pair = solve_pair_full(outer_tau);
        e2 = spectral_pair_jet(pair.eta0,pair.eta1,1,1);
        e2x = spectral_pair_jet(pair.eta0,pair.eta1,1i*kx,1);
        e2y = spectral_pair_jet(pair.eta0,pair.eta1,1i*ky,1);
        p2x = spectral_pair_jet(pair.phi0,pair.phi1,1i*kx,1);
        p2y = spectral_pair_jet(pair.phi0,pair.phi1,1i*ky,1);
        p2z = spectral_pair_jet2(pair.phi0,pair.phi1,pair.phi2,G);
        p2zz = spectral_pair_jet(pair.phi0,pair.phi1,q2,1);
        p2tz = {1i*p2z{2},1i*p2z{3}};

        fk0 = add_values( ...
            product_value(b1.phix{1},e2x{1}),product_value(b1.phiy{1},e2y{1}), ...
            product_value(p2x{1},b1.etax{1}),product_value(p2y{1},b1.etay{1}), ...
            product_value(b1.eta{1},add_values( ...
                product_value(b1.phizx{1},b1.etax{1}), ...
                product_value(b1.phizy{1},b1.etay{1}))), ...
            -product_value(b1.eta{1},p2zz{1}), ...
            -product_value(e2{1},b1.phizz{1}), ...
            -0.5*product_value(product_value(b1.eta{1},b1.eta{1}),b1.phizzz{1}));

        fd = jet_add( ...
            jet_product(b1.eta,p2tz,1),jet_product(e2,b1.phitz,1), ...
            jet_scale(jet_product(jet_square(b1.eta,1),b1.phitzz,1),0.5), ...
            jet_product(b1.phix,p2x,1),jet_product(b1.phiy,p2y,1), ...
            jet_product(b1.phiz,{p2z{1},p2z{2}},1), ...
            jet_product(b1.eta,jet_add( ...
                jet_product(b1.phix,b1.phizx,1), ...
                jet_product(b1.phiy,b1.phizy,1), ...
                jet_product(b1.phiz,b1.phizz,1)),1));
        fd1_hat = fft_counted(fd{2});
        fk0_hat = fft_counted(fk0);
        numerator_hat = 1i*fd1_hat-fk0_hat;
    end

    function pair = solve_pair_full(outer_tau)
        counts.pair_solves_full = counts.pair_solves_full+1;
        pair = zero_pair(true);
        rules = stage_rules(2);
        for ir = 1:numel(rules)
            rule = rules{ir};
            for ii = 1:numel(rule.nodes)
                inner_tau = rule.nodes(ii)/rule.scale;
                counts.pair_nodes_full = counts.pair_nodes_full+1;
                actual_times.(rule.time_field)(end+1) = inner_tau;
                b = first_bank_inner_order3(outer_tau+inner_tau);
                fk = order2_fk(b,2);
                fd = order2_fd(b,3);
                fkh = cellfun_explicit(fk);
                fdh = cellfun_explicit(fd);
                pair.eta0 = pair.eta0+apply(fdh{1}.*G+1i*fkh{2});
                pair.eta1 = pair.eta1+apply(fdh{2}.*G+1i*fkh{3});
                pair.phi0 = pair.phi0+apply(1i*fdh{2}-fkh{1});
                pair.phi1 = pair.phi1+apply(1i*fdh{3}-fkh{2});
                pair.phi2 = pair.phi2+apply(1i*fdh{4}-fkh{3});
            end
        end
        pair = mask_pair(pair,true);

        function out = apply(in)
            out = resolvent_filter(in.*mask2,inner_tau,rule.weights(ii),rule.nodes(ii), ...
                rule.scale,rule.kind);
        end
    end

    function pair = solve_pair_value(outer_tau)
        counts.pair_solves_value = counts.pair_solves_value+1;
        pair = zero_pair(false);
        rules = stage_rules(2);
        for ir = 1:numel(rules)
            rule = rules{ir};
            for ii = 1:numel(rule.nodes)
                inner_tau = rule.nodes(ii)/rule.scale;
                counts.pair_nodes_value = counts.pair_nodes_value+1;
                b = first_bank_inner_order1(outer_tau+inner_tau);
                fk = order2_fk(b,1);
                fd = order2_fd(b,1);
                fd0_hat = fft_counted(fd{1});
                fk1_hat = fft_counted(fk{2});
                phi0_numerator_hat = fft_counted(1i*fd{2}-fk{1});
                pair.eta0 = pair.eta0+apply(fd0_hat.*G+1i*fk1_hat);
                pair.phi0 = pair.phi0+apply(phi0_numerator_hat);
            end
        end
        pair = mask_pair(pair,false);

        function out = apply(in)
            out = resolvent_filter(in.*mask2,inner_tau,rule.weights(ii),rule.nodes(ii), ...
                rule.scale,rule.kind);
        end
    end

    function pair = zero_pair(full)
        pair.eta0 = zeros(size(q)); pair.phi0 = zeros(size(q));
        if full
            pair.eta1 = zeros(size(q)); pair.phi1 = zeros(size(q));
            pair.phi2 = zeros(size(q));
        end
    end

    function pair = mask_pair(pair,full)
        pair.eta0 = pair.eta0.*mask2; pair.phi0 = pair.phi0.*mask2;
        if full
            pair.eta1 = pair.eta1.*mask2; pair.phi1 = pair.phi1.*mask2;
            pair.phi2 = pair.phi2.*mask2;
        end
    end

    function fk = order2_fk(b,order)
        fk = jet_add(jet_product(b.phix,b.etax,order), ...
            jet_product(b.phiy,b.etay,order), ...
            jet_scale(jet_product(b.eta,b.phizz,order),-1));
    end

    function fd = order2_fd(b,order)
        fd = jet_add(jet_product(b.eta,b.phitz,order), ...
            jet_scale(jet_add(jet_square(b.phix,order), ...
            jet_square(b.phiy,order),jet_square(b.phiz,order)),0.5));
    end

    function hats = cellfun_explicit(values)
        hats = cell(size(values));
        for kk = 1:numel(values), hats{kk} = fft_counted(values{kk}); end
    end

    function b = first_bank_order1(tau)
        d = balanced_input(tau);
        a0=I(d); aw1=I(omega.*d); aw2=I(omega.^2.*d); aw3=I(omega.^3.*d);
        xm1=I(kx.*invomega.*d); x0=I(kx.*d); x1=I(kx.*omega.*d); x2=I(kx.*omega.^2.*d);
        ym1=I(ky.*invomega.*d); y0=I(ky.*d); y1=I(ky.*omega.*d); y2=I(ky.*omega.^2.*d);
        qm1=I(q2.*invomega.*d); q0=I(q2.*d); q1=I(q2.*omega.*d); q2w=I(q2.*omega.^2.*d);
        b.eta={a0,-aw1}; b.etax={1i*x0,-1i*x1}; b.etay={1i*y0,-1i*y1};
        b.phix={xm1,-x0}; b.phiy={ym1,-y0}; b.phiz={-1i*aw1,1i*aw2};
        b.phizz={-1i*qm1,1i*q0}; b.phizzz={-1i*q1,1i*q2w};
        b.phizx={x1,-x2}; b.phizy={y1,-y2}; b.phitz={-aw2,aw3};
        b.phitzz={-q0,q1};
    end

    function b = first_bank_inner_order3(tau)
        d = balanced_input(tau);
        a=cell(1,6); for r=0:5, a{r+1}=I(omega.^r.*d); end
        x=cell(1,5); y=cell(1,5);
        x{1}=I(kx.*invomega.*d); y{1}=I(ky.*invomega.*d);
        for r=0:3
            x{r+2}=I(kx.*omega.^r.*d); y{r+2}=I(ky.*omega.^r.*d);
        end
        z=cell(1,4); z{1}=I(q2.*invomega.*d);
        for r=0:2, z{r+2}=I(q2.*omega.^r.*d); end
        b.eta={a{1},-a{2},a{3},-a{4}};
        b.etax={1i*x{2},-1i*x{3},1i*x{4},-1i*x{5}};
        b.etay={1i*y{2},-1i*y{3},1i*y{4},-1i*y{5}};
        b.phix={x{1},-x{2},x{3},-x{4}}; b.phiy={y{1},-y{2},y{3},-y{4}};
        b.phiz={-1i*a{2},1i*a{3},-1i*a{4},1i*a{5}};
        b.phizz={-1i*z{1},1i*z{2},-1i*z{3},1i*z{4}};
        b.phitz={-a{3},a{4},-a{5},a{6}};
    end

    function b = first_bank_inner_order1(tau)
        d = balanced_input(tau);
        a0=I(d); aw1=I(omega.*d); aw2=I(omega.^2.*d); aw3=I(omega.^3.*d);
        xm1=I(kx.*invomega.*d); x0=I(kx.*d); x1=I(kx.*omega.*d);
        ym1=I(ky.*invomega.*d); y0=I(ky.*d); y1=I(ky.*omega.*d);
        qm1=I(q2.*invomega.*d); q0=I(q2.*d);
        b.eta={a0,-aw1}; b.etax={1i*x0,-1i*x1}; b.etay={1i*y0,-1i*y1};
        b.phix={xm1,-x0}; b.phiy={ym1,-y0}; b.phiz={-1i*aw1,1i*aw2};
        b.phizz={-1i*qm1,1i*q0}; b.phitz={-aw2,aw3};
    end

    function b = first_surface_bank()
        d = eta11_hat/2;
        b.eta = I(d); b.phiz = -1i*I(omega.*d);
        b.phizz = -1i*I(q2.*invomega.*d);
    end

    function d = balanced_input(tau)
        d = (eta11_hat/2).*exp(-tau*(omega-balance_velocity*projection));
    end

    function out = I(in), out = ifft_counted(in); end

    function jet = spectral_pair_jet(v0,v1,multiplier,order)
        jet={ifft_counted(multiplier.*v0)};
        if order>=1, jet{2}=ifft_counted(multiplier.*v1); end
    end

    function jet = spectral_pair_jet2(v0,v1,v2,multiplier)
        jet={ifft_counted(multiplier.*v0),ifft_counted(multiplier.*v1), ...
            ifft_counted(multiplier.*v2)};
    end

    function out = jet_product(a,b,order)
        out=cell(1,order+1);
        for r=0:order
            value=zeros(size(q));
            for s=0:r
                counts.products=counts.products+1;
                value=value+nchoosek(r,s)*(a{s+1}.*b{r-s+1});
            end
            out{r+1}=value;
        end
    end

    function out = jet_square(a,order)
        out=cell(1,order+1);
        for r=0:order
            value=zeros(size(q));
            for s=0:floor(r/2)
                other=r-s;
                coefficient=nchoosek(r,s);
                if s~=other, coefficient=2*coefficient; end
                counts.products=counts.products+1;
                value=value+coefficient*(a{s+1}.*a{other+1});
            end
            out{r+1}=value;
        end
    end

    function out = jet_scale(a,factor)
        out=cellfun(@(x) factor*x,a,'UniformOutput',false);
    end

    function out = jet_add(varargin)
        out=cell(size(varargin{1}));
        for r=1:numel(out)
            value=zeros(size(q));
            for aa=1:nargin, value=value+varargin{aa}{r}; end
            out{r}=value;
        end
    end

    function out = product_value(a,b)
        counts.products=counts.products+1; out=a.*b;
    end

    function out = add_values(varargin)
        out=zeros(size(q)); for aa=1:nargin, out=out+varargin{aa}; end
    end

    function out = fft_counted(in), counts.fft=counts.fft+1; out=fft2(in); end
    function out = ifft_counted(in), counts.ifft=counts.ifft+1; out=ifft2(in); end
end

function out = default_field(in,name,value)
out=in; if ~isfield(out,name), out.(name)=value; end
end
