function coeffs = gl_spectral_coefficients(order,g,h,a,b,kx,ky,Ux,Uy,varargin)
%GL_SPECTRAL_COEFFICIENTS Prepare a Green--Laplace spectral reconstruction.
%   The first nine inputs intentionally match mf12_spectral_coefficients.
%   An optional final struct controls the prescribed Green--Laplace graph.

opts = default_options();
if ~isempty(varargin)
    if numel(varargin)~=1 || ~isstruct(varargin{1})
        error('green_laplace:Options', ...
            'The optional input must be one options struct.');
    end
    opts = merge_options(opts,varargin{1});
end

validateattributes(order,{'numeric'},{'scalar','integer','>=',1,'<=',3});
validateattributes(g,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(h,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(Ux,{'numeric'},{'scalar','real','finite'});
validateattributes(Uy,{'numeric'},{'scalar','real','finite'});
if Ux~=0 || Uy~=0
    error('green_laplace:UniformCurrent', ...
        ['The released API preserves the MF12-compatible Ux, Uy input slots ' ...
         'but has only been certified for Ux=Uy=0.']);
end

a = double(a(:).');
b = double(b(:).');
kx = double(kx(:).');
ky = double(ky(:).');
if isempty(a) || ~isequal(numel(a),numel(b),numel(kx),numel(ky))
    error('green_laplace:ComponentSize', ...
        'a, b, kx, and ky must be nonempty vectors of equal length.');
end
if any(~isfinite([a,b,kx,ky])) || any(hypot(kx,ky)<=0)
    error('green_laplace:Components', ...
        'Component amplitudes and nonzero wavenumbers must be finite.');
end
if any(kx<=0)
    error('green_laplace:Support', ...
        'The released GL graphs require strict-forward support kx>0.');
end
if opts.sector~="pure-sum"
    error('green_laplace:Sector', ...
        'The released total-field API supports the positive pure-sum sector only.');
end
validateattributes(opts.eta22_rank,{'numeric'}, ...
    {'scalar','integer','positive','<=',12});

amplitude = hypot(a,b);
if ~any(amplitude>0)
    error('green_laplace:ZeroInput','At least one input amplitude is required.');
end
if isempty(opts.peak_wavenumber)
    [~,peak_index] = max(amplitude.^2);
    peak_wavenumber = hypot(kx(peak_index),ky(peak_index));
else
    validateattributes(opts.peak_wavenumber,{'numeric'}, ...
        {'scalar','real','finite','positive'});
    peak_wavenumber = double(opts.peak_wavenumber);
end

parent_q = h*hypot(kx,ky);
if order>=2 && any(parent_q<0.3) && ~opts.allow_below_model_parent_domain
    error('green_laplace:Order2Domain', ...
        ['Order-two surface-potential reconstruction requires kh>=0.3. ' ...
         'A labelled tail extrapolation can be enabled in options.']);
end
if order>=3 && any(parent_q<=0.5)
    error('green_laplace:Order3Domain', ...
        'Order-three reconstruction requires every parent kh>0.5.');
end

intrinsic_omega = sqrt(g*hypot(kx,ky).*tanh(h*hypot(kx,ky)));
coeffs = struct( ...
    'api','green-laplace-spectral-v1', ...
    'order',double(order), ...
    'g',double(g),'h',double(h), ...
    'a',a,'b',b,'kx',kx,'ky',ky, ...
    'Ux',double(Ux),'Uy',double(Uy), ...
    'intrinsic_omega',intrinsic_omega, ...
    'omega',kx*Ux+ky*Uy+intrinsic_omega, ...
    'peak_wavenumber',peak_wavenumber, ...
    'peak_depth',h*peak_wavenumber, ...
    'options',opts, ...
    'information_boundary','first-order surface elevation only', ...
    'surface_potential_definition','psi=phi(x,z=eta,t)', ...
    'flat_potential_is_public_output',false);
end

function opts = default_options()
opts = struct( ...
    'sector',"pure-sum", ...
    'eta22_rank',6, ...
    'peak_wavenumber',[], ...
    'allow_below_model_parent_domain',false, ...
    'grid_tolerance',1e-10);
end

function target = merge_options(target,source)
names = fieldnames(source);
allowed = fieldnames(target);
for index = 1:numel(names)
    if ~ismember(names{index},allowed)
        error('green_laplace:UnknownOption', ...
            'Unknown Green--Laplace option: %s',names{index});
    end
    target.(names{index}) = source.(names{index});
end
target.sector = string(target.sector);
target.allow_below_model_parent_domain = logical( ...
    target.allow_below_model_parent_domain);
end
