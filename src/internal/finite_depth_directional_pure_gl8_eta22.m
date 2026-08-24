function [eta22_plus,audit] = finite_depth_directional_pure_gl8_eta22( ...
    eta11_spectrum,qx,qy,peak_depth_qp,project_root)
%FINITE_DEPTH_DIRECTIONAL_PURE_GL8_ETA22 Retained pure-GL eta22 graph.
% The graph contains no Stokes-diagonal or empirical angular correction.
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
freeze_path=fullfile(project_root,'symbolic','generated', ...
    'finite_depth_directional_order2_eta22_pure_gl8.json');
freeze=jsondecode(fileread(freeze_path));
if ~freeze.overall_exact_gate_pass || freeze.quadrature_rank~=8 ...
        || freeze.oracle_or_mf12_used_in_formula ...
        || freeze.stokes_correction_used || freeze.angular_correction_used
    error('Invalid pure-GL8 Wolfram freeze.');
end
nodes=reshape(freeze.nodes,1,[]);
weights=reshape(freeze.weights,1,[]);
if numel(nodes)~=8 || numel(weights)~=8
    error('The pure-GL8 freeze must contain eight nodes and weights.');
end

q=hypot(qx,qy);
support=abs(eta11_spectrum)>1e-13*max(1,max(abs(eta11_spectrum(:))));
if ~any(support,'all') || any(qx(support)<=0) || any(q(support)<=0)
    error('Strict-forward nonzero analytic parent support is required.');
end
assert_quadratic_support_is_alias_safe(support);
support_count=real(ifft2(fft2(double(support)).^2));
output_support=support_count>0.5;
nu=sqrt(q.*tanh(q));
A=q.*tanh(q);
sqrt_A=sqrt(A);
lambda=sqrt((2*sqrt(peak_depth_qp*tanh(peak_depth_qp)))^2 ...
    -2*peak_depth_qp*tanh(2*peak_depth_qp));
if ~isfinite(lambda) || lambda<=0
    error('Invalid determinant-root Green--Laplace scale.');
end

candidate_spectrum=complex(zeros(size(q)));
safe_nu=max(nu,realmin);
minimum_parent_nu=min(nu(support));
maximum_pair_a=max(sqrt_A(output_support));
balance_shift=0.5*(minimum_parent_nu+maximum_pair_a/2);
for node_index=1:8
    x=nodes(node_index);
    t=x/lambda;
    damping=exp(-t*(nu-balance_shift));
    v=ifft2(eta11_spectrum.*damping);
    v_nu=ifft2(eta11_spectrum.*damping.*nu);
    v_nu2=ifft2(eta11_spectrum.*damping.*nu.^2);
    hx=ifft2(eta11_spectrum.*damping.*qx./safe_nu);
    hy=ifft2(eta11_spectrum.*damping.*qy./safe_nu);
    radial_over_nu=ifft2(eta11_spectrum.*damping.*q.^2./safe_nu);
    jx=ifft2(eta11_spectrum.*damping.*qx);
    jy=ifft2(eta11_spectrum.*damping.*qy);
    source_d=2*v.*v_nu2+v_nu.^2-hx.^2-hy.^2;
    source_k=2*v.*radial_over_nu+2*(hx.*jx+hy.*jy);
    coefficient=weights(node_index)*exp(x)/lambda/4 ...
        *exp(-2*balance_shift*t);
    sinh_over_a=t*ones(size(q));
    nonzero=sqrt_A>0;
    sinh_over_a(nonzero)=sinh(sqrt_A(nonzero)*t)./sqrt_A(nonzero);
    source_d_spectrum=fft2(source_d);
    source_k_spectrum=fft2(source_k);
    source_d_spectrum(~output_support)=0;
    source_k_spectrum(~output_support)=0;
    source_d_spectrum=suppress_fft_roundoff( ...
        source_d_spectrum,output_support);
    source_k_spectrum=suppress_fft_roundoff( ...
        source_k_spectrum,output_support);
    candidate_spectrum=candidate_spectrum+coefficient*( ...
        -A.*sinh_over_a.*source_d_spectrum ...
        +cosh(sqrt_A*t).*source_k_spectrum);
end
eta22_plus=ifft2(candidate_spectrum);
audit=struct( ...
    'candidate_id',freeze.candidate_id, ...
    'status','retained-bounded-pure-gl8', ...
    'information_boundary','eta11-only', ...
    'quadrature_rank',8, ...
    'quadrature_exact_moment_degree',15, ...
    'stokes_correction_used',false, ...
    'angular_correction_used',false, ...
    'learned_coefficients',0, ...
    'peak_depth_qp',peak_depth_qp, ...
    'quadrature_scale',lambda, ...
    'balance_shift',balance_shift, ...
    'stable_fft_backend',true, ...
    'fft_ifft_count',85, ...
    'pointwise_product_count',57, ...
    'pair_loops',0, ...
    'wolfram_freeze',freeze_path);
end

function assert_quadratic_support_is_alias_safe(support)
[ny,nx]=size(support);
mode_x=[0:(ceil(nx/2)-1),-floor(nx/2):-1];
mode_y=[0:(ceil(ny/2)-1),-floor(ny/2):-1];
[column,row]=meshgrid(mode_x,mode_y);
if max(2*abs(column(support)))>=nx/2 || max(2*abs(row(support)))>=ny/2
    error('Quadratic output support reaches a Nyquist boundary.');
end
end

function spectrum=suppress_fft_roundoff(spectrum,support)
active_values=abs(spectrum(support));
if isempty(active_values),return,end
scale=max(active_values);
if scale==0,return,end
relative_floor=16*eps*max(1,log2(numel(spectrum)));
spectrum(support & abs(spectrum)<relative_floor*scale)=0;
end
