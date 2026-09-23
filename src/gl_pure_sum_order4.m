function [eta44,psi44,audit] = gl_pure_sum_order4(s,qx,qy,qp,rank)
%GL_PURE_SUM_ORDER4 Experimental dimensionless positive pure-sum order-four fields.
% Input s is fft2(eta11_plus/h); qx=kx*h, qy=ky*h, qp=kp*h.
% Outputs are analytic eta44/h and psi44/(h*sqrt(g*h)). Take real parts for
% physical fields. Strict-forward, quartically alias-safe support is required.
% No Stokes correction or externally supplied higher-order field is accepted.
if nargin<5,rank=4;end
root=string(fileparts(fileparts(mfilename('fullpath'))));
[eta44,audit,~,state]=finite_depth_directional_green_laplace_eta44_three_kernels( ...
    s,qx,qy,qp,root,quadrature_rank=rank);
psi44=state.Psi44;
audit.surface_fft_ifft_count=state.fft_ifft_count;
audit.surface_pointwise_product_count=state.pointwise_product_count;
audit.release_status='experimental no-Stokes migration; not a total-field API extension';
end
