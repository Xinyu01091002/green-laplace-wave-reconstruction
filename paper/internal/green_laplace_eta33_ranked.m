function varargout = green_laplace_eta33_ranked(s,qx,qy,qp,rank,root,outer_rank)
% Paper entry point using the same no-Stokes implementation as the public API.
if nargin<7,outer_rank=rank;end
[varargout{1:nargout}] = gl_no_stokes_eta33(s,qx,qy,qp,rank,root,outer_rank);
end
