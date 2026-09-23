function varargout = finite_depth_directional_nested_gl_eta33_stable(s,qx,qy,qp,root)
% Compatibility entry point: pure GL4 with no Stokes correction.
[varargout{1:nargout}] = gl_no_stokes_eta33(s,qx,qy,qp,4,root,4);
end
