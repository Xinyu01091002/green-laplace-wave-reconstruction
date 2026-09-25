function [eta,linear,source] = mf12_eta22_time(A,k,g,h,t)
% Same complex-amplitude convention as gl_mf12_order2_comparison.
A=A(:).'; k=k(:).'; t=t(:);
source=which('mf12_spectral_coefficients'); assert(~isempty(source));
c=mf12_spectral_coefficients(2,g,h,real(A),imag(A),k,zeros(size(k)),0,0, ...
    0,struct('enable_subharmonic',false));
linear=real(exp(-1i*t*c.omega)*A.');
C=complex(c.A_2,c.B_2).*c.G_2;
eta=exp(-1i*t*(2*c.omega))*C.';
plus=1:2:numel(c.G_npm);
C=complex(c.A_npm(plus),c.B_npm(plus)).*c.G_npm(plus);
eta=real(eta+exp(-1i*t*c.omega_npm(plus))*C.');
assert(all(isfinite(eta)));
end
