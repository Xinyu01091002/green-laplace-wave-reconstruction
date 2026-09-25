function eta=mf12_eta20_time(A,k,g,h,t)
% Difference coefficients only. Never computes third-order MF12.
A=A(:).';k=k(:).';t=t(:);
c=mf12_spectral_coefficients(2,g,h,real(A),imag(A),k,zeros(size(k)),0,0, ...
    0,struct('enable_subharmonic',true));
minus=2:2:numel(c.G_npm);
C=complex(c.A_npm(minus),c.B_npm(minus)).*c.G_npm(minus);
eta=real(exp(-1i*t*c.omega_npm(minus))*C.');
assert(all(isfinite(eta)));
end
