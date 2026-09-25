function psi=mf12_psi22_time(A,k,g,h,t)
% Same surface-potential coefficients and convention as gl_mf12_order2_comparison.
A=A(:).';k=k(:).';t=t(:);
c=mf12_spectral_coefficients(2,g,h,real(A),imag(A),k,zeros(size(k)),0,0, ...
    0,struct('enable_subharmonic',false));
coeff=1i*complex(c.A_2,c.B_2).*c.mu_2;
psi=exp(-1i*t*(2*c.omega))*coeff.';
plus=1:2:numel(c.G_npm);
coeff=1i*complex(c.A_npm(plus),c.B_npm(plus)).*c.mu_npm(plus);
psi=real(psi+exp(-1i*t*c.omega_npm(plus))*coeff.');assert(all(isfinite(psi)));
end
