function [psi,audit]=gl_psi22_time_pairs(A,omega,k,g,h,kp,t)
% Released dual-branch GL2+2 true surface potential, physical m^2/s.
A=A(:);omega=omega(:);k=k(:);t=t(:);
assert(all(k*h>=.3) && all(omega>0) && g>0 && h>0 && kp>0);
q=h*k;nu=sqrt(q.*tanh(q));Q=q+q.';a=sqrt(Q.*tanh(Q));s=nu+nu.';
sd=nu.^2+nu*nu.'+(nu.').^2-(q*q.')./(nu*nu.');
sk=(q.^2+q*q.')./nu+((q.').^2+q*q.')./nu.';
qp=h*kp;nup=sqrt(qp*tanh(qp));ap=sqrt(2*qp*tanh(2*qp));
scales=[2*nup-ap,2*nup+ap]; assert(all(scales>0));
[V,D]=eig([1,1;1,3],'vector');[x,ix]=sort(D);weights=V(1,ix).^2;
kernel=complex(zeros(size(Q)));
for branch=1:2
    signBranch=2*branch-3;
    for j=1:2
        factor=weights(j)/scales(branch)*exp(x(j)-(s+signBranch*a)*x(j)/scales(branch));
        kernel=kernel+1i/8*factor.*(sd+signBranch*sk./a);
    end
end
% 0.5 eta1 phi1_z, symmetrized over the two ordered parents.
taylor=-1i*s/4;
coeff=sqrt(g/h)*(kernel+taylor).*(A*A.'); w=omega+omega.';
psi=exp(-1i*t*w(:).')*coeff(:);assert(all(isfinite(psi)));
audit=struct('surface_potential',true,'quadrature_nodes_per_branch',2, ...
    'scales',scales,'stokes_or_angular_correction',false,'taylor_term_retained',true);
end
