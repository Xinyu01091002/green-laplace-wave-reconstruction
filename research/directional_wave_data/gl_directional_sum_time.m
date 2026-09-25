function [eta,psi]=gl_directional_sum_time(A,omega,kx,ky,g,h,kp,t,J)
% Directional eta22 GL-J and true surface psi22 GL2+2; ordered pairs.
A=A(:);omega=omega(:);kx=kx(:);ky=ky(:);t=t(:);
q=h*hypot(kx,ky);nu=sqrt(q.*tanh(q));assert(all(q>=.3) && all(kx>0));
dot12=h^2*(kx*kx.'+ky*ky.');Q=h*hypot(kx+kx.',ky+ky.');a=sqrt(Q.*tanh(Q));s=nu+nu.';
sd=nu.^2+nu*nu.'+(nu.').^2-dot12./(nu*nu.');
sk=(q.^2+dot12)./nu+((q.').^2+dot12)./nu.';
qp=h*kp;vp=sqrt(qp*tanh(qp));ap=sqrt(2*qp*tanh(2*qp));
[S,C]=gl_laplace_terms(s,a,sqrt(4*vp^2-ap^2),J);
etaCoeff=(-a.^2.*S.*sd+C.*sk)/4.*(A*A.')/h;
[V,D]=eig([1,1;1,3],'vector');[x,ix]=sort(D);weights=V(1,ix).^2;
phi=complex(zeros(size(Q)));scales=[2*vp-ap,2*vp+ap];
for branch=1:2
    signBranch=2*branch-3;
    for j=1:2
        factor=weights(j)/scales(branch)*exp(x(j)-(s+signBranch*a)*x(j)/scales(branch));
        phi=phi+1i/8*factor.*(sd+signBranch*sk./a);
    end
end
psiCoeff=sqrt(g/h)*(phi-1i*s/4).*(A*A.');
w=omega+omega.';w=w(:);etaCoeff=etaCoeff(:);psiCoeff=psiCoeff(:);
eta=complex(zeros(size(t)));psi=eta;
for start=1:512:numel(w)
    ids=start:min(start+511,numel(w));E=exp(-1i*t*w(ids).');
    eta=eta+E*etaCoeff(ids);psi=psi+E*psiCoeff(ids);
end
assert(all(isfinite([eta;psi])));
end
