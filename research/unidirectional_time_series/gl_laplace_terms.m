function [S,C]=gl_laplace_terms(s,a,lambda,J)
% Stable prescribed Gauss-Laguerre sinh(a*t)/a and cosh(a*t) integrals.
assert(all(s>a,'all') && all(a>0,'all') && lambda>0);
[V,D]=eig(diag(2*(1:J)-1)+diag(1:J-1,1)+diag(1:J-1,-1),'vector');
[x,idx]=sort(D); w=V(1,idx).^2;
S=zeros(size(s)); C=S;
for j=1:J
    ep=exp(x(j)-(s-a)*x(j)/lambda);
    em=exp(x(j)-(s+a)*x(j)/lambda);
    S=S+w(j)*(ep-em)./(2*lambda*a);
    C=C+w(j)*(ep+em)/(2*lambda);
end
end
