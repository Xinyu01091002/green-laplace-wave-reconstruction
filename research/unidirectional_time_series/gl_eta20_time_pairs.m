function [eta,audit]=gl_eta20_time_pairs(A,omega,k,h,t,J)
% Frozen shared-scale eta20 diagnostic: nonzero differences only, physical m.
A=A(:);omega=omega(:);k=k(:);t=t(:);
assert(all(k>0) && ismember(J,[6,12,16]) && all(isfinite([A;omega;k;t])));
q=h*k; nu=sqrt(q.*tanh(q)); E=abs(A).^2;
q0=sum(q.*E)/sum(E); scale=sqrt(2*sum((q-q0).^2.*E)/sum(E));
assert(scale>0);
Q=q-q.'; sigma=nu-nu.'; outnu=sqrt(abs(Q).*tanh(abs(Q)));
fd=-nu.^2-(nu.').^2+nu*nu.'+(q./nu)*(q./nu).';
fk=1i*Q.*(q./nu+(q./nu).');
active=Q~=0; assert(all(outnu(active)>abs(sigma(active))));
[V,D]=eig(diag(2*(1:J)-1)+diag(1:J-1,1)+diag(1:J-1,-1),'vector');
[x,idx]=sort(D);w=V(1,idx).^2;
kernel=zeros(size(Q));
for j=1:J
    ep=exp(x(j)-(outnu-sigma)*x(j)/scale);
    em=exp(x(j)-(outnu+sigma)*x(j)/scale);
    increment=w(j)/(4*scale)*((1i*fk-outnu.*fd).*ep+(-1i*fk-outnu.*fd).*em);
    increment(~active)=0; kernel=kernel+increment;
end
% The published figure_eta20_gl_rank assembles 2*the raw diagnostic field.
% Keep that declared physical normalization; this is not an error-fitted gain.
coeff=kernel.*(A*A')/(2*h); % u=A/2 and physical assembly factor 2
w=omega-omega.';
eta=real(exp(-1i*t*w(:).')*coeff(:));
assert(all(isfinite(eta)));
audit=struct('rank',J,'scale',scale,'parents',numel(A),'strict_zero','excluded', ...
    'source','diagnostics/eta20/eta20_green_laplace_shared.m', ...
    'physical_assembly_factor',2, ...
    'status','shared-scale GL diagnostic; not production or a complete mean-flow model');
end
