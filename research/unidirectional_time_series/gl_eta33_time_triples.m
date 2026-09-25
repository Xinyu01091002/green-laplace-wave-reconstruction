function [eta,audit]=gl_eta33_time_triples(A,omega,k,h,kp,t,J,bins)
% Ordered triple transfer of src/internal/gl_no_stokes_eta33.m, elevation only.
% Uniform temporal bins enable exact aggregation of equal sum frequencies.
if nargin<8,bins=[];end
A=A(:);omega=omega(:);k=k(:);t=t(:); n=numel(A);
assert(all(h*k>0.5) && ismember(J,[4,6,8]) && all(isfinite([A;omega;k;t])));
q=h*k; nu=sqrt(q.*tanh(q)); qp=h*kp;
lambda2=2*sqrt(qp*tanh(qp))-sqrt(2*qp*tanh(2*qp));
lambda3=3*sqrt(qp*tanh(qp))-sqrt(3*qp*tanh(3*qp));
Q=q+q.'; s=nu+nu.'; G=Q.*tanh(Q);
sd=nu.^2+nu*nu.'+(nu.').^2-(q*q.')./(nu*nu.');
sk=(q.^2+q*q.')./nu+((q.').^2+q*q.')./nu.';
[S,C]=gl_laplace_terms(s,sqrt(G),lambda2,J);
e2=(-G.*S.*sd+C.*sk)/4;
p2=1i*(C.*sd-S.*sk)/4;
eta=complex(zeros(size(t))); temporal=eta;
if ~isempty(bins)
    bins=bins(:); N=numel(t); dt=mean(diff(t));
    assert(max(abs(t-(0:N-1)'*dt))<1e-8 && max(abs(omega-2*pi*bins/(N*dt)))<1e-10);
    assert(3*max(bins)<N/2);
    pairbins=bins+bins.';
end
for a=1:n
    % One single parent a times the ordered lower pair (b,c).
    fkpair=1i*(q(a)*Q+q(a)^2)/nu(a).*e2-(q(a)*Q+Q.^2).*p2;
    fdpair=-1i*s.*G.*p2-nu(a)^2*e2+1i*q(a)*Q/nu(a).*p2-1i*nu(a)*G.*p2;
    fkdirect=1i*(nu.*q)*q.'+0.5i*ones(n,1)*(q.^2.*nu).';
    fddirect=-0.5*ones(n,1)*(q.^2).'+(q./nu)*(nu.*q).'-nu*(q.^2./nu).';
    fk=fkpair/2+fkdirect/4; fd=fdpair/2+fddirect/4;
    outq=q(a)+Q; outG=outq.*tanh(outq);
    [So,Co]=gl_laplace_terms(nu(a)+s,sqrt(outG),lambda3,J);
    coeff=(outG.*So.*fd-1i*Co.*fk).*(A(a)*(A*A.'))/h^2;
    if isempty(bins)
        w=omega(a)+omega+omega.';
        eta=eta+exp(-1i*t*w(:).')*coeff(:);
    else
        index=1+bins(a)+pairbins;
        temporal=temporal+accumarray(index(:),coeff(:),[N,1]);
    end
end
if ~isempty(bins),eta=fft(temporal);end
assert(all(isfinite(eta)));
audit=struct('rank',J,'inner_scale',lambda2,'outer_scale',lambda3, ...
    'parents',n,'ordered_triples',n^3,'no_stokes_correction',true, ...
    'source','src/internal/gl_no_stokes_eta33.m','mean_or_difference_included',false);
end
