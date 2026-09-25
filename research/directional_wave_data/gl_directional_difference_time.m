function [eta,audit]=gl_directional_difference_time(A,omega,kx,ky,h,t,J,bins)
% Physical eta20 shared-scale GL diagnostic, all nonzero spatial differences.
if nargin<8,bins=[];end
A=A(:);omega=omega(:);kx=kx(:);ky=ky(:);t=t(:);
qx=h*kx;qy=h*ky;q=hypot(qx,qy);nu=sqrt(q.*tanh(q));E=abs(A).^2;E=E/sum(E);
scale=sqrt(2*sum(E.*((qx-sum(E.*qx)).^2+(qy-sum(E.*qy)).^2)));assert(scale>0);
Q=hypot(qx-qx.',qy-qy.');s=nu-nu.';a=sqrt(Q.*tanh(Q));dot12=qx*qx.'+qy*qy.';
fd=-nu.^2-(nu.').^2+nu*nu.'+dot12./(nu*nu.');
fk=1i*((q.^2-dot12)./nu+(dot12-(q.').^2)./nu.');
active=Q>0;assert(all(a(active)>abs(s(active))));
[V,D]=eig(diag(2*(1:J)-1)+diag(1:J-1,1)+diag(1:J-1,-1),'vector');[x,ix]=sort(D);w=V(1,ix).^2;
kernel=zeros(size(Q));
for j=1:J
    ep=exp(x(j)-(a-s)*x(j)/scale);em=exp(x(j)-(a+s)*x(j)/scale);
    v=w(j)/(4*scale)*((1i*fk-a.*fd).*ep+(-1i*fk-a.*fd).*em);v(~active)=0;kernel=kernel+v;
end
C=kernel.*(A*A')/(2*h);frequency=omega-omega.';
stationary=frequency==0 & active;
audit=struct('rank',J,'scale',scale,'nonzero_spatial_stationary_pairs',nnz(stationary), ...
    'stationary_nonzero_spatial_value_m',real(sum(C(stationary))), ...
    'strict_spatial_zero_excluded',true,'temporal_zero_not_automatically_excluded',true);
if ~isempty(bins)
    bins=bins(:);N=numel(t);dt=mean(diff(t));
    assert(max(abs(t-(0:N-1)'*dt))<1e-8 && max(abs(omega-2*pi*bins/(N*dt)))<1e-10);
    index=mod(bins-bins.',N)+1;eta=real(fft(accumarray(index(:),C(:),[N,1])));
else
    eta=complex(zeros(size(t)));C=C(:);frequency=frequency(:);
    for start=1:512:numel(C)
        ids=start:min(start+511,numel(C));eta=eta+exp(-1i*t*frequency(ids).')*C(ids);
    end
    eta=real(eta);
end
assert(all(isfinite(eta)));
end
