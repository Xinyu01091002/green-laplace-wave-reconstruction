function [eta, audit] = gl_eta22_time_pairs(A,omega,k,h,kp,t,J)
%GL_ETA22_TIME_PAIRS Ordered-pair evaluation of the released GL eta22 kernel.
% A is physical amplitude in Re sum A exp(-i omega t), k>0; eta is in m.
% No spatial-grid approximation, Stokes repair, or reference-field input.
A=A(:); omega=omega(:); k=k(:); t=t(:);
assert(numel(A)==numel(k) && numel(k)==numel(omega));
assert(all(isfinite([real(A);imag(A);omega;k;t])));
assert(h>0 && kp>0 && all(k*h>=0.3) && all(omega>0));
assert(ismember(J,[4,6,8,12]));
q=h*k; nu=sqrt(q.*tanh(q));
assert(max(abs(omega./nu-omega(1)/nu(1)))<1e-10*max(1,omega(1)/nu(1)));
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
freeze=jsondecode(fileread(fullfile(root,'symbolic','generated', ...
    'finite_depth_directional_order2_eta22_pure_gl8.json')));
assert(freeze.overall_exact_gate_pass && ~freeze.stokes_correction_used ...
    && ~freeze.angular_correction_used && ~freeze.oracle_or_mf12_used_in_formula);
[V,D]=eig(diag(2*(1:J)-1)+diag(1:J-1,1)+diag(1:J-1,-1),'vector');
[nodes,idx]=sort(D); weights=V(1,idx).^2;
qp=h*kp; lambda=sqrt(4*qp*tanh(qp)-2*qp*tanh(2*qp));
assert(lambda>0);
Q=q+q.'; a=sqrt(Q.*tanh(Q)); s=nu+nu.';
sd=nu.^2+nu*nu.'+(nu.').^2-(q*q.')./(nu*nu.');
sk=(q.^2+q*q.')./nu+((q.').^2+q*q.')./nu.';
assert(all(s>a,'all'),'Only regular positive sum interactions are supported.');
kernel=zeros(size(Q));
for j=1:J
    tau=nodes(j)/lambda;
    % Algebraically identical to exp(x-s*tau)*sinh/cosh(a*tau),
    % avoiding overflow before exponential cancellation.
    ep=exp(nodes(j)-(s-a)*tau); em=exp(nodes(j)-(s+a)*tau);
    kernel=kernel+weights(j)/(8*lambda)*((-a.*sd+sk).*ep+(a.*sd+sk).*em);
end
C=(A*A.').*kernel/h; eta=complex(zeros(size(t)));
for j=1:numel(A)
    eta=eta+exp(-1i*t*(omega(j)+omega.'))*C(j,:).';
end
assert(all(isfinite(eta)));
audit=struct('rank',J,'parent_count',numel(A),'ordered_pairs',numel(A)^2, ...
    'lambda',lambda,'minimum_parent_kh',min(q),'no_stokes_correction',true, ...
    'output','physical analytic eta22; take real part','information_boundary','eta1 only');
end
