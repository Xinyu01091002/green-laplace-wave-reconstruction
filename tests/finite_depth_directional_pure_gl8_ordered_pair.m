function eta22_plus = finite_depth_directional_pure_gl8_ordered_pair( ...
    eta11_spectrum,qx,qy,peak_depth_qp,project_root)
%FINITE_DEPTH_DIRECTIONAL_PURE_GL8_ORDERED_PAIR Direct frozen-kernel check.
arguments
    eta11_spectrum (:,:) {mustBeNumeric}
    qx (:,:) double
    qy (:,:) double
    peak_depth_qp (1,1) double {mustBePositive}
    project_root (1,1) string = string(pwd)
end
freeze_path=fullfile(project_root,'symbolic','generated', ...
    'finite_depth_directional_order2_eta22_pure_gl8.json');
freeze=jsondecode(fileread(freeze_path));
if ~freeze.overall_exact_gate_pass || freeze.quadrature_rank~=8 ...
        || freeze.stokes_correction_used || freeze.angular_correction_used
    error('Invalid pure-GL8 Wolfram freeze.');
end
nodes=reshape(freeze.nodes,1,[]);
weights=reshape(freeze.weights,1,[]);
[ny,nx]=size(eta11_spectrum);
scale=nx*ny;
threshold=1e-13*max(1,max(abs(eta11_spectrum(:))));
active=find(abs(eta11_spectrum)>threshold);
[rows,columns]=ind2sub([ny,nx],active);
amplitudes=eta11_spectrum(active)/scale;
px=qx(active);
py=qy(active);
q=hypot(px,py);
nu=sqrt(q.*tanh(q));
lambda=sqrt((2*sqrt(peak_depth_qp*tanh(peak_depth_qp)))^2 ...
    -2*peak_depth_qp*tanh(2*peak_depth_qp));
output=complex(zeros(ny,nx));
for first=1:numel(active)
    for second=1:numel(active)
        dot12=px(first)*px(second)+py(first)*py(second);
        Q=hypot(px(first)+px(second),py(first)+py(second));
        A=Q*tanh(Q);
        a=sqrt(A);
        s=nu(first)+nu(second);
        sd=nu(first)^2+nu(first)*nu(second)+nu(second)^2 ...
            -dot12/(nu(first)*nu(second));
        sk=(q(first)^2+dot12)/nu(first) ...
            +(q(second)^2+dot12)/nu(second);
        kernel=0;
        for index=1:8
            x=nodes(index);
            t=x/lambda;
            if a==0,sinh_over_a=t;else,sinh_over_a=sinh(a*t)/a;end
            coefficient=weights(index)*exp(x)/lambda*exp(-s*t)/4;
            kernel=kernel+coefficient*( ...
                -A*sd*sinh_over_a+sk*cosh(a*t));
        end
        row=mod((rows(first)-1)+(rows(second)-1),ny)+1;
        column=mod((columns(first)-1)+(columns(second)-1),nx)+1;
        output(row,column)=output(row,column) ...
            +scale*amplitudes(first)*amplitudes(second)*kernel;
    end
end
eta22_plus=ifft2(output);
end
