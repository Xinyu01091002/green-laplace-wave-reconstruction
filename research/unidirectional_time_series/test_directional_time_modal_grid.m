function report=test_directional_time_modal_grid
%TEST_DIRECTIONAL_TIME_MODAL_GRID Single-mode and off-grid prototype gates.
q=1;g=1;h=1;kp=1;nu=sqrt(q*tanh(q));
N=128;input_bin=4;dt=2*pi*input_bin/(N*nu);t=(0:N-1)'*dt;
A=0.01+0.003i;
options=struct('Nx',16,'Ny',1,'J',8,'QxLimit',4,'QyLimit',0);
[eta,audit,state]=gl_directional_time_modal_grid( ...
    A,nu,q,0,g,h,kp,t,options);
absolute_options=options;absolute_options.Baseband=false;
[eta_absolute,absolute_audit]=gl_directional_time_modal_grid( ...
    A,nu,q,0,g,h,kp,t,absolute_options);
baseband_relative=norm(eta-eta_absolute)/norm(eta_absolute);
assert(baseband_relative<2e-12);
assert(audit.work_frequency_count==1 && absolute_audit.work_frequency_count==9);

[nodes,weights]=local_gauss_laguerre(options.J);
lambda=sqrt(4*kp*tanh(kp)-2*kp*tanh(2*kp));
Q=2*q;a=sqrt(Q*tanh(Q));s=2*nu;
md=0;mk=0;
for node_index=1:options.J
    tau=nodes(node_index)/lambda;
    ep=weights(node_index)/lambda*exp(nodes(node_index)-(s-a)*tau);
    em=weights(node_index)/lambda*exp(nodes(node_index)-(s+a)*tau);
    md=md-a*(ep-em)/2;
    mk=mk+(ep+em)/2;
end
Sd=(3*nu^2-q^2/nu^2)*A^2;
Sk=4*q^2/nu*A^2;
expected=(md*Sd+mk*Sk)/(4*h);
actual=state.coefficients(2*input_bin+1);
single_mode_relative=abs(actual-expected)/max(abs(expected),realmin);
assert(single_mode_relative<2e-12);
assert(audit.station_input_relative<2e-14);
assert(audit.pair_loops==0 && ~audit.polynomial_resolvent);

q_off=1.07;nu_off=sqrt(q_off*tanh(q_off));
dt_off=2*pi*input_bin/(N*nu_off);t_off=(0:N-1)'*dt_off;
[~,off_audit]=gl_directional_time_modal_grid( ...
    A,nu_off,q_off,0,g,h,kp,t_off,options);
assert(off_audit.station_input_relative<2e-14);
assert(off_audit.wavevector_projection_rms_relative>0);
assert(off_audit.dispersion_residual_rms_relative>0);
assert(all(isfinite(eta)));

frequency_step=.1;multi_bins=(4:6)';multi_omega=frequency_step*multi_bins;
multi_q=arrayfun(@(w)fzero(@(x)x*tanh(x)-w^2,[0,2]),multi_omega);
multi_t=(0:N-1)'*(2*pi/(N*frequency_step));
multi_A=[.01+.002i;-.004+.003i;.006-.001i];
multi_options=struct('Nx',32,'Ny',1,'J',8,'QxLimit',4,'QyLimit',0);
multi_absolute=multi_options;multi_absolute.Baseband=false;
multi_baseband=gl_directional_time_modal_grid( ...
    multi_A,multi_omega,multi_q,zeros(size(multi_q)), ...
    g,h,kp,multi_t,multi_options);
multi_native=gl_directional_time_modal_grid( ...
    multi_A,multi_omega,multi_q,zeros(size(multi_q)), ...
    g,h,kp,multi_t,multi_absolute);
multifrequency_baseband_relative=norm(multi_baseband-multi_native) ...
    /norm(multi_native);
assert(multifrequency_baseband_relative<2e-12);

report=struct( ...
    'status','pass', ...
    'single_mode_relative',single_mode_relative, ...
    'baseband_relative',baseband_relative, ...
    'on_grid_input_relative',audit.station_input_relative, ...
    'off_grid_input_relative',off_audit.station_input_relative, ...
    'off_grid_wavevector_rms_relative', ...
        off_audit.wavevector_projection_rms_relative, ...
    'off_grid_dispersion_rms_relative', ...
        off_audit.dispersion_residual_rms_relative, ...
    'multifrequency_baseband_relative',multifrequency_baseband_relative, ...
    'eta_is_finite',all(isfinite(eta)));
disp(report)
end

function [nodes,weights]=local_gauss_laguerre(J)
diagonal=2*(1:J)-1;off_diagonal=1:(J-1);
[vectors,values]=eig(diag(diagonal)+diag(off_diagonal,1) ...
    +diag(off_diagonal,-1),'vector');
[nodes,order]=sort(values.');weights=vectors(1,order).^2;
end
