function report=test_directional_time_modal_grid_eta33
%TEST_DIRECTIONAL_TIME_MODAL_GRID_ETA33 Compare to the frozen spatial GL graph.
q=1;g=1;h=1;kp=1;nu=sqrt(q*tanh(q));
N=128;input_bin=4;dt=2*pi*input_bin/(N*nu);t=(0:N-1)'*dt;
A=.01;
grid_options=struct('Nx',16,'Ny',1,'J',8,'QxLimit',4,'QyLimit',0);
[candidate,audit]=gl_directional_time_modal_grid_eta33( ...
    A,nu,q,0,g,h,kp,t,grid_options);
absolute_options=grid_options;absolute_options.Baseband=false;
[absolute_candidate,absolute_audit]=gl_directional_time_modal_grid_eta33( ...
    A,nu,q,0,g,h,kp,t,absolute_options);
baseband_relative=norm(candidate-absolute_candidate)/norm(absolute_candidate);
assert(baseband_relative<5e-12);
assert(audit.work_frequency_count==1 && absolute_audit.work_frequency_count==13);
mode=[0:7,-8:-1];qx=.5*mode;qy=zeros(size(qx));
spectrum=complex(zeros(size(qx)));spectrum(3)=A*numel(mode);
project_root=string(fileparts(fileparts(fileparts(mfilename('fullpath')))));
[reference,reference_audit]=gl_no_stokes_eta33( ...
    spectrum,qx,qy,kp,8,project_root,8);
relative=abs(candidate(1)-reference(1))/abs(reference(1));
assert(relative<5e-12);
assert(audit.pair_loops==0 && audit.triple_loops==0);
assert(~audit.polynomial_resolvent && ~audit.stokes_correction);

frequency_step=.1;multi_bins=(4:6)';multi_omega=frequency_step*multi_bins;
multi_q=arrayfun(@(w)fzero(@(x)x*tanh(x)-w^2,[0,2]),multi_omega);
multi_t=(0:N-1)'*(2*pi/(N*frequency_step));
multi_A=[.01+.002i;-.004+.003i;.006-.001i];
multi_options=struct('Nx',32,'Ny',1,'J',8,'QxLimit',4,'QyLimit',0);
multi_absolute=multi_options;multi_absolute.Baseband=false;
multi_baseband=gl_directional_time_modal_grid_eta33( ...
    multi_A,multi_omega,multi_q,zeros(size(multi_q)), ...
    g,h,kp,multi_t,multi_options);
multi_native=gl_directional_time_modal_grid_eta33( ...
    multi_A,multi_omega,multi_q,zeros(size(multi_q)), ...
    g,h,kp,multi_t,multi_absolute);
multifrequency_baseband_relative=norm(multi_baseband-multi_native) ...
    /norm(multi_native);
assert(multifrequency_baseband_relative<5e-12);
report=struct('status','pass','relative',relative, ...
    'baseband_relative',baseband_relative, ...
    'multifrequency_baseband_relative',multifrequency_baseband_relative, ...
    'station_input_relative',audit.station_input_relative, ...
    'reference_pair_loops',reference_audit.pair_loops, ...
    'reference_triple_loops',reference_audit.triple_loops);
disp(report)
end
