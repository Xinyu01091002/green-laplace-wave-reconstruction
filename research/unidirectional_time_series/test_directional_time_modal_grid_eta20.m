function report=test_directional_time_modal_grid_eta20
%TEST_DIRECTIONAL_TIME_MODAL_GRID_ETA20 Compare to frozen spatial shared GL.
ratio=@(q)sqrt(2*q*tanh(2*q))/sqrt(q*tanh(q));
q1=fzero(@(q)ratio(q)-1.5,[.1,10]);q=[q1;2*q1];
nu=sqrt(q.*tanh(q));g=1;h=1;N=128;bins=[2;3];
dt=2*pi*bins(1)/(N*nu(1));t=(0:N-1)'*dt;
A=[.012+.003i;-.006+.004i];limit=8*q1;
options=struct('Nx',16,'Ny',16,'J',6,'QxLimit',limit,'QyLimit',limit);
[candidate,audit]=gl_directional_time_modal_grid_eta20( ...
    A,nu,q,zeros(size(q)),g,h,t,options);

mode=[0:7,-8:-1];[MX,MY]=meshgrid(mode,mode);
qx=MX*q1;qy=MY*q1;spectrum=complex(zeros(16));
normalization=16*16/2;
spectrum(1,2)=A(1)*normalization;spectrum(1,3)=A(2)*normalization;
spectrum(1,16)=conj(A(1))*normalization;
spectrum(1,15)=conj(A(2))*normalization;
project_root=string(fileparts(fileparts(fileparts(mfilename('fullpath')))));
[reference_half,reference_audit]=eta20_green_laplace_shared( ...
    spectrum,qx,qy,h,"shared6",project_root);
reference=2*reference_half;
relative=abs(candidate(1)-reference(1))/max(abs(reference(1)),realmin);
assert(relative<5e-11);
assert(audit.pair_loops==0 && reference_audit.pair_loops==0);
assert(audit.station_input_relative<2e-14);
report=struct('status','pass','relative',relative, ...
    'station_input_relative',audit.station_input_relative, ...
    'delta_q',audit.delta_q_rms_pair);
disp(report)
end
