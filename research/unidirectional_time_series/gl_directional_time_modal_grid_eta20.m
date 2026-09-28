function [eta20,audit,state]=gl_directional_time_modal_grid_eta20( ...
        A,omega,kx,ky,g,h,t,options)
%GL_DIRECTIONAL_TIME_MODAL_GRID_ETA20 Nonzero-K difference-frequency GL eta20.
% Strict spatial Q=0 is excluded. Same-frequency, nonzero-Q directional
% differences are retained. No parent-pair loop is used.
arguments
    A (:,1) {mustBeNumeric}
    omega (:,1) double {mustBePositive,mustBeFinite}
    kx (:,1) double {mustBeFinite}
    ky (:,1) double {mustBeFinite}
    g (1,1) double {mustBePositive,mustBeFinite}
    h (1,1) double {mustBePositive,mustBeFinite}
    t (:,1) double {mustBeFinite}
    options (1,1) struct = struct
end

options=with_defaults(options,struct( ...
    'Nx',256,'Ny',256,'J',16,'QxLimit',0,'QyLimit',0));
mustBeInteger(options.Nx);mustBePositive(options.Nx);
mustBeInteger(options.Ny);mustBePositive(options.Ny);
if ~ismember(options.J,[6,12,16])
    error('The frozen shared-scale eta20 ranks are GL6, GL12 and GL16.');
end
mustBeNonnegative(options.QxLimit);mustBeFinite(options.QxLimit);
mustBeNonnegative(options.QyLimit);mustBeFinite(options.QyLimit);
if ~isequal(numel(A),numel(omega),numel(kx),numel(ky))
    error('A, omega, kx and ky must contain the same number of modes.');
end
if options.Nx<16 || options.Ny<16 ...
        || mod(options.Nx,2)~=0 || mod(options.Ny,2)~=0
    error('Nx and Ny must be even and at least 16.');
end

N=numel(t);dt=mean(diff(t));
if max(abs(t-t(1)-(0:N-1)'*dt))>1e-9*max(1,abs(t(end)-t(1)))
    error('t must be uniformly spaced.');
end
bins=round(omega*N*dt/(2*pi));
if any(bins<1) || max(abs(omega-2*pi*bins/(N*dt)))>1e-10*max(1,max(omega))
    error('omega must lie on positive temporal DFT bins of the record.');
end
base_bin=min(bins);work_bins=bins-base_bin;span=max(work_bins);
work_count=2*span+1;

qx=h*kx;qy=h*ky;q=hypot(qx,qy);nu=omega*sqrt(h/g);
if any(q<=0) || any(qx<=0)
    error('Strict-forward, nonzero parent modes are required.');
end
[qx_axis,dqx]=modal_axis(options.Nx,options.QxLimit,qx,2);
[qy_axis,dqy]=modal_axis(options.Ny,options.QyLimit,qy,2);
[joint_shifted,projection]=deposit_joint( ...
    A,work_bins,qx,qy,nu,qx_axis,qy_axis,dqx,dqy,work_count);
joint=ifftshift(ifftshift(joint_shifted,1),2);
qx_fft=ifftshift(qx_axis);qy_fft=ifftshift(qy_axis);
[QX,QY]=meshgrid(qx_fft,qy_fft);Q=hypot(QX,QY);
a=sqrt(Q.*tanh(Q));QX3=reshape(QX,options.Ny,options.Nx,1);
QY3=reshape(QY,options.Ny,options.Nx,1);Q23=QX3.^2+QY3.^2;
frequency_step=2*pi/(N*dt)*sqrt(h/g);
input_nu_axis=(base_bin+(0:work_count-1))*frequency_step;
NU=reshape(input_nu_axis,1,1,[]);safe_nu=NU;safe_nu(safe_nu==0)=1;

one=joint_field(joint,1,options.Nx,options.Ny);
dno=joint_field(joint,NU.^2,options.Nx,options.Ny);
nu_field=joint_field(joint,NU,options.Nx,options.Ny);
bx=joint_field(joint,QX3./safe_nu,options.Nx,options.Ny);
by=joint_field(joint,QY3./safe_nu,options.Nx,options.Ny);
qx_field=joint_field(joint,QX3,options.Nx,options.Ny);
qy_field=joint_field(joint,QY3,options.Nx,options.Ny);
radial=joint_field(joint,Q23./safe_nu,options.Nx,options.Ny);
forcing_d=-dno.*conj(one)-one.*conj(dno) ...
    +nu_field.*conj(nu_field)+bx.*conj(bx)+by.*conj(by);
forcing_k=1i*radial.*conj(one)-1i*bx.*conj(qx_field) ...
    -1i*by.*conj(qy_field)+1i*qx_field.*conj(bx) ...
    +1i*qy_field.*conj(by)-1i*one.*conj(radial);
clear one dno nu_field bx by qx_field qy_field radial
fd_hat=joint_spectrum(forcing_d,options.Nx,options.Ny);
fk_hat=joint_spectrum(forcing_k,options.Nx,options.Ny);
clear forcing_d forcing_k

energy=abs(A).^2;energy_sum=sum(energy);
mean_qx=sum(qx.*energy)/energy_sum;mean_qy=sum(qy.*energy)/energy_sum;
delta_q=sqrt(2*sum(((qx-mean_qx).^2+(qy-mean_qy).^2).*energy)/energy_sum);
if delta_q<=64*eps(max(1,max(q)))
    eta20=zeros(N,1);audit=degenerate_audit(options,work_count,projection);
    state=struct('coefficients',zeros(N,1),'delta_q',delta_q);return
end
[nodes,weights]=gauss_laguerre_rule(options.J);
signed_bins=[0:span,-span:-1];coefficients=complex(zeros(N,1));
total_source=0;strict_zero_source=0;nonzero_invalid_source=0;
for index=1:work_count
    signed_bin=signed_bins(index);sigma=signed_bin*frequency_step;
    fd=fd_hat(:,:,index);fk=fk_hat(:,:,index);
    source_energy=abs(fd).^2+abs(fk).^2;
    valid=Q>0 & a>abs(sigma);
    total_source=total_source+sum(source_energy,'all');
    strict_zero_source=strict_zero_source+sum(source_energy(Q==0),'all');
    nonzero_invalid=Q>0 & a<=abs(sigma);
    nonzero_invalid_source=nonzero_invalid_source ...
        +sum(source_energy(nonzero_invalid),'all');
    md=zeros(size(Q));mk=md;av=a(valid);
    for node_index=1:options.J
        tau=nodes(node_index)/delta_q;
        ep=weights(node_index)/delta_q ...
            *exp(nodes(node_index)-(av-sigma)*tau);
        em=weights(node_index)/delta_q ...
            *exp(nodes(node_index)-(av+sigma)*tau);
        md(valid)=md(valid)-av.*(ep+em)/4;
        mk(valid)=mk(valid)+1i*(ep-em)/4;
    end
    response=complex(zeros(size(Q)));
    response(valid)=(md(valid).*fd(valid)+mk(valid).*fk(valid))/(2*h);
    value=sum(response,'all');
    if signed_bin>=0
        native_index=signed_bin+1;
    else
        native_index=N+signed_bin+1;
    end
    coefficients(native_index)=value;
end
coefficients=hermitian_project_1d(coefficients);
eta20=real(fft(coefficients));

station_input=squeeze(sum(sum(joint,1),2));
original_input=accumarray(work_bins+1,A,[work_count,1],@sum,0);
audit=struct( ...
    'status','joint-modal-cartesian-shared-scale-green-laplace-eta20-prototype', ...
    'pair_loops',0,'time_snapshot_loops',0,'quadrature_rank',options.J, ...
    'strict_zero_spatial_mode','excluded','temporal_zero_nonzero_Q','retained', ...
    'base_bin',base_bin,'work_frequency_count',work_count, ...
    'difference_bins',signed_bins,'modal_grid',[options.Ny,options.Nx], ...
    'delta_q_rms_pair',delta_q, ...
    'station_input_relative',norm(station_input-original_input) ...
        /max(norm(original_input),realmin), ...
    'wavevector_projection_rms_relative',projection.q_rms_relative, ...
    'dispersion_residual_rms_relative',projection.dispersion_rms_relative, ...
    'strict_zero_source_energy_fraction', ...
        strict_zero_source/max(total_source,realmin), ...
    'nonzero_invalid_laplace_source_energy_fraction', ...
        nonzero_invalid_source/max(total_source,realmin), ...
    'fitted_coefficients',0,'polynomial_resolvent',false);
state=struct('coefficients',coefficients,'delta_q',delta_q, ...
    'qx_axis',qx_axis,'qy_axis',qy_axis,'projection',projection);
end

function audit=degenerate_audit(options,work_count,projection)
audit=struct('status','degenerate-zero-difference-output','pair_loops',0, ...
    'time_snapshot_loops',0,'quadrature_rank',options.J, ...
    'strict_zero_spatial_mode','excluded','temporal_zero_nonzero_Q','retained', ...
    'work_frequency_count',work_count,'modal_grid',[options.Ny,options.Nx], ...
    'wavevector_projection_rms_relative',projection.q_rms_relative, ...
    'dispersion_residual_rms_relative',projection.dispersion_rms_relative, ...
    'fitted_coefficients',0,'polynomial_resolvent',false);
end

function coefficients=hermitian_project_1d(coefficients)
N=numel(coefficients);
for bin=1:floor((N-1)/2)
    positive=bin+1;negative=N-bin+1;
    average=(coefficients(positive)+conj(coefficients(negative)))/2;
    coefficients(positive)=average;coefficients(negative)=conj(average);
end
coefficients(1)=real(coefficients(1));
if mod(N,2)==0,coefficients(N/2+1)=real(coefficients(N/2+1));end
end

function field=joint_field(joint,multiplier,nx,ny)
field=fft(ifft(ifft(joint.*multiplier,[],1),[],2)*(nx*ny),[],3);
end

function spectrum=joint_spectrum(field,nx,ny)
spectrum=fft2(ifft(field,[],3))/(nx*ny);
end

function [axis,dq]=modal_axis(count,requested_limit,parent,order)
if requested_limit==0
    parent_max=max(abs(parent));
    if parent_max==0,error('A nonzero modal-axis limit is required.');end
    requested_limit=(order+.05)*parent_max/(1-2/count);
end
dq=2*requested_limit/count;axis=(-count/2:count/2-1)*dq;
if min(parent)<axis(1) || max(parent)>axis(end)
    error('A parent wavevector lies outside the requested modal axis.');
end
end

function [joint,report]=deposit_joint( ...
        A,bins,qx,qy,nu,qx_axis,qy_axis,dqx,dqy,work_count)
nx=numel(qx_axis);ny=numel(qy_axis);
[ix0,ix1,wx0,wx1]=bracket(qx,qx_axis,dqx);
[iy0,iy1,wy0,wy1]=bracket(qy,qy_axis,dqy);
indices=cell(4,1);values=cell(4,1);mass=cell(4,1);
qerror=cell(4,1);derror=cell(4,1);corners=[0,0;1,0;0,1;1,1];
for corner=1:4
    if corners(corner,1)==0,ix=ix0;wx=wx0;else,ix=ix1;wx=wx1;end
    if corners(corner,2)==0,iy=iy0;wy=wy0;else,iy=iy1;wy=wy1;end
    weight=wx.*wy;indices{corner}=iy+(ix-1)*ny+bins*ny*nx;
    values{corner}=A.*weight;mass{corner}=abs(A).^2.*weight;
    qxnode=reshape(qx_axis(ix),[],1);qynode=reshape(qy_axis(iy),[],1);
    qnode=hypot(qxnode,qynode);qtarget=hypot(qx,qy);
    qerror{corner}=hypot(qxnode-qx,qynode-qy)./max(qtarget,realmin);
    derror{corner}=abs(nu.^2-qnode.*tanh(qnode))./max(nu.^2,realmin);
end
joint=reshape(accumarray(vertcat(indices{:}),vertcat(values{:}), ...
    [ny*nx*work_count,1],@sum,complex(0)),ny,nx,work_count);
all_mass=vertcat(mass{:});all_qerror=vertcat(qerror{:});
all_derror=vertcat(derror{:});active=all_mass>0;total_mass=sum(all_mass(active));
report=struct( ...
    'q_rms_relative',sqrt(sum(all_mass(active).*all_qerror(active).^2) ...
        /max(total_mass,realmin)), ...
    'dispersion_rms_relative',sqrt(sum(all_mass(active).*all_derror(active).^2) ...
        /max(total_mass,realmin)));
end

function [i0,i1,w0,w1]=bracket(value,axis,dq)
position=(value-axis(1))/dq+1;i0=floor(position);
i0=max(1,min(numel(axis)-1,i0));i1=i0+1;
w1=max(0,min(1,position-i0));w0=1-w1;
end

function [nodes,weights]=gauss_laguerre_rule(J)
diagonal=2*(1:J)-1;off_diagonal=1:(J-1);
[vectors,values]=eig(diag(diagonal)+diag(off_diagonal,1) ...
    +diag(off_diagonal,-1),'vector');
[nodes,order]=sort(values.');weights=vectors(1,order).^2;
end

function options=with_defaults(options,defaults)
names=fieldnames(defaults);
for index=1:numel(names)
    name=names{index};if ~isfield(options,name),options.(name)=defaults.(name);end
end
unknown=setdiff(fieldnames(options),names);
if ~isempty(unknown),error('Unknown option: %s.',unknown{1});end
end
