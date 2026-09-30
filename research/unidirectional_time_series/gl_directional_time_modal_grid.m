function [eta,audit,state]=gl_directional_time_modal_grid( ...
        A,omega,kx,ky,g,h,kp,t,options)
%GL_DIRECTIONAL_TIME_MODAL_GRID Joint-modal finite-rank GL eta22 at a probe.
% Non-Cartesian parent wavevectors are conservatively deposited on a
% Cartesian modal grid. The Green--Laplace resolvent itself is not fitted.
arguments
    A (:,1) {mustBeNumeric}
    omega (:,1) double {mustBePositive,mustBeFinite}
    kx (:,1) double {mustBeFinite}
    ky (:,1) double {mustBeFinite}
    g (1,1) double {mustBePositive,mustBeFinite}
    h (1,1) double {mustBePositive,mustBeFinite}
    kp (1,1) double {mustBePositive,mustBeFinite}
    t (:,1) double {mustBeFinite}
    options (1,1) struct = struct
end

options=with_defaults(options,struct( ...
    'Nx',64,'Ny',64,'J',8,'QxLimit',0,'QyLimit',0,'Baseband',true));
mustBeInteger(options.Nx);mustBePositive(options.Nx);
mustBeInteger(options.Ny);mustBePositive(options.Ny);
mustBeInteger(options.J);mustBePositive(options.J);
mustBeNonnegative(options.QxLimit);mustBeFinite(options.QxLimit);
mustBeNonnegative(options.QyLimit);mustBeFinite(options.QyLimit);
mustBeNumericOrLogical(options.Baseband);

if ~isequal(numel(A),numel(omega),numel(kx),numel(ky))
    error('A, omega, kx and ky must contain the same number of modes.');
end
if options.Nx>1 && mod(options.Nx,2)~=0
    error('Nx must be even unless Nx=1.');
end
if options.Ny>1 && mod(options.Ny,2)~=0
    error('Ny must be even unless Ny=1.');
end

N=numel(t);
if N<3,error('At least three uniformly spaced samples are required.');end
dt=mean(diff(t));
if max(abs(t-t(1)-(0:N-1)'*dt))>1e-9*max(1,abs(t(end)-t(1)))
    error('t must be uniformly spaced.');
end
bins=round(omega*N*dt/(2*pi));
if any(bins<1) || max(abs(omega-2*pi*bins/(N*dt)))>1e-10*max(1,max(omega))
    error('omega must lie on positive temporal DFT bins of the record.');
end
if options.Baseband,base_bin=min(bins);else,base_bin=0;end
work_bins=bins-base_bin;
maximum_work_bin=max(work_bins);
work_count=2*maximum_work_bin+1;
if work_count>N
    error('The positive-sum output band does not fit in the record DFT.');
end

qx=h*kx;
qy=h*ky;
q=hypot(qx,qy);
nu=omega*sqrt(h/g);
if any(q<=0) || any(qx<=0)
    error('Strict-forward, nonzero parent modes are required.');
end

[qx_axis,dqx]=modal_axis(options.Nx,options.QxLimit,qx);
[qy_axis,dqy]=modal_axis(options.Ny,options.QyLimit,qy);
[joint_shifted,projection]=deposit_joint( ...
    A,work_bins,qx,qy,nu,qx_axis,qy_axis,dqx,dqy,work_count);
joint=ifftshift(ifftshift(joint_shifted,1),2);
qx_fft=ifftshift(qx_axis);
qy_fft=ifftshift(qy_axis);
[QX,QY]=meshgrid(qx_fft,qy_fft);
Q=hypot(QX,QY);

frequency_step=2*pi/(N*dt)*sqrt(h/g);
nu_axis=(base_bin+(0:work_count-1))*frequency_step;
NU=reshape(nu_axis,1,1,[]);
safe_nu=NU;
safe_nu(safe_nu==0)=1;
QX3=reshape(QX,options.Ny,options.Nx,1);
QY3=reshape(QY,options.Ny,options.Nx,1);
Q23=QX3.^2+QY3.^2;

v=joint_field(joint,ones(size(joint)),options.Nx,options.Ny);
vn2=joint_field(joint,NU.^2,options.Nx,options.Ny);
source_d=2*v.*vn2;
clear vn2
vn=joint_field(joint,NU,options.Nx,options.Ny);
source_d=source_d+vn.^2;
clear vn

hx=joint_field(joint,QX3./safe_nu,options.Nx,options.Ny);
jx=joint_field(joint,QX3,options.Nx,options.Ny);
source_d=source_d-hx.^2;
source_k=2*hx.*jx;
clear hx jx
hy=joint_field(joint,QY3./safe_nu,options.Nx,options.Ny);
jy=joint_field(joint,QY3,options.Nx,options.Ny);
source_d=source_d-hy.^2;
source_k=source_k+2*hy.*jy;
clear hy jy
radial=joint_field(joint,Q23./safe_nu,options.Nx,options.Ny);
source_k=source_k+2*v.*radial;
clear radial v

source_d_hat=fft2(ifft(source_d,[],3))/(options.Nx*options.Ny);
source_k_hat=fft2(ifft(source_k,[],3))/(options.Nx*options.Ny);
clear source_d source_k

[nodes,weights]=gauss_laguerre_rule(options.J);
qp=kp*h;
lambda=sqrt(4*qp*tanh(qp)-2*qp*tanh(2*qp));
if ~isfinite(lambda) || lambda<=0
    error('Invalid determinant-root Green--Laplace scale.');
end
Aq=Q.*tanh(Q);
a=sqrt(Aq);
outbins_work=(2*min(work_bins)):(2*maximum_work_bin);
outbins=2*base_bin+outbins_work;
coefficients=complex(zeros(N,1));
invalid_fraction=zeros(numel(outbins_work),1);
output_source_energy=zeros(numel(outbins_work),1);
invalid_source_energy=zeros(numel(outbins_work),1);
minimum_gap=Inf;
for output_index=1:numel(outbins_work)
    output_work_bin=outbins_work(output_index);
    output_bin=2*base_bin+output_work_bin;
    s=output_bin*frequency_step;
    sd=source_d_hat(:,:,output_work_bin+1);
    sk=source_k_hat(:,:,output_work_bin+1);
    source_energy=abs(sd).^2+abs(sk).^2;
    valid=s>a;
    output_source_energy(output_index)=sum(source_energy,'all');
    invalid_source_energy(output_index)=sum(source_energy(~valid),'all');
    invalid_fraction(output_index)=invalid_source_energy(output_index) ...
        /max(output_source_energy(output_index),realmin);
    if any(valid,'all')
        minimum_gap=min(minimum_gap,min(s-a(valid)));
    end
    md=zeros(size(Q));
    mk=zeros(size(Q));
    av=a(valid);
    for node_index=1:options.J
        tau=nodes(node_index)/lambda;
        ep=weights(node_index)/lambda ...
            *exp(nodes(node_index)-(s-av)*tau);
        em=weights(node_index)/lambda ...
            *exp(nodes(node_index)-(s+av)*tau);
        md(valid)=md(valid)-av.*(ep-em)/2;
        mk(valid)=mk(valid)+(ep+em)/2;
    end
    response=complex(zeros(size(Q)));
    response(valid)=(md(valid).*sd(valid)+mk(valid).*sk(valid))/(4*h);
    coefficients(output_bin+1)=sum(response,'all');
end
eta=fft(coefficients);

station_input=squeeze(sum(sum(joint,1),2));
original_input=accumarray(work_bins+1,A,[work_count,1],@sum,0);
station_input_error=norm(station_input-original_input) ...
    /max(norm(original_input),realmin);
significant=output_source_energy>1e-10*max(output_source_energy);
if any(significant)
    significant_invalid=max(invalid_fraction(significant));
else
    significant_invalid=0;
end
audit=struct( ...
    'status','joint-modal-cartesian-green-laplace-prototype', ...
    'pair_loops',0,'time_snapshot_loops',0, ...
    'resolvent_approximation','Gauss-Laguerre Laplace nodes only', ...
    'wavevector_projection','bilinear conservative Cartesian deposit', ...
    'quadrature_rank',options.J,'quadrature_exact_moment_degree',2*options.J-1, ...
    'input_mode_count',numel(A),'input_bin_count',numel(unique(bins)), ...
    'baseband_enabled',logical(options.Baseband),'base_bin',base_bin, ...
    'work_frequency_count',work_count,'output_bins',outbins, ...
    'modal_grid',[options.Ny,options.Nx], ...
    'qx_limits',[qx_axis(1),qx_axis(end)], ...
    'qy_limits',[qy_axis(1),qy_axis(end)], ...
    'station_input_relative',station_input_error, ...
    'wavevector_projection_rms_relative',projection.q_rms_relative, ...
    'wavevector_projection_max_relative',projection.q_max_relative, ...
    'dispersion_residual_rms_relative',projection.dispersion_rms_relative, ...
    'dispersion_residual_max_relative',projection.dispersion_max_relative, ...
    'maximum_invalid_laplace_source_fraction',max(invalid_fraction), ...
    'significant_bin_invalid_laplace_source_fraction',significant_invalid, ...
    'energy_weighted_invalid_laplace_source_fraction', ...
        sum(invalid_source_energy)/max(sum(output_source_energy),realmin), ...
    'minimum_valid_s_minus_a',minimum_gap, ...
    'aliasing_control','none', ...
    'fitted_coefficients',0,'polynomial_resolvent',false);
state=struct('coefficients',coefficients,'qx_axis',qx_axis, ...
    'qy_axis',qy_axis,'joint_station_spectrum',station_input, ...
    'invalid_laplace_source_fraction',invalid_fraction, ...
    'projection',projection);
end

function [axis,dq]=modal_axis(count,requested_limit,parent)
if count==1
    if any(abs(parent)>100*eps(max(1,max(abs(parent)))))
        error('A one-point modal axis can represent only zero components.');
    end
    axis=0;dq=Inf;return
end
if requested_limit==0
    parent_max=max(abs(parent));
    requested_limit=2.05*parent_max/(1-2/count);
end
dq=2*requested_limit/count;
axis=(-count/2:count/2-1)*dq;
if min(parent)<axis(1) || max(parent)>axis(end)
    error('A parent wavevector lies outside the requested modal axis.');
end
end

function [joint,report]=deposit_joint( ...
        A,bins,qx,qy,nu,qx_axis,qy_axis,dqx,dqy,work_count)
nx=numel(qx_axis);ny=numel(qy_axis);mode_count=numel(A);
[ix0,ix1,wx0,wx1]=bracket(qx,qx_axis,dqx);
[iy0,iy1,wy0,wy1]=bracket(qy,qy_axis,dqy);
indices=cell(4,1);values=cell(4,1);mass=cell(4,1);
qerror=cell(4,1);derror=cell(4,1);
corners=[0,0;1,0;0,1;1,1];
for corner=1:4
    if corners(corner,1)==0,ix=ix0;wx=wx0;else,ix=ix1;wx=wx1;end
    if corners(corner,2)==0,iy=iy0;wy=wy0;else,iy=iy1;wy=wy1;end
    weight=wx.*wy;
    indices{corner}=iy+(ix-1)*ny+bins*ny*nx;
    values{corner}=A.*weight;
    mass{corner}=abs(A).^2.*weight;
    qxnode=reshape(qx_axis(ix),[],1);
    qynode=reshape(qy_axis(iy),[],1);
    qnode=hypot(qxnode,qynode);
    qtarget=hypot(qx,qy);
    qerror{corner}=hypot(qxnode-qx,qynode-qy) ...
        ./max(qtarget,realmin);
    derror{corner}=abs(nu.^2-qnode.*tanh(qnode)) ...
        ./max(nu.^2,realmin);
end
all_indices=vertcat(indices{:});
all_values=vertcat(values{:});
joint=reshape(accumarray(all_indices,all_values,[ny*nx*work_count,1], ...
    @sum,complex(0)),ny,nx,work_count);
all_mass=vertcat(mass{:});
all_qerror=vertcat(qerror{:});
all_derror=vertcat(derror{:});
active=all_mass>0;
total_mass=sum(all_mass(active));
report=struct( ...
    'q_rms_relative',sqrt(sum(all_mass(active).*all_qerror(active).^2) ...
        /max(total_mass,realmin)), ...
    'q_max_relative',max(all_qerror(active)), ...
    'dispersion_rms_relative',sqrt(sum(all_mass(active).*all_derror(active).^2) ...
        /max(total_mass,realmin)), ...
    'dispersion_max_relative',max(all_derror(active)), ...
    'deposited_corner_count',4*mode_count);
end

function [i0,i1,w0,w1]=bracket(value,axis,dq)
if numel(axis)==1
    i0=ones(size(value));i1=i0;w0=ones(size(value));w1=zeros(size(value));
    return
end
position=(value-axis(1))/dq+1;
i0=floor(position);
i0=max(1,min(numel(axis)-1,i0));
i1=i0+1;
w1=position-i0;
w1=max(0,min(1,w1));
w0=1-w1;
end

function field=joint_field(joint,multiplier,nx,ny)
spectrum=joint.*multiplier;
field=fft(ifft(ifft(spectrum,[],1),[],2)*(nx*ny),[],3);
end

function [nodes,weights]=gauss_laguerre_rule(J)
diagonal=2*(1:J)-1;
off_diagonal=1:(J-1);
[vectors,values]=eig(diag(diagonal)+diag(off_diagonal,1) ...
    +diag(off_diagonal,-1),'vector');
[nodes,order]=sort(values.');
weights=vectors(1,order).^2;
for degree=0:(2*J-1)
    relative_error=abs(sum(weights.*nodes.^degree)-factorial(degree)) ...
        /max(1,factorial(degree));
    if relative_error>2e-11
        error('Gauss--Laguerre moment gate failed at degree %d.',degree);
    end
end
end

function options=with_defaults(options,defaults)
names=fieldnames(defaults);
for index=1:numel(names)
    name=names{index};
    if ~isfield(options,name),options.(name)=defaults.(name);end
end
unknown=setdiff(fieldnames(options),names);
if ~isempty(unknown)
    error('Unknown option: %s.',unknown{1});
end
end
