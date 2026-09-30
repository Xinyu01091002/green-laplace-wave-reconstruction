function [eta33,audit,state]=gl_directional_time_modal_grid_eta33( ...
        A,omega,kx,ky,g,h,kp,t,options)
%GL_DIRECTIONAL_TIME_MODAL_GRID_ETA33 Nested joint-modal GL eta33 at a probe.
% The inner eta22/Phi22 state and outer cubic response are generated from
% the declared first-order joint spectrum. No pair/triple loops are used.
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
    'Nx',128,'Ny',128,'J',8,'QxLimit',0,'QyLimit',0,'Baseband',true));
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

N=numel(t);dt=mean(diff(t));
if max(abs(t-t(1)-(0:N-1)'*dt))>1e-9*max(1,abs(t(end)-t(1)))
    error('t must be uniformly spaced.');
end
bins=round(omega*N*dt/(2*pi));
if any(bins<1) || max(abs(omega-2*pi*bins/(N*dt)))>1e-10*max(1,max(omega))
    error('omega must lie on positive temporal DFT bins of the record.');
end
if options.Baseband,base_bin=min(bins);else,base_bin=0;end
work_bins=bins-base_bin;maximum_work_bin=max(work_bins);
work_count=3*maximum_work_bin+1;
if work_count>N
    error('The positive cubic output band does not fit in the record DFT.');
end

qx=h*kx;qy=h*ky;q=hypot(qx,qy);nu=omega*sqrt(h/g);
if any(q<=0) || any(qx<=0)
    error('Strict-forward, nonzero parent modes are required.');
end
[qx_axis,dqx]=modal_axis(options.Nx,options.QxLimit,qx,3);
[qy_axis,dqy]=modal_axis(options.Ny,options.QyLimit,qy,3);
[joint_shifted,projection]=deposit_joint( ...
    A,work_bins,qx,qy,nu,qx_axis,qy_axis,dqx,dqy,work_count);
joint=ifftshift(ifftshift(joint_shifted,1),2);
qx_fft=ifftshift(qx_axis);qy_fft=ifftshift(qy_axis);
[QX,QY]=meshgrid(qx_fft,qy_fft);Q2=QX.^2+QY.^2;
Q=sqrt(Q2);Aq=Q.*tanh(Q);a=sqrt(Aq);
QX3=reshape(QX,options.Ny,options.Nx,1);
QY3=reshape(QY,options.Ny,options.Nx,1);
Q23=reshape(Q2,options.Ny,options.Nx,1);
frequency_step=2*pi/(N*dt)*sqrt(h/g);
input_nu_axis=(base_bin+(0:work_count-1))*frequency_step;
inner_nu_axis=(2*base_bin+(0:work_count-1))*frequency_step;
NU=reshape(input_nu_axis,1,1,[]);safe_nu=NU;safe_nu(safe_nu==0)=1;

first=struct();
first.eta=joint_field(joint,1,options.Nx,options.Ny);
first.etax=joint_field(joint,1i*QX3,options.Nx,options.Ny);
first.etay=joint_field(joint,1i*QY3,options.Nx,options.Ny);
first.phix=joint_field(joint,QX3./safe_nu,options.Nx,options.Ny);
first.phiy=joint_field(joint,QY3./safe_nu,options.Nx,options.Ny);
first.phiz=joint_field(joint,-1i*NU,options.Nx,options.Ny);
first.phixz=joint_field(joint,NU.*QX3,options.Nx,options.Ny);
first.phiyz=joint_field(joint,NU.*QY3,options.Nx,options.Ny);
first.phizz=joint_field(joint,-1i*Q23./safe_nu,options.Nx,options.Ny);
first.phizzz=joint_field(joint,-1i*Q23.*NU,options.Nx,options.Ny);
first.phitz=joint_field(joint,-NU.^2,options.Nx,options.Ny);
first.phitzz=joint_field(joint,-Q23,options.Nx,options.Ny);

vn=1i*first.phiz;
vn2=joint_field(joint,NU.^2,options.Nx,options.Ny);
jx=-1i*first.etax;jy=-1i*first.etay;
radial=1i*first.phizz;
source_d=2*first.eta.*vn2+vn.^2-first.phix.^2-first.phiy.^2;
source_k=2*first.eta.*radial ...
    +2*(first.phix.*jx+first.phiy.*jy);
clear vn vn2 jx jy radial
source_d_hat=joint_spectrum(source_d,options.Nx,options.Ny);
source_k_hat=joint_spectrum(source_k,options.Nx,options.Ny);
clear source_d source_k

[nodes,weights]=gauss_laguerre_rule(options.J);
qp=kp*h;nu_p=sqrt(qp*tanh(qp));
lambda2=2*nu_p-sqrt(2*qp*tanh(2*qp));
lambda3=3*nu_p-sqrt(3*qp*tanh(3*qp));
if min(lambda2,lambda3)<=0
    error('Invalid nested Green--Laplace scales.');
end

outbins2_work=(2*min(work_bins)):(2*maximum_work_bin);
E2=complex(zeros(options.Ny,options.Nx,work_count));
P2=E2;inner_total=0;inner_invalid=0;
for output_index=1:numel(outbins2_work)
    output_work_bin=outbins2_work(output_index);
    s=inner_nu_axis(output_work_bin+1);
    sd=source_d_hat(:,:,output_work_bin+1);
    sk=source_k_hat(:,:,output_work_bin+1);
    energy=abs(sd).^2+abs(sk).^2;valid=s>a;
    inner_total=inner_total+sum(energy,'all');
    inner_invalid=inner_invalid+sum(energy(~valid),'all');
    [S,C]=node_response(s,a,valid,lambda2,nodes,weights);
    e=complex(zeros(size(Q)));p=e;
    e(valid)=(-Aq(valid).*S(valid).*sd(valid)+C(valid).*sk(valid))/4;
    p(valid)=1i*(C(valid).*sd(valid)-S(valid).*sk(valid))/4;
    E2(:,:,output_work_bin+1)=e;P2(:,:,output_work_bin+1)=p;
end
clear source_d_hat source_k_hat

pair=struct();
pair.eta=spectral_field(E2,options.Nx,options.Ny);
pair.etax=spectral_field(E2.*(1i*QX3),options.Nx,options.Ny);
pair.etay=spectral_field(E2.*(1i*QY3),options.Nx,options.Ny);
pair.phix=spectral_field(P2.*(1i*QX3),options.Nx,options.Ny);
pair.phiy=spectral_field(P2.*(1i*QY3),options.Nx,options.Ny);
pair.phiz=spectral_field(P2.*reshape(Aq,options.Ny,options.Nx,1), ...
    options.Nx,options.Ny);
pair.phizz=spectral_field(P2.*Q23,options.Nx,options.Ny);
pair.phitz=spectral_field(P2.*reshape(-1i*inner_nu_axis,1,1,[]) ...
    .*reshape(Aq,options.Ny,options.Nx,1),options.Nx,options.Ny);
clear E2 P2

forcing_k_pair=first.phix.*pair.etax+first.phiy.*pair.etay ...
    +pair.phix.*first.etax+pair.phiy.*first.etay ...
    -first.eta.*pair.phizz-pair.eta.*first.phizz;
forcing_k_direct=first.eta.*(first.phixz.*first.etax ...
    +first.phiyz.*first.etay)-0.5*first.eta.^2.*first.phizzz;
forcing_d_pair=first.eta.*pair.phitz+pair.eta.*first.phitz ...
    +first.phix.*pair.phix+first.phiy.*pair.phiy ...
    +first.phiz.*pair.phiz;
forcing_d_direct=0.5*first.eta.^2.*first.phitzz ...
    +first.eta.*(first.phix.*first.phixz ...
    +first.phiy.*first.phiyz+first.phiz.*first.phizz);
forcing_k=forcing_k_direct/4+forcing_k_pair/2;
forcing_d=forcing_d_direct/4+forcing_d_pair/2;
clear first pair forcing_k_pair forcing_k_direct forcing_d_pair forcing_d_direct
forcing_k_hat=joint_spectrum(forcing_k,options.Nx,options.Ny);
forcing_d_hat=joint_spectrum(forcing_d,options.Nx,options.Ny);
clear forcing_k forcing_d

outbins3_work=(3*min(work_bins)):(3*maximum_work_bin);
outbins3=3*base_bin+outbins3_work;
coefficients=complex(zeros(N,1));outer_total=0;outer_invalid=0;
for output_index=1:numel(outbins3_work)
    output_work_bin=outbins3_work(output_index);
    output_bin=3*base_bin+output_work_bin;
    s=output_bin*frequency_step;
    fk=forcing_k_hat(:,:,output_work_bin+1);
    fd=forcing_d_hat(:,:,output_work_bin+1);
    energy=abs(fk).^2+abs(fd).^2;valid=s>a;
    outer_total=outer_total+sum(energy,'all');
    outer_invalid=outer_invalid+sum(energy(~valid),'all');
    [S,C]=node_response(s,a,valid,lambda3,nodes,weights);
    response=complex(zeros(size(Q)));
    response(valid)=(Aq(valid).*S(valid).*fd(valid) ...
        -1i*C(valid).*fk(valid))/h^2;
    coefficients(output_bin+1)=sum(response,'all');
end
eta33=fft(coefficients);

station_input=squeeze(sum(sum(joint,1),2));
original_input=accumarray(work_bins+1,A,[work_count,1],@sum,0);
station_input_error=norm(station_input-original_input) ...
    /max(norm(original_input),realmin);
audit=struct( ...
    'status','joint-modal-cartesian-nested-green-laplace-eta33-prototype', ...
    'pair_loops',0,'triple_loops',0,'time_snapshot_loops',0, ...
    'inner_quadrature_rank',options.J,'outer_quadrature_rank',options.J, ...
    'inner_scale',lambda2,'outer_scale',lambda3, ...
    'modal_grid',[options.Ny,options.Nx], ...
    'baseband_enabled',logical(options.Baseband),'base_bin',base_bin, ...
    'work_frequency_count',work_count,'input_bin_count',numel(unique(bins)), ...
    'eta33_output_bins',outbins3,'station_input_relative',station_input_error, ...
    'wavevector_projection_rms_relative',projection.q_rms_relative, ...
    'dispersion_residual_rms_relative',projection.dispersion_rms_relative, ...
    'inner_invalid_source_energy_fraction',inner_invalid/max(inner_total,realmin), ...
    'outer_invalid_source_energy_fraction',outer_invalid/max(outer_total,realmin), ...
    'aliasing_control','none','fitted_coefficients',0, ...
    'polynomial_resolvent',false,'stokes_correction',false);
state=struct('coefficients',coefficients,'qx_axis',qx_axis, ...
    'qy_axis',qy_axis,'projection',projection);
end

function [S,C]=node_response(s,a,valid,lambda,nodes,weights)
S=zeros(size(a));C=S;nonzero=valid&a>0;zero=valid&a==0;
av=a(nonzero);
for index=1:numel(nodes)
    tau=nodes(index)/lambda;
    ep=weights(index)/lambda*exp(nodes(index)-(s-av)*tau);
    em=weights(index)/lambda*exp(nodes(index)-(s+av)*tau);
    S(nonzero)=S(nonzero)+(ep-em)./(2*av);
    C(nonzero)=C(nonzero)+(ep+em)/2;
    if any(zero,'all')
        base=weights(index)/lambda*exp(nodes(index)-s*tau);
        S(zero)=S(zero)+base*tau;C(zero)=C(zero)+base;
    end
end
end

function field=joint_field(joint,multiplier,nx,ny)
field=spectral_field(joint.*multiplier,nx,ny);
end

function field=spectral_field(spectrum,nx,ny)
field=fft(ifft(ifft(spectrum,[],1),[],2)*(nx*ny),[],3);
end

function spectrum=joint_spectrum(field,nx,ny)
spectrum=fft2(ifft(field,[],3))/(nx*ny);
end

function [axis,dq]=modal_axis(count,requested_limit,parent,order)
if count==1
    if any(abs(parent)>100*eps(max(1,max(abs(parent)))))
        error('A one-point modal axis can represent only zero components.');
    end
    axis=0;dq=Inf;return
end
if requested_limit==0
    requested_limit=(order+.05)*max(abs(parent))/(1-2/count);
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
    weight=wx.*wy;
    indices{corner}=iy+(ix-1)*ny+bins*ny*nx;
    values{corner}=A.*weight;mass{corner}=abs(A).^2.*weight;
    qxnode=reshape(qx_axis(ix),[],1);qynode=reshape(qy_axis(iy),[],1);
    qnode=hypot(qxnode,qynode);qtarget=hypot(qx,qy);
    qerror{corner}=hypot(qxnode-qx,qynode-qy)./max(qtarget,realmin);
    derror{corner}=abs(nu.^2-qnode.*tanh(qnode))./max(nu.^2,realmin);
end
joint=reshape(accumarray(vertcat(indices{:}),vertcat(values{:}), ...
    [ny*nx*work_count,1],@sum,complex(0)),ny,nx,work_count);
all_mass=vertcat(mass{:});all_qerror=vertcat(qerror{:});
all_derror=vertcat(derror{:});active=all_mass>0;
total_mass=sum(all_mass(active));
report=struct( ...
    'q_rms_relative',sqrt(sum(all_mass(active).*all_qerror(active).^2) ...
        /max(total_mass,realmin)), ...
    'dispersion_rms_relative',sqrt(sum(all_mass(active).*all_derror(active).^2) ...
        /max(total_mass,realmin)));
end

function [i0,i1,w0,w1]=bracket(value,axis,dq)
if numel(axis)==1
    i0=ones(size(value));i1=i0;w0=ones(size(value));w1=zeros(size(value));
    return
end
position=(value-axis(1))/dq+1;i0=floor(position);
i0=max(1,min(numel(axis)-1,i0));i1=i0+1;
w1=max(0,min(1,position-i0));w0=1-w1;
end

function [nodes,weights]=gauss_laguerre_rule(J)
diagonal=2*(1:J)-1;off_diagonal=1:(J-1);
[vectors,values]=eig(diag(diagonal)+diag(off_diagonal,1) ...
    +diag(off_diagonal,-1),'vector');
[nodes,order]=sort(values.');weights=vectors(1,order).^2;
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
if ~isempty(unknown),error('Unknown option: %s.',unknown{1});end
end
