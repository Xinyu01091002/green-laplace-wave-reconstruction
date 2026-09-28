function [eta,audit]=gl_native_joint_eta22(initial,first,t,s,probe,J)
%GL_NATIVE_JOINT_ETA22 Native-k-grid, observed-bin-conditioned GL eta22.
if nargin<6,J=8;end
N=numel(t);dt=mean(diff(t));df=2*pi/(N*dt);tr=t-t(1);
C=initial.C(:);kx=initial.kx(:);ky=initial.ky(:);om=initial.om(:);
use=abs(C)>0&kx>0;C=C(use);kx=kx(use);ky=ky(use);om=om(use);
nx=s.grid(1);ny=s.grid(2);Lx=s.domain_m(1);Ly=s.domain_m(2);
lo=ceil(min(om)/df);hi=floor(max(om)/df);bins=(lo:hi)';K=numel(bins);W=2*K-1;
F=fft(first)/N;observed=2*conj(F(bins+1));

mx=round(kx*Lx/(2*pi));my=round(ky*Ly/(2*pi));
grid_kx=mx*2*pi/Lx;grid_ky=my*2*pi/Ly;
k_error=max(hypot(grid_kx-kx,grid_ky-ky));
ix=mod(mx,nx)+1;iy=mod(my,ny)+1;
position=om/df;lower=floor(position);upper=lower+1;upper_weight=position-lower;lower_weight=1-upper_weight;
xy=s.probes_xy_m(probe,:);phase=C.*exp(1i*(kx*xy(1)+ky*xy(2)));
indices=[];values=[];mass=[];
for side=1:2
    if side==1,m=lower;weight=lower_weight;else,m=upper;weight=upper_weight;end
    active=m>=lo&m<=hi&weight>0;local=m(active)-lo;
    indices=[indices;iy(active)+(ix(active)-1)*ny+local*ny*nx]; %#ok<AGROW>
    values=[values;phase(active).*weight(active)]; %#ok<AGROW>
    mass=[mass;abs(phase(active)).*weight(active)]; %#ok<AGROW>
end
joint=reshape(accumarray(indices,values,[ny*nx*W,1],@sum,complex(0)),ny,nx,W);
prior=squeeze(sum(sum(joint(:,:,1:K),1),2));
condition=zeros(K,1);
for local=1:K
    active=indices>(local-1)*ny*nx&indices<=local*ny*nx;
    condition(local)=sum(mass(active))/max(abs(prior(local)),realmin);
end
for local=1:K
    slice=joint(:,:,local);weights=abs(slice).^2;weight_sum=sum(weights,'all');
    if weight_sum==0,error('A native temporal bin has zero prior support.');end
    slice=slice+(observed(local)-prior(local))*weights/weight_sum;
    joint(:,:,local)=slice;
end
recovered=real(exp(-1i*tr*(df*bins).')*observed);

mode_x=[0:nx/2-1,-nx/2:-1];mode_y=[0:ny/2-1,-ny/2:-1];
[KX,KY]=meshgrid(mode_x*2*pi/Lx,mode_y*2*pi/Ly);QX=hreshape(s.h*KX);QY=hreshape(s.h*KY);Q2=QX.^2+QY.^2;
nu=(df*(lo+(0:W-1)))*sqrt(s.h/s.g);NU=reshape(nu,1,1,[]);safe=NU;safe(safe==0)=1;
v=field(1);vn2=field(NU.^2);sd=2*v.*vn2;clear vn2
vn=field(NU);sd=sd+vn.^2;clear vn
hx=field(QX./safe);hy=field(QY./safe);sd=sd-hx.^2-hy.^2;
jx=field(QX);jy=field(QY);radial=field(Q2./safe);sk=2*v.*radial+2*(hx.*jx+hy.*jy);
clear v hx hy jx jy radial
ds=fft2(ifft(sd,[],3))/(nx*ny);ks=fft2(ifft(sk,[],3))/(nx*ny);clear sd sk

Q=hypot(QX(:,:,1),QY(:,:,1));a=sqrt(Q.*tanh(Q));[nodes,weights]=laguerre(J);
qp=s.kp*s.h;lambda=sqrt(4*qp*tanh(qp)-2*qp*tanh(2*qp));coeff=complex(zeros(N,1));
invalid=0;total=0;
for local=0:2*(K-1)
    output_bin=2*lo+local;sout=output_bin*df*sqrt(s.h/s.g);
    d=ds(:,:,local+1);k=ks(:,:,local+1);energy=abs(d).^2+abs(k).^2;valid=sout>a;
    total=total+sum(energy,'all');invalid=invalid+sum(energy(~valid),'all');md=zeros(size(Q));mk=md;av=a(valid);
    for j=1:J,tau=nodes(j)/lambda;ep=weights(j)/lambda*exp(nodes(j)-(sout-av)*tau);em=weights(j)/lambda*exp(nodes(j)-(sout+av)*tau);md(valid)=md(valid)-av.*(ep-em)/2;mk(valid)=mk(valid)+(ep+em)/2;end
    response=complex(zeros(size(Q)));response(valid)=(md(valid).*d(valid)+mk(valid).*k(valid))/(4*s.h);coeff(output_bin+1)=sum(response,'all');
end
eta=fft(coeff);
audit=struct('status','native-k-grid-minimum-energy-update-eta22', ...
    'native_grid',[ny,nx],'input_modes',numel(C),'input_bins',bins.', ...
    'work_frequency_count',W,'first_projection_relative',norm(recovered-first)/norm(first), ...
    'retained_first_energy',sum(2*abs(F(bins+1)).^2)/sum(abs(F).^2), ...
    'maximum_temporal_bin_condition',max(condition),'energy_weighted_condition', ...
    sum(condition.*abs(observed).^2)/sum(abs(observed).^2), ...
    'maximum_native_k_rounding_error',k_error,'invalid_source_energy_fraction',invalid/max(total,realmin), ...
    'pair_loops',0,'spatial_wavevector_projection',false,'temporal_frequency_linear_deposit',true, ...
    'observed_bin_update','minimum weighted-norm additive correction; exact station-bin sum');
    function u=field(multiplier)
        u=fft(ifft(ifft(joint.*multiplier,[],1),[],2)*(nx*ny),[],3);
    end
end
function value=hreshape(value),value=reshape(value,size(value,1),size(value,2),1);end
function [nodes,weights]=laguerre(J),[V,D]=eig(diag(2*(1:J)-1)+diag(1:J-1,1)+diag(1:J-1,-1),'vector');[nodes,index]=sort(D.');weights=V(1,index).^2;end
