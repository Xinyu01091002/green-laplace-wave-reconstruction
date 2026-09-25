function report=design_compact_ow3d(out)
% First-order geometry/spectrum design only; no nonlinear solver execution.
% Execute remotely. The random prototype has unit eta RMS until normalized.
if nargin<1,error('Specify an explicit remote output directory.');end
if isfolder(out),error('Output directory already exists; choose a new design ID.');end
mkdir(out);
g=9.81;kp=.0279;h=1/kp;Akp=.12;A=Akp/kp;lambda=2*pi/kp;
nx=1024;ny=256;nz=17;Lx=50*lambda;Ly=20*lambda;
wp=sqrt(g*kp*tanh(kp*h));Tp=2*pi/wp;
[mx,my]=meshgrid([0:nx/2-1,-nx/2:-1],[0:ny/2-1,-ny/2:-1]);
kx=mx*(2*pi/Lx);ky=my*(2*pi/Ly);k=hypot(kx,ky);theta=atan2(ky,kx);
omega=sqrt(g*k.*tanh(h*k));kw=.004606;alpha=8;
width=kw*ones(size(k));width(k>=kp)=kp/sqrt(2*log(10^alpha));
polarWeight=exp(-.5*((k-kp)./width).^2).*exp(-.5*(theta/deg2rad(25)).^2);
% Cartesian quadrature of a polar modal-amplitude density dk dtheta.
weights=zeros(size(k));forward=kx>0;
weights(forward)=polarWeight(forward)./k(forward)*(2*pi/Lx)*(2*pi/Ly);
% Prescribed numerical sparsification, chosen before any OW3D/GL comparison.
thresholds=[1e-4,1e-5,1e-6,1e-7,1e-8];counts=zeros(size(thresholds));lost=counts;lostEnergy=counts;
for j=1:numel(thresholds)
    active=weights>=max(weights,[],'all')*thresholds(j);
    counts(j)=nnz(active);lost(j)=sum(weights(~active))/sum(weights,'all');
    lostEnergy(j)=sum(weights(~active).^2)/sum(weights.^2,'all');
end
tailTable=table(thresholds(:),counts(:),lost(:),lostEnergy(:), ...
    'VariableNames',{'relative_amplitude_threshold','parents','omitted_L1_fraction','omitted_squared_amplitude_fraction'});
% Retain the full sampled first-order spectrum in this geometry audit.
weights=weights/sum(weights,'all');xf=Lx/2;yf=Ly/2;
leads=[60,90,110,130];records=cell(numel(leads),1);
edge=false(ny,nx);edge([1,end],:)=true;edge(:,[1,end])=true;
[X,Y]=meshgrid((0:nx-1)*Lx/nx,(0:ny-1)*Ly/ny);
nearEdge=X<2*lambda|X>Lx-2*lambda|Y<2*lambda|Y>Ly-2*lambda;
for j=1:numel(leads)
    C=A*weights.*exp(-1i*(kx*xf+ky*yf)+1i*omega*leads(j));
    metrics=zeros(0,7);
    for t=0:10:2*leads(j)
        z=C.*exp(-1i*omega*t);eta=ifft2(nx*ny*z);
        p=zeros(size(z));p(forward)=-1i*g./omega(forward).*z(forward);psi=ifft2(nx*ny*p);
        metrics(end+1,:)=[t,max(abs(eta(edge)))/A,max(abs(psi(edge)))/(g/wp*A), ...
            sum(abs(eta(nearEdge)).^2)/sum(abs(eta).^2,'all'), ...
            max(abs(eta),[],'all'),max(abs(real(eta)),[],'all'),max(abs(real(psi)),[],'all')]; %#ok<AGROW>
    end
    records{j}=struct('focus_time_s',leads(j),'metrics',metrics, ...
        'metric_columns',{{'time_s','edge_analytic_eta_over_A','edge_analytic_psi_over_gAwp', ...
        'two_lambda_edge_band_squared_eta_fraction','peak_analytic_eta_m','peak_real_eta_m','peak_real_psi_m2s'}}, ...
        'max_edge_eta_over_A',max(metrics(:,2)),'max_edge_psi_over_gAwp',max(metrics(:,3)));
end
rng(20260925,'twister');randomPhase=2*pi*rand(size(weights));
randomUnit=sqrt(2)*weights/sqrt(sum(weights.^2,'all')).*exp(1i*randomPhase);
eRandom=real(ifft2(nx*ny*randomUnit));
randomInfo=struct('seed',20260925,'normalization','unit spatial eta RMS; physical steepness pending', ...
    'measured_rms',sqrt(mean(eRandom.^2,'all')), ...
    'edge_rms',sqrt(mean(eRandom(edge).^2)), ...
    'boundary','Nonzero at tank walls. Not a ready OW3D initial state; boundary design required.');
energy=weights.^2;meanW=sum(omega.*energy,'all')/sum(energy,'all');
beta=sqrt(sum((omega-meanW).^2.*energy,'all')/sum(energy,'all'))/meanW;
report=struct('status','FIRST_ORDER_DESIGN_ONLY_NOT_RUN_READY','g',g,'kp',kp,'h',h, ...
    'kpd',1,'Akp_group',Akp,'group_linear_focus_amplitude_m',A,'lambda_p_m',lambda, ...
    'Tp_s',Tp,'OW3D_nodes',[nx+1,ny+1,nz],'unique_FFT_nodes',[nx,ny], ...
    'domain_m',[Lx,Ly],'focus_xy_m',[xf,yf],'dx_dy_m',[Lx/nx,Ly/ny], ...
    'output_dt_s',.2,'samples_at_3fp',Tp/(3*.2),'radial_amplitude_width_below_kp',kw, ...
    'radial_high_tail_alpha',alpha,'angular_amplitude_Gaussian_width_deg',25, ...
    'angle_support','strict kx>0; no additional angular truncation', ...
    'frequency_energy_std_over_mean',beta,'leads',[records{:}],'sparsification_candidates',table2struct(tailTable), ...
    'random_prototype',randomInfo,'phase_convention','Re(C exp(i k dot x - i omega t)); psi1 coefficient=-i g C/omega', ...
    'nonlinear_initialization','pending independent MF12 order 2, eta11+eta20+eta22 and true surface potential', ...
    'no_Stokes_repair',true,'reference_or_GL_errors_used_to_select_input',false);
save(fullfile(out,'first_order_design.mat'),'report','weights','kx','ky','omega','randomUnit','X','Y','-v7.3');
fid=fopen(fullfile(out,'design.json'),'w');assert(fid>=0);c=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(report));writetable(tailTable,fullfile(out,'spectrum_size.csv'));
disp(tailTable);disp(struct2table(rmfield([records{:}],{'metrics','metric_columns'})));disp(report);
end
