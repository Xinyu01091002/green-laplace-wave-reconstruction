function report=prepare_compact_wavegroup(designFile,mf12Root,out,phaseRandomizationSource,randomSeed)
% Independent MF12 order-2 initial fields; no GL formula or OW3D run here.
% Read an immutable first-order design snapshot and write a new remote case set.
% Optional fourth argument: preserve that family's modal amplitudes and
% randomize phases only. Export periodic fields, not closed-wall OW3D inputs.
randomized=nargin>=4;
if nargin<5,randomSeed=20260925;end
validateattributes(randomSeed,{'numeric'},{'scalar','integer','nonnegative','<',2^32});
if isfolder(out),error('Case directory already exists.');end
addpath(mf12Root);mkdir(out);d=load(designFile);r=d.report;
g=r.g;h=r.h;kp=r.kp;nx=r.unique_FFT_nodes(1);ny=r.unique_FFT_nodes(2);
Lx=r.domain_m(1);Ly=r.domain_m(2);A=r.group_linear_focus_amplitude_m;
tf=110;duration=220;dt=.2;threshold=1e-4;
active=d.weights>=max(d.weights,[],'all')*threshold;
omitted=sum(d.weights(~active))/sum(d.weights,'all');
w=d.weights(active);w=w/sum(w);kx=d.kx(active).';ky=d.ky(active).';om=d.omega(active).';
C=A*w.'.*exp(-1i*(kx*Lx/2+ky*Ly/2)+1i*om*tf);
family='wavegroup';randomization=struct();
if randomized
    parent=load(phaseRandomizationSource,'C','kx','ky','om');
    assert(isequal(kx,parent.kx) && isequal(ky,parent.ky) && isequal(om,parent.om));
    rng(randomSeed,'twister');phaseIncrements=2*pi*rand(size(parent.C));
    C=parent.C.*exp(1i*phaseIncrements);family='randomphase';tf=NaN;
    relativeAmplitudeChange=max(abs(abs(C)-abs(parent.C)))/max(abs(parent.C));
    assert(relativeAmplitudeChange<1e-14);
    sigma=sqrt(sum(abs(C).^2)/2);
    randomization=struct('rule','C_random=C_group*exp(i*independent_uniform_phase); no amplitude change', ...
        'seed',randomSeed,'phase_increments',phaseIncrements,'source',phaseRandomizationSource, ...
        'maximum_relative_amplitude_change',relativeAmplitudeChange, ...
        'linear_spatial_rms_m',sigma,'linear_Hs_4sigma_m',4*sigma, ...
        'kp_Hs_over_2',kp*2*sigma,'focusing_Akp_label',kp*sum(abs(C)), ...
        'physical_space_taper',false,'amplitude_renormalization',false);
    fprintf('Random phases only: linear Hs=%.9g m, kp*Hs/2=%.9g.\n',4*sigma,kp*2*sigma);
end
fprintf('Preparing %d declared first-order parents; omitted modal L1 fraction %.6g.\n',numel(C),omitted);
% Verify the external preallocated spectral path includes both order-2 sectors.
% Its superharmonic_only flag controls third-order retention, not these pairs.
verify_external_order2(g,h,kp,Lx,Ly,nx,ny);
timer=tic;
c=mf12_spectral_coefficients(2,g,h,real(C),imag(C),kx,ky,0,0);
coefficientSeconds=toc(timer);
assert(~isfield(c,'G_3') && all(c.muStar==0));
assert(numel(c.G_npm)==numel(C)*(numel(C)-1));
assert(all(isfinite(c.G_npm)) && all(isfinite(c.mu_npm)));
assert(all(hypot(c.kx_npm,c.ky_npm)>0),'Duplicate parent or strict-zero pair.');
assert(all(abs(c.kx_npm)<pi*nx/Lx) && all(abs(c.ky_npm)<pi*ny/Ly));
assert(all(abs(2*kx)<pi*nx/Lx) && all(abs(2*ky)<pi*ny/Ly));
fprintf('MF12 order-2 coefficients prepared in %.3f s.\n',coefficientSeconds);
% Deposit declared first-order coefficients exactly on their native grid.
lin=complex(zeros(ny,nx));ind=find(active);lin(ind)=nx*ny*C;
linP=complex(zeros(ny,nx));linP(ind)=nx*ny*(-1i*g./om.*C);
eta1=ifft2(lin);psi1=ifft2(linP);
phases=[0,90,180,270];E=zeros(ny,nx,4);P=E;metrics=zeros(4,7);
edge=false(ny,nx);edge([1,end],:)=true;edge(:,[1,end])=true;
plus=1:2:numel(c.G_npm);self=complex(c.A_2,c.B_2);pairs=complex(c.A_npm,c.B_npm);
for ip=1:4
    ph=deg2rad(phases(ip));cc=c;
    z=C*exp(1i*ph);cc.a=real(z);cc.b=imag(z);
    z=self*exp(2i*ph);cc.A_2=real(z);cc.B_2=imag(z);
    z=pairs;z(plus)=z(plus)*exp(2i*ph);cc.A_npm=real(z);cc.B_npm=imag(z);
    [eta,psi]=mf12_spectral_surface(cc,Lx,Ly,nx,ny,0);
    assert(all(isfinite(eta),'all') && all(isfinite(psi),'all'));
    assert(abs(mean(eta,'all'))<1e-10 && abs(mean(psi,'all'))<1e-9);
    E(:,:,ip)=eta;P(:,:,ip)=psi;
    name=sprintf('%s_kpd1_akp012_phi%03d',family,phases(ip));folder=fullfile(out,name);mkdir(folder);
    if randomized
        % Same unique periodic grid for future HOS import. No taper changes
        % the user-prescribed amplitudes; closed-wall OW3D is not configured.
        save(fullfile(folder,'initial_surface.mat'),'eta','psi','Lx','Ly','g','h','-v7.3');
    else
        write_initial(fullfile(folder,'OceanWave3D.init'),eta,psi,Lx,Ly,0);
        write_input(fullfile(folder,'OceanWave3D.inp'),Lx,Ly,h,nx+1,ny+1,17,g,dt,duration);
    end
    metrics(ip,:)=[phases(ip),max(abs(eta),[],'all'),max(abs(psi),[],'all'), ...
        max(abs(eta(edge)))/A,max(abs(psi(edge)))/(g/(2*pi/r.Tp_s)*A), ...
        min(eta,[],'all'),max(eta,[],'all')];
    fprintf('Wrote %s; eta max %.6g m.\n',name,metrics(ip,2));
end
etaFirst=(E(:,:,1)-1i*E(:,:,2)-E(:,:,3)+1i*E(:,:,4))/2;
psiFirst=(P(:,:,1)-1i*P(:,:,2)-P(:,:,3)+1i*P(:,:,4))/2;
assert(norm(etaFirst-eta1,'fro')/norm(eta1,'fro')<1e-10);
assert(norm(psiFirst-psi1,'fro')/norm(psi1,'fro')<1e-10);
eta20=mean(E,3);psi20=mean(P,3);
eta22=(E(:,:,1)-E(:,:,2)+E(:,:,3)-E(:,:,4))/4;
psi22=(P(:,:,1)-P(:,:,2)+P(:,:,3)-P(:,:,4))/4;
report=struct('status','INITIAL_FIELDS_PREPARED_NOT_PROPAGATION_VALIDATED','design',r, ...
    'focus_time_s',tf,'duration_s',duration,'integration_dt_s',dt,'output_dt_s',dt, ...
    'initialization_order',2,'initialization_source','external MF12 order 2, complete nonzero sum and difference pairs', ...
    'strict_zero_mode','eta and psi mean zero; no imposed return current', ...
    'parent_count',numel(C),'sparsification_relative_amplitude',threshold,'omitted_modal_L1_fraction',omitted, ...
    'coefficient_seconds',coefficientSeconds,'phase_metrics',metrics, ...
    'phase_metric_columns',{{'phase_degrees','peak_abs_eta_m','peak_abs_psi_m2s','edge_eta_over_focus_A', ...
    'edge_psi_over_gAwp','min_eta_m','max_eta_m'}}, ...
    'no_GL_high_order_initialization',true,'third_order_MF12_computed',false, ...
    'boundary','Native straight-wall boundary; packet edge audit required. No periodic OW3D boundary is claimed.', ...
    'initial_output_source','Use EP_00000.bin for eta/psi at t=0; kinematics begins at dt', ...
    'kinematics_physical_node_range',[1,nx+1,1,ny/2-3,ny/2+5,1,2,round(duration/dt)+1,1]);
if randomized
    report.status='PERIODIC_RANDOM_INITIAL_FIELDS_PREPARED_SOLVER_IMPORT_PENDING';
    report.randomization=randomization;
    report.boundary='Periodic realization; no taper or closed-wall OW3D run configuration';
    report.initial_output_source='initial_surface.mat eta and true surface psi on the unique periodic grid';
    report.kinematics_physical_node_range=[];
    report.design=rmfield(report.design,{'random_prototype','nonlinear_initialization','status'});
    measuredSigma=sqrt(mean(real(eta1).^2,'all'));
    assert(abs(measuredSigma/randomization.linear_spatial_rms_m-1)<1e-12);
end
save(fullfile(out,'initial_fields.mat'),'report','C','kx','ky','om','E','P','eta1','psi1','eta20','psi20','eta22','psi22','-v7.3');
fid=fopen(fullfile(out,'initialization.json'),'w');assert(fid>=0);guard=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(report));disp(metrics);disp(report.status);
end

function verify_external_order2(g,h,kp,Lx,Ly,nx,ny)
kx=(2*pi/Lx)*[42,50,54];ky=(2*pi/Ly)*[-5,0,5];
a=[.03,-.02,.01]/kp;b=[.01,.005,-.02]/kp;
c=mf12_spectral_coefficients(2,g,h,a,b,kx,ky,0,0);
ref=mf12_direct_coefficients(2,g,h,a,b,kx,ky,0,0);
for f={'G_2','mu_2','G_npm','mu_npm'}
    z=c.(f{1});v=ref.(f{1});assert(norm(z-v)/max(norm(v),realmin)<1e-11);
end
[e,p,X,Y]=mf12_spectral_surface(c,Lx,Ly,nx,ny,0);
for ind=[1,129*ny+9,511*ny+128]
    [er,pr]=mf12_direct_surface(2,ref,X(ind),Y(ind),0);
    assert(abs(e(ind)-er)<1e-10 && abs(p(ind)-pr)<1e-9);
end
fprintf('External MF12 order-2 sum/difference and surface-potential parity passed.\n');
end

function write_initial(name,eta,psi,Lx,Ly,t0)
% OW3D writes x fastest, then y; repeat the periodic construction endpoints.
e=eta([1:end,1],[1:end,1]).';p=psi([1:end,1],[1:end,1]).';
fid=fopen(name,'w');assert(fid>=0);guard=onCleanup(@()fclose(fid));
fprintf(fid,'Independent MF12 second-order initialization from declared first-order parents\n');
fprintf(fid,'%.17g %.17g %d %d %.17g\n',Lx,Ly,size(e,1),size(e,2),t0);
fprintf(fid,'%.17g %.17g\n',[e(:),p(:)].');
end

function write_input(name,Lx,Ly,h,nx,ny,nz,g,dt,duration)
steps=round(duration/dt)+1;saveEP=round(55/dt);
fid=fopen(name,'w');assert(fid>=0);guard=onCleanup(@()fclose(fid));
fprintf(fid,'Compact directional wavegroup kpd=1 Akp=0.12; independent order-2 initialization\n');
% Supply the optional third value: list-directed Fortran reads otherwise
% continue onto the grid line when there is no trailing nonnumeric comment.
fprintf(fid,'-1 0 1000\n');
fprintf(fid,'%.17g %.17g %.17g %d %d %d 0 0 1 1 1 1\n',Lx,Ly,h,nx,ny,nz);
fprintf(fid,'4 4 4 1 1 1\n%d %.17g 1 0.0 1 0.0\n%.17g 1000\n',steps,dt,g);
fprintf(fid,'1 3 0 55 1e-6 1e-6 1 V 1 1 20\n');
fprintf(fid,'0.05 1.00 1.84 2 0 0 1 6 32\n');
fprintf(fid,'%d 20 1 1\n',saveEP);
fprintf(fid,'1 %d 1 %d %d 1 2 %d 1\n',nx,(ny+1)/2-4,(ny+1)/2+4,steps);
fprintf(fid,'1 0\n0 6 10 0.08 0.08 0.4\n0 0 0 0 0 0 0\n0 0.0 0 X 0.0\n0 0\n0 2.0 2 0 0 1 0\n0\n');
end
