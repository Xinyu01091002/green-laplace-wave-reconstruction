function prepare_jonswap_rescaled(root,source,target,durationTp)
% Exact order-1/order-2 amplitude homogeneity of the frozen R4-GL input.
% The source fields are generated initial conditions, not evolved references.
if nargin<4,durationTp=80;end
assert(ismember(durationTp,[20,80]));
run=fullfile(root,'medium');assert(~isfolder(run));
s=jsondecode(fileread(fullfile(source,'high','settings.json')));
d=load(fullfile(source,'high','inputs','initial_fields.mat'));
low=load(fullfile(source,'low','inputs','initial_fields.mat'),'E','P','C');
assert(s.kpHs_over_2==.12&&target==.06);
assert(contains(s.initialization,'no initial31/33'));
nx=s.grid(1);ny=s.grid(2);Lx=s.domain_m(1);Ly=s.domain_m(2);
kx=d.kx;ky=d.ky;om=d.om;a=target/s.kpHs_over_2;
mx=kx*Lx/(2*pi);my=ky*Ly/(2*pi);
assert(max(abs(mx-round(mx)))<1e-9&&max(abs(my-round(my)))<1e-9);
ids=sub2ind([ny,nx],mod(round(my),ny)+1,mod(round(mx),nx)+1);
assert(numel(unique(ids))==numel(ids));
F=complex(zeros(ny,nx));F(ids)=nx*ny*d.C;
eta1=ifft2(F);F(ids)=-1i*s.g./om(:).*F(ids);psi1=ifft2(F);
E=zeros(size(d.E));P=zeros(size(d.P));E2=E;P2=P;
lowErrors=zeros(4,2);firstErrors=zeros(4,2);
for j=1:4
    phase=(j-1)*pi/2;e1=real(eta1*exp(1i*phase));p1=real(psi1*exp(1i*phase));
    opposite=mod(j+1,4)+1;
    firstErrors(j,:)=[rel((d.E(:,:,j)-d.E(:,:,opposite))/2,e1), ...
        rel((d.P(:,:,j)-d.P(:,:,opposite))/2,p1)];
    E2(:,:,j)=d.E(:,:,j)-e1;P2(:,:,j)=d.P(:,:,j)-p1;
    lowErrors(j,:)=[rel(e1/6+E2(:,:,j)/36,low.E(:,:,j)), ...
        rel(p1/6+P2(:,:,j)/36,low.P(:,:,j))];
    E(:,:,j)=a*e1+a^2*E2(:,:,j);P(:,:,j)=a*p1+a^2*P2(:,:,j);
end
assert(max(firstErrors,[],'all')<1e-11&&max(lowErrors,[],'all')<1e-11);
assert(rel(d.C/6,low.C)<1e-12);
assert(all(isfinite(E),'all')&&all(isfinite(P),'all'));
assert(max(abs(mean(E,[1,2])),[],'all')<1e-9);
C=a*d.C;Hs=4*std(real(a*eta1(:)),1);
assert(abs(Hs-2*target/s.kp)/(2*target/s.kp)<1e-11);
audit=d.audit;
audit.amplitude_scaling='kpHs/2 .12 to .06: identical parents/phases, linear x1/2 and quadratic x1/4; no total-field scaling';
audit.source_campaign=source;audit.source_field_role='Frozen declared-input R4-GL initial conditions, not HOS output';
audit.first_order_from_declared_C_relative_errors=firstErrors;
audit.existing_low_input_recovery_relative_errors=lowErrors;
audit.Hs_linear_measured_m=Hs;
audit.source_initialization_audit_file=fullfile(source,'initialization-audit.json');
KX=repmat(2*pi/Lx*[0:nx/2-1,-nx/2:-1],ny,1);stats=cell(4,1);
for j=1:4
    e=E(:,:,j);ef=real(interpft(e,4*nx,2));
    ex=real(ifft2(1i*KX.*fft2(e)));ex=real(interpft(interpft(ex,2*ny,1),2*nx,2));
    stats{j}=struct('phase_deg',(j-1)*90,'max_eta_native_m',max(e,[],'all'), ...
        'min_eta_native_m',min(e,[],'all'),'max_abs_x_slope_2x',max(abs(ex),[],'all'), ...
        'max_x_upcrossing_height_native_m',wave_height(e), ...
        'max_x_upcrossing_height_x4_m',wave_height(ef), ...
        'eta2_over_eta1_L2',norm(a^2*E2(:,:,j),'fro')/norm(real(a*eta1*exp(1i*(j-1)*pi/2)),'fro'));
end
audit.initial_statistics=[stats{:}];
mkdir(run);mkdir(fullfile(run,'inputs'));
save(fullfile(run,'inputs','initial_fields.mat'),'C','kx','ky','om','E','P','audit','-v7.3');
s.kpHs_over_2=target;s.Hs_linear=Hs;s.source_campaign=source;s.amplitude_scaling=audit.amplitude_scaling;
s.requested_duration_Tp=durationTp;
s.duration_s=.4*round(durationTp*s.Tp/.4);
s.duration_Tp=s.duration_s/s.Tp;
s.expected_samples=round(s.duration_s/s.output_dt_s)+1;
if durationTp==20,s.scoring_window_s=[3,17]*s.Tp;end
assert(mod(s.expected_samples,2)==1);
for j=1:4
    for p=1:size(s.probe_indices,1)
        s.expected_initial_probe_eta(j,p)=E(s.probe_indices(p,2),s.probe_indices(p,1),j);
    end
end
writejson(fullfile(run,'settings.json'),s);export_jonswap_hos(run);
exportError=zeros(1,2);
for j=1:4
    for rank=0:s.MPI_ranks_per_phase-1
        file=fullfile(run,'cases',sprintf('phi%03d',(j-1)*90),'Results',sprintf('3d_ini_%03d.dat',rank));
        fid=fopen(file,'r');assert(fid>=0);guard=onCleanup(@()fclose(fid));line=fgetl(fid);
        while ischar(line)&&~startsWith(strtrim(line),'ZONE'),line=fgetl(fid);end
        assert(ischar(line));v=fscanf(fid,'%f',[2,Inf]);clear guard;
        rows=rank*(ny/s.MPI_ranks_per_phase)+(1:ny/s.MPI_ranks_per_phase);
        ee=E(rows,:,j).';pp=P(rows,:,j).';assert(size(v,2)==numel(ee));
        exportError=max(exportError,[max(abs(v(1,:).'-ee(:))),max(abs(v(2,:).'-pp(:)))]);
    end
end
assert(all(exportError==0));audit.export_eta_psi_max_abs_error=exportError;
writejson(fullfile(root,'initialization-audit.json'),audit);
disp(struct2table([stats{:}]));fprintf('Hs=%.12g m; maximum input/export discrepancy=%.3g\n',Hs,max(exportError));
end

function v=rel(a,b)
v=norm(a(:)-b(:))/norm(b(:));
end

function height=wave_height(e)
% Complete periodic +x spatial waves; no temporal-height claim or demeaning.
height=0;n=size(e,2);
for row=1:size(e,1)
    v=e(row,:);cross=find(v<=0 & circshift(v,-1)>0);
    if isempty(cross),continue;end
    edges=[cross,cross(1)+n];periodic=[v,v];
    for j=1:numel(edges)-1
        segment=periodic(edges(j)+1:edges(j+1));height=max(height,max(segment)-min(segment));
    end
end
end

function writejson(file,value)
fid=fopen(file,'w');assert(fid>=0);guard=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(value,PrettyPrint=true));
end
