function report=validate_ow3d_campaign_case(caseDir)
% Validate native output structure and extract eta/phi_s remotely, no GL fit.
c=jsondecode(fileread(fullfile(caseDir,'expected.json')));
assert(contains(fileread(fullfile(caseDir,'ow3d.log')),'JOB IS COMPLETE'));
fid=fopen(fullfile(caseDir,'OceanWave3D.end'),'r');assert(fid>=0);guard=onCleanup(@()fclose(fid));
fgetl(fid);ending=sscanf(fgetl(fid),'%f');clear guard;
assert(numel(ending)>=5 && abs(ending(5)-c.expected_final_time_s)<1e-8);
[X,Y,E0,P0]=read_ep(fullfile(caseDir,'EP_00000.bin'));
fid=fopen(fullfile(caseDir,'OceanWave3D.init'),'r');assert(fid>=0);guard=onCleanup(@()fclose(fid));
fgetl(fid);inputHeader=sscanf(fgetl(fid),'%f');pairs=fscanf(fid,'%f',[2,Inf]);clear guard;
assert(isequal(inputHeader(3:4).',[1025,257]) && size(pairs,2)==1025*257);
e=reshape(pairs(1,:),1025,257);p=reshape(pairs(2,:),1025,257);
assert(max(abs(E0(2:1026,2:258)-e),[],'all')<1e-11);
assert(max(abs(P0(2:1026,2:258)-p),[],'all')<1e-10);
files=dir(fullfile(caseDir,'Kinematics*.bin'));assert(isscalar(files));
fid=fopen(fullfile(files(1).folder,files(1).name),'r','ieee-le');assert(fid>=0);guard=onCleanup(@()fclose(fid));
assert(fread(fid,1,'int32')==48);header=fread(fid,9,'int32').';
dt=fread(fid,1,'double');nz=fread(fid,1,'int32');assert(fread(fid,1,'int32')==48);
assert(isequal(header,reshape(c.expected_kinematics_header,1,[])) && nz==18);
assert(abs(dt-c.output_dt_s)<1e-14);
nx=1025;ny=9;n=nx*ny;
grid=reshape(read_record(fid,5*n),5,nx,ny);sigma=read_record(fid,nz);
assert(abs(sigma(end)-1)<1e-14);
ix=2:1026;iy=126:134;x=squeeze(grid(1,:,:));y=squeeze(grid(2,:,:));
assert(max(abs(x-X(ix,iy)),[],'all')<1e-10 && max(abs(y-Y(ix,iy)),[],'all')<1e-10);
steps=(header(7)-1):header(9):(header(8)-1);t=[0,steps*dt].';
assert(numel(t)==c.expected_saved_samples_with_initial && abs(t(end)-220)<1e-8);
eta=zeros(numel(t),nx,ny);psi=eta;
eta(1,:,:)=E0(ix,iy);psi(1,:,:)=P0(ix,iy);checkpoints=zeros(0,4);
for it=1:numel(steps)
    eta(it+1,:,:)=reshape(read_record(fid,n),nx,ny);
    skip_record(fid,n);skip_record(fid,n);
    phi=reshape(read_record(fid,nz*n),nz,nx,ny);psi(it+1,:,:)=squeeze(phi(end,:,:));
    for j=1:12,skip_record(fid,nz*n);end
    if mod(steps(it),275)==0
        [~,~,ec,pc]=read_ep(fullfile(caseDir,sprintf('EP_%05d.bin',steps(it))));
        de=max(abs(squeeze(eta(it+1,:,:))-ec(ix,iy)),[],'all');
        dp=max(abs(squeeze(psi(it+1,:,:))-pc(ix,iy)),[],'all');
        assert(de<1e-11);
        checkpoints(end+1,:)=[t(it+1),de,dp,dp/max(max(abs(pc(ix,iy)),[],'all'),realmin)]; %#ok<AGROW>
    end
end
assert(isempty(fread(fid,1,'uint8')));clear guard;
assert(all(isfinite(eta),'all') && all(isfinite(psi),'all'));
assert(size(checkpoints,1)==4);
out=fullfile(caseDir,'processed');assert(~isfolder(out));mkdir(out);
report=struct('status','RAW_OUTPUT_INTEGRITY_VERIFIED','case_id',c.case_id, ...
    'completed_time_s',ending(5),'samples',numel(t),'kinematics_header',header, ...
    'initial_sample_source','EP_00000.bin; native t=0 volume phi is not used', ...
    'checkpoint_columns',{{'time_s','eta_max_abs_difference_m','psi_max_abs_difference_m2s','psi_relative_Linf'}}, ...
    'kinematics_vs_EP_checkpoints',checkpoints, ...
    'psi_checkpoint_relative_tolerance_warning',any(checkpoints(:,4)>1e-4), ...
    'boundary',c.boundary,'physical_accuracy_certified',false,'GL_comparison_performed',false, ...
    'max_abs_eta_m',max(abs(eta),[],'all'),'max_abs_psi_m2s',max(abs(psi),[],'all'));
save(fullfile(out,'surface_strip.mat'),'t','x','y','eta','psi','report','-v7.3');
fid=fopen(fullfile(out,'validation.json'),'w');assert(fid>=0);guard=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(report));disp(report);
end

function a=read_record(fid,n)
assert(fread(fid,1,'int32')==8*n);a=fread(fid,n,'double');
assert(numel(a)==n && fread(fid,1,'int32')==8*n && all(isfinite(a)));
end

function skip_record(fid,n)
assert(fread(fid,1,'int32')==8*n);assert(fseek(fid,8*n,'cof')==0);
assert(fread(fid,1,'int32')==8*n);
end

function [x,y,e,p]=read_ep(name)
fid=fopen(name,'r','ieee-le');assert(fid>=0);guard=onCleanup(@()fclose(fid));
assert(fread(fid,1,'int32')==8);sz=fread(fid,2,'int32').';assert(fread(fid,1,'int32')==8);
n=prod(sz);xy=read_record(fid,2*n);ep=read_record(fid,2*n);
x=reshape(xy(1:n),sz);y=reshape(xy(n+1:end),sz);e=reshape(ep(1:n),sz);p=reshape(ep(n+1:end),sz);
end
