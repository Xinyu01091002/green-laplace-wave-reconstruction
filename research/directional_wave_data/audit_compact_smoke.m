function report=audit_compact_smoke(smoke,initialFile,designFile,initialOnly)
% Native binary IO parity and a surface-gradient diagnostic, not Euler certification.
if nargin<4,initialOnly=false;end
resource=jsondecode(fileread(fullfile(smoke,'resources.json')));
frames=3;finalTime=NaN;status='TWO_STEP_IO_SMOKE_PASS';
if initialOnly
    frames=1;status='INITIAL_OUTPUT_ONLY_NOT_A_COMPLETED_TIME_STEP';
else
    assert(resource.exit_code==0 && ~resource.timed_out);
    assert(contains(fileread(fullfile(smoke,'ow3d.log')),'JOB IS COMPLETE'));
    fid=fopen(fullfile(smoke,'OceanWave3D.end'),'r');assert(fid>=0);guard=onCleanup(@()fclose(fid));
    fgetl(fid);ending=sscanf(fgetl(fid),'%f');clear guard;finalTime=ending(5);
    assert(abs(finalTime-.4)<1e-10);
end
files=dir(fullfile(smoke,'Kinematics*.bin'));assert(isscalar(files));
fid=fopen(fullfile(files(1).folder,files(1).name),'r','ieee-le');assert(fid>=0);guard=onCleanup(@()fclose(fid));
assert(fread(fid,1,'int32')==48);header=fread(fid,9,'int32').';
dt=fread(fid,1,'double');nz=fread(fid,1,'int32');assert(fread(fid,1,'int32')==48);
nx=(header(2)-header(1))/header(3)+1;ny=(header(5)-header(4))/header(6)+1;
assert(isequal(header,[1,1025,1,125,133,1,1,3,1]) && nz==18 && dt==.2);
grid=read_record(fid,5*nx*ny);sigma=read_record(fid,nz);assert(abs(sigma(end)-1)<1e-14);
grid=reshape(grid,5,nx,ny);errors=zeros(frames,2);
for it=1:frames
    eta=reshape(read_record(fid,nx*ny),nx,ny);
    read_record(fid,nx*ny);read_record(fid,nx*ny);
    volume=reshape(read_record(fid,nz*nx*ny),nz,nx,ny);psi=squeeze(volume(end,:,:));
    for j=1:12,read_record(fid,nz*nx*ny);end
    [X,Y,E,P]=read_ep(fullfile(smoke,sprintf('EP_%05d.bin',it-1)));
    ix=2:1026;iy=126:134;
    xx=squeeze(grid(1,:,:));yy=squeeze(grid(2,:,:));
    assert(max(abs(xx-X(ix,iy)),[],'all')<1e-10 && max(abs(yy-Y(ix,iy)),[],'all')<1e-10);
    errors(it,:)=[max(abs(eta-E(ix,iy)),[],'all'),max(abs(psi-P(ix,iy)),[],'all')];
end
if ~initialOnly,assert(isempty(fread(fid,1,'uint8')));end
clear guard;
disp(errors);
assert(all(isfinite(errors),'all') && max(errors(:,1))<1e-11);
% This native build writes its unsolved volume phi at t=0. Use initial EP P
% for that time; only compare kinematics phi after completed advances.
if ~initialOnly,assert(max(errors(2:end,2))<1e-8);end
d=load(designFile,'kx','ky','report');s=load(initialFile,'P');r=d.report;
velocityScale=r.g*r.kp*r.group_linear_focus_amplitude_m/(2*pi/r.Tp_s);
gradients=zeros(4,2);
for j=1:4
    F=fft2(s.P(:,:,j));px=real(ifft2(1i*d.kx.*F));py=real(ifft2(1i*d.ky.*F));
    gradients(j,:)=[max(abs(px(:,[1,end])),[],'all'),max(abs(py([1,end],:)),[],'all')]/velocityScale;
end
report=struct('status',status,'final_time_s',finalTime,'kinematics_header',header, ...
    'checked_times_s',(0:frames-1)*dt,'kinematics_minus_EP_max_abs_eta_psi',errors, ...
    'initial_kinematics_phi_matches_EP',errors(1,2)<1e-8, ...
    'initial_surface_potential_source','EP_00000.bin or the prepared OceanWave3D.init, not unsolved t=0 volume phi', ...
    'top_sigma',sigma(end),'OW3D_resources',resource, ...
    'initial_normal_coordinate_surface_psi_gradient_over_gkpAwp',gradients, ...
    'boundary_scope','Surface psi gradients only, not exact volume normal-velocity or reduced-domain propagation certification', ...
    'production_run_started',false);
name='audit.json';if initialOnly,name='initial-output-audit.json';end
fid=fopen(fullfile(smoke,name),'w');assert(fid>=0);guard=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(report));disp(report);disp(errors);disp(gradients);
end

function data=read_record(fid,n)
assert(fread(fid,1,'int32')==8*n);data=fread(fid,n,'double');
assert(numel(data)==n && fread(fid,1,'int32')==8*n && all(isfinite(data)));
end

function [x,y,e,p]=read_ep(name)
fid=fopen(name,'r','ieee-le');assert(fid>=0);guard=onCleanup(@()fclose(fid));
assert(fread(fid,1,'int32')==8);sz=fread(fid,2,'int32').';assert(fread(fid,1,'int32')==8);
n=prod(sz);xy=read_record(fid,2*n);ep=read_record(fid,2*n);
x=reshape(xy(1:n),sz);y=reshape(xy(n+1:end),sz);e=reshape(ep(1:n),sz);p=reshape(ep(n+1:end),sz);
end
