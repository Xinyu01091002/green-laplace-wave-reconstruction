function output=extract_directional_joint_input(dataRoot)
% Read-only directional OW3D pilot; future fields supply only probe eta1.
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
out=fullfile(root,'results','directional_joint_input');if ~isfolder(out),mkdir(out);end
steps=(0:10:900)';phases=[0,90,180,270];N=numel(steps);
probe1=complex(zeros(N,1));eta2=zeros(N,1);psi2=eta2;
files=cell(N*4,1);hashes=files;counter=0;
for it=1:N
    first=[];second=[];potential=[];
    for p=1:4
        folder=fullfile(dataRoot,sprintf('kd1.0_spread_25_Akp_0.02_phi_shift_%d',phases(p)));
        if it==1
            lines=splitlines(string(fileread(fullfile(folder,'OceanWave3D.inp'))));
            grid=sscanf(lines(3),'%f');clock=sscanf(lines(5),'%f');grav=sscanf(lines(6),'%f');
            values=[grid(1:5).',clock(2),grav(1)];
            if p==1,parameters=values;else,assert(isequal(parameters,values));end
        end
        file=fullfile(folder,sprintf('EP_%05d.bin',steps(it)));
        [X,Y,e,phi]=read_ep(file);
        if it==1 && p==1
            x=X(:,1);y=Y(1,:).';dx=median(diff(x));dy=median(diff(y));
            ix=find(x>=-1e-7 & x<grid(1)-dx/2);iy=find(y>=-1e-7 & y<grid(2)-dy/2);
            assert(numel(ix)==grid(4)-1 && numel(iy)==grid(5)-1);
            x=x(ix);y=y(iy);nx=numel(x);ny=numel(y);
            [~,px]=min(abs(x-grid(1)/2));[~,py]=min(abs(y-grid(2)/2));
            [kx,ky]=meshgrid(2*pi/grid(1)*[0:nx/2-1,-nx/2:-1], ...
                2*pi/grid(2)*[0:ny/2-1,-ny/2:-1]);
            forward=kx>0;
        end
        assert(all(isfinite([e(:);phi(:)])));
        e=e(ix,iy).';phi=phi(ix,iy).';
        if p==1,first=complex(zeros(ny,nx));second=zeros(ny,nx);potential=second;end
        first=first+.5*exp(-1i*deg2rad(phases(p)))*e;
        second=second+(-1)^(p-1)*e/4;
        potential=potential+(-1)^(p-1)*phi/4;
        counter=counter+1;files{counter}=file;hashes{counter}=hash_file(file);
    end
    spectrum=fft2(first);spectrum(~forward)=0;
    f=ifft2(spectrum);probe1(it)=f(py,px);
    eta2(it)=second(py,px);psi2(it)=potential(py,px);
    if it==1
        initialSpectrum=spectrum/(nx*ny);
        initialBoundaryRatio=norm([f(1,:),f(end,:),f(:,1).',f(:,end).'])/norm(f(:));
    end
    if mod(it-1,15)==0,fprintf('Read directional snapshot %d/%d\n',it,N);end
end
g=parameters(7);h=parameters(3);kp=.0279;t=steps*parameters(6);
metadata=struct('h',h,'g',g,'kp',kp,'kph',h*kp,'alpha','not certified from historical generator', ...
    'spread_label_degrees',25,'Akp',.02,'dt',10*parameters(6),'probe',[x(px),y(py)], ...
    'grid',[ny,nx],'domain',[parameters(1),parameters(2)],'initial_boundary_ratio',initialBoundaryRatio, ...
    'first_input','positive-kx projection of four-phase first sector, not exact eta11 certification', ...
    'snapshot_reference','raw second phase sector at probe; no temporal filtering');
output=fullfile(out,'extracted.mat');save(output,'initialSpectrum','kx','ky','probe1','eta2','psi2','t','metadata','-v7.3');
writetable(table(files,hashes,'VariableNames',{'source','sha256'}),fullfile(out,'source_files.csv'));
disp(metadata);
end

function [x,y,e,p]=read_ep(file)
fid=fopen(file,'r','ieee-le');assert(fid>=0);c=onCleanup(@()fclose(fid));
assert(fread(fid,1,'int32')==8);dims=fread(fid,2,'int32').';assert(fread(fid,1,'int32')==8);
n=prod(dims);assert(fread(fid,1,'int32')==16*n);
x=fread(fid,dims,'double');y=fread(fid,dims,'double');assert(fread(fid,1,'int32')==16*n);
assert(fread(fid,1,'int32')==16*n);e=fread(fid,dims,'double');p=fread(fid,dims,'double');
assert(numel(p)==n);last=fread(fid,1,'int32');assert(isempty(last)||last==16*n);
end

function value=hash_file(file)
fid=fopen(file,'r');assert(fid>=0);c=onCleanup(@()fclose(fid));
d=java.security.MessageDigest.getInstance('SHA-256');while ~feof(fid),d.update(fread(fid,1024*1024,'*uint8'));end
value=lower(reshape(dec2hex(typecast(d.digest(),'uint8'),2).',1,[]));
end
