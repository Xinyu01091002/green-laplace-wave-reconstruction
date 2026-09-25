function audit_directional_initial_convention(dataRoot)
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));out=fullfile(root,'results','directional_joint_input');
d=load(fullfile(out,'extracted.mat'));m=d.metadata;raw=[];
for phase=[0,90,180,270]
    file=fullfile(dataRoot,sprintf('kd1.0_spread_25_Akp_0.02_phi_shift_%d',phase),'EP_00000.bin');
    fid=fopen(file,'r','ieee-le');assert(fid>=0);assert(fread(fid,1,'int32')==8);dims=fread(fid,2,'int32').';
    assert(fread(fid,1,'int32')==8);n=prod(dims);assert(fread(fid,1,'int32')==16*n);
    X=fread(fid,dims,'double');Y=fread(fid,dims,'double');assert(fread(fid,1,'int32')==16*n);
    assert(fread(fid,1,'int32')==16*n);fread(fid,n,'double');p=fread(fid,dims,'double');fclose(fid);
    x=X(:,1);y=Y(1,:).';dx=median(diff(x));dy=median(diff(y));
    ix=x>=-1e-7 & x<m.domain(1)-dx/2;iy=y>=-1e-7 & y<m.domain(2)-dy/2;
    p=p(ix,iy).';if isempty(raw),raw=complex(zeros(size(p)));end
    raw=raw+.5*exp(-1i*deg2rad(phase))*p;
end
P=fft2(raw)/numel(raw);r=hypot(d.kx,d.ky);w=sqrt(m.g*r.*tanh(r*m.h));
use=d.kx>0&r*m.h>=.3&w<pi/m.dt;
eta=d.initialSpectrum(use);psi=P(use);factor=m.g./w(use);
minus=norm(psi+1i*factor.*eta)/norm(psi);plus=norm(psi-1i*factor.*eta)/norm(psi);
% This diagnoses the stored physical propagation convention from first-order
% initialization only. No high-order target or shifted time record is used.
assert(min(plus,minus)<.05 && max(plus,minus)>1.9,'Initial propagation convention is unresolved.');
storedPositiveKHasPositiveTime=(plus<.05);
report=struct('minus_time_polarization_relative_error',minus, ...
    'plus_time_polarization_relative_error',plus, ...
    'stored_positive_k_has_positive_time',storedPositiveKHasPositiveTime, ...
    'basis','initial eta1/psi1 only; no higher-order reference fit');
fid=fopen(fullfile(out,'initial_convention.json'),'w');fprintf(fid,'%s',jsonencode(report));fclose(fid);disp(report);
end
