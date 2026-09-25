function report = run_ow3d_eta22_pilot(dataRoot,mf12Root,kph)
% Read-only four-phase OW3D pilot. All methods receive identical parent bins.
% Selection is fixed by GL domain and native temporal Nyquist, not reference error.
if nargin<3, kph=1; end
assert(ismember(kph,[0.5,0.6,0.8,1,2,5]));
caseName=['kh',strrep(sprintf('%g',kph),'.','p'),'_alpha1_akp002'];
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root); setup_green_laplace('MF12Root',mf12Root);
out=fullfile(root,'results','unidirectional_time_series',['ow3d_',caseName]);
if ~isfolder(out), mkdir(out); end
steps=(2800:4:4200)'; probe=3800; phases=[0,90,180,270];
raw=zeros(numel(steps),4); rawpsi=raw;
files=cell(numel(steps)*4,1); hashes=files; count=0;
metaFiles={}; metaHashes={}; xref=[]; yref=[];
for p=1:4
    folder=fullfile(dataRoot,sprintf('T_init-40_Tp_Alpha_1.0_Akp_002_kd%.1f_phi_%d',kph,phases(p)));
    inp=splitlines(string(fileread(fullfile(folder,'OceanWave3D.inp'))));
    grid=sscanf(inp(3),'%f'); timing=sscanf(inp(5),'%f'); gravity=sscanf(inp(6),'%f');
    rd=fileread(fullfile(folder,'OW_readme.txt'));
    kpToken=regexp(rd,'(?<![A-Za-z])kp=([\d.eE+-]+)','tokens','once');
    tpToken=regexp(rd,'(?<![A-Za-z])Tp=([\d.eE+-]+)','tokens','once');
    assert(~isempty(kpToken) && ~isempty(tpToken),'Missing kp/Tp metadata.');
    values=[grid(3),timing(2),gravity(1),str2double(kpToken{1}),str2double(tpToken{1})];
    if p==1, parameters=values; else, assert(isequal(values,parameters)); end
    for name={'OceanWave3D.inp','OW_readme.txt'}
        file=fullfile(folder,name{1}); metaFiles{end+1,1}=file; metaHashes{end+1,1}=sha256(file); %#ok<AGROW>
    end
    for it=1:numel(steps)
        file=fullfile(folder,sprintf('EP_%05d.bin',steps(it)));
        [x,y,e,psi]=read_ep(file);
        if isempty(xref), xref=x; yref=y; end
        assert(isequal(x,xref) && isequal(y,yref),'Spatial grids differ.');
        raw(it,p)=e(probe); rawpsi(it,p)=psi(probe);
        count=count+1; files{count}=file; hashes{count}=sha256(file);
    end
    fprintf('Read phase %d: %d complete snapshots\n',phases(p),numel(steps));
end
h=parameters(1); dt=4*parameters(2); g=parameters(3); kp=parameters(4); Tp=parameters(5);
assert(abs(kp*h-kph)<1e-5,'Depth label and input metadata disagree.');
t=steps*parameters(2); trel=t-t(1); N=numel(t);
% Temporal Hilbert sign matches paperplot_VWA_time_series.m.
H=imag(analytic_time(raw));
eta1phase=(raw(:,1)-raw(:,3)-H(:,2)+H(:,4))/4;
reference=(raw(:,1)-raw(:,2)+raw(:,3)-raw(:,4))/4;
% These are phase sectors, not strictly isolated perturbation orders.
F=fft(eta1phase)/N;
bins=(1:floor((N-1)/2))'; omegaAll=2*pi*bins/(N*dt);
kAll=zeros(size(omegaAll));
for j=1:numel(kAll)
    w=omegaAll(j); kAll(j)=fzero(@(k) g*k*tanh(k*h)-w^2,[0,max(1,2*w^2/g+2/h)]);
end
nyquist=pi/dt;
keep=(kAll*h>=0.3) & (2*omegaAll<nyquist);
% No energy-ranked truncation, gain fitting, or interpolation in k.
A=2*conj(F(bins(keep)+1)); omega=omegaAll(keep); k=kAll(keep);
E=exp(-1i*trel*omega.'); z=E*A; eta1=real(z);
inputFraction=sum(abs(F(bins(keep)+1)).^2)/sum(abs(F(bins+1)).^2);
assert(numel(A)>1 && all(isfinite(raw),'all'));
ranks=[6,8,12]; gl=zeros(N,numel(ranks)); audits=cell(1,numel(ranks));
for j=1:numel(ranks)
    [v,audits{j}]=gl_eta22_time_pairs(A,omega,k,h,kp,trel,ranks(j));
    gl(:,j)=real(v);
end
[mf12,mflinear,mfsource]=mf12_eta22_time(A,k,g,h,trel);
linearParity=norm(mflinear-eta1)/norm(eta1); assert(linearParity<1e-11);
B=(3-tanh(k*h).^2)./(4*tanh(k*h).^3);
vwa=real(z.*(E*(A.*k.*B)));
bp=(3-tanh(kp*h)^2)/(4*tanh(kp*h)^3);
walker=kp*bp*real(z.^2);
pred=[gl,mf12,vwa,walker];
names=["GL6","GL8","GL12","Spectral MF12","VWA","Walker"];
assert(all(isfinite(pred),'all'));
% Full window and fixed middle half; no peak search or alignment.
central=trel>=trel(end)/4 & trel<=3*trel(end)/4;
rows=cell(0,8);
for window=1:2
    mask=true(N,1); label="full";
    if window==2, mask=central; label="middle_half"; end
    for j=1:numel(names)
        v=pred(mask,j); r=reference(mask); diff=v-r;
        rows(end+1,:)={names(j),label,norm(diff)/norm(r),max(abs(diff))/max(abs(r)), ...
            norm(diff)/(norm(v)+norm(r)),norm(v)/norm(r),max(v),min(v)}; %#ok<AGROW>
    end
end
metrics=cell2table(rows,'VariableNames',{'method','window','relative_L2','relative_Linf','Q','norm_ratio','max_m','min_m'});
writetable(metrics,fullfile(out,'metrics.csv'));
writetable(table(files,hashes,'VariableNames',{'source_path','sha256'}),fullfile(out,'source_files.csv'));
writetable(table(metaFiles,metaHashes,'VariableNames',{'source_path','sha256'}),fullfile(out,'metadata_files.csv'));
trace=table(t,eta1phase,eta1,reference,gl(:,1),gl(:,2),gl(:,3),mf12,vwa,walker, ...
    'VariableNames',{'simulation_time_s','first_phase_sector_m','common_eta1_m','ow3d_second_phase_sector_m', ...
    'gl6_m','gl8_m','gl12_m','spectral_mf12_m','vwa_m','walker_m'});
writetable(trace,fullfile(out,'time_series.csv'));
report=struct('case',caseName,'depth_m',h,'kp_rad_m',kp,'Tp_s',Tp, ...
    'integration_dt_s',parameters(2),'sample_dt_s',dt,'sample_count',N, ...
    'probe_matlab_index',probe,'probe_x_m',xref(probe),'probe_y_m',yref(probe), ...
    'first_step',steps(1),'last_step',steps(end),'parent_count',numel(A), ...
    'retained_positive_frequency_energy_fraction',inputFraction, ...
    'excluded_dc_m',real(F(1)), ...
    'input_projection_relative_L2',norm(eta1-eta1phase)/norm(eta1phase), ...
    'minimum_parent_kh',min(k*h),'maximum_parent_kh',max(k*h), ...
    'omega_range_rad_s',[min(omega),max(omega)],'nyquist_rad_s',nyquist, ...
    'mf12_linear_parity',linearParity,'mf12_source',mfsource,'mf12_source_sha256',sha256(mfsource), ...
    'gl6_vs_mf12_relative_L2',norm(gl(:,1)-mf12)/norm(mf12), ...
    'gl8_vs_mf12_relative_L2',norm(gl(:,2)-mf12)/norm(mf12), ...
    'gl12_vs_mf12_relative_L2',norm(gl(:,3)-mf12)/norm(mf12), ...
    'reference_note','OW3D second phase sector, not pure perturbation order 2', ...
    'input_note','first phase sector projected to common declared parent support', ...
    'primary_rank',6,'comparison_tuned',false);
save(fullfile(out,'pilot.mat'),'report','metrics','t','raw','rawpsi','eta1phase','eta1', ...
    'reference','pred','names','A','omega','k','audits','steps');
fid=fopen(fullfile(out,'report.json'),'w'); assert(fid>=0);
fprintf(fid,'%s\n',jsonencode(report)); fclose(fid);
plot_ow3d_eta22_pilot(out);
disp(report); disp(metrics);
end

function z=analytic_time(x)
N=size(x,1); mask=zeros(N,1); mask(1)=1;
if mod(N,2)==0,mask(2:N/2)=2;mask(N/2+1)=1;else,mask(2:(N+1)/2)=2;end
z=ifft(fft(x).*mask);
end

function [x,y,e,p]=read_ep(file)
fid=fopen(file,'r','ieee-le'); assert(fid>=0,'Missing file: %s',file);
cleanup=onCleanup(@()fclose(fid));
assert(fread(fid,1,'int32')==8); dims=fread(fid,2,'int32');
assert(numel(dims)==2 && dims(2)==1 && fread(fid,1,'int32')==8);
n=prod(dims); marker=fread(fid,1,'int32'); assert(marker==16*n);
x=fread(fid,n,'double'); y=fread(fid,n,'double'); assert(fread(fid,1,'int32')==marker);
marker=fread(fid,1,'int32'); assert(marker==16*n);
e=fread(fid,n,'double'); p=fread(fid,n,'double');
assert(numel(p)==n && all(isfinite([x;y;e;p])));
closing=fread(fid,1,'int32'); assert(isempty(closing)||closing==marker);
assert(isempty(fread(fid,1,'uint8')));
end

function value=sha256(file)
fid=fopen(file,'r'); assert(fid>=0); cleanup=onCleanup(@()fclose(fid));
d=java.security.MessageDigest.getInstance('SHA-256');
while ~feof(fid), d.update(fread(fid,1024*1024,'*uint8')); end
value=lower(reshape(dec2hex(typecast(d.digest(),'uint8'),2).',1,[]));
end
