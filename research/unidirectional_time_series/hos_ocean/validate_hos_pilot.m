function validate_hos_pilot(out)
s=load(fullfile(out,'import_reference.mat')); nx=size(s.eta,1);ny=size(s.eta,2);N=nx*ny;
f=fopen(fullfile(out,'Results','3d.dat')); assert(f>=0); guard=onCleanup(@()fclose(f));
t=[];peaks=[]; etaStrip=[]; psiStrip=[];
while ~feof(f)
    line=fgetl(f); if ~ischar(line),break;end
    if ~startsWith(strtrim(line),'ZONE'),continue;end
    tok=regexp(line,'SOLUTIONTIME\s*=\s*([+\-\d.Ee]+)','tokens','once'); assert(~isempty(tok));
    t(end+1)=str2double(tok{1});
    if numel(t)==1,cols=4;else,cols=2;end
    a=fscanf(f,'%f',[cols,N]); assert(isequal(size(a),[cols,N])); assert(all(isfinite(a),'all'));
    if numel(t)==1
        x=reshape(a(1,:),nx,ny); y=reshape(a(2,:),nx,ny);
        xr=(0:nx-1)'*s.r.domain_m(1)/nx;yr=(0:ny-1)*s.r.domain_m(2)/ny;
        xerr=max(abs(x-xr),[],'all');yerr=max(abs(y-yr),[],'all');
        assert(xerr<.1 && yerr<.1,'Unexpected output coordinates');
        e0=a(3,:).';p0=a(4,:).';
        qeta=norm(e0-s.q(1,:).')/norm(s.q(1,:));qpsi=norm(p0-s.q(2,:).')/norm(s.q(2,:));
        raweta=norm(e0-s.eta(:))/norm(s.eta(:));rawpsi=norm(p0-s.psi(:))/norm(s.psi(:));
        fprintf('t0 versus quantized input: eta=%.17g psi=%.17g\n',qeta,qpsi);
        % Pinned HOS-ocean.f90 lines 386-390 zero both mean modes when eta mean is nonzero.
        % This is a source-prescribed import check, not a fitted correction to output.
        means=mean(s.q,2);
        expectedE=s.q(1,:).'-means(1);expectedP=s.q(2,:).'-means(2);
        halfUnitE=.5*10.^(floor(log10(max(abs(e0),realmin)))-5);
        halfUnitP=.5*10.^(floor(log10(max(abs(p0),realmin)))-5);
        boundE=halfUnitE+1e-12*max(abs(s.q(1,:)));
        boundP=halfUnitP+1e-12*max(abs(s.q(2,:)));
        assert(all(abs(e0-expectedE)<=boundE),'eta import exceeds documented zero-mode and output rounding bound');
        assert(all(abs(p0-expectedP)<=boundP),'psi import exceeds documented zero-mode and output rounding bound');
        assert(raweta<1e-5 && rawpsi<1e-5,'Raw import exceeds two ES12.5 rounding operations');
    end
    e=reshape(a(end-1,:),nx,ny);p=reshape(a(end,:),nx,ny);
    peaks(end+1,:)=[max(abs(e),[],'all'),max(abs(p),[],'all')];
    etaStrip(:,:,numel(t))=e(:,125:133);psiStrip(:,:,numel(t))=p(:,125:133);
end
assert(numel(t)==11 && max(abs(t-(0:.2:2)))<1e-10,'Incorrect output times');
r=struct('status','SHORT_PILOT_IO_AND_FINITE_OUTPUT_PASSED_NOT_PHYSICAL_ACCURACY_VALIDATION', ...
    'times_s',t,'max_eta_psi_by_time',peaks,'t0_relative_L2_vs_quantized_input',[qeta,qpsi], ...
    't0_relative_L2_vs_original_input',[raweta,rawpsi],'max_coordinate_text_rounding_m',[xerr,yerr], ...
    'input_quantized_means_eta_psi',means.', ...
    'zero_mode_rule','Pinned HOS source zeros eta and psi mean at import; raw output retained', ...
    'surface_file_numeric_format','ES12.5; about six significant digits', ...
    'no_gain_phase_shift_or_bias_fit',true,'output_columns','eta[m], true surface psi[m2/s]');
% Remove the preallocation-free empty first slice created by end+1 on [].
if size(etaStrip,3)==numel(t)+1,etaStrip=etaStrip(:,:,2:end);psiStrip=psiStrip(:,:,2:end);end
save(fullfile(out,'surface_strip.mat'),'t','etaStrip','psiStrip','r','-v7.3');
fid=fopen(fullfile(out,'validation.json'),'w');fprintf(fid,'%s',jsonencode(r,PrettyPrint=true));fclose(fid);disp(r);
end