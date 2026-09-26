function validate_hos_precision(root)
% Preregistered checks: full-input t0 <1e-12 relative L2; coordinates <1e-9 m;
% same-input dynamics agree within old output rounding plus 1e-11 peak.
[a,ta,ca]=read_surface(fullfile(root,'case-baseline-low','Results','3d.dat'));
[b,tb,cb]=read_surface(fullfile(root,'case-precise-low','Results','3d.dat'));
[c,tc,cc]=read_surface(fullfile(root,'case-precise-full','Results','3d.dat'));
assert(max(abs(ta-tb))<1e-12 && max(abs(ta-tc))<1e-12);
s=load(fullfile(root,'case-precise-full','import_reference.mat'));
ref=[s.eta(:),s.psi(:)];raw=relative(c(:,:,1)-ref,ref);
assert(all(raw<1e-12),'Double-precision initial import failed');
nx=size(s.eta,1);ny=size(s.eta,2);
x=repmat((0:nx-1)'*s.r.domain_m(1)/nx,1,ny);y=repmat((0:ny-1)*s.r.domain_m(2)/ny,nx,1);
coord=max(abs(cc-[x(:),y(:)]),[],1);assert(all(coord<1e-9));
lowComparison=zeros(numel(ta),2);inputEffect=lowComparison;
for it=1:numel(ta)
    lowComparison(it,:)=relative(a(:,:,it)-b(:,:,it),b(:,:,it));
    inputEffect(it,:)=relative(c(:,:,it)-b(:,:,it),c(:,:,it));
    for k=1:2
        halfUnit=.5*10.^(floor(log10(max(abs(a(:,k,it)),realmin)))-5);
        bound=halfUnit+1e-11*max(abs(b(:,k,it)));
        assert(all(abs(a(:,k,it)-b(:,k,it))<=bound),'Matched build differs beyond legacy rounding bound');
    end
end
r=struct('status','PRECISION_IMPORT_AND_MATCHED_BUILD_IO_CHECKS_PASSED', ...
    'times_s',tc,'t0_relative_L2_vs_original_eta_psi',raw, ...
    'input_quantization_relative_L2',[s.r.quantization_eta_relative_L2,s.r.quantization_psi_relative_L2], ...
    't0_max_abs_error_eta_psi',max(abs(c(:,:,1)-ref),[],1), ...
    'coordinate_max_abs_error_m',coord, ...
    'same_quantized_input_baseline_vs_precise_relative_L2_by_time',lowComparison, ...
    'precise_output_low_vs_full_input_relative_L2_by_time',inputEffect, ...
    'same_input_check','Pointwise legacy ES12.5 half-unit rounding plus 1e-11 peak bound', ...
    'all_frames_finite',true,'no_alignment_or_fitted_correction',true, ...
    'scope','2 s wavegroup phase 0, not long-time or physical-accuracy validation');
eta=reshape(c(:,1,:),nx,ny,[]);psi=reshape(c(:,2,:),nx,ny,[]);etaStrip=eta(:,125:133,:);psiStrip=psi(:,125:133,:);t=tc;
save(fullfile(root,'case-precise-full','surface_strip.mat'),'etaStrip','psiStrip','t','r','-v7.3');
f=fopen(fullfile(root,'precision-validation.json'),'w');fprintf(f,'%s',jsonencode(r,PrettyPrint=true));fclose(f);disp(r);
end
function z=relative(a,b)
z=sqrt(sum(abs(a).^2,1))./sqrt(sum(abs(b).^2,1));
end
function [data,t,coord]=read_surface(path)
f=fopen(path);assert(f>=0);guard=onCleanup(@()fclose(f));N=1024*256;data=zeros(N,2,11);t=[];coord=[];
while ~feof(f)
    line=fgetl(f);if ~ischar(line),break;end
    if ~startsWith(strtrim(line),'ZONE'),continue;end
    tok=regexp(line,'SOLUTIONTIME\s*=\s*([+\-\d.Ee]+)','tokens','once');assert(~isempty(tok));
    t(end+1)=str2double(tok{1});it=numel(t);assert(it<=11);
    if it==1,cols=4;else,cols=2;end
    values=fscanf(f,'%f',[cols,N]);assert(isequal(size(values),[cols,N]) && all(isfinite(values),'all'));
    if it==1,coord=values(1:2,:).';end
    data(:,:,it)=values(end-1:end,:).';
end
assert(numel(t)==11 && max(abs(t-(0:.2:2)))<1e-12);
end