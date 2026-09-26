function validate_short(root)
s=load(fullfile(root,'inputs','initial_fields.mat'),'E','P');settings=jsondecode(fileread(fullfile(root,'staging4','settings.json')));N=1024*256;ref=[reshape(s.E(:,:,1).',N,1),reshape(s.P(:,:,1).',N,1)];report=struct();
for np=[4 8]
 field=zeros(1024,256,2,11);ny=256/np;
 for rank=0:np-1
  file=fullfile(root,sprintf('bench_np%d',np),'Results',sprintf('3d_%03d.dat',rank));[a,t,xy]=readfield(file,1024*ny);ys=rank*ny+(1:ny);assert(max(abs(reshape(xy(:,2),1024,ny)-(ys-1)*settings.domain_m(2)/256),[],'all')<1e-9);field(:,ys,:,:)=reshape(a,1024,ny,2,11);
 end
 field=reshape(field,N,2,11);err=sqrt(sum((field(:,:,1)-ref).^2))./sqrt(sum(ref.^2));assert(max(err)<1e-10);
 f=fopen(fullfile(root,sprintf('bench_np%d',np),'Results','probes.dat'));while ~feof(f),line=fgetl(f);if startsWith(line,'VARIABLES'),break;end;end;probes=fscanf(f,'%f',[6 Inf]).';fclose(f);assert(isequal(size(probes),[11 6]));
 for j=1:5,index=settings.probe_indices(j,1)+(settings.probe_indices(j,2)-1)*1024;assert(max(abs(squeeze(field(index,1,:))-probes(:,j+1)))<1e-10);end
 report.(sprintf('np%d_initial_relative_L2',np))=err;if np==4,baseline=field;else,parallel=field;end
end
parity=zeros(11,2);for it=1:11,parity(it,:)=sqrt(sum((parallel(:,:,it)-baseline(:,:,it)).^2))./sqrt(sum(baseline(:,:,it).^2));end
assert(max(parity,[],'all')<1e-10);report.status='PASSED_2_SECOND_MPI4_MPI8_AND_PROBE_CHECK';report.max_relative_L2=max(parity,[],1);report.times_s=t;
f=fopen(fullfile(root,'benchmark-validation.json'),'w');fprintf(f,'%s',jsonencode(report,PrettyPrint=true));fclose(f);disp(report);
end
function [a,t,xy]=readfield(file,N)
f=fopen(file);assert(f>=0);cl=onCleanup(@()fclose(f));a=zeros(N,2,11);t=[];xy=[];
while ~feof(f)
 line=fgetl(f);if ~ischar(line),break;end;if ~startsWith(line,'ZONE'),continue;end;tok=regexp(line,'SOLUTIONTIME\s*=\s*([+\-\d.Ee]+)','tokens','once');t(end+1)=str2double(tok{1});it=numel(t);assert(it<=11);if it==1,nc=4;else,nc=2;end;v=fscanf(f,'%f',[nc,N]);assert(isequal(size(v),[nc N])&&all(isfinite(v),'all'));if it==1,xy=v(1:2,:).';end;a(:,:,it)=v(end-1:end,:).';
end
assert(numel(t)==11&&max(abs(t-(0:.2:2)))<1e-12);
end