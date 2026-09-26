function validate_mpi_check(root)
base='/home/lxy/green-laplace-unidirectional-time-series/artifacts/hos-ocean-v2.1.0/precision-v1/case-precise-full';
[ref,tr]=readfile(fullfile(base,'Results','3d.dat'),1024*256);report=struct();
for np=[1 4]
    field=zeros(1024,256,2,11);ny=256/np;
    for rank=0:np-1
        [a,t,xy]=readfile(fullfile(root,sprintf('np%d',np),'Results',sprintf('3d_%03d.dat',rank)),1024*ny);
        assert(max(abs(t-tr))<1e-12);ys=rank*ny+(1:ny);
        expectedY=(ys-1)*4504.075489017624/256;Y=reshape(xy(:,2),1024,ny);assert(max(abs(Y-expectedY),[],'all')<1e-9,'Unexpected MPI slab coordinate order');
        field(:,ys,:,:)=reshape(a,[1024,ny,2,11]);
    end
    data=reshape(field,1024*256,2,11);delta=zeros(11,2);maxabs=delta;
    for it=1:11
        delta(it,:)=sqrt(sum((data(:,:,it)-ref(:,:,it)).^2))./sqrt(sum(ref(:,:,it).^2));maxabs(it,:)=max(abs(data(:,:,it)-ref(:,:,it)),[],1);
    end
    assert(max(delta,[],'all')<1e-10,'MPI result differs materially from serial reference');
    report.(sprintf('np%d_vs_serial',np))=struct('relative_L2_by_frame',delta,'max_abs_by_frame',maxabs);
    if np==1,one=data;else,four=data;end
end
parity=zeros(11,2);
for it=1:11,parity(it,:)=sqrt(sum((four(:,:,it)-one(:,:,it)).^2))./sqrt(sum(one(:,:,it).^2));end
assert(max(parity,[],'all')<1e-10);
report.np4_vs_np1_relative_L2_by_frame=parity;report.status='MPI_1_AND_4_AGREE_WITH_SERIAL';report.times_s=tr;report.no_alignment_or_fit=true;
f=fopen(fullfile(root,'validation.json'),'w');fprintf(f,'%s',jsonencode(report,PrettyPrint=true));fclose(f);disp(report);disp(max(parity,[],1));
end
function [a,t,xy]=readfile(file,N)
f=fopen(file);assert(f>=0,file);cleanup=onCleanup(@()fclose(f));t=[];a=zeros(N,2,11);xy=[];
while ~feof(f)
    line=fgetl(f);if ~ischar(line),break;end;if ~startsWith(strtrim(line),'ZONE'),continue;end
    tok=regexp(line,'SOLUTIONTIME\s*=\s*([+\-\d.Ee]+)','tokens','once');assert(~isempty(tok));t(end+1)=str2double(tok{1});it=numel(t);assert(it<=11);
    if it==1,nc=4;else,nc=2;end;v=fscanf(f,'%f',[nc,N]);assert(isequal(size(v),[nc,N])&&all(isfinite(v),'all'));
    if it==1,xy=v(1:2,:).';end;a(:,:,it)=v(end-1:end,:).';
end
assert(numel(t)==11&&max(abs(t-(0:.2:2)))<1e-12);
end