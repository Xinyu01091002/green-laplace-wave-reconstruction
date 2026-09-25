function report=check_strip_projection()
% Check existing positive-kx extraction on nine full-x rows, using old data.
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
base=fullfile(root,'results','directional_amplitude_gate_akp012');
d=load(fullfile(base,'extracted.mat'),'metadata','probe1','t');
sources=readtable(fullfile(base,'source_files.csv'),'TextType','string', ...
    'Delimiter',',','ReadVariableNames',true,'NumHeaderLines',0);
assert(width(sources)==2);
assert(endsWith(string(sources{1,1}),"phi_shift_0"+filesep+"EP_00000.bin"));
times=[120,240,360];errors=zeros(size(times));parity=errors;absolute=errors;roundoffBound=errors;
for it=1:numel(times)
    first=[];
    for p=1:4
        stem=fileparts(sources{p,1});
        file=fullfile(stem,sprintf('EP_%05d.bin',round(times(it)/.4)));
        [X,Y,E]=read_e(file);
        if p==1
            x=X(:,1);y=Y(1,:);dx=median(diff(x));dy=median(diff(y));
            ix=find(x>=-1e-7 & x<d.metadata.domain(1)-dx/2);
            iy=find(y>=-1e-7 & y<d.metadata.domain(2)-dy/2);
            first=complex(zeros(numel(iy),numel(ix)));
        end
        first=first+.5*exp(-1i*(p-1)*pi/2)*E(ix,iy).';
    end
    nx=size(first,2);keep=[false,true(1,nx/2-1),false(1,nx/2)];
    F=fft2(first);F(:,~keep)=0;full=ifft2(F);
    rows=253:261;S=fft(first(rows,:),[],2);S(:,~keep)=0;strip=ifft(S,[],2);
    reference=full(rows,:);errors(it)=norm(strip-reference,'fro')/norm(reference,'fro');
    values=zeros(1,size(d.metadata.probes,1));
    for q=1:numel(values)
        [~,px]=min(abs(x(ix)-d.metadata.probes(q,1)));
        [~,py]=min(abs(y(iy)-d.metadata.probes(q,2)));
        values(q)=strip(py-rows(1)+1,px);
    end
    ref=d.probe1(d.t==times(it),:);assert(numel(ref)==numel(values));
    parity(it)=norm(values-ref)/norm(ref);
    absolute(it)=max(abs(values-ref));
    roundoffBound(it)=128*eps(max(abs(first(:))));
end
disp(table(times(:),errors(:),parity(:),absolute(:), ...
    'VariableNames',{'time_s','strip_error','saved_probe_error','absolute_error_m'}));
% Tail probes have tiny signal: use an input-amplitude roundoff bound there.
assert(max(errors)<1e-12 && all(absolute<=roundoffBound), ...
    'Strip projection did not match existing extraction within roundoff.');
report=struct('times_s',times,'strip_full_relative_L2',errors, ...
    'saved_probe_relative_L2',parity,'saved_probe_absolute_m',absolute, ...
    'roundoff_bound_m',roundoffBound,'positive_kx_only',true, ...
    'scope','Existing spatial first-sector extraction only; not a new temporal separator or OW3D run.');
out=fullfile(root,'artifacts','remote_campaign_design');if ~isfolder(out),mkdir(out);end
fid=fopen(fullfile(out,'strip_projection.json'),'w');assert(fid>=0);
c=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));disp(report);
end

function [x,y,e]=read_e(file)
fid=fopen(file,'r','ieee-le');assert(fid>=0);c=onCleanup(@()fclose(fid));
assert(fread(fid,1,'int32')==8);dims=fread(fid,2,'int32').';assert(fread(fid,1,'int32')==8);
n=prod(dims);assert(fread(fid,1,'int32')==16*n);
x=fread(fid,dims,'double');y=fread(fid,dims,'double');assert(fread(fid,1,'int32')==16*n);
assert(fread(fid,1,'int32')==16*n);e=fread(fid,dims,'double');assert(numel(e)==n);
end
