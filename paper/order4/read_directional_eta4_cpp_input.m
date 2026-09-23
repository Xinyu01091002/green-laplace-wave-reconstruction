function input=read_directional_eta4_cpp_input(path)
%READ_DIRECTIONAL_ETA4_CPP_INPUT Read DIR4IN01 modes and analytic amplitudes.
arguments
    path (1,1) string
end
fid=fopen(path,'rb','ieee-le');assert(fid>=0,'Unable to open directional input.');
cleanup=onCleanup(@()fclose(fid));
assert(strcmp(char(fread(fid,8,'*char').'),'DIR4IN01'),'Directional input magic changed.');
modeCount=double(fread(fid,1,'*uint64'));nx=double(fread(fid,1,'*uint64'));
ny=double(fread(fid,1,'*uint64'));caseCount=double(fread(fid,1,'*uint64'));
gravity=fread(fid,1,'*double');depth=fread(fid,1,'*double');
deltaK=[fread(fid,1,'*double'),fread(fid,1,'*double')];
threshold=fread(fid,1,'*double');
modes=double(fread(fid,[2,modeCount],'*int32').');
raw=fread(fid,[2,modeCount*caseCount],'*double');
assert(size(raw,2)==modeCount*caseCount,'Directional input amplitude table is truncated.');
amplitudes=reshape(complex(raw(1,:),raw(2,:)),modeCount,caseCount);
assert(isempty(fread(fid,1,'*uint8')),'Directional input has trailing bytes.');
input=struct('modeCount',modeCount,'gridSize',[nx,ny],'caseCount',caseCount, ...
    'gravity',gravity,'depth',depth,'deltaK',deltaK,'threshold',threshold, ...
    'modeIndices',modes,'eta11Analytic',amplitudes);
end
