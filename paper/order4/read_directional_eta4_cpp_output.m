function result=read_directional_eta4_cpp_output(path)
%READ_DIRECTIONAL_ETA4_CPP_OUTPUT Read paired directional eta4 sectors.
arguments
    path (1,1) string
end
fid=fopen(path,'r','ieee-le');assert(fid>=0,'Unable to open C++ output.');
cleanup=onCleanup(@()fclose(fid));
magic=char(fread(fid,8,'*uint8').');
allowed={'DIR40O01','DIR42O01','DIR44O01'};
assert(ismember(magic,allowed),'Directional C++ output magic changed.');
nx=double(fread(fid,1,'*uint64'));ny=double(fread(fid,1,'*uint64'));
caseCount=double(fread(fid,1,'*uint64'));
total=double(fread(fid,1,'*uint64'));kept=double(fread(fid,1,'*uint64'));
near=double(fread(fid,1,'*uint64'));
strictZeroOutput=double(fread(fid,1,'*uint64'));
strictZeroIntermediate=double(fread(fid,1,'*uint64'));
seconds=fread(fid,1,'double');maximumResidual=fread(fid,1,'double');
minimumDetuning=fread(fid,1,'double');
maximumStrictZeroSource=fread(fid,1,'double');
count=nx*ny*caseCount;
arrays=cell(5,1);
for index=1:5
    values=fread(fid,[2,count],'double');
    assert(size(values,2)==count,'Directional C++ output is truncated.');
    arrays{index}=reshape(complex(values(1,:),values(2,:)),nx,ny,caseCount);
end
assert(isempty(fread(fid,1,'uint8')),'Directional C++ output has trailing bytes.');
result=struct('magic',magic,'gridSize',[nx,ny],'caseCount',caseCount, ...
    'totalCandidates',total,'retainedCount',kept,'nearResonantCount',near, ...
    'strictZeroOutputCount',strictZeroOutput, ...
    'strictZeroIntermediateCount',strictZeroIntermediate, ...
    'elapsedSeconds',seconds,'maximumPairedResidual',maximumResidual, ...
    'minimumNormalizedDetuning',minimumDetuning, ...
    'maximumStrictZeroSource',maximumStrictZeroSource, ...
    'etaFixedB',arrays{1},'psiFixedB',arrays{2}, ...
    'etaCoordinate',arrays{3},'psiCoordinate',arrays{4}, ...
    'retainedNear',arrays{5},'allFinite',all(isfinite(vertcat(arrays{:})),'all'));
end
