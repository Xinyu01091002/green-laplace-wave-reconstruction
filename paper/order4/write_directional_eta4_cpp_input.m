function manifest=write_directional_eta4_cpp_input(path,modeIndices, ...
        amplitudes,gravity,depth,deltaK,gridSize,resonanceThreshold)
%WRITE_DIRECTIONAL_ETA4_CPP_INPUT Frozen directional C++ binary adapter.
arguments
    path (1,1) string
    modeIndices (:,2) double {mustBeInteger,mustBeFinite}
    amplitudes (:,:) double {mustBeFinite}
    gravity (1,1) double {mustBePositive,mustBeFinite}
    depth (1,1) double {mustBePositive,mustBeFinite}
    deltaK (1,2) double {mustBePositive,mustBeFinite}
    gridSize (1,2) double {mustBeInteger,mustBePositive}
    resonanceThreshold (1,1) double {mustBePositive,mustBeFinite}=1e-6
end
assert(size(amplitudes,1)==size(modeIndices,1), ...
    'Each C++ case must share the directional mode table.');
fid=fopen(path,'w','ieee-le');assert(fid>=0,'Unable to open C++ input.');
cleanup=onCleanup(@()fclose(fid));
fwrite(fid,uint8('DIR4IN01'),'uint8');
fwrite(fid,uint64(size(modeIndices,1)),'uint64');
fwrite(fid,uint64(gridSize(1)),'uint64');
fwrite(fid,uint64(gridSize(2)),'uint64');
fwrite(fid,uint64(size(amplitudes,2)),'uint64');
fwrite(fid,[gravity,depth,deltaK,resonanceThreshold],'double');
fwrite(fid,int32(modeIndices.'),'int32');
for caseIndex=1:size(amplitudes,2)
    values=amplitudes(:,caseIndex);
    fwrite(fid,[real(values).';imag(values).'],'double');
end
manifest=struct('magic','DIR4IN01','modeCount',size(modeIndices,1), ...
    'gridSize',gridSize,'caseCount',size(amplitudes,2), ...
    'gravity',gravity,'depth',depth,'deltaK',deltaK, ...
    'resonanceThreshold',resonanceThreshold,'path',char(path));
end
