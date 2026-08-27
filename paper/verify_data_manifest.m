function result = verify_data_manifest(data_root)
%VERIFY_DATA_MANIFEST Verify external matched-field inputs before plotting.
arguments
    data_root (1,1) string = string(getenv('GL_PSI33_DATA_ROOT'))
end
project_root=string(fileparts(fileparts(mfilename('fullpath'))));
manifest=jsondecode(fileread(fullfile(project_root,'paper','data_manifest.json')));
records=repmat(struct('case_id','','path','','bytes',0,'sha256','', ...
    'pass',false),numel(manifest.files),1);
for index=1:numel(manifest.files)
    expected=manifest.files(index);
    path=fullfile(data_root,strrep(string(expected.relative_path),'/',filesep));
    if ~isfile(path),error('Missing external field input: %s',path);end
    info=dir(path);digest=sha256_file(path);
    pass=info.bytes==expected.bytes && strcmpi(digest,expected.sha256);
    records(index)=struct('case_id',char(expected.case_id), ...
        'path',char(path),'bytes',info.bytes,'sha256',digest,'pass',pass);
    if ~pass,error('External field hash or size mismatch for %s.',path);end
end
result=struct('dataset',manifest.dataset,'records',records, ...
    'overall_pass',all([records.pass]));
end

function digest=sha256_file(path)
file_id=fopen(path,'rb');
if file_id<0,error('Could not read %s.',path);end
cleanup=onCleanup(@()fclose(file_id));
bytes=fread(file_id,Inf,'*uint8');
engine=java.security.MessageDigest.getInstance('SHA-256');
engine.update(typecast(bytes,'int8'));
raw=typecast(engine.digest(),'uint8');
digest=lower(reshape(dec2hex(raw,2).',1,[]));
end
