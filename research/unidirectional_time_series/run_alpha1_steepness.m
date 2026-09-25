function run_alpha1_steepness(escRoot,mf12Root,boundaryCorrected)
% Raw compact records only; no ESC reconstruction or correction is imported.
if nargin<3,boundaryCorrected=true;end
data=fullfile(escRoot,'outputs','remote_runs','esc-ow3d-full-20260827T015024Z','compact-results','compact');
inputFamily='vwa_general_mf12_fourphase'; prefix='compact';
if boundaryCorrected
    data=fullfile(escRoot,'outputs','remote_runs','esc-ow3d-kd1-boundary-corrected-20260901T013632Z','compact');
    inputFamily='vwa_general_mf12_fourphase_kd1_boundary_corrected'; prefix='boundary';
end
manifestFile=fullfile(data,'probe_manifest.json');
manifest=jsondecode(fileread(manifestFile));
configFile=fullfile(escRoot,'configs','ow3d_esc_input_plan.json');
config=jsondecode(fileread(configFile));
kp=config.peak_wavenumber;
for steepness=[.02,.12]
    compact=struct(); compact.steepness=steepness;
    compact.caseName=sprintf('%s_kh1_alpha1_akp%03d',prefix,round(100*steepness));
    compact.files=cell(4,1); compact.hashes=cell(4,1);
    compact.metaFiles={manifestFile;configFile};
    compact.metaHashes={hash_file(manifestFile);hash_file(configFile)};
    if boundaryCorrected
        plan=fullfile(escRoot,'configs','ow3d_esc_kd1_boundary_corrected_plan.json');
        compact.metaFiles{end+1,1}=plan; compact.metaHashes{end+1,1}=hash_file(plan);
    end
    for p=1:4
        name=sprintf('T_init-40_Tp_Alpha_1.0_Akp_%03d_kd1.0_phi_%d',round(100*steepness),(p-1)*90);
        idx=find(strcmp({manifest.cases.case_name},name)); assert(isscalar(idx));
        entry=manifest.cases(idx); file=fullfile(data,entry.csv_relative_path);
        digest=hash_file(file); assert(strcmpi(digest,entry.csv_sha256));
        tab=readtable(file); assert(height(tab)==entry.sample_count);
        inpFile=fullfile(escRoot,'outputs','ow3d_inputs',inputFamily,name,'OceanWave3D.inp');
        inp=splitlines(string(fileread(inpFile)));
        grid=sscanf(inp(3),'%f'); timing=sscanf(inp(5),'%f'); gravity=sscanf(inp(6),'%f');
        params=[grid(3),timing(2),gravity(1),kp,entry.peak_period];
        assert(abs(timing(2)-entry.dt)<1e-10 && abs(kp*grid(3)-1)<1e-5);
        assert(max(abs(tab.time_relative_s-entry.initial_time-tab.step*entry.dt))<1e-8);
        if p==1
            compact.parameters=params; compact.steps=tab.step;
            compact.probe=entry.probe_index_zero_based+1; compact.x=entry.probe_x;
            compact.raw=zeros(height(tab),4); compact.rawpsi=compact.raw;
        else
            assert(isequal(compact.parameters,params) && isequal(compact.steps,tab.step));
            assert(compact.x==entry.probe_x && compact.probe==entry.probe_index_zero_based+1);
        end
        compact.raw(:,p)=tab.eta; compact.rawpsi(:,p)=tab.psi;
        compact.files{p}=file; compact.hashes{p}=digest;
        compact.metaFiles{end+1,1}=inpFile; compact.metaHashes{end+1,1}=hash_file(inpFile);
    end
    assert(all(isfinite(compact.raw),'all') && all(isfinite(compact.rawpsi),'all'));
    run_ow3d_eta22_pilot('',mf12Root,1,compact);
    root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
    plot_ow3d_main_group(fullfile(root,'results','unidirectional_time_series',['ow3d_',compact.caseName]));
end
end

function value=hash_file(file)
fid=fopen(file,'r'); assert(fid>=0); cleanup=onCleanup(@()fclose(fid));
digest=java.security.MessageDigest.getInstance('SHA-256');
while ~feof(fid),digest.update(fread(fid,1024*1024,'*uint8'));end
value=lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2).',1,[]));
end
