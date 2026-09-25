function matrix=export_ow3d16_inputs(base,out)
% Materialize sixteen distinct native-wall OW3D cases from frozen eta/psi.
% No taper, renormalization, new transfer coefficient or initial-field change.
assert(~isfolder(out),'Refusing to overwrite a prior export.');mkdir(out);
families={'wavegroup','random_s20260925','random_s20260926','random_s20260927'};
sources={'wavegroup-v2','random-phase-only-v1','random-seed20260926','random-seed20260927'};
template=splitlines(string(fileread(fullfile(base,'wavegroup-v4','wavegroup_kpd1_akp012_phi000','OceanWave3D.inp'))));
assert(startsWith(template(2),'-1 0 1000'));
assert(startsWith(template(10),'1 1025 1 125 133 1 2 1101 1'));
rows=cell(16,1);n=0;phases=[0,90,180,270];C0=[];records=cell(4,1);
for f=1:4
    source=fullfile(base,sources{f},'initial_fields.mat');d=load(source);
    if f==1,C0=d.C;else,assert(max(abs(abs(d.C)-abs(C0)))/max(abs(C0))<1e-14);end
    assert(isequal(size(d.E),[256,1024,4]) && isequal(size(d.P),size(d.E)));
    assert(all(isfinite(d.E),'all') && all(isfinite(d.P),'all'));
    r=d.report.design;assert(abs(r.kpd-1)<1e-12 && r.OW3D_nodes(3)==17);
    seed=[];if f>1,seed=d.report.randomization.seed;end
    records{f}=struct('family',families{f},'seed',seed,'source',source, ...
        'parent_count',numel(d.C),'linear_Hs_4sigma_m',4*sqrt(sum(abs(d.C).^2)/2));
    for p=1:4
        n=n+1;id=sprintf('%s_phi%03d',families{f},phases(p));folder=fullfile(out,id);mkdir(folder);
        E=d.E(:,:,p);P=d.P(:,:,p);e=E([1:end,1],[1:end,1]).';v=P([1:end,1],[1:end,1]).';
        fid=fopen(fullfile(folder,'OceanWave3D.init'),'w');assert(fid>=0);guard=onCleanup(@()fclose(fid));
        fprintf(fid,'Independent MF12 order-2; %s; native-wall OW3D finite-time experiment\n',id);
        fprintf(fid,'%.17g %.17g 1025 257 0\n',r.domain_m(1),r.domain_m(2));
        fprintf(fid,'%.17g %.17g\n',[e(:),v(:)].');clear guard;
        lines=template;lines(1)="GL compact "+id+"; native walls; no taper or rescaling";
        fid=fopen(fullfile(folder,'OceanWave3D.inp'),'w');assert(fid>=0);guard=onCleanup(@()fclose(fid));
        fprintf(fid,'%s',join(lines,newline));clear guard;
        fid=fopen(fullfile(folder,'OceanWave3D.init'),'r');assert(fid>=0);guard=onCleanup(@()fclose(fid));
        fgetl(fid);header=sscanf(fgetl(fid),'%f');pairs=fscanf(fid,'%f',[2,Inf]);clear guard;
        assert(isequal(header(3:5).',[1025,257,0]) && size(pairs,2)==numel(e));
        assert(max(abs(pairs(1,:).'-e(:)))<1e-13 && max(abs(pairs(2,:).'-v(:)))<1e-12);
        rows{n}=struct('case_id',id,'family',families{f},'phase_degrees',phases(p),'seed',seed, ...
            'source_fields',source,'source_phase_index',p,'directory',folder, ...
            'expected_final_time_s',220,'output_dt_s',.2,'expected_saved_samples_with_initial',1101, ...
            'expected_kinematics_header',[1,1025,1,125,133,1,2,1101,1], ...
            'native_nodes',[1025,257,17],'g',r.g,'depth_m',r.h, ...
            'boundary','native OW3D walls; not periodic-HOS-equivalent; boundary influence remains to be assessed');
    end
end
assert(n==16);matrix=struct('cases',[rows{:}],'families',[records{:}], ...
    'scope','One wavegroup and three phase-only random realizations, four phases each', ...
    'no_new_depth_steepness_dt_cases',true,'no_taper_or_amplitude_fit',true, ...
    'complete_propagation_validated',false,'initial_eta_psi_source','EP_00000.bin; kinematics starts at t=.2 s');
fid=fopen(fullfile(out,'case_matrix.json'),'w');assert(fid>=0);guard=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(matrix));fprintf('Exported %d distinct phase case definitions.\n',n);
end
