function run_directional_multiprobe(dataRoot)
% Fixed cross of probes chosen before target scoring, offsets in peak wavelengths.
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
offsets=[0,0;-1,0;1,0;0,-.5;0,.5];labels=["center","x_minus","x_plus","y_minus","y_plus"];
base=fullfile(root,'results','directional_joint_input_multiprobe');
extract_directional_joint_input(dataRoot,offsets,'directional_joint_input_multiprobe');
oldBase=fullfile(root,'results','directional_joint_input');
old=load(fullfile(oldBase,'extracted.mat'));d=load(fullfile(base,'extracted.mat'));
assert(isequal(old.initialSpectrum,d.initialSpectrum) && isequal(old.kx,d.kx) && isequal(old.ky,d.ky));
assert(isequal(old.probe1,d.probe1(:,1)) && isequal(old.eta2,d.eta2(:,1)) && isequal(old.psi2,d.psi2(:,1)));
assert(isequal(old.t,d.t) && old.metadata.h==d.metadata.h && old.metadata.g==d.metadata.g && old.metadata.kp==d.metadata.kp);
% The initial convention audit is global and probe-independent. Reuse only
% after the same source spectrum and center extraction have been verified.
copyfile(fullfile(oldBase,'initial_convention.json'),fullfile(base,'initial_convention.json'));
folders=strings(5,1);folders(1)=string(fullfile(oldBase,'verified_convention'));
for ip=2:5
    fprintf('RUN PROBE %s, x=%.6f, y=%.6f\n',labels(ip),d.metadata.probes(ip,1),d.metadata.probes(ip,2));
    run_directional_joint_pilot(base,char(labels(ip)),ip);
    folders(ip)=string(fullfile(base,char(labels(ip))));
end
rows=cell(0,10);
for ip=1:5
    f=load(fullfile(folders(ip),'joint_pilot.mat'));r=f.report;mt=f.metrics;
    for variable=["eta22","psi22"]
        match=string(mt.variable)==variable & string(mt.window)=="main_group";
        m=mt(match,:);get=@(name)m.relative_L2(string(m.method)==name);
        rows(end+1,:)={labels(ip),r.metadata.probe(1),r.metadata.probe(2),variable, ...
            get("Joint input, 15 deg"),get("Joint input, 7.5 deg"), ...
            get("Initial spectrum only"),get("Eta1 only, all directions zero"), ...
            r.conditioning{2}.observed_energy_weighted,nnz(f.focus)}; %#ok<AGROW>
    end
end
summary=cell2table(rows,'VariableNames',{'probe','x_m','y_m','variable','joint15_L2', ...
    'joint7p5_L2','initial_only_L2','single_direction_L2','weighted_condition7p5','main_samples'});
writetable(summary,fullfile(base,'summary.csv'));
save(fullfile(base,'summary.mat'),'summary','folders','labels','offsets');
plot_directional_multiprobe(base);disp(summary);
end
