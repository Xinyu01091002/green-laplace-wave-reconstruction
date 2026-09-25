function run_directional_amplitude_gate(dataRoot,Akp)
% User-defined eligibility: sampled OW3D eta22 peak >= same-x centerline / 3.
% This selects a useful signal region, not GL coefficients or error minima.
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
if nargin<2,Akp=.02;end
assert(ismember(Akp,[.02,.12]));
cfg=struct('kph',1,'spread',25,'Akp',Akp,'lastStep',900);
if Akp==.02
    oldBase=fullfile(root,'results','directional_joint_input');centerFolder=fullfile(oldBase,'verified_convention');
    outputName='directional_amplitude_gate';
else
    oldBase=fullfile(root,'results','directional_sweep','kh1_s25_a012');centerFolder=fullfile(oldBase,'center');
    outputName='directional_amplitude_gate_akp012';cfg.lastStep=1200;
end
center=load(fullfile(centerFolder,'joint_pilot.mat'));assert(center.report.metadata.Akp==Akp);
referencePeak=max(abs(center.d.eta2));referenceRange=range(center.d.eta2);
centerMainPeak=max(abs(center.d.eta2(center.focus)));
folders=string(centerFolder);labels="center";
if Akp==.02
    near=load(fullfile(root,'results','directional_joint_input_multiprobe','summary.mat'));
    farNames=["x_m3","x_p3","y_m1","y_p1","y_m1p5","y_p1p5"];
    folders=[string(near.folders(:));strings(numel(farNames),1)];labels=[string(near.labels(:));farNames(:)];
    farBase=fullfile(root,'results','directional_sweep','kh1_s25_a002');
    for j=1:numel(farNames)
        folders(numel(near.folders)+j)=string(fullfile(farBase,char(farNames(j))));
    end
end
rows=cell(0,13);
for j=1:numel(folders)
    d=load(fullfile(folders(j),'joint_pilot.mat'));
    rows(end+1,:)=record(labels(j),d,referencePeak,referenceRange,centerMainPeak,center.report.metadata.probe); %#ok<AGROW>
end
base=fullfile(root,'results',outputName);if ~isfolder(base),mkdir(base);end
oldTable=cell2table(rows,'VariableNames',columns());writetable(oldTable,fullfile(base,'existing_points.csv'));disp(oldTable);
% Fixed closer lateral locations, not a search for minimum high-order error.
offsets=[0,0;0,-.25;0,.25;0,-.35;0,.35];newLabels=["center","y_m0p25","y_p0p25","y_m0p35","y_p0p35"];
extract_directional_joint_input(dataRoot,offsets,outputName,cfg);
old=load(fullfile(oldBase,'extracted.mat'));raw=load(fullfile(base,'extracted.mat'));
assert(isequal(old.initialSpectrum,raw.initialSpectrum) && isequal(old.t,raw.t));
assert(isequal(old.eta2,raw.eta2(:,1)) && isequal(old.probe1,raw.probe1(:,1)));
assert(old.metadata.h==raw.metadata.h && old.metadata.g==raw.metadata.g && raw.metadata.Akp==Akp);
copyfile(fullfile(oldBase,'initial_convention.json'),fullfile(base,'initial_convention.json'));
newFolders=strings(5,1);newFolders(1)=string(centerFolder);ran=false(5,1);ran(1)=true;
for j=2:5
    ratio=max(abs(raw.eta2(:,j)))/referencePeak;
    fprintf('Eligibility %s: eta22 peak ratio %.8f (threshold 1/3)\n',newLabels(j),ratio);
    if ratio>=1/3
        run_directional_joint_pilot(base,char(newLabels(j)),j);
        newFolders(j)=string(fullfile(base,char(newLabels(j))));ran(j)=true;
        d=load(fullfile(newFolders(j),'joint_pilot.mat'));
        rows(end+1,:)=record(newLabels(j),d,referencePeak,referenceRange,centerMainPeak,center.report.metadata.probe); %#ok<AGROW>
    else
        rows(end+1,:)={newLabels(j),raw.metadata.probes(j,1),raw.metadata.probes(j,2),raw.metadata.Akp, ...
            max(abs(raw.eta2(:,j))),ratio,range(raw.eta2(:,j))/referenceRange,NaN,false,false,NaN,NaN,"not run: below amplitude threshold"}; %#ok<AGROW>
    end
end
summary=cell2table(rows,'VariableNames',columns());writetable(summary,fullfile(base,'summary.csv'));
criteria=struct('primary','max(abs(OW3D eta22)) over full saved record, same-x centerline baseline', ...
    'threshold',1/3,'center_peak_m',referencePeak,'center_peak_to_trough_m',referenceRange, ...
    'center_main_peak_m',centerMainPeak,'Akp',center.report.metadata.Akp, ...
    'selection_requested_after_previous_results',true,'criterion_defined_before_new_probe_errors',true, ...
    'nominal_probe_offsets_wavelengths',offsets,'actual_probe_coordinates_m',raw.metadata.probes,'kernel_tuned',false);
save(fullfile(base,'summary.mat'),'summary','criteria','newFolders','newLabels','ran','offsets');
fid=fopen(fullfile(base,'criteria.json'),'w');fprintf(fid,'%s',jsonencode(criteria));fclose(fid);
plot_directional_amplitude_gate(base);disp(summary);
end

function row=record(label,d,centerPeak,centerRange,centerMainPeak,centerXY)
r=d.d.eta2;xy=d.report.metadata.probe;
onLine=abs(xy(2)-centerXY(2))<1e-8;assert(onLine || abs(xy(1)-centerXY(1))<1e-8);
ratio=max(abs(r))/centerPeak;rangeRatio=range(r)/centerRange;
mainRatio=max(abs(r(d.focus)))/centerMainPeak;
pass=onLine || ratio>=1/3;
mt=d.metrics;sel=string(mt.window)=="main_group" & string(mt.method)=="Joint input, 7.5 deg";
eta=mt.relative_L2(sel & string(mt.variable)=="eta22");psi=mt.relative_L2(sel & string(mt.variable)=="psi22");
row={label,xy(1),xy(2),d.report.metadata.Akp,max(abs(r)),ratio,rangeRatio,mainRatio, ...
    onLine,pass,eta,psi,"evaluated"};
end

function names=columns()
names={'probe','x_m','y_m','Akp','eta22_peak_m','peak_ratio_to_focus_center', ...
    'range_ratio_to_focus_center','main_peak_ratio_to_focus_center','on_centerline', ...
    'passes_same_x_amplitude_gate','eta22_main_L2','psi22_main_L2','status'};
end
