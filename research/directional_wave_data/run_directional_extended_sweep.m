function run_directional_extended_sweep(dataRoot)
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));out=fullfile(root,'results','directional_sweep');
if ~isfolder(out),mkdir(out);end
% Prespecified cases; no selection from high-order scores.
specs={ 'kh1_s25_a002','test1',1,25,.02,900; ...
        'kh2_s25_a002','test1',2,25,.02,900; ...
        'kh5_s25_a002','test1',5,25,.02,900; ...
        'kh1_s5_a002','test1',1,5,.02,900; ...
        'kh1_s25_a012','test6',1,25,.12,1200; ...
        'kh5_s15_a002','test1',5,15,.02,900; ...
        'kh5_s15_a012','test6',5,15,.12,1200};
rows=cell(0,10);diffrows=cell(0,9);failures=cell(0,3);folders={};
for ic=1:size(specs,1)
    id=specs{ic,1};base=fullfile(out,id);
    cfg=struct('kph',specs{ic,3},'spread',specs{ic,4},'Akp',specs{ic,5},'lastStep',specs{ic,6});
    offsets=[0,0];labels="center";
    if ic==1,offsets=[0,0;-3,0;3,0;0,-1;0,1;0,-1.5;0,1.5];labels=["center","x_m3","x_p3","y_m1","y_p1","y_m1p5","y_p1p5"];end
    fprintf('EXTRACT FAMILY %s\n',id);
    extract_directional_joint_input(fullfile(dataRoot,specs{ic,2}),offsets,fullfile('directional_sweep',id),cfg);
    d=load(fullfile(base,'extracted.mat'));m=d.metadata;
    r=hypot(d.kx,d.ky);w=sqrt(m.g*r.*tanh(r*m.h));use=d.kx>0&r*m.h>=.3&w<pi/m.dt;
    e=d.initialSpectrum(use);p=d.initialPsiSpectrum(use);factor=m.g./w(use);
    minus=norm(p+1i*factor.*e)/norm(p);plus=norm(p-1i*factor.*e)/norm(p);
    convention=struct('minus_time_polarization_relative_error',minus,'plus_time_polarization_relative_error',plus, ...
        'stored_positive_k_has_positive_time',plus<minus,'basis','initial fields only');
    fid=fopen(fullfile(base,'initial_convention.json'),'w');fprintf(fid,'%s',jsonencode(convention));fclose(fid);
    if min(plus,minus)>=.05 || max(plus,minus)<=1.9
        failures(end+1,:)={id,'all','Initial propagation convention unresolved'};continue; %#ok<AGROW>
    end
    for ip=1:size(offsets,1)
        try
            fprintf('RUN %s / %s\n',id,labels(ip));
            run_directional_joint_pilot(base,char(labels(ip)),ip,true);
            folder=fullfile(base,char(labels(ip)));folders{end+1,1}=folder; %#ok<AGROW>
            f=load(fullfile(folder,'joint_pilot.mat'));mt=f.metrics;rr=f.report;
            for variable=["eta22","psi22"]
                ix=string(mt.variable)==variable&string(mt.method)=="Joint input, 7.5 deg"&string(mt.window)=="main_group";
                rows(end+1,:)={id,labels(ip),m.kph,m.spread_label_degrees,m.Akp,variable,mt.relative_L2(ix),rr.conditioning{2}.observed_energy_weighted,rr.main_window_complete,rr.initial_frequency_std_over_mean}; %#ok<AGROW>
            end
            a=readtable(fullfile(folder,'eta20_metrics.csv'),'TextType','string');
            for band=[.5,3]
                ix=a.rank==16&a.cutoff_ratio==band&a.window=="main_group";
                diffrows(end+1,:)={id,labels(ip),m.kph,m.spread_label_degrees,m.Akp,band,a.relative_L2(ix),a.nonzero_bins(ix),rr.main_window_complete}; %#ok<AGROW>
            end
        catch exception
            failures(end+1,:)={id,char(labels(ip)),exception.message}; %#ok<AGROW>
            fprintf(2,'FAILED %s / %s: %s\n',id,labels(ip),exception.message);
        end
        save(fullfile(out,'progress.mat'),'rows','diffrows','failures','folders','specs');
    end
end
summary=cell2table(rows,'VariableNames',{'case_id','probe','kph','spread_deg','Akp','variable','main_L2','weighted_condition','complete_main_window','frequency_std_over_mean'});
eta20summary=cell2table(diffrows,'VariableNames',{'case_id','probe','kph','spread_deg','Akp','cutoff_ratio','main_L2','retained_bins','complete_main_window'});
failureTable=cell2table(failures,'VariableNames',{'case_id','probe','reason'});
writetable(summary,fullfile(out,'summary.csv'));writetable(eta20summary,fullfile(out,'eta20_summary.csv'));writetable(failureTable,fullfile(out,'failures.csv'));
save(fullfile(out,'summary.mat'),'summary','eta20summary','failureTable','folders','specs');disp(summary);disp(eta20summary);disp(failureTable);
plot_directional_extended_sweep(out);
end
