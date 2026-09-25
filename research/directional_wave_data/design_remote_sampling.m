function report=design_remote_sampling()
% Planning from already-local frozen initial spectra; no new remote field fetch.
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
out=fullfile(root,'artifacts','remote_campaign_design');if ~isfolder(out),mkdir(out);end
ids={'kh1_s25_a012','kh5_s15_a012'};percentiles=[.99,.999,.9999,.99999];entries=cell(1,numel(ids));
for j=1:numel(ids)
    source=fullfile(root,'results','directional_sweep',ids{j},'extracted.mat');d=load(source);
    m=d.metadata;r=hypot(d.kx,d.ky);f=sqrt(m.g*r.*tanh(r*m.h))/(2*pi);
    use=d.kx>0;energy=abs(d.initialSpectrum(use)).^2;frequencies=f(use);
    assert(all(isfinite(energy))&&sum(energy)>0,'Invalid initial spectrum energy.');
    [frequencies,ix]=sort(frequencies);energy=energy(ix);cdf=cumsum(energy)/sum(energy);
    quantiles=zeros(size(percentiles));for k=1:numel(percentiles),quantiles(k)=frequencies(find(cdf>=percentiles(k),1));end
    fp=sqrt(m.g*m.kp*tanh(m.kp*m.h))/(2*pi);
    entries{j}=struct('case_id',ids{j},'source',source,'peak_frequency_Hz',fp,'peak_period_s',1/fp, ...
        'energy_quantiles',percentiles,'frequency_quantiles_Hz',quantiles, ...
        'dt_for_twenty_samples_at_3fp_s',1/(20*3*fp), ...
        'samples_at_3fp_for_dt0p2',1/(.2*3*fp), ...
        'samples_per_fifth_harmonic_at_dt_0p1',1/(.1*5*quantiles(3)), ...
        'samples_per_third_harmonic_at_dt_0p2',1/(.2*3*quantiles(3)));
end
records=[entries{:}];
dur=480;headerNx=2051;headerNy=515;fieldBytes=32*headerNx*headerNy+32;
strides=[4,.4,.2,.1,.05,2,5];rows=zeros(numel(strides),3);
for j=1:numel(strides)
    count=round(dur/strides(j))+1;rows(j,:)=[strides(j),count,4*count*fieldBytes/2^30];
end
storage=array2table(rows,'VariableNames',{'output_dt_s','samples_per_phase','four_phase_EP_GiB'});
% Native Cartesian format 20: three surface and thirteen Nz-volume records.
% Nz=9 physical nodes plus one bottom ghost. All sixteen records have markers.
nz=10;names=["full_x_nine_rows";"local_81_by_9";"center_x_nine_probes"];
pointCounts=[2049*9;81*9;9];nt=1201;
kinBytes=4*((40*pointCounts+72+8*nz)+nt*(8*(3+13*nz)*pointCounts+16*8));
kinStorage=table(names,pointCounts,kinBytes/2^30, ...
    'VariableNames',{'layout','horizontal_points','four_phase_120_to_360_dt0p2_GiB'});
report=struct('initial_spectra',records,'strong_EP_bytes',fieldBytes,'duration_s',dur, ...
    'old_output_Nyquist_Hz',1/(2*4), ...
    'four_fine_phases_plus_one_coarse_phase_GiB',(4*4801+2401)*fieldBytes/2^30, ...
    'kinematics_window_s',[120,360],'kinematics_output_dt_s',.2, ...
    'kinematics_storage',table2struct(kinStorage), ...
    'note',['EP alternative stores coordinates, eta and surface potential. ' ...
    'Native local kinematics also writes volume phi, velocities and gradients. ' ...
    'Initial positive-kx spectral energy quantiles are planning proxies, not hard cutoffs. ' ...
    'Uncompressed size estimates exclude checkpoints, logs and derived data; no runtime estimate.']);
writetable(storage,fullfile(out,'storage_estimates.csv'));
writetable(kinStorage,fullfile(out,'local_kinematics_storage.csv'));
fid=fopen(fullfile(out,'sampling_report.json'),'w');fprintf(fid,'%s',jsonencode(report));fclose(fid);
disp(records);disp(storage);disp(kinStorage);disp(report);
end
