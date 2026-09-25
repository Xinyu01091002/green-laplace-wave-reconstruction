function reproject_eta20_band(cutoffRatio)
% User-requested output-band experiment; inputs and kernels are unchanged.
if nargin<1,cutoffRatio=3;end
validateattributes(cutoffRatio,{'double'},{'scalar','positive','finite'});
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
for akp=[.02,.12]
    prior=fullfile(root,'results','unidirectional_time_series', ...
        sprintf('eta20_eta33_boundary_alpha1_akp%03d',round(100*akp)));
    source=fullfile(prior,'fields.mat'); d=load(source);
    out=[prior,'_eta20_band',strrep(sprintf('%g',cutoffRatio),'.','p')];
    if ~isfolder(out),mkdir(out);end
    N=numel(d.t);dt=mean(diff(d.t));
    omega=2*pi*[0:ceil(N/2)-1,-floor(N/2):-1]'/(N*dt);
    wp=d.report.eta20_nonzero_cutoff_rad_s/.5;
    assert(cutoffRatio*wp<pi/dt,'Requested band exceeds temporal Nyquist.');
    band=abs(omega)>0 & abs(omega)<cutoffRatio*wp;
    d.pred20=real(ifft(fft(d.raw20).*band));
    d.obs20=real(ifft(fft(d.obs20raw).*band));
    assert(all(isfinite([d.pred20,d.obs20]),'all'));
    d.report.eta20_nonzero_cutoff_rad_s=cutoffRatio*wp;
    d.report.eta20_cutoff_ratio=cutoffRatio;
    d.report.eta20_lowpass_bins=nnz(band);
    d.report.projection_source=source;
    d.report.recomputed_kernels=false;
    d.report.parent_input_unchanged=true;
    rows=cell(0,5);
    for region=1:2
        mask=true(N,1);label="full";
        if region==2,mask=d.focus;label="main_group";end
        r=d.obs20(mask);
        for j=1:numel(d.names20)
            v=d.pred20(mask,j);e=v-r;
            rows(end+1,:)={d.names20(j),label,norm(e)/norm(r),max(abs(e))/max(abs(r)),norm(v)/norm(r)}; %#ok<AGROW>
        end
    end
    metrics=cell2table(rows,'VariableNames',{'method','window','relative_L2','relative_Linf','norm_ratio'});
    writetable(metrics,fullfile(out,'eta20_metrics.csv'));
    save(fullfile(out,'fields.mat'),'-struct','d');
    fid=fopen(fullfile(out,'report.json'),'w');assert(fid>=0);fprintf(fid,'%s',jsonencode(d.report));fclose(fid);
    plot_eta20_eta33(out,20);
    old=readtable(fullfile(prior,'eta20_metrics.csv'),'TextType','string');
    baseline=old(old.window=="main_group",:);current=metrics(metrics.window=="main_group",:);
    assert(isequal(string(baseline.method),string(current.method)));
    comparison=table(current.method,baseline.relative_L2,current.relative_L2, ...
        'VariableNames',{'method','relative_L2_band0p5','relative_L2_new_band'});
    writetable(comparison,fullfile(out,'cutoff_comparison.csv'));
    fprintf('Akp %.2f: %.1f omega_p, %d nonzero bins\n',akp,cutoffRatio,nnz(band));disp(comparison);
end
end
