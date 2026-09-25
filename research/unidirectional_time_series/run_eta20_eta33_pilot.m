function run_eta20_eta33_pilot(mf12Root)
% Independent nonzero eta20 / positive eta33 trials. MF12 is order TWO only.
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));addpath(root);
setup_green_laplace('MF12Root',mf12Root);
for akp=[.02,.12]
    source=fullfile(root,'results','unidirectional_time_series', ...
        sprintf('ow3d_boundary_kh1_alpha1_akp%03d',round(100*akp)),'pilot.mat');
    d=load(source); h=d.report.depth_m; kp=d.report.kp_rad_m; g=9.81;
    t=d.t;trel=t-t(1);N=numel(t);dt=mean(diff(t));Tp=d.report.Tp_s;
    F=fft(d.eta1phase)/N;allbins=(1:floor((N-1)/2))';
    allw=2*pi*allbins/(N*dt); allk=zeros(size(allw));
    for j=1:numel(allw)
        allk(j)=fzero(@(k)g*k*tanh(k*h)-allw(j)^2,[0,max(1,2*allw(j)^2/g+2/h)]);
    end
    keep=allk*h>0.5 & 3*allw<pi/dt;
    bins=allbins(keep); w=allw(keep); k=allk(keep); A=2*conj(F(bins+1));
    z=exp(-1i*trel*w.')*A; eta1=real(z);
    raw=d.raw; H=imag(analytic(raw));
    obs20=mean(raw,2);
    obs33=(raw(:,1)-raw(:,3)+H(:,2)-H(:,4))/4;
    frequencies=2*pi*[0:floor(N/2),-floor(N/2):-1]'/(N*dt);
    wp=sqrt(g*kp*tanh(kp*h));
    lowpass=abs(frequencies)>0 & abs(frequencies)<.5*wp;
    names20=["GL6 diagnostic","GL12 diagnostic","GL16 diagnostic","Spectral MF12 (order 2)","Walker Eq.14"];
    pred20=zeros(N,5);audit20=cell(1,3);
    for j=1:3
        ranks=[6,12,16];[pred20(:,j),audit20{j}]=gl_eta20_time_pairs(A,w,k,h,trel,ranks(j));
    end
    pred20(:,4)=mf12_eta20_time(A,k,g,h,trel);
    pred20(:,5)=-kp/(2*sinh(2*kp*h))*abs(z).^2;
    raw20=pred20; obs20raw=obs20;
    pred20=real(ifft(fft(pred20).*lowpass));obs20=real(ifft(fft(obs20).*lowpass));
    names33=["GL4","GL6","GL8","VWA","Walker"];
    pred33=zeros(N,5);audit33=cell(1,3);
    for j=1:3
        ranks=[4,6,8];[v,audit33{j}]=gl_eta33_time_triples(A,w,k,h,kp,trel,ranks(j),bins);
        pred33(:,j)=real(v);
    end
    sigma=tanh(k*h); B33=(27-9*sigma.^2+9*sigma.^4-3*sigma.^6)./(64*sigma.^6);
    E=exp(-1i*trel*w.');
    % Preserve the established VWA time-series product of two first-order multipliers.
    pred33(:,4)=real(z.*(E*(A.*B33.*k)).*(E*(A.*k)));
    cp=sech(2*kp*h); bp=3*(1+3*cp+3*cp^2+2*cp^3)/(8*(1-cp)^3);
    pred33(:,5)=kp^2*bp*real(z.^3);
    assert(all(isfinite([pred20,pred33]),'all'));
    [~,index]=max(abs(z)); displayLimits=t(index)+[-2,2]*Tp;
    focus=t>=displayLimits(1)&t<=displayLimits(2);
    out=fullfile(root,'results','unidirectional_time_series',sprintf('eta20_eta33_boundary_alpha1_akp%03d',round(100*akp)));
    if ~isfolder(out),mkdir(out);end
    metrics20=metrics(names20,pred20,obs20,focus);metrics33=metrics(names33,pred33,obs33,focus);
    writetable(metrics20,fullfile(out,'eta20_metrics.csv'));writetable(metrics33,fullfile(out,'eta33_metrics.csv'));
    report=struct('Akp',akp,'kph',kp*h,'alpha',1,'h',h,'Tp',Tp, ...
        'parent_count',numel(A),'minimum_parent_kh',min(k*h), ...
        'retained_positive_frequency_energy',sum(abs(F(bins+1)).^2)/sum(abs(F(allbins+1)).^2), ...
        'input_projection_relative_L2',norm(eta1-d.eta1phase)/norm(d.eta1phase), ...
        'eta20_nonzero_cutoff_rad_s',.5*wp,'eta20_lowpass_bins',nnz(lowpass), ...
        'eta20_excluded_observed_dc_m',mean(obs20raw),'source',source, ...
        'eta33_mf12_computed',false,'display_limits_s',displayLimits);
    save(fullfile(out,'fields.mat'),'report','t','eta1','pred20','pred33','obs20','obs33', ...
        'raw20','obs20raw','names20','names33','A','w','k','bins','focus','audit20','audit33');
    fid=fopen(fullfile(out,'report.json'),'w');assert(fid>=0);fprintf(fid,'%s',jsonencode(report));fclose(fid);
    plot_eta20_eta33(out);disp(report);disp(metrics20);disp(metrics33);
end
end

function a=analytic(x)
N=size(x,1);mask=zeros(N,1);mask(1)=1;mask(2:(N+1)/2)=2;a=ifft(fft(x).*mask);
end

function tab=metrics(names,pred,obs,focus)
rows=cell(0,5);
for region=1:2
    mask=true(size(focus));label="full";
    if region==2,mask=focus;label="main_group";end
    r=obs(mask);
    for j=1:numel(names)
        v=pred(mask,j); delta=v-r;
        rows(end+1,:)={names(j),label,norm(delta)/norm(r),max(abs(delta))/max(abs(r)),norm(v)/norm(r)}; %#ok<AGROW>
    end
end
tab=cell2table(rows,'VariableNames',{'method','window','relative_L2','relative_Linf','norm_ratio'});
end
