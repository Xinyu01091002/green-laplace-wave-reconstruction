function run_directional_joint_pilot()
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));addpath(root);setup_green_laplace;
addpath(fullfile(root,'research','unidirectional_time_series'));
base=fullfile(root,'results','directional_joint_input');
d=load(fullfile(base,'extracted.mat'));m=d.metadata;t=d.t;tr=t-t(1);N=numel(t);dt=m.dt;
convention=jsondecode(fileread(fullfile(base,'initial_convention.json')));
out=fullfile(base,'verified_convention');if ~isfolder(out),mkdir(out);end
g=m.g;h=m.h;kp=m.kp;wp=sqrt(g*kp*tanh(kp*h));Tp=2*pi/wp;
realInput=real(d.probe1);F=fft(realInput)/N;allbins=(1:floor((N-1)/2))';
allw=2*pi*allbins/(N*dt);allk=zeros(size(allw));
for j=1:numel(allw),allk(j)=fzero(@(k)g*k*tanh(k*h)-allw(j)^2,[0,1]);end
keep=allk*h>=.3;bins=allbins(keep);w=allw(keep);k=allk(keep);observed=2*conj(F(bins+1));
input=real(exp(-1i*tr*w.')*observed);
radial=hypot(d.kx,d.ky);frequency=sqrt(g*radial.*tanh(radial*h));
initialEligible=d.kx>0 & h*radial>=.3 & frequency<pi/dt;
initialMass=abs(d.initialSpectrum).^2;
initialRetained=sum(initialMass(initialEligible))/sum(initialMass,'all');
ids=find(initialEligible);kx0=d.kx(ids);ky0=d.ky(ids);w0=frequency(ids);
A0=d.initialSpectrum(ids).*exp(1i*(kx0*m.probe(1)+ky0*m.probe(2)));
if convention.stored_positive_k_has_positive_time
    % Conjugate the stored positive-k/+omega analytic representation to the
    % physical negative-k/-omega representation; rotate both axes by pi so
    % the released strict-forward kernel retains positive kx. Dot products
    % and sum magnitudes are invariant; no target-selected sign is applied.
    A0=conj(A0);
end
angles=atan2(ky0,kx0)*180/pi;
eta=zeros(N,4);psi=eta;allocation=cell(1,2);condition=cell(1,2);predictedInput=zeros(N,2);
angleWidths=[15,7.5];
for resolution=1:2
    width=angleWidths(resolution);centers=(-90+width/2:width:90-width/2)';
    group=1+floor((angles+90)/width);assert(all(group>=1 & group<=numel(centers)));
    directionalPriorTime=complex(zeros(N,numel(centers)));
    for direction=1:numel(centers)
        active=find(group==direction);
        for start=1:256:numel(active)
            ix=active(start:min(start+255,numel(active)));
            directionalPriorTime(:,direction)=directionalPriorTime(:,direction)+exp(-1i*tr*w0(ix).')*A0(ix);
        end
    end
    % Negative temporal exponent convention: ifft returns the complex coefficients.
    prior=ifft(directionalPriorTime,[],1);prior=prior(bins+1,:);
    [joint,condition{resolution}]=allocate_directional_record(prior,observed);
    predictedInput(:,resolution)=real(exp(-1i*tr*w.')*sum(prior,2));
    recovered=real(exp(-1i*tr*w.')*sum(joint,2));
    assert(norm(recovered-input)/norm(input)<1e-11);
    [KX,TH]=ndgrid(k,deg2rad(centers));KY=KX.*sin(TH);KX=KX.*cos(TH);
    OM=repmat(w,1,numel(centers));
    fprintf('Directional width %.1f deg, parents %d, max conditioning %.4g\n',width,numel(joint),condition{resolution}.max);
    [e,p]=gl_directional_sum_time(joint(:),OM(:),KX(:),KY(:),g,h,kp,tr,8);
    eta(:,resolution)=real(e);psi(:,resolution)=real(p);
    if resolution==1
        [e,p]=gl_directional_sum_time(prior(:),OM(:),KX(:),KY(:),g,h,kp,tr,8);
        eta(:,3)=real(e);psi(:,3)=real(p);
    end
    allocation{resolution}=struct('angles_degrees',centers,'prior',prior,'joint',joint, ...
        'kx',KX,'ky',KY,'omega',OM);
end
% Declared direction-blind control; not selected to improve an error metric.
[e,p]=gl_directional_sum_time(observed,w,k,zeros(size(k)),g,h,kp,tr,8);
eta(:,4)=real(e);psi(:,4)=real(p);
names=["Joint input, 15 deg","Joint input, 7.5 deg","Initial spectrum only","Eta1 only, all directions zero"];
[~,peak]=max(abs(exp(-1i*tr*w.')*observed));limits=t(peak)+[-2,2]*Tp;
focus=t>=limits(1)&t<=limits(2);
rows=cell(0,6);
for variable=1:2
    ref=d.eta2;pred=eta;varname="eta22";
    if variable==2,ref=d.psi2;pred=psi;varname="psi22";end
    for region=1:2
        use=true(N,1);label="full";if region==2,use=focus;label="main_group";end
        for j=1:4
            r=ref(use);v=pred(use,j);
            rows(end+1,:)={varname,names(j),label,norm(v-r)/norm(r),max(abs(v-r))/max(abs(r)),norm(v)/norm(r)}; %#ok<AGROW>
        end
    end
end
metrics=cell2table(rows,'VariableNames',{'variable','method','window','relative_L2','relative_Linf','norm_ratio'});
writetable(metrics,fullfile(out,'metrics.csv'));
report=struct('metadata',m,'initial_convention',convention, ...
    'angle_widths_degrees',angleWidths,'temporal_bins',bins.', ...
    'initial_spectrum_domain_retained_energy',initialRetained, ...
    'temporal_input_retained_energy',sum(abs(F(bins+1)).^2)/sum(abs(F(allbins+1)).^2), ...
    'input_projection_relative_L2',norm(input-realInput)/norm(realInput), ...
    'linear_prior_probe_relative_L2',norm(predictedInput(:,1)-input)/norm(input), ...
    'angular_refinement_eta22_relative_L2',norm(eta(focus,2)-eta(focus,1))/norm(eta(focus,2)), ...
    'angular_refinement_psi22_relative_L2',norm(psi(focus,2)-psi(focus,1))/norm(psi(focus,2)), ...
    'conditioning',{condition},'display_limits_s',limits, ...
    'reference_interpolated',false,'output_sums_above_nyquist','evaluated at saved times, not claimed temporally resolved', ...
    'initial_first_order_exact',false,'high_order_reference_used_for_allocation',false);
save(fullfile(out,'joint_pilot.mat'),'report','metrics','allocation','t','realInput','input', ...
    'predictedInput','eta','psi','names','d','focus','-v7.3');
fid=fopen(fullfile(out,'report.json'),'w');fprintf(fid,'%s',jsonencode(report));fclose(fid);
plot_directional_joint_pilot(out);disp(report);disp(metrics);
end
