function plot_random_space(root)
run=fullfile(root,'akp012');out=fullfile(root,'spatial-view-v3');assert(~isfolder(out));mkdir(out);
s=load(fullfile(run,'inputs','initial_fields.mat'));d=s.report.design;
cfg=jsondecode(fileread(fullfile(run,'settings.json')));lambda=2*pi/d.kp;
[ny,nx,~]=size(s.E);x=(0:nx-1)*d.domain_m(1)/nx;y=(0:ny-1)*d.domain_m(2)/ny;
eta=s.E(:,:,1);linear=real(s.eta1);assert(all(isfinite(eta),'all'));
% Audit the actual random input against its recorded focused parent and seed.
parentPath=s.report.randomization.source;if ~startsWith(parentPath,'/'),parentPath=fullfile('/home/lxy/green-laplace-unidirectional-time-series',parentPath);end;parent=load(parentPath,'C');rng(s.report.randomization.seed,'twister');increments=2*pi*rand(size(parent.C));
phaseParity=norm(s.C-parent.C.*exp(1i*increments))/norm(s.C);assert(phaseParity<1e-12);
energy=abs(s.C(:)).^2;om=s.om(:);meanOm=sum(energy.*om)/sum(energy);sigmaOm=sqrt(sum(energy.*(om-meanOm).^2)/sum(energy));
report=struct('field','Actual saved HOS initial elevation, phase0, t=0; MF12 order2 eta11+eta20+eta22', ...
 'evolved_spatial_fields_available',false,'seed',s.report.randomization.seed,'random_phase_reproduction_relative_L2',phaseParity, ...
 'modal_amplitude_parity',norm(abs(s.C)-abs(parent.C))/norm(abs(parent.C)), ...
 'Hs_linear_spatial',4*std(linear(:),1),'eta_min',min(eta,[],'all'),'eta_max',max(eta,[],'all'), ...
 'frequency_mean_Hz',meanOm/(2*pi),'frequency_std_Hz',sigmaOm/(2*pi),'relative_frequency_bandwidth',sigmaOm/meanOm, ...
 'record_duration_s',cfg.duration_s,'record_duration_Tp',cfg.duration_s/cfg.Tp, ...
 'spectrum_note','Random phases preserve modal amplitudes. Variance weights are abs(C)^2, not the original amplitude weights. No spatial taper or refocusing.');
f=figure('Visible','off','Position',[40 40 1450 1000]);tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
nexttile;imagesc(x/lambda,y/lambda,eta);set(gca,'YDir','normal');axis image;colorbar;clim([-1 1]*max(abs(eta),[],'all'));colormap(turbo);hold on;
plot(cfg.probes_xy_m(:,1)/lambda,cfg.probes_xy_m(:,2)/lambda,'wo','MarkerFaceColor','k','MarkerSize',5);
xlabel('x / lambda_p');ylabel('y / lambda_p');title('Actual initial random surface eta(x,y,0), phase 0 (m); circles: probe locations');
nexttile;idx=cfg.probe_indices(1,2);plot(x/lambda,eta(idx,:),'k',x/lambda,linear(idx,:),'r--','LineWidth',.9);grid on;xlim([0 50]);xlabel('x / lambda_p');ylabel('Elevation (m)');title('Initial along-propagation section at y = 10 lambda_p');legend('Total initial elevation','First-order elevation','Location','best');
nexttile;hold on;for j=[1 2 4],iy=cfg.probe_indices(j,2);plot(x/lambda,eta(iy,:),'LineWidth',1);end;grid on;xlim([20 30]);xlabel('x / lambda_p');ylabel('Elevation (m)');title('Same initial field: magnified sections near the probes');legend('y offset 0 m','y offset -52.78 m','y offset -70.38 m','Location','best');
exportgraphics(f,fullfile(out,'random_initial_space.png'),'Resolution',160);exportgraphics(f,fullfile(out,'random_initial_space.pdf'),'ContentType','vector');close(f);
a=load(fullfile(run,'random-gl-comparison-v1','comparison.mat'),'t','raw','first','reference');
f=figure('Visible','off','Position',[40 40 1400 900]);tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
nexttile;plot(a.t,a.raw(:,1,1),'k');hold on;plot(a.t,a.first(:,1),'r--');grid on;xlim([0 220]);ylabel('Elevation (m)');title('Center probe: complete HOS record, not only the scoring window');legend('Total eta, phase0','Separated first harmonic','Location','northwest','FontSize',10);
nexttile;plot(a.t,a.first(:,1),'b');hold on;env=abs(hilbert(a.first(:,1)));plot(a.t,env,'k--',a.t,-env,'k--');grid on;xlim([0 220]);ylabel('Elevation (m)');title('First harmonic and its finite-record envelope');
nexttile;plot(a.t,a.reference(:,1),'k');grid on;xlim([0 220]);ylabel('Second harmonic (m)');xlabel('Time (s)');title('Four-phase second-harmonic sector: this small component was shown in the GL comparison');
exportgraphics(f,fullfile(out,'random_full_timeseries.png'),'Resolution',160);close(f);
f=fopen(fullfile(out,'report.json'),'w');fprintf(f,'%s',jsonencode(report,PrettyPrint=true));fclose(f);disp(report);
end
