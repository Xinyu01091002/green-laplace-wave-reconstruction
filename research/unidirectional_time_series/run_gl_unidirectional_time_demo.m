function run_gl_unidirectional_time_demo(refined)
% Use measured OW3D phase records only; never load GL reference predictions.
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));addpath(fileparts(mfilename('fullpath')));
if nargin<1,refined=false;end
if refined,tag='unidirectional-gl-demo-refined';else,tag='unidirectional-gl-demo';end
out=fullfile(root,'artifacts','unidirectional_time_series',tag);
assert(~isfolder(out));mkdir(out);profile clear;profile on;
rows=cell(0,9);reports=cell(2,1);traces=cell(2,1);
for icase=1:2
    amplitudes=[.02,.12];akp=amplitudes(icase);
    source=fullfile(root,'results','unidirectional_time_series',sprintf('ow3d_boundary_kh1_alpha1_akp%03d',round(100*akp)),'pilot.mat');
    d=load(source,'t','raw','report');g=9.81;h=d.report.depth_m;kp=d.report.kp_rad_m;Tp=d.report.Tp_s;
    ht=imag(analytic_signal(d.raw));eta1=(d.raw(:,1)-d.raw(:,3)-ht(:,2)+ht(:,4))/4;
    observed=(d.raw(:,1)-d.raw(:,2)+d.raw(:,3)-d.raw(:,4))/4;
    opts=struct('omega_max',4*sqrt(g*kp*tanh(kp*h)),'peak_wavenumber',kp,'quadrature_rank',8, ...
        'domain_lengths',[25 50 100 200 400],'relative_tolerance',.005);
    if refined,opts.domain_lengths=[25 50 100 200 400 800];opts.relative_tolerance=.0005;end
    start=tic;[result,audit]=gl_unidirectional_time_series(eta1,d.t,g,h,opts);total=toc(start);
    [~,peak]=max(abs(analytic_signal(result.eta1_used)));main=abs(d.t-d.t(peak))<=2*Tp;
    errMain=norm(result.eta22(main)-observed(main))/norm(observed(main));
    errFull=norm(result.eta22-observed)/norm(observed);
    rows(end+1,:)={akp,numel(eta1),numel(audit.parent_bins),audit.levels(end).L_over_h, ...
        audit.levels(end).Nx,audit.converged,audit.levels(end).relative_change,errMain,total}; %#ok<AGROW>
    reports{icase}=struct('source',source,'audit',audit,'main_relative',errMain,'full_relative',errFull,'total_seconds',total);
    traces{icase}=struct('result',result,'observed',observed,'main',main,'Tp',Tp);
    fprintf('Akp %.2f: converged %d, grid change %.6g, observed main %.6g, full %.6g, total %.3fs\n', ...
        akp,audit.converged,audit.levels(end).relative_change,errMain,errFull,total);
end
profile off;p=profile('info');names=string({p.FunctionTable.FunctionName});
assert(~any(contains(names,["time_pairs","time_triples","ordered_pair","pair_reference","integral_reference"])));
metrics=cell2table(rows,'VariableNames',{'akp','samples','parent_bins','L_over_h','Nx','converged','grid_change','OW3D_main_relative','total_seconds'});
writetable(metrics,fullfile(out,'metrics.csv'));save(fullfile(out,'demo.mat'),'reports','traces','metrics');
fid=fopen(fullfile(out,'reports.json'),'w');fprintf(fid,'%s\n',jsonencode(reports));fclose(fid);
fid=fopen(fullfile(out,'called_functions.json'),'w');fprintf(fid,'%s\n',jsonencode(names));fclose(fid);
fig=figure('Visible','off','Color','w','Position',[100,100,1200,700]);tiledlayout(2,2,'TileSpacing','compact');
for icase=1:2
    z=traces{icase};tx=(z.result.t-z.result.t(1))/z.Tp;
    nexttile(icase);plot(tx,z.observed,'k-','DisplayName','OW3D second phase sector');hold on;
    plot(tx,z.result.eta22,'--','DisplayName','GL8 operator time series');xlim([min(tx(z.main)),max(tx(z.main))]);
    grid on;legend('Location','best');ylabel('Elevation (m)');title(sprintf('Akp %.2f',metrics.akp(icase)));
    nexttile(icase+2);lev=reports{icase}.audit.levels;
    semilogy([lev.L_over_h],[lev.relative_change],'-o');grid on;xlabel('Auxiliary length L/h');ylabel('Successive-grid relative change');
end
exportgraphics(fig,fullfile(out,'demo.png'),'Resolution',150);exportgraphics(fig,fullfile(out,'demo.pdf'),'ContentType','vector');close(fig);
disp(metrics);
end
function z=analytic_signal(y)
N=size(y,1);mask=zeros(N,1);mask(1)=1;
if mod(N,2),mask(2:(N+1)/2)=2;else,mask(2:N/2)=2;mask(N/2+1)=1;end
z=ifft(fft(y).*mask);
end
