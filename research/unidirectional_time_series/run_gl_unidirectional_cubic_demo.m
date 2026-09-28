function run_gl_unidirectional_cubic_demo
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));addpath(fileparts(mfilename('fullpath')));
out=fullfile(root,'artifacts','unidirectional_time_series','unidirectional-gl-cubic-demo');
assert(~isfolder(out));mkdir(out);profile clear;profile on;
reports=cell(2,1);traces=cell(2,1);
for icase=1:2
    akp=[.02,.12];source=fullfile(root,'results','unidirectional_time_series',sprintf('ow3d_boundary_kh1_alpha1_akp%03d',round(100*akp(icase))),'pilot.mat');
    d=load(source,'t','raw','report');g=9.81;h=d.report.depth_m;kp=d.report.kp_rad_m;Tp=d.report.Tp_s;
    N=numel(d.t);mask=zeros(N,1);mask(1)=1;mask(2:(N+1)/2)=2;
    ht=imag(ifft(fft(d.raw).*mask));eta1=(d.raw(:,1)-d.raw(:,3)-ht(:,2)+ht(:,4))/4;
    observed2=(d.raw(:,1)-d.raw(:,2)+d.raw(:,3)-d.raw(:,4))/4;
    observed3=(d.raw(:,1)-d.raw(:,3)+ht(:,2)-ht(:,4))/4;
    opts=struct('order',3,'omega_max',4*sqrt(g*kp*tanh(kp*h)),'peak_wavenumber',kp, ...
        'quadrature_rank',8,'memory_budget_MiB',2048);
    timer=tic;[result,audit]=gl_unidirectional_time_series(eta1,d.t,g,h,opts);seconds=toc(timer);
    F=fft(result.eta1_used);analyticMask=zeros(N,1);analyticMask(1)=1;analyticMask(2:(N+1)/2)=2;
    [~,peak]=max(abs(ifft(F.*analyticMask)));main=abs(d.t-d.t(peak))<=2*Tp;
    reports{icase}=struct('akp',akp(icase),'audit',audit,'profiled_seconds',seconds, ...
        'eta22_main_relative',norm(result.eta22(main)-observed2(main))/norm(observed2(main)), ...
        'eta33_main_relative',norm(result.eta33(main)-observed3(main))/norm(observed3(main)), ...
        'eta33_full_relative',norm(result.eta33-observed3)/norm(observed3));
    traces{icase}=struct('result',result,'reference2',observed2,'reference3',observed3,'main',main,'Tp',Tp);
    fprintf('Akp %.2f: grid change %.6g, converged %d; eta22 %.6g eta33 %.6g\n',akp(icase),audit.levels(end).relative_change,audit.converged,reports{icase}.eta22_main_relative,reports{icase}.eta33_main_relative);
end
profile off;p=profile('info');names=string({p.FunctionTable.FunctionName});
assert(~any(contains(names,["time_pairs","time_triples","ordered_pair","pair_reference","integral_reference"])));
save(fullfile(out,'demo.mat'),'reports','traces');fid=fopen(fullfile(out,'reports.json'),'w');fprintf(fid,'%s\n',jsonencode(reports));fclose(fid);
fid=fopen(fullfile(out,'called_functions.json'),'w');fprintf(fid,'%s\n',jsonencode(names));fclose(fid);
fig=figure('Visible','off','Color','w','Position',[100,100,1100,650]);tiledlayout(2,2,'TileSpacing','compact');
for j=1:2
    z=traces{j};tx=(z.result.t-z.result.t(1))/z.Tp;
    nexttile(j);plot(tx,z.reference2,'k-',tx,z.result.eta22,'--');grid on;xlim([min(tx(z.main)),max(tx(z.main))]);title(sprintf('Akp %.2f, eta22',reports{j}.akp));ylabel('Elevation (m)');
    nexttile(j+2);plot(tx,z.reference3,'k-',tx,z.result.eta33,'--');grid on;xlim([min(tx(z.main)),max(tx(z.main))]);title('eta33');xlabel('Time from record start / Tp');ylabel('Elevation (m)');legend('OW3D phase sector','GL8 operator','Location','best');
end
exportgraphics(fig,fullfile(out,'demo.png'),'Resolution',150);exportgraphics(fig,fullfile(out,'demo.pdf'),'ContentType','vector');close(fig);
end
