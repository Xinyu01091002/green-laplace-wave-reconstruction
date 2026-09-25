function run_surface_potential_pilot(mf12Root)
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));addpath(root);
setup_green_laplace('MF12Root',mf12Root);set(groot,'defaultFigureVisible','off');
for akp=[.02,.12]
    source=fullfile(root,'results','unidirectional_time_series',sprintf('ow3d_boundary_kh1_alpha1_akp%03d',round(100*akp)),'pilot.mat');
    d=load(source);t=d.t;tr=t-t(1);N=numel(t);h=d.report.depth_m;kp=d.report.kp_rad_m;g=9.81;
    A=d.A;w=d.omega;k=d.k; E=exp(-1i*tr*w.');z=E*A;
    psi11=real(E*(-1i*g*A./w));
    [psi22,audit]=gl_psi22_time_pairs(A,w,k,g,h,kp,tr);psi22=real(psi22);
    mf22=mf12_psi22_time(A,k,g,h,tr);
    % OW3D psi=phi(x,eta,t); used only as held-out output, never executor input.
    raw=d.rawpsi;mask=zeros(N,1);mask(1)=1;mask(2:(N+1)/2)=2;
    H=imag(ifft(fft(raw).*mask));
    obs11=(raw(:,1)-raw(:,3)-H(:,2)+H(:,4))/4;
    obs22=(raw(:,1)-raw(:,2)+raw(:,3)-raw(:,4))/4;
    [~,i]=max(abs(z));lims=t(i)+[-2,2]*d.report.Tp_s;focus=t>=lims(1)&t<=lims(2);
    pred=[psi11,psi22,mf22];obs=[obs11,obs22,obs22];
    names=["linear psi11 from eta1","GL2+2 psi22","spectral MF12 psi22"];
    rows=cell(0,5);
    for region=1:2
        use=true(N,1);label="full";if region==2,use=focus;label="main_group";end
        for j=1:3
            v=pred(use,j);r=obs(use,j);delta=v-r;
            rows(end+1,:)={names(j),label,norm(delta)/norm(r),max(abs(delta))/max(abs(r)),norm(v)/norm(r)}; %#ok<AGROW>
        end
    end
    metrics=cell2table(rows,'VariableNames',{'method','window','relative_L2','relative_Linf','norm_ratio'});
    out=fullfile(root,'results','unidirectional_time_series',sprintf('surface_potential_boundary_alpha1_akp%03d',round(100*akp)));
    if ~isfolder(out),mkdir(out);end
    report=struct('Akp',akp,'source',source,'definition','psi=phi(x,eta,t)', ...
        'units','m^2/s','input','same declared eta1 as eta22; OW3D psi is not an input', ...
        'parent_count',numel(A),'gauge_offset_fit',false,'reference_filter',false, ...
        'gl_vs_mf12_psi22_main_L2',norm(psi22(focus)-mf22(focus))/norm(mf22(focus)), ...
        'display_limits_s',lims);
    save(fullfile(out,'fields.mat'),'report','metrics','t','psi11','psi22','mf22','obs11','obs22','focus','audit');
    writetable(metrics,fullfile(out,'metrics.csv'));
    fid=fopen(fullfile(out,'report.json'),'w');fprintf(fid,'%s',jsonencode(report));fclose(fid);
    f=figure('Color','w','Position',[100,100,1200,850]);tiledlayout(3,1,'Padding','compact','TileSpacing','compact');
    nexttile;plot(t,obs11,'k',t,psi11,'--','LineWidth',1.3);xlim(lims);grid on;ylabel('\psi_{11} (m^2/s)');
    legend('OW3D first phase sector','Linear prediction from eta1','Location','eastoutside');
    nexttile;plot(t,obs22,'k','LineWidth',1.7);hold on;plot(t,psi22,'-','Color',[0,.4,.7],'LineWidth',1.3);
    plot(t,mf22,'--','Color',[.85,.25,.1],'LineWidth',1.3);xlim(lims);grid on;ylabel('\psi_{22} (m^2/s)');
    legend('OW3D second phase sector','GL2+2','Spectral MF12 (order 2)','Location','eastoutside');
    nexttile;plot(t,psi22-obs22,'-','Color',[0,.4,.7],'LineWidth',1.3);hold on;
    plot(t,mf22-obs22,'--','Color',[.85,.25,.1],'LineWidth',1.3);xlim(lims);grid on;
    ylabel('Prediction - OW3D (m^2/s)');xlabel('Simulation elapsed time (s)');
    sgtitle(sprintf('Surface potential: k_p h = 1, Alpha = 1, Ak_p = %.2f',akp));
    exportgraphics(f,fullfile(out,'psi11_psi22_main_group.png'),'Resolution',180);
    exportgraphics(f,fullfile(out,'psi11_psi22_main_group.pdf'),'ContentType','vector');close(f);
    disp(report);disp(metrics);
end
end
