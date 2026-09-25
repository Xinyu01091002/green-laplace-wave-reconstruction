function run_directional_bandwidth_benchmark(mf12Root)
% Controlled directional synthetic test, not OW3D evidence or nonlinear evolution.
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));addpath(root);setup_green_laplace('MF12Root',mf12Root);
addpath(fullfile(root,'research','unidirectional_time_series'));
out=fullfile(root,'results','directional_bandwidth');if ~isfolder(out),mkdir(out);end
g=9.81;kp=.0279;N=256;bin=(8:32)';angles=deg2rad(-60:20:60);rows=cell(0,8);cases=cell(2,3);
for depthIndex=1:2
    depths=[1,5];kph=depths(depthIndex);h=kph/kp;wp=sqrt(g*kp*tanh(kph));
    dt=(20*2*pi/wp)/N;t=(0:N-1)'*dt;w=bin*wp/20;k=zeros(size(w));
    for j=1:numel(w),k(j)=fzero(@(x)g*x*tanh(x*h)-w(j)^2,[0,1]);end
    assert(all(k*h>=.3));[K,TH]=ndgrid(k,angles);kx=K.*cos(TH);ky=K.*sin(TH);
    W=repmat(w,1,numel(angles));ib=repmat(bin,1,numel(angles));
    for ibeta=1:3
        betas=[.08,.16,.32];beta=betas(ibeta);
        radial=exp(-((w/wp-1).^2)/(4*beta^2));direction=exp(-.5*(angles/deg2rad(25)).^2);
        initial=radial*direction.*exp(1i*(.2*sin(w/wp)+.35*sin(angles)));
        initial=initial*(.02/kp)/sum(abs(initial),'all');
        multiplier=(.9+.2*cos(pi*w/wp)).*exp(.4i*(w/wp-1).^2);
        truth=initial.*multiplier;observed=sum(truth,2);
        [joint,condition]=allocate_directional_record(initial,observed);
        assert(norm(joint-truth,'fro')/norm(truth,'fro')<1e-12);
        [eta22,psi22]=gl_directional_sum_time(joint(:),W(:),kx(:),ky(:),g,h,kp,t,8,ib(:));
        eta20=gl_directional_difference_time(joint(:),W(:),kx(:),ky(:),h,t,16,ib(:));
        A=truth(:).';x=kx(:).';y=ky(:).';
        c=mf12_spectral_coefficients(2,g,h,real(A),imag(A),x,y,0,0,0,struct('enable_subharmonic',true));
        plus=1:2:numel(c.G_npm);minus=2:2:numel(c.G_npm);
        self=complex(c.A_2,c.B_2);cross=complex(c.A_npm,c.B_npm);
        E=exp(-1i*t*(2*c.omega));Ep=exp(-1i*t*c.omega_npm(plus));Em=exp(-1i*t*c.omega_npm(minus));
        ref22=real(E*(self.*c.G_2).'+Ep*(cross(plus).*c.G_npm(plus)).');
        refpsi=real(E*(1i*self.*c.mu_2).'+Ep*(1i*cross(plus).*c.mu_npm(plus)).');
        ref20=real(Em*(cross(minus).*c.G_npm(minus)).');
        energy=sum(abs(truth).^2,2);meanW=sum(w.*energy)/sum(energy);actualBeta=sqrt(sum((w-meanW).^2.*energy)/sum(energy))/meanW;
        err=[norm(real(eta22)-ref22)/norm(ref22),norm(real(psi22)-refpsi)/norm(refpsi),norm(eta20-ref20)/norm(ref20)];
        rows(end+1,:)={kph,beta,actualBeta,err(1),err(2),err(3),condition.max,numel(A)}; %#ok<AGROW>
        cases{depthIndex,ibeta}=struct('t',t,'eta22',real(eta22),'psi22',real(psi22),'eta20',eta20, ...
            'reference',[ref22,refpsi,ref20],'truth',truth,'initial',initial,'joint',joint,'w',w,'angles',angles,'h',h);
    end
end
summary=cell2table(rows,'VariableNames',{'kph','nominal_beta','measured_frequency_std_over_mean','eta22_L2','psi22_L2','eta20_L2','allocation_condition_max','parents'});
writetable(summary,fullfile(out,'summary.csv'));save(fullfile(out,'fields.mat'),'summary','cases','-v7.3');disp(summary);
set(groot,'defaultFigureVisible','off');f=figure('Color','w','Position',[50,50,1250,550]);tiledlayout(1,3,'Padding','compact');
for j=1:3
    nexttile;hold on;
    for dep=[1,5]
        use=summary.kph==dep;vars={'eta22_L2','psi22_L2','eta20_L2'};
        plot(summary.measured_frequency_std_over_mean(use),100*summary.(vars{j})(use),'-o','LineWidth',1.4);
    end
    grid on;xlabel('Actual relative frequency bandwidth');ylabel('Relative L2 vs spectral MF12 (%)');
    names={'eta22: GL8','psi22: GL2+2','eta20: GL16 nonzero spatial K'};title(names{j});legend('kh=1','kh=5','Location','best');
end
sgtitle('Controlled synthetic bandwidth benchmark; not OW3D validation');
exportgraphics(f,fullfile(out,'bandwidth_comparison.png'),'Resolution',180);exportgraphics(f,fullfile(out,'bandwidth_comparison.pdf'),'ContentType','vector');close(f);
end
