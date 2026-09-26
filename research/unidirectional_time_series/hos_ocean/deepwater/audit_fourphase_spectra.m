function audit_fourphase_spectra(root)
out=fullfile(root,'fourphase-spectrum-audit');assert(~isfolder(out));mkdir(out);
coef=[1 0 -1 0 0 -1 0 1;1 -1 1 -1 0 0 0 0;1 0 -1 0 0 1 0 -1;1 1 1 1 0 0 0 0]/4;
fig=figure('Visible','off','Position',[50 50 1650 850]);tl=tiledlayout(2,3,'TileSpacing','compact','Padding','compact');colors=[.05 .35 .72;.9 .36 .05;.49 .15 .65];records=[];
for row=1:2
    if row==1,akp=.02;folder=fullfile(root,'tol16','akp002');else,akp=.12;folder=fullfile(root,'tol14','akp012');end
    s=load(fullfile(folder,'comparison.mat'));r=s.report.settings;N=r.N;E=zeros(N,4,3);P=E;targets=[0 3 20]*r.Tp;
    for j=1:4
        f=fopen(fullfile(folder,sprintf('phase%03d',(j-1)*90),'Results','3d.dat'));assert(f>=0);frames=0;hits=false(1,3);
        while ~feof(f)
            line=fgetl(f);if ~ischar(line),break;end;if ~startsWith(strtrim(line),'ZONE'),continue;end
            tok=regexp(line,'SOLUTIONTIME\s*=\s*([+\-\d.Ee]+)','tokens','once');t=str2double(tok{1});frames=frames+1;if frames==1,nc=4;else,nc=2;end
            v=fscanf(f,'%f',[nc,N]);assert(isequal(size(v),[nc,N]));it=find(abs(t-targets)<1e-8*r.Tp);
            if ~isempty(it),E(:,j,it)=v(end-1,:).';P(:,j,it)=v(end,:).';hits(it)=true;end
        end
        fclose(f);assert(all(hits));
    end
    k=(0:N/2)'*2*pi/r.L;krel=k/r.kp;
    for it=1:3
        H=-imag(hilbert(E(:,:,it)));legacy=[E(:,:,it),H]*coef.';
        HP=-imag(hilbert(P(:,:,it)));legacyP=[P(:,:,it),HP]*coef.';
        old=[s.first(:,it),s.obs2(:,it),s.obs3(:,it)];delta=sqrt(sum((legacy(:,1:3)-old).^2))./sqrt(sum(old.^2));assert(all(delta<1e-9));
        oldP=old_operator(P(:,:,it));deltaP=sqrt(sum((legacyP(:,1:3)-oldP).^2))./sqrt(sum(oldP.^2));assert(all(deltaP<1e-9));
        record=struct('Akp',akp,'time_Tp',targets(it)/r.Tp,'eta_operator_relative_L2',delta,'psi_operator_relative_L2',deltaP);records=[records;record];
        A=abs(fft(legacy(:,1:3)))/N;A=A(1:N/2+1,:);A(2:end-1,:)=2*A(2:end-1,:);
        tab=array2table([krel,A],'VariableNames',{'k_over_kp','first_amplitude_m','second_amplitude_m','third_amplitude_m'});writetable(tab,fullfile(out,sprintf('spectrum_akp%03d_t%02dTp.csv',round(100*akp),round(targets(it)/r.Tp))));
        nexttile;hold on;for h=1:3,semilogy(krel,max(A(:,h),realmin),'Color',colors(h,:),'LineWidth',1.3);end
        set(gca,'YScale','log');xlim([0 6]);ylim([1e-12 .2]);grid on;for n=1:3,xline(n,':','Color',[.65 .65 .65],'HandleVisibility','off');end
        xlabel('k / k_p');ylabel('One-sided modal amplitude (m)');title(sprintf('Akp = %.2f, t = %g T_p',akp,targets(it)/r.Tp));
        if row==1&&it==1,legend('First harmonic','Second harmonic','Third harmonic','Location','southwest');end
    end
end
title(tl,'HOS spatial spectra after Hilbert four-phase separation — before band-pass filtering');
exportgraphics(fig,fullfile(out,'fourphase_spatial_spectra.png'),'Resolution',170);exportgraphics(fig,fullfile(out,'fourphase_spatial_spectra.pdf'),'ContentType','vector');close(fig);
f=fopen(fullfile(out,'operator-audit.json'),'w');fprintf(f,'%s',jsonencode(records,PrettyPrint=true));fclose(f);disp(records);
end
function o=old_operator(E)
N=size(E,1);c1=(E(:,1)-1i*E(:,2)-E(:,3)+1i*E(:,4))/4;c3=conj(c1);
f=fft(c1);z=zeros(N,1);z(2:N/2)=2*f(2:N/2);o(:,1)=real(ifft(z));
o(:,2)=(E(:,1)-E(:,2)+E(:,3)-E(:,4))/4;
f=fft(c3);z(:)=0;z(2:N/2)=2*f(2:N/2);o(:,3)=real(ifft(z));
end