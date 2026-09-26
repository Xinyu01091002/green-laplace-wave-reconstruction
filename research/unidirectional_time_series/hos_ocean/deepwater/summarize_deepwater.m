function summarize_deepwater(root)
rows=[];
for akp=[.02,.12]
    name=sprintf('akp%03d',round(100*akp));if akp==.02,final='tol16';previous='tol14';else,final='tol14';previous='tol12';end;s=load(fullfile(root,final,name,'comparison.mat'));p=load(fullfile(root,previous,name,'comparison.mat'));ini=load(fullfile(root,name,'initial.mat'),'bins');N=numel(s.x);signed=[0:N/2,-N/2+1:-1]';below=abs(signed)<3*min(ini.bins);
    tvals=[0,3,20];for it=1:3
        f=fft(s.obs3(:,it));err=fft(s.obs3(:,it)-s.pred3(:,it));
        rows(end+1,:)=[akp,tvals(it),norm(f(below))/norm(f),norm(err(below))/norm(err),norm(err(~below))/norm(f(~below)),norm(s.obs2(:,it)-p.obs2(:,it))/norm(s.obs2(:,it)),norm(s.obs3(:,it)-p.obs3(:,it))/norm(s.obs3(:,it))];
    end
    fig=figure('Visible','off','Position',[50,50,1300,750]);tiledlayout(2,2);lambda=2*pi/s.report.settings.kp;
    for h=2:3,for it=2:3
        z=fft(s.first(:,it));z(N/2+1:end)=0;z(2:N/2)=2*z(2:N/2);[~,center]=max(abs(ifft(z)));
        dx=mod(s.x-s.x(center)+s.report.settings.L/2,s.report.settings.L)-s.report.settings.L/2;[xx,order]=sort(dx/lambda);mask=abs(xx)<=4;
        if h==2,o=s.obs2(:,it);pred=s.pred2(:,it);else,o=s.obs3(:,it);pred=s.pred3(:,it);end
        nexttile;plot(xx(mask),o(order(mask)),'k-',xx(mask),pred(order(mask)),'r--','LineWidth',1);grid on;xlabel('(x - x_{envelope}) / lambda_p');ylabel(sprintf('harmonic %d, eta (m)',h));title(sprintf('Akp=%.2f, t=%g Tp',akp,s.report.metrics(2*it-1).time_Tp));if h==2&&it==2,legend('HOS four-phase','MF12 reconstruction','Location','northwest','FontSize',9);end
    end,end
    exportgraphics(fig,fullfile(root,final,name,'comparison_zoom.png'),'Resolution',150);close(fig);
end
tab=array2table(rows,'VariableNames',{'Akp','time_Tp','third_below_triple_support_relative_norm','third_error_below_support_fraction','third_error_within_support_relative_L2','previous_vs_final_observed_second_relative_L2','previous_vs_final_observed_third_relative_L2'});writetable(tab,fullfile(root,'support_and_tolerance_audit_final.csv'));disp(tab);
end