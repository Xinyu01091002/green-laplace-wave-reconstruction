function replot_spectra(root)
out=fullfile(root,'fourphase-spectrum-audit');fig=figure('Visible','off','Position',[50 50 1650 850]);tl=tiledlayout(2,3,'TileSpacing','compact','Padding','compact');colors=[.05 .35 .72;.9 .36 .05;.49 .15 .65];
for akp=[.02 .12],for tp=[0 3 20]
    tab=readtable(fullfile(out,sprintf('spectrum_akp%03d_t%02dTp.csv',round(100*akp),tp)));a=nexttile;hold(a,'on');set(a,'YScale','log');
    for h=1:3,plot(a,tab{:,1},max(tab{:,h+1},realmin),'Color',colors(h,:),'LineWidth',1.4);end
    xlim([0 6]);ylim([1e-12 .2]);yticks(10.^(-12:2:-2));grid on;for n=1:3,xline(n,':','Color',[.65 .65 .65],'HandleVisibility','off');end
    xlabel('k / k_p');ylabel('One-sided modal amplitude (m)');title(sprintf('Akp = %.2f, t = %g T_p',akp,tp));if akp==.02&&tp==0,legend('First harmonic','Second harmonic','Third harmonic','Location','southwest');end
end,end
title(tl,'HOS spatial spectra after Hilbert four-phase separation (no band-pass filter)');
exportgraphics(fig,fullfile(out,'fourphase_spatial_spectra_final.png'),'Resolution',170);exportgraphics(fig,fullfile(out,'fourphase_spatial_spectra_final.pdf'),'ContentType','vector');close(fig);
end