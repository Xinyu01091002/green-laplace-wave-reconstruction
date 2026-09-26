function plot_full_gl_record(root)
d=load(fullfile(root,'record-direction-audit-v2','audit.mat'),'T','cases');T=d.T;
f=figure('Visible','off','Position',[30 30 1400 1150]);tiledlayout(5,1,'TileSpacing','compact','Padding','compact');
offsets=[0 -52.7823 52.7823 -70.3762 70.3762];
for p=1:5
 ix=find(T.probe==p&T.N==1101&T.angle_deg==1.875&T.taper_s==0&~T.full_Hilbert_prefix,1);a=d.cases{ix};
 nexttile;plot(a.t,a.filtered(:,1),'k-',a.t,a.filtered(:,2),'r--','LineWidth',1);hold on;xline(41.285991,':');xline(178.714009,':');grid on;xlim([0 220]);xticks(0:20:220);ylabel('eta_2 (m)');title(sprintf('y offset %.2f m | full 0-220 s record',offsets(p)));
 if p==1,legend('HOS phase sector','GL 1.875 deg','Location','northwest','FontSize',9);end
end
xlabel('Time (s)');exportgraphics(f,fullfile(root,'spatial-view-v3','random_full_GL_comparison.png'),'Resolution',160);close(f);
end