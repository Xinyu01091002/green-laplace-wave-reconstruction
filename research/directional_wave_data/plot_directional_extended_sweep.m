function plot_directional_extended_sweep(out)
s=load(fullfile(out,'summary.mat'));set(groot,'defaultFigureVisible','off');
for kind=1:2
    folders=s.folders;
    chosen=contains(string(folders),'kh1_s25_a002');
    if kind==2,chosen=endsWith(string(folders),[filesep,'center']);end
    folders=folders(chosen);if isempty(folders),continue;end
    f=figure('Color','w','Position',[30,30,1550,220*numel(folders)]);tl=tiledlayout(numel(folders),3,'Padding','compact','TileSpacing','compact');
    for j=1:numel(folders)
        d=load(fullfile(folders{j},'joint_pilot.mat'));a=load(fullfile(folders{j},'eta20_fields.mat'));
        for col=1:3
            ref=d.d.eta2;pred=d.eta(:,2);label='eta22 (m)';
            if col==2,ref=d.d.psi2;pred=d.psi(:,2);label='psi22 (m^2/s)';end
            if col==3,ref=a.reference20{2};pred=a.varying{2}(:,3);label='eta20, 3 wp/Nyquist (m)';end
            nexttile;plot(d.t,ref,'ko','MarkerSize',3);hold on;plot(d.t,pred,'-','Color',[0,.4,.7],'LineWidth',1.1);
            xlim(d.report.display_limits_s);grid on;ylabel(label);
            if col==1
                if kind==1,title(sprintf('x=%.0f, y=%.0f m',d.report.metadata.probe));
                else,title(sprintf('kh=%g, spread=%g, Akp=%.2f',round(d.report.metadata.kph),d.report.metadata.spread_label_degrees,d.report.metadata.Akp));end
            end
        end
    end
    title(tl,'Directional joint input: OW3D saved samples (circles), GL (blue)');
    name='far_probes';if kind==2,name='case_matrix';end
    exportgraphics(f,fullfile(out,[name,'.png']),'Resolution',160);exportgraphics(f,fullfile(out,[name,'.pdf']),'ContentType','vector');close(f);
end
end
