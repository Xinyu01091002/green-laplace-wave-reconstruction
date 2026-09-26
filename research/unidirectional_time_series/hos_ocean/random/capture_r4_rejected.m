function capture_r4_rejected(root)
addpath(fullfile(root,'spark'));addpath(fullfile(root,'spark_diagnostic'),'-begin');s=load(fullfile(root,'input.mat'));
reports=struct;
for kind={'eta','psi'}
 kind=kind{1};target=fullfile(root,['r4_' kind '_rejected_capture.mat']);setenv('R4_CAPTURE_PATH',target);tic;
 try
  if strcmp(kind,'eta'),spark_eta20_r_series(s.F,s.kx,s.ky,s.h,4);else,spark_psi20_r_series(s.F,s.kx,s.ky,s.h,4);end
  error('Capture unexpectedly bypassed original rejection');
 catch ME
  assert(contains(ME.identifier,'RSeriesStability'));assert(isfile(target));
  reports.(kind)=struct('status','captured_before_original_rejection','identifier',ME.identifier,'message',ME.message,'operator_seconds',toc);
 end
end
setenv('R4_CAPTURE_PATH','');
ref=load(fullfile(root,'mf12.mat'));rows={};
f=figure('Visible','off','Position',[30 30 1300 850]);tiledlayout(2,2);
for kind={'eta','psi'}
 kind=kind{1};a=load(fullfile(root,['r4_' kind '_rejected_capture.mat']));
 if strcmp(kind,'eta'),field=2*a.eta20;else,field=2*sqrt(s.g)*a.psi20;end
 b=ref.([kind '20']);delta=field-b;assert(all(isfinite(field),'all'));
 rows(end+1,:)={kind,norm(delta(:))/norm(b(:)),max(abs(delta),[],'all')/max(abs(b),[],'all'),sqrt(mean(delta.^2,'all')),norm(field(:))/norm(b(:)),a.defect,a.outsideDefect,a.repairAudit.repaired_output_count,a.repairAudit.evaluated_oriented_pairs,a.repairAudit.seconds};
 nexttile;plot((0:s.nx-1)*s.Lx/s.nx,b(s.ny/2+1,:),'k', (0:s.nx-1)*s.Lx/s.nx,field(s.ny/2+1,:),'r--');grid on;xlabel('x (m)');ylabel([kind '20']);title('Centerline: R4 field before rejection');legend('MF12','R4 diagnostic','Location','best');
 nexttile;imagesc(delta);axis image;colorbar;title([kind '20: raw R4 diagnostic minus MF12']);
end
T=cell2table(rows,'VariableNames',{'field','relative_L2','relative_Linf','absolute_RMS','norm_ratio','preprojection_defect','outside_repair_defect','repair_output_count','repair_oriented_pairs','repair_seconds'});writetable(T,fullfile(root,'rejected_diagnostic_comparison.csv'));
exportgraphics(f,fullfile(root,'rejected_diagnostic_comparison.png'),'Resolution',160);close(f);
report=struct('status','DIAGNOSTIC_ONLY_ORIGINAL_R4_REJECTED','captures',reports,'metrics',table2struct(T),'note','Only save-before-check instrumentation. Original formula, low-output repair, threshold and rejection unchanged. These scores characterize the rejected intermediate fields; not production acceptance or HOS readiness.');
fid=fopen(fullfile(root,'rejected_diagnostic_comparison.json'),'w');fprintf(fid,'%s',jsonencode(report,PrettyPrint=true));fclose(fid);disp(T);
end
