function prepare_low_run(root)
geometry='/home/lxy/green-laplace-unidirectional-time-series/results/ow3d_redesign/20260925-compact/geometry-v1/first_order_design.mat';
high='/home/lxy/green-laplace-unidirectional-time-series-runs/hos-directional-fourphase-20260926-v2/inputs/initial_fields.mat';
d=load(geometry);d.report.Akp_group=.02;d.report.group_linear_focus_amplitude_m=.02/d.report.kp;d.report.leads=[];d.report.random_prototype=[];d.report.status='AKP002_FROM_FROZEN_GEOMETRY_NEW_MF12_COEFFICIENTS';
mkdir(fullfile(root,'inputs'));save(fullfile(root,'inputs','design.mat'),'-struct','d','-v7.3');
prepare_compact_low(fullfile(root,'inputs','design.mat'),fullfile(root,'mf12'),fullfile(root,'generated-initial'));
low=load(fullfile(root,'generated-initial','initial_fields.mat'));old=load(high);a=.02/.12;assert(norm(low.C-a*old.C)/norm(low.C)<1e-12);
checkE=zeros(size(low.E));checkP=checkE;
for j=1:4
 le=real(old.eta1*exp(1i*(j-1)*pi/2));lp=real(old.psi1*exp(1i*(j-1)*pi/2));
 checkE(:,:,j)=a*le+a^2*(old.E(:,:,j)-le);checkP(:,:,j)=a*lp+a^2*(old.P(:,:,j)-lp);
end
audit=struct('Akp',.02,'coefficient_recomputed',true,'total_field_not_linearly_scaled',true,'eta_homogeneity_relative_L2',norm(low.E(:)-checkE(:))/norm(low.E(:)),'psi_homogeneity_relative_L2',norm(low.P(:)-checkP(:))/norm(low.P(:)));
assert(audit.eta_homogeneity_relative_L2<1e-11&&audit.psi_homogeneity_relative_L2<1e-11);
copyfile(fullfile(root,'generated-initial','initial_fields.mat'),fullfile(root,'inputs','initial_fields.mat'));
mkdir(fullfile(root,'staging4'));prepare_directional_hos(fullfile(root,'staging4'),fullfile(root,'inputs','initial_fields.mat'));
f=fopen(fullfile(root,'initialization-audit.json'),'w');fprintf(f,'%s',jsonencode(audit,PrettyPrint=true));fclose(f);disp(audit);
end