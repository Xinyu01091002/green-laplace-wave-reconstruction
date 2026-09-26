function prepare_random_run(run,label,mf12)
if label==.02 && isfile(fullfile(fileparts(run),'cancel-low.json'))
 error('HOS:UserCancelledLow','USER_CANCELLED_LOW_BEFORE_INITIALIZATION');
end
base='/home/lxy/green-laplace-unidirectional-time-series/results/ow3d_redesign/20260925-compact';
mkdir(fullfile(run,'inputs'));
if label==.12
 source=fullfile(base,'random-phase-only-v1','initial_fields.mat');
 copyfile(source,fullfile(run,'inputs','initial_fields.mat'));
else
 assert(label==.02);
 focused='/home/lxy/green-laplace-unidirectional-time-series-runs/hos-directional-akp002-20260926-v1/inputs/initial_fields.mat';
 geometry='/home/lxy/green-laplace-unidirectional-time-series-runs/hos-directional-akp002-20260926-v1/inputs/design.mat';
 prepare_compact_wavegroup(geometry,mf12,fullfile(run,'generated-initial'),focused,20260925);
 copyfile(fullfile(run,'generated-initial','initial_fields.mat'),fullfile(run,'inputs','initial_fields.mat'));
end
s=load(fullfile(run,'inputs','initial_fields.mat'));
assert(s.report.randomization.seed==20260925);
assert(abs(s.report.design.kp*sum(abs(s.C))-label)<1e-10);
if label==.02
 high=load(fullfile(base,'random-phase-only-v1','initial_fields.mat'),'C','kx','ky');
 assert(isequal(s.kx,high.kx)&&isequal(s.ky,high.ky));
 assert(norm(s.C-high.C/6)/norm(s.C)<1e-12,'High/low random phases differ');
end
mkdir(fullfile(run,'staging4'));
prepare_directional_hos(fullfile(run,'staging4'),fullfile(run,'inputs','initial_fields.mat'));
p=fullfile(run,'staging4','settings.json');d=jsondecode(fileread(p));
d.family='phase_only_random';d.seed=20260925;d.nominal_focus_s=[];d.nominal_focus_Tp=[];d.nominal_focus_xy=[];
d.linear_Hs_4sigma_m=4*sqrt(sum(abs(s.C).^2)/2);d.kp_Hs_over_2=d.kp*d.linear_Hs_4sigma_m/2;
d.amplitude_rule='Preserve focused modal amplitudes; Akp is potential focusing label, not sea-state steepness';
d.scoring_window_s=[3*d.Tp,d.duration_s-3*d.Tp];
f=fopen(p,'w');assert(f>=0);fprintf(f,'%s',jsonencode(d,PrettyPrint=true));fclose(f);
disp(d);
end
