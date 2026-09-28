function screen_hos_high_xy_probes(out)
%SCREEN_HOS_HIGH_XY_PROBES Select off-axis probes before nonlinear scoring.
source='/home/lxy/green-laplace-unidirectional-time-series-runs/hos-directional-fourphase-20260926-v2';
s=jsondecode(fileread(fullfile(source,'settings.json')));
initial=load(fullfile(source,'inputs','initial_fields.mat'),'C','kx','ky','om');
use=abs(initial.C)>0 & initial.kx>0;
C=initial.C(use);kx=initial.kx(use);ky=initial.ky(use);om=initial.om(use);
C=C(:);kx=kx(:);ky=ky(:);om=om(:);
t=(0:s.output_dt_s:s.duration_s)';lambda=s.lambda_p;
dx_levels=[-1,-.5,.5,1];dy_levels=[-.75,-.5,-.25,.25,.5,.75];
[DX,DY]=ndgrid(dx_levels,dy_levels);
candidate_offsets=[DX(:),DY(:)];
centre_offsets=[dx_levels(:),zeros(numel(dx_levels),1)];
offsets=[centre_offsets;candidate_offsets];
xy=s.nominal_focus_xy(:).'+lambda*offsets;
phase=exp(1i*(kx*xy(:,1).'+ky*xy(:,2).'));
analytic=complex(zeros(numel(t),size(xy,1)));
for start=1:256:numel(C)
    ids=start:min(start+255,numel(C));
    analytic=analytic+exp(-1i*t*om(ids).')*(C(ids).*phase(ids,:));
end
physical=real(analytic);
max_abs=max(abs(physical),[],1);range=max(physical,[],1)-min(physical,[],1);
envelope=max(abs(analytic),[],1);
[~,peak_index]=max(abs(physical),[],1);peak_time=t(peak_index);
rows=zeros(size(candidate_offsets,1),14);
for index=1:size(candidate_offsets,1)
    dx=candidate_offsets(index,1);dy=candidate_offsets(index,2);
    centre=find(dx_levels==dx,1);
    point=numel(dx_levels)+index;
    rows(index,:)=[dx,dy,xy(point,:),max_abs(point),max_abs(centre), ...
        max_abs(point)/max_abs(centre),range(point),range(centre), ...
        range(point)/range(centre),envelope(point),envelope(point)/envelope(centre), ...
        peak_time(point),max_abs(point)/max_abs(1)];
end
metrics=array2table(rows,'VariableNames',{ ...
    'dx_lambda','dy_lambda','x_m','y_m','max_abs_m','same_x_centre_max_abs_m', ...
    'same_x_max_abs_ratio','peak_to_trough_m','same_x_centre_peak_to_trough_m', ...
    'same_x_peak_to_trough_ratio','max_envelope_m','same_x_envelope_ratio', ...
    'peak_time_s','nominal_focus_centre_max_abs_ratio'});
metrics.passes_one_third_gate=metrics.same_x_max_abs_ratio>=1/3;
writetable(metrics,fullfile(out,'screen.csv'));
report=struct('source',source,'Akp',s.Akp,'lambda_p_m',lambda, ...
    'candidate_dx_lambda',dx_levels,'candidate_dy_lambda',dy_levels, ...
    'gate','linear first-harmonic full-record max(abs eta) >= one third of same-x centreline max', ...
    'eligible_count',nnz(metrics.passes_one_third_gate), ...
    'candidate_count',height(metrics),'metrics',table2struct(metrics), ...
    'no_GL_executed',true,'no_HOS_reference_error_used_for_selection',true);
writejson(fullfile(out,'report.json'),report);
save(fullfile(out,'screen.mat'),'report','metrics','t','analytic','xy','offsets');
figure_handle=figure('Visible','off','Color','w','Position',[100,100,1000,650]);
scatter(metrics.dx_lambda,metrics.dy_lambda,90,metrics.same_x_max_abs_ratio,'filled');
hold on;failed=~metrics.passes_one_third_gate;
plot(metrics.dx_lambda(failed),metrics.dy_lambda(failed),'kx','MarkerSize',12,'LineWidth',1.5);
colorbar;grid on;xlabel('Delta x / lambda_p');ylabel('Delta y / lambda_p');
title('Linear wave-group amplitude / same-x centreline amplitude');
exportgraphics(figure_handle,fullfile(out,'screen.png'),'Resolution',160);
close(figure_handle);
disp(metrics)
end

function writejson(file,value)
fid=fopen(file,'w');assert(fid>=0);fprintf(fid,'%s\n',jsonencode(value));fclose(fid);
end
