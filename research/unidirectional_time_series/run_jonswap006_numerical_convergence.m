function run_jonswap006_numerical_convergence(out)
%RUN_JONSWAP006_NUMERICAL_CONVERGENCE Fixed-input grid/rank audit at probe 1.
run='/home/lxy/green-laplace-unidirectional-time-series-runs/hos-jonswap-kphs006-20tp-20260927-v1/medium';
s=jsondecode(fileread(fullfile(run,'settings.json')));initial=load(fullfile(run,'inputs','initial_fields.mat'),'C','kx','ky','om');
phases=[0,90,180,270];v=read_probe_file(fullfile(run,'cases','phi000','Results','probes.dat'));
t=v(:,1);raw=zeros(numel(t),5,4);raw(:,:,1)=v(:,2:end);
for p=2:4
    v=read_probe_file(fullfile(run,'cases',sprintf('phi%03d',phases(p)),'Results','probes.dat'));
    assert(isequal(t,v(:,1)));
    raw(:,:,p)=v(:,2:end);
end
H=zeros(size(raw));for p=1:4,H(:,:,p)=imag(hilbert(raw(:,:,p)));end
first=(raw(:,:,1)-raw(:,:,3)-H(:,:,2)+H(:,:,4))/4;
ref22=(raw(:,:,1)-raw(:,:,2)+raw(:,:,3)-raw(:,:,4))/4;
ref33=(raw(:,:,1)-raw(:,:,3)+H(:,:,2)-H(:,:,4))/4;
score=t>=s.scoring_window_s(1)&t<=s.scoring_window_s(2);
[model,info]=prepare_directional_joint_wavegroup(initial,first(:,1),t,s,1,7.5,.999,20);

[qx2,qy2]=limits(model,s.h,2,256);[qx3,qy3]=limits(model,s.h,3,256);
grid22=[256,512];grid33=[256,512];
[values22,rows22]=grid_runs(@eta22,grid22,qx2,qy2,ref22(:,1),score);
[values33,rows33]=grid_runs(@eta33,grid33,qx3,qy3,ref33(:,1),score);
grid_metrics=[cell2table(rows22,'VariableNames',names());cell2table(rows33,'VariableNames',names())];
writetable(grid_metrics,fullfile(out,'grid_metrics.csv'));

ranks22=[6,8,10,12];ranks33=[6,8,10,12];
rank_rows=[rank_runs(@eta22,ranks22,256,qx2,qy2,ref22(:,1),score); ...
    rank_runs(@eta33,ranks33,256,qx3,qy3,ref33(:,1),score)];
rank_metrics=cell2table(rank_rows,'VariableNames',{'component','rank','grid_size','seconds','change_to_highest_rank','HOS_relative','norm_ratio'});
writetable(rank_metrics,fullfile(out,'rank_metrics.csv'));
report=struct('source_run',run,'probe',1,'input_audit',rmfield(info,'input'), ...
    'grid_metrics',table2struct(grid_metrics),'rank_metrics',table2struct(rank_metrics), ...
    'selection_rule','fixed observed input and fixed physical modal limits; numerical changes do not use HOS error for selection', ...
    'pair_loops',0,'triple_loops',0);
writejson(fullfile(out,'report.json'),report);save(fullfile(out,'convergence.mat'),'report','grid_metrics','rank_metrics','values22','values33','t','score');
disp(grid_metrics);disp(rank_metrics)

    function value=eta22(grid,rank,qxlim,qylim)
        value=real(gl_directional_time_modal_grid(model.A,model.omega,model.kx,model.ky,s.g,s.h,s.kp,t-t(1),struct('Nx',grid,'Ny',grid,'J',rank,'QxLimit',qxlim,'QyLimit',qylim,'Baseband',true)));
    end
    function value=eta33(grid,rank,qxlim,qylim)
        value=real(gl_directional_time_modal_grid_eta33(model.A,model.omega,model.kx,model.ky,s.g,s.h,s.kp,t-t(1),struct('Nx',grid,'Ny',grid,'J',rank,'QxLimit',qxlim,'QyLimit',qylim,'Baseband',true)));
    end
end

function [values,rows]=grid_runs(executor,grids,qx,qy,reference,mask)
values=cell(1,numel(grids));rows=cell(numel(grids),8);
for i=1:numel(grids),clock=tic;values{i}=executor(grids(i),default_rank(executor),qx,qy);seconds=toc(clock);rows(i,:)={component_name(executor),grids(i),default_rank(executor),seconds,NaN,relative(values{i},reference,mask),norm(values{i}(mask))/norm(reference(mask)),NaN};end
for i=1:numel(grids),rows{i,5}=relative(values{i},values{end},mask);rows{i,8}=norm(values{i}(mask)-values{end}(mask))/norm(values{end}(mask));end
end
function rows=rank_runs(executor,ranks,grid,qx,qy,reference,mask)
values=cell(1,numel(ranks));rows=cell(numel(ranks),7);
for i=1:numel(ranks),clock=tic;values{i}=executor(grid,ranks(i),qx,qy);seconds=toc(clock);rows(i,:)={component_name(executor),ranks(i),grid,seconds,NaN,relative(values{i},reference,mask),norm(values{i}(mask))/norm(reference(mask))};end
for i=1:numel(ranks),rows{i,5}=relative(values{i},values{end},mask);end
end
function rank=default_rank(executor),name=func2str(executor);if contains(name,'eta22'),rank=8;else,rank=10;end,end
function name=component_name(executor),text=func2str(executor);if contains(text,'eta22'),name='eta22';else,name='eta33';end,end
function n=names(),n={'component','grid_size','rank','seconds','change_to_finest_grid','HOS_relative','norm_ratio','duplicate_change_check'};end
function [qx,qy]=limits(model,h,order,N),qx=(order+.05)*max(abs(h*model.kx))/(1-2/N);qy=(order+.05)*max(abs(h*model.ky))/(1-2/N);end
function value=relative(a,b,mask),value=norm(a(mask)-b(mask))/norm(b(mask));end
function v=read_probe_file(file),lines=readlines(file);header=find(startsWith(strtrim(lines),'VARIABLES'),1,'last');v=readmatrix(file,'NumHeaderLines',header);v=v(all(isfinite(v),2),:);end
function writejson(file,value),fid=fopen(file,'w');assert(fid>=0);fprintf(fid,'%s\n',jsonencode(value));fclose(fid);end
