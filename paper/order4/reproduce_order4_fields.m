function metrics=reproduce_order4_fields(recompute,ranks)
%REPRODUCE_ORDER4_FIELDS Replot archived fields or rerun the no-Stokes GL graph.
if nargin<1,recompute=false;end
if nargin<2,ranks=[6,8,10];end
root=setup_green_laplace();here=fileparts(mfilename('fullpath'));addpath(here);
out=fullfile(root,'results','order4');if ~isfolder(out),mkdir(out);end
saved=load(fullfile(here,'data','eta44_no_stokes_rank_fields.mat'));
manifest=jsondecode(fileread(fullfile(here,'inputs','field_n0382','manifest.json')));
wit=read_directional_eta4_cpp_output(fullfile(here,'data','wit_field_n0382.bin'));
reference=real(ifft2(wit.etaFixedB(:,:,1).'*manifest.nx*manifest.ny));
assert(norm(reference(:)-saved.reference(:))/norm(reference(:))<1e-13);
assert(wit.nearResonantCount==0 && wit.allFinite);
if recompute
    input=prepare_order4_inputs(fullfile(here,'inputs','field_n0382'), ...
        fullfile(out,'field_input'));
    n=manifest.nx;axis=[0:n/2-1,-n/2:-1];[mx,my]=meshgrid(axis,axis);
end
fields=cell(1,numel(ranks));l2=zeros(numel(ranks),1);linf=l2;Q=l2;parity=l2;
for k=1:numel(ranks)
    old=saved.fields{find([6,8,10]==ranks(k),1)};
    if recompute
        [plus,audit]=finite_depth_directional_green_laplace_eta44_three_kernels( ...
            input.eta11_spectrum,manifest.h*manifest.dk*mx,manifest.h*manifest.dk*my, ...
            manifest.peakDepth,root,quadrature_rank=ranks(k));
        assert(~audit.all_nested_stokes_corrections_used);
        fields{k}=manifest.h*real(plus);
    else
        fields{k}=old;
    end
    difference=fields{k}-reference;
    l2(k)=norm(difference(:))/norm(reference(:));
    linf(k)=max(abs(difference),[],'all')/max(abs(reference),[],'all');
    Q(k)=norm(difference(:))/(norm(reference(:))+norm(fields{k}(:)));
    parity(k)=norm(fields{k}(:)-old(:))/norm(old(:));
    if recompute,assert(parity(k)<2e-10,'GL field changed from paper source.');end
end
metrics=table(ranks(:),l2,linf,Q,parity,'VariableNames', ...
    {'rank','raw_relative_l2','raw_relative_linf','Q','archived_field_relative_l2'});
writetable(metrics,fullfile(out,'field_metrics.csv'));
n=manifest.nx;x=((-n/2):(n/2-1))*manifest.domainPeakWavelengths/n;
ref=fftshift(reference);colors=lines(numel(ranks));
fig=figure('Color','w','Position',[100 100 1000 430]);
tiledlayout(1,2,'TileSpacing','compact');
nexttile;plot(x,ref(n/2+1,:),'k-','LineWidth',1.4,'DisplayName','WIT');hold on;
for k=1:numel(ranks)
    f=fftshift(fields{k});plot(x,f(n/2+1,:),'Color',colors(k,:), ...
        'LineWidth',1.2,'DisplayName',sprintf('GL%d, Q=%.4g',ranks(k),Q(k)));
end
xlim([-0.6 0.6]);xlabel('x / wavelength');ylabel('Elevation');legend('Location','best');grid on;
nexttile;hold on;
for k=1:numel(ranks)
    d=fftshift(fields{k}-reference);plot(x,d(n/2+1,:),'Color',colors(k,:), ...
        'LineWidth',1.2,'DisplayName',sprintf('GL%d - WIT',ranks(k)));
end
xlim([-0.6 0.6]);xlabel('x / wavelength');ylabel('Elevation difference');legend('Location','best');grid on;
exportgraphics(fig,fullfile(out,'eta44_gl_wit_fields.png'),'Resolution',200);close(fig);
disp(metrics);
end
