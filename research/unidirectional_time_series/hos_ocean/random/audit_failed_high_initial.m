function report=audit_failed_high_initial(run,out)
% Read-only initial-field and saved-log diagnostics. No solver restart/change.
assert(~isfolder(out));mkdir(out);levels={'low','high'};rows=cell(0,12);summaries=cell(2,1);
for il=1:2
    folder=fullfile(run,levels{il});s=jsondecode(fileread(fullfile(folder,'settings.json')));
    d=load(fullfile(folder,'inputs','initial_fields.mat'),'E','P','C','kx','ky','om');
    nx=s.grid(1);ny=s.grid(2);Lx=s.domain_m(1);Ly=s.domain_m(2);
    [KX,KY]=meshgrid(2*pi/Lx*[0:nx/2-1,-nx/2:-1],2*pi/Ly*[0:ny/2-1,-ny/2:-1]);
    K=hypot(KX,KY);first=(d.E(:,:,1)-1i*d.E(:,:,2)-d.E(:,:,3)+1i*d.E(:,:,4))/2;
    for ip=1:4
        e=d.E(:,:,ip);p=d.P(:,:,ip);F=fft2(e);
        ex=real(ifft2(1i*KX.*F));ey=real(ifft2(1i*KY.*F));
        exFine=real(interpft(interpft(ex,2*ny,1),2*nx,2));
        [peakSlope,idx]=max(abs(exFine),[],'all','linear');[iy,ix]=ind2sub(size(exFine),idx);
        e1=real(first*exp(1i*(ip-1)*pi/2));
        power=abs(F).^2;slopePower=KX.^2.*power;high=K>3*s.kp;
        fileError=zeros(1,2);
        for rank=0:s.MPI_ranks_per_phase-1
            file=fullfile(folder,'cases',sprintf('phi%03d',(ip-1)*90),'Results',sprintf('3d_ini_%03d.dat',rank));
            fid=fopen(file,'r');assert(fid>=0);guard=onCleanup(@()fclose(fid));
            line=fgetl(fid);while ischar(line)&&~startsWith(strtrim(line),'ZONE'),line=fgetl(fid);end
            assert(ischar(line));a=fscanf(fid,'%f',[2,Inf]);clear guard;
            cols=rank*(ny/s.MPI_ranks_per_phase)+(1:ny/s.MPI_ranks_per_phase);
            ee=e(cols,:).';pp=p(cols,:).';assert(size(a,2)==numel(ee));
            fileError=max(fileError,[max(abs(a(1,:).'-ee(:))),max(abs(a(2,:).'-pp(:)))]);
        end
        rows(end+1,:)={levels{il},(ip-1)*90,max(abs(e),[],'all'), ...
            max(abs(ex),[],'all'),peakSlope,max(abs(ey),[],'all'), ...
            norm(e-e1,'fro')/norm(e1,'fro'),sum(power(high))/sum(power,'all'), ...
            sum(slopePower(high))/sum(slopePower,'all'), ...
            [(ix-1)*Lx/(2*nx),(iy-1)*Ly/(2*ny)],fileError(1),fileError(2)}; %#ok<AGROW>
    end
    kapp=hypot(d.kx,d.ky);summaries{il}=struct('level',levels{il}, ...
        'kpHs_over_2',s.kpHs_over_2,'Hs_m',s.Hs_linear, ...
        'max_parent_k_over_kp',max(kapp)/s.kp,'max_parent_frequency_over_fp',max(d.om)/(2*pi/s.Tp), ...
        'x_Nyquist_over_kp',pi*nx/Lx/s.kp,'y_Nyquist_over_kp',pi*ny/Ly/s.kp);
end
metrics=cell2table(rows,'VariableNames',{'family','phase_deg','max_initial_abs_eta_m', ...
    'max_initial_abs_eta_x','interpolated2x_max_abs_eta_x','max_initial_abs_eta_y', ...
    'second_order_eta_over_first_L2','eta_power_fraction_k_above_3kp', ...
    'x_slope_power_fraction_k_above_3kp','initial_max_x_slope_xy_m','export_eta_max_abs_error','export_psi_max_abs_error'});
report=struct('scope','Initial base/refined-grid diagnostics and ASCII export equality; not a failed-time field reconstruction', ...
    'source_run',run,'settings',[summaries{:}],'metrics',table2struct(metrics), ...
    'failure_guard_threshold',10,'logged_failure_slope_x',10.08, ...
    'last_saved_high_phi270_time_s',5.6,'failed_field_saved',false, ...
    'new_HOS_runs',false,'physical_input_or_threshold_changed',false);
writetable(metrics,fullfile(out,'initial_diagnostics.csv'));
fid=fopen(fullfile(out,'report.json'),'w');assert(fid>=0);guard=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(report));disp(metrics(:,[1:9,11:12]));disp([summaries{:}]);
end
