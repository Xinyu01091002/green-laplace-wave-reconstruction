function count_energy_support(root)
src=fullfile(root,'jonswap-design-v1','linear_design.mat');s=load(src);out=fullfile(root,'jonswap-energy-support-v1');assert(~isfolder(out));mkdir(out);
C=s.C(:);energy=abs(C).^2/2;total=sum(energy);[sorted,order]=sort(energy,'descend');cum=cumsum(sorted)/total;
kx=s.kx(s.use);ky=s.ky(s.use);k=hypot(kx,ky);h=s.report.h;g=9.81;omega=sqrt(g*k.*tanh(k*h));freq=omega/(2*pi);theta=atan2(ky,kx);fp=s.report.fp;
rows={};details=cell(1,3);levels=[.99 .995 .999];
for a=1:numel(levels)
 n=find(cum>=levels(a),1);keep=order(1:n);retained=cum(n);v=energy(keep)/sum(energy(keep));mu=sum(v.*freq(keep));beta=sqrt(sum(v.*(freq(keep)-mu).^2))/mu;
 rows(end+1,:)={levels(a),n,n/numel(C),retained,1/sqrt(retained),n*(n-1),8*n*(n-1)/2^30,beta,min(freq(keep))/fp,max(freq(keep))/fp};
 bands=[.5 1;1 1.5;1.5 2;2 2.5];bandRows=zeros(4,4);mask=false(size(C));mask(keep)=true;
 for j=1:4,b=freq>=bands(j,1)*fp&freq<=bands(j,2)*fp;bandRows(j,:)=[bands(j,:),sum(energy(b))/total,sum(energy(b&mask))/max(sum(energy(b)),realmin)];end
 details{a}=struct('target_energy',levels(a),'frequency_bands_columns',{{'fmin_over_fp','fmax_over_fp','original_total_energy_fraction','fraction_of_band_retained'}},'frequency_bands',bandRows,'direction_rms_deg',sqrt(sum(v.*theta(keep).^2))*180/pi,'max_abs_direction_deg',max(abs(theta(keep)))*180/pi,'max_self_second_Kx_over_kp',2*max(abs(kx(keep)))/s.report.kp,'max_self_second_Ky_over_kp',2*max(abs(ky(keep)))/s.report.kp);
end
T=cell2table(rows,'VariableNames',{'target_energy','retained_modes','fraction_of_modes','actual_retained_energy','amplitude_renormalization','ordered_cross_pairs','one_pair_real_array_GiB','relative_bandwidth','min_f_over_fp','max_f_over_fp'});writetable(T,fullfile(out,'counts.csv'));
r=struct('status','ENERGY_COUNT_ONLY_NO_MF12_NO_HOS','source_design',src,'steepness_definition','kp Hs / 2', 'steepness_values',[.02 .12],'Hs_values',2*[.02 .12]/s.report.kp,'original_modes',numel(C),'original_ordered_pairs',numel(C)*(numel(C)-1),'energy_definition','abs(C_j)^2/2 on the actual Cartesian discretization','details',{details},'note','Counts independent of uniform amplitude rescaling. Threshold is relative to already truncated .5--2.5 fp directional spectrum. Ranked modal-energy pruning is not a bound on second-harmonic error.');
f=fopen(fullfile(out,'report.json'),'w');fprintf(f,'%s',jsonencode(r,PrettyPrint=true));fclose(f);save(fullfile(out,'retained_indices.mat'),'order','cum','levels','T','-v7.3');disp(T);disp(r);
end
