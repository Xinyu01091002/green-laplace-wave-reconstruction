function compare_directional_partial(run)
out=fullfile(run,'partial-comparison-20260926-v1');addpath(fullfile(out,'source','research','directional_wave_data'));addpath(fullfile(out,'source','research','unidirectional_time_series'));
settings=jsondecode(fileread(fullfile(run,'settings.json')));initial=load(fullfile(run,'inputs','initial_fields.mat'),'C','kx','ky','om');snapshot=jsondecode(fileread(fullfile(out,'snapshot.json')));
N=snapshot.sample_count;raw=zeros(N,5,4);
for j=1:4,a=readmatrix(fullfile(out,sprintf('phi%03d.csv',(j-1)*90)));assert(isequal(size(a),[N,6]));t=a(:,1);raw(:,:,j)=a(:,2:end);end
HH=imag(hilbert(raw));first=(raw(:,:,1)-raw(:,:,3)-HH(:,:,2)+HH(:,:,4))/4;reference=(raw(:,:,1)-raw(:,:,2)+raw(:,:,3)-raw(:,:,4))/4;
ratio=max(abs(reference),[],1)/max(abs(reference(:,1)));eligible=ratio>=1/3;assert(eligible(1));allMetrics=table();results=cell(1,5);records=cell(1,5);timers=zeros(1,5);
for probe=1:5
 if ~eligible(probe),continue;end
 timer=tic;[results{probe},records{probe},metrics]=one(first(:,probe),reference(:,probe),t,initial,settings,probe);timers(probe)=toc(timer);allMetrics=[allMetrics;metrics];
end
% Fixed 10 s shorter prefix diagnoses record-length/frequency-grid sensitivity.
ns=N-50;shortRaw=raw(1:ns,1,:);hh=imag(hilbert(shortRaw));shortFirst=(shortRaw(:,1,1)-shortRaw(:,1,3)-hh(:,1,2)+hh(:,1,4))/4;
[shortResult,shortReport,shortMetrics]=one(shortFirst,reference(1:ns,1),t(1:ns),initial,settings,1);
full=results{1};mask=t(1:ns)>=full.limits(1)&t(1:ns)<=full.limits(2);assert(full.limits(2)<t(ns));
windowAudit=struct('full_end_s',t(end),'short_end_s',t(ns),'main_window_s',full.limits,'first_input_relative_L2',norm(shortResult.input(mask)-full.input(find(mask)))/norm(full.input(find(mask))), ...
 'GL7p5_prediction_relative_L2',norm(shortResult.pred(mask,2)-full.pred(find(mask),2))/norm(full.pred(find(mask),2)), ...
 'common_band_reference_relative_L2',norm(shortResult.filtered(mask,1)-full.filtered(find(mask),1))/norm(full.filtered(find(mask),1)), ...
 'note','Sensitivity includes Hilbert endpoints, frequency grid and input-derived support; no target-selected window');
writetable(allMetrics,fullfile(out,'metrics.csv'));writetable(shortMetrics,fullfile(out,'short_center_metrics.csv'));
report=struct('status','PARTIAL_RECORD_PREVIEW_NOT_FINAL_220S_RESULT','snapshot',snapshot,'source_settings',settings,'amplitude_ratios',ratio,'eligible',eligible,'probe_reports',{records},'window_sensitivity',windowAudit,'processing_seconds_by_probe',timers, ...
 'parent_band_rule','Frequency bins within original eligible first-order frequency range; kh>=.3 and sums below Nyquist. Same band for observed first harmonic and directional prior.', ...
 'target_filter_rule','Raw results retained; common_sum_band applies exactly the same Fourier sum-support mask to reference and all predictions.', ...
 'directional_assumption','Initial directional complex distribution retained per frequency, constrained by observed first harmonic; no higher-order reference used in allocation');
writejson(fullfile(out,'report.json'),report);save(fullfile(out,'comparison.mat'),'t','raw','first','reference','results','records','report','shortResult','shortReport','-v7.3');
f=figure('Visible','off','Position',[40 40 1500 1000]);tiledlayout(3,2,'TileSpacing','compact','Padding','compact');
nexttile;plot(t,first(:,1),'Color',[.6 .6 .6]);hold on;plot(t,full.input,'b');xline(full.limits(1),':');xline(full.limits(2),':');grid on;xlabel('Time since startup (s)');ylabel('First harmonic (m)');title(sprintf('Frozen partial record: 0--%.1f s',t(end)));legend('Separated','Common input','Location','northwest','Box','off');
for probe=1:5
 nexttile;if ~eligible(probe),text(.1,.5,'Below one-third amplitude gate');axis off;continue;end
 d=results{probe};plot(t,d.filtered(:,1),'k-',t,d.filtered(:,3),'r--',t,d.filtered(:,2),'b:','LineWidth',1.1);grid on;xlim(d.limits);bound=1.1*max(abs(d.filtered(:,1:3)),[],'all');ylim([-bound bound]);xlabel('Time since startup (s)');ylabel('Second harmonic (m)');
 selected=allMetrics.probe==probe&strcmp(allMetrics.filter,'common_sum_band')&strcmp(allMetrics.window,'main_group')&strcmp(allMetrics.method,'Joint 7.5 deg');err=allMetrics.relative_L2(selected);
 title(sprintf('y offset %.2f m | GL 7.5 deg L2 %.2f%%',settings.probes_xy_m(probe,2)-settings.probes_xy_m(1,2),100*err));if probe==1,legend('HOS reference','GL 7.5 deg','GL 15 deg','Location','northwest','Box','off');end
end
sgtitle('Directional GL time-series reconstruction — partial HOS record, common output filtering');exportgraphics(f,fullfile(out,'partial_comparison.png'),'Resolution',160);exportgraphics(f,fullfile(out,'partial_comparison.pdf'),'ContentType','vector');close(f);disp(allMetrics);disp(windowAudit);
end
function [d,r,metrics]=one(realInput,ref,t,initial,s,probe)
N=numel(t);dt=mean(diff(t));tr=t-t(1);g=s.g;h=s.h;kp=s.kp;
k0=hypot(initial.kx(:),initial.ky(:));eligible=k0*h>=.3&initial.kx(:)>0;w0=initial.om(eligible);w0=w0(:);kx0=initial.kx(eligible);kx0=kx0(:);ky0=initial.ky(eligible);ky0=ky0(:);C=initial.C(eligible);C=C(:);xy=s.probes_xy_m(probe,:);A0=C.*exp(1i*(kx0*xy(1)+ky0*xy(2)));
allbins=(1:floor((N-1)/2))';allw=2*pi*allbins/(N*dt);keep=allw>=min(w0)&allw<=max(w0)&2*allw<pi/dt;bins=allbins(keep);w=allw(keep);k=zeros(size(w));
for j=1:numel(w),k(j)=fzero(@(z)g*z*tanh(z*h)-w(j)^2,[0,max(1,2*w(j)^2/g+2/h)]);end
keep=k*h>=.3;bins=bins(keep);w=w(keep);k=k(keep);F=fft(realInput)/N;observed=2*conj(F(bins+1));input=real(exp(-1i*tr*w.')*observed);angles=atan2(ky0,kx0)*180/pi;pred=zeros(N,4);condition=cell(1,2);
widths=[15 7.5];
for resolution=1:2
 width=widths(resolution);centers=(-90+width/2:width:90-width/2)';group=1+floor((angles+90)/width);priorTime=complex(zeros(N,numel(centers)));
 for direction=1:numel(centers)
  ids=find(group==direction);
  for start=1:256:numel(ids),ix=ids(start:min(start+255,numel(ids)));priorTime(:,direction)=priorTime(:,direction)+exp(-1i*tr*w0(ix).')*A0(ix);end
 end
 prior=ifft(priorTime,[],1);prior=prior(bins+1,:);[joint,condition{resolution}]=allocate_directional_record(prior,observed);
 assert(norm(real(exp(-1i*tr*w.')*sum(joint,2))-input)/norm(input)<1e-11);
 [K,TH]=ndgrid(k,deg2rad(centers));OM=repmat(w,1,numel(centers));pairBins=repmat(bins,1,numel(centers));
 [e,~]=gl_directional_sum_time(joint(:),OM(:),reshape(K.*cos(TH),[],1),reshape(K.*sin(TH),[],1),g,h,kp,tr,8,pairBins(:));pred(:,resolution)=real(e);
 if resolution==1,[e,~]=gl_directional_sum_time(prior(:),OM(:),reshape(K.*cos(TH),[],1),reshape(K.*sin(TH),[],1),g,h,kp,tr,8,pairBins(:));pred(:,3)=real(e);end
end
[e,~]=gl_directional_sum_time(observed,w,k,zeros(size(k)),g,h,kp,tr,8,bins);pred(:,4)=real(e);
[~,peak]=max(abs(hilbert(input)));limits=t(peak)+[-2 2]*s.Tp;main=t>=limits(1)&t<=limits(2);assert(limits(1)>=t(1)&&limits(2)<=t(end),'Main group incomplete');
signed=[0:floor(N/2),-floor(N/2):-1]';band=abs(signed)>=2*min(bins)&abs(signed)<=2*max(bins);combined=[ref,pred];filtered=real(ifft(fft(combined).*band));rows={};names={'Joint 15 deg','Joint 7.5 deg','Initial spectrum only','Direction blind'};
for filt=1:2
 if filt==1,a=combined;kind='raw';else,a=filtered;kind='common_sum_band';end
 for window=1:2
  if window==1,mask=true(N,1);label='available_full';else,mask=main;label='main_group';end
  for method=1:4,res=a(mask,method+1)-a(mask,1);rows(end+1,:)={probe,kind,label,names{method},norm(res)/norm(a(mask,1)),max(abs(res))/max(abs(a(mask,1)))};end
 end
end
metrics=cell2table(rows,'VariableNames',{'probe','filter','window','method','relative_L2','relative_Linf'});
r=struct('probe',probe,'parent_frequency_bins',bins.','frequency_band_rad_s',[min(w),max(w)],'original_frequency_range_rad_s',[min(w0),max(w0)],'first_projection_relative_L2',norm(input-realInput)/norm(realInput),'retained_positive_frequency_energy',sum(abs(F(bins+1)).^2)/sum(abs(F(allbins+1)).^2),'conditioning',{condition},'main_window_s',limits,'main_window_complete',true,'angle_refinement_relative_L2',norm(pred(main,2)-pred(main,1))/norm(pred(main,2)),'end_first_amplitude_over_peak',abs(realInput(end))/max(abs(realInput)));
d=struct('t',t,'input',input,'reference',ref,'pred',pred,'filtered',filtered,'limits',limits);
end
function writejson(file,r)
f=fopen(file,'w');fprintf(f,'%s',jsonencode(r,PrettyPrint=true));fclose(f);
end