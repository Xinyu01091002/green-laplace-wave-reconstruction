function audit_record_direction(run,out)
assert(~isfolder(out));mkdir(out);
src=fullfile(run,'random-gl-comparison-v1');addpath(fullfile(src,'source','research','directional_wave_data'));addpath(fullfile(src,'source','research','unidirectional_time_series'));
s=jsondecode(fileread(fullfile(run,'settings.json')));ini=load(fullfile(run,'inputs','initial_fields.mat'),'C','kx','ky','om');old=load(fullfile(src,'comparison.mat'),'raw','t','results');
t=old.t;raw=old.raw;fullFirst=separate(raw);rows={};cases={};
% Predeclared angular refinements at unchanged temporal resolution, all probes.
for p=1:5
 for width=[15 7.5 3.75 1.875]
  add(p,1101,width,0,false);
 end
end
% Center-only temporal diagnostics; same absolute scoring interval throughout.
for n=[1051 1001]
 add(1,n,1.875,0,false);add(1,n,1.875,0,true);
end
% Cosine endpoint taper of first-order input and directional prior. Since the
% reconstruction is quadratic, reference is multiplied by w^2. In the fixed
% scoring interval w=1. These are diagnostics, not an exact window covariance.
for taper=[10 20]
 for n=[1101 1051 1001],add(1,n,1.875,taper,false);end
end
T=cell2table(rows,'VariableNames',{'probe','N','end_s','angle_deg','taper_s','full_Hilbert_prefix','raw_L2','common_band_L2','vs_full_fine_prediction_L2','condition_weighted','condition_max','seconds'});
writetable(T,fullfile(out,'metrics.csv'));
baseline=cases{find(T.probe==1&T.N==1101&T.angle_deg==7.5&T.taper_s==0&~T.full_Hilbert_prefix,1)};
parity=norm(baseline.pred-old.results{1}.pred(:,2))/norm(old.results{1}.pred(:,2));assert(parity<1e-10,'Original baseline changed');
summary=struct('status','completed','baseline_prediction_parity',parity,'fixed_window_s',s.scoring_window_s,'taper_rule','First order input and directional prior w; quadratic reference w^2; w=1 in scored interval. Same Fourier output mask applied to reference and prediction. Diagnostic only, not a physical correction.','notes','Short records change Fourier grid and support. Full-Hilbert-prefix cases isolate redoing phase separation from later frequency-grid changes. No fitting or HOS rerun.');
f=fopen(fullfile(out,'report.json'),'w');fprintf(f,'%s',jsonencode(summary,PrettyPrint=true));fclose(f);
save(fullfile(out,'audit.mat'),'T','cases','summary','-v7.3');
f=figure('Visible','off','Position',[30 30 1300 800]);tiledlayout(2,2);
nexttile;hold on;for p=1:5,ix=T.probe==p&T.N==1101&T.taper_s==0;plot(T.angle_deg(ix),100*T.common_band_L2(ix),'-o');end;set(gca,'XDir','reverse');grid on;xlabel('Direction bin width (deg)');ylabel('Relative L2 error (%)');legend('center','y-52.8','y+52.8','y-70.4','y+70.4','Location','best');
nexttile;hold on;for taper=[0 10 20],ix=T.probe==1&T.angle_deg==1.875&T.taper_s==taper&~T.full_Hilbert_prefix;plot(T.end_s(ix),100*T.common_band_L2(ix),'-o');end;grid on;xlabel('Record end (s)');ylabel('Relative L2 error (%)');legend('Untapered','10 s taper','20 s taper','Location','best');
nexttile;hold on;for width=[7.5 1.875],ix=find(T.probe==1&T.N==1101&T.angle_deg==width&T.taper_s==0);d=cases{ix};plot(d.t,d.filtered(:,2));end;plot(d.t,d.filtered(:,1),'k--');xlim(s.scoring_window_s);grid on;legend('GL 7.5 deg','GL 1.875 deg','HOS sector');xlabel('Time (s)');ylabel('Second harmonic (m)');
nexttile;hold on;for taper=[0 10 20],ix=T.probe==1&T.angle_deg==1.875&T.taper_s==taper&~T.full_Hilbert_prefix;plot(T.end_s(ix),100*T.vs_full_fine_prediction_L2(ix),'-o');end;grid on;xlabel('Record end (s)');ylabel('Change from full fine prediction (%)');legend('Untapered','10 s taper','20 s taper','Location','best');
exportgraphics(f,fullfile(out,'audit.png'),'Resolution',160);close(f);disp(T);
 function add(p,n,width,taper,prefix)
  timer=tic;d=one(raw(1:n,p,:),t(1:n),ini,s,p,width,taper,prefix,fullFirst(1:n,p));
  match=find(cellfun(@(x)x.p==p&&x.n==1101&&x.width==1.875&&x.taper==0&&~x.prefix,cases),1);
  delta=NaN;if ~isempty(match),b=cases{match};delta=norm(d.pred(d.mask)-b.pred(find(d.mask)))/norm(b.pred(find(d.mask)));end
  if n==1101&&width==1.875&&taper==0&&~prefix,delta=0;end
  d.p=p;d.n=n;d.width=width;d.taper=taper;d.prefix=prefix;cases{end+1}=d;
  rows(end+1,:)={p,n,t(n),width,taper,prefix,d.rawError,d.filteredError,delta,d.condition.observed_energy_weighted,d.condition.max,toc(timer)};
  fprintf('DONE p%d N%d angle%g taper%g prefix%d error %.4f%%\n',p,n,width,taper,prefix,100*d.filteredError);
 end
end
function first=separate(raw)
H=imag(hilbert(raw));first=(raw(:,:,1)-raw(:,:,3)-H(:,:,2)+H(:,:,4))/4;
end
function d=one(raw,t,ini,s,p,width,taper,prefix,fullFirst)
N=numel(t);dt=mean(diff(t));tr=t-t(1);win=ones(N,1);
if taper>0,u=min(t-t(1),t(end)-t);sel=u<taper;win(sel)=.5*(1-cos(pi*u(sel)/taper));end
first=separate(raw.*reshape(win,[],1,1));if prefix,assert(taper==0);first=fullFirst;end
ref=(raw(:,1,1)-raw(:,1,2)+raw(:,1,3)-raw(:,1,4))/4.*win.^2;
use=hypot(ini.kx(:),ini.ky(:))*s.h>=.3&ini.kx(:)>0;
w0=ini.om(use);w0=w0(:);kx0=ini.kx(use);kx0=kx0(:);ky0=ini.ky(use);ky0=ky0(:);C=ini.C(use);C=C(:);xy=s.probes_xy_m(p,:);A=C.*exp(1i*(kx0*xy(1)+ky0*xy(2)));
bins=(1:floor((N-1)/2))';w=2*pi*bins/(N*dt);keep=w>=min(w0)&w<=max(w0)&2*w<pi/dt;bins=bins(keep);w=w(keep);k=zeros(size(w));
for j=1:numel(w),k(j)=fzero(@(z)s.g*z*tanh(z*s.h)-w(j)^2,[0,max(1,2*w(j)^2/s.g+2/s.h)]);end
keep=k*s.h>=.3;bins=bins(keep);w=w(keep);k=k(keep);F=fft(first)/N;observed=2*conj(F(bins+1));
centers=(-90+width/2:width:90-width/2)';group=1+floor((atan2(ky0,kx0)*180/pi+90)/width);priorTime=complex(zeros(N,numel(centers)));
for direction=1:numel(centers)
 ids=find(group==direction);for start=1:256:numel(ids),ix=ids(start:min(start+255,numel(ids)));priorTime(:,direction)=priorTime(:,direction)+exp(-1i*tr*w0(ix).')*A(ix);end
end
prior=ifft(priorTime.*win,[],1);prior=prior(bins+1,:);[joint,condition]=allocate_directional_record(prior,observed);
[K,TH]=ndgrid(k,deg2rad(centers));OM=repmat(w,1,numel(centers));pairBins=repmat(bins,1,numel(centers));
[e,~]=gl_directional_sum_time(joint(:),OM(:),reshape(K.*cos(TH),[],1),reshape(K.*sin(TH),[],1),s.g,s.h,s.kp,tr,8,pairBins(:));pred=real(e);
mask=t>=s.scoring_window_s(1)&t<=s.scoring_window_s(2);assert(all(win(mask)==1));signed=[0:floor(N/2),-floor(N/2):-1]';band=abs(signed)>=2*min(bins)&abs(signed)<=2*max(bins);filtered=real(ifft(fft([ref,pred]).*band));
d=struct('t',t,'pred',pred,'reference',ref,'filtered',filtered,'mask',mask,'condition',condition,'rawError',norm(pred(mask)-ref(mask))/norm(ref(mask)),'filteredError',norm(filtered(mask,2)-filtered(mask,1))/norm(filtered(mask,1)));
end
