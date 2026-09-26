function prepare_jonswap_hybrid(root)
addpath(fullfile(root,'spark'));
source='/home/lxy/green-laplace-unidirectional-time-series-runs/r4-mf12-jonswap11128-20260926-v1/input.mat';s=load(source);
assert(s.input.parents==11128);F=s.F;h=s.h;g=s.g;kp=s.kp;nx=s.nx;ny=s.ny;
tic;[half,aeta]=spark_eta20_r_series(F,s.kx,s.ky,h,4);eta20=2*half;tEta=toc;
tic;[half,apsi]=spark_psi20_r_series(F,s.kx,s.ky,h,4);psi20=2*sqrt(g)*half;tPsi=toc;
u=zeros(size(F));u(s.kx>0)=2*F(s.kx>0)/h;prior=[];rankAudit={};
for rank=[12 16 24 32]
 timer=tic;v=spark_native_superharmonics(u,h*s.kx,h*s.ky,h*kp,rank,6,4,2,true);
 if isempty(prior),delta=[NaN NaN];else,delta=[norm(v.eta22(:)-prior.eta22(:))/norm(v.eta22(:)),norm(v.psi22(:)-prior.psi22(:))/norm(v.psi22(:))];end
 rankAudit{end+1}=struct('rank',rank,'seconds',toc(timer),'relative_change_eta_psi',delta);disp(rankAudit{end});
 if rank>=16 && all(delta<.005),break;end
 prior=v;
end
assert(all(delta<.005),'GL initialization rank changes exceed 0.5 percent');
analytic=2*F;analytic(s.kx<=0)=0;eta1=ifft2(analytic);omFull=sqrt(g*hypot(s.kx,s.ky).*tanh(h*hypot(s.kx,s.ky)));pF=zeros(size(F));pF(s.ids)=-1i*g./omFull(s.ids).*analytic(s.ids);psi1=ifft2(pF);
eta22=h*v.eta22;psi22=h*sqrt(g*h)*v.psi22;
Tp=2*pi/sqrt(g*kp*tanh(kp*h));duration=.4*round(80*Tp/.4);dt=.2;count=round(duration/dt)+1;
audit=struct('initialization','11 analytic + pure R4 20 + pure GL 22; no initial31/33','R4_eta_seconds',tEta,'R4_psi_seconds',tPsi,'R4_eta_audit',aeta,'R4_psi_audit',apsi,'GL_rank_audit',{rankAudit},'GL_selected_rank',rank,'R4_numerical_check_policy','Unchanged source formulas and repair; defects 1e-8..1e-7 warn, greater than1e-7 stop. Explicit run-local policy for previously inspected ~1e-8 residuals; original SPARK untouched.','amplitude_scaling','High and low share C phases, linear ratio1/6 and quadratic ratio1/36');
for level=[.02 .12]
 if level==.02,name='low';else,name='high';end
 run=fullfile(root,name);mkdir(run);mkdir(fullfile(run,'inputs'));a=level/.12;C=a*s.C;kx=s.kx(s.ids);ky=s.ky(s.ids);om=omFull(s.ids);
 E=zeros(ny,nx,4);P=E;for j=1:4,phase=(j-1)*pi/2;E(:,:,j)=real(a*eta1*exp(1i*phase)+a*a*eta22*exp(2i*phase))+a*a*eta20;P(:,:,j)=real(a*psi1*exp(1i*phase)+a*a*psi22*exp(2i*phase))+a*a*psi20;end
 assert(all(isfinite(E),'all')&&all(isfinite(P),'all'));assert(max(abs(mean(E,[1 2])),[],'all')<1e-9);
 first=(E(:,:,1)-1i*E(:,:,2)-E(:,:,3)+1i*E(:,:,4))/2;assert(norm(first-a*eta1,'fro')/norm(a*eta1,'fro')<1e-11);
 save(fullfile(run,'inputs','initial_fields.mat'),'C','kx','ky','om','E','P','audit','-v7.3');
 probeIndex=[513 257;513 251;513 263;513 249;513 265];probes=(probeIndex-1).*[s.Lx/nx,s.Ly/ny];expected=zeros(4,5);
 for j=1:4,for p=1:5,expected(j,p)=E(probeIndex(p,2),probeIndex(p,1),j);end,end
 settings=struct('g',g,'kp',kp,'h',h,'kph',kp*h,'Tp',Tp,'lambda_p',2*pi/kp,'grid',[nx ny],'domain_m',[s.Lx s.Ly],'domain_lambda_p',[50 20],'kpHs_over_2',level,'Hs_linear',2*level/kp,'parents',numel(C),'seed',20260925,'jonswap_gamma',3.3,'retained_energy',s.input.selected_energy,'direction_energy_sigma_deg',25/sqrt(2),'direction_mean_deg',0,'initialization',audit.initialization,'GL_initialization_rank',rank,'duration_s',duration,'duration_Tp',duration/Tp,'output_dt_s',dt,'expected_samples',count,'MPI_ranks_per_phase',8,'concurrent_phases',4,'probe_indices',probeIndex,'probes_xy_m',probes,'expected_initial_probe_eta',expected,'M',5,'dealiasing',[3 3],'absolute_tolerance',1e-12,'ramp_Ta_s',0,'scoring_window_s',[10 70]*Tp,'boundary','periodic, no absorbing layer','field_output','probes only; full initial eta/psi retained','all_probe_vars','time plus five eta probes');
 writejson(fullfile(run,'settings.json'),settings);export_jonswap_hos(run);
end
writejson(fullfile(root,'initialization-audit.json'),audit);
end
function writejson(p,data)
f=fopen(p,'w');assert(f>=0);fprintf(f,'%s',jsonencode(data,PrettyPrint=true));fclose(f);
end
