function benchmark_r4_mf12(root,stage)
addpath(fullfile(root,'spark'));addpath(fullfile(root,'mf12'));
switch stage
 case 'prepare'
  source='/home/lxy/green-laplace-unidirectional-time-series-runs/hos-random-fourphase-20260926-v1';
  s=load(fullfile(source,'jonswap-design-v1','linear_design.mat'));sel=load(fullfile(source,'jonswap-energy-support-v1','retained_indices.mat'));
  n=find(sel.cum>=.99,1);assert(n==11128);ii=sel.order(1:n);ids=find(s.use);ids=ids(ii);C=s.C(ii);kp=s.report.kp;h=s.report.h;g=9.81;
  targetHs=2*.12/kp;C=C*targetHs/(4*sqrt(sum(abs(C).^2)/2));kx=s.kx;ky=s.ky;[ny,nx]=size(kx);Lx=s.report.domain_m(1);Ly=s.report.domain_m(2);
  % Real FFT input to SPARK is half the analytic amplitude on each conjugate side.
  F=complex(zeros(ny,nx));F(ids)=nx*ny*C/2;F=F+conj(F([1 ny:-1:2],[1 nx:-1:2]));eta11=real(ifft2(F));
  assert(abs(4*std(eta11(:),1)-targetHs)/targetHs<1e-12);assert(nnz(F)==2*n);
  c1=mf12_spectral_coefficients(1,g,h,real(C),imag(C),kx(ids),ky(ids),0,0);
  [e1,~]=mf12_spectral_surface(c1,Lx,Ly,nx,ny,0);parity=norm(e1(:)-eta11(:))/norm(eta11(:));assert(parity<1e-12);
  input=struct('parents',n,'seed',20260925,'kpHs_over_2',.12,'Hs',targetHs,'grid',[nx ny],'domain_m',[Lx Ly],'h',h,'g',g,'MF12_linear_input_parity',parity,'selected_energy',sel.cum(n),'scope','Same 99 percent energy support, dimensional raw nonzero-difference eta20 and true-surface psi20; no fitting/filtering');
  save(fullfile(root,'input.mat'),'F','C','ids','kx','ky','kp','h','g','nx','ny','Lx','Ly','input','-v7.3');writejson(fullfile(root,'input.json'),input);
 case {'r4_eta','r4_psi'}
  s=load(fullfile(root,'input.mat'));tic;
  try
   if strcmp(stage,'r4_eta'),[half,audit]=spark_eta20_r_series(s.F,s.kx,s.ky,s.h,4);field=2*half;
   else,[half,audit]=spark_psi20_r_series(s.F,s.kx,s.ky,s.h,4);field=2*sqrt(s.g)*half;end
   seconds=toc;report=struct('status','returned_successfully','operator_seconds',seconds,'audit',audit);
   save(fullfile(root,[stage '.mat']),'field','report','-v7.3');
  catch ME
   report=struct('status','rejected_by_original_operator','operator_seconds',toc,'identifier',ME.identifier,'message',ME.message);
  end
  writejson(fullfile(root,[stage '.json']),report);disp(report);
 case 'mf12'
  s=load(fullfile(root,'input.mat'));tic;
  c=mf12_spectral_coefficients(2,s.g,s.h,real(s.C),imag(s.C),s.kx(s.ids),s.ky(s.ids),0,0);coefficient_seconds=toc;
  fprintf('Coefficient generation completed: %.3f s\n',coefficient_seconds);
  % Published MF12 decomposition: suppress 11/self/sum, retain all difference slots.
  tic;c.a(:)=0;c.b(:)=0;c.muStar(:)=0;c.G_2(:)=0;c.mu_2(:)=0;c.G_npm(1:2:end)=0;c.mu_npm(1:2:end)=0;selection_seconds=toc;
  tic;[eta20,psi20]=mf12_spectral_surface(c,s.Lx,s.Ly,s.nx,s.ny,0);surface_seconds=toc;
  assert(all(isfinite(eta20),'all')&&all(isfinite(psi20),'all'));
  report=struct('status','completed','coefficient_seconds',coefficient_seconds,'selection_seconds',selection_seconds,'surface_seconds',surface_seconds,'total_operator_seconds',coefficient_seconds+selection_seconds+surface_seconds,'timing_scope','Unmodified MF12 order2 coefficient API computes sum and difference; surface called once for difference eta/psi only. No MF12 third order.');
  clear c;save(fullfile(root,'mf12.mat'),'eta20','psi20','report','-v7.3');writejson(fullfile(root,'mf12.json'),report);disp(report);
 case 'summary'
  s=load(fullfile(root,'input.mat'));ref=load(fullfile(root,'mf12.mat'));rows={};
  for name={'eta','psi'}
   name=name{1};label=['r4_' name];info=jsondecode(fileread(fullfile(root,[label '.json'])));
   if strcmp(info.status,'returned_successfully')
    a=load(fullfile(root,[label '.mat']));b=ref.([name '20']);delta=a.field-b;
    rows(end+1,:)={name,info.status,info.operator_seconds,norm(delta(:))/norm(b(:)),max(abs(delta),[],'all')/max(abs(b),[],'all'),sqrt(mean(delta.^2,'all')),norm(a.field(:))/norm(b(:))};
   else,rows(end+1,:)={name,info.status,info.operator_seconds,NaN,NaN,NaN,NaN};end
  end
  T=cell2table(rows,'VariableNames',{'field','R4_status','R4_operator_seconds','relative_L2','relative_Linf','absolute_RMS','norm_ratio'});writetable(T,fullfile(root,'comparison.csv'));
  writejson(fullfile(root,'comparison.json'),struct('input',s.input,'metrics',table2struct(T),'MF12',ref.report,'note','A rejected R4 result is not a validated field and has no accepted accuracy score. Individual process wall/RSS recorded separately.'));disp(T);
end
end
function writejson(path,data)
f=fopen(path,'w');assert(f>=0);fprintf(f,'%s',jsonencode(data,PrettyPrint=true));fclose(f);
end
