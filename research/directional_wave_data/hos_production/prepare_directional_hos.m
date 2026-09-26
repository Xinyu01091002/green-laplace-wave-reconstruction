function prepare_directional_hos(root,source)
assert(~isfolder(fullfile(root,'cases')));mkdir(fullfile(root,'cases'));s=load(source);d=s.report.design;
assert(isequal(size(s.E),[256 1024 4])&&isequal(size(s.P),size(s.E)));assert(abs(d.kp*d.h-1)<1e-12);
etaFirst=(s.E(:,:,1)-1i*s.E(:,:,2)-s.E(:,:,3)+1i*s.E(:,:,4))/2;psiFirst=(s.P(:,:,1)-1i*s.P(:,:,2)-s.P(:,:,3)+1i*s.P(:,:,4))/2;
assert(norm(etaFirst-s.eta1,'fro')/norm(s.eta1,'fro')<1e-10);assert(norm(psiFirst-s.psi1,'fro')/norm(s.psi1,'fro')<1e-10);
nx=1024;ny=256;Lx=d.domain_m(1);Ly=d.domain_m(2);Tp=2*pi/sqrt(d.g*d.kp*tanh(d.kp*d.h));
probeIndex=[513 129;513 126;513 132;513 125;513 133];probes=(probeIndex-1).*[Lx/nx,Ly/ny];expected=zeros(4,5);phases=[0 90 180 270];
for p=1:4
    folder=fullfile(root,'cases',sprintf('phi%03d',phases(p)));mkdir(folder);mkdir(fullfile(folder,'Results'));e=s.E(:,:,p).';ps=s.P(:,:,p).';assert(all(isfinite(e),'all')&&all(isfinite(ps),'all'));
    for rank=0:3
        cols=rank*64+(1:64);f=fopen(fullfile(folder,'Results',sprintf('3d_ini_%03d.dat',rank)),'w');assert(f>=0);
        for n=1:67,fprintf(f,'# Frozen independent MF12 order-2 input, phase %d, header %d\n',phases(p),n);end
        fprintf(f,'%-20s%12.5E%4s%5d%4s%5d\n','ZONE SOLUTIONTIME = ',0,', I=',nx,', J=',64);
        er=e(:,cols);pr=ps(:,cols);fprintf(f,'%.17e %.17e\n',[er(:),pr(:)].');fclose(f);
    end
    f=fopen(fullfile(folder,'prob.inp'),'w');fprintf(f,'%.17g %.17g\n',probes.');fclose(f);
    for j=1:5,expected(p,j)=e(probeIndex(j,1),probeIndex(j,2));end
    f=fopen(fullfile(folder,'input.yml'),'w');
    fprintf(f,'restart: disable\ngravity: %.17g\ndomain size:\n  x: %.17g\n  y: %.17g\ncomputation case:\n  3D: true\n  duration: 220.0\n  type: Initial surface quantities\n',d.g,Lx,Ly);
    fprintf(f,'numerical parameters:\n  time:\n    integration tolerance: 1.e-12\n    relative tolerance: false\n    Dommermuth initialisation:\n      n: 4\n      Ta: 0.0\n  discretization:\n    x: 1024\n    y: 256\n  surface nonlinearity order: 5\n  dealiasing:\n    x: 3\n    y: 3\nbathymetry:\n  depth: %.17g\noutput:\n  directory: Results\n  dimensional: true\n  frequency: 5.0\n  max number: 5\n  free surface:\n    physical space: false\n    modal space: false\n  probes:\n    activate: true\n',d.h);fclose(f);
end
r=struct('g',d.g,'kp',d.kp,'h',d.h,'kph',d.kp*d.h,'lambda_p',2*pi/d.kp,'Tp',Tp,'domain_m',[Lx Ly],'domain_lambda_p',[50 20],'grid',[nx ny],'dx_dy_m',[Lx/nx Ly/ny],'Akp',d.kp*sum(abs(s.C)),'parents',numel(s.C),'initialization','frozen MF12 order2: eta/true-surface psi 11+20+22; no33 or31','direction','positive-x mean; Gaussian AMPLITUDE angular sigma 25 degrees','radial_spectrum','semi-Gaussian, lower width .004606, upper kp/sqrt(2*log(10^8)); original frozen support','duration_s',220,'duration_Tp',220/Tp,'nominal_focus_s',110,'nominal_focus_Tp',110/Tp,'nominal_focus_xy',[Lx/2 Ly/2],'output_dt_s',.2,'expected_samples',1101,'time_integrator','Cash-Karp embedded Runge-Kutta 5(4), 6 stages, adaptive; linear part integrated analytically','absolute_tolerance',1e-12,'M',5,'dealiasing',[3 3],'padded_FFT_grid',[2048 512],'full_dealiasing',false,'ramp_Ta_s',0,'breaking',false,'added_dissipation',false,'boundary','periodic in x and y, constant finite depth','phases_degrees',phases,'MPI_ranks_per_phase',4,'concurrent_phases',4,'probe_indices',probeIndex,'probes_xy_m',probes,'expected_initial_probe_eta',expected,'all_probe_vars','time [s], eta at 5 probes [m]; no psi columns','first_sector_input_parity_eta',norm(etaFirst-s.eta1,'fro')/norm(s.eta1,'fro'));
f=fopen(fullfile(root,'settings.json'),'w');fprintf(f,'%s',jsonencode(r,PrettyPrint=true));fclose(f);disp(r);
end