function prepare_hos_random_double_domain(sourceRoot,targetRoot)
% Repeat the exact 68-lambda initial field twice on a 136-lambda domain.
source=fullfile(sourceRoot,'seed20260927','akp012');d=load(fullfile(source,'initial.mat'),'E','P','report');r=d.report;
assert(r.N==4096);N=2*r.N;L=2*r.L;E=[d.E;d.E];P=[d.P;d.P];assert(size(E,1)==N);
mkdir(targetRoot);
for phase=1:4
    folder=fullfile(targetRoot,sprintf('phase%03d',(phase-1)*90));mkdir(folder);mkdir(fullfile(folder,'Results'));
    fid=fopen(fullfile(folder,'Results','3d_ini.dat'),'w');
    for line=1:67,fprintf(fid,'# repeated random-phase initial field on doubled periodic domain; line %d\n',line);end
    fprintf(fid,'%-20s%12.5E%4s%5d%4s%5d\n','ZONE SOLUTIONTIME = ',0,', I=',N,', J=',1);
    fprintf(fid,'%.17e %.17e\n',[E(:,phase),P(:,phase)].');fclose(fid);
    fid=fopen(fullfile(folder,'input.yml'),'w');
    fprintf(fid,'restart: disable\ngravity: %.17g\ndomain size:\n  x: %.17g\ncomputation case:\n  3D: false\n  duration: %.17g\n  type: Initial surface quantities\n',r.g,L,50*r.Tp);
    fprintf(fid,'numerical parameters:\n  time:\n    integration tolerance: 1.e-10\n    relative tolerance: true\n    Dommermuth initialisation:\n      n: 4\n      Ta: 0.0\n  discretization:\n    x: %d\n  surface nonlinearity order: 5\n  dealiasing:\n    x: 5\nbathymetry:\n  depth: %.17g\noutput:\n  directory: Results\n  dimensional: true\n  frequency: %.17g\n  free surface:\n    physical space: false\n  probes:\n    activate: true\n  max number: 2\n',N,r.h,40/r.Tp);fclose(fid);
    fid=fopen(fullfile(folder,'prob.inp'),'w');fprintf(fid,'%.17g\n%.17g\n',r.probe_x,r.probe_x+r.L);fclose(fid);
end
report=struct('source',source,'seed',r.random_seed,'nominal_Akp',r.nominal_Akp,'actual_kp_Hs_over_2',r.actual_kp_Hs_over_2, ...
    'original_domain_m',r.L,'doubled_domain_m',L,'original_N',r.N,'doubled_N',N, ...
    'same_dx',true,'construction','exact repeated halves','probe_x',[r.probe_x,r.probe_x+r.L], ...
    'duration_Tp',50,'output_dt',r.output_dt);
fid=fopen(fullfile(targetRoot,'design.json'),'w');fprintf(fid,'%s\n',jsonencode(report));fclose(fid);disp(report);
end
