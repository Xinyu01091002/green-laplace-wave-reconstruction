function prepare_hos_precision(inputFile,out)
assert(~isfolder(out),'Output already exists');
s=load(inputFile,'E','P','report'); d=s.report.design;
assert(isequal(size(s.E),[256,1024,4]));
mkdir(out); mkdir(fullfile(out,'Results'));
eta=s.E(:,:,1).'; psi=s.P(:,:,1).'; nx=size(eta,1); ny=size(eta,2);
f=fopen(fullfile(out,'Results','3d_ini.dat'),'w'); assert(f>=0);
for j=1:67, fprintf(f,'# External MF12 order-2 initial surface; header line %d\n',j); end
fprintf(f,'%-20s%12.5E%4s%5d%4s%5d\n','ZONE SOLUTIONTIME = ',0,', I=',nx,', J=',ny);
fprintf(f,'%.17e %.17e\n',[eta(:),psi(:)].'); fclose(f);
f=fopen(fullfile(out,'Results','3d_ini.dat')); for j=1:68,fgetl(f);end
q=fscanf(f,'%f',[2,Inf]); fclose(f); assert(isequal(size(q),[2,nx*ny]));
r=struct('input_file',inputFile,'phase_degrees',0,'grid',[nx,ny], ...
    'domain_m',d.domain_m,'depth_m',d.h,'g',d.g, ...
    'quantization_eta_relative_L2',norm(q(1,:).'-eta(:))/norm(eta(:)), ...
    'quantization_psi_relative_L2',norm(q(2,:).'-psi(:))/norm(psi(:)), ...
    'quantization_eta_max_abs_m',max(abs(q(1,:).'-eta(:))), ...
    'quantization_psi_max_abs_m2s',max(abs(q(2,:).'-psi(:))), ...
    'M',5,'qx',3,'qy',3,'duration_s',2,'output_interval_s',.2,'tolerance',1e-8, ...
    'ramp_s',0,'breaking',false,'dissipation',false,'periodic',true);
f=fopen(fullfile(out,'input_HOS-Ocean.yml'),'w');
fprintf(f,'restart: disable\ngravity: %.17g\ndomain size:\n  x: %.17g\n  y: %.17g\n',d.g,d.domain_m(1),d.domain_m(2));
fprintf(f,'computation case:\n  3D: true\n  duration: 2.0\n  type: Initial surface quantities\n');
fprintf(f,'numerical parameters:\n  time:\n    integration tolerance: 1.e-8\n    Dommermuth initialisation:\n      n: 4\n      Ta: 0.0\n  discretization:\n    x: %d\n    y: %d\n  surface nonlinearity order: 5\n  dealiasing:\n    x: 3\n    y: 3\n',nx,ny);
fprintf(f,'bathymetry:\n  depth: %.17g\noutput:\n  directory: Results\n  dimensional: true\n  frequency: 5.0\n  HDF5: false\n  free surface:\n    physical space: true\n    modal space: false\n',d.h); fclose(f);
assert(r.quantization_eta_relative_L2<1e-15 && r.quantization_psi_relative_L2<1e-15);
r.input_format='list-directed 18 significant digits';
f=fopen(fullfile(out,'input-audit.json'),'w');fprintf(f,'%s',jsonencode(r,PrettyPrint=true));fclose(f);
save(fullfile(out,'import_reference.mat'),'eta','psi','q','r','-v7.3');disp(r);
end