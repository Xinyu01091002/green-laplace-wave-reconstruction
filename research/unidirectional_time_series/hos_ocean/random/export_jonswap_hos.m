function export_jonswap_hos(run)
s=jsondecode(fileread(fullfile(run,'settings.json')));d=load(fullfile(run,'inputs','initial_fields.mat'),'E','P');nx=s.grid(1);ny=s.grid(2);nr=s.MPI_ranks_per_phase;assert(mod(ny,nr)==0);
for phase=1:4
 folder=fullfile(run,'cases',sprintf('phi%03d',(phase-1)*90));assert(~isfolder(folder));mkdir(folder);mkdir(fullfile(folder,'Results'));e=d.E(:,:,phase).';p=d.P(:,:,phase).';
 for rank=0:nr-1
  f=fopen(fullfile(folder,'Results',sprintf('3d_ini_%03d.dat',rank)),'w');assert(f>=0);for j=1:67,fprintf(f,'# R4 GL JONSWAP initial field\n');end
  fprintf(f,'%-20s%12.5E%4s%5d%4s%5d\n','ZONE SOLUTIONTIME = ',0,', I=',nx,', J=',ny/nr);cols=rank*(ny/nr)+(1:ny/nr);ee=e(:,cols);pp=p(:,cols);fprintf(f,'%.17e %.17e\n',[ee(:),pp(:)].');fclose(f);
 end
 f=fopen(fullfile(folder,'prob.inp'),'w');fprintf(f,'%.17g %.17g\n',s.probes_xy_m.');fclose(f);
 f=fopen(fullfile(folder,'input.yml'),'w');
 fprintf(f,'restart: disable\ngravity: %.17g\ndomain size:\n  x: %.17g\n  y: %.17g\ncomputation case:\n  3D: true\n  duration: %.17g\n  type: Initial surface quantities\n',s.g,s.domain_m(1),s.domain_m(2),s.duration_s);
 fprintf(f,'numerical parameters:\n  time:\n    integration tolerance: 1.e-12\n    relative tolerance: false\n    Dommermuth initialisation:\n      n: 4\n      Ta: 0.0\n  discretization:\n    x: %d\n    y: %d\n  surface nonlinearity order: 5\n  dealiasing:\n    x: 3\n    y: 3\nbathymetry:\n  depth: %.17g\noutput:\n  directory: Results\n  dimensional: true\n  frequency: 5.0\n  max number: 5\n  free surface:\n    physical space: false\n    modal space: false\n  probes:\n    activate: true\n',nx,ny,s.h);fclose(f);
end
end
