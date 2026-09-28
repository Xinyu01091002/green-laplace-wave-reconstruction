function prepare_hos_high_xy_probe_run(root,source_root)
%PREPARE_HOS_HIGH_XY_PROBE_RUN Freeze Akp=.12 off-axis probe replay.
assert(~isfolder(fullfile(root,'cases')));mkdir(fullfile(root,'cases'));
source=fullfile(source_root,'inputs','initial_fields.mat');s=load(source);d=s.report.design;
assert(isequal(size(s.E),[256 1024 4])&&isequal(size(s.P),size(s.E)));
nx=1024;ny=256;Lx=d.domain_m(1);Ly=d.domain_m(2);
centre=[493,129;503,129;523,129;533,129];
offaxis=[493,119;493,135;503,123;503,132;523,126;523,135;533,123;533,139];
probe_index=[centre;offaxis];probes=(probe_index-1).*[Lx/nx,Ly/ny];
focus_xy=[Lx/2,Ly/2];offset_lambda=(probes-focus_xy)/(2*pi/d.kp);
centre_for_offaxis=[1;1;2;2;3;3;4;4];
phases=[0,90,180,270];expected=zeros(4,size(probe_index,1));
for phase_index=1:4
    phase=phases(phase_index);folder=fullfile(root,'cases',sprintf('phi%03d',phase));
    mkdir(folder);mkdir(fullfile(folder,'Results'));
    eta=s.E(:,:,phase_index).';psi=s.P(:,:,phase_index).';
    for rank=0:3
        columns=rank*64+(1:64);
        file=fopen(fullfile(folder,'Results',sprintf('3d_ini_%03d.dat',rank)),'w');assert(file>=0);
        for header=1:67
            fprintf(file,'# Frozen Akp=.12 MF12 order-2 input, phase %d, header %d\n',phase,header);
        end
        fprintf(file,'%-20s%12.5E%4s%5d%4s%5d\n','ZONE SOLUTIONTIME = ',0,', I=',nx,', J=',64);
        eta_rank=eta(:,columns);psi_rank=psi(:,columns);
        fprintf(file,'%.17e %.17e\n',[eta_rank(:),psi_rank(:)].');fclose(file);
    end
    file=fopen(fullfile(folder,'prob.inp'),'w');assert(file>=0);
    fprintf(file,'%.17g %.17g\n',probes.');fclose(file);
    for probe=1:size(probe_index,1)
        expected(phase_index,probe)=eta(probe_index(probe,1),probe_index(probe,2));
    end
    file=fopen(fullfile(folder,'input.yml'),'w');assert(file>=0);
    fprintf(file,'restart: disable\ngravity: %.17g\ndomain size:\n  x: %.17g\n  y: %.17g\ncomputation case:\n  3D: true\n  duration: 220.0\n  type: Initial surface quantities\n',d.g,Lx,Ly);
    fprintf(file,'numerical parameters:\n  time:\n    integration tolerance: 1.e-12\n    relative tolerance: false\n    Dommermuth initialisation:\n      n: 4\n      Ta: 0.0\n  discretization:\n    x: 1024\n    y: 256\n  surface nonlinearity order: 5\n  dealiasing:\n    x: 3\n    y: 3\nbathymetry:\n  depth: %.17g\noutput:\n  directory: Results\n  dimensional: true\n  frequency: 5.0\n  max number: %d\n  free surface:\n    physical space: false\n    modal space: false\n  probes:\n    activate: true\n',d.h,size(probe_index,1));
    fclose(file);
end

t=(0:.2:220)';use=abs(s.C)>0&s.kx>0;
C=s.C(use);kx=s.kx(use);ky=s.ky(use);om=s.om(use);
C=C(:);kx=kx(:);ky=ky(:);om=om(:);
phase=exp(1i*(kx*probes(:,1).'+ky*probes(:,2).'));
analytic=complex(zeros(numel(t),size(probes,1)));
for start=1:256:numel(C)
    ids=start:min(start+255,numel(C));
    analytic=analytic+exp(-1i*t*om(ids).')*(C(ids).*phase(ids,:));
end
linear_max=max(abs(real(analytic)),[],1);
same_x_ratio=ones(1,size(probe_index,1));
for index=1:8
    same_x_ratio(4+index)=linear_max(4+index)/linear_max(centre_for_offaxis(index));
end
assert(all(same_x_ratio(5:end)>=1/3));
Tp=2*pi/sqrt(d.g*d.kp*tanh(d.kp*d.h));
settings=struct('g',d.g,'kp',d.kp,'h',d.h,'kph',d.kp*d.h, ...
    'lambda_p',2*pi/d.kp,'Tp',Tp,'domain_m',[Lx,Ly],'grid',[nx,ny], ...
    'Akp',d.kp*sum(abs(s.C)),'duration_s',220,'output_dt_s',.2, ...
    'expected_samples',1101,'nominal_focus_xy',focus_xy, ...
    'phases_degrees',phases,'MPI_ranks_per_phase',4,'concurrent_phases',4, ...
    'cpu_sets',{{'4-7','8-11','12-15','36-39'}}, ...
    'probe_indices',probe_index,'probes_xy_m',probes, ...
    'probe_offset_lambda',offset_lambda,'centre_probe_rows',(1:4)', ...
    'centre_for_offaxis_rows',centre_for_offaxis+4*0, ...
    'offaxis_probe_rows',(5:12)','linear_first_max_abs_m',linear_max, ...
    'linear_same_x_centre_ratio',same_x_ratio, ...
    'amplitude_gate','linear full-record max(abs eta1) >= one third same-x centreline', ...
    'expected_initial_probe_eta',expected, ...
    'initialization','same frozen Akp=.12 MF12 11+20+22 initial fields as source run', ...
    'source_run',source_root);
writejson(fullfile(root,'settings.json'),settings);
save(fullfile(root,'selection.mat'),'settings','analytic','t');
disp(settings)
end

function writejson(file,value)
fid=fopen(file,'w');assert(fid>=0);fprintf(fid,'%s\n',jsonencode(value,PrettyPrint=true));fclose(fid);
end
