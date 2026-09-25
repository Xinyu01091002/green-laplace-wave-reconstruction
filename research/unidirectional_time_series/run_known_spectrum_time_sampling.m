function report = run_known_spectrum_time_sampling()
% Existing spatial GL API sampled in time; not a time-series inversion.
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root); setup_green_laplace;
g=9.81; h=1; a=[0.010,0.006]; b=[0,0.002]; k=[2,3];
L=2*pi; nx=32; ny=4; times=[0,0.125,0.25,0.5];
opts=struct('eta22_rank',6);
c=gl_spectral_coefficients(3,g,h,a,b,k,[0,0],0,0,opts);
[~,~,~,~,baseline]=gl_spectral_surface(c,L,L,nx,ny,0);
names={'eta11','psi11','eta22','psi22','eta33','psi33'};
relL2=zeros(numel(times),6); relLinf=relL2; trace=relL2;
linear_error=zeros(numel(times),2); y_error=relL2;
for it=1:numel(times)
    t=times(it);
    [~,~,X,~,parts,audit]=gl_spectral_surface(c,L,L,nx,ny,t);
    A=complex(a,b).*exp(-1i*c.omega*t);
    shifted=gl_spectral_coefficients(3,g,h,real(A),imag(A),k,[0,0],0,0,opts);
    [~,~,~,~,reference]=gl_spectral_surface(shifted,L,L,nx,ny,0);
    for j=1:numel(names)
        v=parts.(names{j}); r=reference.(names{j});
        if t==0, assert(isequal(v,baseline.(names{j}))); end
        assert(all(isfinite(v),'all'));
        relL2(it,j)=norm(v(:)-r(:))/max(norm(r(:)),realmin);
        relLinf(it,j)=max(abs(v(:)-r(:)))/max(max(abs(r(:))),realmin);
        y_error(it,j)=max(abs(v-repmat(v(1,:),ny,1)),[],'all');
        trace(it,j)=v(1,1);
    end
    e=zeros(size(X)); p=e;
    for j=1:numel(k)
        wave=A(j)*exp(1i*k(j)*X);
        e=e+real(wave); p=p+real(-1i*g/c.omega(j)*wave);
    end
    linear_error(it,:)=[norm(parts.eta11(:)-e(:))/norm(e(:)), ...
        norm(parts.psi11(:)-p(:))/norm(p(:))];
    assert(audit.time==t && ~audit.rescaling_alignment_or_gain);
end
assert(max(relL2,[],'all')<1e-12 && max(relLinf,[],'all')<1e-12);
assert(max(linear_error,[],'all')<1e-12 && max(y_error,[],'all')<1e-14);
report=struct('pass',true,'matlab_release',version('-release'), ...
    'times_s',times,'component_names',{names},'relative_L2',relL2, ...
    'relative_Linf',relLinf,'linear_relative_L2',linear_error, ...
    'y_absolute_error',y_error,'t0_exact',true,'fixed_point_trace',trace, ...
    'coefficients',c,'grid',[ny,nx],'domain',[L,L],'position',[0,0]);
out=fullfile(root,'results','unidirectional_time_series');
if ~isfolder(out), mkdir(out); end
save(fullfile(out,'known_spectrum_baseline.mat'),'report');
fid=fopen(fullfile(out,'known_spectrum_baseline.json'),'w');
assert(fid>=0); cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(report));
fprintf('PASS time sampling: t0 exact; max L2 %.3e; Linf %.3e; linear %.3e\n', ...
    max(relL2,[],'all'),max(relLinf,[],'all'),max(linear_error,[],'all'));
end
