function [result,audit]=gl_unidirectional_time_series(eta1,t,g,h,options)
% Single-direction eta1(t) -> eta22(t), optionally eta33(t), GL nodes and FFTs.
% No parent-pair enumeration and no pair-based reference dependency.
% Input is a declared first-order record, NOT total measured elevation.
% Spatial auxiliary-domain error is assessed by domain refinement.
arguments
    eta1 (:,1) double {mustBeFinite,mustBeReal}
    t (:,1) double {mustBeFinite,mustBeReal}
    g (1,1) double {mustBePositive}
    h (1,1) double {mustBePositive}
    options struct
end
allowed={'omega_max','peak_wavenumber','quadrature_rank','domain_lengths', ...
    'relative_tolerance','spatial_points','memory_budget_MiB','order'};
assert(all(ismember(fieldnames(options),allowed)),'Unknown option.');
assert(isfield(options,'omega_max') && options.omega_max>0, ...
    'Declare omega_max explicitly; no hidden parent-frequency cutoff.');
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
freezePath=fullfile(root,'symbolic','generated','finite_depth_directional_order2_eta22_pure_gl8.json');
freeze=jsondecode(fileread(freezePath));
assert(freeze.overall_exact_gate_pass && ~freeze.oracle_or_mf12_used_in_formula ...
    && ~freeze.stokes_correction_used && ~freeze.angular_correction_used);
if ~isfield(options,'quadrature_rank'),options.quadrature_rank=8;end
if ~isfield(options,'order'),options.order=2;end
assert(ismember(options.order,[2 3]));
if ~isfield(options,'domain_lengths'),options.domain_lengths=[25 50 100 200 400 800];end
if ~isfield(options,'relative_tolerance'),options.relative_tolerance=.0005;end
if ~isfield(options,'spatial_points'),options.spatial_points=[];end
if ~isfield(options,'memory_budget_MiB'),options.memory_budget_MiB=1024;end
J=options.quadrature_rank;assert(ismember(J,[4 6 8 12 16]));
if options.order==3,assert(ismember(J,[4 6 8]),'Nested original-GL execution currently supports shared inner/outer ranks 4, 6 or 8.');end
lengths=options.domain_lengths(:)';assert(all(lengths>0) && all(diff(lengths)>0));
assert(options.relative_tolerance>0 && options.memory_budget_MiB>0);
assert(isempty(options.spatial_points) || numel(options.spatial_points)==numel(lengths));
N=numel(t);assert(N==numel(eta1) && N>=9 && all(diff(t)>0));
dt=mean(diff(t));tr=t-t(1);period=N*dt;
timeDeviation=max(abs(tr-(0:N-1)'*dt));timeTolerance=max(1e-10,1024*eps(max(1,tr(end))));
assert(timeDeviation<=timeTolerance,'A uniform time record is required.');
positiveMax=floor((N-1)/2);allBins=(1:positiveMax)';allOmega=2*pi*allBins/period;
bins=allBins(allOmega<=options.omega_max);assert(~isempty(bins));omega=2*pi*bins/period;
F=fft(eta1)/N;A=2*conj(F(bins+1));assert(norm(A)>0,'The declared input band is empty.');
q=zeros(size(omega));nu=omega*sqrt(h/g);
for n=1:numel(q)
    q(n)=fzero(@(v)v*tanh(v)-nu(n)^2,[0,max(1,2*nu(n)^2+2)]);
end
assert(all(diff(q,2)>-1e-12),'Convex dispersion is required.');
if isfield(options,'peak_wavenumber')
    qp=h*options.peak_wavenumber;assert(qp>0);
else
    [~,peak]=max(abs(A));qp=q(peak);
end
lambda=sqrt(4*qp*tanh(qp)-2*qp*tanh(2*qp));assert(lambda>0);
[nodes,weights]=laguerre(J);
if options.order==3
    assert(3*max(bins)<=positiveMax,'Cubic execution currently requires all declared cubic sums inside the native output band. No parent cutoff was applied.');
    freeze3=jsondecode(fileread(fullfile(root,'symbolic','generated','finite_depth_directional_order3_nested_green_laplace.json')));
    assert(freeze3.overall_exact_gate_pass && ~freeze3.oracle_or_mf12_used && ~freeze3.input_chi_route_used);
end
workCount=min(N,2^nextpow2(options.order*max(bins)+1));assert(workCount>options.order*max(bins));
outbins=(2*min(bins):min(2*max(bins),positiveMax));assert(~isempty(outbins));
base=zeros(N,1);base(bins+1)=A;used=real(fft(base));
lower=max(min(bins),outbins-max(bins));upper=outbins-lower;
Qmax=reshape(q(lower-min(bins)+1)+q(upper-min(bins)+1),1,[]);
Qmin=reshape(q(floor(outbins/2)-min(bins)+1)+q(ceil(outbins/2)-min(bins)+1),1,[]);
s=2*pi*outbins/period*sqrt(h/g);
template=struct('L_over_h',0,'Nx',0,'seconds',0,'estimated_MiB',0,'relative_change',NaN,'eta22_change',NaN,'eta33_change',NaN);
levels=repmat(template,1,numel(lengths));previous=[];converged=false;timer=tic;
for level=1:numel(lengths)
    L=lengths(level);
    if isempty(options.spatial_points)
        Nx=2^nextpow2(max(16,floor(options.order*max(q)*L/pi+2)+1));
    else
        Nx=options.spatial_points(level);
    end
    assert(mod(Nx,2)==0 && Nx>=8);
    workFields=14;if options.order==3,workFields=24;end
    estimateMiB=16*Nx*(numel(A)+workFields*workCount+8*numel(outbins))/2^20;
    assert(estimateMiB<=options.memory_budget_MiB,'Requested auxiliary grid exceeds the declared memory budget. No input frequencies were removed.');
    mode=[0:Nx/2-1,-Nx/2:-1]';x=mode*L/Nx;K=mode*2*pi/L;
    assert(options.order*max(q)<max(K),'Spatial grid cannot represent the declared products.');
    ticLevel=tic;
    inputSpectrum=complex(zeros(Nx,workCount));inputSpectrum(:,bins+1)=exp(1i*x*q.').*A.';
    filters=zeros(6,workCount);
    filters(:,bins+1)=[ones(size(q))';nu';(nu.^2)';(q./nu)';(q.^2./nu)';q'];
    v=fft(inputSpectrum.*filters(1,:),[],2);
    vn=fft(inputSpectrum.*filters(2,:),[],2);vn2=fft(inputSpectrum.*filters(3,:),[],2);
    hx=fft(inputSpectrum.*filters(4,:),[],2);radial=fft(inputSpectrum.*filters(5,:),[],2);
    jx=fft(inputSpectrum.*filters(6,:),[],2);
    if options.order==3
        extra=zeros(3,workCount);extra(:,bins+1)=[(nu.*q)';(q.^2.*nu)';(q.^2)'];
        nujx=fft(inputSpectrum.*extra(1,:),[],2);
        q2nu=fft(inputSpectrum.*extra(2,:),[],2);q2field=fft(inputSpectrum.*extra(3,:),[],2);
    end
    clear inputSpectrum;
    sd=2*v.*vn2+vn.^2-hx.^2;sk=2*v.*radial+2*hx.*jx;
    if options.order==2,clear v vn vn2 hx radial jx;end
    ds=ifft(sd,[],2);ks=ifft(sk,[],2);clear sd sk;
    ds=fft(ds(:,outbins+1),[],1);ks=fft(ks(:,outbins+1),[],1);
    % Same prescribed GL nodes as the spatial graph. Parent damping is
    % transferred to the OUTPUT temporal frequency before spatial response.
    Keval=min(max(K,Qmin),Qmax);a=sqrt(Keval.*tanh(Keval));
    assert(all(s>a,'all'));
    md=zeros(size(a));mk=md;
    for node=1:J
        tau=nodes(node)/lambda;
        ep=weights(node)/lambda*exp(nodes(node)-(s-a)*tau);
        em=weights(node)/lambda*exp(nodes(node)-(s+a)*tau);
        md=md-a.*(ep-em)/2;mk=mk+(ep+em)/2;
    end
    output=complex(zeros(N,1));output(outbins+1)=sum(md.*ds+mk.*ks,1).'/(4*h*Nx);
    analytic=fft(output);analytic3=[];outbins3=[];
    if options.order==3
        % Original nested GL graph: its inner/outer scales differ from the
        % standalone eta22 determinant-root scale and are preserved here.
        lambda2=2*sqrt(qp*tanh(qp))-sqrt(2*qp*tanh(2*qp));
        lambda3=3*sqrt(qp*tanh(qp))-sqrt(3*qp*tanh(3*qp));
        [S2,C2]=node_response(s,a,lambda2,nodes,weights);
        E2=(-a.^2.*S2.*ds+C2.*ks)/4;
        P2=1i*(C2.*ds-S2.*ks)/4;clear ds ks S2 C2;
        lower=restore(E2,outbins,workCount);fk=1i*lower.*radial;fd=-lower.*vn2;
        lower=restore(1i*Keval.*E2,outbins,workCount);fk=fk+hx.*lower;clear E2;
        lower=restore(1i*Keval.*P2,outbins,workCount);fk=fk+1i*lower.*jx;fd=fd+hx.*lower;
        lower=restore(Keval.^2.*P2,outbins,workCount);fk=fk-v.*lower;
        lower=restore(-1i*s.*a.^2.*P2,outbins,workCount);fd=fd+v.*lower;
        lower=restore(a.^2.*P2,outbins,workCount);fd=fd-1i*vn.*lower;clear lower P2;
        fk=fk/2+(1i*v.*nujx.*jx+.5i*v.^2.*q2nu)/4;
        fd=fd/2+(-.5*v.^2.*q2field+v.*(hx.*nujx-vn.*radial))/4;
        clear v vn vn2 hx radial jx nujx q2nu q2field;
        outbins3=3*min(bins):3*max(bins);
        ds3=ifft(fd,[],2);ks3=ifft(fk,[],2);clear fd fk;
        ds3=fft(ds3(:,outbins3+1),[],1);ks3=fft(ks3(:,outbins3+1),[],1);
        [lo3,hi3]=cubic_support(outbins3,bins,q);
        Q3=min(max(K,lo3),hi3);a3=sqrt(Q3.*tanh(Q3));
        s3=2*pi*outbins3/period*sqrt(h/g);
        [S3,C3]=node_response(s3,a3,lambda3,nodes,weights);
        output3=complex(zeros(N,1));
        output3(outbins3+1)=sum(a3.^2.*S3.*ds3-1i*C3.*ks3,1).'/(Nx*h^2);
        analytic3=fft(output3);
    end
    seconds=toc(ticLevel);current=[analytic,analytic3];changes=NaN(1,size(current,2));
    if ~isempty(previous),changes=sqrt(sum(abs(current-previous).^2,1))./max(sqrt(sum(abs(current).^2,1)),realmin);end
    change=max(changes);change3=NaN;if options.order==3,change3=changes(2);end
    levels(level)=struct('L_over_h',L,'Nx',Nx,'seconds',seconds,'estimated_MiB',estimateMiB,'relative_change',change,'eta22_change',changes(1),'eta33_change',change3);
    previous=current;
    if isfinite(change) && change<=options.relative_tolerance,converged=true;break;end
end
levels=levels(1:level);
result=struct('t',t,'eta1_used',used,'eta22',real(analytic),'eta22_analytic',analytic);
if options.order==3,result.eta33=real(analytic3);result.eta33_analytic=analytic3;end
nyquistCoefficient=0;if mod(N,2)==0,nyquistCoefficient=real(F(N/2+1));end
audit=struct('implementation','prescribed-GL-quadrature-space-time-FFT', ...
    'order',options.order,'output_component','positive pure-sum elevation', ...
    'formula_freeze',freeze.candidate_id,'formula_freeze_path',freezePath, ...
    'quadrature_rank',J,'peak_wavenumber',qp/h,'lambda',lambda, ...
    'parent_bins',bins,'parent_omega',omega,'parent_k',q/h, ...
    'output_bins',outbins,'output_projection','native representable positive temporal sums', ...
    'omitted_output_bins',max(0,2*max(bins)-positiveMax), ...
    'input_projection_relative',norm(used-eta1)/norm(eta1),'minimum_parent_kh',min(q), ...
    'parents_below_released_kh_point3',nnz(q<.3), ...
    'excluded_DC_m',real(F(1)),'excluded_Nyquist_coefficient_m',nyquistCoefficient, ...
    'output_time_count',N,'work_time_count',workCount, ...
    'pair_enumeration',false,'triple_enumeration',false,'stokes_correction',false, ...
    'time_grid_deviation_s',timeDeviation,'levels',levels,'converged',converged, ...
    'relative_tolerance',options.relative_tolerance,'solve_seconds',toc(timer), ...
    'response_extension','constant outside the declared positive-sum spatial support', ...
    'convergence_note','Successive-domain change is a diagnostic, not an absolute error bound.');
if options.order==3
    audit.order3_formula_freeze=freeze3.candidate_id;audit.inner_scale=lambda2;audit.outer_scale=lambda3;
    audit.eta33_output_bins=outbins3;
    audit.parents_at_or_below_released_kh_point5=nnz(q<=.5);
end
end
function [S,C]=node_response(s,a,lambda,nodes,weights)
assert(lambda>0 && all(s>a,'all'));
S=zeros(size(a));C=S;
for j=1:numel(nodes)
    tau=nodes(j)/lambda;
    ep=weights(j)/lambda*exp(nodes(j)-(s-a)*tau);
    em=weights(j)/lambda*exp(nodes(j)-(s+a)*tau);
    S=S+(ep-em)./(2*a);C=C+(ep+em)/2;
end
end
function value=restore(coefficients,bins,N)
spectrum=complex(zeros(size(coefficients,1),N));spectrum(:,bins+1)=ifft(coefficients,[],1);
value=fft(spectrum,[],2);
end
function [low,high]=cubic_support(out,bins,q)
base=floor(out/3);extra=mod(out,3);
low=(3-extra).*reshape(q(base-min(bins)+1),1,[]) ...
    +extra.*reshape(q(min(base+1,max(bins))-min(bins)+1),1,[]);
remaining=out;high=zeros(size(out));
for leaf=1:3
    index=min(max(bins),remaining-(3-leaf)*min(bins));
    high=high+reshape(q(index-min(bins)+1),1,[]);remaining=remaining-index;
end
end
function [nodes,weights]=laguerre(J)
[vectors,values]=eig(diag(2*(1:J)-1)+diag(1:J-1,1)+diag(1:J-1,-1),'vector');
[nodes,index]=sort(values);weights=vectors(1,index).^2;
for degree=0:2*J-1
    assert(abs(sum(weights.*nodes'.^degree)-factorial(degree))/max(1,factorial(degree))<3e-11);
end
end
