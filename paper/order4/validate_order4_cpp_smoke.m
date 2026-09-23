function report=validate_order4_cpp_smoke()
root=setup_green_laplace();here=fileparts(mfilename('fullpath'));addpath(here);
out=fullfile(root,'artifacts','order4-smoke');saved=load(fullfile(out,'input.mat'));
input=saved.input;n=size(input.eta11_spectrum,1);axis=[0:n/2-1,-n/2:-1];
[mx,my]=meshgrid(axis,axis);ranks=[6 8 10];errors=zeros(3,1);wit_errors=errors;
wit=read_directional_eta4_cpp_output(fullfile(out,'wit.bin'));
materialized=read_directional_eta4_cpp_output(fullfile(out,'wit_materialized.bin'));
assert(wit.allFinite && wit.nearResonantCount==0 && wit.maximumPairedResidual<1e-10);
wit_parity=norm(wit.etaFixedB(:)-materialized.etaFixedB(:))/norm(materialized.etaFixedB(:));
assert(wit_parity<1e-12);
ref=real(ifft2(wit.etaFixedB(:,:,1).'*n*n));
for k=1:numel(ranks)
    fid=fopen(fullfile(out,sprintf('gl%d.bin',ranks(k))),'rb');
    assert(strcmp(char(fread(fid,8,'*uint8').'),'GLE44O01'));
    header=fread(fid,3,'uint64');elapsed=fread(fid,1,'double'); %#ok<NASGU>
    raw=fread(fid,[2 n*n],'double');assert(isempty(fread(fid,1,'uint8')));fclose(fid);
    assert(header(1)==n && header(2)==ranks(k));
    actual=reshape(complex(raw(1,:),raw(2,:)),n,n);
    [expected,a]=finite_depth_directional_green_laplace_eta44_three_kernels( ...
        input.eta11_spectrum,.5*mx,.5*my,1,root,quadrature_rank=ranks(k));
    assert(~a.all_nested_stokes_corrections_used);
    errors(k)=norm(actual(:)-expected(:))/norm(expected(:));
    assert(errors(k)<2e-8,'C++ / MATLAB GL mismatch.');
    candidate=real(actual);wit_errors(k)=norm(candidate(:)-ref(:))/norm(ref(:));
end
fid=fopen(fullfile(out,'gl6_split.bin'),'rb');fread(fid,8,'*uint8');fread(fid,3,'uint64');fread(fid,1,'double');
raw=fread(fid,[2 n*n],'double');fclose(fid);split=reshape(complex(raw(1,:),raw(2,:)),n,n);
[full,~]=finite_depth_directional_green_laplace_eta44_three_kernels(input.eta11_spectrum,.5*mx,.5*my,1,root,quadrature_rank=6);
split_error=norm(split(:)-full(:))/norm(full(:));assert(split_error<2e-8);
report=struct('matlab_release',version('-release'),'precision','double', ...
    'grid',[n n],'ranks',ranks,'cpp_matlab_gl_relative_l2',errors, ...
    'wit_materialized_streaming_relative_l2',wit_parity, ...
    'gl6_split_matlab_relative_l2',split_error,'gl_wit_raw_relative_l2',wit_errors, ...
    'stokes_correction_used',false, ...
    'scope','migration and implementation parity; GL-WIT errors are non-selecting diagnostics');
fid=fopen(fullfile(out,'validation.json'),'w');fprintf(fid,'%s\n',jsonencode(report,PrettyPrint=true));fclose(fid);
disp(report);
end
