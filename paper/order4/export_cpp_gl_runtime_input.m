function report=export_cpp_gl_runtime_input(input_path,output_path)
%EXPORT_CPP_GL_RUNTIME_INPUT Export frozen spectrum and GL rules to C++.
arguments
    input_path (1,1) string
    output_path (1,1) string
end
loaded=load(input_path,'input');input=loaded.input;
output_root=fileparts(output_path);
if ~isfolder(output_root),mkdir(output_root);end
fid=fopen(output_path,'wb');if fid<0,error('Cannot open C++ GL input.');end
cleanup=onCleanup(@()fclose(fid));
fwrite(fid,uint8('GLRTI001'),'uint8');
n=size(input.eta11_spectrum,1);
depth=input.config.peak_kh/input.config.peak_wavenumber;
dk=input.config.delta_k_over_kp*input.config.peak_wavenumber;
fwrite(fid,uint64(n),'uint64');
fwrite(fid,[dk,depth,input.config.peak_kh],'double');
ranks=[6,8,10];
fwrite(fid,uint64(numel(ranks)),'uint64');
for rank=ranks
    [nodes,weights]=gauss_laguerre_rule(rank);
    fwrite(fid,uint64(rank),'uint64');
    fwrite(fid,nodes,'double');
    fwrite(fid,weights,'double');
end
interleaved=zeros(2,numel(input.eta11_spectrum));
interleaved(1,:)=real(input.eta11_spectrum(:));
interleaved(2,:)=imag(input.eta11_spectrum(:));
fwrite(fid,interleaved,'double');
report=struct('status','PASS','grid_size',n,'dk',dk,'depth',depth, ...
    'peak_depth',input.config.peak_kh,'ranks',ranks, ...
    'input_path',input_path,'output_path',output_path);
fprintf('CPP_GL_RUNTIME_INPUT=PASS N=%d %s\n',n,output_path);
end

function [nodes,weights]=gauss_laguerre_rule(rank)
indices=1:rank;
jacobi=diag(2*indices-1)+diag(1:rank-1,1)+diag(1:rank-1,-1);
[vectors,values]=eig(jacobi,'vector');
[nodes,order]=sort(real(values));vectors=vectors(:,order);
weights=real(vectors(1,:)).^2;
nodes=nodes(:).';weights=weights(:).';
end

