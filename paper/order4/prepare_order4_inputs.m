function input=prepare_order4_inputs(case_directory,output_directory)
%PREPARE_ORDER4_INPUTS Deposit the frozen components without renormalization.
if ~isfolder(output_directory),mkdir(output_directory);end
manifest=jsondecode(fileread(fullfile(case_directory,'manifest.json')));
frozen=read_directional_eta4_cpp_input(fullfile(case_directory,'directional_eta44_input.bin'));
assert(frozen.caseCount==1);
c=frozen.eta11Analytic;
t=table(frozen.modeIndices(:,1),frozen.modeIndices(:,2),'VariableNames',{'mx','my'});
n=manifest.nx;assert(n==manifest.ny,'Square FFT grids required.');
assert(all(4*abs(t.mx)<n/2) && all(4*abs(t.my)<n/2),'Quartic aliasing.');
assert(all(t.mx>0),'Strict-forward pure-sum support required.');
s=complex(zeros(n));
indices=sub2ind([n,n],mod(double(t.my),n)+1,mod(double(t.mx),n)+1);
assert(numel(unique(indices))==height(t));
s(indices)=n*n*c/manifest.h;
input=struct('eta11_spectrum',s,'config',struct( ...
    'gravity',manifest.g,'peak_kh',manifest.peakDepth, ...
    'peak_wavenumber',manifest.kp,'delta_k_over_kp',manifest.dk/manifest.kp), ...
    'reference_or_higher_order_field_included',false);
save(fullfile(output_directory,'input.mat'),'input');
export_cpp_gl_runtime_input(fullfile(output_directory,'input.mat'), ...
    fullfile(output_directory,'gl_input.bin'));
write_directional_eta4_cpp_input(fullfile(output_directory,'wit_input.bin'), ...
    [t.mx,t.my],c,frozen.gravity,frozen.depth,frozen.deltaK,[n,n],frozen.threshold);
% The archive's WIT binary must be exactly reproduced, including normalization.
fid=fopen(fullfile(case_directory,'directional_eta44_input.bin'),'rb');
expected=fread(fid,Inf,'*uint8');fclose(fid);
fid=fopen(fullfile(output_directory,'wit_input.bin'),'rb');
actual=fread(fid,Inf,'*uint8');fclose(fid);
assert(isequal(actual,expected),'order4:FrozenInput','Frozen WIT input differs.');
end
