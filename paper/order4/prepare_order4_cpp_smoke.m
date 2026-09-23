function prepare_order4_cpp_smoke()
root=setup_green_laplace();here=fileparts(mfilename('fullpath'));addpath(here);
out=fullfile(root,'artifacts','order4-smoke');if ~isfolder(out),mkdir(out);end
n=32;modes=[2,0;3,1;3,-1];c=[.012;.007*exp(.31i);.005*exp(-.47i)];
s=complex(zeros(n));s(sub2ind([n,n],mod(modes(:,2),n)+1,modes(:,1)+1))=n^2*c;
input=struct('eta11_spectrum',s,'config',struct('gravity',1,'peak_kh',1, ...
    'peak_wavenumber',1,'delta_k_over_kp',0.5));
save(fullfile(out,'input.mat'),'input');
export_cpp_gl_runtime_input(fullfile(out,'input.mat'),fullfile(out,'gl_input.bin'));
write_directional_eta4_cpp_input(fullfile(out,'wit_input.bin'),modes,c,1,1,[.5 .5],[n n]);
end
