function report=validate_order4_cpp_paper(output_path)
%VALIDATE_ORDER4_CPP_PAPER Compare a C++ output to the archived paper GL field.
root=setup_green_laplace();here=fileparts(mfilename('fullpath'));
if nargin<1,output_path=fullfile(root,'artifacts','order4-smoke','paper_gl6.bin');end
fid=fopen(output_path,'rb');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));
assert(strcmp(char(fread(fid,8,'*uint8').'),'GLE44O01'));
header=fread(fid,3,'uint64');fread(fid,1,'double');n=header(1);rank=header(2);
raw=fread(fid,[2 n*n],'double');assert(isempty(fread(fid,1,'uint8')));
field=real(reshape(complex(raw(1,:),raw(2,:)),n,n));
saved=load(fullfile(here,'data','eta44_no_stokes_rank_fields.mat'));
expected=saved.fields{find([6 8 10]==rank,1)};
assert(isequal(size(field),size(expected)));
relative_error=norm(field(:)-expected(:))/norm(expected(:));
assert(relative_error<2e-8,'C++ paper field differs from archived GL.');
report=struct('grid',[n n],'rank',rank,'raw_relative_l2',relative_error, ...
    'stokes_correction_used',false,'reference','archived same-rank GL paper field');
fid2=fopen(fullfile(root,'artifacts','order4-smoke','paper_cpp_validation.json'),'w');
fprintf(fid2,'%s\n',jsonencode(report,PrettyPrint=true));fclose(fid2);disp(report);
end
