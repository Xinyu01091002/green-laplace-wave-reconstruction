function test_directional_joint_model()
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));addpath(root);
setup_green_laplace;addpath(fullfile(root,'research','unidirectional_time_series'));
A=[.009+.002i,.006-.001i,.003+.003i];kx=[2,2,3];ky=[1,-1,1];h=.7;g=9.81;kp=hypot(kx(1),ky(1));
w=sqrt(g*hypot(kx,ky).*tanh(h*hypot(kx,ky)));
c=gl_spectral_coefficients(2,g,h,real(A),imag(A),kx,ky,0,0,struct('eta22_rank',8,'peak_wavenumber',kp));
errors=[];
for t=[0,.21]
    [~,~,X,Y,parts]=gl_spectral_surface(c,2*pi,2*pi,32,32,t);
    for pos=[1,7,19]
        shifted=A.*exp(1i*(kx*X(3,pos)+ky*Y(3,pos)));
        [e,p]=gl_directional_sum_time(shifted,w,kx,ky,g,h,kp,t,8);
        errors(end+1,:)=[abs(e-parts.eta22_plus(3,pos)),abs(p-parts.psi22_plus(3,pos))]; %#ok<AGROW>
    end
end
prior=[1+.2i,.4-.1i,.2+.3i;.7-.1i,.2+.4i,.5-.2i];
observed=sum(prior,2).*[1.3*exp(.4i);.8*exp(-.3i)];
[coeff,condition]=allocate_directional_record(prior,observed);
assert(norm(sum(coeff,2)-observed)<1e-13 && max(errors,[],'all')<1e-11);
[unchanged,~]=allocate_directional_record(prior,sum(prior,2));assert(norm(unchanged-prior)<1e-13);
blocked=false;
try
    allocate_directional_record([1,-1],1);
catch
    blocked=true;
end
assert(blocked);
report=struct('maximum_absolute_eta22_difference',max(errors(:,1)), ...
    'maximum_absolute_psi22_difference',max(errors(:,2)), ...
    'allocation_sum_error',condition.sum_relative_error,'exact_cancellation_rejected',blocked);
out=fullfile(root,'artifacts','directional_joint_input');if ~isfolder(out),mkdir(out);end
fid=fopen(fullfile(out,'implementation_checks.json'),'w');c=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));disp(report);
end
