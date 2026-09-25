function report=test_surface_potential_time(mf12Root)
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));addpath(root);
setup_green_laplace('MF12Root',mf12Root);
g=9.81;h=.7;kp=2;k=[1,2,3];A=[.009+.002i,.006-.001i,.003+.003i];
w=sqrt(g*k.*tanh(h*k));errors=[];mferrors=[];
c=gl_spectral_coefficients(2,g,h,real(A),imag(A),k,[0,0,0],0,0,struct('peak_wavenumber',kp));
for t=[0,.17,.43]
    [~,~,X,~,p]=gl_spectral_surface(c,2*pi,2*pi,32,4,t);
    comparison=gl_mf12_order2_comparison(c,2*pi,2*pi,32,4,t);
    for col=[1,5,11]
        shifted=A.*exp(1i*k*X(1,col));
        pred=gl_psi22_time_pairs(shifted,w,k,g,h,kp,t);
        errors(end+1)=abs(pred-p.psi22_plus(1,col)); %#ok<AGROW>
        mf=mf12_psi22_time(shifted,k,g,h,t);
        mferrors(end+1)=abs(mf-comparison.mf12_psi22(1,col)); %#ok<AGROW>
    end
end
report=struct('gl_max_absolute_error',max(errors),'mf12_max_absolute_error',max(mferrors));
disp(report);assert(max(errors)<1e-11 && max(mferrors)<1e-11);
fid=fopen(fullfile(root,'artifacts','unidirectional_time_series','psi22_checks.json'),'w');
cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));
end
