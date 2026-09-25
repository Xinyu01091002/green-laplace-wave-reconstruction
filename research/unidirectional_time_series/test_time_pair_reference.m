function report = test_time_pair_reference(mf12root)
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root); setup_green_laplace('MF12Root',mf12root);
g=9.81; h=1; k=[1,2,3]; A=[.009+.002i,.006-.001i,.003+.003i];
w=sqrt(g*k.*tanh(k*h)); times=[0,.13,.47]; L=2*pi;
errors=[]; absolute_errors=[]; mferrors=[];
for rank=[4,6,8,12]
    c=gl_spectral_coefficients(2,g,h,real(A),imag(A),k,[0,0,0],0,0, ...
        struct('eta22_rank',rank,'peak_wavenumber',1));
    for t=times
        [~,~,X,~,p]=gl_spectral_surface(c,L,L,32,4,t);
        % Several points, nonzero complex phases and times exercise conventions.
        for col=[1,4,9]
            Ap=A.*exp(1i*k*X(1,col));
            y=gl_eta22_time_pairs(Ap,w,k,h,1,t,rank);
            errors(end+1)=abs(y-p.eta22_plus(1,col))/max(abs(p.eta22_plus(1,col)),1e-12); %#ok<AGROW>
            absolute_errors(end+1)=abs(y-p.eta22_plus(1,col)); %#ok<AGROW>
        end
        comparison=gl_mf12_order2_comparison(c,L,L,32,4,t);
        y=mf12_eta22_time(A,k,g,h,t);
        mferrors(end+1)=abs(y-comparison.mf12_eta22(1,1))/max(abs(y),1e-12); %#ok<AGROW>
    end
end
% Uniform temporal bins: recover the original complex amplitudes exactly.
N=128; dt=.2; bins=[3,5,8]; omega=2*pi*bins/(N*dt);
tt=(0:N-1)'*dt; signal=real(exp(-1i*tt*omega)*A.');
F=fft(signal)/N; recovered=2*conj(F(bins+1)).';
fft_error=norm(recovered-A)/norm(A);
fprintf('GL spatial max %.16g; MF12 spatial max %.16g; FFT %.16g\n',max(errors),max(mferrors),fft_error);
% The existing balanced FFT graph loses small absolute precision for this
% broad support. Also check the independently frozen ordered GL8 evaluator.
addpath(fullfile(root,'tests'));
[qx,qy]=meshgrid([0:15,-16:-1],[0:1,-2:-1]);
spec=complex(zeros(4,32)); spec(1,k+1)=128*A;
ordered=finite_depth_directional_pure_gl8_ordered_pair(spec,qx,qy,1,string(root));
direct=gl_eta22_time_pairs(A,w,k,h,1,0,8);
frozen_error=abs(direct-ordered(1,1))/abs(ordered(1,1));
assert(frozen_error<1e-12 && max(absolute_errors)<1e-12 && max(errors)<1e-8);
assert(max(mferrors)<1e-10 && fft_error<1e-13);
report=struct('pass',true,'max_gl_spatial_relative_error',max(errors), ...
    'max_gl_spatial_absolute_error_m',max(absolute_errors),'frozen_gl8_relative_error',frozen_error, ...
    'max_mf12_spatial_relative_error',max(mferrors),'fft_amplitude_relative_error',fft_error);
out=fullfile(root,'artifacts','unidirectional_time_series');
if ~isfolder(out),mkdir(out);end
fid=fopen(fullfile(out,'implementation_checks.json'),'w'); cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(report)); disp(report);
end
