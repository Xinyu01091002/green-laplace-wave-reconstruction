function report=test_eta20_eta33_time()
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root);setup_green_laplace;
addpath(fullfile(root,'diagnostics','eta20'));addpath(fullfile(root,'diagnostics','eta20','generated'));
n=32;ny=16;h=.5;g=9.81;k=[2,3,4];A=[.009+.002i,.006-.001i,.003+.003i];
omega=sqrt(g*k.*tanh(h*k));kp=3;
[kx,ky]=meshgrid([0:n/2-1,-n/2:-1],[0:ny/2-1,-ny/2:-1]);
err20=[];err33=[];
for time=[0,.27]
    amplitude=A.*exp(-1i*omega*time);
    sp=complex(zeros(ny,n)); sp(1,k+1)=n*ny*amplitude/h;
    for rank=[4,6,8]
        field=gl_no_stokes_eta33(sp,h*kx,h*ky,h*kp,rank,string(root));
        for col=[1,4,9]
            shifted=A.*exp(1i*k*(col-1)*2*pi/n);
            y=gl_eta33_time_triples(shifted,omega,k,h,kp,time,rank);
            err33(end+1)=abs(y-h*field(1,col)); %#ok<AGROW>
        end
    end
    sp=complex(zeros(ny,n));sp(1,k+1)=n*ny*amplitude/2;
    sp(1,mod(-k,n)+1)=n*ny*conj(amplitude)/2;
    for rank=[6,12,16]
        field=eta20_green_laplace_shared(sp,kx,ky,h,"shared"+string(rank),string(root));
        for col=[1,4,9]
            shifted=A.*exp(1i*k*(col-1)*2*pi/n);
            y=gl_eta20_time_pairs(shifted,omega,k,h,time,rank);
            err20(end+1)=abs(y-2*field(1,col)); %#ok<AGROW>
        end
    end
end
tt=(0:63)'*.2;bins=[5;7;9];ww=2*pi*bins/(64*.2);kk=zeros(3,1);
for j=1:3,kk(j)=fzero(@(x)g*x*tanh(x)-ww(j)^2,[0,10]);end
y1=gl_eta33_time_triples(A,ww,kk,1,kk(2),tt,6);
y2=gl_eta33_time_triples(A,ww,kk,1,kk(2),tt,6,bins);
aggregation_error=norm(y1-y2)/norm(y1);assert(aggregation_error<1e-12);
report=struct('eta20_max_absolute_difference_m',max(err20),'eta33_max_absolute_difference_m',max(err33), ...
    'eta33_uniform_bin_aggregation_relative_error',aggregation_error);
disp(report);assert(max(err20)<1e-11 && max(err33)<1e-11);
out=fullfile(root,'artifacts','unidirectional_time_series','eta20_eta33_checks.json');
fid=fopen(out,'w');cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));
end
