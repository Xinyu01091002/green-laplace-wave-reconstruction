function report=test_directional_difference_time()
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));addpath(root);setup_green_laplace;
addpath(fullfile(root,'research','unidirectional_time_series'));addpath(fullfile(root,'diagnostics','eta20'));addpath(fullfile(root,'diagnostics','eta20','generated'));
N=32;h=.7;g=9.81;kx=[2,2,3];ky=[1,-1,1];A=[.009+.002i,.006-.001i,.003+.003i];
w=sqrt(g*hypot(kx,ky).*tanh(h*hypot(kx,ky)));[KX,KY]=meshgrid([0:15,-16:-1]);err=[];
for time=[0,.27]
    sp=complex(zeros(N));a=A.*exp(-1i*w*time);
    for j=1:3,sp(mod(ky(j),N)+1,kx(j)+1)=N*N*a(j)/2;sp(mod(-ky(j),N)+1,mod(-kx(j),N)+1)=N*N*conj(a(j))/2;end
    for rank=[6,12,16]
        ref=2*eta20_green_laplace_shared(sp,KX,KY,h,"shared"+string(rank),string(root));
        for col=[1,5,9]
            shifted=A.*exp(1i*(kx*(col-1)*2*pi/N+ky*3*2*pi/N));
            [candidate,audit]=gl_directional_difference_time(shifted,w,kx,ky,h,time,rank);
            assert(audit.nonzero_spatial_stationary_pairs==2);
            err(end+1)=abs(candidate-ref(4,col)); %#ok<AGROW>
        end
    end
end
tt=(0:63)'*.2;bins=[18;18;23];ww=2*pi*bins/(64*.2);kk=zeros(3,1);
for j=1:3,kk(j)=fzero(@(k)g*k*tanh(k)-ww(j)^2,[0,30]);end
ang=[.3;-.3;.1];xx=kk.*cos(ang);yy=kk.*sin(ang);
[e,p]=gl_directional_sum_time(A,ww,xx,yy,g,1,kk(1),tt,8);
[ea,pa]=gl_directional_sum_time(A,ww,xx,yy,g,1,kk(1),tt,8,bins);
v=gl_directional_difference_time(A,ww,xx,yy,1,tt,16);
va=gl_directional_difference_time(A,ww,xx,yy,1,tt,16,bins);
aggregation=max([norm(e-ea)/norm(e),norm(p-pa)/norm(p),norm(v-va)/norm(v)]);
report=struct('spatial_max_absolute_error_m',max(err),'sample_aggregation_relative_error',aggregation);
disp(report);assert(max(err)<1e-10 && aggregation<1e-11);
out=fullfile(root,'artifacts','directional_sweep');if ~isfolder(out),mkdir(out);end
fid=fopen(fullfile(out,'eta20_checks.json'),'w');c=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));
end
