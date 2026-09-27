function make_time_reference(outdir)
% Direct coefficient convolution + RK4, independent of GPU/FFT implementation.
assert(~isfolder(outdir),'Use a fresh directory');mkdir(outdir);
M=5;bx=3;by=2;nx=64;ny=32;lx=9;ly=7;h=1.3;g=9.81;dt=.005;steps=40;
sx=2*M*bx+1;sy=2*M*by+1;cx=M*bx+1;cy=M*by+1;
[mx,my]=meshgrid(-M*bx:M*bx,-M*by:M*by);kx=2*pi*mx/lx;ky=2*pi*my/ly;k=hypot(kx,ky);
mask=abs(mx)<=bx & abs(my)<=by;
eta=complex(zeros(sy,sx));psi=eta;eta(cy,cx)=.007;psi(cy,cx)=.03;
parent=[1 2 .031 .012 -.017 .027;3 -1 -.014 .009 .024 -.011];
for row=parent.'
    x=row(1);y=row(2);e=complex(row(3),row(4));p=complex(row(5),row(6));
    eta(cy+y,cx+x)=e;eta(cy-y,cx-x)=conj(e);psi(cy+y,cx+x)=p;psi(cy-y,cx-x)=conj(p);
end
D=cell(M+1,1);for j=1:M+1,D{j}=k.^j.*(mod(j,2)*tanh(k*h)+(1-mod(j,2)));end
one=complex(zeros(sy,sx));one(cy,cx)=1;
[ef,pf]=integrate(eta,psi,dt,steps);[er,pr]=integrate(eta,psi,dt/2,2*steps);
refinement=struct('eta_coefficient_max_abs',max(abs(ef-er),[],'all'),...
    'psi_coefficient_max_abs',max(abs(pf-pr),[],'all'),'dt',dt,'steps',steps,'duration',dt*steps);
assert(refinement.eta_coefficient_max_abs<1e-8 && refinement.psi_coefficient_max_abs<1e-8);
assert(abs(ef(cy,cx)-eta(cy,cx))<1e-12,'Mean eta changed');
file=fullfile(outdir,'hos_time.bin');fid=fopen(file,'w','ieee-le');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));
fwrite(fid,[nx ny M bx by sx sy steps],'int32');fwrite(fid,[lx ly h g dt],'double');
for c={eta,psi,ef,pf}
    value=c{1};assert(max(abs(value-conj(rot90(value,2))),[],'all')<1e-12);
    fwrite(fid,real(value).','double');fwrite(fid,imag(value).','double');
end
clear cleanup
fid=fopen(fullfile(outdir,'time_refinement.json'),'w');fprintf(fid,'%s\n',jsonencode(refinement));fclose(fid);
disp(refinement);
    function [et,pt]=rhs(e,p)
        mul=@(a,b)conv2(a,b,'same');ep=cell(M,1);ep{1}=one;
        for j=1:M-1,ep{j+1}=mul(ep{j},e)/j;end
        phi=cell(M,1);W=cell(M,1);phi{1}=p;
        for n=1:M
            if n>1
                phi{n}=complex(zeros(sy,sx));for j=1:n-1,phi{n}=phi{n}-mul(ep{j+1},D{j}.*phi{n-j});end
            end
            W{n}=complex(zeros(sy,sx));for j=0:n-1,W{n}=W{n}+mul(ep{j+1},D{j+1}.*phi{n-j});end
        end
        ex=1i*kx.*e;ey=1i*ky.*e;px=1i*kx.*p;py=1i*ky.*p;Q=mul(ex,ex)+mul(ey,ey);
        et=-mul(ex,px)-mul(ey,py);pt=-g*e-.5*(mul(px,px)+mul(py,py));
        for n=1:M
            et=et+W{n};if n>=3,et=et+mul(Q,W{n-2});end
            for j=1:n-1,pt=pt+.5*mul(W{j},W{n-j});end
            for j=1:n-3,pt=pt+.5*mul(Q,mul(W{j},W{n-2-j}));end
        end
        et=et.*mask;pt=pt.*mask;
    end
    function [e,p]=integrate(e,p,step,count)
        for it=1:count
            [a,b]=rhs(e,p);[c,d]=rhs(e+.5*step*a,p+.5*step*b);
            [v,w]=rhs(e+.5*step*c,p+.5*step*d);[u,z]=rhs(e+step*v,p+step*w);
            e=e+step/6*(a+2*c+2*v+u);p=p+step/6*(b+2*d+2*w+z);
        end
    end
end
