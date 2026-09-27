function make_hos_reference(outdir)
% Independent Fourier-convolution HOS reference. No FFT, no CUDA code reuse.
% Coefficients multiply exp(i k.x). Constant bed and true surface potential.
assert(~isfolder(outdir),'Use a fresh reference directory');mkdir(outdir);
M=5;bx=3;by=2;nx=64;ny=32;lx=9;ly=7;g=9.81;
sx=2*M*bx+1;sy=2*M*by+1;cx=M*bx+1;cy=M*by+1;
[mx,my]=meshgrid(-M*bx:M*bx,-M*by:M*by);
kx=2*pi*mx/lx;ky=2*pi*my/ly;k=hypot(kx,ky);
eta=complex(zeros(sy,sx));psi=eta;one=eta;one(cy,cx)=1;
eta(cy,cx)=.007;psi(cy,cx)=.03;
parent=[1 2 .031 .012 -.017 .027;3 -1 -.014 .009 .024 -.011];
for row=parent.'
    x=row(1);y=row(2);e=complex(row(3),row(4));p=complex(row(5),row(6));
    eta(cy+y,cx+x)=e;eta(cy-y,cx-x)=conj(e);
    psi(cy+y,cx+x)=p;psi(cy-y,cx-x)=conj(p);
end
for h=[.15 1.3 20]
    D=@(j) k.^j.*(mod(j,2)*tanh(k*h)+(1-mod(j,2)));
    multiply=@(a,b) conv2(a,b,'same');
    powers=cell(M,1);powers{1}=one;
    for j=1:M-1,powers{j+1}=multiply(powers{j},eta)/j;end
    phi=cell(M,1);W=cell(M,1);phi{1}=psi;
    for n=1:M
        if n>1
            phi{n}=complex(zeros(sy,sx));
            for j=1:n-1,phi{n}=phi{n}-multiply(powers{j+1},D(j).*phi{n-j});end
        end
        W{n}=complex(zeros(sy,sx));
        for j=0:n-1,W{n}=W{n}+multiply(powers{j+1},D(j+1).*phi{n-j});end
    end
    ex=1i*kx.*eta;ey=1i*ky.*eta;px=1i*kx.*psi;py=1i*ky.*psi;
    grad2=multiply(ex,ex)+multiply(ey,ey);
    E=cell(M,1);P=cell(M,1);
    for n=1:M
        E{n}=W{n};P{n}=complex(zeros(sy,sx));
        if n==1,P{n}=-g*eta;end
        if n==2
            E{n}=E{n}-multiply(ex,px)-multiply(ey,py);
            P{n}=P{n}-.5*(multiply(px,px)+multiply(py,py));
        end
        if n>=3,E{n}=E{n}+multiply(grad2,W{n-2});end
        for j=1:n-1,P{n}=P{n}+.5*multiply(W{j},W{n-j});end
        for j=1:n-3,P{n}=P{n}+.5*multiply(grad2,multiply(W{j},W{n-2-j}));end
    end
    % Independent order-two identity and linear limit, before exporting.
    G=k.*tanh(k*h);et2=-G.*multiply(eta,G.*psi)-1i*kx.*multiply(eta,px)-1i*ky.*multiply(eta,py);
    assert(max(abs(E{2}-et2),[],'all')<1e-13);
    assert(max(abs(E{1}-G.*psi),[],'all')<1e-13);
    fields={eta,psi};for n=1:M,fields=[fields,{W{n},E{n},P{n}}];end %#ok<AGROW>
    for j=1:numel(fields)
        assert(max(abs(fields{j}-conj(rot90(fields{j},2))),[],'all')<1e-12);
    end
    path=fullfile(outdir,sprintf('hos_h%g.bin',h));
    fid=fopen(path,'w','ieee-le');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));
    fwrite(fid,[nx ny M bx by sx sy],'int32');fwrite(fid,[lx ly h g],'double');
    for j=1:numel(fields)
        fwrite(fid,real(fields{j}).','double');fwrite(fid,imag(fields{j}).','double');
    end
    clear cleanup
    fprintf('Reference written: %s\n',path);
end
end
