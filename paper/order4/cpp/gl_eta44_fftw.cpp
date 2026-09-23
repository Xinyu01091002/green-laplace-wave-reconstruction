#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>
#include <fftw3.h>
#include <omp.h>
using C=std::complex<double>;
using V=std::vector<C>;
using R=std::vector<double>;
using M=std::vector<uint8_t>;

template<class T>void rd(std::ifstream&f,T&v){
f.read(reinterpret_cast<char*>(&v),sizeof(T));
if(!f)throw std::runtime_error("input eof");

}

template<class T>void wr(std::ofstream&f,const T&v){
f.write(reinterpret_cast<const char*>(&v),sizeof(T));
if(!f)throw std::runtime_error("output failed");

}

class FFT{
int n;
size_t z;
fftw_complex*a,*b;
fftw_plan pf,pb;
public:FFT(int n_,int nt):n(n_),z(size_t(n_)*n_){
a=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*z);
b=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*z);
if(!a||!b||!fftw_init_threads())throw std::runtime_error("fftw init");
fftw_plan_with_nthreads(nt);
pf=fftw_plan_dft_2d(n,n,a,b,FFTW_FORWARD,FFTW_ESTIMATE);
pb=fftw_plan_dft_2d(n,n,a,b,FFTW_BACKWARD,FFTW_ESTIMATE);
if(!pf||!pb)throw std::runtime_error("fftw plan");

}
~FFT(){
fftw_destroy_plan(pf);
fftw_destroy_plan(pb);
fftw_cleanup_threads();
fftw_free(a);
fftw_free(b);

}
V go(const V&x,bool inv){
for(size_t i=0;
i<z;
++i){
a[i][0]=x[i].real();
a[i][1]=x[i].imag();

}
fftw_execute(inv?pb:pf);
V y(z);
double s=inv?1.0/double(z):1.;
for(size_t i=0;
i<z;
++i)y[i]=s*C(b[i][0],b[i][1]);
return y;

}

}
;

struct Grid{
int n,nt,rank;
size_t z;
double dk,depth,peak,l2,l3,l4,shift2,shift3,shift4;
R qx,qy,q,nu,A,sqrtA,safe;
M support,pairMask,tripleMask,quarticMask;
std::vector<int>axis;
FFT fft;
std::vector<double>nodes,weights;
Grid(int n_,int nt_,double dk_,double d,double p,V const&eta,std::vector<double>no,std::vector<double>we):n(n_),nt(nt_),rank(int(no.size())),z(size_t(n_)*n_),dk(dk_),depth(d),peak(p),qx(z),qy(z),q(z),nu(z),A(z),sqrtA(z),safe(z),support(z),axis(n),fft(n,nt),nodes(std::move(no)),weights(std::move(we)){
for(int i=0;
i<n;
++i)axis[i]=i<n/2?i:i-n;
double sc=0;
for(auto&v:eta)sc=std::max(sc,std::abs(v));
double th=1e-13*std::max(1.,sc);
for(size_t k=0;
k<z;
++k){
int ix=int(k/n),iy=int(k%n);
qx[k]=depth*dk*axis[ix];
qy[k]=depth*dk*axis[iy];
q[k]=std::hypot(qx[k],qy[k]);
nu[k]=std::sqrt(q[k]*std::tanh(q[k]));
A[k]=q[k]*std::tanh(q[k]);
sqrtA[k]=std::sqrt(A[k]);
safe[k]=std::max(nu[k],std::numeric_limits<double>::min());
if(std::abs(eta[k])>th){
if(qx[k]<=0||q[k]<=0)throw std::runtime_error("support");
support[k]=1;

}

}
pairMask=powerMask(2);
tripleMask=powerMask(3);
quarticMask=powerMask(4);
double mn=1e300,m2=0,m3=0,m4=0;
for(size_t k=0;
k<z;
++k){
if(support[k])mn=std::min(mn,nu[k]);
if(pairMask[k])m2=std::max(m2,sqrtA[k]);
if(tripleMask[k])m3=std::max(m3,sqrtA[k]);
if(quarticMask[k])m4=std::max(m4,sqrtA[k]);

}
shift2=.5*(mn+m2/2);
shift3=.5*(mn+m3/3);
shift4=.5*(mn+m4/4);
l2=2*std::sqrt(peak*std::tanh(peak))-std::sqrt(2*peak*std::tanh(2*peak));
l3=3*std::sqrt(peak*std::tanh(peak))-std::sqrt(3*peak*std::tanh(3*peak));
l4=4*std::sqrt(peak*std::tanh(peak))-std::sqrt(4*peak*std::tanh(4*peak));

}
M powerMask(int p){
V s(z);
for(size_t k=0;
k<z;
++k)s[k]=double(support[k]);
auto f=fft.go(s,false);
for(auto&v:f){
C b=v,o=1.;
for(int i=0;
i<p;
++i)o*=b;
v=o;

}
auto c=fft.go(f,true);
M m(z);
for(size_t k=0;
k<z;
++k)m[k]=c[k].real()>.5;
return m;

}

}
;

void zeroSuppress(V&x,M const&m){
double s=0;
for(size_t k=0;
k<x.size();
++k)if(m[k])s=std::max(s,std::abs(x[k]));
double f=16*std::numeric_limits<double>::epsilon()*std::max(1.,std::log2(double(x.size())));
for(size_t k=0;
k<x.size();
++k){
if(!m[k])x[k]=0.;
else if(s&&std::abs(x[k])<f*s)x[k]=0.;

}

}

template<class F>V invMul(Grid&g,V const&s,F fn){
V a(g.z);
for(size_t k=0;
k<g.z;
++k)a[k]=s[k]*fn(k);
return g.fft.go(a,true);

}

struct First{
V eta,ex,ey,px,py,pz,pxz,pyz,pzz,pzzz,ptz,ptzz,pxz_t,pyz_t,pzzz_t,ptz_t,pxzz,pyzz,pzzzz;

}
;

First first(Grid&g,V const&s,bool order4=false){
First f;
f.eta=g.fft.go(s,true);
f.ex=invMul(g,s,[&](size_t k){
return C(0,g.qx[k]);

}
);
f.ey=invMul(g,s,[&](size_t k){
return C(0,g.qy[k]);

}
);
f.px=invMul(g,s,[&](size_t k){
return C(g.qx[k]/g.safe[k],0);

}
);
f.py=invMul(g,s,[&](size_t k){
return C(g.qy[k]/g.safe[k],0);

}
);
f.pz=invMul(g,s,[&](size_t k){
return C(0,-g.nu[k]);

}
);
f.pxz=invMul(g,s,[&](size_t k){
return C(g.nu[k]*g.qx[k],0);

}
);
f.pyz=invMul(g,s,[&](size_t k){
return C(g.nu[k]*g.qy[k],0);

}
);
f.pzz=invMul(g,s,[&](size_t k){
return C(0,-g.q[k]*g.q[k]/g.safe[k]);

}
);
f.pzzz=invMul(g,s,[&](size_t k){
return C(0,-g.q[k]*g.q[k]*g.nu[k]);

}
);
f.ptz=invMul(g,s,[&](size_t k){
return C(-g.A[k],0);

}
);
f.ptzz=invMul(g,s,[&](size_t k){
return C(-g.q[k]*g.q[k],0);

}
);
f.pxz_t=invMul(g,s,[&](size_t k){
return C(0,-g.nu[k]*g.nu[k]*g.qx[k]);

}
);
f.pyz_t=invMul(g,s,[&](size_t k){
return C(0,-g.nu[k]*g.nu[k]*g.qy[k]);

}
);
f.pzzz_t=invMul(g,s,[&](size_t k){
return C(-g.q[k]*g.q[k]*g.A[k],0);

}
);
f.ptz_t=invMul(g,s,[&](size_t k){
return C(0,g.A[k]*g.nu[k]);

}
);
if(order4){
V phi(g.z);
for(size_t k=0;
k<g.z;
++k)if(g.nu[k]>0)phi[k]=C(0,-1)*s[k]/g.nu[k];
f.pxzz=invMul(g,phi,[&](size_t k){
return C(0,g.qx[k]*g.q[k]*g.q[k]);

}
);
f.pyzz=invMul(g,phi,[&](size_t k){
return C(0,g.qy[k]*g.q[k]*g.q[k]);

}
);
f.pzzzz=invMul(g,phi,[&](size_t k){
return C(g.q[k]*g.q[k]*g.q[k]*g.q[k],0);

}
);

}
return f;

}

struct Pair{
V eta,ex,ey,eta_t,ex_t,ey_t,px,py,pz,pzz,ptz,px_t,py_t,pzz_t,ptz_t,pxz,pyz,pzzz;
V phiSpec;

}
;

Pair pair(Grid&g,V const&s,bool order4=false){
V es(g.z),ets(g.z),ps(g.z),pts(g.z),ptts(g.z);
for(size_t ni=0;
ni<g.nodes.size();
++ni){
double x=g.nodes[ni],t=x/g.l2;
V damp(g.z),tmp(g.z);
for(size_t k=0;
k<g.z;
++k)damp[k]=s[k]*std::exp(-t*(g.nu[k]-g.shift2));
auto fil=[&](auto fn){
for(size_t k=0;
k<g.z;
++k)tmp[k]=damp[k]*fn(k);
return g.fft.go(tmp,true);

}
;
auto v=fil([&](size_t){
return C(1,0);

}
);
auto vn=fil([&](size_t k){
return C(g.nu[k],0);

}
);
auto vn2=fil([&](size_t k){
return C(g.nu[k]*g.nu[k],0);

}
);
auto vn3=fil([&](size_t k){
return C(g.nu[k]*g.nu[k]*g.nu[k],0);

}
);
auto vn4=fil([&](size_t k){
return C(g.nu[k]*g.nu[k]*g.nu[k]*g.nu[k],0);

}
);
auto hx=fil([&](size_t k){
return C(g.qx[k]/g.safe[k],0);

}
);
auto hy=fil([&](size_t k){
return C(g.qy[k]/g.safe[k],0);

}
);
auto jx=fil([&](size_t k){
return C(g.qx[k],0);

}
);
auto jy=fil([&](size_t k){
return C(g.qy[k],0);

}
);
auto ro=fil([&](size_t k){
return C(g.q[k]*g.q[k]/g.safe[k],0);

}
);
auto q2=fil([&](size_t k){
return C(g.q[k]*g.q[k],0);

}
);
auto q2n=fil([&](size_t k){
return C(g.q[k]*g.q[k]*g.nu[k],0);

}
);
auto njx=fil([&](size_t k){
return C(g.nu[k]*g.qx[k],0);

}
);
auto njy=fil([&](size_t k){
return C(g.nu[k]*g.qy[k],0);

}
);
auto n2jx=fil([&](size_t k){
return C(g.nu[k]*g.nu[k]*g.qx[k],0);

}
);
auto n2jy=fil([&](size_t k){
return C(g.nu[k]*g.nu[k]*g.qy[k],0);

}
);
V sd(g.z),sk(g.z),sds(g.z),sks(g.z),sdss(g.z),skss(g.z);
for(size_t k=0;
k<g.z;
++k){
sd[k]=2.*v[k]*vn2[k]+vn[k]*vn[k]-hx[k]*hx[k]-hy[k]*hy[k];
sk[k]=2.*v[k]*ro[k]+2.*(hx[k]*jx[k]+hy[k]*jy[k]);
sds[k]=2.*v[k]*vn3[k]+4.*vn[k]*vn2[k]-2.*hx[k]*jx[k]-2.*hy[k]*jy[k];
sks[k]=2.*(vn[k]*ro[k]+v[k]*q2[k])+2.*(jx[k]*jx[k]+hx[k]*njx[k]+jy[k]*jy[k]+hy[k]*njy[k]);
sdss[k]=2.*v[k]*vn4[k]+6.*vn[k]*vn3[k]+4.*vn2[k]*vn2[k]-2.*(jx[k]*jx[k]+hx[k]*njx[k]+jy[k]*jy[k]+hy[k]*njy[k]);
skss[k]=2.*(vn2[k]*ro[k]+2.*vn[k]*q2[k]+v[k]*q2n[k])+6.*(jx[k]*njx[k]+jy[k]*njy[k])+2.*(hx[k]*n2jx[k]+hy[k]*n2jy[k]);

}
sd=g.fft.go(sd,false);
sk=g.fft.go(sk,false);
sds=g.fft.go(sds,false);
sks=g.fft.go(sks,false);
sdss=g.fft.go(sdss,false);
skss=g.fft.go(skss,false);
zeroSuppress(sd,g.pairMask);
zeroSuppress(sk,g.pairMask);
zeroSuppress(sds,g.pairMask);
zeroSuppress(sks,g.pairMask);
zeroSuppress(sdss,g.pairMask);
zeroSuppress(skss,g.pairMask);
double co=g.weights[ni]*std::exp(x)/g.l2*std::exp(-2*g.shift2*t)/4;
for(size_t k=0;
k<g.z;
++k){
double sh=g.sqrtA[k]>0?std::sinh(g.sqrtA[k]*t)/g.sqrtA[k]:t,ch=std::cosh(g.sqrtA[k]*t);
es[k]+=co*(-g.A[k]*sh*sd[k]+ch*sk[k]);
ets[k]-=C(0,1)*co*(-g.A[k]*sh*sds[k]+ch*sks[k]);
ps[k]+=C(0,1)*co*(ch*sd[k]-sh*sk[k]);
pts[k]+=co*(ch*sds[k]-sh*sks[k]);
ptts[k]-=C(0,1)*co*(ch*sdss[k]-sh*skss[k]);

}

}
zeroSuppress(es,g.pairMask);
zeroSuppress(ets,g.pairMask);
zeroSuppress(ps,g.pairMask);
zeroSuppress(pts,g.pairMask);
zeroSuppress(ptts,g.pairMask);
Pair p;
p.eta=g.fft.go(es,true);
p.ex=invMul(g,es,[&](size_t k){
return C(0,g.qx[k]);

}
);
p.ey=invMul(g,es,[&](size_t k){
return C(0,g.qy[k]);

}
);
p.eta_t=g.fft.go(ets,true);
p.ex_t=invMul(g,ets,[&](size_t k){
return C(0,g.qx[k]);

}
);
p.ey_t=invMul(g,ets,[&](size_t k){
return C(0,g.qy[k]);

}
);
p.px=invMul(g,ps,[&](size_t k){
return C(0,g.qx[k]);

}
);
p.py=invMul(g,ps,[&](size_t k){
return C(0,g.qy[k]);

}
);
p.pz=invMul(g,ps,[&](size_t k){
return C(g.A[k],0);

}
);
p.pzz=invMul(g,ps,[&](size_t k){
return C(g.q[k]*g.q[k],0);

}
);
p.ptz=invMul(g,pts,[&](size_t k){
return C(g.A[k],0);

}
);
p.px_t=invMul(g,pts,[&](size_t k){
return C(0,g.qx[k]);

}
);
p.py_t=invMul(g,pts,[&](size_t k){
return C(0,g.qy[k]);

}
);
p.pzz_t=invMul(g,pts,[&](size_t k){
return C(g.q[k]*g.q[k],0);

}
);
p.ptz_t=invMul(g,ptts,[&](size_t k){
return C(g.A[k],0);

}
);
p.phiSpec=ps;
if(order4){
p.pxz=invMul(g,ps,[&](size_t k){
return C(0,g.qx[k]*g.A[k]);

}
);
p.pyz=invMul(g,ps,[&](size_t k){
return C(0,g.qy[k]*g.A[k]);

}
);
p.pzzz=invMul(g,ps,[&](size_t k){
return C(g.q[k]*g.q[k]*g.A[k],0);

}
);

}
return p;

}

struct Third{
V eta,ex,ey,px,py,pz,pzz,pzt;

}
;

Third third(Grid&g,V const&s,int requested=-1){
V es(g.z),ets(g.z),ps(g.z),pts(g.z);
int begin=requested>=0?requested:0;
int end=requested>=0?requested+1:g.rank;
for(int oi=begin;
oi<end;
++oi){
double x=g.nodes[oi],t=x/g.l3;
V damp(g.z);
for(size_t k=0;
k<g.z;
++k)damp[k]=s[k]*std::exp(-t*(g.nu[k]-g.shift3));
First f=first(g,damp);
Pair p=pair(g,damp);
V fk(g.z),fd(g.z),fkt(g.z),fdt(g.z);
for(size_t k=0;
k<g.z;
++k){
C kp=f.px[k]*p.ex[k]+f.py[k]*p.ey[k]+p.px[k]*f.ex[k]+p.py[k]*f.ey[k]-f.eta[k]*p.pzz[k]-p.eta[k]*f.pzz[k];
C kd=f.eta[k]*(f.pxz[k]*f.ex[k]+f.pyz[k]*f.ey[k])-.5*f.eta[k]*f.eta[k]*f.pzzz[k];
fk[k]=kd/4.+kp/2.;
C dp=f.eta[k]*p.ptz[k]+p.eta[k]*f.ptz[k]+f.px[k]*p.px[k]+f.py[k]*p.py[k]+f.pz[k]*p.pz[k];
C dd=.5*f.eta[k]*f.eta[k]*f.ptzz[k]+f.eta[k]*(f.px[k]*f.pxz[k]+f.py[k]*f.pyz[k]+f.pz[k]*f.pzz[k]);
fd[k]=dd/4.+dp/2.;
C kpt=-f.ex[k]*p.ex[k]+f.px[k]*p.ex_t[k]
    -f.ey[k]*p.ey[k]+f.py[k]*p.ey_t[k]
    +p.px_t[k]*f.ex[k]+p.px[k]*f.pxz[k]
    +p.py_t[k]*f.ey[k]+p.py[k]*f.pyz[k]
    -f.pz[k]*p.pzz[k]-f.eta[k]*p.pzz_t[k]
    -p.eta_t[k]*f.pzz[k]-p.eta[k]*f.ptzz[k];
C kdt=f.pz[k]*(f.pxz[k]*f.ex[k]+f.pyz[k]*f.ey[k])+f.eta[k]*(f.pxz_t[k]*f.ex[k]+f.pxz[k]*f.pxz[k]+f.pyz_t[k]*f.ey[k]+f.pyz[k]*f.pyz[k])-f.eta[k]*f.pz[k]*f.pzzz[k]-.5*f.eta[k]*f.eta[k]*f.pzzz_t[k];
fkt[k]=kdt/4.+kpt/2.;
C dpt=f.pz[k]*p.ptz[k]+f.eta[k]*p.ptz_t[k]+p.eta_t[k]*f.ptz[k]+p.eta[k]*f.ptz_t[k]-f.ex[k]*p.px[k]+f.px[k]*p.px_t[k]-f.ey[k]*p.py[k]+f.py[k]*p.py_t[k]+f.ptz[k]*p.pz[k]+f.pz[k]*p.ptz[k];
C ddt=f.eta[k]*f.pz[k]*f.ptzz[k]-.5*f.eta[k]*f.eta[k]*f.pzzz[k]+f.pz[k]*(f.px[k]*f.pxz[k]+f.py[k]*f.pyz[k]+f.pz[k]*f.pzz[k])+f.eta[k]*(-f.ex[k]*f.pxz[k]+f.px[k]*f.pxz_t[k]-f.ey[k]*f.pyz[k]+f.py[k]*f.pyz_t[k]+f.ptz[k]*f.pzz[k]+f.pz[k]*f.ptzz[k]);
fdt[k]=ddt/4.+dpt/2.;

}
fk=g.fft.go(fk,false);
fd=g.fft.go(fd,false);
fkt=g.fft.go(fkt,false);
fdt=g.fft.go(fdt,false);
zeroSuppress(fk,g.tripleMask);
zeroSuppress(fd,g.tripleMask);
zeroSuppress(fkt,g.tripleMask);
zeroSuppress(fdt,g.tripleMask);
double co=g.weights[oi]*std::exp(x)/g.l3*std::exp(-3*g.shift3*t);
for(size_t k=0;
k<g.z;
++k){
double sh=g.sqrtA[k]>0?std::sinh(g.sqrtA[k]*t)/g.sqrtA[k]:t,ch=std::cosh(g.sqrtA[k]*t);
es[k]+=co*(g.A[k]*sh*fd[k]-C(0,1)*ch*fk[k]);
ets[k]+=co*(g.A[k]*sh*fdt[k]-C(0,1)*ch*fkt[k]);
ps[k]+=co*(-C(0,1)*ch*fd[k]-sh*fk[k]);
pts[k]+=co*(-C(0,1)*ch*fdt[k]-sh*fkt[k]);

}

}
zeroSuppress(es,g.tripleMask);
zeroSuppress(ps,g.tripleMask);
zeroSuppress(pts,g.tripleMask);
Third t;
t.eta=g.fft.go(es,true);
t.ex=invMul(g,es,[&](size_t k){
return C(0,g.qx[k]);

}
);
t.ey=invMul(g,es,[&](size_t k){
return C(0,g.qy[k]);

}
);
t.px=invMul(g,ps,[&](size_t k){
return C(0,g.qx[k]);

}
);
t.py=invMul(g,ps,[&](size_t k){
return C(0,g.qy[k]);

}
);
t.pz=invMul(g,ps,[&](size_t k){
return C(g.A[k],0);

}
);
t.pzz=invMul(g,ps,[&](size_t k){
return C(g.q[k]*g.q[k],0);

}
);
t.pzt=invMul(g,pts,[&](size_t k){
return C(g.A[k],0);

}
);
return t;

}

int main(int argc,char**argv){
try{
if(argc!=5&&argc!=6&&argc!=7)throw std::runtime_error("usage INPUT OUTPUT RANK THREADS [OUTER4_INDEX [PART_INDEX]]");
int rank=std::stoi(argv[3]),nt=std::stoi(argv[4]),requested=argc==6?std::stoi(argv[5]):-1;
if(argc==7)requested=std::stoi(argv[5]);
int requestedPart=argc==7?std::stoi(argv[6]):-2;
if(requested>=rank)throw std::runtime_error("outer index");
if(requestedPart>=rank||requestedPart<-2)throw std::runtime_error("part index");
if(argc==7&&requested<0)throw std::runtime_error("split mode requires one outer4 index");
std::ifstream in(argv[1],std::ios::binary);
char mg[8]{

}
;
in.read(mg,8);
if(std::string(mg,8)!="GLRTI001")throw std::runtime_error("magic");
uint64_t n64=0;
rd(in,n64);
int n=int(n64);
double dk=0,d=0,peak=0;
rd(in,dk);
rd(in,d);
rd(in,peak);
uint64_t nr=0;
rd(in,nr);
std::map<int,std::pair<std::vector<double>,std::vector<double>>>rules;
for(uint64_t r=0;
r<nr;
++r){
uint64_t j=0;
rd(in,j);
std::vector<double>no(j),we(j);
for(auto&v:no)rd(in,v);
for(auto&v:we)rd(in,v);
rules[int(j)]={
no,we
}
;

}
if(!rules.count(rank))throw std::runtime_error("rank");
size_t z=size_t(n)*n;
V eta(z);
for(auto&v:eta){
double re,im;
rd(in,re);
rd(in,im);
v={
re,im
}
;

}
auto start=std::chrono::steady_clock::now();
Grid g(n,nt,dk,d,peak,eta,rules[rank].first,rules[rank].second);
int b=requested>=0?requested:0,e=requested>=0?requested+1:rank;
V outSpec(z);
for(int oi=b;
oi<e;
++oi){
double x=g.nodes[oi],tm=x/g.l4;
V damp(z);
for(size_t k=0;
k<z;
++k)damp[k]=eta[k]*std::exp(-tm*(g.nu[k]-g.shift4));
bool include31=requestedPart==-2||requestedPart>=0;
bool includeBase=requestedPart==-2||requestedPart==-1;
Third t;
if(include31)t=third(g,damp,requestedPart>=0?requestedPart:-1);
First f=first(g,damp,includeBase);
Pair p;
if(includeBase)p=pair(g,damp,true);
V fk(z),fd(z);
for(size_t k=0;
k<z;
++k){
C ksum=0.,dsum=0.;
if(include31){
C k31=f.px[k]*t.ex[k]+f.py[k]*t.ey[k]+t.px[k]*f.ex[k]+t.py[k]*f.ey[k]-f.eta[k]*t.pzz[k]-t.eta[k]*f.pzz[k];
C d31=f.eta[k]*t.pzt[k]+t.eta[k]*f.ptz[k]+f.px[k]*t.px[k]+f.py[k]*t.py[k]+f.pz[k]*t.pz[k];
ksum+=k31/2.;
dsum+=d31/2.;
}
if(includeBase){
C k22=p.px[k]*p.ex[k]+p.py[k]*p.ey[k]-p.eta[k]*p.pzz[k];
C d22=p.eta[k]*p.ptz[k]+.5*(p.px[k]*p.px[k]+p.py[k]*p.py[k]+p.pz[k]*p.pz[k]);
C k211=f.eta[k]*(f.pxz[k]*p.ex[k]+f.pyz[k]*p.ey[k])+(f.eta[k]*p.pxz[k]+p.eta[k]*f.pxz[k])*f.ex[k]+(f.eta[k]*p.pyz[k]+p.eta[k]*f.pyz[k])*f.ey[k]-.5*f.eta[k]*f.eta[k]*p.pzzz[k]-f.eta[k]*p.eta[k]*f.pzzz[k];
C d211=.5*f.eta[k]*f.eta[k]*p.pzz_t[k]+f.eta[k]*p.eta[k]*f.ptzz[k]+f.px[k]*(f.eta[k]*p.pxz[k]+p.eta[k]*f.pxz[k])+f.py[k]*(f.eta[k]*p.pyz[k]+p.eta[k]*f.pyz[k])+f.pz[k]*(f.eta[k]*p.pzz[k]+p.eta[k]*f.pzz[k])+p.px[k]*(f.eta[k]*f.pxz[k])+p.py[k]*(f.eta[k]*f.pyz[k])+p.pz[k]*(f.eta[k]*f.pzz[k]);
C k1111=.5*f.eta[k]*f.eta[k]*(f.pxzz[k]*f.ex[k]+f.pyzz[k]*f.ey[k])-f.eta[k]*f.eta[k]*f.eta[k]*f.pzzzz[k]/6.;
C d1111=f.eta[k]*f.eta[k]*f.eta[k]*f.pzzz_t[k]/6.+.5*f.eta[k]*f.eta[k]*(f.px[k]*f.pxzz[k]+f.py[k]*f.pyzz[k]+f.pz[k]*f.pzzz[k])+.5*f.eta[k]*f.eta[k]*(f.pxz[k]*f.pxz[k]+f.pyz[k]*f.pyz[k]+f.pzz[k]*f.pzz[k]);
ksum+=k1111/8.+k211/4.+k22/2.;
dsum+=d1111/8.+d211/4.+d22/2.;
}
fk[k]=ksum;
fd[k]=dsum;

}
fk=g.fft.go(fk,false);
fd=g.fft.go(fd,false);
zeroSuppress(fk,g.quarticMask);
zeroSuppress(fd,g.quarticMask);
double co=g.weights[oi]*std::exp(x)/g.l4*std::exp(-4*g.shift4*tm);
for(size_t k=0;
k<z;
++k){
double sh=g.sqrtA[k]>0?std::sinh(g.sqrtA[k]*tm)/g.sqrtA[k]:tm,ch=std::cosh(g.sqrtA[k]*tm);
outSpec[k]+=co*(g.A[k]*sh*fd[k]-C(0,1)*ch*fk[k]);

}

}
zeroSuppress(outSpec,g.quarticMask);
V field=g.fft.go(outSpec,true);
double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
std::ofstream out(argv[2],std::ios::binary);
out.write("GLE44O01",8);
wr(out,n64);
wr(out,uint64_t(rank));
wr(out,uint64_t(nt));
wr(out,elapsed);
for(auto&v:field){
wr(out,v.real());
wr(out,v.imag());

}
std::cerr<<"CPP_GL_ETA44_DONE n="<<n<<" rank="<<rank<<" threads="<<nt<<" outer="<<requested<<" part="<<requestedPart<<" elapsed="<<elapsed<<"\n";
return 0;

}
catch(std::exception const&e){
std::cerr<<"CPP_GL_ETA44_ERROR "<<e.what()<<"\n";
return 1;

}

}
