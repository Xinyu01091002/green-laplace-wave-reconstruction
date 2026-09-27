#include "second_order_rhs.cuh"
#include <array>
#include <complex>
#include <map>
#include <vector>
#include <iostream>
#include <iomanip>
#include <algorithm>

using Z=std::complex<double>;
using Coeff=std::array<Z,8>; // eta,psi,W,A,Bx,By,eta_t,psi_t
using Fields=std::array<std::vector<double>,8>;
struct Mode {int x,y;Z eta,psi;};
constexpr double tau=6.2831853071795864769;
constexpr double lx=9,ly=7,g=9.81;

// Sparse DIRECT convolution; phase tables avoid recomputing trigonometry
// for every pair at every point. No FFT or GPU operator is used here.
Fields reference(int nx,int ny,double h) {
    std::vector<Mode> m={{0,0,Z(.007),Z(.03)},
        {1,ny==1?0:2,Z(.031,.012),Z(-.017,.027)},
        {3,ny==1?0:-1,Z(-.014,.009),Z(.024,-.011)}};
    for(int j=1;j<=2;++j){auto q=m[j];m.push_back({-q.x,-q.y,std::conj(q.eta),std::conj(q.psi)});}
    auto d=[&](int a,int b){double k=std::hypot(tau*a/lx,tau*b/ly);return k*std::tanh(k*h);};
    std::map<std::pair<int,int>,Coeff> c;
    for(auto q:m) {
        auto& v=c[{q.x,q.y}];v[0]+=q.eta;v[1]+=q.psi;v[2]+=d(q.x,q.y)*q.psi;
        v[6]+=d(q.x,q.y)*q.psi;v[7]-=g*q.eta;
        for(auto p:m) {
            int a=p.x+q.x,b=p.y+q.y;auto& w=c[{a,b}];
            Z A=d(a,b)*d(q.x,q.y)*p.eta*q.psi;
            Z Bx=-tau*tau*a*q.x/(lx*lx)*p.eta*q.psi;
            Z By=-tau*tau*b*q.y/(ly*ly)*p.eta*q.psi;
            w[3]+=A;w[4]+=Bx;w[5]+=By;w[6]-=A+Bx+By;
            w[7]+=0.5*(tau*tau*(p.x*q.x/(lx*lx)+p.y*q.y/(ly*ly))+d(p.x,p.y)*d(q.x,q.y))*p.psi*q.psi;
        }
    }
    Fields f;for(auto& v:f)v.assign(size_t(nx)*ny,0);
    for(auto [mode,v]:c) {
        // Check conjugate symmetry of the independent coefficient reference.
        auto partner=c.find({-mode.first,-mode.second});
        if(partner==c.end())throw std::runtime_error("Missing conjugate reference");
        for(int j=0;j<8;++j)if(std::abs(v[j]-std::conj(partner->second[j]))>1e-13)
            throw std::runtime_error("Non-Hermitian reference");
        std::vector<Z> ex(nx),ey(ny);
        for(int x=0;x<nx;++x)ex[x]=std::exp(Z(0,tau*mode.first*x/nx));
        for(int y=0;y<ny;++y)ey[y]=std::exp(Z(0,tau*mode.second*y/ny));
        for(int y=0;y<ny;++y)for(int x=0;x<nx;++x) {
            Z phase=ex[x]*ey[y];size_t i=size_t(y)*nx+x;
            for(int j=0;j<8;++j)f[j][i]+=(v[j]*phase).real();
        }
    }
    return f;
}
template<class T> double error(const std::vector<T>& a,const std::vector<double>& b) {
    double max=0;for(size_t i=0;i<a.size();++i){
        if(!std::isfinite(a[i]))throw std::runtime_error("Nonfinite field");
        max=std::max(max,std::abs(double(a[i])-b[i]));
    }return max;
}
template<class T> int run(int nx,int ny,double h) {
    auto f=reference(nx,ny,h);size_t n=size_t(nx)*ny;
    std::vector<T> eta(f[0].begin(),f[0].end()),psi(f[1].begin(),f[1].end()),et(n),pt(n),term(n);
    hos::SecondOrderRHS<T> rhs(nx,ny,3,ny==1?0:2,T(lx),T(ly),T(h),T(g));
    rhs.upload(eta.data(),psi.data());
    const char* precision=(std::is_same<T,float>::value)?"fp32":"fp64";
    double tol=(std::is_same<T,float>::value)?2e-5:2e-11;
    bool both_pass=false;
    for(int variant=0;variant<4;++variant) {
        bool input=variant&1,products=variant&2;
        auto prefix=[&](const char* type){
            std::cout<<"{\"type\":\""<<type<<"\",\"precision\":\""<<precision<<"\",\"nx\":"<<nx
              <<",\"ny\":"<<ny<<",\"depth\":"<<h<<",\"variant\":"<<variant;
        };
        auto spectral=[&](const char* label,const typename hos::DeviceFFT<T>::Complex* device,int op,int bx,int by) {
            using C=typename hos::DeviceFFT<T>::Complex;std::vector<C> a(size_t(ny)*(nx/2+1));
            hos::check(cudaMemcpy(a.data(),device,a.size()*sizeof(C),cudaMemcpyDeviceToHost));
            double peak=0,weighted_peak=0,sum=0,weighted_sum=0;
            for(size_t i=0;i<a.size();++i){
                int x=int(i%(nx/2+1)),y=int(i/(nx/2+1));if(y>ny/2)y-=ny;
                if(!std::isfinite(a[i].x)||!std::isfinite(a[i].y))throw std::runtime_error("Nonfinite spectrum");
                if(x<=bx && std::abs(y)<=by)continue;
                double value=std::hypot(double(a[i].x),double(a[i].y))/n;
                double kx=tau*x/lx,ky=tau*y/ly,k=std::hypot(kx,ky);
                double factor=op==0?k*std::tanh(k*h):(op==1?(x==nx/2?0:kx):((ny>1&&std::abs(y)==ny/2)?0:ky));
                double weight=(x==0||x==nx/2)?1:2;
                peak=std::max(peak,value);weighted_peak=std::max(weighted_peak,std::abs(factor)*value);
                sum+=weight*value*value;weighted_sum+=weight*value*value*factor*factor;
            }
            prefix("spectrum");std::cout<<",\"stage\":\""<<label<<"\",\"outside_peak\":"<<peak
              <<",\"outside_rms\":"<<std::sqrt(sum)<<",\"after_operator_outside_peak\":"<<weighted_peak
              <<",\"after_operator_outside_rms\":"<<std::sqrt(weighted_sum)<<"}\n";
        };
        auto field=[&](const char* label,const T* device){
            int j=std::string(label)=="W"?2:(std::string(label)=="A"?3:(std::string(label)=="Bx"?4:5));
            hos::check(cudaMemcpy(term.data(),device,n*sizeof(T),cudaMemcpyDeviceToHost));
            prefix("term");std::cout<<",\"stage\":\""<<label<<"\",\"max_abs\":"<<error(term,f[j])<<"}\n";
        };
        rhs.evaluate(input,products,spectral,field);rhs.download(et.data(),pt.data());
        double ee=error(et,f[6]),pe=error(pt,f[7]);bool pass=ee<tol&&pe<tol;
        prefix("rhs");std::cout<<",\"eta_max_abs\":"<<ee<<",\"psi_max_abs\":"<<pe
          <<",\"pass\":"<<(pass?"true":"false")<<"}\n";
        if(variant==3)both_pass=pass;
        // Debug callbacks must not alter arithmetic, and projected eta must not
        // accumulate inverse-transform roundoff in the original input state.
        auto old_et=et,old_pt=pt;rhs.evaluate(input,products);rhs.download(et.data(),pt.data());
        if(et!=old_et||pt!=old_pt)throw std::runtime_error("Repeat/callback consistency failure");
    }
    return both_pass?0:1; // all variant pass/fail values remain in JSON output
}
int main(int argc,char** argv){
    try {
        if(argc!=5)throw std::invalid_argument("Usage: hos_band_check float|double nx ny h");
        std::string p=argv[1];int nx=std::stoi(argv[2]),ny=std::stoi(argv[3]);double h=std::stod(argv[4]);
        if((p!="float"&&p!="double")||nx<16||(ny!=1&&ny<16)||h<=0||!std::isfinite(h))throw std::invalid_argument("Invalid input");
        int n=0;auto e=cudaGetDeviceCount(&n);if(e!=cudaSuccess||!n){std::cerr<<"SKIP: no GPU\n";return 77;}
        std::cout<<std::setprecision(17);
        return p=="float"?run<float>(nx,ny,h):run<double>(nx,ny,h);
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
