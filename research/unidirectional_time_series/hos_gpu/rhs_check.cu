#include "second_order_rhs.cuh"
#include <complex>
#include <vector>
#include <iostream>
#include <iomanip>
#include <algorithm>

using Z=std::complex<double>;
struct Mode {int x,y; Z eta,psi;};
// Independent sparse Fourier convolution reference: no FFT, no GPU kernels.
std::pair<double,double> reference(const std::vector<Mode>& modes,double x,double y,
                                  double lx,double ly,double h,double g) {
    constexpr double tau=6.2831853071795864769;
    auto d=[&](int a,int b) {double k=std::hypot(tau*a/lx,tau*b/ly);return k*std::tanh(k*h);};
    auto phase=[&](int a,int b) {return std::exp(Z(0,tau*(a*x/lx+b*y/ly)));};
    Z et=0,pt=0;
    for(const auto& q:modes) {
        et+=d(q.x,q.y)*q.psi*phase(q.x,q.y);pt-=g*q.eta*phase(q.x,q.y);
        for(const auto& p:modes) {
            int a=p.x+q.x,b=p.y+q.y;
            double kdotq=tau*tau*(double(a)*q.x/(lx*lx)+double(b)*q.y/(ly*ly));
            double pdotq=tau*tau*(double(p.x)*q.x/(lx*lx)+double(p.y)*q.y/(ly*ly));
            et+=(kdotq-d(a,b)*d(q.x,q.y))*p.eta*q.psi*phase(a,b);
            pt+=0.5*(pdotq+d(p.x,p.y)*d(q.x,q.y))*p.psi*q.psi*phase(a,b);
        }
    }
    if(std::abs(et.imag())>1e-10 || std::abs(pt.imag())>1e-10) throw std::runtime_error("Reference not real");
    return {et.real(),pt.real()};
}
template<class T> int run(int nx,int ny) {
    constexpr double tau=6.2831853071795864769;
    const double lx=9.0,ly=7.0,g=9.81;
    for(double h:{0.15,1.3,20.0}) {
        std::vector<Mode> modes={{0,0,Z(.007),Z(.03)},
          {1,ny==1?0:2,Z(.031,.012),Z(-.017,.027)},
          {3,ny==1?0:-1,Z(-.014,.009),Z(.024,-.011)}};
        for(int j=1;j<=2;++j) {auto m=modes[j];modes.push_back({-m.x,-m.y,std::conj(m.eta),std::conj(m.psi)});}
        size_t n=size_t(nx)*ny;std::vector<T> eta(n),psi(n),et(n),pt(n);
        for(int y=0;y<ny;++y)for(int x=0;x<nx;++x) {
            Z e=0,p=0;for(auto m:modes) {Z z=std::exp(Z(0,tau*(double(m.x)*x/nx+double(m.y)*y/ny)));e+=m.eta*z;p+=m.psi*z;}
            eta[size_t(y)*nx+x]=T(e.real());psi[size_t(y)*nx+x]=T(p.real());
        }
        hos::SecondOrderRHS<T> rhs(nx,ny,3,ny==1?0:2,T(lx),T(ly),T(h),T(g));
        rhs.upload(eta.data(),psi.data());rhs.evaluate();rhs.download(et.data(),pt.data());
        double ee=0,pe=0;
        for(int y=0;y<ny;++y)for(int x=0;x<nx;++x) {
            auto r=reference(modes,lx*x/nx,ly*y/ny,lx,ly,h,g);size_t i=size_t(y)*nx+x;
            if(!std::isfinite(et[i]) || !std::isfinite(pt[i])) throw std::runtime_error("Nonfinite GPU RHS");
            ee=std::max(ee,std::abs(double(et[i])-r.first));pe=std::max(pe,std::abs(double(pt[i])-r.second));
        }
        const double tol=(std::is_same<T,float>::value)?2e-5:2e-11;
        bool pass=ee<tol && pe<tol;
        std::cout<<std::setprecision(17)<<"{\"precision\":\""<<((std::is_same<T,float>::value)?"fp32":"fp64")
          <<"\",\"nx\":"<<nx<<",\"ny\":"<<ny<<",\"depth\":"<<h<<",\"eta_rhs_max_abs\":"<<ee
          <<",\"psi_rhs_max_abs\":"<<pe<<",\"allocated_bytes\":"<<rhs.device_bytes()<<",\"pass\":"<<(pass?"true":"false")<<"}\n";
        if(!pass)return 1;
        // Same inputs evaluated again with the same buffers/plans: catch accidental
        // destruction of psi_hat, stale outputs and accumulation between calls.
        std::vector<T> first_et=et,first_pt=pt;rhs.evaluate();rhs.download(et.data(),pt.data());
        if(et!=first_et || pt!=first_pt)throw std::runtime_error("Repeated evaluation changed result");
    }
    return 0;
}
int main(int argc,char** argv) {
    try {
        if(argc!=4)throw std::invalid_argument("Usage: hos_rhs_check float|double nx ny");
        std::string p=argv[1];if(p!="float" && p!="double")throw std::invalid_argument("Unknown precision");
        int nx=std::stoi(argv[2]),ny=std::stoi(argv[3]);
        if(nx<16 || (ny!=1 && ny<16))throw std::invalid_argument("Test grid too small");
        int count=0;auto s=cudaGetDeviceCount(&count);
        if(s!=cudaSuccess || !count){std::cerr<<"SKIP: no usable GPU\n";return 77;}
        return p=="float"?run<float>(nx,ny):run<double>(nx,ny);
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
