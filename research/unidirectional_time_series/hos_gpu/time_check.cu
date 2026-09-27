#include "rk4.cuh"
#include <complex>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
using Z=std::complex<double>;
template<class T> int run(const char* path){
    std::ifstream file(path,std::ios::binary);int head[8];double cfg[5];
    file.read(reinterpret_cast<char*>(head),sizeof(head));file.read(reinterpret_cast<char*>(cfg),sizeof(cfg));
    if(!file)throw std::runtime_error("Invalid time reference");
    int nx=head[0],ny=head[1],m=head[2],bx=head[3],by=head[4],sx=head[5],sy=head[6],steps=head[7];
    if(nx<2||ny<1||sx<1||sy<1||sx>1001||sy>1001||steps<1||steps>10000)throw std::runtime_error("Invalid sizes");
    size_t n=size_t(nx)*ny;std::vector<std::vector<double>> fields;
    for(int j=0;j<4;++j){
        std::vector<double> re(size_t(sx)*sy),im(re.size()),values(n,0);
        file.read(reinterpret_cast<char*>(re.data()),re.size()*sizeof(double));file.read(reinterpret_cast<char*>(im.data()),im.size()*sizeof(double));
        if(!file)throw std::runtime_error("Truncated time reference");
        for(int iy=0;iy<sy;++iy)for(int ix=0;ix<sx;++ix){size_t q=size_t(iy)*sx+ix;
            if(re[q]==0&&im[q]==0)continue;
            for(int y=0;y<ny;++y)for(int x=0;x<nx;++x){
                double phase=6.2831853071795864769*(double(ix-sx/2)*x/nx+double(iy-sy/2)*y/ny);
                values[size_t(y)*nx+x]+=(Z(re[q],im[q])*std::exp(Z(0,phase))).real();
            }
        }fields.push_back(std::move(values));
    }
    std::vector<T> eta(fields[0].begin(),fields[0].end()),psi(fields[1].begin(),fields[1].end());
    hos::HOSRHS<T> rhs(nx,ny,m,bx,by,T(cfg[0]),T(cfg[1]),T(cfg[2]),T(cfg[3]));hos::RK4<T> rk(rhs);
    rk.upload(eta.data(),psi.data());for(int j=0;j<steps;++j)rk.step(T(cfg[4]));rk.download(eta.data(),psi.data());
    double em=0,pm=0,e2=0,p2=0,en=0,pn=0,mean0=0,mean1=0;
    for(size_t i=0;i<n;++i){
        if(!std::isfinite(eta[i])||!std::isfinite(psi[i]))throw std::runtime_error("Nonfinite trajectory");
        double e=double(eta[i])-fields[2][i],p=double(psi[i])-fields[3][i];em=std::max(em,std::abs(e));pm=std::max(pm,std::abs(p));
        e2+=e*e;p2+=p*p;en+=fields[2][i]*fields[2][i];pn+=fields[3][i]*fields[3][i];mean0+=fields[0][i]/n;mean1+=double(eta[i])/n;
    }
    double tol=(std::is_same<T,float>::value)?2e-5:2e-11;bool pass=em<tol&&pm<tol;
    std::cout<<std::setprecision(17)<<"{\"precision\":\""<<((std::is_same<T,float>::value)?"fp32":"fp64")
      <<"\",\"M\":"<<m<<",\"nx\":"<<nx<<",\"ny\":"<<ny<<",\"dt\":"<<cfg[4]<<",\"steps\":"<<steps
      <<",\"eta_max_abs\":"<<em<<",\"psi_max_abs\":"<<pm<<",\"eta_relative_l2\":"<<sqrt(e2/en)
      <<",\"psi_relative_l2\":"<<sqrt(p2/pn)<<",\"eta_mean_change\":"<<mean1-mean0
      <<",\"device_bytes\":"<<rk.device_bytes()<<",\"pass\":"<<(pass?"true":"false")<<"}\n";return pass?0:1;
}
int main(int argc,char** argv){try{
    if(argc!=3)throw std::invalid_argument("Usage: hos_time_check float|double hos_time.bin");
    std::string p=argv[1];if(p!="float"&&p!="double")throw std::invalid_argument("Unknown precision");
    int count=0;auto e=cudaGetDeviceCount(&count);if(e!=cudaSuccess||!count){std::cerr<<"SKIP: no GPU\n";return 77;}
    return p=="float"?run<float>(argv[2]):run<double>(argv[2]);
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
