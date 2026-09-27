#include "hos_rhs.cuh"
#include <complex>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <array>
#include <algorithm>
using Z=std::complex<double>;
struct Term {int x,y;Z c;};
struct Reference {
    int nx,ny,m,bx,by,sx,sy;double lx,ly,h,g;
    std::vector<std::vector<Term>> fields;
    explicit Reference(const char* path){
        std::ifstream f(path,std::ios::binary);if(!f)throw std::runtime_error("Missing reference");
        int header[7];double params[4];f.read(reinterpret_cast<char*>(header),sizeof(header));
        f.read(reinterpret_cast<char*>(params),sizeof(params));
        nx=header[0];ny=header[1];m=header[2];bx=header[3];by=header[4];sx=header[5];sy=header[6];
        lx=params[0];ly=params[1];h=params[2];g=params[3];
        if(!f||m<1||m>5||sx!=2*m*bx+1||sy!=2*m*by+1||sx>1001||sy>1001)throw std::runtime_error("Invalid reference header");
        for(int j=0;j<2+3*m;++j){
            std::vector<double> re(size_t(sx)*sy),im(re.size());
            f.read(reinterpret_cast<char*>(re.data()),re.size()*sizeof(double));
            f.read(reinterpret_cast<char*>(im.data()),im.size()*sizeof(double));
            if(!f)throw std::runtime_error("Truncated reference");
            std::vector<Term> terms;
            for(int y=0;y<sy;++y)for(int x=0;x<sx;++x){size_t i=size_t(y)*sx+x;
                if(!std::isfinite(re[i])||!std::isfinite(im[i]))throw std::runtime_error("Nonfinite reference");
                if(re[i]!=0||im[i]!=0)terms.push_back({x-sx/2,y-sy/2,Z(re[i],im[i])});
            }fields.push_back(std::move(terms));
        }
    }
};
double synth(const std::vector<Term>& terms,int x,int y,int nx,int ny){
    Z z=0;for(auto t:terms)z+=t.c*std::exp(Z(0,6.2831853071795864769*(double(t.x)*x/nx+double(t.y)*y/ny)));
    if(std::abs(z.imag())>1e-11)throw std::runtime_error("Reference not real");return z.real();
}
template<class T> int run(const Reference& r,int nx,int ny){
    size_t n=size_t(nx)*ny;std::vector<T> eta(n),psi(n),got(n);
    for(int y=0;y<ny;++y)for(int x=0;x<nx;++x){size_t i=size_t(y)*nx+x;
        eta[i]=T(synth(r.fields[0],x,y,nx,ny));psi[i]=T(synth(r.fields[1],x,y,nx,ny));}
    std::vector<size_t> points;
    if(n<=32768){for(size_t i=0;i<n;++i)points.push_back(i);}
    else {for(size_t j=0;j<257;++j)points.push_back((j*7919)%n);}
    hos::HOSRHS<T> solver(nx,ny,r.m,r.bx,r.by,T(r.lx),T(r.ly),T(r.h),T(r.g));
    solver.upload(eta.data(),psi.data());bool all=true;
    std::vector<double> esum(points.size(),0),psum(points.size(),0);
    auto observer=[&](const char* field,int degree,const T* data){
        hos::check(cudaMemcpy(got.data(),data,n*sizeof(T),cudaMemcpyDeviceToHost));
        for(T v:got)if(!std::isfinite(v))throw std::runtime_error("Nonfinite GPU field anywhere on grid");
        int offset=std::string(field)=="W"?0:(std::string(field)=="eta"?1:2);
        double max=0,norm=0,err2=0,refmax=0;
        for(size_t j=0;j<points.size();++j){size_t i=points[j];
            double truth=synth(r.fields[2+3*(degree-1)+offset],int(i%nx),int(i/nx),nx,ny);
            double delta=double(got[i])-truth;max=std::max(max,std::abs(delta));norm+=truth*truth;err2+=delta*delta;refmax=std::max(refmax,std::abs(truth));
            if(offset==1)esum[j]+=truth;if(offset==2)psum[j]+=truth;
        }
        const double atol=(std::is_same<T,float>::value)?1e-10:1e-12;
        const double rtol=(std::is_same<T,float>::value)?5e-4:2e-9;
        bool pass=max<=atol+rtol*refmax;all=all&&pass;
        std::cout<<std::setprecision(17)<<"{\"precision\":\""<<((std::is_same<T,float>::value)?"fp32":"fp64")
          <<"\",\"nx\":"<<nx<<",\"ny\":"<<ny<<",\"depth\":"<<r.h<<",\"degree\":"<<degree
          <<",\"field\":\""<<field<<"\",\"checked_points\":"<<points.size()<<",\"max_abs\":"<<max
          <<",\"relative_l2\":"<<(norm?std::sqrt(err2/norm):0)<<",\"reference_max\":"<<refmax
          <<",\"pass\":"<<(pass?"true":"false")<<"}\n";
    };
    solver.evaluate(observer);std::vector<T> et(n),pt(n);solver.download(et.data(),pt.data());
    double e=0,p=0,emax=0,pmax=0;for(size_t j=0;j<points.size();++j){e=std::max(e,std::abs(double(et[points[j]])-esum[j]));p=std::max(p,std::abs(double(pt[points[j]])-psum[j]));emax=std::max(emax,std::abs(esum[j]));pmax=std::max(pmax,std::abs(psum[j]));}
    double atol=(std::is_same<T,float>::value)?5e-10:5e-12,rtol=(std::is_same<T,float>::value)?5e-4:2e-9;
    all=all&&(e<=atol+rtol*emax)&&(p<=atol+rtol*pmax);
    solver.evaluate();solver.download(eta.data(),psi.data());
    if(eta!=et||psi!=pt)throw std::runtime_error("Repeated HOS evaluation changed output");
    std::cout<<"{\"summary\":true,\"precision\":\""<<((std::is_same<T,float>::value)?"fp32":"fp64")
      <<"\",\"nx\":"<<nx<<",\"ny\":"<<ny<<",\"depth\":"<<r.h<<",\"eta_sum_max_abs\":"<<e
      <<",\"psi_sum_max_abs\":"<<p<<",\"device_bytes\":"<<solver.device_bytes()<<",\"pass\":"<<(all?"true":"false")<<"}\n";
    return all?0:1;
}
int main(int argc,char** argv){try{
    if(argc!=3&&argc!=5)throw std::invalid_argument("Usage: hos_order_check float|double reference.bin [nx ny]");
    std::string p=argv[1];if(p!="float"&&p!="double")throw std::invalid_argument("Unknown precision");
    Reference r(argv[2]);int nx=argc==5?std::stoi(argv[3]):r.nx,ny=argc==5?std::stoi(argv[4]):r.ny;
    if(nx<2||ny<1||nx%2||(ny>1&&ny%2))throw std::invalid_argument("Invalid grid");
    int n=0;auto e=cudaGetDeviceCount(&n);if(e!=cudaSuccess||!n){std::cerr<<"SKIP: no GPU\n";return 77;}
    return p=="float"?run<float>(r,nx,ny):run<double>(r,nx,ny);
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
