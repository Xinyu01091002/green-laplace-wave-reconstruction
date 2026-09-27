#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <cstdint>
#include <stdexcept>
int main(int argc,char** argv){try{
    if(argc!=4)throw std::runtime_error("Usage: compare reference.bin candidate.bin relative_tolerance");
    std::ifstream ref(argv[1],std::ios::binary),test(argv[2],std::ios::binary);uint64_t n=0,m=0;
    ref.read(reinterpret_cast<char*>(&n),8);test.read(reinterpret_cast<char*>(&m),8);
    if(!ref||!test||n!=m||n==0)throw std::runtime_error("Mismatched files");double tol=std::stod(argv[3]);bool all=true;
    const char* names[]={"eta_rhs","psi_rhs","eta_after_step","psi_after_step"};
    for(int j=0;j<4;++j){double sum=0,norm=0,peak=0;
        for(uint64_t i=0;i<n;++i){double a,b;ref.read(reinterpret_cast<char*>(&a),8);test.read(reinterpret_cast<char*>(&b),8);
            if(!ref||!test||!std::isfinite(a)||!std::isfinite(b))throw std::runtime_error("Invalid field");
            double d=b-a;sum+=d*d;norm+=a*a;peak=std::max(peak,std::abs(d));}
        double rel=std::sqrt(sum/norm);bool pass=rel<=tol;all=all&&pass;
        std::cout<<std::setprecision(17)<<"{\"field\":\""<<names[j]<<"\",\"relative_l2\":"<<rel<<",\"max_abs\":"<<peak<<",\"pass\":"<<(pass?"true":"false")<<"}\n";
    }return all?0:1;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
