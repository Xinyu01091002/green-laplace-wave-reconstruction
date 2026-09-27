// SPDX-License-Identifier: GPL-3.0-or-later
#include "ocean_q3.cuh"
#include <fstream>
#include <vector>
#include <iostream>
#include <iomanip>
#include <algorithm>
template<class T> int run(const char* path){
    std::ifstream f(path,std::ios::binary);int header[6];double p[4];f.read(reinterpret_cast<char*>(header),sizeof(header));f.read(reinterpret_cast<char*>(p),sizeof(p));
    if(!f||header[2]!=2*header[0]||header[3]!=2*header[1]||header[4]!=5)throw std::runtime_error("Invalid official reference");
    int nx=header[0],ny=header[1];size_t n=size_t(nx)*ny,ne=4*n;
    std::vector<std::vector<double>> ref;
    for(int j=0;j<8;++j){ref.emplace_back(j<2?ne:n);f.read(reinterpret_cast<char*>(ref.back().data()),ref.back().size()*sizeof(double));}
    if(!f)throw std::runtime_error("Truncated reference");
    std::vector<T> input_eta(ref[2].begin(),ref[2].end()),input_psi(ref[3].begin(),ref[3].end());
    hos::OceanQ3<T> solver(nx,ny,T(p[0]),T(p[1]),T(p[2]),T(p[3]));solver.upload(input_eta.data(),input_psi.data());
    std::vector<std::vector<T>> got(6);got[0].resize(ne);got[1].resize(ne);for(int j=2;j<6;++j)got[j].resize(n);
    solver.expanded(got[0].data(),got[1].data());solver.evaluate();solver.download(got[2].data(),got[3].data(),got[4].data(),got[5].data());
    const char* labels[]={"expanded_eta","expanded_psi","eta_nonlinear","psi_nonlinear","eta_full","psi_full"};bool all=true;
    for(int j=0;j<6;++j){auto& r=ref[j<2?j:j+2];double err=0,norm=0,peak=0,rmax=0;
        for(size_t i=0;i<r.size();++i){if(!std::isfinite(got[j][i])||!std::isfinite(r[i]))throw std::runtime_error("Nonfinite field");
            double d=got[j][i]-r[i];err+=d*d;norm+=r[i]*r[i];peak=std::max(peak,std::abs(d));rmax=std::max(rmax,std::abs(r[i]));}
        double rtol=(std::is_same<T,float>::value)?5e-4:1e-8,atol=(std::is_same<T,float>::value)?2e-6:2e-12;
        double rel=norm?std::sqrt(err/norm):0;
        bool pass=peak<=atol+rtol*rmax&&(norm==0||rel<=rtol);all=all&&pass;
        std::cout<<std::setprecision(17)<<"{\"precision\":\""<<((std::is_same<T,float>::value)?"fp32":"fp64")<<"\",\"nx\":"<<nx<<",\"ny\":"<<ny<<",\"case\":"<<header[5]<<",\"field\":\""<<labels[j]<<"\",\"relative_l2\":"<<rel<<",\"max_abs\":"<<peak<<",\"reference_max\":"<<rmax<<",\"pass\":"<<(pass?"true":"false")<<"}\n";
    }return all?0:1;
}
int main(int argc,char** argv){try{
    if(argc==2)return run<double>(argv[1]);
    if(argc!=3)throw std::runtime_error("Usage: hos_ocean_check [float|double] official_reference.bin");
    std::string precision=argv[1];if(precision!="float"&&precision!="double")throw std::runtime_error("Unknown precision");
    return precision=="float"?run<float>(argv[2]):run<double>(argv[2]);
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
