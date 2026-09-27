// SPDX-License-Identifier: GPL-3.0-or-later
#include "ocean_cash_karp.cuh"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
int main(int argc,char** argv){try{
    if(argc!=2)throw std::runtime_error("Usage: hos_ocean_ck_check reference.bin");
    std::ifstream f(argv[1],std::ios::binary);int head[5];double par[7];f.read(reinterpret_cast<char*>(head),sizeof(head));f.read(reinterpret_cast<char*>(par),sizeof(par));
    if(!f||head[0]<4||head[1]<4||head[2]!=5)throw std::runtime_error("Invalid reference header");
    size_t n=size_t(head[1])*(head[0]/2+1);using C=cufftDoubleComplex;std::vector<C> ref[4];
    for(auto& a:ref){a.resize(n);f.read(reinterpret_cast<char*>(a.data()),n*sizeof(C));}if(!f)throw std::runtime_error("Truncated reference");
    hos::OceanCashKarp<double> solver(head[0],head[1],par[0],par[1],par[2],par[3]);solver.upload(ref[0].data(),ref[1].data());
    double estimate=solver.attempt(par[4],par[5]);std::vector<C> eta(n),psi(n);solver.download(eta.data(),psi.data());
    bool all=true;double error_delta=std::abs(estimate-par[6]);bool error_pass=error_delta<=5e-15+1e-5*std::abs(par[6]);
    bool decision_match=(estimate<=1e-12)==(par[6]<=1e-12);all=error_pass&&decision_match;
    for(int j=0;j<2;++j){auto& a=j?psi:eta;auto& b=ref[j+2];double sum=0,norm=0,peak=0,refpeak=0;
        for(size_t i=0;i<n;++i){if(!std::isfinite(a[i].x)||!std::isfinite(a[i].y))throw std::runtime_error("Nonfinite step");
            double d=std::hypot(a[i].x-b[i].x,a[i].y-b[i].y),r=std::hypot(b[i].x,b[i].y);sum+=d*d;norm+=r*r;peak=std::max(peak,d);refpeak=std::max(refpeak,r);}
        bool pass=peak<=2e-12+1e-8*refpeak;all=all&&pass;
        std::cout<<std::setprecision(17)<<"{\"nx\":"<<head[0]<<",\"ny\":"<<head[1]<<",\"case\":"<<head[3]<<",\"trial\":"<<head[4]
          <<",\"field\":\""<<(j?"psi":"eta")<<"\",\"relative_l2\":"<<sqrt(sum/norm)<<",\"max_abs\":"<<peak<<",\"pass\":"<<(pass?"true":"false")<<"}\n";
    }
    std::cout<<"{\"nx\":"<<head[0]<<",\"ny\":"<<head[1]<<",\"case\":"<<head[3]<<",\"trial\":"<<head[4]
      <<",\"t0\":"<<par[4]<<",\"step\":"<<par[5]<<",\"gpu_error_estimate\":"<<estimate<<",\"cpu_error_estimate\":"<<par[6]
      <<",\"error_abs_difference\":"<<error_delta<<",\"error_pass\":"<<(error_pass?"true":"false")
      <<",\"decision_match_at_1e_minus12\":"<<(decision_match?"true":"false")<<",\"all_pass\":"<<(all?"true":"false")<<"}\n";
    if(head[4]==2 && par[6]>1e-12){
        // An unaccepted attempt must leave initial state intact, even after
        // a very large rejected candidate. Retry against the original small-step reference.
        std::string path=argv[1];auto pos=path.rfind("_trial2.bin");
        if(pos==std::string::npos)throw std::runtime_error("Cannot locate retry reference");
        path.replace(pos,11,"_trial1.bin");std::ifstream retry_file(path,std::ios::binary);
        int rh[5];double rp[7];retry_file.read(reinterpret_cast<char*>(rh),sizeof(rh));retry_file.read(reinterpret_cast<char*>(rp),sizeof(rp));
        if(!retry_file||rh[0]!=head[0]||rh[1]!=head[1]||rh[3]!=head[3])throw std::runtime_error("Invalid retry reference");
        std::vector<C> rr[4];for(auto& a:rr){a.resize(n);retry_file.read(reinterpret_cast<char*>(a.data()),n*sizeof(C));}
        if(!retry_file)throw std::runtime_error("Truncated retry reference");
        double re=solver.attempt(rp[4],rp[5]);solver.download(eta.data(),psi.data());double max_delta=0,refmax=0;
        for(size_t i=0;i<n;++i){
            max_delta=std::max(max_delta,std::hypot(eta[i].x-rr[2][i].x,eta[i].y-rr[2][i].y));
            max_delta=std::max(max_delta,std::hypot(psi[i].x-rr[3][i].x,psi[i].y-rr[3][i].y));
            refmax=std::max(refmax,std::hypot(rr[2][i].x,rr[2][i].y));refmax=std::max(refmax,std::hypot(rr[3][i].x,rr[3][i].y));
        }
        bool retry_pass=max_delta<=2e-12+1e-8*refmax&&std::abs(re-rp[6])<=5e-15+1e-5*std::abs(rp[6]);all=all&&retry_pass;
        std::cout<<"{\"nx\":"<<head[0]<<",\"ny\":"<<head[1]<<",\"case\":"<<head[3]<<",\"retry_after_rejection\":true,\"max_abs\":"<<max_delta<<",\"error_abs_difference\":"<<std::abs(re-rp[6])<<",\"pass\":"<<(retry_pass?"true":"false")<<"}\n";
    }
    return all?0:1;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
