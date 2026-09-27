// SPDX-License-Identifier: GPL-3.0-or-later
#include "ocean_adaptive.cuh"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
template<class T> int run(const char* path,double tolerance_override){
    using C=typename hos::OceanQ3<T>::Complex;using D=cufftDoubleComplex;
    std::ifstream f(path,std::ios::binary);int h[5];double p[11];f.read(reinterpret_cast<char*>(h),sizeof(h));f.read(reinterpret_cast<char*>(p),sizeof(p));
    if(!f||h[0]<4||h[1]<4||h[2]!=5||h[4]<1||h[4]>100)throw std::runtime_error("Invalid trajectory header");
    size_t n=size_t(h[1])*(h[0]/2+1);std::vector<D> ref_e(n),ref_p(n);std::vector<C> e(n),q(n);
    f.read(reinterpret_cast<char*>(ref_e.data()),n*sizeof(D));f.read(reinterpret_cast<char*>(ref_p.data()),n*sizeof(D));if(!f)throw std::runtime_error("Missing initial state");
    for(size_t i=0;i<n;++i){e[i]={T(ref_e[i].x),T(ref_e[i].y)};q[i]={T(ref_p[i].x),T(ref_p[i].y)};}
    hos::OceanCashKarp<T> solver(h[0],h[1],T(p[0]),T(p[1]),T(p[2]),T(p[3]));solver.upload(e.data(),q.data());
    double tol=tolerance_override>0?tolerance_override:p[8];hos::OceanAdaptive<T> controller(solver,p[4],p[5],p[6],p[7],tol);
    std::vector<hos::AdaptiveTrace> traces;bool all=true;double largest=0;
    const char* precision=(std::is_same<T,float>::value)?"fp32":"fp64";
    for(int output=0;output<h[4];++output){
        double target;f.read(reinterpret_cast<char*>(&target),sizeof(target));
        f.read(reinterpret_cast<char*>(ref_e.data()),n*sizeof(D));f.read(reinterpret_cast<char*>(ref_p.data()),n*sizeof(D));if(!f)throw std::runtime_error("Truncated trajectory");
        controller.advance_output([&](const hos::AdaptiveTrace& t){traces.push_back(t);});solver.download_state(e.data(),q.data());
        if(std::abs(controller.time()-target)>1e-12)throw std::runtime_error("Output time mismatch");
        for(int field=0;field<2;++field){auto& a=field?q:e;auto& b=field?ref_p:ref_e;double err=0,norm=0,peak=0,rmax=0;
            for(size_t i=0;i<n;++i){
                if(!std::isfinite(a[i].x)||!std::isfinite(a[i].y)||!std::isfinite(b[i].x)||!std::isfinite(b[i].y))throw std::runtime_error("Nonfinite trajectory");
                double d=std::hypot(double(a[i].x)-b[i].x,double(a[i].y)-b[i].y),r=std::hypot(b[i].x,b[i].y);
                int x=int(i%(h[0]/2+1));double weight=(x==0||x==h[0]/2)?1.:.5;
                err+=weight*d*d;norm+=weight*r*r;peak=std::max(peak,d);rmax=std::max(rmax,r);
            }
            double rel=std::sqrt(err/norm),atol=(std::is_same<T,float>::value)?2e-6:2e-12,rtol=(std::is_same<T,float>::value)?1e-3:1e-8;
            bool pass=rel<=rtol&&peak<=atol+rtol*rmax;all=all&&pass;largest=std::max(largest,rel);
            std::cout<<std::setprecision(17)<<"{\"type\":\"state\",\"precision\":\""<<precision<<"\",\"case\":"<<h[3]<<",\"nx\":"<<h[0]<<",\"ny\":"<<h[1]
              <<",\"output\":"<<output+1<<",\"time\":"<<target<<",\"physical_time\":"<<target*p[10]<<",\"tolerance\":"<<tol<<",\"field\":\""<<(field?"psi":"eta")
              <<"\",\"physical_l2_relative\":"<<rel<<",\"modal_max_abs\":"<<peak<<",\"pass\":"<<(pass?"true":"false")<<"}\n";
        }
    }
    if(!controller.finished())throw std::runtime_error("Trajectory did not reach requested endpoint");
    std::string trace_path=path;auto dot=trace_path.rfind(".bin");trace_path.replace(dot,4,".trace");std::ifstream ref_trace(trace_path);
    std::string line;std::getline(ref_trace,line);std::vector<hos::AdaptiveTrace> reference;hos::AdaptiveTrace t{};int accept;
    while(ref_trace>>t.attempt>>t.time_before>>t.proposal>>t.actual>>t.error>>t.time_after>>accept>>t.next){t.accepted=accept!=0;reference.push_back(t);}
    if(reference.empty())throw std::runtime_error("Missing CPU step trace");
    bool require_trace=std::is_same<T,double>::value&&tol==p[8];bool match=reference.size()==traces.size();double trace_delta=0;
    if(match)for(size_t i=0;i<traces.size();++i){
        const auto& a=traces[i];const auto& b=reference[i];match=match&&(a.accepted==b.accepted);
        double av[5]={a.time_before,a.proposal,a.actual,a.time_after,a.next},bv[5]={b.time_before,b.proposal,b.actual,b.time_after,b.next};
        for(int j=0;j<5;++j){trace_delta=std::max(trace_delta,std::abs(av[j]-bv[j]));match=match&&std::abs(av[j]-bv[j])<=1e-10+1e-5*std::abs(bv[j]);}
        match=match&&std::abs(a.error-b.error)<=5e-15+1e-5*std::abs(b.error);
    }
    if(require_trace)all=all&&match;
    for(const auto& a:traces)std::cout<<"{\"type\":\"trace\",\"attempt\":"<<a.attempt<<",\"t\":"<<a.time_before<<",\"proposal\":"<<a.proposal<<",\"h\":"<<a.actual<<",\"error\":"<<a.error<<",\"accepted\":"<<(a.accepted?"true":"false")<<",\"after\":"<<a.time_after<<",\"next\":"<<a.next<<"}\n";
    std::cout<<"{\"type\":\"summary\",\"precision\":\""<<precision<<"\",\"case\":"<<h[3]<<",\"nx\":"<<h[0]<<",\"ny\":"<<h[1]<<",\"tolerance\":"<<tol
      <<",\"accepted\":"<<controller.accepted()<<",\"rejected\":"<<controller.rejected()<<",\"reference_attempts\":"<<reference.size()<<",\"trace_required\":"<<(require_trace?"true":"false")
      <<",\"trace_match\":"<<(match?"true":"false")<<",\"max_trace_time_difference\":"<<trace_delta<<",\"max_state_relative_l2\":"<<largest<<",\"pass\":"<<(all?"true":"false")<<"}\n";
    return all?0:1;
}
int main(int argc,char** argv){try{
    if(argc<3||argc>4)throw std::runtime_error("Usage: hos_adaptive_check float|double reference.bin [tolerance]");
    std::string precision=argv[1];double tol=argc==4?std::stod(argv[3]):-1;
    if(precision!="float"&&precision!="double")throw std::runtime_error("Invalid precision");
    return precision=="float"?run<float>(argv[2],tol):run<double>(argv[2],tol);
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
