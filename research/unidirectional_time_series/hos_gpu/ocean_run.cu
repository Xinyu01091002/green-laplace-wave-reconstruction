// SPDX-License-Identifier: GPL-3.0-or-later
#include "ocean_adaptive.cuh"
#include "ocean_diagnostics.cuh"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <sys/resource.h>
using Clock=std::chrono::steady_clock;
double seconds(Clock::time_point t){return std::chrono::duration<double>(Clock::now()-t).count();}
template<class T> int run(const char* input,const char* output,double duration,double tolerance){
    using C=typename hos::OceanQ3<T>::Complex;using D=cufftDoubleComplex;
    auto start=Clock::now();std::ifstream f(input,std::ios::binary);int h[5];double p[11];
    f.read(reinterpret_cast<char*>(h),sizeof(h));f.read(reinterpret_cast<char*>(p),sizeof(p));
    if(!f||h[0]!=1024||h[1]!=512||h[2]!=5)throw std::runtime_error("Expected verified production-state input");
    size_t n=size_t(h[1])*(h[0]/2+1);std::vector<D> ed(n),pd(n);std::vector<C> e(n),q(n);
    f.read(reinterpret_cast<char*>(ed.data()),n*sizeof(D));f.read(reinterpret_cast<char*>(pd.data()),n*sizeof(D));if(!f)throw std::runtime_error("Truncated input");
    for(size_t i=0;i<n;++i){e[i]={T(ed[i].x),T(ed[i].y)};q[i]={T(pd[i].x),T(pd[i].y)};}
    std::filesystem::path out(output);if(std::filesystem::exists(out))throw std::runtime_error("Never overwrite a run");std::filesystem::create_directory(out);
    hos::OceanCashKarp<T> solver(h[0],h[1],T(p[0]),T(p[1]),T(p[2]),T(p[3]));solver.upload(e.data(),q.data());
    hos::OceanAdaptive<T> controller(solver,duration/p[10],p[5],p[6],p[7],tolerance);
    hos::OceanDiagnostics<T> diagnostics(h[0],h[1],p[3]);double setup=seconds(start),checkpoint_seconds=0,initial_energy=0;
    std::ofstream probes(out/"diagnostics.csv"),trace(out/"steps.csv");probes<<std::setprecision(17);trace<<std::setprecision(17);
    probes<<"time,p1,p2,p3,p4,p5,potential,kinetic,total,relative_energy_change\n";
    trace<<"attempt,t_before,h_proposed,h_actual,error,t_after,accepted,h_next\n";
    auto checkpoint=[&](int index,double time){
        auto t=Clock::now();solver.download_state(e.data(),q.data());
        for(size_t i=0;i<n;++i){ed[i]={double(e[i].x),double(e[i].y)};pd[i]={double(q[i].x),double(q[i].y)};
            if(!std::isfinite(ed[i].x)||!std::isfinite(ed[i].y)||!std::isfinite(pd[i].x)||!std::isfinite(pd[i].y))throw std::runtime_error("Nonfinite state");}
        std::ofstream file(out/("state-"+std::to_string(index)+".bin"),std::ios::binary);uint64_t count=n;int dims[2]={h[0],h[1]};
        file.write(reinterpret_cast<char*>(dims),sizeof(dims));file.write(reinterpret_cast<char*>(&time),8);file.write(reinterpret_cast<char*>(&count),8);
        file.write(reinterpret_cast<char*>(ed.data()),n*sizeof(D));file.write(reinterpret_cast<char*>(pd.data()),n*sizeof(D));file.close();if(!file)throw std::runtime_error("Checkpoint write failed");
        checkpoint_seconds+=seconds(t);
    };
    auto integration_start=Clock::now();int index=0;
    while(true){
        double t=controller.time()*p[10],diag[7];diagnostics.evaluate(solver,diag);
        for(double v:diag)if(!std::isfinite(v))throw std::runtime_error("Nonfinite diagnostic");
        double energy=diag[5]+diag[6];if(index==0)initial_energy=energy;
        double energy_scale=p[9]*p[9]*p[9]/(p[10]*p[10]);
        probes<<t;for(int j=0;j<5;++j)probes<<','<<diag[j]*p[9];
        probes<<','<<diag[5]*energy_scale<<','<<diag[6]*energy_scale<<','<<energy*energy_scale<<','<<(energy-initial_energy)/initial_energy<<'\n';probes.flush();
        if(index==0||std::abs(t-2.0)<1e-7||std::abs(t-13.8)<1e-7||controller.finished())checkpoint(index,t);
        {std::ofstream status(out/"status.json");status<<std::setprecision(17)<<"{\"time\":"<<t<<",\"accepted\":"<<controller.accepted()<<",\"rejected\":"<<controller.rejected()<<",\"wall_seconds\":"<<seconds(integration_start)<<"}\n";}
        if(controller.finished())break;
        controller.advance_output([&](const hos::AdaptiveTrace& a){trace<<a.attempt<<','<<a.time_before<<','<<a.proposal<<','<<a.actual<<','<<a.error<<','<<a.time_after<<','<<a.accepted<<','<<a.next<<'\n';});trace.flush();++index;
    }
    hos::check(cudaDeviceSynchronize());double total=seconds(integration_start),advance=total-checkpoint_seconds;rusage usage{};getrusage(RUSAGE_SELF,&usage);
    std::ofstream report(out/"summary.json");report<<std::setprecision(17)<<"{\"precision\":\""<<((std::is_same<T,float>::value)?"fp32":"fp64")<<"\",\"duration_s\":"<<duration
      <<",\"tolerance\":"<<tolerance<<",\"accepted\":"<<controller.accepted()<<",\"rejected\":"<<controller.rejected()<<",\"setup_seconds\":"<<setup
      <<",\"advance_with_diagnostics_seconds\":"<<advance<<",\"checkpoint_seconds\":"<<checkpoint_seconds<<",\"integration_wall_seconds\":"<<total<<",\"process_peak_rss_kib\":"<<usage.ru_maxrss<<"}\n";
    std::cout<<"Completed "<<output<<" duration="<<duration<<" advance_seconds="<<advance<<'\n';return 0;
}
int main(int argc,char** argv){try{
    if(argc!=6)throw std::runtime_error("Usage: hos_run float|double initial.bin output_dir duration_s tolerance");
    std::string precision=argv[1];double duration=std::stod(argv[4]),tol=std::stod(argv[5]);
    if((precision!="float"&&precision!="double")||duration<=0||tol<=0||!std::isfinite(duration)||!std::isfinite(tol))throw std::runtime_error("Invalid parameters");
    return precision=="float"?run<float>(argv[2],argv[3],duration,tol):run<double>(argv[2],argv[3],duration,tol);
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
