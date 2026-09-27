#ifdef HOS_CPU
#include "cpu_generated/rk4.cuh"
#else
#include "rk4.cuh"
#endif
#include <fftw3.h>
#include <chrono>
#include <sys/resource.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <complex>
#include <cstring>
using Clock=std::chrono::steady_clock;
double elapsed(Clock::time_point t){return std::chrono::duration<double>(Clock::now()-t).count();}
double median(std::vector<double> v){std::sort(v.begin(),v.end());return v[v.size()/2];}

// Deterministic 273-positive-mode directional input, with conjugate partners.
// Generated once in FP64 before timing. Exact same sampled arrays for both backends.
template<class T> void input(int nx,int ny,std::vector<T>& eta,std::vector<T>& psi){
    size_t n=size_t(nx)*ny,ns=size_t(ny)*(nx/2+1);double lx=90,ly=60,h=1.3,g=9.81;
    auto* spectrum=fftw_alloc_complex(ns);auto* physical=fftw_alloc_real(n);
    auto plan=fftw_plan_dft_c2r_2d(ny,nx,spectrum,physical,FFTW_ESTIMATE);
    if(!spectrum||!physical||!plan)throw std::runtime_error("Input FFT allocation failure");
    double variance=0;
    for(int x=8;x<=20;++x)for(int y=-10;y<=10;++y){double a=std::exp(-std::pow((x-14)/5.,2)-std::pow(y/8.,2));variance+=2*a*a;}
    double scale=.01/std::sqrt(variance);
    for(int field=0;field<2;++field){
        std::memset(spectrum,0,ns*sizeof(fftw_complex));
        for(int x=8;x<=20;++x)for(int y=-10;y<=10;++y){
            double a=scale*std::exp(-std::pow((x-14)/5.,2)-std::pow(y/8.,2));
            double phase=std::sin(x*17.13+y*9.71)*100.;std::complex<double> c=a*std::exp(std::complex<double>(0,phase));
            if(field){double k=std::hypot(6.283185307179586*x/lx,6.283185307179586*y/ly);double omega=std::sqrt(g*k*std::tanh(k*h));c*=std::complex<double>(0,-g/omega);}
            size_t i=size_t(y<0?ny+y:y)*(nx/2+1)+x;spectrum[i][0]=c.real();spectrum[i][1]=c.imag();
        }
        fftw_execute_dft_c2r(plan,spectrum,physical);
        auto& out=field?psi:eta;for(size_t i=0;i<n;++i)out[i]=T(physical[i]);
    }
    fftw_destroy_plan(plan);fftw_free(spectrum);fftw_free(physical);
}
template<class T> int run(int nx,int ny,int threads,int reps,const char* output){
#ifdef HOS_CPU
    cpu_threads=threads;omp_set_dynamic(0);omp_set_num_threads(threads);const char* backend="cpu_fftw";
#else
    const char* backend="gpu_cufft";
    int devices=0;hos::check(cudaGetDeviceCount(&devices));if(!devices)throw std::runtime_error("No CUDA device");
#endif
    size_t n=size_t(nx)*ny;std::vector<T> eta(n),psi(n),et(n),pt(n),ef(n),pf(n);input(nx,ny,eta,psi);
    auto start=Clock::now();
    hos::HOSRHS<T> rhs(nx,ny,5,20,10,T(90),T(60),T(1.3),T(9.81));
    hos::RK4<T> rk(rhs);
    rhs.upload(eta.data(),psi.data());rhs.evaluate();rhs.download(et.data(),pt.data());hos::check(cudaDeviceSynchronize());
    double cold=elapsed(start);
    rhs.evaluate();hos::check(cudaDeviceSynchronize());
    std::vector<double> times,step_times;
    for(int j=0;j<reps;++j){auto t=Clock::now();rhs.evaluate();hos::check(cudaDeviceSynchronize());times.push_back(elapsed(t));}
    rk.upload(eta.data(),psi.data());rk.step(T(.005));hos::check(cudaDeviceSynchronize());
    for(int j=0;j<reps;++j){
        rk.upload(eta.data(),psi.data());hos::check(cudaDeviceSynchronize());auto t=Clock::now();
        rk.step(T(.005));hos::check(cudaDeviceSynchronize());step_times.push_back(elapsed(t));
    }
    rk.download(ef.data(),pf.data());
    for(auto* v:{&et,&pt,&ef,&pf})for(T x:*v)if(!std::isfinite(x))throw std::runtime_error("Nonfinite output");
    std::ofstream file(output,std::ios::binary);uint64_t count=n;file.write(reinterpret_cast<char*>(&count),sizeof(count));
    // Standard double output permits cross-precision comparisons without scaling.
    for(auto* v:{&et,&pt,&ef,&pf})for(T x:*v){double d=x;file.write(reinterpret_cast<char*>(&d),sizeof(d));}
    file.close();if(!file)throw std::runtime_error("Failed to save benchmark output");
    rusage usage{};getrusage(RUSAGE_SELF,&usage);
    std::cout<<std::setprecision(17)<<"{\"backend\":\""<<backend<<"\",\"precision\":\""<<((std::is_same<T,float>::value)?"fp32":"fp64")
      <<"\",\"nx\":"<<nx<<",\"ny\":"<<ny<<",\"cpu_threads\":"<<threads<<",\"repetitions\":"<<reps
      <<",\"setup_first_rhs_roundtrip_seconds\":"<<cold<<",\"rhs_median_seconds\":"<<median(times)
      <<",\"rk4_step_median_seconds\":"<<median(step_times)<<",\"declared_backend_bytes\":"<<rk.device_bytes()
      <<",\"process_peak_rss_kib\":"<<usage.ru_maxrss<<",\"rhs_samples\":[";
    for(size_t i=0;i<times.size();++i)std::cout<<(i?",":"")<<times[i];std::cout<<"],\"step_samples\":[";
    for(size_t i=0;i<step_times.size();++i)std::cout<<(i?",":"")<<step_times[i];std::cout<<"]}\n";
    return 0;
}
int main(int argc,char** argv){try{
    if(argc!=7)throw std::invalid_argument("Usage: benchmark float|double nx ny threads repetitions output.bin");
    std::string p=argv[1];int nx=std::stoi(argv[2]),ny=std::stoi(argv[3]),threads=std::stoi(argv[4]),reps=std::stoi(argv[5]);
    if((p!="float"&&p!="double")||nx<202||ny<102||nx%2||ny%2||threads<1||reps<3)throw std::invalid_argument("Invalid benchmark inputs");
    return p=="float"?run<float>(nx,ny,threads,reps,argv[6]):run<double>(nx,ny,threads,reps,argv[6]);
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
