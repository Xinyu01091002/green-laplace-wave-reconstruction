#include "device_fft.cuh"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

template<class T> int run(int nx,int ny) {
    hos::DeviceFFT<T> fft(nx,ny);
    std::vector<T> input(fft.cells()), output(fft.cells());
    constexpr double pi=3.14159265358979323846;
    // Includes DC and a diagonal mode: catches ordering/sign/normalisation.
    for(int y=0;y<ny;++y) for(int x=0;x<nx;++x)
        input[size_t(y)*nx+x]=T(0.25+0.5*std::cos(2*pi*(double(x)/nx+double(y)/ny)));
    hos::check(cudaMemcpy(fft.real(),input.data(),input.size()*sizeof(T),cudaMemcpyHostToDevice));
    fft.forward();
    std::vector<typename hos::DeviceFFT<T>::Complex> spectrum(size_t(ny)*(nx/2+1));
    hos::check(cudaMemcpy(spectrum.data(),fft.spectral(),spectrum.size()*sizeof(spectrum[0]),cudaMemcpyDeviceToHost));
    double spectral_error=0;
    for(size_t i=0;i<spectrum.size();++i) {
        if(!std::isfinite(spectrum[i].x) || !std::isfinite(spectrum[i].y))
            throw std::runtime_error("Nonfinite spectrum");
        double expected=i==0?0.25: (i==size_t(ny==1?0:1)*(nx/2+1)+1?0.25:0.0);
        spectral_error=std::max(spectral_error,std::hypot(double(spectrum[i].x)/fft.cells()-expected,
                                                       double(spectrum[i].y)/fft.cells()));
    }
    fft.inverse();
    hos::check(cudaMemcpy(output.data(),fft.real(),output.size()*sizeof(T),cudaMemcpyDeviceToHost));
    double max_error=0;
    for(size_t i=0;i<input.size();++i) {
        if(!std::isfinite(output[i])) throw std::runtime_error("Nonfinite inverse output");
        max_error=std::max(max_error,std::abs(double(output[i])-double(input[i])));
    }
    const double tolerance=(std::is_same<T, float>::value) ? 2e-6:2e-13;
    bool ok=max_error<tolerance && spectral_error<tolerance;
    std::cout<<std::setprecision(17)<<"{\"precision\":\""<<((std::is_same<T, float>::value) ? "fp32":"fp64")
        <<"\",\"nx\":"<<nx<<",\"ny\":"<<ny<<",\"array_bytes\":"<<fft.array_bytes()
        <<",\"fft_workspace_bytes\":"<<fft.workspace_bytes()<<",\"roundtrip_max_abs\":"<<max_error
        <<",\"spectrum_max_abs\":"<<spectral_error<<",\"pass\":"<<(ok?"true":"false")<<"}\n";
    return ok?0:1;
}
int main(int argc,char** argv) {
    try {
        if(argc!=4) throw std::invalid_argument("Usage: hos_fft_check float|double nx ny");
        std::string precision=argv[1];
        if(precision!="float" && precision!="double") throw std::invalid_argument("Unknown precision");
        int nx=std::stoi(argv[2]),ny=std::stoi(argv[3]);
        if(nx<8 || (ny!=1 && ny<8)) throw std::invalid_argument("Analytic test requires nx>=8, ny=1 or ny>=8");
        int devices=0; auto status=cudaGetDeviceCount(&devices);
        if(status!=cudaSuccess || devices==0) {
            std::cerr<<"SKIP: no usable CUDA device ("<<cudaGetErrorString(status)<<")\n";
            return 77;
        }
        return precision=="float"?run<float>(nx,ny):run<double>(nx,ny);
    } catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
