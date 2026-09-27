#pragma once
#include <cuda_runtime.h>
#include <cufft.h>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>

namespace hos {
inline void check(cudaError_t s) {
    if (s != cudaSuccess) throw std::runtime_error(cudaGetErrorString(s));
}
inline void check(cufftResult s) {
    if (s != CUFFT_SUCCESS)
        throw std::runtime_error("cuFFT status " + std::to_string(static_cast<int>(s)));
}
template<class T> class Buffer {
    T* p_ = nullptr;
    size_t count_;
public:
    explicit Buffer(size_t n) : count_(n) {
        if (!n || n > std::numeric_limits<size_t>::max()/sizeof(T))
            throw std::invalid_argument("Invalid device allocation size");
        check(cudaMalloc(reinterpret_cast<void**>(&p_), bytes()));
    }
    ~Buffer() { if (p_) cudaFree(p_); }
    Buffer(const Buffer&) = delete;
    Buffer& operator=(const Buffer&) = delete;
    T* get() const { return p_; }
    size_t bytes() const { return count_*sizeof(T); }
};
class Plan {
    cufftHandle p_ = 0;
public:
    Plan() { check(cufftCreate(&p_)); }
    ~Plan() { if (p_) cufftDestroy(p_); }
    Plan(const Plan&) = delete;
    Plan& operator=(const Plan&) = delete;
    operator cufftHandle() const { return p_; }
};
template<class T> __global__ void normalize(T* x, size_t n) {
    for (size_t i = size_t(blockIdx.x)*blockDim.x + threadIdx.x;
         i < n; i += size_t(blockDim.x)*gridDim.x) x[i] /= T(n);
}

// Row-major physical layout [ny][nx]; R2C half-spectrum [ny][nx/2+1].
// Forward coefficients are unnormalised; inverse applies 1/(nx*ny).
// One instance is sequential on the default stream, not concurrently reusable.
// C2R may overwrite its spectral input. Callers must preserve state separately.
template<class T> class DeviceFFT {
    static_assert(std::is_same_v<T,float> || std::is_same_v<T,double>);
public:
    using Complex = std::conditional_t<std::is_same_v<T,float>,cufftComplex,cufftDoubleComplex>;
private:
    static size_t cells(int nx, int ny) {
        if (nx < 2 || ny < 1 || nx%2 || (ny != 1 && ny%2))
            throw std::invalid_argument("Require even nx>=2, and ny=1 or even ny>=2");
        return size_t(nx)*size_t(ny);
    }
    size_t n_, work_bytes_ = 0;
    Buffer<T> real_;
    Buffer<Complex> spectral_;
    // Plans destroyed before their workspace.
    std::unique_ptr<Buffer<unsigned char>> work_;
    Plan forward_, inverse_;
public:
    DeviceFFT(int nx, int ny) : n_(cells(nx,ny)), real_(n_),
        spectral_(size_t(ny)*(nx/2+1)) {
        size_t a=0,b=0;
        check(cufftSetAutoAllocation(forward_,0));
        check(cufftSetAutoAllocation(inverse_,0));
        check(cufftMakePlan2d(forward_,ny,nx,
              (std::is_same<T, float>::value) ? CUFFT_R2C:CUFFT_D2Z,&a));
        check(cufftMakePlan2d(inverse_,ny,nx,
              (std::is_same<T, float>::value) ? CUFFT_C2R:CUFFT_Z2D,&b));
        work_bytes_ = a>b?a:b;
        work_ = std::make_unique<Buffer<unsigned char>>(work_bytes_?work_bytes_:1);
        check(cufftSetWorkArea(forward_,work_->get()));
        check(cufftSetWorkArea(inverse_,work_->get()));
    }
    ~DeviceFFT() { cudaDeviceSynchronize(); }
    DeviceFFT(const DeviceFFT&) = delete;
    DeviceFFT& operator=(const DeviceFFT&) = delete;
    T* real() const { return real_.get(); }
    Complex* spectral() const { return spectral_.get(); }
    size_t cells() const { return n_; }
    size_t array_bytes() const { return real_.bytes()+spectral_.bytes(); }
    size_t workspace_bytes() const { return work_bytes_; }
    void forward() {
        if constexpr(std::is_same_v<T,float>)
            check(cufftExecR2C(forward_,real(),spectral()));
        else check(cufftExecD2Z(forward_,real(),spectral()));
    }
    void inverse() {
        if constexpr(std::is_same_v<T,float>)
            check(cufftExecC2R(inverse_,spectral(),real()));
        else check(cufftExecZ2D(inverse_,spectral(),real()));
        normalize<<<256,256>>>(real(),n_);
        check(cudaGetLastError());
    }
};
}
