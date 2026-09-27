#pragma once
// CPU compatibility layer for mechanically translated, unchanged HOS kernels.
// FFTW MEASURE plans + OpenMP loops; no CUDA runtime is linked into CPU binary.
#include <fftw3.h>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <omp.h>
inline int cpu_threads=1;
struct Dim {int x;};
inline constexpr Dim blockIdx{0},threadIdx{0},blockDim{1},gridDim{1};
enum cudaError_t {cudaSuccess=0,cudaErrorMemoryAllocation=2};
enum cudaMemcpyKind {cudaMemcpyHostToDevice,cudaMemcpyDeviceToHost,cudaMemcpyDeviceToDevice};
inline const char* cudaGetErrorString(cudaError_t x){return x==cudaSuccess?"Success":"CPU allocation failed";}
inline cudaError_t cudaMalloc(void** p,size_t n){*p=fftw_malloc(n);return *p?cudaSuccess:cudaErrorMemoryAllocation;}
inline cudaError_t cudaFree(void* p){fftw_free(p);return cudaSuccess;}
inline cudaError_t cudaMemcpy(void* d,const void* s,size_t n,cudaMemcpyKind){std::memcpy(d,s,n);return cudaSuccess;}
inline cudaError_t cudaMemcpyAsync(void* d,const void* s,size_t n,cudaMemcpyKind k){return cudaMemcpy(d,s,n,k);}
inline cudaError_t cudaMemsetAsync(void* d,int c,size_t n){std::memset(d,c,n);return cudaSuccess;}
inline cudaError_t cudaDeviceSynchronize(){return cudaSuccess;}
inline cudaError_t cudaGetLastError(){return cudaSuccess;}
inline cudaError_t cudaGetDeviceCount(int* n){*n=1;return cudaSuccess;}
enum cufftResult {CUFFT_SUCCESS=0,CUFFT_ALLOC_FAILED=2};
enum cufftType {CUFFT_R2C,CUFFT_C2R,CUFFT_D2Z,CUFFT_Z2D};
struct cufftComplex {float x,y;};struct cufftDoubleComplex {double x,y;};
struct CPUPlan {fftw_plan d=nullptr;fftwf_plan f=nullptr;};
using cufftHandle=CPUPlan*;
inline cufftResult cufftCreate(cufftHandle* p){*p=new CPUPlan;return CUFFT_SUCCESS;}
inline cufftResult cufftDestroy(cufftHandle p){if(p->d)fftw_destroy_plan(p->d);if(p->f)fftwf_destroy_plan(p->f);delete p;return CUFFT_SUCCESS;}
inline cufftResult cufftSetAutoAllocation(cufftHandle,int){return CUFFT_SUCCESS;}
inline cufftResult cufftSetWorkArea(cufftHandle,void*){return CUFFT_SUCCESS;}
inline cufftResult cufftMakePlan2d(cufftHandle p,int ny,int nx,cufftType type,size_t* workspace){
    *workspace=0;size_t n=size_t(nx)*ny,s=size_t(ny)*(nx/2+1);
    if(type==CUFFT_R2C||type==CUFFT_C2R){
        fftwf_init_threads();fftwf_plan_with_nthreads(cpu_threads);
        auto* r=fftwf_alloc_real(n);auto* c=fftwf_alloc_complex(s);
        if(!r||!c){fftwf_free(r);fftwf_free(c);return CUFFT_ALLOC_FAILED;}
        p->f=type==CUFFT_R2C?fftwf_plan_dft_r2c_2d(ny,nx,r,c,FFTW_MEASURE):fftwf_plan_dft_c2r_2d(ny,nx,c,r,FFTW_MEASURE);
        fftwf_free(r);fftwf_free(c);return p->f?CUFFT_SUCCESS:CUFFT_ALLOC_FAILED;
    }
    fftw_init_threads();fftw_plan_with_nthreads(cpu_threads);
    auto* r=fftw_alloc_real(n);auto* c=fftw_alloc_complex(s);
    if(!r||!c){fftw_free(r);fftw_free(c);return CUFFT_ALLOC_FAILED;}
    p->d=type==CUFFT_D2Z?fftw_plan_dft_r2c_2d(ny,nx,r,c,FFTW_MEASURE):fftw_plan_dft_c2r_2d(ny,nx,c,r,FFTW_MEASURE);
    fftw_free(r);fftw_free(c);return p->d?CUFFT_SUCCESS:CUFFT_ALLOC_FAILED;
}
inline cufftResult cufftExecR2C(cufftHandle p,float* r,cufftComplex* c){fftwf_execute_dft_r2c(p->f,r,reinterpret_cast<fftwf_complex*>(c));return CUFFT_SUCCESS;}
inline cufftResult cufftExecC2R(cufftHandle p,cufftComplex* c,float* r){fftwf_execute_dft_c2r(p->f,reinterpret_cast<fftwf_complex*>(c),r);return CUFFT_SUCCESS;}
inline cufftResult cufftExecD2Z(cufftHandle p,double* r,cufftDoubleComplex* c){fftw_execute_dft_r2c(p->d,r,reinterpret_cast<fftw_complex*>(c));return CUFFT_SUCCESS;}
inline cufftResult cufftExecZ2D(cufftHandle p,cufftDoubleComplex* c,double* r){fftw_execute_dft_c2r(p->d,reinterpret_cast<fftw_complex*>(c),r);return CUFFT_SUCCESS;}
