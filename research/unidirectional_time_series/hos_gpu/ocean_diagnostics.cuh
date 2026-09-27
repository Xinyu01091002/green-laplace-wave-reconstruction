// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "ocean_cash_karp.cuh"
namespace hos {
// Diagnostics accumulate in FP64 regardless of state precision. They do not
// feed back into evolution. Energy uses the original c=1 staged eta derivative.
template<class C> __global__ void ocean_diag_blocks(double* out,const C* eta,const C* psi,const C* deta,int nx,int ny,double gravity){
    __shared__ double cache[256];int channel=blockIdx.y;double sum=0;
    const size_t n=size_t(ny)*(nx/2+1);
    const int offsets[5]={0,-6,6,-8,8};
    for(size_t i=size_t(blockIdx.x)*blockDim.x+threadIdx.x;i<n;i+=size_t(blockDim.x)*gridDim.x){
        int x=int(i%(nx/2+1)),y=int(i/(nx/2+1));if(y>ny/2)y-=ny;
        double er=eta[i].x,ei=eta[i].y;
        if(channel<5){
            double phase=6.2831853071795864769*(.5*x+double(ny/2+offsets[channel])*y/ny);
            sum+=er*cos(phase)-ei*sin(phase);
        }else{
            double weight=x==0?1.:.5;if(y==ny/2)weight*=.5;
            if(channel==5)sum+=.5*gravity*weight*(er*er+ei*ei);
            else sum+=.5*weight*(double(psi[i].x)*deta[i].x+double(psi[i].y)*deta[i].y);
        }
    }
    cache[threadIdx.x]=sum;__syncthreads();
    for(int offset=128;offset>0;offset/=2){if(threadIdx.x<offset)cache[threadIdx.x]+=cache[threadIdx.x+offset];__syncthreads();}
    if(threadIdx.x==0)out[channel*256+blockIdx.x]=cache[0];
}
__global__ void ocean_diag_reduce(double* out,const double* blocks){
    __shared__ double cache[256];int channel=blockIdx.x;cache[threadIdx.x]=blocks[channel*256+threadIdx.x];__syncthreads();
    for(int offset=128;offset>0;offset/=2){if(threadIdx.x<offset)cache[threadIdx.x]+=cache[threadIdx.x+offset];__syncthreads();}
    if(threadIdx.x==0)out[channel]=cache[0];
}
template<class T> class OceanDiagnostics {
    Buffer<double> blocks_{7*256},out_{7};int nx_,ny_;double g_;
public:
    OceanDiagnostics(int nx,int ny,double g):nx_(nx),ny_(ny),g_(g){}
    void evaluate(const OceanCashKarp<T>& solver,double* out){
        ocean_diag_blocks<<<dim3(256,7),256>>>(blocks_.get(),solver.eta_modes(),solver.psi_modes(),solver.diagnostic_eta_dot(),nx_,ny_,g_);check(cudaGetLastError());
        ocean_diag_reduce<<<7,256>>>(out_.get(),blocks_.get());check(cudaGetLastError());check(cudaMemcpy(out,out_.get(),7*sizeof(double),cudaMemcpyDeviceToHost));
    }
};
}
