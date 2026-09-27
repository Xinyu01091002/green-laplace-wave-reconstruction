#pragma once
#include "hos_rhs.cuh"
namespace hos {
template<class T> __global__ void stage_state(T* out,const T* state,const T* rate,T dt,size_t n){
    for(size_t i=size_t(blockIdx.x)*blockDim.x+threadIdx.x;i<n;i+=size_t(blockDim.x)*gridDim.x)out[i]=state[i]+dt*rate[i];
}
// Fixed-step RK4 of the explicitly retained-input-band Galerkin system.
// Validation integrator only; not Cash-Karp, not integrating-factor adaptive RK.
template<class T> class RK4 {
    HOSRHS<T>& rhs_;size_t n_;
    Buffer<T> eta_,psi_,stage_eta_,stage_psi_,sum_eta_,sum_psi_;
public:
    explicit RK4(HOSRHS<T>& rhs):rhs_(rhs),n_(rhs.cells()),eta_(n_),psi_(n_),stage_eta_(n_),stage_psi_(n_),sum_eta_(n_),sum_psi_(n_){}
    void upload(const T* eta,const T* psi){
        check(cudaMemcpy(eta_.get(),eta,n_*sizeof(T),cudaMemcpyHostToDevice));
        check(cudaMemcpy(psi_.get(),psi,n_*sizeof(T),cudaMemcpyHostToDevice));
        rhs_.project_state_to_input_band(eta_.get(),psi_.get());
    }
    void step(T dt){
        if(dt<=0||!std::isfinite(dt))throw std::invalid_argument("Invalid RK time step");
        check(cudaMemsetAsync(sum_eta_.get(),0,n_*sizeof(T)));check(cudaMemsetAsync(sum_psi_.get(),0,n_*sizeof(T)));
        rhs_.set_device_state(eta_.get(),psi_.get());
        for(int s=0;s<4;++s){
            rhs_.evaluate();rhs_.project_rhs_to_input_band();
            T weight=(s==0||s==3)?T(1)/T(6):T(1)/T(3);
            axpy<<<256,256>>>(sum_eta_.get(),rhs_.eta_rhs(),weight,n_);
            axpy<<<256,256>>>(sum_psi_.get(),rhs_.psi_rhs(),weight,n_);
            if(s<3){
                T c=s==2?dt:dt/T(2);
                stage_state<<<256,256>>>(stage_eta_.get(),eta_.get(),rhs_.eta_rhs(),c,n_);
                stage_state<<<256,256>>>(stage_psi_.get(),psi_.get(),rhs_.psi_rhs(),c,n_);
                check(cudaGetLastError());rhs_.set_device_state(stage_eta_.get(),stage_psi_.get());
            }
        }
        axpy<<<256,256>>>(eta_.get(),sum_eta_.get(),dt,n_);axpy<<<256,256>>>(psi_.get(),sum_psi_.get(),dt,n_);
        check(cudaGetLastError());rhs_.project_state_to_input_band(eta_.get(),psi_.get());
    }
    void download(T* eta,T* psi){
        check(cudaMemcpy(eta,eta_.get(),n_*sizeof(T),cudaMemcpyDeviceToHost));
        check(cudaMemcpy(psi,psi_.get(),n_*sizeof(T),cudaMemcpyDeviceToHost));
    }
    size_t device_bytes()const{return rhs_.device_bytes()+6*n_*sizeof(T);}
};
}
