#pragma once
// SPDX-License-Identifier: GPL-3.0-or-later
// Matches HOS-Ocean runge_kutta.f90 linear-split Cash-Karp stages/error formula.
#include "ocean_q3.cuh"
namespace hos {
template<class T,class C> __device__ void ocean_rotate(C eta_scaled,C psi,T angle,C& e,C& p){
    T c=cos(angle),s=sin(angle);
    e.x=c*eta_scaled.x+s*psi.x;e.y=c*eta_scaled.y+s*psi.y;
    p.x=-s*eta_scaled.x+c*psi.x;p.y=-s*eta_scaled.y+c*psi.y;
}
template<class T,class C> __global__ void ocean_ck_stage(C* out_e,C* out_p,const C* e,const C* p,
        const C* ke,const C* kp,const T* A,const T* b,const T* c,int stage,int nx,int ny,T lx,T ly,T depth,T g,T h){
    size_t n=size_t(ny)*(nx/2+1);
    for(size_t i=size_t(blockIdx.x)*blockDim.x+threadIdx.x;i<n;i+=size_t(blockDim.x)*gridDim.x){
        int x=int(i%(nx/2+1)),y=int(i/(nx/2+1));if(y>ny/2)y-=ny;
        T k=hypot(T(6.2831853071795864769)*x/lx,T(6.2831853071795864769)*y/ly);
        T omega=sqrt(g*k*tanh(k*depth)),go=omega==0?T(1):g/omega;
        C se{0,0},sp{0,0};
        for(int j=0;j<stage;++j){
            C u=ke[size_t(j)*n+i],v=kp[size_t(j)*n+i],a,d;u.x*=go;u.y*=go;
            ocean_rotate(u,v,-omega*h*c[j],a,d);T weight=stage==6?b[j]:A[stage*6+j];
            se.x+=weight*a.x;se.y+=weight*a.y;sp.x+=weight*d.x;sp.y+=weight*d.y;
        }
        C u{e[i].x*go+h*se.x,e[i].y*go+h*se.y},v{p[i].x+h*sp.x,p[i].y+h*sp.y},a,d;
        ocean_rotate(u,v,omega*h*(stage==6?T(1):c[stage]),a,d);
        a.x/=go;a.y/=go;out_e[i]=a;out_p[i]=d;
    }
}
template<class T,class C> __global__ void ocean_ck_error(T* blocks,const C* ke,const C* kp,
        const T* b,const T* low,const T* c,int nx,int ny,T lx,T ly,T depth,T g,T t0,T h,T scale1,T scale2){
    __shared__ T cache[256];size_t n=size_t(ny)*(nx/2+1);T largest=0;
    for(size_t i=size_t(blockIdx.x)*blockDim.x+threadIdx.x;i<n;i+=size_t(blockDim.x)*gridDim.x){
        int x=int(i%(nx/2+1)),y=int(i/(nx/2+1));if(y>ny/2)y-=ny;
        T k=hypot(T(6.2831853071795864769)*x/lx,T(6.2831853071795864769)*y/ly);
        T omega=sqrt(g*k*tanh(k*depth)),go=omega==0?T(1):g/omega;C se{0,0},sp{0,0};
        for(int j=0;j<6;++j){
            C u=ke[size_t(j)*n+i],v=kp[size_t(j)*n+i],a,d;u.x*=go;u.y*=go;
            // Preserve the original error-estimator rotation, including t0.
            ocean_rotate(u,v,omega*(t0+h*c[j]),a,d);T weight=b[j]-low[j];
            se.x+=weight*a.x;se.y+=weight*a.y;sp.x+=weight*d.x;sp.y+=weight*d.y;
        }
        T e1=h*hypot(se.x,se.y)/scale1,e2=h*hypot(sp.x,sp.y)/scale2;
        if(!isfinite(e1)||!isfinite(e2))largest=T(INFINITY);
        else {largest=fmax(largest,e1);largest=fmax(largest,e2);}
    }
    cache[threadIdx.x]=largest;__syncthreads();
    for(int offset=128;offset>0;offset/=2){if(threadIdx.x<offset)cache[threadIdx.x]=fmax(cache[threadIdx.x],cache[threadIdx.x+offset]);__syncthreads();}
    if(threadIdx.x==0)blocks[blockIdx.x]=cache[0];
}
template<class T> __global__ void ocean_ck_max(T* out,const T* blocks){
    __shared__ T cache[256];cache[threadIdx.x]=blocks[threadIdx.x];__syncthreads();
    for(int offset=128;offset>0;offset/=2){if(threadIdx.x<offset)cache[threadIdx.x]=fmax(cache[threadIdx.x],cache[threadIdx.x+offset]);__syncthreads();}
    if(threadIdx.x==0)out[0]=cache[0];
}
template<class T> __global__ void ocean_slope_max(T* blocks,const T* x,const T* y,size_t n){
    __shared__ T cache[256];T largest=0;
    for(size_t i=size_t(blockIdx.x)*blockDim.x+threadIdx.x;i<n;i+=size_t(blockDim.x)*gridDim.x){
        if(!isfinite(x[i])||!isfinite(y[i]))largest=T(INFINITY);
        else {largest=fmax(largest,fabs(x[i]));largest=fmax(largest,fabs(y[i]));}
    }
    cache[threadIdx.x]=largest;__syncthreads();
    for(int offset=128;offset>0;offset/=2){if(threadIdx.x<offset)cache[threadIdx.x]=fmax(cache[threadIdx.x],cache[threadIdx.x+offset]);__syncthreads();}
    if(threadIdx.x==0)blocks[blockIdx.x]=cache[0];
}
template<class T> class OceanCashKarp {
    using C=typename OceanQ3<T>::Complex;
    int nx_,ny_;size_t n_;T lx_,ly_,h_,g_;
    OceanQ3<T> rhs_;
    Buffer<C> eta_,psi_,stage_eta_,stage_psi_,candidate_eta_,candidate_psi_,ke_,kp_,diagnostic_eta_dot_;
    Buffer<T> a_,b_,c_,low_,blocks_,error_;
public:
    OceanCashKarp(int nx,int ny,T lx,T ly,T depth,T gravity):nx_(nx),ny_(ny),n_(size_t(ny)*(nx/2+1)),lx_(lx),ly_(ly),h_(depth),g_(gravity),
        rhs_(nx,ny,lx,ly,depth,gravity),eta_(n_),psi_(n_),stage_eta_(n_),stage_psi_(n_),candidate_eta_(n_),candidate_psi_(n_),ke_(6*n_),kp_(6*n_),diagnostic_eta_dot_(n_),
        a_(36),b_(6),c_(6),low_(6),blocks_(256),error_(1){
        T a[36]={},b[6]={T(37./378),0,T(250./621),T(125./594),0,T(512./1771)};
        T c[6]={0,T(1./5),T(3./10),T(3./5),1,T(7./8)};
        T low[6]={T(2825./27648),0,T(18575./48384),T(13525./55296),T(277./14336),T(1./4)};
        a[6]=T(1./5);a[12]=T(3./40);a[13]=T(9./40);a[18]=T(3./10);a[19]=T(-9./10);a[20]=T(6./5);
        a[24]=T(-11./54);a[25]=T(5./2);a[26]=T(-70./27);a[27]=T(35./27);
        a[30]=T(1631./55296);a[31]=T(175./512);a[32]=T(575./13824);a[33]=T(44275./110592);a[34]=T(253./4096);
        check(cudaMemcpy(a_.get(),a,sizeof(a),cudaMemcpyHostToDevice));check(cudaMemcpy(b_.get(),b,sizeof(b),cudaMemcpyHostToDevice));
        check(cudaMemcpy(c_.get(),c,sizeof(c),cudaMemcpyHostToDevice));check(cudaMemcpy(low_.get(),low,sizeof(low),cudaMemcpyHostToDevice));
    }
    void upload(const C* eta,const C* psi){check(cudaMemcpy(eta_.get(),eta,n_*sizeof(C),cudaMemcpyHostToDevice));check(cudaMemcpy(psi_.get(),psi,n_*sizeof(C),cudaMemcpyHostToDevice));}
    T attempt(T t0,T step,T scale1=T(1),T scale2=T(1)){
        if(step<=0||scale1<=0||scale2<=0||!std::isfinite(t0)||!std::isfinite(step))throw std::invalid_argument("Invalid CK step");
        rhs_.set_modal_state(eta_.get(),psi_.get());rhs_.modal_nonlinear(ke_.get(),kp_.get(),eta_.get());
        for(int stage=1;stage<6;++stage){
            ocean_ck_stage<<<256,256>>>(stage_eta_.get(),stage_psi_.get(),eta_.get(),psi_.get(),ke_.get(),kp_.get(),a_.get(),b_.get(),c_.get(),stage,nx_,ny_,lx_,ly_,h_,g_,step);check(cudaGetLastError());
            rhs_.set_modal_state(stage_eta_.get(),stage_psi_.get());rhs_.modal_nonlinear(ke_.get()+size_t(stage)*n_,kp_.get()+size_t(stage)*n_,stage_eta_.get());
            if(stage==4)rhs_.full_eta_modal(diagnostic_eta_dot_.get());
        }
        ocean_ck_error<<<256,256>>>(blocks_.get(),ke_.get(),kp_.get(),b_.get(),low_.get(),c_.get(),nx_,ny_,lx_,ly_,h_,g_,t0,step,scale1,scale2);check(cudaGetLastError());
        ocean_ck_max<<<1,256>>>(error_.get(),blocks_.get());check(cudaGetLastError());
        ocean_ck_stage<<<256,256>>>(candidate_eta_.get(),candidate_psi_.get(),eta_.get(),psi_.get(),ke_.get(),kp_.get(),a_.get(),b_.get(),c_.get(),6,nx_,ny_,lx_,ly_,h_,g_,step);check(cudaGetLastError());
        T err;check(cudaMemcpy(&err,error_.get(),sizeof(T),cudaMemcpyDeviceToHost));return err;
    }
    void accept(){check(cudaMemcpyAsync(eta_.get(),candidate_eta_.get(),n_*sizeof(C),cudaMemcpyDeviceToDevice));check(cudaMemcpyAsync(psi_.get(),candidate_psi_.get(),n_*sizeof(C),cudaMemcpyDeviceToDevice));}
    void prime(){rhs_.set_modal_state(eta_.get(),psi_.get());rhs_.modal_nonlinear(ke_.get(),kp_.get(),eta_.get());rhs_.full_eta_modal(diagnostic_eta_dot_.get());}
    const C* eta_modes() const {return eta_.get();}
    const C* psi_modes() const {return psi_.get();}
    const C* diagnostic_eta_dot() const {return diagnostic_eta_dot_.get();}
    T current_slope(){
        ocean_slope_max<<<256,256>>>(blocks_.get(),rhs_.slope_x(),rhs_.slope_y(),rhs_.expanded_cells());check(cudaGetLastError());
        ocean_ck_max<<<1,256>>>(error_.get(),blocks_.get());check(cudaGetLastError());
        T value;check(cudaMemcpy(&value,error_.get(),sizeof(T),cudaMemcpyDeviceToHost));return value;
    }
    void download_state(C* eta,C* psi){check(cudaMemcpy(eta,eta_.get(),n_*sizeof(C),cudaMemcpyDeviceToHost));check(cudaMemcpy(psi,psi_.get(),n_*sizeof(C),cudaMemcpyDeviceToHost));}
    void download(C* eta,C* psi){check(cudaMemcpy(eta,candidate_eta_.get(),n_*sizeof(C),cudaMemcpyDeviceToHost));check(cudaMemcpy(psi,candidate_psi_.get(),n_*sizeof(C),cudaMemcpyDeviceToHost));}
};
}
