#pragma once
#include "second_order_rhs.cuh"
#include <vector>

namespace hos {
template<class T> __global__ void fill(T* a,T value,size_t n) {
    for(size_t i=size_t(blockIdx.x)*blockDim.x+threadIdx.x;i<n;i+=size_t(blockDim.x)*gridDim.x)a[i]=value;
}
template<class T> __global__ void axpy(T* a,const T* b,T scale,size_t n) {
    for(size_t i=size_t(blockIdx.x)*blockDim.x+threadIdx.x;i<n;i+=size_t(blockDim.x)*gridDim.x)a[i]+=scale*b[i];
}
template<class T> __global__ void multiply_add(T* a,const T* b,const T* c,T scale,size_t n) {
    for(size_t i=size_t(blockIdx.x)*blockDim.x+threadIdx.x;i<n;i+=size_t(blockDim.x)*gridDim.x)a[i]+=scale*b[i]*c[i];
}
template<class T> __global__ void power_step(T* a,const T* previous,const T* eta,int j,size_t n) {
    for(size_t i=size_t(blockIdx.x)*blockDim.x+threadIdx.x;i<n;i+=size_t(blockDim.x)*gridDim.x)a[i]=previous[i]*eta[i]/T(j);
}
template<class T,class C> __global__ void vertical_derivative(C* a,int nx,int ny,T lx,T ly,T h,int j) {
    const size_t ns=size_t(ny)*(nx/2+1);
    for(size_t i=size_t(blockIdx.x)*blockDim.x+threadIdx.x;i<ns;i+=size_t(blockDim.x)*gridDim.x){
        int x=int(i%(nx/2+1)),y=int(i/(nx/2+1));if(y>ny/2)y-=ny;
        T kx=T(6.2831853071795864769)*T(x)/lx,ky=T(6.2831853071795864769)*T(y)/ly;
        T k=sqrt(kx*kx+ky*ky),v=T(1);for(int p=0;p<j;++p)v*=k;
        if(j%2)v*=tanh(k*h);
        a[i].x*=v;a[i].y*=v;
    }
}

// Taylor HOS, constant depth, complete homogeneous RHS degrees 1..M.
// All generated order-M modes must fit STRICTLY below the working Nyquist.
// This is a full-support validation path, NOT HOS-Ocean's truncated q=3 masks.
template<class T> class HOSRHS {
    using C=typename DeviceFFT<T>::Complex;
    int nx_,ny_,m_,bx_,by_;T lx_,ly_,h_,g_;size_t n_,ns_;
    DeviceFFT<T> fft_;
    Buffer<T> eta_,psi_,eta_p_,ex_,ey_,px_,py_,grad2_,a_,b_,lower_,et_,pt_;
    Buffer<C> eta_hat_,psi_hat_,phi_;
    Buffer<T> powers_,w_;
    void copy(T* to,const T* from){check(cudaMemcpyAsync(to,from,n_*sizeof(T),cudaMemcpyDeviceToDevice));}
    void copy_spectrum(C* to,const C* from){check(cudaMemcpyAsync(to,from,ns_*sizeof(C),cudaMemcpyDeviceToDevice));}
    void zero(T* a){check(cudaMemsetAsync(a,0,n_*sizeof(T)));}
    void mask(int degree){project_band<<<256,256>>>(fft_.spectral(),nx_,ny_,degree*bx_,degree*by_);check(cudaGetLastError());}
    void project(T* a,int degree){copy(fft_.real(),a);fft_.forward();mask(degree);fft_.inverse();copy(a,fft_.real());}
    C* phi(int degree){return phi_.get()+size_t(degree-1)*ns_;}
    T* w(int degree){return w_.get()+size_t(degree-1)*n_;}
    T* power(int degree){return powers_.get()+size_t(degree)*n_;}
    void derivative(const C* input,int j){
        copy_spectrum(fft_.spectral(),input);
        vertical_derivative<<<256,256>>>(fft_.spectral(),nx_,ny_,lx_,ly_,h_,j);
        check(cudaGetLastError());fft_.inverse();
    }
    void horizontal(const C* input,T* out,int op){
        copy_spectrum(fft_.spectral(),input);
        multiplier<<<256,256>>>(fft_.spectral(),nx_,ny_,lx_,ly_,h_,op);
        check(cudaGetLastError());fft_.inverse();copy(out,fft_.real());
    }
    void w_square(T* out,int degree){
        zero(out);
        for(int p=1;p<degree;++p)multiply_add<<<256,256>>>(out,w(p),w(degree-p),T(1),n_);
        check(cudaGetLastError());project(out,degree);
    }
public:
    using Observer=std::function<void(const char*,int,const T*)>;
    HOSRHS(int nx,int ny,int order,int bx,int by,T lx,T ly,T depth,T gravity)
      :nx_(nx),ny_(ny),m_(order),bx_(bx),by_(by),lx_(lx),ly_(ly),h_(depth),g_(gravity),
       n_(size_t(nx)*ny),ns_(size_t(ny)*(nx/2+1)),fft_(nx,ny),
       eta_(n_),psi_(n_),eta_p_(n_),ex_(n_),ey_(n_),px_(n_),py_(n_),grad2_(n_),
       a_(n_),b_(n_),lower_(n_),et_(n_),pt_(n_),eta_hat_(ns_),psi_hat_(ns_),
       phi_(ns_*size_t(order>0?order:1)),powers_(n_*size_t(order>0?order:1)),w_(n_*size_t(order>0?order:1)) {
        if(order<1||order>5||bx<0||by<0||2LL*order*bx>=nx||(ny==1?by!=0:2LL*order*by>=ny))
            throw std::invalid_argument("Require 1<=M<=5 and complete order-M support below Nyquist");
        if(lx<=0||ly<=0||depth<=0||gravity<=0||!std::isfinite(lx)||!std::isfinite(ly)||!std::isfinite(depth)||!std::isfinite(gravity))
            throw std::invalid_argument("Invalid physical parameters");
    }
    void upload(const T* eta,const T* psi){
        check(cudaMemcpy(eta_.get(),eta,n_*sizeof(T),cudaMemcpyHostToDevice));
        check(cudaMemcpy(psi_.get(),psi,n_*sizeof(T),cudaMemcpyHostToDevice));
    }
    void evaluate(const Observer& observe={}){
        copy(fft_.real(),eta_.get());fft_.forward();mask(1);copy_spectrum(eta_hat_.get(),fft_.spectral());
        fft_.inverse();copy(eta_p_.get(),fft_.real());
        copy(fft_.real(),psi_.get());fft_.forward();mask(1);copy_spectrum(psi_hat_.get(),fft_.spectral());
        horizontal(eta_hat_.get(),ex_.get(),1);horizontal(eta_hat_.get(),ey_.get(),2);
        horizontal(psi_hat_.get(),px_.get(),1);horizontal(psi_hat_.get(),py_.get(),2);
        zero(grad2_.get());
        multiply_add<<<256,256>>>(grad2_.get(),ex_.get(),ex_.get(),T(1),n_);
        multiply_add<<<256,256>>>(grad2_.get(),ey_.get(),ey_.get(),T(1),n_);
        check(cudaGetLastError());if(m_>=2)project(grad2_.get(),2);
        fill<<<256,256>>>(power(0),T(1),n_);check(cudaGetLastError());
        for(int j=1;j<m_;++j){
            power_step<<<256,256>>>(power(j),power(j-1),eta_p_.get(),j,n_);
            check(cudaGetLastError());project(power(j),j);
        }
        copy_spectrum(phi(1),psi_hat_.get());
        for(int degree=1;degree<=m_;++degree){
            if(degree>1){
                zero(a_.get());
                for(int j=1;j<degree;++j){
                    derivative(phi(degree-j),j);
                    multiply_add<<<256,256>>>(a_.get(),power(j),fft_.real(),T(-1),n_);
                    check(cudaGetLastError());
                }
                copy(fft_.real(),a_.get());fft_.forward();mask(degree);copy_spectrum(phi(degree),fft_.spectral());
            }
            zero(w(degree));
            for(int j=0;j<degree;++j){
                derivative(phi(degree-j),j+1);
                multiply_add<<<256,256>>>(w(degree),power(j),fft_.real(),T(1),n_);
                check(cudaGetLastError());
            }
            project(w(degree),degree);if(observe)observe("W",degree,w(degree));
        }
        zero(et_.get());zero(pt_.get());
        for(int degree=1;degree<=m_;++degree){
            copy(a_.get(),w(degree));zero(b_.get());
            if(degree==1)axpy<<<256,256>>>(b_.get(),eta_p_.get(),-g_,n_);
            if(degree==2){
                multiply_add<<<256,256>>>(a_.get(),ex_.get(),px_.get(),T(-1),n_);
                multiply_add<<<256,256>>>(a_.get(),ey_.get(),py_.get(),T(-1),n_);
                multiply_add<<<256,256>>>(b_.get(),px_.get(),px_.get(),T(-0.5),n_);
                multiply_add<<<256,256>>>(b_.get(),py_.get(),py_.get(),T(-0.5),n_);
            }
            if(degree>=3)multiply_add<<<256,256>>>(a_.get(),grad2_.get(),w(degree-2),T(1),n_);
            if(degree>=2){
                w_square(lower_.get(),degree);axpy<<<256,256>>>(b_.get(),lower_.get(),T(0.5),n_);
            }
            if(degree>=4){
                w_square(lower_.get(),degree-2);
                multiply_add<<<256,256>>>(b_.get(),grad2_.get(),lower_.get(),T(0.5),n_);
            }
            check(cudaGetLastError());project(a_.get(),degree);project(b_.get(),degree);
            if(observe){observe("eta",degree,a_.get());observe("psi",degree,b_.get());}
            axpy<<<256,256>>>(et_.get(),a_.get(),T(1),n_);axpy<<<256,256>>>(pt_.get(),b_.get(),T(1),n_);
            check(cudaGetLastError());
        }
    }
    void download(T* eta_t,T* psi_t){
        check(cudaMemcpy(eta_t,et_.get(),n_*sizeof(T),cudaMemcpyDeviceToHost));
        check(cudaMemcpy(psi_t,pt_.get(),n_*sizeof(T),cudaMemcpyDeviceToHost));
    }
    size_t cells() const {return n_;}
    void set_device_state(const T* eta,const T* psi){copy(eta_.get(),eta);copy(psi_.get(),psi);}
    const T* eta_rhs() const {return et_.get();}
    const T* psi_rhs() const {return pt_.get();}
    void project_rhs_to_input_band(){project(et_.get(),1);project(pt_.get(),1);}
    void project_state_to_input_band(T* eta,T* psi){project(eta,1);project(psi,1);}
    size_t device_bytes() const {
        return fft_.array_bytes()+fft_.workspace_bytes()+(13+2*size_t(m_))*n_*sizeof(T)+(2+size_t(m_))*ns_*sizeof(C);
    }
};
}
