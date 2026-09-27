#pragma once
#include "device_fft.cuh"
#include <cmath>
#include <functional>

namespace hos {
template<class C> __global__ void project_band(C* a,int nx,int ny,int bx,int by) {
    size_t n=size_t(ny)*(nx/2+1);
    for(size_t i=size_t(blockIdx.x)*blockDim.x+threadIdx.x;i<n;i+=size_t(blockDim.x)*gridDim.x) {
        int x=int(i%(nx/2+1)),y=int(i/(nx/2+1));if(y>ny/2)y-=ny;
        if(x>bx || y>by || y < -by) {a[i].x=0;a[i].y=0;}
    }
}
// Flat-bed finite-depth G0 and horizontal derivatives, physical wavenumbers.
template<class T, class C> __global__ void multiplier(C* a, int nx, int ny,
        T lx,T ly,T depth,int operation) {
    const size_t ns=size_t(ny)*(nx/2+1);
    for(size_t i=size_t(blockIdx.x)*blockDim.x+threadIdx.x;i<ns;i+=size_t(blockDim.x)*gridDim.x) {
        int ix=int(i%(nx/2+1)),iy=int(i/(nx/2+1));
        int sy=iy<=ny/2?iy:iy-ny;
        T kx=T(6.2831853071795864769)*T(ix)/lx;
        T ky=T(6.2831853071795864769)*T(sy)/ly;
        if(operation==0) {
            T k=sqrt(kx*kx+ky*ky),d=k*tanh(k*depth);
            a[i].x*=d;a[i].y*=d;
        } else {
            T k=operation==1?(ix==nx/2?T(0):kx):((ny>1 && iy==ny/2)?T(0):ky);
            T re=a[i].x;a[i].x=-k*a[i].y;a[i].y=k*re;
        }
    }
}
template<class T> __global__ void product(T* out,const T* a,const T* b,size_t n) {
    for(size_t i=size_t(blockIdx.x)*blockDim.x+threadIdx.x;i<n;i+=size_t(blockDim.x)*gridDim.x)
        out[i]=a[i]*b[i];
}
template<class T> __global__ void subtract(T* out,const T* a,size_t n) {
    for(size_t i=size_t(blockIdx.x)*blockDim.x+threadIdx.x;i<n;i+=size_t(blockDim.x)*gridDim.x)
        out[i]-=a[i];
}
template<class T> __global__ void dynamic_rhs(T* out,const T* eta,const T* w,
        const T* px,const T* py,T gravity,size_t n) {
    for(size_t i=size_t(blockIdx.x)*blockDim.x+threadIdx.x;i<n;i+=size_t(blockDim.x)*gridDim.x)
        out[i]=-gravity*eta[i]+T(0.5)*(w[i]*w[i]-px[i]*px[i]-py[i]*py[i]);
}

// Complete quadratic RHS in eta and TRUE surface potential psi.
// Caller supplies band-limited fields on an expanded grid: 4*bx<nx,
// 4*by<ny (ny=1 allowed). Optional projections use only declared support.
// This is a standalone order-2 gate, not the production M=5 HOS solver.
template<class T> class SecondOrderRHS {
    using C=typename DeviceFFT<T>::Complex;
    int nx_,ny_,bx_,by_; T lx_,ly_,depth_,g_;
    DeviceFFT<T> fft_;
    Buffer<C> psi_hat_;
    Buffer<T> eta_,psi_,w_,px_,py_,eta_t_,psi_t_,eta_projected_;
public:
    // Optional host diagnostic callbacks; absent in the ordinary device-only path.
    using SpectralObserver=std::function<void(const char*,const C*,int,int,int)>;
    using FieldObserver=std::function<void(const char*,const T*)>;
private:
    void copy(T* to,const T* from) { check(cudaMemcpyAsync(to,from,eta_.bytes(),cudaMemcpyDeviceToDevice)); }
    void apply(int op) {
        multiplier<<<256,256>>>(fft_.spectral(),nx_,ny_,lx_,ly_,depth_,op);
        check(cudaGetLastError());fft_.inverse();
    }
    void psi_operator(T* target,int op) {
        check(cudaMemcpyAsync(fft_.spectral(),psi_hat_.get(),psi_hat_.bytes(),cudaMemcpyDeviceToDevice));
        apply(op);copy(target,fft_.real());
    }
    void subtract_product_derivative(const T* eta,const T* field,int op,bool project,
                                    const SpectralObserver& spectral,const FieldObserver& real) {
        product<<<256,256>>>(fft_.real(),eta,field,fft_.cells());
        check(cudaGetLastError());fft_.forward();
        const char* name=op==0?"A":(op==1?"Bx":"By");
        if(spectral)spectral(name,fft_.spectral(),op,2*bx_,2*by_);
        if(project) {project_band<<<256,256>>>(fft_.spectral(),nx_,ny_,2*bx_,2*by_);check(cudaGetLastError());}
        apply(op);
        if(real)real(name,fft_.real());
        subtract<<<256,256>>>(eta_t_.get(),fft_.real(),fft_.cells());check(cudaGetLastError());
    }
public:
    SecondOrderRHS(int nx,int ny,int bx,int by,T lx,T ly,T depth,T gravity)
        :nx_(nx),ny_(ny),bx_(bx),by_(by),lx_(lx),ly_(ly),depth_(depth),g_(gravity),fft_(nx,ny),
        psi_hat_(size_t(ny)*(nx/2+1)),eta_(fft_.cells()),psi_(fft_.cells()),
        w_(fft_.cells()),px_(fft_.cells()),py_(fft_.cells()),eta_t_(fft_.cells()),psi_t_(fft_.cells()),eta_projected_(fft_.cells()) {
        if(bx<0 || by<0 || 4LL*bx>=nx || (ny==1?by!=0:4LL*by>=ny))
            throw std::invalid_argument("Quadratic support must fit strictly below Nyquist");
        if(!std::isfinite(lx) || !std::isfinite(ly) || !std::isfinite(depth) || !std::isfinite(gravity)
            || lx<=0 || ly<=0 || depth<=0 || gravity<=0) throw std::invalid_argument("Invalid physical parameters");
    }
    void upload(const T* eta,const T* psi) {
        check(cudaMemcpy(eta_.get(),eta,eta_.bytes(),cudaMemcpyHostToDevice));
        check(cudaMemcpy(psi_.get(),psi,psi_.bytes(),cudaMemcpyHostToDevice));
    }
    void evaluate(bool project_input=false,bool project_products=false,
                  const SpectralObserver& spectral={},const FieldObserver& real={}) {
        const T* eta=eta_.get();
        if(project_input) {
            copy(fft_.real(),eta_.get());fft_.forward();
            project_band<<<256,256>>>(fft_.spectral(),nx_,ny_,bx_,by_);check(cudaGetLastError());
            fft_.inverse();copy(eta_projected_.get(),fft_.real());eta=eta_projected_.get();
        }
        copy(fft_.real(),psi_.get());fft_.forward();
        if(spectral)spectral("psi_input",fft_.spectral(),0,bx_,by_);
        if(project_input) {project_band<<<256,256>>>(fft_.spectral(),nx_,ny_,bx_,by_);check(cudaGetLastError());}
        check(cudaMemcpyAsync(psi_hat_.get(),fft_.spectral(),psi_hat_.bytes(),cudaMemcpyDeviceToDevice));
        psi_operator(w_.get(),0);psi_operator(px_.get(),1);psi_operator(py_.get(),2);
        if(real)real("W",w_.get());
        copy(eta_t_.get(),w_.get());
        subtract_product_derivative(eta,w_.get(),0,project_products,spectral,real);
        subtract_product_derivative(eta,px_.get(),1,project_products,spectral,real);
        subtract_product_derivative(eta,py_.get(),2,project_products,spectral,real);
        dynamic_rhs<<<256,256>>>(psi_t_.get(),eta,w_.get(),px_.get(),py_.get(),g_,fft_.cells());
        check(cudaGetLastError());
    }
    void download(T* eta_t,T* psi_t) {
        check(cudaMemcpy(eta_t,eta_t_.get(),eta_t_.bytes(),cudaMemcpyDeviceToHost));
        check(cudaMemcpy(psi_t,psi_t_.get(),psi_t_.bytes(),cudaMemcpyDeviceToHost));
    }
    size_t device_bytes() const {
        return fft_.array_bytes()+fft_.workspace_bytes()+psi_hat_.bytes()+8*eta_.bytes();
    }
};
}
