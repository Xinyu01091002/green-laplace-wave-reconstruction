#pragma once
// SPDX-License-Identifier: GPL-3.0-or-later
// CUDA adaptation of the HOS-Ocean v2.1.0 M=5/q=3 flat-bed algorithm.
// Original HOS-Ocean: Copyright (C) 2014 LHEEA Lab., Ecole Centrale de Nantes.
// This adapter follows its MPI normalization, intermediate masks and reduction.
#include "hos_rhs.cuh"
namespace hos {
template<class C> __global__ void ocean_modal_scale(C* a,int nx,int ny,bool to_raw){
    size_t n=size_t(ny)*(nx/2+1);
    for(size_t i=size_t(blockIdx.x)*blockDim.x+threadIdx.x;i<n;i+=size_t(blockDim.x)*gridDim.x){
        int x=int(i%(nx/2+1));double factor=(x==0||x==nx/2)?1.:2.;
        factor=to_raw?double(nx)*ny/factor:factor/(double(nx)*ny);a[i].x*=factor;a[i].y*=factor;
    }
}
template<class T,class C> __global__ void ocean_mean_gravity(C* out,const C* eta,T g){
    if(blockIdx.x==0&&threadIdx.x==0){out[0].x-=g*eta[0].x;out[0].y-=g*eta[0].y;}
}
template<class T> __global__ void scale_half(T* p,size_t n){for(size_t i=size_t(blockIdx.x)*blockDim.x+threadIdx.x;i<n;i+=size_t(blockDim.x)*gridDim.x)p[i]*=T(.5);}
template<class C> __global__ void ocean_extend(C* out,const C* in,int nx,int ny){
    int ex=2*nx,ey=2*ny;size_t n=size_t(ey)*(ex/2+1);
    for(size_t i=size_t(blockIdx.x)*blockDim.x+threadIdx.x;i<n;i+=size_t(blockDim.x)*gridDim.x){
        int x=int(i%(ex/2+1)),y=int(i/(ex/2+1));if(y>ey/2)y-=ey;
        if(x>nx/2||y>ny/2||y< -ny/2){out[i].x=0;out[i].y=0;continue;}
        int iy=y<0?ny+y:y;double scale=4.;if(x==nx/2)scale*=.5;if(y==ny/2||y== -ny/2)scale*=.5;
        auto z=in[size_t(iy)*(nx/2+1)+x];out[i].x=scale*z.x;out[i].y=scale*z.y;
    }
}
template<class C> __global__ void ocean_reduce(C* out,const C* in,int nx,int ny){
    size_t n=size_t(ny)*(nx/2+1);
    for(size_t i=size_t(blockIdx.x)*blockDim.x+threadIdx.x;i<n;i+=size_t(blockDim.x)*gridDim.x){
        int x=int(i%(nx/2+1)),y=int(i/(nx/2+1)),iy=y<=ny/2?y:y+ny;
        double scale=.25;if(x==nx/2)scale*=2;if(y==ny/2)scale*=2;
        auto z=in[size_t(iy)*(nx+1)+x];out[i].x=scale*z.x;out[i].y=scale*z.y;
    }
}
template<class T> class OceanQ3 {
    using C=typename DeviceFFT<T>::Complex;
    int nx_,ny_,ex_,ey_;size_t n_,ne_,se_;T lx_,ly_,h_,g_;
    DeviceFFT<T> base_,ext_;
    Buffer<T> eta_,psi_,eta_base_,psi_base_,px_,py_,gx_,gy_,grad2_,powers_,w_,work_,sum_,a_,b_,cross_,en_,pn_,ef_,pf_;
    Buffer<C> eta_hat_,psi_hat_,phi_;
    void copy(T* dst,const T* src,size_t n){check(cudaMemcpyAsync(dst,src,n*sizeof(T),cudaMemcpyDeviceToDevice));}
    void copy_spectrum(C* dst,const C* src){check(cudaMemcpyAsync(dst,src,se_*sizeof(C),cudaMemcpyDeviceToDevice));}
    void zero(T* p,size_t n){check(cudaMemsetAsync(p,0,n*sizeof(T)));}
    T* power(int j){return powers_.get()+size_t(j)*ne_;}
    T* w(int j){return w_.get()+size_t(j-1)*ne_;}
    C* phi(int j){return phi_.get()+size_t(j-1)*se_;}
    void dealias(T* p,int degree){
        if((degree-1)%2!=0)return;
        copy(ext_.real(),p,ne_);ext_.forward();
        project_band<<<256,256>>>(ext_.spectral(),ex_,ey_,nx_/2,ny_/2);check(cudaGetLastError());
        ext_.inverse();copy(p,ext_.real(),ne_);
    }
    void reduce(const T* p,T* target){
        copy(ext_.real(),p,ne_);ext_.forward();
        ocean_reduce<<<256,256>>>(base_.spectral(),ext_.spectral(),nx_,ny_);check(cudaGetLastError());
        base_.inverse();copy(target,base_.real(),n_);
    }
    void derivative(const C* p,int j){
        copy_spectrum(ext_.spectral(),p);
        vertical_derivative<<<256,256>>>(ext_.spectral(),ex_,ey_,lx_,ly_,h_,j);check(cudaGetLastError());ext_.inverse();
    }
    void horizontal(const C* p,int op,T* target){
        copy_spectrum(ext_.spectral(),p);
        multiplier<<<256,256>>>(ext_.spectral(),ex_,ey_,lx_,ly_,h_,op);check(cudaGetLastError());ext_.inverse();copy(target,ext_.real(),ne_);
    }
    void expand(T* source,T* target,C* spectrum){
        copy(base_.real(),source,n_);base_.forward();
        ocean_extend<<<256,256>>>(ext_.spectral(),base_.spectral(),nx_,ny_);check(cudaGetLastError());
        copy_spectrum(spectrum,ext_.spectral());ext_.inverse();copy(target,ext_.real(),ne_);
    }
public:
    using Complex=C;
    size_t expanded_cells() const {return ne_;}
    const T* slope_x() const {return gx_.get();}
    const T* slope_y() const {return gy_.get();}
    OceanQ3(int nx,int ny,T lx,T ly,T depth,T gravity):nx_(nx),ny_(ny),ex_(2*nx),ey_(2*ny),
       n_(size_t(nx)*ny),ne_(4*n_),se_(size_t(2*ny)*(nx+1)),lx_(lx),ly_(ly),h_(depth),g_(gravity),base_(nx,ny),ext_(2*nx,2*ny),
       eta_(ne_),psi_(ne_),eta_base_(n_),psi_base_(n_),px_(ne_),py_(ne_),gx_(ne_),gy_(ne_),grad2_(ne_),powers_(6*ne_),w_(5*ne_),
       work_(ne_),sum_(ne_),a_(ne_),b_(ne_),cross_(ne_),en_(n_),pn_(n_),ef_(n_),pf_(n_),eta_hat_(se_),psi_hat_(se_),phi_(5*se_){
        if(nx<4||ny<4||nx%2||ny%2||lx<=0||ly<=0||depth<=0||gravity<=0)throw std::invalid_argument("Invalid q3 grid/parameters");
    }
    void upload(const T* eta,const T* psi){
        check(cudaMemcpy(eta_base_.get(),eta,n_*sizeof(T),cudaMemcpyHostToDevice));
        check(cudaMemcpy(psi_base_.get(),psi,n_*sizeof(T),cudaMemcpyHostToDevice));
        expand(eta_base_.get(),eta_.get(),eta_hat_.get());expand(psi_base_.get(),psi_.get(),psi_hat_.get());
    }
    void expanded(T* eta,T* psi){check(cudaMemcpy(eta,eta_.get(),ne_*sizeof(T),cudaMemcpyDeviceToHost));check(cudaMemcpy(psi,psi_.get(),ne_*sizeof(T),cudaMemcpyDeviceToHost));}
    // HOS-Ocean modal convention. Device arrays are not modified by this call.
    void set_modal_state(const C* eta,const C* psi){
        const C* inputs[2]={eta,psi};T* base_fields[2]={eta_base_.get(),psi_base_.get()};
        T* ext_fields[2]={eta_.get(),psi_.get()};C* hats[2]={eta_hat_.get(),psi_hat_.get()};
        size_t ns=size_t(ny_)*(nx_/2+1);
        for(int j=0;j<2;++j){
            check(cudaMemcpyAsync(base_.spectral(),inputs[j],ns*sizeof(C),cudaMemcpyDeviceToDevice));
            ocean_modal_scale<<<256,256>>>(base_.spectral(),nx_,ny_,true);check(cudaGetLastError());
            ocean_extend<<<256,256>>>(ext_.spectral(),base_.spectral(),nx_,ny_);check(cudaGetLastError());
            copy_spectrum(hats[j],ext_.spectral());ext_.inverse();copy(ext_fields[j],ext_.real(),ne_);
            base_.inverse();copy(base_fields[j],base_.real(),n_);
        }
    }
    void modal_nonlinear(C* eta_rhs,C* psi_rhs,const C* eta_state){
        evaluate();size_t ns=size_t(ny_)*(nx_/2+1);
        T* fields[2]={en_.get(),pn_.get()};C* outputs[2]={eta_rhs,psi_rhs};
        for(int j=0;j<2;++j){
            copy(base_.real(),fields[j],n_);base_.forward();
            ocean_modal_scale<<<256,256>>>(base_.spectral(),nx_,ny_,false);check(cudaGetLastError());
            check(cudaMemcpyAsync(outputs[j],base_.spectral(),ns*sizeof(C),cudaMemcpyDeviceToDevice));
        }
        ocean_mean_gravity<<<1,1>>>(psi_rhs,eta_state,g_);check(cudaGetLastError());
    }
    void full_eta_modal(C* target){
        copy(base_.real(),ef_.get(),n_);base_.forward();
        ocean_modal_scale<<<256,256>>>(base_.spectral(),nx_,ny_,false);check(cudaGetLastError());
        check(cudaMemcpyAsync(target,base_.spectral(),size_t(ny_)*(nx_/2+1)*sizeof(C),cudaMemcpyDeviceToDevice));
    }
    void evaluate(){
        fill<<<256,256>>>(power(0),T(1),ne_);copy(power(1),eta_.get(),ne_);
        for(int j=2;j<=5;++j){power_step<<<256,256>>>(power(j),power(j-1),eta_.get(),j,ne_);check(cudaGetLastError());dealias(power(j),j);}
        horizontal(eta_hat_.get(),1,gx_.get());horizontal(eta_hat_.get(),2,gy_.get());
        horizontal(psi_hat_.get(),1,px_.get());horizontal(psi_hat_.get(),2,py_.get());
        zero(grad2_.get(),ne_);multiply_add<<<256,256>>>(grad2_.get(),gx_.get(),gx_.get(),T(1),ne_);
        multiply_add<<<256,256>>>(grad2_.get(),gy_.get(),gy_.get(),T(1),ne_);check(cudaGetLastError());
        copy_spectrum(phi(1),psi_hat_.get());
        for(int degree=1;degree<=5;++degree){
            zero(work_.get(),ne_);zero(w(degree),ne_);
            for(int j=degree;j>=1;--j){
                derivative(phi(degree-j+1),j);
                multiply_add<<<256,256>>>(work_.get(),ext_.real(),power(j),T(-1),ne_);
                multiply_add<<<256,256>>>(w(degree),ext_.real(),power(j-1),T(1),ne_);check(cudaGetLastError());
            }
            if(degree>1)dealias(w(degree),degree);
            if(degree<5){dealias(work_.get(),degree+1);copy(ext_.real(),work_.get(),ne_);ext_.forward();copy_spectrum(phi(degree+1),ext_.spectral());}
        }
        zero(sum_.get(),ne_);for(int j=1;j<=5;++j)axpy<<<256,256>>>(sum_.get(),w(j),T(1),ne_);
        check(cudaGetLastError());reduce(sum_.get(),en_.get());reduce(w(1),ef_.get());
        subtract<<<256,256>>>(en_.get(),ef_.get(),n_);check(cudaGetLastError());
        zero(a_.get(),ne_);multiply_add<<<256,256>>>(a_.get(),gx_.get(),px_.get(),T(-1),ne_);
        multiply_add<<<256,256>>>(a_.get(),gy_.get(),py_.get(),T(-1),ne_);check(cudaGetLastError());
        for(int j=1;j<=3;++j){
            product<<<256,256>>>(work_.get(),grad2_.get(),w(j),ne_);check(cudaGetLastError());dealias(work_.get(),j+2);
            axpy<<<256,256>>>(a_.get(),work_.get(),T(1),ne_);
        }
        check(cudaGetLastError());reduce(a_.get(),ef_.get());axpy<<<256,256>>>(en_.get(),ef_.get(),T(1),n_);
        zero(b_.get(),ne_);multiply_add<<<256,256>>>(b_.get(),px_.get(),px_.get(),T(-1),ne_);
        multiply_add<<<256,256>>>(b_.get(),py_.get(),py_.get(),T(-1),ne_);check(cudaGetLastError());
        for(int p=1;p<=4;++p)for(int q=1;q<=p&&p+q<=5;++q){
            zero(work_.get(),ne_);multiply_add<<<256,256>>>(work_.get(),w(p),w(q),T(p==q?1:2),ne_);check(cudaGetLastError());
            dealias(work_.get(),p+q);axpy<<<256,256>>>(b_.get(),work_.get(),T(1),ne_);
            if(p+q<=3){
                product<<<256,256>>>(cross_.get(),grad2_.get(),work_.get(),ne_);check(cudaGetLastError());dealias(cross_.get(),p+q+2);
                axpy<<<256,256>>>(b_.get(),cross_.get(),T(1),ne_);
            }
        }
        check(cudaGetLastError());reduce(b_.get(),pn_.get());
        scale_half<<<256,256>>>(pn_.get(),n_);check(cudaGetLastError());
        // Full RHS: restore the linear part handled separately by HOS-Ocean RK.
        copy(base_.real(),psi_base_.get(),n_);base_.forward();
        multiplier<<<256,256>>>(base_.spectral(),nx_,ny_,lx_,ly_,h_,0);check(cudaGetLastError());base_.inverse();
        copy(ef_.get(),base_.real(),n_);axpy<<<256,256>>>(ef_.get(),en_.get(),T(1),n_);
        copy(pf_.get(),pn_.get(),n_);axpy<<<256,256>>>(pf_.get(),eta_base_.get(),-g_,n_);check(cudaGetLastError());
    }
    void download(T* en,T* pn,T* ef,T* pf){
        check(cudaMemcpy(en,en_.get(),n_*sizeof(T),cudaMemcpyDeviceToHost));check(cudaMemcpy(pn,pn_.get(),n_*sizeof(T),cudaMemcpyDeviceToHost));
        check(cudaMemcpy(ef,ef_.get(),n_*sizeof(T),cudaMemcpyDeviceToHost));check(cudaMemcpy(pf,pf_.get(),n_*sizeof(T),cudaMemcpyDeviceToHost));
    }
};
}
