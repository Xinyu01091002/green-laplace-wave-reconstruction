#include <algorithm>
#include <complex>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
using C=std::complex<double>;
template<class T>void rd(std::ifstream&f,T&v){f.read(reinterpret_cast<char*>(&v),sizeof(T));if(!f)throw std::runtime_error("input eof");}
template<class T>void wr(std::ofstream&f,const T&v){f.write(reinterpret_cast<const char*>(&v),sizeof(T));if(!f)throw std::runtime_error("output failed");}
int main(int argc,char**argv){try{if(argc<4)throw std::runtime_error("usage OUTPUT INPUT1 INPUT2 ...");uint64_t n=0,rank=0;double wall=0;std::vector<C>sum;for(int j=2;j<argc;++j){std::ifstream in(argv[j],std::ios::binary);char magic[8]{};in.read(magic,8);if(std::string(magic,8)!="GLE44O01")throw std::runtime_error("magic");uint64_t cn=0,cr=0,ct=0;double elapsed=0;rd(in,cn);rd(in,cr);rd(in,ct);rd(in,elapsed);if(j==2){n=cn;rank=cr;sum.assign(size_t(n)*n,C(0));}else if(cn!=n||cr!=rank)throw std::runtime_error("incompatible partial outputs");wall=std::max(wall,elapsed);for(auto&v:sum){double re=0,im=0;rd(in,re);rd(in,im);v+=C(re,im);}}std::ofstream out(argv[1],std::ios::binary);out.write("GLE44O01",8);wr(out,n);wr(out,rank);wr(out,uint64_t(argc-2));wr(out,wall);for(auto&v:sum){wr(out,v.real());wr(out,v.imag());}std::cerr<<"CPP_GL_ETA44_MERGE rank="<<rank<<" partials="<<(argc-2)<<" max_node_seconds="<<wall<<"\n";return 0;}catch(std::exception const&e){std::cerr<<"CPP_GL_ETA44_MERGE_ERROR "<<e.what()<<"\n";return 1;}}
