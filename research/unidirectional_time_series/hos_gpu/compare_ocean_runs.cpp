#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>
#include <stdexcept>
using Row=std::array<double,10>;struct Z{double x,y;};
std::vector<Row> csv(const std::filesystem::path& path){
    std::ifstream f(path);std::string line;std::getline(f,line);std::vector<Row> rows;
    while(std::getline(f,line)){if(line.empty())continue;std::replace(line.begin(),line.end(),',',' ');std::istringstream s(line);Row r{};
        for(auto& v:r){if(!(s>>v)||!std::isfinite(v))throw std::runtime_error("Bad diagnostic row");}rows.push_back(r);}
    if(rows.empty())throw std::runtime_error("Missing diagnostics");return rows;
}
int main(int argc,char** argv){try{
    if(argc!=5)throw std::runtime_error("Usage: compare_runs probes|fields reference candidate_directory tolerance");
    bool fields=std::string(argv[1])=="fields";std::filesystem::path ref=argv[2],cand=argv[3];double tol=std::stod(argv[4]);
    auto a=csv(fields?ref/"diagnostics.csv":ref),b=csv(cand/"diagnostics.csv");if(b.size()<a.size())throw std::runtime_error("Candidate record too short");bool all=true;
    for(size_t i=0;i<a.size();++i)if(std::abs(a[i][0]-b[i][0])>1e-7)throw std::runtime_error("Time mismatch");
    for(int j=1;j<=5;++j){double sum=0,norm=0,peak=0;
        for(size_t i=0;i<a.size();++i){double d=b[i][j]-a[i][j];sum+=d*d;norm+=a[i][j]*a[i][j];peak=std::max(peak,std::abs(d));}
        double rel=std::sqrt(sum/norm);bool pass=rel<=tol;all=all&&pass;
        std::cout<<std::setprecision(17)<<"{\"type\":\"probe\",\"probe\":"<<j<<",\"relative_l2\":"<<rel<<",\"max_abs_m\":"<<peak<<",\"end_time_s\":"<<a.back()[0]<<",\"pass\":"<<(pass?"true":"false")<<"}\n";
    }
    double drift_difference=0,candidate_drift=0,reference_drift=0;
    for(size_t i=0;i<a.size();++i){
        // Official output.f90 writes ABS(E-E0)/E0; GPU files retain signed drift.
        double delta=fields?std::abs(a[i][9]-b[i][9]):std::abs(std::abs(a[i][9])-std::abs(b[i][9]));
        drift_difference=std::max(drift_difference,delta);candidate_drift=std::max(candidate_drift,std::abs(b[i][9]));reference_drift=std::max(reference_drift,std::abs(a[i][9]));}
    // Absolute drift is reported, not assumed to be purely temporal error.
    bool energy_pass=drift_difference<=std::max(tol,1e-8);all=all&&energy_pass;
    std::cout<<"{\"type\":\"energy\",\"max_drift_difference\":"<<drift_difference<<",\"reference_max_relative_drift\":"<<reference_drift<<",\"candidate_max_relative_drift\":"<<candidate_drift<<",\"end_time_s\":"<<a.back()[0]<<",\"pass\":"<<(energy_pass?"true":"false")<<"}\n";
    if(fields)for(const auto& entry:std::filesystem::directory_iterator(ref)){
        auto name=entry.path().filename();if(name.string().rfind("state-",0)!=0||entry.path().extension()!=".bin")continue;
        std::ifstream x(entry.path(),std::ios::binary),y(cand/name,std::ios::binary);int dx[2],dy[2];double tx,ty;uint64_t nx,ny;
        x.read(reinterpret_cast<char*>(dx),8);y.read(reinterpret_cast<char*>(dy),8);x.read(reinterpret_cast<char*>(&tx),8);y.read(reinterpret_cast<char*>(&ty),8);
        x.read(reinterpret_cast<char*>(&nx),8);y.read(reinterpret_cast<char*>(&ny),8);
        if(!x||!y||nx!=ny||dx[0]!=dy[0]||dx[1]!=dy[1]||std::abs(tx-ty)>1e-7)throw std::runtime_error("Checkpoint mismatch");
        for(int field=0;field<2;++field){double err=0,norm=0,peak=0,maxref=0,phase=0;size_t dominant=0;
            for(size_t i=0;i<nx;++i){Z ar,br;x.read(reinterpret_cast<char*>(&ar),16);y.read(reinterpret_cast<char*>(&br),16);
                if(!x||!y||!std::isfinite(ar.x)||!std::isfinite(ar.y)||!std::isfinite(br.x)||!std::isfinite(br.y))throw std::runtime_error("Nonfinite checkpoint");
                double d=std::hypot(ar.x-br.x,ar.y-br.y),r=std::hypot(ar.x,ar.y);int kx=int(i%(dx[0]/2+1));double w=(kx==0||kx==dx[0]/2)?1.:.5;
                err+=w*d*d;norm+=w*r*r;peak=std::max(peak,d);
                if(i!=0&&r>maxref){maxref=r;dominant=i;phase=std::atan2(br.y*ar.x-br.x*ar.y,br.x*ar.x+br.y*ar.y);}
            }
            double rel=std::sqrt(err/norm);bool pass=rel<=tol&&std::abs(phase)<=1e-3;all=all&&pass;
            std::cout<<"{\"type\":\"checkpoint\",\"time_s\":"<<tx<<",\"field\":\""<<(field?"psi":"eta")<<"\",\"physical_relative_l2\":"<<rel<<",\"modal_max_abs\":"<<peak<<",\"dominant_mode_index\":"<<dominant<<",\"dominant_phase_error_rad\":"<<phase<<",\"pass\":"<<(pass?"true":"false")<<"}\n";
        }
    }
    return all?0:1;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
