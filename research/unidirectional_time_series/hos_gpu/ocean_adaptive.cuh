#pragma once
// SPDX-License-Identifier: GPL-3.0-or-later
// Flat-bed/no-breaking/no-relaxation branch of the HOS-Ocean outer controller.
#include "ocean_cash_karp.cuh"
#include <limits>
#include <algorithm>
namespace hos {
struct AdaptiveTrace {
    int attempt;double time_before,proposal,actual,error,time_after;bool accepted;double next;
};
template<class T> class OceanAdaptive {
    OceanCashKarp<T>& solver_;
    double time_=0,proposal_,stop_,output_,linear_limit_,tolerance_;
    int attempts_=0,accepted_=0,rejected_=0;
public:
    using Trace=std::function<void(const AdaptiveTrace&)>;
    OceanAdaptive(OceanCashKarp<T>& solver,double stop,double output,double initial,double linear_limit,double tolerance)
        :solver_(solver),proposal_(initial),stop_(stop),output_(output),linear_limit_(linear_limit),tolerance_(tolerance){
        if(stop<=0||output<=0||initial<=0||linear_limit<=0||tolerance<=0)throw std::invalid_argument("Invalid adaptive parameters");
        solver_.prime();
    }
    bool finished() const {return std::abs(time_-stop_)<=std::numeric_limits<double>::epsilon();}
    double time()const{return time_;}int accepted()const{return accepted_;}int rejected()const{return rejected_;}
    void advance_output(const Trace& trace={}){
        if(finished())return;
        double target=time_+output_;
        // Exact source uses epsilon(0.0), a default-REAL constant, at this boundary.
        if(target-stop_ > -double(std::numeric_limits<float>::epsilon()))target=stop_;
        while(time_<target){
            if(++attempts_>10000)throw std::runtime_error("Diagnostic attempt budget exceeded");
            double slope=double(solver_.current_slope());
            if(!std::isfinite(slope)||slope>10)throw std::runtime_error("HOS slope check failed");
            double before=time_,old_proposal=proposal_,actual=std::min(proposal_,target-time_);
            if(actual<=0||time_+actual==time_)throw std::runtime_error("Time-step underflow");
            double error=double(solver_.attempt(T(time_),T(actual)));
            if(std::isnan(error)||error<0)throw std::runtime_error("Invalid temporal error estimate");
            bool accept=error<=tolerance_;double correction;
            if(!accept){correction=std::pow(error/tolerance_,-.25);++rejected_;}
            else {
                solver_.accept();time_+=actual;++accepted_;
                correction=error>std::numeric_limits<double>::epsilon()?std::pow(error/tolerance_,-.2):4.;
            }
            correction=.98*std::max(.125,std::min(4.,correction));
            proposal_=std::min({proposal_*correction,output_,linear_limit_});
            if(trace)trace({attempts_,before,old_proposal,actual,error,time_,accept,proposal_});
        }
    }
};
}
