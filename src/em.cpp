#include "pmbfa_types.h"

namespace pm {
List fit(const NumericVector& x,const List& input,const List& options,bool normal){
  Controls control(options);
  using Clock=std::chrono::steady_clock;auto began=Clock::now();V seconds(7,0);
  struct Timer{bool enabled;double& elapsed;Clock::time_point start;
    Timer(bool on,double& value):enabled(on),elapsed(value){if(enabled)start=Clock::now();}
    ~Timer(){if(enabled)elapsed+=std::chrono::duration<double>(Clock::now()-start).count();}};
  Params accepted=parse(x,input,control,normal);project(accepted);
  // The default is calibrated once. An explicit tol is never overwritten.
  if(!options.containsElementNamed("mag")||Rf_isNull(options["mag"]))control.mag=std::log10(std::sqrt((double)accepted.N)*accepted.n*accepted.p)+3;
  State state;V curve,pcg_curve,rejected_sizes;int iterations=0,pcg_iterations=0,retained=0,reductions=0;bool converged=false;std::string status="initialization_failed",message="";
  double tolerance=control.tol;
  try{
    {Timer timer(control.profile,seconds[0]);state=expectation(x,accepted,control.maxit>0);}curve.push_back(state.ll);
    if(!control.explicit_tol)tolerance=std::max(1.0,std::abs(state.ll))*std::pow(10.0,-control.mag);
    status=control.maxit==0?"evaluated":"max_iter";
    auto valid_sizes=[&](const State& candidate){
      return std::all_of(candidate.ng.begin(),candidate.ng.end(),[&](double ng){return ng>control.mincomp;});
    };
    for(int iter=1;iter<=control.maxit;++iter){checkUserInterrupt();
      if(!valid_sizes(state)){rejected_sizes=state.ng;status="insufficient_component_size";message="Effective component sizes must be strictly greater than min.component.";break;}
      Params trial=accepted;int pcg_before=pcg_iterations;
      {Timer timer(control.profile,seconds[1]);update_location(x,trial,state,control,pcg_iterations);}
      {Timer timer(control.profile,seconds[2]);update_family(trial,state,control,retained);}
      {Timer timer(control.profile,seconds[3]);update_scale(x,trial,state,control,false);}
      {Timer timer(control.profile,seconds[4]);update_scale(x,trial,state,control,true);}
      State next;{Timer timer(control.profile,seconds[5]);next=expectation(x,trial,iter<control.maxit);}
      if(!valid_sizes(next)){rejected_sizes=next.ng;status="insufficient_component_size";message="Proposed effective component size is not strictly greater than min.component; retained last accepted state.";break;}
      if(next.ll<state.ll){bool found=false;
        for(int b=1;b<=control.backtracks;++b){Timer timer(control.profile,seconds[6]);Params damped=blend(accepted,trial,std::ldexp(1.0,-b));State proposed=expectation(x,damped,iter<control.maxit);++reductions;
          if(proposed.ll>=state.ll&&valid_sizes(proposed)){trial=std::move(damped);next=std::move(proposed);found=true;break;}}
        if(!found){double drop=state.ll-next.ll;
          if(drop<=std::min(tolerance,1e-8*(1+std::abs(state.ll)))&&iterations>=control.minit){converged=true;status="converged";message="Stationary within floating-point resolution; retained last accepted state.";}
          else{status="nonmonotone";message="A proposed AECM cycle decreased likelihood; retained last accepted state.";}break;
        }
      }
      double gain=next.ll-state.ll;accepted=std::move(trial);state=std::move(next);curve.push_back(state.ll);pcg_curve.push_back(pcg_iterations-pcg_before);++iterations;
      if(control.verbose){std::ostringstream line;line<<std::fixed<<std::setprecision(2)<<"Iteration "<<iterations<<": log likelihood = "<<state.ll<<", change = "<<gain<<"\n";Rcout<<line.str();}
      bool stop_now=false;
      if(iterations>=control.minit){
        if(!control.aitken)stop_now=gain<=tolerance;
        else if(curve.size()>=3){double previous=curve[curve.size()-2]-curve[curve.size()-3];
          if(gain==0)stop_now=true;else if(previous>0){double acceleration=gain/previous;
            if(acceleration>=0&&acceleration<1){double gap=gain*acceleration/(1-acceleration);stop_now=finite(gap)&&gap>=0&&gap<=tolerance;}}}
      }
      if(stop_now){converged=true;status="converged";break;}
    }
  }catch(const Rcpp::internal::InterruptedException&){throw;}
  catch(const std::exception& ex){status=curve.empty()?"initialization_failed":"numerical_error";message=ex.what();}
  double k=count(accepted),bic=curve.empty()?INF:-2*curve.back()+k*std::log((double)accepted.N),icl=bic;
  IntegerVector cluster(accepted.N,NA_INTEGER);NumericMatrix posterior(accepted.N,accepted.G);std::fill(posterior.begin(),posterior.end(),NA_REAL);
  if(!curve.empty())for(int i=0;i<accepted.N;++i){double best=-1;int group=0;
    for(int g=0;g<accepted.G;++g){posterior(i,g)=state.z(i,g);if(state.z(i,g)>best){best=state.z(i,g);group=g;}}
    cluster[i]=group+1;if(best>0)icl-=2*std::log(best);
  }
  std::string criterion=option(options,"criterion",option(accepted.metadata,"icl",false)?std::string("ICL"):std::string("BIC"));
  List out=List::create(_["params"]=serialize(accepted),_["loglik"]=wrap(curve),_["iter"]=iterations,_["converged"]=converged,_["status"]=status,_["message"]=message,
    _["posterior"]=posterior,_["cluster"]=cluster,_["npar"]=k,_["bic"]=bic,_["icl"]=icl,_["BIC"]=bic,_["ICL"]=icl,_["criterion"]=criterion,_["score"]=converged?(criterion=="ICL"?icl:bic):INF,
    _["control"]=clone(options),_["diagnostics"]=List::create(_["tolerance"]=tolerance,_["stopping"]=control.aitken?"aitken":"increment",_["pcg_iterations"]=pcg_iterations,_["shape_updates_retained"]=retained,_["backtracking_evaluations"]=reductions,_["effective_sizes"]=wrap(state.ng),_["rejected_effective_sizes"]=wrap(rejected_sizes)));
  if(control.profile){NumericVector timing=wrap(seconds);timing.names()=CharacterVector::create("initial_expectation","location","family","row","column","expectation","backtracking");
    List diagnostic=out["diagnostics"];diagnostic["timing"]=timing;diagnostic["pcg_iterations_by_cycle"]=wrap(pcg_curve);diagnostic["total_native_seconds"]=std::chrono::duration<double>(Clock::now()-began).count();out["diagnostics"]=diagnostic;}
  out.attr("class")=accepted.family=="Normal"?CharacterVector::create("pmbfa"):CharacterVector::create("pmsbfa","pmbfa");return out;
}

} // namespace pm
