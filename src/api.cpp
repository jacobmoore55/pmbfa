#include "pmbfa_types.h"
using namespace Rcpp;

// The entire normal/skewed algorithm remains inside these single calls.
// [[Rcpp::export]]
List pmsbfaEM(const NumericVector& x,const List& params,const List& aargs){return pm::fit(x,params,aargs,false);}
// [[Rcpp::export]]
List pmbfaEM(const NumericVector& x,const List& params,const List& aargs){return pm::fit(x,params,aargs,true);}

// Low-level likelihood evaluation uses the same implementation as fitting.
// [[Rcpp::export]]
List pmbfa_evaluate(const NumericVector& x,const List& params,Nullable<List> aargs=R_NilValue){List a=aargs.isNotNull()?clone(as<List>(aargs)):List::create();a["max.iter"]=0;return pm::fit(x,params,a,false);}
// [[Rcpp::export]]
std::string pmbfa_core_version(){return "1.2.0";}

// Count a candidate before fitting it. Used by optional nested-model screening.
// [[Rcpp::export]]
double pmbfa_model_npar(int n,int p,int G,int q,int r,std::string model,std::string family="ST"){
  if(n<2||p<2||G<1||q<1||q>=n||r<1||r>=p)stop("Invalid model dimensions.");
  if(family!="Normal"&&family!="ST"&&family!="GH"&&family!="NIG"&&family!="VG")stop("Unknown family.");
  pm::check_name(model,family=="Normal");
  if((model[2]=='U'&&(n-q)*(n-q)<n+q)||(model[6]=='U'&&(p-r)*(p-r)<p+r))stop("Factor rank exceeds the identifiable diagonal factor-analysis range.");
  pm::Params metadata;metadata.n=n;metadata.p=p;metadata.G=G;metadata.q=q;metadata.r=r;metadata.name=model;metadata.family=family;return pm::count(metadata);
}

// Post-fit diagnostics, deliberately outside the hot fitting loop.
// [[Rcpp::export]]
List pmbfa_posterior_moments(const NumericVector& x,const List& params,Nullable<List> aargs=R_NilValue){
  List controls=aargs.isNotNull()?as<List>(aargs):List::create();pm::Controls settings(controls);
  pm::Params model=pm::parse(x,params,settings,false);pm::project(model);
  try{pm::State e=pm::expectation(x,model,true);
    return List::create(_["z"]=pm::to_r(e.z),_["a"]=pm::to_r(e.a),_["b"]=pm::to_r(e.b),_["c"]=pm::to_r(e.c),_["ng"]=wrap(e.ng),_["sa"]=wrap(e.sa),_["sb"]=wrap(e.sb),_["sc"]=wrap(e.sc),_["ll"]=e.ll,_["status"]="ok",_["message"]="");
  }catch(const Rcpp::internal::InterruptedException&){throw;}catch(const std::exception& ex){return List::create(_["status"]="numerical_error",_["message"]=ex.what());}
}

// [[Rcpp::export]]
Rcpp::List pmbfa_initializations(const Rcpp::NumericVector& x,const int G,const int q,const int r,
  const std::string family="ST",Rcpp::List aargs=R_NilValue,
  const Rcpp::NumericVector labels=Rcpp::NumericVector::create()) {
  return pmbfa_init::initialize(x,G,q,r,family,aargs,labels);
}
