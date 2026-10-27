#ifndef PMBFA_TYPES_H
#define PMBFA_TYPES_H
// R headers can name diagnostic groups absent from the host Clang version.
// Localize this compatibility suppression to external headers only.
#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunknown-warning-option"
#endif
#include <Rcpp.h>
#if defined(__clang__)
#pragma clang diagnostic pop
#endif
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <limits>
#include <stdexcept>
#include <functional>
#include <chrono>
#include <sstream>
#include <iomanip>
namespace pm {
using namespace Rcpp;
using V = std::vector<double>;
constexpr double PI=3.1415926535897932384626433832795;
constexpr double INF=std::numeric_limits<double>::infinity();
inline bool finite(double x) {return std::isfinite(x);}
inline void need(bool b,const char* s) {if(!b) throw std::runtime_error(s);}
template<class T> T option(const List& a,const char* name,T fallback) {
  return a.containsElementNamed(name)&&!Rf_isNull(a[name]) ? as<T>(a[name]):fallback;
}
struct Mat {
  int n=0,p=0; V v;
  Mat()=default;
  Mat(int nr,int nc):n(nr),p(nc),v((size_t)nr*nc,0){}
  double& operator()(int i,int j){return v[i+(size_t)n*j];}
  double operator()(int i,int j)const{return v[i+(size_t)n*j];}
};
Mat transpose(const Mat& x);
Mat multiply(const Mat& a,const Mat& b);
double dot(const V& a,const V& b);
Mat identity(int n);
Mat chol(const Mat& a);
Mat forward(const Mat& l,const Mat& b);
Mat backward(const Mat& l,const Mat& b);
Mat solve(const Mat& a,const Mat& b);
NumericMatrix to_r(const Mat& m);
Mat from_r(SEXP x,int n,int p,const char* label);
struct Precision {
  V inv;Mat K,Kt,B,W;double logdet=0;
  Precision()=default;
  Precision(const V& d,const Mat& loading);
  Mat left(const Mat& a) const;
  Mat right(const Mat& a) const;
};
struct Controls {
  int maxit,minit,backtracks,pcgmax;double tol,mag,floor,numin,numax,shapemin,mincomp,pcgtol;
  bool explicit_tol,track_likelihood,aitken,profile;std::string pcg_preconditioner;
  explicit Controls(const List& a);
};
struct Component {Mat M,A,L,D;V sig,psi;double pi=1,nu=10,omega=50,lambda=2,gam=1,gamma=1;};
struct Params {int n,p,N,G,q,r;std::string family,name;List metadata;V labels;std::vector<Component> groups;};
void check_name(const std::string& s,bool normal);
Params parse(const NumericVector& x,const List& input,const Controls& c,bool normal);
void project(Params& p);
struct Moments {double a,b,c;};
double logbk(double x,double nu);
double dlogbk(double x,double nu);
Moments gig(double lam,double chi,double psi,double cached_logk=NA_REAL);
double log_integral(double lam,double chi,double psi,double* cached_logk=nullptr);
double shape_root(double previous,double target,double lower,double upper,double scale);
double maximize(const std::function<double(double)>& f,double old,double step,double lower,double upper);
struct State {Mat z,a,b,c;V ng,sa,sb,sc;double ll=-INF;};
double quadratic_distance(const Mat& y,const Precision& s,const Precision& t);
Mat residual(const NumericVector& x,const Component& a,int obs,int n,int p,bool trans=false);
State expectation(const NumericVector& x,const Params& p,bool moments);
Mat weighted_operator(const Mat& x,const std::vector<Precision>& s,const std::vector<Precision>& t,const V& w);
Mat pcg(const Mat& start,const Mat& rhs,const std::vector<Precision>& s,const std::vector<Precision>& t,const V& w,const Controls& c,int& iterations,const Params& parameters);
void update_location(const NumericVector& x,Params& p,const State& e,const Controls& c,int& pcg_iterations);
void update_family(Params& p,const State& e,const Controls& c,int& retained);
struct Sufficient {V diagonal;Mat cross,second;};
Sufficient sufficient(const NumericVector& x,const Params& p,const Component& a,const State& e,int g,bool columns);
void update_scale(const NumericVector& x,Params& p,const State& e,const Controls& c,bool columns);
double count(const Params& p);
List serialize(const Params& p);
Params blend(const Params& a,const Params& b,double step);
List fit(const NumericVector& x,const List& input,const List& options,bool normal);
} // namespace pm
namespace pmbfa_init {
Rcpp::List initialize(const Rcpp::NumericVector& x,const int G,const int q,const int r,
                     const std::string family,Rcpp::List aargs,const Rcpp::NumericVector labels);
}
#endif
