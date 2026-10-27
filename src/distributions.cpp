#include "pmbfa_types.h"

namespace pm {
// Scaled log K: the leading uniform asymptotic plus first Debye correction
// is only a fallback for orders outside the finite range of R's implementation.
double logbk(double x,double nu){
  need(x>0&&finite(x)&&finite(nu),"Invalid Bessel K arguments.");double v=std::abs(nu);
  long double z=(v>0)?(long double)x/v:1,t=std::hypot(1.0L,z),ap=0;
  if(v>=10){long double eta=t+std::log(z)-std::log1p(t),u=1/t;
    ap=(long double)x-v*eta+.5L*std::log((long double)PI/(2*v))-.5L*std::log(t);
    long double correction=(3*u-5*u*u*u)/(24*v);ap+=std::log1p(-correction);
  }
  if(v<10||ap<690){double b=Rf_bessel_k(x,v,2.0);if(b>0&&finite(b))return std::log(b);}
  if(v>=10&&finite((double)ap))return (double)ap;
  throw std::runtime_error("Bessel K evaluation failed.");
}
double dlogbk(double x,double nu){
  double h=1e-4*(1+std::abs(nu));
  return (logbk(x,nu-2*h)-8*logbk(x,nu-h)+8*logbk(x,nu+h)-logbk(x,nu+2*h))/(12*h);
}
Moments gig(double lam,double chi,double psi,double cached_logk){
  need(chi>=0&&psi>=0,"Invalid GIG parameters.");Moments m;
  if(psi==0){need(chi>0&&lam< -1,"Inverse gamma conditional moments are undefined.");m={chi/(-2*lam-2),-2*lam/chi,std::log(chi/2)-R::digamma(-lam)};}
  else if(chi==0){need(lam>1,"Gamma conditional inverse moment is undefined.");m={2*lam/psi,psi/(2*(lam-1)),R::digamma(lam)-std::log(psi/2)};}
  else {double arg=std::sqrt(chi)*std::sqrt(psi),base=finite(cached_logk)?cached_logk:logbk(arg,lam);
    m={std::sqrt(chi/psi)*std::exp(logbk(arg,lam+1)-base),std::sqrt(psi/chi)*std::exp(logbk(arg,lam-1)-base),.5*std::log(chi/psi)+dlogbk(arg,lam)};}
  need(m.a>0&&m.b>0&&finite(m.a)&&finite(m.b)&&finite(m.c),"Conditional mixing moments are nonfinite.");return m;
}
double log_integral(double lam,double chi,double psi,double* cached_logk){
  if(psi==0){need(chi>0&&lam<0,"Undefined mixing integral.");return R::lgammafn(-lam)+lam*std::log(chi/2);}
  if(chi==0){need(psi>0&&lam>0,"Unbounded variance-gamma density at an observation.");return R::lgammafn(lam)-lam*std::log(psi/2);}
  double arg=std::sqrt(chi)*std::sqrt(psi),base=logbk(arg,lam);
  if(cached_logk)*cached_logk=base;
  return std::log(2.0)+.5*lam*std::log(chi/psi)+base-arg;
}
// Adaptively bracketed, safeguarded Newton/bisection as in mshmm.cpp.
double shape_root(double previous,double target,double lower,double upper,double scale){
  auto f=[&](double x){return std::log(x/scale)+1-R::digamma(x/scale)-target;};
  double lo=lower,fl=f(lo);if(!finite(fl))return previous;if(fl<=0)return lo;
  double hi=std::min(upper,std::max(lo*2,previous)),fh=f(hi);
  while(finite(fh)&&fh>0&&hi<upper){hi=std::min(upper,2*hi);fh=f(hi);}
  if(!finite(fh)||fh>0)return previous;
  double v=std::max(lo,std::min(hi,previous));
  for(int i=0;i<100;++i){double fv=f(v);if(!finite(fv))return previous;
    if(std::abs(fv)<1e-11)return v;if(fv>0)lo=v;else hi=v;
    if(hi-lo<1e-10*(1+v))return (hi+lo)/2;
    double df=1/v-R::trigamma(v/scale)/scale,next=v-fv/df;
    v=finite(next)&&next>lo&&next<hi?next:(lo+hi)/2;
  }return (lo+hi)/2;
}
// Bracket and maximize a one-dimensional conditional objective; retain the
// previous value if the search cannot produce an improvement.
double maximize(const std::function<double(double)>& f,double old,double step,double lower,double upper){
  double a=std::max(lower,old-step),b=std::min(upper,old+step),fo=f(old);
  for(int j=0;j<30;++j){double fa=f(a),fb=f(b);
    if(fa<=fo&&fb<=fo)break;
    if(fa>fo&&fa>=fb){b=old;old=a;fo=fa;a=std::max(lower,a-2*step);}
    else{a=old;old=b;fo=fb;b=std::min(upper,b+2*step);}step*=2;
  }
  constexpr double ratio=.6180339887498948482;double l=b-ratio*(b-a),r=a+ratio*(b-a),fl=f(l),fr=f(r);
  for(int i=0;i<70&&b-a>1e-8*(1+std::abs(old));++i){
    if(fl>fr){b=r;r=l;fr=fl;l=b-ratio*(b-a);fl=f(l);}else{a=l;l=r;fl=fr;r=a+ratio*(b-a);fr=f(r);}
  }double v=fl>fr?l:r;return finite(f(v))&&f(v)>fo?v:old;
}
void update_family(Params& p,const State& e,const Controls& c,int& retained){
  if(p.family=="Normal")return;bool common=p.name[10]=='C';
  if(p.family=="GH"){
    // Conditional coordinate maximizations replace the unstable multiplicative
    // lambda update. Common concentration means common omega only.
    for(int g=0;g<p.G;++g){Component& a=p.groups[g];double cb=e.sc[g]/e.ng[g],om=a.omega,old=a.lambda;
      a.lambda=maximize([&](double l){return l*cb-logbk(om,l);},old,1,-1e6,1e6);if(a.lambda==old)++retained;}
    for(int g=0;g<(common?1:p.G);++g){double old=std::log(p.groups[g].omega);
      auto obj=[&](double logomega){double om=std::exp(logomega),v=0;for(int h=common?0:g;h<(common?p.G:g+1);++h)v+=e.ng[h]*(om-logbk(om,p.groups[h].lambda))-.5*om*(e.sa[h]+e.sb[h]);return v;};
      double value=std::exp(maximize(obj,old,.5,std::log(1e-8),std::log(1e8)));
      if(value==p.groups[g].omega)++retained;for(int h=common?0:g;h<(common?p.G:g+1);++h)p.groups[h].omega=value;
    }return;
  }
  for(int g=0;g<(common?1:p.G);++g){double ng=0,sa=0,sb=0,sc=0;
    for(int h=common?0:g;h<(common?p.G:g+1);++h){ng+=e.ng[h];sa+=e.sa[h];sb+=e.sb[h];sc+=e.sc[h];}
    Component& a=p.groups[g];double value=0,old=0;
    if(p.family=="ST"){old=a.nu;value=shape_root(old,(sb+sc)/ng,c.numin,c.numax,2);}
    else if(p.family=="VG"){old=a.gamma;value=shape_root(old,(sa-sc)/ng,c.shapemin,1e12,1);}
    else{old=a.gam;value=std::max(c.shapemin,ng/sa);}
    if(value==old)++retained;
    for(int h=common?0:g;h<(common?p.G:g+1);++h){if(p.family=="ST")p.groups[h].nu=value;else if(p.family=="VG")p.groups[h].gamma=value;else p.groups[h].gam=value;}
  }
}

} // namespace pm
