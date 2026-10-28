#include "pmbfa_types.h"

namespace pm {
// Scalar Woodbury expansion avoids constructing S^-1*Y*P^-1 just for its
// Frobenius inner product. Accumulate its four terms in extended precision.
double quadratic_distance(const Mat& y,const Precision& s,const Precision& t){
  Mat sy=multiply(s.K,y),yt=multiply(y,t.Kt);
  Mat syt=y.n<=y.p?multiply(s.K,yt):multiply(sy,t.Kt);
  long double base=0,row=0,column=0,correction=0;
  for(int j=0;j<y.p;++j)for(int i=0;i<y.n;++i)base+=(long double)s.inv[i]*t.inv[j]*y(i,j)*y(i,j);
  for(int j=0;j<sy.p;++j)for(int u=0;u<sy.n;++u)row+=(long double)t.inv[j]*sy(u,j)*sy(u,j);
  for(int u=0;u<yt.p;++u)for(int i=0;i<yt.n;++i)column+=(long double)s.inv[i]*yt(i,u)*yt(i,u);
  for(double z:syt.v)correction+=(long double)z*z;
  double delta=(double)(base-row-column+correction);
  if(delta<=1e-10*(double)base){Mat full=t.right(s.left(y));delta=dot(y.v,full.v);}
  return delta;
}
Mat residual(const NumericVector& x,const Component& a,int obs,int n,int p,bool trans){
  Mat y(trans?p:n,trans?n:p);size_t offset=(size_t)obs*n*p;
  for(int j=0;j<p;++j)for(int i=0;i<n;++i){double v=x[offset+i+(size_t)n*j]-a.M(i,j);if(trans)y(j,i)=v;else y(i,j)=v;}return y;
}
State expectation(const NumericVector& x,const Params& p,bool moments,BesselWorkspace& workspace){
  State e;e.z=Mat(p.N,p.G);e.a=Mat(p.N,p.G);e.b=Mat(p.N,p.G);e.c=Mat(p.N,p.G);
  e.ng.assign(p.G,0);e.sa.assign(p.G,0);e.sb.assign(p.G,0);e.sc.assign(p.G,0);
  for(int g=0;g<p.G;++g){const Component& a=p.groups[g];Precision s(a.sig,a.L),t(a.psi,a.D);
    double base=-.5*p.n*p.p*std::log(2*PI)-.5*p.p*s.logdet-.5*p.n*t.logdet+std::log(a.pi);
    double lam=0,chi0=0,psi0=0,prior=0,rho=0;Mat norm;
    if(p.family!="Normal"){
      norm=t.right(s.left(a.A));rho=dot(a.A.v,norm.v);need(rho>= -1e-8&&finite(rho),"Invalid skewness quadratic form.");rho=std::max(0.0,rho);
      if(p.family=="ST"){lam=-a.nu/2;chi0=a.nu;prior=.5*a.nu*std::log(a.nu/2)-R::lgammafn(a.nu/2);}
      if(p.family=="GH"){lam=a.lambda;chi0=psi0=a.omega;prior=-std::log(2.0)-logbk(a.omega,a.lambda,workspace)+a.omega;}
      if(p.family=="NIG"){lam=-.5;chi0=1;psi0=a.gam*a.gam;prior=a.gam-.5*std::log(2*PI);}
      if(p.family=="VG"){lam=a.gamma;psi0=2*a.gamma;prior=a.gamma*std::log(a.gamma)-R::lgammafn(a.gamma);}
      lam-=.5*p.n*p.p;
    }
    for(int i=0;i<p.N;++i){if((i&255)==0)checkUserInterrupt();Mat y=residual(x,a,i,p.n,p.p);
      double delta=quadratic_distance(y,s,t);need(finite(delta)&&delta>= -1e-7,"Invalid Mahalanobis distance.");delta=std::max(0.0,delta);
      if(p.family=="Normal"){e.z(i,g)=base-.5*delta;e.a(i,g)=e.b(i,g)=1;}
      else{double chi=chi0+delta,psi=psi0+rho,cached_logk=NA_REAL;
        e.z(i,g)=base+dot(y.v,norm.v)+prior+log_integral(lam,chi,psi,workspace,&cached_logk);
        if(moments){Moments w=gig(lam,chi,psi,workspace,cached_logk);e.a(i,g)=w.a;e.b(i,g)=w.b;e.c(i,g)=w.c;}}
      need(finite(e.z(i,g)),"Observed component density is nonfinite.");
    }
  }
  e.ll=0;
  for(int i=0;i<p.N;++i){
    double mx=-INF;for(int g=0;g<p.G;++g)mx=std::max(mx,e.z(i,g));double sum=0;for(int g=0;g<p.G;++g)sum+=std::exp(e.z(i,g)-mx);
    double ll=mx+std::log(sum);int known=(int)p.labels[i]-1;
    e.ll+=known>=0?e.z(i,known):ll;
    for(int g=0;g<p.G;++g){e.z(i,g)=known>=0?(g==known?1:0):std::exp(e.z(i,g)-ll);
      e.ng[g]+=e.z(i,g);e.sa[g]+=e.z(i,g)*e.a(i,g);e.sb[g]+=e.z(i,g)*e.b(i,g);e.sc[g]+=e.z(i,g)*e.c(i,g);}
  }need(finite(e.ll),"Observed likelihood is nonfinite.");return e;
}


} // namespace pm
