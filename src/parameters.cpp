#include "pmbfa_types.h"

namespace pm {
Controls::Controls(const List& a){
    maxit=option(a,"max.iter",10000);minit=option(a,"min.iter",2);
    backtracks=option(a,"backtracks",12);pcgmax=option(a,"pcg.max.iter",0);
    explicit_tol=a.containsElementNamed("tol")&&!Rf_isNull(a["tol"]);
    tol=option(a,"tol",1e-6);mag=option(a,"mag",7.0);
    floor=option(a,"variance.floor",.1);numin=option(a,"nu.min",4.1);numax=option(a,"nu.max",1e12);
    shapemin=option(a,"shape.min",.5);mincomp=option(a,"min.component",6.0);pcgtol=option(a,"pcg.tol",1e-9);
    track_likelihood=option(a,"track_likelihood",false);profile=option(a,"profile",false);std::string s=option(a,"stopping",std::string("increment"));
    aitken=s=="aitken";
    pcg_preconditioner=option(a,"pcg.preconditioner",std::string("diagonal"));
    if(pcg_preconditioner!="component"&&pcg_preconditioner!="diagonal")stop("pcg.preconditioner must be 'component' or 'diagonal'.");
    if(s!="increment"&&s!="aitken")stop("stopping must be 'increment' or 'aitken'.");
    if(maxit<0||minit<1||backtracks<0||pcgmax<0||!(tol>0)||!finite(tol)||!finite(mag)||!(floor>0)||!finite(floor)||!(numin>0)||numax<numin||!finite(numax)||!(shapemin>0)||!finite(shapemin)||mincomp<0||!finite(mincomp)||!(pcgtol>0)||!finite(pcgtol))stop("Invalid AECM controls.");
  }
void check_name(const std::string& s,bool normal){
  if(s.size()!=(normal?7:11)||s[3]!='-'||(!normal&&(s[7]!='-'||s[9]!='-')))stop("Invalid constraint model string.");
  for(size_t i=0;i<s.size();++i)if(i!=3&&i!=7&&i!=9&&s[i]!='C'&&s[i]!='U')stop("Constraint flags must be C or U.");
}
Params parse(const NumericVector& x,const List& input,const Controls& c,bool normal){
  if(!input.containsElementNamed("moddat")||!input.containsElementNamed("gmp")||!input.containsElementNamed("args"))stop("params must contain moddat, gmp and args.");
  if(!x.hasAttribute("dim"))stop("x must be an n by p by N numeric array.");
  IntegerVector dim=x.attr("dim");if(dim.size()!=3||dim[0]<2||dim[1]<2||dim[2]<=10)stop("x must have matrix dimensions >=2 and N > 10.");
  Params o; o.n=dim[0];o.p=dim[1];o.N=dim[2];o.metadata=clone(as<List>(input["moddat"]));List md=o.metadata,gm=input["gmp"],args=input["args"];
  o.G=option(md,"G",0);o.q=option(md,"q",0);o.r=option(md,"r",0);
  o.family=option(md,"family",std::string(normal?"Normal":"ST"));o.name=option(md,"name",std::string(normal?"UUU-UUU":"UUU-UUU-U-U"));
  if(normal&&o.family!="Normal")stop("pmbfaEM requires family='Normal'.");
  if(o.family!="Normal"&&o.family!="ST"&&o.family!="GH"&&o.family!="NIG"&&o.family!="VG")stop("Unknown family.");
  check_name(o.name,o.family=="Normal");
  if(o.G<1||o.G>o.N||o.q<1||o.q>=o.n||o.r<1||o.r>=o.p||args.size()!=o.G)stop("Invalid G, q, r or component parameter count.");
  if((o.name[2]=='U'&&(o.n-o.q)*(o.n-o.q)<o.n+o.q)||
     (o.name[6]=='U'&&(o.p-o.r)*(o.p-o.r)<o.p+o.r))
    stop("Factor rank exceeds the identifiable diagonal factor-analysis range.");
  for(double z:x)if(!finite(z))stop("x must contain only finite values.");
  for(const char* s:{"Sigmak","Psik","Lambdak","Deltak","pik"})if(!gm.containsElementNamed(s))stop("gmp is missing %s.",s);
  Mat sig=from_r(gm["Sigmak"],o.n,o.G,"Sigmak"),psi=from_r(gm["Psik"],o.p,o.G,"Psik");
  NumericVector lv=as<NumericVector>(gm["Lambdak"]),dv=as<NumericVector>(gm["Deltak"]),pv=as<NumericVector>(gm["pik"]);
  auto array_dims=[](const NumericVector& v,int g,int n,int k){if(!v.hasAttribute("dim"))return false;IntegerVector d=v.attr("dim");return d.size()==3&&d[0]==g&&d[1]==n&&d[2]==k;};
  if(!array_dims(lv,o.G,o.n,o.q)||!array_dims(dv,o.G,o.p,o.r)||pv.size()!=o.G)stop("Incorrect loading array or mixing weight dimensions.");
  o.labels.assign(o.N,0);
  if(md.containsElementNamed("class")&&!Rf_isNull(md["class"])){
    NumericVector labels=as<NumericVector>(md["class"]);if(labels.size()!=o.N)stop("class must have length N.");
    for(int i=0;i<o.N;++i){double z=labels[i];if(!finite(z)||z<0||z>o.G||z!=std::floor(z))stop("class must use 0 for unknown labels and 1..G for known labels.");o.labels[i]=z;}
  }
  o.groups.resize(o.G);double sumpi=0;
  for(int g=0;g<o.G;++g){Component& a=o.groups[g];List ag=args[g];
    if(!ag.containsElementNamed("M"))stop("Each component must contain M.");a.M=from_r(ag["M"],o.n,o.p,"M");
    a.A=o.family=="Normal"?Mat(o.n,o.p):from_r(ag["A"],o.n,o.p,"A");a.L=Mat(o.n,o.q);a.D=Mat(o.p,o.r);
    a.sig.resize(o.n);a.psi.resize(o.p);
    for(int i=0;i<o.n;++i){a.sig[i]=std::max(c.floor,sig(i,g));for(int u=0;u<o.q;++u)a.L(i,u)=lv[g+o.G*(i+o.n*u)];}
    for(int j=0;j<o.p;++j){a.psi[j]=std::max(c.floor,psi(j,g));for(int u=0;u<o.r;++u)a.D(j,u)=dv[g+o.G*(j+o.p*u)];}
    for(double v:a.L.v)if(!finite(v))stop("Loadings must be finite.");for(double v:a.D.v)if(!finite(v))stop("Loadings must be finite.");
    a.pi=pv[g];if(!(a.pi>0)||!finite(a.pi))stop("Mixing weights must be positive and finite.");sumpi+=a.pi;
    a.nu=option(ag,"nu",10.0);a.omega=option(ag,"omega",50.0);a.lambda=option(ag,"lambda",2.0);a.gam=option(ag,"gam",1.0);a.gamma=option(ag,"gamma",1.0);
    for(double z:{a.nu,a.omega,a.lambda,a.gam,a.gamma})if(!finite(z))stop("Family parameters must be finite.");
    a.nu=std::min(c.numax,std::max(c.numin,a.nu));a.gam=std::max(c.shapemin,a.gam);a.gamma=std::max(c.shapemin,a.gamma);
    if(o.family=="GH"&&!(a.omega>0))stop("GH omega must be positive.");
  }for(auto& a:o.groups)a.pi/=sumpi;return o;
}
// Project each fresh shared start into the candidate parameter space.
void project(Params& p){
  auto matrices=[&](int flag,Mat Component::*member){if(p.name[flag]!='C')return;Mat mean=p.groups[0].*member;
    std::fill(mean.v.begin(),mean.v.end(),0);for(const auto& a:p.groups)for(size_t i=0;i<mean.v.size();++i)mean.v[i]+=(a.*member).v[i]/p.G;
    for(auto& a:p.groups)a.*member=mean;};
  matrices(0,&Component::L);matrices(4,&Component::D);
  auto diagonal=[&](int start,V Component::*member){
    if(p.name[start+1]=='C'){V mean=(p.groups[0].*member);std::fill(mean.begin(),mean.end(),0);
      for(const auto& a:p.groups)for(size_t i=0;i<mean.size();++i)mean[i]+=(a.*member)[i]/p.G;
      for(auto& a:p.groups)a.*member=mean;}
    if(p.name[start+2]=='C')for(auto& a:p.groups){V& d=a.*member;double s=std::accumulate(d.begin(),d.end(),0.0)/d.size();std::fill(d.begin(),d.end(),s);}
  };diagonal(0,&Component::sig);diagonal(4,&Component::psi);
  if(p.family!="Normal"){
    matrices(8,&Component::A);
    if(p.name[10]=='C'){
      double Component::*key=p.family=="ST"?&Component::nu:p.family=="GH"?&Component::omega:p.family=="NIG"?&Component::gam:&Component::gamma;
      double s=0;for(const auto& a:p.groups)s+=a.*key/p.G;for(auto& a:p.groups)a.*key=s;
    }
  }
}

double count(const Params& p){
  auto block=[&](int d,int k,int offset){double loading=d*(double)k-k*(k-1)/2.0;
    return loading*(p.name[offset]=='C'?1:p.G)+(p.name[offset+2]=='C'?1:d)*(p.name[offset+1]=='C'?1:p.G);};
  double shapes=0,skew=0;
  if(p.family!="Normal"){skew=p.n*(double)p.p*(p.name[8]=='C'?1:p.G);shapes=p.family=="GH"?(p.name[10]=='C'?p.G+1:2*p.G):(p.name[10]=='C'?1:p.G);}
  int scale=p.name[0]=='U'&&p.name[1]=='U'&&p.name[4]=='U'&&p.name[5]=='U'?p.G:1;
  return p.G-1+p.G*(double)p.n*p.p+skew+block(p.n,p.q,0)+block(p.p,p.r,4)+shapes-scale;
}
List serialize(const Params& p){
  NumericMatrix sig(p.n,p.G),psi(p.p,p.G);NumericVector l((R_xlen_t)p.G*p.n*p.q),d((R_xlen_t)p.G*p.p*p.r),pi(p.G);List args(p.G);
  l.attr("dim")=IntegerVector::create(p.G,p.n,p.q);d.attr("dim")=IntegerVector::create(p.G,p.p,p.r);
  for(int g=0;g<p.G;++g){const Component& a=p.groups[g];pi[g]=a.pi;
    for(int i=0;i<p.n;++i){sig(i,g)=a.sig[i];for(int u=0;u<p.q;++u)l[g+p.G*(i+p.n*u)]=a.L(i,u);}
    for(int j=0;j<p.p;++j){psi(j,g)=a.psi[j];for(int u=0;u<p.r;++u)d[g+p.G*(j+p.p*u)]=a.D(j,u);}
    List ag=List::create(_["M"]=to_r(a.M));if(p.family!="Normal")ag["A"]=to_r(a.A);
    if(p.family=="ST")ag["nu"]=a.nu;if(p.family=="GH"){ag["omega"]=a.omega;ag["lambda"]=a.lambda;}if(p.family=="NIG")ag["gam"]=a.gam;if(p.family=="VG")ag["gamma"]=a.gamma;args[g]=ag;
  }
  List md=clone(p.metadata);md["n"]=p.n;md["p"]=p.p;md["N"]=p.N;md["G"]=p.G;md["q"]=p.q;md["r"]=p.r;md["family"]=p.family;md["name"]=p.name;md["class"]=wrap(p.labels);
  return List::create(_["moddat"]=md,_["gmp"]=List::create(_["Sigmak"]=sig,_["Psik"]=psi,_["Lambdak"]=l,_["Deltak"]=d,_["pik"]=pi),_["args"]=args);
}
Params blend(const Params& a,const Params& b,double step){
  Params o=a;
  for(int g=0;g<a.G;++g){Component& z=o.groups[g];const Component& x=a.groups[g];const Component& y=b.groups[g];
    for(Mat Component::*key:{&Component::M,&Component::A,&Component::L,&Component::D})for(size_t i=0;i<(z.*key).v.size();++i)(z.*key).v[i]=(x.*key).v[i]+step*((y.*key).v[i]-(x.*key).v[i]);
    for(V Component::*key:{&Component::sig,&Component::psi})for(size_t i=0;i<(z.*key).size();++i)(z.*key)[i]=(x.*key)[i]+step*((y.*key)[i]-(x.*key)[i]);
    for(double Component::*key:{&Component::pi,&Component::nu,&Component::omega,&Component::lambda,&Component::gam,&Component::gamma})z.*key=x.*key+step*(y.*key-x.*key);
  }return o;
}

} // namespace pm
