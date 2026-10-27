#include "pmbfa_types.h"

namespace pm {
// Stream diag(C), C*B' and B*C*B'. The positive decomposition
// H=sqrt(b)*(X-M)-A/sqrt(b) avoids cancellation of four large scatters.
Sufficient sufficient(const NumericVector& x,const Params& p,const Component& a,const State& e,int g,bool columns){
  const Mat& loading=columns?a.D:a.L;const V& diagonal=columns?a.psi:a.sig;
  Precision own(diagonal,loading),other(columns?a.sig:a.psi,columns?a.L:a.D);
  int d=loading.n,k=loading.p,otherdim=columns?p.n:p.p,otherk=other.K.n;
  Sufficient out{V(d,0),Mat(d,k),Mat(k,k)};double kappa=0;
  Mat A=columns?transpose(a.A):a.A,Bt=transpose(own.B);
  bool normal=p.family=="Normal";
  // Reuse these buffers for every observation; never form the full residual
  // matrix or its product with the opposite precision.
  V h(d),bh(k),diagonal_part(d);Mat hk(d,otherk);
  auto accumulate=[&](int obs,double weight,double sqrtb,double inverse_sqrtb,bool only_A){
    std::fill(hk.v.begin(),hk.v.end(),0);std::fill(diagonal_part.begin(),diagonal_part.end(),0);
    size_t offset=(size_t)obs*p.n*p.p;
    for(int j=0;j<otherdim;++j){
      if(only_A)for(int i=0;i<d;++i)h[i]=A(i,j);
      else if(normal)for(int i=0;i<d;++i)h[i]=x[offset+(columns?j+(size_t)p.n*i:i+(size_t)p.n*j)]-(columns?a.M(j,i):a.M(i,j));
      else for(int i=0;i<d;++i)h[i]=sqrtb*(x[offset+(columns?j+(size_t)p.n*i:i+(size_t)p.n*j)]-(columns?a.M(j,i):a.M(i,j)))-inverse_sqrtb*A(i,j);
      for(int i=0;i<d;++i)diagonal_part[i]+=other.inv[j]*h[i]*h[i];
      for(int u=0;u<k;++u){
        // Four independent reduction chains permit vectorization without
        // changing the compiler's floating-point flags.
        double s0=0,s1=0,s2=0,s3=0;int i=0;
        for(;i+3<d;i+=4){s0+=Bt(i,u)*h[i];s1+=Bt(i+1,u)*h[i+1];s2+=Bt(i+2,u)*h[i+2];s3+=Bt(i+3,u)*h[i+3];}
        double sum=(s0+s1)+(s2+s3);for(;i<d;++i)sum+=Bt(i,u)*h[i];bh[u]=sum;
      }
      for(int u=0;u<k;++u){double v=weight*other.inv[j]*bh[u];for(int i=0;i<d;++i)out.cross(i,u)+=v*h[i];}
      for(int u=0;u<otherk;++u){double v=other.K(u,j);for(int i=0;i<d;++i)hk(i,u)+=v*h[i];}
    }
    Mat bhk=multiply(own.B,hk);
    for(int i=0;i<d;++i){long double v=diagonal_part[i];for(int u=0;u<otherk;++u)v-=(long double)hk(i,u)*hk(i,u);out.diagonal[i]+=weight*(double)v;}
    for(int u=0;u<k;++u)for(int v=0;v<otherk;++v){double z=weight*bhk(u,v);for(int i=0;i<d;++i)out.cross(i,u)-=z*hk(i,v);}
  };
  for(int i=0;i<p.N;++i){if((i&255)==0)checkUserInterrupt();double z=e.z(i,g);if(z==0)continue;
    double sqrtb=1,inverse_sqrtb=0;
    if(p.family!="Normal"){
      double b=e.b(i,g);sqrtb=std::sqrt(b);inverse_sqrtb=1/sqrtb;
      double remainder=e.a(i,g)-1/b;
      need(remainder>= -1e-7*e.a(i,g),"Conditional GIG moments violate the moment inequality.");kappa+=z*std::max(0.0,remainder);
    }accumulate(i,z,sqrtb,inverse_sqrtb,false);
  }
  if(kappa>0)accumulate(0,kappa,0,0,true);
  out.second=multiply(own.B,out.cross);
  for(int j=0;j<k;++j)for(int i=0;i<k;++i)out.second(i,j)+=e.ng[g]*otherdim*own.W(i,j);
  for(int j=0;j<k;++j)for(int i=0;i<j;++i)out.second(i,j)=out.second(j,i)=.5*(out.second(i,j)+out.second(j,i));
  return out;
}

void update_scale(const NumericVector& x,Params& p,const State& e,const Controls& c,bool columns){
  int start=columns?4:0,d=columns?p.p:p.n,k=columns?p.r:p.q,otherdim=columns?p.n:p.p;
  bool shared_loading=p.name[start]=='C',shared_diagonal=p.name[start+1]=='C',isotropic=p.name[start+2]=='C';
  std::vector<Sufficient> stat;stat.reserve(p.G);
  for(int g=0;g<p.G;++g)stat.push_back(sufficient(x,p,p.groups[g],e,g,columns));
  auto load=[&](int g)->Mat&{return columns?p.groups[g].D:p.groups[g].L;};
  auto diag=[&](int g)->V&{return columns?p.groups[g].psi:p.groups[g].sig;};
  if(!shared_loading){for(int g=0;g<p.G;++g)load(g)=transpose(solve(stat[g].second,transpose(stat[g].cross)));}
  else{
    Mat common(d,k);
    // One solve suffices when diagonal variances share their row weights.
    if(shared_diagonal||isotropic){Mat q(k,k),t(d,k);
      for(int g=0;g<p.G;++g){double w=shared_diagonal?1:1/diag(g)[0];for(size_t j=0;j<q.v.size();++j)q.v[j]+=w*stat[g].second.v[j];for(size_t j=0;j<t.v.size();++j)t.v[j]+=w*stat[g].cross.v[j];}
      common=transpose(solve(q,transpose(t)));
    }else for(int i=0;i<d;++i){Mat q(k,k),t(k,1);
      for(int g=0;g<p.G;++g){double w=1/diag(g)[i];for(size_t j=0;j<q.v.size();++j)q.v[j]+=w*stat[g].second.v[j];for(int u=0;u<k;++u)t(u,0)+=w*stat[g].cross(i,u);}
      Mat solution=solve(q,t);for(int u=0;u<k;++u)common(i,u)=solution(u,0);
    }for(int g=0;g<p.G;++g)load(g)=common;
  }
  std::vector<V> residuals(p.G,V(d));
  for(int g=0;g<p.G;++g){const Mat& l=load(g);Mat lq=multiply(l,stat[g].second);
    for(int i=0;i<d;++i){double v=stat[g].diagonal[i];for(int u=0;u<k;++u)v+=l(i,u)*(lq(i,u)-2*stat[g].cross(i,u));
      need(finite(v)&&v>= -1e-6*std::max(1.0,stat[g].diagonal[i]),"Negative residual variance numerator.");residuals[g][i]=std::max(0.0,v);}
  }
  if(shared_diagonal){V values(d,0);for(int g=0;g<p.G;++g)for(int i=0;i<d;++i)values[i]+=residuals[g][i]/(p.N*(double)otherdim);
    if(isotropic){double v=std::accumulate(values.begin(),values.end(),0.0)/d;std::fill(values.begin(),values.end(),v);}
    for(double& v:values)v=std::max(c.floor,v);for(int g=0;g<p.G;++g)diag(g)=values;
  }else for(int g=0;g<p.G;++g){V values=residuals[g];for(double& v:values)v/=e.ng[g]*otherdim;
    if(isotropic){double v=std::accumulate(values.begin(),values.end(),0.0)/d;std::fill(values.begin(),values.end(),v);}
    for(double& v:values)v=std::max(c.floor,v);diag(g)=values;
  }
}

} // namespace pm
