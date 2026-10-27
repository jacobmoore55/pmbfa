#include "pmbfa_types.h"

namespace pm {
Mat weighted_operator(const Mat& x,const std::vector<Precision>& s,const std::vector<Precision>& t,const V& w){
  Mat out(x.n,x.p);
  for(size_t g=0;g<w.size();++g)if(w[g]>0){Mat v=t[g].right(s[g].left(x));for(size_t i=0;i<out.v.size();++i)out.v[i]+=w[g]*v.v[i];}return out;
}
Mat pcg_attempt(const Mat& start,const Mat& rhs,const std::vector<Precision>& s,const std::vector<Precision>& t,const V& w,const Controls& c,int& iterations,const Params& parameters,bool component){
  size_t reference=std::distance(w.begin(),std::max_element(w.begin(),w.end()));
  const Component& pre=parameters.groups[reference];
  const Mat lt=component?transpose(pre.L):Mat(),dt=component?transpose(pre.D):Mat();
  double total_weight=std::accumulate(w.begin(),w.end(),0.0);
  auto precondition=[&](const Mat& residual){
    Mat left=multiply(pre.L,multiply(lt,residual));
    for(int j=0;j<residual.p;++j)for(int i=0;i<residual.n;++i)left(i,j)+=pre.sig[i]*residual(i,j);
    Mat out=multiply(multiply(left,pre.D),dt);
    for(int j=0;j<residual.p;++j)for(int i=0;i<residual.n;++i)out(i,j)=(out(i,j)+pre.psi[j]*left(i,j))/total_weight;
    return out;
  };
  Mat x=start,diag(x.n,x.p),ax=weighted_operator(x,s,t,w),r=rhs,z(x.n,x.p);
  for(size_t g=0;g<w.size();++g){V sd=s[g].inv,pd=t[g].inv;
    for(int i=0;i<x.n;++i)for(int u=0;u<s[g].K.n;++u)sd[i]-=s[g].K(u,i)*s[g].K(u,i);
    for(int j=0;j<x.p;++j)for(int u=0;u<t[g].K.n;++u)pd[j]-=t[g].K(u,j)*t[g].K(u,j);
    for(int j=0;j<x.p;++j)for(int i=0;i<x.n;++i)diag(i,j)+=w[g]*sd[i]*pd[j];
  }
  for(size_t i=0;i<r.v.size();++i){r.v[i]-=ax.v[i];need(diag.v[i]>0&&finite(diag.v[i]),"Invalid common-skewness preconditioner.");}
  if(component)z=precondition(r);else for(size_t i=0;i<r.v.size();++i)z.v[i]=r.v[i]/diag.v[i];
  Mat direction=z;double rz=dot(r.v,z.v),r0=std::sqrt(dot(r.v,r.v)),target=c.pcgtol*std::max(1.0,std::sqrt(dot(rhs.v,rhs.v)));
  if(r0<=target)return x;
  int maxit=c.pcgmax>0?c.pcgmax:(int)x.v.size();
  for(int it=0;it<maxit;++it){if((it&31)==0)checkUserInterrupt();Mat ad=weighted_operator(direction,s,t,w);double denom=dot(direction.v,ad.v);
    need(denom>0&&finite(denom)&&rz>0,"Common-skewness conjugate gradient broke down.");double alpha=rz/denom;
    for(size_t i=0;i<r.v.size();++i){x.v[i]+=alpha*direction.v[i];r.v[i]-=alpha*ad.v[i];}++iterations;
    if(std::sqrt(dot(r.v,r.v))<=target){Mat truth=weighted_operator(x,s,t,w);for(size_t i=0;i<truth.v.size();++i)truth.v[i]=rhs.v[i]-truth.v[i];if(std::sqrt(dot(truth.v,truth.v))<=2*target)return x;}
    if(component)z=precondition(r);else for(size_t i=0;i<r.v.size();++i)z.v[i]=r.v[i]/diag.v[i];double next=dot(r.v,z.v),beta=next/rz;
    for(size_t i=0;i<r.v.size();++i)direction.v[i]=z.v[i]+beta*direction.v[i];rz=next;
  }throw std::runtime_error("Common-skewness conjugate gradient did not converge.");
}
// A representative covariance supplies an SPD inverse approximation without
// forming a Kronecker matrix. Keep the original diagonal solver as a fallback.
Mat pcg(const Mat& start,const Mat& rhs,const std::vector<Precision>& s,const std::vector<Precision>& t,const V& w,const Controls& c,int& iterations,const Params& parameters){
  if(c.pcg_preconditioner=="diagonal")return pcg_attempt(start,rhs,s,t,w,c,iterations,parameters,false);
  try{return pcg_attempt(start,rhs,s,t,w,c,iterations,parameters,true);}
  catch(const Rcpp::internal::InterruptedException&){throw;}
  catch(const std::exception&){return pcg_attempt(start,rhs,s,t,w,c,iterations,parameters,false);}
}
void update_location(const NumericVector& x,Params& p,const State& e,const Controls& c,int& pcg_iterations){
  std::vector<Mat> sum(p.G,Mat(p.n,p.p)),bsum=sum;size_t d=(size_t)p.n*p.p;
  for(int g=0;g<p.G;++g){need(e.ng[g]>0&&e.sb[g]>0,"Empty component during location update.");
    p.groups[g].pi=e.ng[g]/p.N;
    for(int i=0;i<p.N;++i){double z=e.z(i,g),zb=z*e.b(i,g);for(size_t j=0;j<d;++j){double v=x[(size_t)i*d+j];sum[g].v[j]+=z*v;bsum[g].v[j]+=zb*v;}}
  }
  if(p.family=="Normal"){for(int g=0;g<p.G;++g)for(size_t j=0;j<d;++j)p.groups[g].M.v[j]=sum[g].v[j]/e.ng[g];return;}
  if(p.name[8]=='U'||p.G==1){
    for(int g=0;g<p.G;++g){Component& a=p.groups[g];double denom=e.sa[g]*e.sb[g]-e.ng[g]*e.ng[g];
      if(denom>1e-12*e.sa[g]*e.sb[g])for(size_t j=0;j<d;++j){a.M.v[j]=(e.sa[g]*bsum[g].v[j]-e.ng[g]*sum[g].v[j])/denom;a.A.v[j]=(e.sb[g]*sum[g].v[j]-e.ng[g]*bsum[g].v[j])/denom;}
      else for(size_t j=0;j<d;++j)a.M.v[j]=(bsum[g].v[j]-e.ng[g]*a.A.v[j])/e.sb[g];
    }
  }else if(p.name[0]=='C'&&p.name[1]=='C'&&p.name[4]=='C'&&p.name[5]=='C'){
    // Both covariance matrices are common: their precisions cancel from the
    // common-skewness equation. Solve directly instead of iterating with PCG.
    Mat numerator(p.n,p.p),A=p.groups[0].A;double weight=0;
    for(int g=0;g<p.G;++g){
      weight+=std::max(0.0,e.sa[g]-e.ng[g]*e.ng[g]/e.sb[g]);
      for(size_t j=0;j<d;++j)numerator.v[j]+=sum[g].v[j]-e.ng[g]*bsum[g].v[j]/e.sb[g];
    }
    if(weight>1e-12)for(size_t j=0;j<d;++j)A.v[j]=numerator.v[j]/weight;
    for(int g=0;g<p.G;++g){p.groups[g].A=A;for(size_t j=0;j<d;++j)p.groups[g].M.v[j]=(bsum[g].v[j]-e.ng[g]*A.v[j])/e.sb[g];}
  }else{
    std::vector<Precision> s,t;V weights(p.G);Mat rhs(p.n,p.p);
    for(int g=0;g<p.G;++g){const Component& a=p.groups[g];s.emplace_back(a.sig,a.L);t.emplace_back(a.psi,a.D);
      weights[g]=std::max(0.0,e.sa[g]-e.ng[g]*e.ng[g]/e.sb[g]);Mat v=sum[g];
      for(size_t j=0;j<d;++j)v.v[j]-=e.ng[g]*bsum[g].v[j]/e.sb[g];
      Mat add=t[g].right(s[g].left(v));for(size_t j=0;j<d;++j)rhs.v[j]+=add.v[j];
    }
    Mat A=p.groups[0].A;
    if(std::accumulate(weights.begin(),weights.end(),0.0)>1e-12)A=pcg(A,rhs,s,t,weights,c,pcg_iterations,p);
    for(int g=0;g<p.G;++g){p.groups[g].A=A;for(size_t j=0;j<d;++j)p.groups[g].M.v[j]=(bsum[g].v[j]-e.ng[g]*A.v[j])/e.sb[g];}
  }
}

} // namespace pm
