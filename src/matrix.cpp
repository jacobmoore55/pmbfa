#include "pmbfa_types.h"

namespace pm {
Mat transpose(const Mat& x) {
  Mat o(x.p,x.n);for(int j=0;j<x.p;++j)for(int i=0;i<x.n;++i)o(j,i)=x(i,j);return o;
}
Mat multiply(const Mat& a,const Mat& b) {
  need(a.p==b.n,"Internal matrix dimensions do not conform.");
  Mat out(a.n,b.p);
  for(int j=0;j<b.p;++j) {
    for(int k=0;k<a.p;++k) {
      const double value=b(k,j);
      if(value==0)continue;
      for(int i=0;i<a.n;++i)out(i,j)+=a(i,k)*value;
    }
  }
  return out;
}
double dot(const V& a,const V& b) {
  double s=0;for(size_t k=0;k<a.size();++k)s+=a[k]*b[k];return s;
}
Mat identity(int n){Mat a(n,n);for(int j=0;j<n;++j)a(j,j)=1;return a;}
Mat chol(const Mat& a) {
  need(a.n==a.p,"Internal Cholesky dimensions."); Mat l(a.n,a.n);
  for(int j=0;j<a.n;++j){
    double s=a(j,j);for(int k=0;k<j;++k)s-=l(j,k)*l(j,k);
    need(s>0&&finite(s),"A conditional second moment is not positive definite.");
    l(j,j)=std::sqrt(s);
    for(int i=j+1;i<a.n;++i){
      double t=.5*(a(i,j)+a(j,i));
      for(int k=0;k<j;++k)t-=l(i,k)*l(j,k);
      l(i,j)=t/l(j,j);
    }
  }return l;
}
Mat forward(const Mat& l,const Mat& b) {
  Mat out=b;
  for(int j=0;j<b.p;++j) {
    for(int i=0;i<l.n;++i) {
      for(int k=0;k<i;++k)out(i,j)-=l(i,k)*out(k,j);
      out(i,j)/=l(i,i);
    }
  }
  return out;
}
Mat backward(const Mat& l,const Mat& b) {
  Mat out=b;
  for(int j=0;j<b.p;++j) {
    for(int i=l.n-1;i>=0;--i) {
      for(int k=i+1;k<l.n;++k)out(i,j)-=l(k,i)*out(k,j);
      out(i,j)/=l(i,i);
    }
  }
  return out;
}
Mat solve(const Mat& a,const Mat& b){Mat l=chol(a);return backward(l,forward(l,b));}
NumericMatrix to_r(const Mat& m){NumericMatrix o(m.n,m.p);std::copy(m.v.begin(),m.v.end(),o.begin());return o;}
Mat from_r(SEXP x,int n,int p,const char* label){
  if(!Rf_isMatrix(x)||!(TYPEOF(x)==REALSXP||TYPEOF(x)==INTSXP))stop("%s must be a numeric matrix.",label);
  NumericMatrix m=as<NumericMatrix>(x);
  if(m.nrow()!=n||m.ncol()!=p)stop("%s has incorrect dimensions.",label);
  Mat o(n,p);std::copy(m.begin(),m.end(),o.v.begin());
  for(double z:o.v)if(!finite(z))stop("%s must be finite.",label);return o;
}
Precision::Precision(const V& d,const Mat& loading){
    int n=d.size(),k=loading.p;inv.resize(n);Mat t(k,n),c=identity(k);
    for(int i=0;i<n;++i){need(d[i]>0&&finite(d[i]),"Nonpositive diagonal variance.");inv[i]=1/d[i];logdet+=std::log(d[i]);
      for(int u=0;u<k;++u)t(u,i)=loading(i,u)*inv[i];}
    for(int j=0;j<k;++j)for(int i=0;i<k;++i)for(int h=0;h<n;++h)c(i,j)+=t(i,h)*loading(h,j);
    Mat l=chol(c);for(int i=0;i<k;++i)logdet+=2*std::log(l(i,i));
    K=forward(l,t);Kt=transpose(K);B=backward(l,K);W=backward(l,forward(l,identity(k)));
  }
Mat Precision::left(const Mat& a)const{
    need((int)inv.size()==a.n,"Internal precision dimensions.");
    Mat out=a,kx=multiply(K,a);
    for(int j=0;j<a.p;++j){
      for(int i=0;i<a.n;++i)out(i,j)*=inv[i];
      for(int u=0;u<K.n;++u){double z=kx(u,j);for(int i=0;i<a.n;++i)out(i,j)-=Kt(i,u)*z;}
    }return out;
  }
Mat Precision::right(const Mat& a)const{
    need((int)inv.size()==a.p,"Internal precision dimensions.");Mat out=a,ak=multiply(a,Kt);
    for(int j=0;j<a.p;++j){
      for(int i=0;i<a.n;++i)out(i,j)*=inv[j];
      for(int u=0;u<K.n;++u){double z=K(u,j);for(int i=0;i<a.n;++i)out(i,j)-=ak(i,u)*z;}
    }return out;
  }

} // namespace pm
