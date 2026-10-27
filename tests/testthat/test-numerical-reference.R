test_that("independent densities, moments, all constraints and native contracts hold", {
  checks <- 0L
  check <- function(condition,label) {
    checks <<- checks+1L
    expect_true(isTRUE(condition),info=label)
    invisible(condition)
  }
  near <- function(actual,expected,tolerance=1e-7)
    isTRUE(all.equal(actual,expected,tolerance=tolerance,check.attributes=FALSE))
  fixture <- numerical_fixture
  control <- native_control
  pmbfaEM <- native("pmbfaEM");pmsbfaEM <- native("pmsbfaEM")
  pmbfa_initializations <- native("pmbfa_initializations")
  fit_one <- function(x,params,opts=control())
    if(params$moddat$family=="Normal")pmbfaEM(x,params,opts)else pmsbfaEM(x,params,opts)
  invisible(capture.output({
# Marginal densities independently integrate the matrix-normal density over
# W, using the original prior density on W and a log-W integration variable.
# Precision is formed directly from full covariance matrices with solve().
log_density_reference <- function(obs, params, g = 1L, moment = NULL) {
  md <- params$moddat; gm <- params$gmp; ag <- params$args[[g]]
  L <- matrix(gm$Lambdak[g,,],md$n,md$q)
  D <- matrix(gm$Deltak[g,,],md$p,md$r)
  Sigma <- diag(gm$Sigmak[,g],md$n) + tcrossprod(L)
  Psi <- diag(gm$Psik[,g],md$p) + tcrossprod(D)
  covariance <- kronecker(Psi,Sigma)
  precision <- solve(covariance)
  v <- as.vector(obs-ag$M); d <- length(v)
  A <- if(md$family=='Normal') rep(0,d) else as.vector(ag$A)
  constant <- -.5*d*log(2*pi) - .5*as.numeric(determinant(covariance,logarithm=TRUE)$modulus)
  delta <- drop(crossprod(v,precision%*%v))
  tilt <- drop(crossprod(v,precision%*%A))
  rho <- drop(crossprod(A,precision%*%A))
  if(md$family=='Normal') return(constant-.5*delta)
  logprior <- switch(md$family,
    ST = function(w) (ag$nu/2)*log(ag$nu/2)-lgamma(ag$nu/2) -
      (ag$nu/2+1)*log(w)-ag$nu/(2*w),
    GH = function(w) {
      z<-ag$omega
      logK<-log(besselK(z,ag$lambda,expon.scaled=TRUE))-z
      -log(2)-logK+(ag$lambda-1)*log(w)-.5*z*(w+1/w)
    },
    NIG = function(w) ag$gam-.5*log(2*pi)-1.5*log(w)-.5*(1/w+ag$gam^2*w),
    VG = function(w) dgamma(w,shape=ag$gamma,rate=ag$gamma,log=TRUE))
  integrand <- function(t) {
    w<-exp(t)
    exp(constant+tilt-.5*delta/w-.5*rho*w-.5*d*t+logprior(w)+t)
  }
  mass <- integrate(integrand,-30,30,subdivisions=1000L,rel.tol=1e-10,abs.tol=1e-100)$value
  if (is.null(moment)) return(log(mass))
  weighted <- function(t) integrand(t) * switch(moment,a=exp(t),b=exp(-t),c=t)
  integrate(weighted,-30,30,subdivisions=1000L,rel.tol=1e-9,abs.tol=1e-12)$value/mass
}
for (family in c('Normal','ST','GH','NIG','VG')) {
  d<-fixture(family,n=3L,p=3L,N=12L,G=1L)
  # Check several nontrivial deviations at fixed parameters, retaining N>10.
  d$x[,,]<-rep(c(.2,-.5,.7,1.1),length.out=length(d$x))
  d$params$args[[1]]$M<-matrix(rep(c(.1,.3,-.1,.2),length.out=9),3)
  evaluated<-fit_one(d$x,d$params,control(max.iter=0L))
  expected<-sum(vapply(seq_len(dim(d$x)[3]),function(i)
    log_density_reference(d$x[,,i],evaluated$params),numeric(1)))
  check(near(tail(evaluated$loglik,1),expected,1e-7),paste(family,'density disagrees with independent W integration'))
}
cat('Independent densities: all five families checked.\n')

# Numerically integrate all three conditional W moments, then independently
# solve the first conditional maximization and the scalar shape equations.
for (family in c('ST','GH','NIG','VG')) {
  d<-fixture(family,n=3L,p=3L,N=12L,G=1L)
  initial<-fit_one(d$x,d$params,control(max.iter=0L))
  moments<-vapply(c('a','b','c'),function(moment) vapply(seq_len(12L),
    function(i)log_density_reference(d$x[,,i],initial$params,moment=moment),numeric(1)),numeric(12L))
  sa<-sum(moments[,'a']);sb<-sum(moments[,'b']);sc<-sum(moments[,'c'])
  sx<-apply(d$x,1:2,sum)
  sbx<-apply(sweep(d$x,3,moments[,'b'],'*'),1:2,sum)
  denominator<-sa*sb-12^2
  expected_M<-(sa*sbx-12*sx)/denominator
  expected_A<-(sb*sx-12*sbx)/denominator
  result<-fit_one(d$x,initial$params,control(max.iter=1L,backtracks=0L))
  check(result$iter==1L,paste(family,'independent conditional moments test accepted cycle'))
  check(near(result$params$args[[1]]$M,expected_M,1e-6),paste(family,'location update from integrated conditional moments'))
  check(near(result$params$args[[1]]$A,expected_A,1e-6),paste(family,'skewness update from integrated conditional moments'))
  if(family %in% c('ST','VG')) {
    scale<-if(family=='ST')2 else 1
    lower<-if(family=='ST')4.1 else .5
    target<-if(family=='ST')(sb+sc)/12 else (sa-sc)/12
    equation<-function(z)log(z/scale)+1-digamma(z/scale)-target
    expected<-if(equation(lower)<=0)lower else uniroot(equation,c(lower,1e8),tol=1e-9)$root
    actual<-if(family=='ST')result$params$args[[1]]$nu else result$params$args[[1]]$gamma
    check(near(actual,expected,1e-6),paste(family,'shape update independently solved'))
  }
  if(family=='NIG')check(near(result$params$args[[1]]$gam,max(.5,12/sa),1e-6),'NIG shape update from integrated moments')
}
cat('Conditional GIG moments: all four skew families verified through independent location/shape updates.\n')

# The symmetric t boundary rho=0 has a closed multivariate-t density and an
# inverse-gamma conditional. A full cycle also exercises its moments.
symmetric<-fixture('ST',n=3L,p=3L,N=12L,G=1L)
symmetric$params$args[[1]]$A[]<-0
sym_eval<-fit_one(symmetric$x,symmetric$params,control(max.iter=0L))
check(near(tail(sym_eval$loglik,1),sum(vapply(seq_len(12L),function(i)
  log_density_reference(symmetric$x[,,i],sym_eval$params),numeric(1))),1e-7),'Symmetric t density rho=0')
sym_fit<-fit_one(symmetric$x,symmetric$params,control(max.iter=2L))
check(length(sym_fit$loglik)>1L && all(is.finite(sym_fit$loglik)), 'Symmetric t conditional moments permit fitting')
check(all(diff(sym_fit$loglik)>=-1e-9),'Symmetric t likelihood is monotone')

# Exhaust all constraint strings on one deterministic data/initialization pair.
constraint_models <- function(skew) {
  flags<-as.matrix(expand.grid(rep(list(c('C','U')),if(skew)8 else 6),stringsAsFactors=FALSE))
  apply(flags,1,function(z) paste0(paste0(z[1:3],collapse=''),'-',
    paste0(z[4:6],collapse=''),if(skew)paste0('-',z[7],'-',z[8])else''))
}
check_constraints <- function(params,name) {
  gm<-params$gmp;G<-params$moddat$G
  equal_groups<-function(a) all(vapply(seq_len(G),function(g) near(a[g,,],a[1,,],1e-8),logical(1)))
  if(substr(name,1,1)=='C')check(equal_groups(gm$Lambdak),paste(name,'shared row loadings'))
  if(substr(name,5,5)=='C')check(equal_groups(gm$Deltak),paste(name,'shared column loadings'))
  for(side in 1:2) {
    offset<-if(side==1)0 else 4;v<-if(side==1)gm$Sigmak else gm$Psik
    if(substr(name,offset+2,offset+2)=='C')check(all(vapply(seq_len(G),function(g)near(v[,g],v[,1],1e-8),logical(1))),paste(name,'shared uniqueness',side))
    if(substr(name,offset+3,offset+3)=='C')check(all(apply(v,2,function(z)max(z)-min(z))<1e-8),paste(name,'isotropic uniqueness',side))
  }
  if(nchar(name)==11L) {
    if(substr(name,9,9)=='C')check(all(vapply(params$args,function(a)near(a$A,params$args[[1]]$A,1e-8),logical(1))),paste(name,'shared skewness'))
    if(substr(name,11,11)=='C')check(length(unique(vapply(params$args,function(a)a$nu,numeric(1))))==1L,paste(name,'shared degrees of freedom'))
  }
}
status_tables<-list()
for(family in c('Normal','ST')) {
  d<-fixture(family); statuses<-character(); rejected<-character()
  for(name in constraint_models(family=='ST')) {
    params<-d$params;params$moddat$name<-name
    result<-fit_one(d$x,params,control(max.iter=2L))
    check(inherits(result,'pmbfa'),paste(family,name,'S3 contract'))
    check(length(result$loglik)>1L && result$iter>=1L,paste(family,name,'did not accept a full AECM cycle'))
    check(all(diff(result$loglik)>=-1e-9),paste(family,name,'likelihood decreased'))
    check(all(result$params$gmp$Sigmak>=.05-1e-12)&&all(result$params$gmp$Psik>=.05-1e-12),paste(family,name,'variance floor'))
    check_constraints(result$params,name)
    statuses<-c(statuses,result$status)
    if(result$iter==0L)rejected<-c(rejected,name)
  }
  status_tables[[family]]<-table(statuses)
  cat(family,':',length(statuses),'constraint models,',paste(names(table(statuses)),as.integer(table(statuses)),collapse=', '),'\n')
  if(length(rejected))cat('Rejected from start:',paste(rejected,collapse=', '),'\n')
}

# A returned snapshot must reproduce both likelihood and posterior independently
# of the fit's history, including a deliberately triggered numerical stop.
for(family in c('Normal','ST','GH','NIG','VG')) {
  d<-fixture(family); original<-serialize(d$params,NULL)
  result<-fit_one(d$x,d$params,control(max.iter=3L))
  reevaluated<-fit_one(d$x,result$params,control(max.iter=0L))
  check(near(tail(result$loglik,1),tail(reevaluated$loglik,1),1e-9),paste(family,'returned state likelihood'))
  check(near(result$posterior,reevaluated$posterior,1e-9),paste(family,'returned state posterior'))
  check(identical(original,serialize(d$params,NULL)),paste(family,'mutated caller parameters'))
  check(near(rowSums(result$posterior),rep(1,24),1e-10),paste(family,'posterior normalization'))
}
d<-fixture('ST')
failed<-fit_one(d$x,d$params,control(min.component=100))
initial<-fit_one(d$x,d$params,control(max.iter=0L))
check(identical(failed$status,'insufficient_component_size')&&!failed$converged&&failed$iter==0L,'Numerical failure status')
check(identical(failed$params,initial$params)&&near(failed$loglik,initial$loglik)&&near(failed$posterior,initial$posterior),'Numerical failure retains matching snapshot')
check(nzchar(failed$message)&&is.finite(failed$bic),'Numerical failure diagnostics and BIC')
output<-capture.output(invisible(fit_one(d$x,d$params,control(min.component=100))))
check(length(output)==0L,'Numerical failure is silent')
bad<-d$params;bad$gmp$Lambdak<-array(1,c(1,1,1))
check(inherits(try(fit_one(d$x,bad),silent=TRUE),'try-error'),'Malformed dimensions raise input error')
check(inherits(try(fit_one(d$x,d$params,control(tol=-1)),silent=TRUE),'try-error'),'Malformed controls raise input error')

# Known memberships stay fixed and participate in the classification likelihood.
labelled<-d$params;labelled$moddat$class<-rep(0,24);labelled$moddat$class[c(1,2,5,6)]<-d$labels[c(1,2,5,6)]
lfit<-fit_one(d$x,labelled,control(max.iter=2L))
for(i in c(1,2,5,6))check(lfit$posterior[i,labelled$moddat$class[i]]==1&&sum(lfit$posterior[i,])==1,paste('fixed label',i))

# Tolerance is an immutable option; Aitken changes the stopping test only.
normal<-fixture('Normal')
long<-fit_one(normal$x,normal$params,control(max.iter=14L,min.iter=12L,tol=1e4))
check(identical(long$diagnostics$tolerance,1e4),'Explicit tolerance remains unchanged after iteration ten')
check(long$iter>=12L,'Iteration-ten tolerance test reached twelve iterations')
limited<-fit_one(normal$x,normal$params,control(max.iter=1L,min.iter=2L))
check(limited$status=='max_iter'&&!limited$converged&&limited$iter==1L,'Maximum iteration stop is not convergence')
afit<-fit_one(normal$x,normal$params,control(max.iter=50L,tol=.01,stopping='aitken'))
check(afit$diagnostics$stopping=='aitken'&&afit$diagnostics$tolerance==.01,'Aitken setting and tolerance')
if(afit$converged&&length(afit$loglik)>=3L&&afit$message=='') {
  gains<-tail(diff(afit$loglik),2);acceleration<-gains[2]/gains[1]
  gap<-gains[2]*acceleration/(1-acceleration)
  check(gains[2]==0||(is.finite(gap)&&gap>=0&&gap<=.01+1e-10),'Aitken convergence has an admissible estimated remaining gain')
}

# Shared starts are deterministic given the R seed and do not depend on a
# previous fitted candidate; test nearest, optional k-means, and supplied labels.
init_control<-control(ncandid=4L,nstart=2L,nrep=2L,init.iter=1L)
set.seed(31);bank1<-pmbfa_initializations(normal$x,2,1,1,'Normal',init_control)
set.seed(31);bank2<-pmbfa_initializations(normal$x,2,1,1,'Normal',init_control)
check(identical(bank1,bank2),'Initializer is deterministic under set.seed')
check(length(bank1)==2L&&all(vapply(bank1,function(p)identical(p$moddat$name,'CCC-CCC'),logical(1))),'Initializer bank structure')
set.seed(31);kbank<-pmbfa_initializations(normal$x,2,1,1,'Normal',modifyList(init_control,list(init_method='kmeans')))
check(length(kbank)==2L&&attr(kbank,'initialization')$kmeans_failures==0L,'Optional kmeans initialization')
set.seed(31);labelbank<-pmbfa_initializations(normal$x,2,1,1,'Normal',init_control,normal$labels)
check(all(vapply(labelbank,function(p)identical(p$moddat$class,as.double(normal$labels)),logical(1))),'Initializer retains fully specified labels')

  }))
  expect_gte(checks,2500L)
})
