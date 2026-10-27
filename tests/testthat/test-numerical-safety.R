test_that("the legacy effective-size bound is strictly greater than six", {
  d<-numerical_fixture("Normal",N=60L)
  d$params$moddat$class<-c(rep(1,6),rep(2,54))
  f<-native("pmbfaEM")(d$x,d$params,list(max.iter=1L))
  initial<-native("pmbfa_evaluate")(d$x,d$params)
  expect_identical(f$status,"insufficient_component_size")
  expect_false(f$converged);expect_identical(f$score,Inf)
  expect_equal(f$diagnostics$effective_sizes,c(6,54))
  expect_identical(f$params,initial$params)
  expect_identical(f$loglik,initial$loglik)
  expect_identical(f$posterior,initial$posterior)
  d$params$moddat$class<-c(rep(1,7),rep(2,53))
  allowed<-native("pmbfaEM")(d$x,d$params,list(max.iter=1L))
  expect_equal(allowed$iter,1L)
  expect_true(all(allowed$diagnostics$effective_sizes>6))
})

test_that("a final proposed collapse is rejected before convergence or acceptance", {
  set.seed(84);x<-array(rnorm(4*4*60),c(4,4,60))
  d<-numerical_fixture("Normal",N=60L);p<-d$params
  p$args[[1]]$M<-matrix(.3,4,4);p$args[[2]]$M<-matrix(0,4,4)
  p$gmp$pik<-c(.2,.8)
  initial<-native("pmbfa_evaluate")(x,p)
  one<-native("pmbfaEM")(x,p,list(max.iter=1L,min.component=0))
  before<-min(initial$diagnostics$effective_sizes)
  after<-min(one$diagnostics$effective_sizes)
  expect_lt(after,before)
  threshold<-(before+after)/2
  guarded<-native("pmbfaEM")(x,p,list(max.iter=1L,min.iter=1L,
     tol=1e10,min.component=threshold))
  expect_identical(guarded$status,"insufficient_component_size")
  expect_false(guarded$converged);expect_equal(guarded$iter,0L)
  expect_identical(guarded$score,Inf)
  expect_true(any(guarded$diagnostics$rejected_effective_sizes<=threshold))
  expect_true(all(guarded$diagnostics$effective_sizes>threshold))
  expect_identical(guarded$params,initial$params)
  expect_identical(guarded$posterior,initial$posterior)
  expect_identical(guarded$loglik,initial$loglik)
})

test_that("extreme finite loadings return a silent diagnostic failure", {
  d<-numerical_fixture("ST",N=60L);d$params$gmp$Lambdak[]<-1e200
  expect_output(f<-native("pmsbfaEM")(d$x,d$params,list(max.iter=2L)),NA)
  expect_s3_class(f,"pmbfa");expect_false(f$converged)
  expect_identical(f$status,"initialization_failed")
  expect_length(f$loglik,0L);expect_identical(f$BIC,Inf)
  expect_identical(f$ICL,Inf);expect_identical(f$score,Inf)
  expect_match(f$message,"positive definite")
  expect_true(all(is.finite(f$params$gmp$Lambdak)))
})

test_that("finite raw information criteria never make a failed fit selectable", {
  d<-numerical_fixture("Normal",N=60L)
  f<-native("pmbfaEM")(d$x,d$params,list(max.iter=1L,min.component=100))
  expect_true(is.finite(f$bic));expect_identical(f$score,Inf)
  valid<-f;valid$converged<-TRUE;valid$status<-"converged";valid$bic<-100
  f$bic<--1e12
  expect_equal(native(".pmbfa_best")(list(f,valid),"BIC")$bic,100)
  failure<-native(".pmbfa_failure")(d$x,2,1,1,"Normal","UUU-UUU","test")
  expect_equal(failure$BIC,failure$bic);expect_equal(failure$ICL,failure$icl)
})

test_that("agreed variance and degrees-of-freedom bounds remain in force", {
  d<-numerical_fixture("ST",N=60L)
  d$params$gmp$Sigmak[]<-1e-12;d$params$gmp$Psik[]<-1e-12
  for(g in 1:2)d$params$args[[g]]$nu<-2
  f<-native("pmsbfaEM")(d$x,d$params,list(max.iter=1L,min.component=0))
  expect_true(all(f$params$gmp$Sigmak>=.1))
  expect_true(all(f$params$gmp$Psik>=.1))
  expect_true(all(vapply(f$params$args,function(a)a$nu,numeric(1))>=4.1))
})
