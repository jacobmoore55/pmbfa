test_that("the installed namespace exposes exactly the four fitting functions", {
  expect_setequal(getNamespaceExports("pmbfa"),c("pmbfa","pmsbfa","search_pmbfa","search_pmsbfa"))
  expect_false(exists(".pmbfa_require_core",asNamespace("pmbfa"),inherits=FALSE))
  routines <- getDLLRegisteredRoutines(getLoadedDLLs()[["pmbfa"]])$.Call
  expect_length(routines,7L)
  expect_false(getLoadedDLLs()[["pmbfa"]][["dynamicLookup"]])
})

test_that("both fitting and factor-search APIs work through registered native code", {
  d <- numerical_fixture("Normal",N=60L)
  a <- pmbfa(d$x,2,1,1,modelSubset="CUC-CUC",control=quick_control())
  b <- pmsbfa(d$x,2,1,1,family="ST",modelSubset="CUC-CUC-C-U",control=quick_control())
  c <- search_pmbfa(d$x,2,max_q=1,max_r=1,modelSubset="CUC-CUC",control=quick_control())
  e <- search_pmsbfa(d$x,2,max_q=1,max_r=1,modelSubset="CUC-CUC-C-U",control=quick_control())
  for(f in list(a,b,c,e)) {
    expect_s3_class(f,"pmbfa")
    expect_length(f$loglik,f$iter+1L)
    expect_true(all(is.finite(f$loglik)))
    expect_true(all(diff(f$loglik)>=0))
    expect_equal(rowSums(f$posterior),rep(1,60),tolerance=1e-10)
    expect_equal(BIC(f),f$bic)
    expect_equal(attr(logLik(f),"df"),f$npar)
    expect_equal(nobs(f),60L)
    if(!f$converged)expect_identical(f$score,Inf)
    expect_equal(f$BIC,f$bic);expect_equal(f$ICL,f$icl)
  }
  expect_s3_class(b,"pmsbfa")
  expect_s3_class(e,"pmsbfa")
  expect_equal(c$factor.search$q,1L)
  expect_equal(e$factor.search$r,1L)
})

test_that("invalid public inputs fail before native numerical fitting", {
  d<-numerical_fixture("Normal",N=60L);x<-d$x; ctl<-quick_control(max.iter=0L)
  expect_error(pmbfa(x[,,1:10],2,1,1,control=ctl),"more than 10")
  expect_error(pmbfa(x,G=2.5,q=1,r=1,control=ctl),"integer")
  expect_error(pmbfa(x,2,2,1,control=ctl),"q must")
  expect_error(pmbfa(x,2,1,1,modelSubset="bad",control=ctl),"invalid constraint")
  expect_error(pmbfa(x,2,1,1,control=list(tol=-1)),"tol")
  expect_error(pmbfa(x,2,1,1,control=list(wrong=1)),"Unknown control")
  expect_error(pmsbfa(x,2,1,1,family="unknown",control=ctl),"arg")
  expect_error(pmbfa(x,2,1,1,class=1:2,control=ctl),"one label")
  expect_error(pmbfa(x,2,1,1,inits=list(d$params),control=ctl,unused=1),"Unused")
  x[1]<-NA_real_;expect_error(pmbfa(x,2,1,1,control=ctl),"finite")
})

test_that("seeds, grids, supplied starts and known labels remain reproducible", {
  d<-numerical_fixture("Normal",N=60L);x<-d$x;ctl<-quick_control(max.iter=3L)
  set.seed(31);before<-.Random.seed
  a<-pmbfa(x,2,1,1,modelSubset="UUU-UUU",control=ctl)
  expect_identical(.Random.seed,before)
  b<-pmbfa(x,2,1,1,modelSubset="UUU-UUU",control=ctl)
  expect_identical(a$params,b$params);expect_identical(a$loglik,b$loglik)
  f<-pmbfa(x,2,1,1,class=d$labels,modelSubset="UUU-UUU",inits=d$params,control=ctl)
  expect_identical(as.integer(f$cluster),as.integer(d$labels))
  g<-pmbfa(x,G=1:2,q=1,r=1,modelSubset="UUU-UUU",control=ctl)
  expect_s3_class(g,"pmbfa_fits");expect_length(g,2L)
  expect_false(is.unsorted(vapply(g,function(z)z$score,numeric(1))))
  expect_identical(capture.output(print(a))[3],sprintf(
    "  logLik: %.2f | BIC: %.2f | ICL: %.2f | parameters: %s",
    tail(a$loglik,1),a$bic,a$icl,a$npar))
})

test_that("fork and socket workers return the same fixed-seed fits", {
  d<-numerical_fixture("Normal",N=60L);ctl<-quick_control(max.iter=2L)
  serial<-pmbfa(d$x,G=1:2,q=1,r=1,modelSubset="UUU-UUU",control=ctl)
  if(.Platform$OS.type=="unix") {
    parallel<-pmbfa(d$x,G=1:2,q=1,r=1,modelSubset="UUU-UUU",control=ctl,parallelize=TRUE,ncores=2)
    expect_identical(lapply(serial,function(f)f$loglik),lapply(parallel,function(f)f$loglik))
  }
  labels<-native(".pmbfa_labels")(NULL,60,1:2)
  controls<-native(".pmbfa_controls")(ctl,list(),FALSE,60,4,4)
  arguments<-list(x=d$x,family="Normal",labels=labels,modelSubset="UUU-UUU",adaptive=TRUE,
    adaptive.method="both",criterion="BIC",seed=123456L,inits=NULL,control=controls,screening="none")
  socket<-native(".pmbfa_psock")(data.frame(G=1:2,q=1L,r=1L),arguments,2L)
  expected<-serial[order(vapply(serial,function(f)f$params$moddat$G,numeric(1)))]
  expect_identical(lapply(expected,function(f)f$loglik),lapply(socket,function(f)f$loglik))
})
