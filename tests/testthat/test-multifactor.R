test_that("multiple-factor shared loading updates agree with a dense R reference", {
  set.seed(581)
  x<-array(rnorm(6*7*90),c(6,7,90));x[,,46:90]<-x[,,46:90]+2
  start<-native("pmbfa_initializations")(x,2,2,2,"Normal",
    list(ncandid=3,nstart=1,nrep=1,init.iter=0,min.component=0))[[1]]
  start$gmp$Sigmak[]<-runif(12,.4,1.8)
  start$gmp$Psik[]<-runif(14,.4,1.8)
  for(model in c("CUU-CUU","CUU-UUU","UUU-CUU")) {
    start$moddat$name<-model
    fit<-native("pmbfaEM")(x,start,list(max.iter=2,tol=1e-14,min.component=0))
    reference<-dense_normal(x,start,cycles=2L)
    expect_equal(fit$iter,2L)
    expect_equal(fit$loglik,reference$loglik,tolerance=1e-10,info=model)
    expect_equal(fit$posterior,reference$posterior,tolerance=1e-10,info=model)
    for(key in c("Lambdak","Deltak","Sigmak","Psik"))
      expect_equal(fit$params$gmp[[key]],reference$params$gmp[[key]],tolerance=1e-9,info=paste(model,key))
  }
})
