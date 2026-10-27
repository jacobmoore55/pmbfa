test_that("normal and skewed factor searches default to both constraint directions", {
  d<-numerical_fixture("Normal",N=60L)
  ctl<-quick_control(max.iter=3L,tol=1e4,min.component=0)
  normal<-search_pmbfa(d$x,G=2,max_q=1,max_r=1,control=ctl)
  skewed<-search_pmsbfa(d$x,G=2,max_q=1,max_r=1,family="ST",control=ctl)
  for(fit in list(normal,skewed)) {
    expect_identical(fit$selection$method,"both")
    expect_identical(names(fit$selection$paths),c("forward","backward"))
  }
  backward<-search_pmbfa(d$x,G=2,max_q=1,max_r=1,adaptive.method="backward",control=ctl)
  expect_identical(names(backward$selection$paths),"backward")
  fixed<-search_pmbfa(d$x,G=2,max_q=1,max_r=1,modelSubset="CCC-CCC",control=ctl)
  expect_identical(fixed$selection$method,"exhaustive")
  expect_null(fixed$selection$paths)
})
