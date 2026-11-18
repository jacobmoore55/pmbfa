test_that("direct fits use one and two components and print their selected criterion", {
  x <- numerical_fixture(N=60L)$x
  for (entry in list(pmbfa,pmsbfa)) {
    fits <- entry(x, q=1, r=1, control=quick_control())
    expect_setequal(vapply(fits,function(f) f$params$moddat$G,numeric(1)),1:2)
    fit <- fits[[1L]]
    expect_identical(fit$criterion,"BIC")
    summary <- capture.output(print(fit))
    expect_true(any(grepl("BIC: ",summary,fixed=TRUE)))
    expect_false(any(grepl("ICL: ",summary,fixed=TRUE)))
    expect_true(any(grepl("iterations: ",summary,fixed=TRUE)))
    fit <- entry(x, G=1, q=1, r=1, criterion="ICL", control=quick_control())
    summary <- capture.output(print(fit))
    expect_true(any(grepl(sprintf("ICL: %.2f",fit$icl),summary,fixed=TRUE)))
    expect_false(any(grepl("BIC: ",summary,fixed=TRUE)))
  }
})

test_that("normal and skewed factor searches default to both constraint directions", {
  d<-numerical_fixture("Normal",N=60L)
  ctl<-quick_control(max.iter=3L,tol=1e4,min.component=0)
  normal<-search_pmbfa(d$x,G=2,max_q=1,max_r=1,control=ctl)
  skewed<-search_pmsbfa(d$x,G=2,max_q=1,max_r=1,family="ST",control=ctl)
  for(fit in list(normal,skewed)) {
    expect_identical(fit$factor.search.method,"adaptive")
    expect_identical(fit$selection$method,"both")
    expect_identical(names(fit$selection$paths),c("forward","backward"))
  }
  backward<-search_pmbfa(d$x,G=2,max_q=1,max_r=1,adaptive.method="backward",control=ctl)
  expect_identical(names(backward$selection$paths),"backward")
  fixed<-search_pmbfa(d$x,G=2,max_q=1,max_r=1,modelSubset="CCC-CCC",control=ctl)
  expect_identical(fixed$selection$method,"exhaustive")
  expect_null(fixed$selection$paths)
})
