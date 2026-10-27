check_elapsed <- function(value) {
  expect_type(value,"double")
  expect_length(value,1L)
  expect_true(is.finite(value)&&value>=0)
}

test_that("scalar fit calls retain the actual user expression and elapsed seconds", {
  d<-numerical_fixture("Normal",N=60L);x<-d$x;ctl<-quick_control(max.iter=0L,init.iter=0L)
  a<-pmbfa(x=x,G=2,q=1,r=1,modelSubset="UUU-UUU",control=ctl)
  expect_identical(a$call,quote(pmbfa(x=x,G=2,q=1,r=1,modelSubset="UUU-UUU",control=ctl)))
  check_elapsed(a$elapsed)
  b<-pmsbfa(x=x,G=2,q=1,r=1,modelSubset="UUU-UUU-U-U",control=ctl)
  expect_identical(b$call,quote(pmsbfa(x=x,G=2,q=1,r=1,modelSubset="UUU-UUU-U-U",control=ctl)))
  check_elapsed(b$elapsed)
  expect_output(print(a),"Elapsed: [0-9]+\\.[0-9]{2} seconds")
  expect_output(print(a),"Call: pmbfa")
})

test_that("grid metadata uses attributes without adding fake fit entries", {
  d<-numerical_fixture("Normal",N=60L);x<-d$x;ctl<-quick_control(max.iter=0L,init.iter=0L)
  fits<-pmbfa(x=x,G=1:2,q=1,r=1,modelSubset="UUU-UUU",control=ctl)
  expected<-quote(pmbfa(x=x,G=1:2,q=1,r=1,modelSubset="UUU-UUU",control=ctl))
  expect_length(fits,2L)
  expect_identical(attr(fits,"call"),expected)
  check_elapsed(attr(fits,"elapsed"))
  for(fit in fits) {
    expect_s3_class(fit,"pmbfa")
    expect_identical(fit$call,expected)
    check_elapsed(fit$elapsed)
  }
})

test_that("factor searches retain the outer call without literal data or closures", {
  d<-numerical_fixture("Normal",N=60L);x<-d$x;ctl<-quick_control(max.iter=0L,init.iter=0L)
  fits<-search_pmbfa(x=x,G=1:2,max_q=1,max_r=1,modelSubset="UUU-UUU",control=ctl)
  expected<-quote(search_pmbfa(x=x,G=1:2,max_q=1,max_r=1,modelSubset="UUU-UUU",control=ctl))
  expect_identical(attr(fits,"call"),expected)
  check_elapsed(attr(fits,"elapsed"))
  for(fit in fits) {
    expect_identical(fit$call,expected)
    expect_lt(as.numeric(object.size(fit$call)),5000)
    check_elapsed(fit$elapsed)
  }
  one<-search_pmsbfa(x=x,G=2,max_q=1,max_r=1,modelSubset="UUU-UUU-U-U",control=ctl)
  expect_identical(one$call,quote(search_pmsbfa(x=x,G=2,max_q=1,max_r=1,modelSubset="UUU-UUU-U-U",control=ctl)))
  check_elapsed(one$elapsed)
})

test_that("diagnostic failure objects also retain call and elapsed metadata", {
  d<-numerical_fixture("Normal",N=60L);x<-d$x;ctl<-quick_control(min.component=100)
  fit<-pmbfa(x=x,G=2,q=1,r=1,modelSubset="UUU-UUU",control=ctl)
  expect_false(fit$converged)
  expect_identical(fit$score,Inf)
  expect_identical(fit$call,quote(pmbfa(x=x,G=2,q=1,r=1,modelSubset="UUU-UUU",control=ctl)))
  check_elapsed(fit$elapsed)
})

test_that("the private search adapter preserves the normal-family API", {
  d<-numerical_fixture("Normal",N=60L)
  expect_error(search_pmbfa(d$x,2,max_q=1,max_r=1,family="ST"),"search_pmsbfa")
  expect_error(search_pmbfa(d$x,2,max_q=1,max_r=1,family=NULL),"search_pmsbfa")
})
