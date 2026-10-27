test_that("common covariance location and skewness updates match independent moments", {
pmbfa_initializations <- native("pmbfa_initializations")
pmbfa_evaluate <- native("pmbfa_evaluate")
pmbfa_posterior_moments <- native("pmbfa_posterior_moments")
pmsbfaEM <- native("pmsbfaEM")
set.seed(819);x <- array(rnorm(5*6*60),c(5,6,60))
flat <- matrix(x,nrow=30)
for(family in c('ST','GH','NIG','VG')) {
 p0 <- pmbfa_initializations(x,2,1,1,family,
   list(ncandid=4,nstart=2,nrep=1,init.iter=0,min.component=0))[[1]]
 for(row in c('CCC','CCU'))for(column in c('CCC','CCU')) {
  p0$moddat$name <- paste(row,column,'C','U',sep='-')
  p <- pmbfa_evaluate(x,p0)$params
  e <- pmbfa_posterior_moments(x,p)
  numerator <- numeric(30);weight <- 0; bsums <- list()
  for(g in 1:2) {
   bsums[[g]] <- drop(flat%*%(e$z[,g]*e$b[,g]))
   numerator <- numerator+drop(flat%*%e$z[,g])-e$ng[g]*bsums[[g]]/e$sb[g]
   weight <- weight+max(0,e$sa[g]-e$ng[g]^2/e$sb[g])
  }
  expected_A <- numerator/weight
  f <- pmsbfaEM(x,p,list(max.iter=1,min.component=0))
  stopifnot(f$iter==1,f$diagnostics$pcg_iterations==0)
  for(g in 1:2)stopifnot(
   max(abs(f$params$args[[g]]$A-expected_A))<1e-9,
   max(abs(f$params$args[[g]]$M-(bsums[[g]]-e$ng[g]*expected_A)/e$sb[g]))<1e-9)
 }
}
cat('Direct shared-covariance location/skewness update: 16 family/model cases passed.\n')

expect_true(TRUE)
})
