test_that("mixing means and finite-difference likelihood scores agree", {
pmbfa_initializations <- native("pmbfa_initializations")
pmbfa_evaluate <- native("pmbfa_evaluate")
pmbfa_posterior_moments <- native("pmbfa_posterior_moments")
pmbfa_component_mean <- native("pmbfa_component_mean")
pmbfa_mean_diagnostics <- native("pmbfa_mean_diagnostics")
.pmbfa_prior_mixing_mean <- native(".pmbfa_prior_mixing_mean")
invisible(capture.output({
set.seed(7341)
x <- array(rnorm(4L*4L*60L),c(4L,4L,60L))
x[,,31:60] <- x[,,31:60]+2
checks <- 0L
check <- function(value) { expect_true(isTRUE(value)); checks <<- checks+1L }
throws <- function(expr,pattern) {
  message <- tryCatch({force(expr);""},error=conditionMessage)
  check(grepl(pattern,message,fixed=TRUE))
}
control <- list(ncandid=3L,nstart=2L,nrep=1L,init.iter=0L,
                min.component=0,tol=1e-8)
fitted <- list()
for(family in c("Normal","ST","GH","NIG","VG")) {
  params <- pmbfa_initializations(x,2L,1L,1L,family,control)[[1L]]
  params$moddat$name <- if(family=="Normal")"UUU-UUU"else"UUU-UUU-U-U"
  if(family!="Normal") for(g in 1:2) params$args[[g]]$A <- matrix(.1*g,4,4)
  fit <- pmbfa_evaluate(x,params,control)
  moments <- pmbfa_posterior_moments(x,fit$params,control)
  diagnostics <- pmbfa_mean_diagnostics(fit,x,moments)
  fitted[[family]] <- list(fit=fit,moments=moments,diagnostics=diagnostics)
  check(abs(diagnostics$posterior_loglik-tail(fit$loglik,1L))<1e-9)
  check(max(abs(rowSums(moments$z)-1))<1e-12)
  for(g in 1:2) {
    check(max(abs(pmbfa_component_mean(fit,g)-diagnostics$matrices[[g]]$population))<1e-12)
    check(abs(pmbfa_component_mean(fit,g,TRUE)-diagnostics$components$population_mean[g])<1e-12)
    for(parameter in if(family=="Normal")"M"else c("M","A")) {
      h <- 1e-5
      plus <- minus <- fit$params
      plus$args[[g]][[parameter]][1,2] <- plus$args[[g]][[parameter]][1,2]+h
      minus$args[[g]][[parameter]][1,2] <- minus$args[[g]][[parameter]][1,2]-h
      finite_difference <- (tail(pmbfa_evaluate(x,plus,control)$loglik,1L)-
                           tail(pmbfa_evaluate(x,minus,control)$loglik,1L))/(2*h)
      analytic <- diagnostics$matrices[[g]][[if(parameter=="M")"location_gradient"else"skew_gradient"]][1,2]*
        diagnostics$components$posterior_n[g]
      check(abs(finite_difference-analytic)<1e-4*(1+abs(analytic)))
    }
  }
  if(family=="Normal") {
    check(all(is.na(diagnostics$components$skew_score_rmse)))
    check(is.na(diagnostics$pooled_skew_gradient_rmse))
    check(grepl("Normal",diagnostics$interpretation))
  } else {
    # Perturb all A matrices together for the common-A score.
    common <- fit$params
    common$moddat$name <- "UUU-UUU-C-U"
    common$args[[2]]$A <- common$args[[1]]$A
    common <- pmbfa_evaluate(x,common,control)
    d <- pmbfa_mean_diagnostics(common,x)
    plus <- minus <- common$params
    h <- 1e-5
    for(g in 1:2) {
      plus$args[[g]]$A[1,2] <- plus$args[[g]]$A[1,2]+h
      minus$args[[g]]$A[1,2] <- minus$args[[g]]$A[1,2]-h
    }
    finite_difference <- (tail(pmbfa_evaluate(x,plus,control)$loglik,1L)-
                         tail(pmbfa_evaluate(x,minus,control)$loglik,1L))/(2*h)
    analytic <- sum(vapply(1:2,function(g)d$matrices[[g]]$skew_gradient[1,2]*
      d$components$posterior_n[g],numeric(1)))
    check(abs(finite_difference-analytic)<1e-4*(1+abs(analytic)))
    check(all(d$components$shared_skewness))
  }
}
check(abs(.pmbfa_prior_mixing_mean("ST",list(nu=5))-5/3)<1e-15)
check(abs(.pmbfa_prior_mixing_mean("NIG",list(gam=.7))-1/.7)<1e-15)
check(.pmbfa_prior_mixing_mean("VG",list())==1)
for(lambda in c(-4,0,2,20)) {
  direct <- besselK(3,lambda+1,expon.scaled=TRUE)/besselK(3,lambda,expon.scaled=TRUE)
  check(abs(.pmbfa_prior_mixing_mean("GH",list(omega=3,lambda=lambda))/direct-1)<1e-13)
}
throws(.pmbfa_prior_mixing_mean("ST",list(nu=2)),"nu > 2")
throws(.pmbfa_prior_mixing_mean("NIG",list(gam=0)),"gam > 0")
throws(.pmbfa_prior_mixing_mean("GH",list(omega=NA,lambda=2)),"omega must be positive")
throws(.pmbfa_prior_mixing_mean("other",list()),"Unknown distribution")
throws(pmbfa_component_mean(fitted$ST$fit,0),"g must identify")
throws(pmbfa_mean_diagnostics(fitted$ST$fit,x,list(status="numerical_error",message="failed")),
       "Posterior moment evaluation failed: failed")
bad <- fitted$ST$moments;bad$a <- NULL
throws(pmbfa_mean_diagnostics(fitted$ST$fit,x,bad),"Posterior moments must be finite")
bad <- fitted$ST$moments;bad$z[,1] <- 0
throws(pmbfa_mean_diagnostics(fitted$ST$fit,x,bad),"Posterior moments must be finite")

# Solve the unrestricted location/skewness CM equations for fixed moments.
# This tests the mean identity, not convergence of those moments with the new params.
fit <- fitted$ST$fit;moments <- fitted$ST$moments
xx <- matrix(x,16L,60L)
for(g in 1:2) {
  ng <- sum(moments$z[,g]);sa <- sum(moments$z[,g]*moments$a[,g])
  sb <- sum(moments$z[,g]*moments$b[,g])
  sx <- matrix(xx%*%moments$z[,g],4L,4L)
  sbx <- matrix(xx%*%(moments$z[,g]*moments$b[,g]),4L,4L)
  fit$params$args[[g]]$M <- (sbx-ng*sx/sa)/(sb-ng^2/sa)
  fit$params$args[[g]]$A <- (sx-ng*fit$params$args[[g]]$M)/sa
}
d <- pmbfa_mean_diagnostics(fit,x,moments)
check(max(d$components$conditional_vs_posterior_rmse)<1e-12)
check(max(d$components$location_score_rmse)<1e-12)
check(max(d$components$skew_score_rmse)<1e-12)
check(max(d$components$population_vs_posterior_rmse)>1e-6)
cat(checks,"mean diagnostic checks passed.\n")

}))
})
