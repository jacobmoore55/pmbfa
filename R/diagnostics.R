# Mean, score, and conditioning diagnostics for the fitted variance-mean mixture.
# Population component means are M + E_prior[W] A. They need not equal hard or
# posterior-weighted sample cluster means in a finite sample, especially when A
# is shared. Native posterior moments avoid a second Bessel implementation.

# Evaluate a stable scaled log Bessel K value for mean diagnostics.
.pmbfa_diag_logK <- function(x,nu) {
  if(length(x)!=1L || length(nu)!=1L || !is.finite(x) || x<=0 || !is.finite(nu))
    stop("GH omega must be positive and lambda must be finite.",call.=FALSE)
  value <- suppressWarnings(log(besselK(x,nu,expon.scaled=TRUE)))
  if(is.finite(value))return(value)
  v <- abs(nu)
  if(v<10 || x<=0)stop("Cannot evaluate the GH population mixing mean.")
  z <- x/v;t <- sqrt(1+z*z);u <- 1/t
  eta <- t+log(z)-log1p(t)
  x-v*eta+.5*log(pi/(2*v))-.5*log(t)+
    log1p(-(3*u-5*u^3)/(24*v))
}

# Compute the family-specific prior expectation of the mixing variable.
.pmbfa_prior_mixing_mean <- function(family,arg) {
  value <- switch(family,Normal=1,
    ST=if(length(arg$nu)==1L && is.finite(arg$nu) && arg$nu>2)
      arg$nu/(arg$nu-2) else stop("The ST population mean requires nu > 2.",call.=FALSE),
    GH=exp(.pmbfa_diag_logK(arg$omega,arg$lambda+1)-
             .pmbfa_diag_logK(arg$omega,arg$lambda)),
    NIG=if(length(arg$gam)==1L && is.finite(arg$gam) && arg$gam>0)
      1/arg$gam else stop("The NIG population mean requires gam > 0.",call.=FALSE),
    VG=1,stop("Unknown distribution family.",call.=FALSE))
  if(!is.finite(value) || value<=0)
    stop("The population mixing mean could not be evaluated reliably.",call.=FALSE)
  value
}

# Calculate a fitted population component mean matrix or its scalar pixel average.
pmbfa_component_mean <- function(fit,g,scalar=FALSE) {
  md <- fit$params$moddat
  if(length(g)!=1L || !is.finite(g) || g!=floor(g) || g<1 || g>md$G)
    stop("g must identify one component.")
  arg <- fit$params$args[[g]]
  mixing_mean <- .pmbfa_prior_mixing_mean(md$family,arg)
  A <- if(md$family=="Normal" || is.null(arg$A))matrix(0,nrow(arg$M),ncol(arg$M))else arg$A
  result <- arg$M+mixing_mean*A
  if(scalar)mean(result)else result
}

# Compare population, posterior-based and empirical means, scores and conditioning.
pmbfa_mean_diagnostics <- function(fit,x,moments=NULL) {
  md <- fit$params$moddat;gm <- fit$params$gmp
  if(!identical(as.integer(dim(x)),as.integer(c(md$n,md$p,md$N))))
    stop("x dimensions do not match the fit.")
  if(is.null(moments)){
    if(!exists("pmbfa_posterior_moments",mode="function"))
      stop("Compile the updated helpers.cpp providing pmbfa_posterior_moments.")
    moments <- pmbfa_posterior_moments(x,fit$params,fit$control)
  }
  if(!is.null(moments$status) && !identical(moments$status,"ok"))
    stop("Posterior moment evaluation failed: ",moments$message,call.=FALSE)
  z <- moments$z;a <- moments$a;b <- moments$b
  expected_dim <- c(as.integer(md$N),as.integer(md$G))
  if(!all(vapply(list(z,a,b),function(v)
    is.numeric(v) && identical(dim(v),expected_dim) && all(is.finite(v)),logical(1))) ||
    any(z<0) || any(abs(rowSums(z)-1)>1e-7) || any(a<=0) || any(b<=0) ||
    any(colSums(z)<=0))
    stop("Posterior moments must be finite N by G matrices with positive component sizes.",call.=FALSE)
  cluster <- max.col(z,ties.method="first")
  wsum <- function(weights)matrix(drop(matrix(x,nrow=md$n*md$p)%*%weights),md$n,md$p)
  norm <- function(value)sqrt(mean(value^2))
  floor <- if(is.null(fit$control$variance.floor)).1 else fit$control$variance.floor
  skewed <- md$family!="Normal"
  shared_A <- skewed&&substr(md$name,9,9)=="C"
  pooled_A_gradient <- matrix(0,md$n,md$p)
  pooled_M_gradient <- matrix(0,md$n,md$p)
  table <- vector("list",md$G);matrices <- vector("list",md$G)
  for(g in seq_len(md$G)){
    arg <- fit$params$args[[g]];M <- arg$M
    A <- if(!skewed || is.null(arg$A))matrix(0,md$n,md$p)else arg$A
    ng <- sum(z[,g]);hard_n <- sum(cluster==g)
    sa <- sum(z[,g]*a[,g]);sb <- sum(z[,g]*b[,g]);abar <- sa/ng;bbar <- sb/ng
    empirical <- wsum(z[,g])/ng
    hard_empirical <- if(hard_n)wsum(as.numeric(cluster==g))/hard_n else matrix(NA_real_,md$n,md$p)
    conditional <- M+abar*A
    prior_w <- .pmbfa_prior_mixing_mean(md$family,arg)
    population <- M+prior_w*A
    score_M <- wsum(z[,g]*b[,g])-sb*M-ng*A
    score_A <- wsum(z[,g])-ng*M-sa*A
    Lambda <- matrix(gm$Lambdak[g,,],md$n,md$q)
    Delta <- matrix(gm$Deltak[g,,],md$p,md$r)
    S <- tcrossprod(Lambda)+diag(gm$Sigmak[,g],md$n)
    P <- tcrossprod(Delta)+diag(gm$Psik[,g],md$p)
    Sinv <- solve(S);Pinv <- solve(P)
    grad_A <- Sinv%*%score_A%*%Pinv
    grad_M <- Sinv%*%score_M%*%Pinv
    pooled_A_gradient <- pooled_A_gradient+grad_A
    pooled_M_gradient <- pooled_M_gradient+grad_M
    Seigen <- eigen(S,symmetric=TRUE,only.values=TRUE)$values
    Peigen <- eigen(P,symmetric=TRUE,only.values=TRUE)$values
    decomposition <- (prior_w-abar)*A-score_A/ng
    stopifnot(norm(population-empirical-decomposition)<1e-7*
                max(1,norm(population),norm(empirical)))
    table[[g]] <- data.frame(
      component=g,pi=gm$pik[g],posterior_n=ng,hard_n=hard_n,
      population_mean=mean(population),posterior_empirical_mean=mean(empirical),
      hard_empirical_mean=mean(hard_empirical),conditional_mean=mean(conditional),
      M_mean=mean(M),A_mean=mean(A),A_rms=norm(A),
      prior_EW=prior_w,posterior_EW=abar,posterior_E_inverseW=bbar,
      mixing_mean_gap=prior_w-abar,
      population_vs_posterior_rmse=norm(population-empirical),
      population_vs_hard_rmse=norm(population-hard_empirical),
      conditional_vs_posterior_rmse=norm(conditional-empirical),
      prior_conditional_A_gap_rmse=norm((prior_w-abar)*A),
      location_score_rmse=norm(score_M/ng),skew_score_rmse=if(skewed)norm(score_A/ng)else NA_real_,
      location_gradient_rmse=norm(grad_M/ng),skew_gradient_rmse=if(skewed)norm(grad_A/ng)else NA_real_,
      nu=if(md$family=="ST")arg$nu else NA_real_,
      determinant_factor=abar*bbar-1,
      sig_min=min(gm$Sigmak[,g]),sig_mean=mean(gm$Sigmak[,g]),
      psi_min=min(gm$Psik[,g]),psi_mean=mean(gm$Psik[,g]),
      sig_floor_fraction=mean(gm$Sigmak[,g]<=floor*(1+1e-8)),
      psi_floor_fraction=mean(gm$Psik[,g]<=floor*(1+1e-8)),
      row_scale_condition=max(Seigen)/min(Seigen),
      column_scale_condition=max(Peigen)/min(Peigen),
      shared_skewness=shared_A)
    matrices[[g]] <- list(population=population,posterior_empirical=empirical,
      hard_empirical=hard_empirical,conditional=conditional,
      prior_conditional_A_gap=(prior_w-abar)*A,
      location_score=score_M/ng,skew_score=if(skewed)score_A/ng else NULL,
      location_gradient=grad_M/ng,skew_gradient=if(skewed)grad_A/ng else NULL)
  }
  list(components=do.call(rbind,table),matrices=matrices,
    pooled_skew_gradient_rmse=if(skewed)norm(pooled_A_gradient/md$N)else NA_real_,
    pooled_location_gradient_rmse=norm(pooled_M_gradient/md$N),
    posterior_loglik=moments$ll,fit_loglik=tail(fit$loglik,1L),
    interpretation=if(!skewed)
      "At a stationary Normal fit, each component location equals its posterior-weighted empirical mean; there is no skewness parameter."
      else if(shared_A)
      "Common A has a pooled precision-weighted score equation; componentwise conditional means need not match posterior empirical means."
      else "At a stationary fit with separate A, posterior empirical means equal M + posterior E[W] A; population E[W] is not constrained to equal posterior E[W].")
}
