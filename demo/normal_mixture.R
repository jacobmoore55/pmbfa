# Normal PMBFA data: G=2, q=r=2, covariance model CCC-CUC.
# Each component contributes 100 full 8 by 8 matrix observations.
library(pmbfa)
set.seed(2026)


# Simulated Data ----------------------------------------------------------

true_class <- rep(1:2,each=100L)
true_locations <- list(
  matrix(0,8L,8L),
  matrix(1,8L,8L)
  )

# Common rank-two row/column loadings. Their columns are independent.
true_Lambda <- .7*cbind(rep(c(1,-1),4L),rep(c(1,1,-1,-1),2L))
true_Delta <- .6*cbind(seq(-1,1,length.out=8L),rep(c(1,-1,-1,1),2L))
true_sigma <- .8
true_psi <- c(.5,1.2)
# CCC: common row loadings and one common isotropic residual variance.
true_row_scale <- tcrossprod(true_Lambda)+diag(true_sigma,8L)
# CUC: common column loadings, component-specific isotropic residuals.
true_column_scales <- lapply(true_psi,function(psi)
  tcrossprod(true_Delta)+diag(psi,8L))
true_model <- list(family="Normal",G=2L,q=2L,r=2L,model="CCC-CUC",
  Lambda=true_Lambda,Delta=true_Delta,sigma=true_sigma,psi=true_psi,
  row_scale=true_row_scale,column_scales=true_column_scales,M=true_locations)

# If Z has independent standard normal entries, L_row %*% Z %*% U_col
# has matrix-normal row scale Sigma and column scale Psi_g.
row_chol <- t(chol(true_row_scale))
column_chols <- lapply(true_column_scales,chol)
x <- array(0,c(8L,8L,200L))
for(i in seq_len(200L)) {
  g <- true_class[i]
  Z <- matrix(rnorm(64L),8L,8L)
  x[,,i] <- true_locations[[g]]+row_chol%*%Z%*%column_chols[[g]]
}
ordering <- sample.int(200L)
x <- x[,,ordering,drop=FALSE]
true_class <- true_class[ordering]


# Parameter Estimation on Simulated Data ----------------------------------

# Adaptive model selection is used by default to cut down on computation time.
models <- search_pmbfa(x, G = 1:4, verbose = TRUE)

comparison <- do.call(rbind,lapply(models,function(fit)data.frame(
  G=fit$params$moddat$G,q=fit$params$moddat$q,r=fit$params$moddat$r,
  model=fit$params$moddat$name,loglik=as.numeric(logLik(fit)),
  BIC=BIC(fit),score=fit$score,converged=fit$converged,status=fit$status,elapsed=fit$elapsed)))
comparison <- comparison[order(comparison$G),]

# Model selections for PMBFA models with 1 to 4 components (G)
comparison

# Output is ordered by best score (lowest BIC)
selected_model <- models[[1L]]


table(true_class, 
      selected_model$cluster)

selected_model

