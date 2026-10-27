dense_normal <- function(x, params, cycles=2L, variance_floor=.1) {
  dims <- dim(x); n <- dims[1L]; p <- dims[2L]; N <- dims[3L]
  G <- params$moddat$G; q <- params$moddat$q; r <- params$moddat$r
  name <- params$moddat$name
  load <- function(side,g) {
    if(side==1L) matrix(params$gmp$Lambdak[g,,],n,q) else
      matrix(params$gmp$Deltak[g,,],p,r)
  }
  diagonal <- function(side,g) {
    if(side==1L) params$gmp$Sigmak[,g] else params$gmp$Psik[,g]
  }
  covariance <- function(side,g) tcrossprod(load(side,g)) + diag(diagonal(side,g))
  # Project each supplied fresh start before evaluating its first likelihood.
  for(side in 1:2) {
    offset <- if(side==1L)0L else 4L
    a <- if(side==1L)params$gmp$Lambdak else params$gmp$Deltak
    d <- if(side==1L)params$gmp$Sigmak else params$gmp$Psik
    if(substr(name,offset+1L,offset+1L)=="C") {
      common <- apply(a,c(2L,3L),mean)
      for(g in seq_len(G))a[g,,] <- common
    }
    d[] <- pmax(variance_floor,d)
    if(substr(name,offset+2L,offset+2L)=="C")d[] <- rowMeans(d)
    if(substr(name,offset+3L,offset+3L)=="C")
      for(g in seq_len(G))d[,g] <- mean(d[,g])
    if(side==1L) {
      params$gmp$Lambdak <- a;params$gmp$Sigmak <- d
    }else {
      params$gmp$Deltak <- a;params$gmp$Psik <- d
    }
  }
  expectation <- function() {
    logs <- matrix(0,N,G)
    for(g in seq_len(G)) {
      S <- covariance(1L,g);P <- covariance(2L,g)
      IS <- solve(S);IP <- solve(P)
      constant <- log(params$gmp$pik[g])-n*p/2*log(2*pi)-
        p/2*as.numeric(determinant(S,TRUE)$modulus)-
        n/2*as.numeric(determinant(P,TRUE)$modulus)
      for(i in seq_len(N)) {
        R <- x[,,i]-params$args[[g]]$M
        logs[i,g] <- constant-.5*sum(R*(IS%*%R%*%IP))
      }
    }
    m <- apply(logs,1L,max)
    observation_ll <- m+log(rowSums(exp(logs-m)))
    list(z=exp(logs-observation_ll),ll=sum(observation_ll))
  }
  update_scale <- function(side,z) {
    offset <- if(side==1L)0L else 4L
    d <- if(side==1L)n else p;k <- if(side==1L)q else r
    opposite_d <- if(side==1L)p else n
    old_diagonal <- lapply(seq_len(G),function(g)diagonal(side,g))
    C <- T <- Q <- vector("list",G)
    for(g in seq_len(G)) {
      L <- load(side,g)
      V <- solve(diag(k)+crossprod(L,L/old_diagonal[[g]]))
      B <- V%*%t(L/old_diagonal[[g]])
      opposite_precision <- solve(covariance(3L-side,g))
      C[[g]] <- matrix(0,d,d)
      for(i in seq_len(N)) {
        R <- x[,,i]-params$args[[g]]$M
        if(side==2L)R <- t(R)
        C[[g]] <- C[[g]]+z[i,g]*R%*%opposite_precision%*%t(R)
      }
      T[[g]] <- C[[g]]%*%t(B)
      Q[[g]] <- sum(z[,g])*opposite_d*V+B%*%C[[g]]%*%t(B)
    }
    fresh_loading <- vector("list",G)
    if(substr(name,offset+1L,offset+1L)=="U") {
      for(g in seq_len(G))fresh_loading[[g]] <- T[[g]]%*%solve(Q[[g]])
    }else {
      common <- matrix(0,d,k)
      for(j in seq_len(d)) {
        rhs <- numeric(k);lhs <- matrix(0,k,k)
        for(g in seq_len(G)) {
          rhs <- rhs+T[[g]][j,]/old_diagonal[[g]][j]
          lhs <- lhs+Q[[g]]/old_diagonal[[g]][j]
        }
        common[j,] <- solve(lhs,rhs)
      }
      for(g in seq_len(G))fresh_loading[[g]] <- common
    }
    residual <- matrix(0,d,G)
    for(g in seq_len(G)) {
      L <- fresh_loading[[g]]
      full_residual <- C[[g]]-L%*%t(T[[g]])-T[[g]]%*%t(L)+L%*%Q[[g]]%*%t(L)
      residual[,g] <- diag(full_residual)
    }
    if(substr(name,offset+2L,offset+2L)=="C") {
      common <- rowSums(residual)/(N*opposite_d)
      for(g in seq_len(G))residual[,g] <- common
    }else {
      for(g in seq_len(G))residual[,g] <- residual[,g]/(sum(z[,g])*opposite_d)
    }
    if(substr(name,offset+3L,offset+3L)=="C")
      for(g in seq_len(G))residual[,g] <- mean(residual[,g])
    residual[] <- pmax(variance_floor,residual)
    if(side==1L) {
      for(g in seq_len(G))params$gmp$Lambdak[g,,] <<- fresh_loading[[g]]
      params$gmp$Sigmak <<- residual
    }else {
      for(g in seq_len(G))params$gmp$Deltak[g,,] <<- fresh_loading[[g]]
      params$gmp$Psik <<- residual
    }
  }
  state <- expectation();curve <- state$ll
  for(iter in seq_len(cycles)) {
    z <- state$z;ng <- colSums(z)
    params$gmp$pik <- ng/N
    for(g in seq_len(G)) {
      center <- matrix(0,n,p)
      for(i in seq_len(N))center <- center+z[i,g]*x[,,i]
      params$args[[g]]$M <- center/ng[g]
    }
    update_scale(1L,z)
    update_scale(2L,z)
    state <- expectation();curve <- c(curve,state$ll)
  }
  list(params=params,loglik=curve,posterior=state$z)
}

