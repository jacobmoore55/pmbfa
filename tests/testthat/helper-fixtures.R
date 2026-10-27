native <- function(name) getFromNamespace(name, "pmbfa")
native_control <- function(...) modifyList(list(max.iter=1L,min.iter=1L,
  min.component=0,variance.floor=.05,tol=1e-10,verbose=FALSE),list(...))
numerical_fixture <- function(family = 'Normal', n = 4L, p = 4L, N = 24L, G = 2L) {
  set.seed(871)
  labels <- rep(seq_len(G), length.out = N)
  x <- array(rnorm(n*p*N, sd = .7), c(n,p,N))
  for (i in seq_len(N)) x[,,i] <- x[,,i] + 1.8*(labels[i] - (G+1)/2)
  args <- lapply(seq_len(G), function(g) {
    z <- list(M = apply(x[,,labels == g,drop=FALSE], 1:2, mean))
    if (family != 'Normal') z$A <- matrix(.08 + .03*g, n, p)
    if (family == 'ST') z$nu <- 6 + g
    if (family == 'GH') { z$omega <- 2 + g/2; z$lambda <- .6 + g/2 }
    if (family == 'NIG') z$gam <- .8 + g/3
    if (family == 'VG') z$gamma <- 3 + g/2
    z
  })
  params <- list(moddat = list(n=n,p=p,N=N,G=G,q=1L,r=1L,family=family,
    name=if(family == 'Normal') 'UUU-UUU' else 'UUU-UUU-U-U',class=rep(0,N)),
    gmp = list(Sigmak = matrix(seq(.7,1.3,length.out=n*G),n,G),
      Psik = matrix(seq(.9,1.5,length.out=p*G),p,G),
      Lambdak = array(runif(G*n,.1,.4),c(G,n,1L)),
      Deltak = array(runif(G*p,.1,.4),c(G,p,1L)), pik=rep(1/G,G)), args=args)
  list(x=x,params=params,labels=labels)
}

quick_control <- function(...) modifyList(list(ncandid=4L,nstart=2L,nrep=1L,
  init.iter=2L,max.iter=100L,tol=1e-4,min.component=1),list(...))
clone_engine <- function() {
  ns <- asNamespace("pmbfa"); env <- new.env(parent=ns)
  for(name in ls(ns,all.names=TRUE)) {
    value <- get(name,ns,inherits=FALSE)
    if(is.function(value) && !startsWith(name,"_pmbfa_")) {
      environment(value) <- env
      assign(name,value,env)
    }
  }
  env
}
