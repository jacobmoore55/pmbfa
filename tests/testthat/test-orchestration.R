test_that("shared starts, adaptive stages and factor search retain their contracts", {
env <- clone_engine()
evalq({
  audit <- new.env(parent=emptyenv())
  audit$calls <- list(); audit$initialize <- 0L; audit$fail <- FALSE
  pmbfa_initializations <- function(x,G,q,r,family="ST",aargs=list(),labels=numeric()) {
    audit$initialize <- audit$initialize+1L
    marker <- runif(1)
    lapply(seq_len(aargs$nrep),function(i) list(
      moddat=list(n=dim(x)[1],p=dim(x)[2],N=dim(x)[3],G=G,q=q,r=r,
                   family=family,name="UUU-UUU-U-U"),
      gmp=list(id=i,marker=marker),args=rep(list(list()),G)))
  }
  pmsbfaEM <- function(x,params,aargs=list()) {
    if (audit$fail) stop("deliberate numerical failure")
    m <- params$moddat
    audit$calls[[length(audit$calls)+1L]] <- list(
      model=m$name,id=params$gmp$id,marker=params$gmp$marker,labels=m$class,
      tol=aargs$tol,stopping=aargs$stopping)
    score <- (m$q-2)^2+(m$r-1)^2+sum(strsplit(m$name,"")[[1]] == "U")
    ll <- 10*params$gmp$id
    structure(list(params=params,loglik=c(ll-1,ll),iter=1L,
      bic=score,icl=-score,score=score,npar=10L,criterion="BIC",
      converged=TRUE,status="converged",message="",
      posterior=matrix(1/m$G,m$N,m$G),cluster=rep(1L,m$N)),
      class=if(m$family == "Normal")"pmbfa" else c("pmsbfa","pmbfa"))
  }
  pmbfaEM <- pmsbfaEM
  pmbfa_core_version <- function() "1.1.0"
  pmbfa_model_npar <- function(...)stop("Counting is not expected in unscreened orchestration tests")
  pmbfa_posterior_moments <- function(...)stop("Moments are not expected in orchestration tests")
  x <- array(seq_len(6*7*60)/100,dim=c(6,7,60))
  controls <- list(ncandid=2L,nstart=2L,nrep=2L,init.iter=1L,
                    tol=1e-7,stopping="aitken")
  set.seed(7); before <- .Random.seed
  fit <- pmsbfa(x,2,1,1,modelSubset="UUU-UUU-U-U",control=controls)
  stopifnot(inherits(fit,"pmbfa"),!inherits(fit,"pmbfa_fits"),
             fit$params$gmp$id == 2L,identical(.Random.seed,before),
             fit$iter == 1L,tail(fit$loglik,1L) == 20,
             all(vapply(audit$calls,function(z)z$tol == 1e-7,logical(1))),
             all(vapply(audit$calls,function(z)z$stopping == "aitken",logical(1))))
  audit$calls <- list(); audit$initialize <- 0L
  selected <- pmsbfa(x,2,1,1,adaptive.method="both",control=controls)
  stopifnot(audit$initialize == 1L,
    identical(names(selected$selection$paths$forward),
      c("concentration","variance.sharing","isotropy","loadings","skewness")),
    identical(names(selected$selection$paths$backward),
      c("skewness","loadings","isotropy","variance.sharing","concentration")),
    identical(selected$loglik,c(19,20)))
  starts <- vapply(audit$calls,function(z)z$marker,numeric(1))
  keys <- vapply(audit$calls,function(z)paste(z$model,z$id),character(1))
  stopifnot(length(unique(starts)) == 1L,!anyDuplicated(keys),
    all(selected$selection$paths$forward$concentration$model %in%
      c("CCC-CCC-C-C","CCC-CCC-C-U")))
  normal <- pmbfa(x,2,1,1,adaptive.method="both",control=controls)
  stopifnot(all(nchar(normal$model.criteria$model) == 7L),
    identical(names(normal$selection$paths$forward),
      c("variance.sharing","isotropy","loadings")))
  grid <- pmbfa(x,G=1:2,q=1:2,r=1,modelSubset="CCC-CCC",control=controls)
  stopifnot(inherits(grid,"pmbfa_fits"),length(grid) == 4L,
    !is.unsorted(vapply(grid,function(z)z$score,numeric(1))))
  aliased <- pmbfa(x,rG=2,rq=1,rr=1,modelSubset="CCC-CCC",control=controls)
  stopifnot(aliased$params$moddat$G == 2L)
  labels <- c(rep("alpha",20),rep("beta",20),rep(NA,20))
  labelled <- pmbfa(x,2,1,1,class=labels,modelSubset="CCC-CCC",control=controls)
  stopifnot(identical(labelled$params$moddat$class,c(rep(1,20),rep(2,20),rep(0,20))))
  full <- pmbfa(x,2,1,1,class=rep(1:2,each=30),
                modelSubset="CCC-CCC",control=controls)
  stopifnot(inherits(full,"pmbfa"))
  icl <- pmbfa(x,2,1,1,modelSubset=c("CCC-CCC","UUU-UUU"),
               criterion="ICL",control=controls)
  stopifnot(icl$params$moddat$name == "UUU-UUU",icl$criterion == "ICL")
  search <- search_pmbfa(x,G=2,max_q=2,max_r=1,modelSubset="CCC-CCC",
                          control=controls)
  stopifnot(search$params$moddat$q == 2L,search$params$moddat$r == 1L)
  audit$fail <- TRUE
  failed <- pmbfa(x,2,1,1,modelSubset="CCC-CCC",control=controls)
  stopifnot(inherits(failed,"pmbfa"),!failed$converged,is.infinite(failed$score),
             failed$status == "numerical_failure",length(failed$loglik) == 0L,
             !is.null(failed$params))
  invalid <- try(pmbfa(x,2,1,1,control=list(tol=-1)),silent=TRUE)
  stopifnot(inherits(invalid,"try-error"))
},env)
expect_true(TRUE)
})
