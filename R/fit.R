# Internal organization: see inst/DEVELOPMENT.md.
# Construct a consistent diagnostic fit after an initialization or numerical failure.
.pmbfa_failure <- function(x, G, q, r, family, name, message,
                           params=NULL, status="initialization_failed") {
  if (is.null(params)) params <- list(
    moddat=list(n=dim(x)[1], p=dim(x)[2], N=dim(x)[3], G=G, q=q, r=r,
                family=family, name=name), gmp=list(), args=list())
  structure(list(params=params, loglik=numeric(), iter=0L,
                 posterior=matrix(NA_real_, dim(x)[3], G),
                 cluster=rep(NA_integer_, dim(x)[3]), bic=Inf, icl=Inf, BIC=Inf, ICL=Inf,
                 criterion="BIC", score=Inf, npar=NA_integer_,
                 converged=FALSE, status=status, message=message),
            class=if (family == "Normal") "pmbfa" else c("pmsbfa", "pmbfa"))
}

# Retrieve the last finite accepted likelihood, or minus infinity if none exists.
.pmbfa_final_ll <- function(fit) {
  values <- fit$loglik[is.finite(fit$loglik)]
  if (length(values)) tail(values, 1L) else -Inf
}

# Choose by the requested criterion while excluding nonconverged fits.
.pmbfa_best <- function(fits, criterion) {
  values <- vapply(fits, function(fit) {
    value <- fit[[tolower(criterion)]]
    if (isTRUE(fit$converged) && length(value) == 1L && is.finite(value))
      value else Inf
  }, numeric(1))
  # If all fits failed, retain the most informative last accepted state.
  index <- if (any(is.finite(values))) which.min(values) else
    which.max(vapply(fits, .pmbfa_final_ll, numeric(1)))
  out <- fits[[index]]
  out$criterion <- criterion
  out$score <- unname(values[index])
  out
}

# Validate a public fitting request and orchestrate its component/factor grid.
.pmbfa_fit <- function(x, G, q, r, class, modelSubset, family, seed, adaptive,
                       adaptive.method, criterion, verbose, parallelize,
                       ncores, inits, control, dots, screening="none",
                       report=TRUE) {
  for (alias in c("rG", "rq", "rr")) {
    if (!is.null(dots[[alias]])) {
      if (alias == "rG") G <- dots[[alias]]
      if (alias == "rq") q <- dots[[alias]]
      if (alias == "rr") r <- dots[[alias]]
      dots[[alias]] <- NULL
    }
  }
  if (!is.null(dots$icl)) criterion <- if (isTRUE(dots$icl)) "ICL" else "BIC"
  if (!is.null(dots[["track"]])) verbose <- isTRUE(dots[["track"]])
  if (!is.null(dots$aitken)) control$aitken <- dots$aitken
  recognized <- c("icl", "track", "quiet", "aitken", "max.iter", "tol",
                   "mag", "stopping", "ncandid", "nstart", "nrep",
                   "init.iter", "min.component", "variance.floor", "nu.min",
                   "shape.min", "verbose", "track_likelihood", "min.iter", "backtracks",
                   "pcg.max.iter", "pcg.tol", "pcg.preconditioner", "nu.max", "init_method", "profile")
  unknown <- setdiff(names(dots), recognized)
  if (length(unknown)) stop("Unused argument(s): ",
    paste(unknown, collapse=", "), call.=FALSE)
  if (!is.numeric(x) || length(dim(x)) != 3L || any(!is.finite(x)) ||
      any(dim(x) < 1L))
    stop("x must be a finite numeric array with dimensions c(n,p,N).",
         call.=FALSE)
  storage.mode(x) <- "double"
  n <- dim(x)[1]; p <- dim(x)[2]; N <- dim(x)[3]
  if (N <= 10L) stop("x must contain more than 10 observations.",call.=FALSE)
  G <- .pmbfa_integer(G,"G", upper=N)
  if (any(N/G < 30)) warning(
    "Some requested component counts have fewer than 30 observations per component on average; fits may be unstable.",
    call.=FALSE)
  q <- .pmbfa_integer(q,"q", upper=.pmbfa_max_factors(n))
  r <- .pmbfa_integer(r,"r", upper=.pmbfa_max_factors(p))
  family <- match.arg(family, c("ST", "GH", "NIG", "VG", "Normal"))
  criterion <- match.arg(toupper(criterion), c("BIC", "ICL"))
  if (is.character(adaptive) && length(adaptive) == 1L) {
    adaptive.method <- adaptive
    adaptive <- TRUE
  }
  if (!is.logical(adaptive) || length(adaptive) != 1L || is.na(adaptive))
    stop("adaptive must be TRUE, FALSE, or a search direction.", call.=FALSE)
  adaptive.method <- match.arg(adaptive.method, c("both", "forward", "backward"))
  screening <- match.arg(screening,c("none","nested"))
  if (!is.null(modelSubset) && (!is.character(modelSubset) ||
      !length(modelSubset) || anyNA(modelSubset) ||
      !all(modelSubset %in% .pmbfa_models(family))))
    stop("modelSubset contains invalid constraint models.", call.=FALSE)
  if (!is.null(seed)) {
    seed <- .pmbfa_integer(seed,"seed",lower=0L,
                           upper=.Machine$integer.max-1L)
    if (length(seed) != 1L) stop("seed must be scalar.", call.=FALSE)
  }
  labels <- .pmbfa_labels(class,N,G)
  control <- .pmbfa_controls(control,dots,verbose,N,n,p)
  tuples <- expand.grid(G=G,q=q,r=r, KEEP.OUT.ATTRS=FALSE)
  if (!is.null(inits) && nrow(tuples) != 1L)
    stop("Supply inits only for a single (G,q,r) combination.", call.=FALSE)
  if (!is.logical(parallelize) || length(parallelize) != 1L ||
      is.na(parallelize)) stop("parallelize must be TRUE or FALSE.", call.=FALSE)
  if (!is.null(ncores)) {
    ncores <- .pmbfa_integer(ncores,"ncores")
    if (length(ncores) != 1L) stop("ncores must be scalar.",call.=FALSE)
  }
  run_parallel <- parallelize && nrow(tuples) > 1L
  runner <- function(i) {
    fit <- .pmbfa_timed_tuple(x,tuples$G[i],tuples$q[i],tuples$r[i],
      family,labels,modelSubset,adaptive,adaptive.method,criterion,seed,inits,control,
      screening=screening)
    if (report && !run_parallel) .pmbfa_report_fit(fit,control$verbose)
    fit
  }
  if (run_parallel) {
    if (is.null(ncores)) ncores <- min(nrow(tuples),
      max(1L, parallel::detectCores(logical=FALSE)-1L, na.rm=TRUE))
    ncores <- min(ncores,nrow(tuples))
    if (.Platform$OS.type == "unix") {
      out <- parallel::mclapply(seq_len(nrow(tuples)),runner,
                                 mc.cores=ncores,mc.set.seed=is.null(seed))
    } else {
      arguments <- list(x=x,family=family,labels=labels,modelSubset=modelSubset,
        adaptive=adaptive,adaptive.method=adaptive.method,criterion=criterion,
        seed=seed,inits=inits,control=control,screening=screening)
      out <- .pmbfa_psock(tuples,arguments,ncores)
    }
  } else out <- lapply(seq_len(nrow(tuples)),runner)
  if (report && run_parallel)
    for (fit in out) .pmbfa_report_fit(fit,control$verbose)
  if (length(out) == 1L) return(out[[1]])
  out <- out[order(vapply(out,function(fit)fit$score,numeric(1)))]
  class(out) <- c("pmbfa_fits","list")
  out
}

