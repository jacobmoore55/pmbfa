# Internal organization: see inst/DEVELOPMENT.md.
# Search factor dimensions, cache comparisons and rank results for each G.
.pmbfa_factor_search <- function(x,G,max_q,max_r,family,search.method,dots) {
  if (!is.numeric(x) || length(dim(x)) != 3L)
    stop("x must be a numeric three-dimensional array.",call.=FALSE)
  if ("family" %in% names(dots))
    stop("Use search_pmsbfa() to choose a distribution family.",call.=FALSE)
  G <- .pmbfa_integer(G,"G",upper=dim(x)[3])
  if (is.null(max_q)) max_q <- .pmbfa_max_factors(dim(x)[1])
  if (is.null(max_r)) max_r <- .pmbfa_max_factors(dim(x)[2])
  max_q <- .pmbfa_integer(max_q,"max_q",upper=.pmbfa_max_factors(dim(x)[1]))
  max_r <- .pmbfa_integer(max_r,"max_r",upper=.pmbfa_max_factors(dim(x)[2]))
  if (length(max_q) != 1L || length(max_r) != 1L)
    stop("max_q and max_r must be scalar.",call.=FALSE)
  if (any(c("q","r","rq","rr","rG","G") %in% names(dots)))
    stop("Use G, max_q, and max_r to specify the factor search.",call.=FALSE)
  search.method <- match.arg(search.method,c("adaptive","grid"))
  criterion <- if (isTRUE(dots$icl)) "ICL" else
    if (!is.null(dots$criterion)) toupper(dots$criterion[1]) else "BIC"
  warned_small_components <- FALSE
  results <- lapply(G,function(g) {
    started_g <- proc.time()[["elapsed"]]
    cache <- new.env(parent=emptyenv())
    evaluate <- function(q,r) {
      key <- paste(q,r,sep=",")
      if (!exists(key,cache,inherits=FALSE)) {
        arguments <- c(list(x=x,G=g,q=q,r=r),dots)
        if (family != "Normal") arguments$family <- family
        fit <- withCallingHandlers(
          do.call(.pmbfa_search_fit,arguments),
          warning=function(condition) {
            if (startsWith(conditionMessage(condition),
              "Some requested component counts have fewer than 30")) {
              if (warned_small_components) invokeRestart("muffleWarning")
              warned_small_components <<- TRUE
            }
          })
        assign(key,fit,cache)
      }
      get(key,cache)
    }
    if (search.method == "grid") {
      for (q in seq_len(max_q)) for (r in seq_len(max_r)) evaluate(q,r)
    } else {
      # Coarse-to-fine coordinate search, followed by unit-neighbor refinement.
      # This is a heuristic and does not claim the global grid optimum.
      q <- max(1L,ceiling(max_q/2)); r <- max(1L,ceiling(max_r/2))
      step_q <- max(1L,ceiling(max_q/2)); step_r <- max(1L,ceiling(max_r/2))
      current <- evaluate(q,r)
      repeat {
        pairs <- unique(rbind(c(q,r),c(max(1L,q-step_q),r),
          c(min(max_q,q+step_q),r),c(q,max(1L,r-step_r)),
          c(q,min(max_r,r+step_r))))
        fits <- lapply(seq_len(nrow(pairs)),function(i)
          evaluate(pairs[i,1],pairs[i,2]))
        best <- .pmbfa_best(fits,criterion)
        if (best$score < current$score) {
          q <- best$params$moddat$q; r <- best$params$moddat$r
          current <- best
        } else if (step_q == 1L && step_r == 1L) break else {
          step_q <- max(1L,ceiling(step_q/2)); step_r <- max(1L,ceiling(step_r/2))
        }
      }
    }
    fits <- mget(ls(cache),cache)
    best <- .pmbfa_best(fits,criterion)
    best$factor.search <- data.frame(
      G=g,q=vapply(fits,function(fit)as.integer(fit$params$moddat$q),integer(1)),
      r=vapply(fits,function(fit)as.integer(fit$params$moddat$r),integer(1)),
      score=vapply(fits,function(fit)fit$score,numeric(1)),
      converged=vapply(fits,function(fit)isTRUE(fit$converged),logical(1)))
    best$factor.search.method <- search.method
    best$elapsed <- .pmbfa_elapsed(started_g)
    .pmbfa_report_fit(best)
    best
  })
  if (length(results) == 1L) return(results[[1]])
  results <- results[order(vapply(results,function(fit)fit$score,numeric(1)))]
  class(results) <- c("pmbfa_fits","list")
  results
}

# Search skewed mixture component counts and factor dimensions.
search_pmsbfa <- function(x,G,max_q=NULL,max_r=NULL,family="ST",
                           search.method=c("adaptive","grid"),...) {
  requested_call <- match.call()
  started <- proc.time()[["elapsed"]]
  family <- match.arg(family,c("ST","GH","NIG","VG","Normal"))
  result <- .pmbfa_factor_search(x,G,max_q,max_r,family,search.method[1L],list(...))
  .pmbfa_record_result(result,requested_call,.pmbfa_elapsed(started))
}

# Search normal mixture component counts and factor dimensions.
search_pmbfa <- function(x,G,max_q=NULL,max_r=NULL,
                          search.method=c("adaptive","grid"),...) {
  requested_call <- match.call()
  started <- proc.time()[["elapsed"]]
  result <- .pmbfa_factor_search(x,G,max_q,max_r,"Normal",search.method[1L],list(...))
  .pmbfa_record_result(result,requested_call,.pmbfa_elapsed(started))
}


# Search candidates use the fitting engine without recording expanded do.call data.
# Fit an internal search candidate without capturing expanded literal-data calls.
.pmbfa_search_fit <- function(x,G,q,r,class=NULL,modelSubset=NULL,
    family="Normal",seed=123456L,adaptive=TRUE,
    adaptive.method=c("both","forward","backward"),criterion=c("BIC","ICL"),
    verbose=FALSE,parallelize=FALSE,ncores=NULL,inits=NULL,control=list(),
    screening=c("none","nested"),...) {
  .pmbfa_fit(x,G,q,r,class,modelSubset,family,seed,adaptive,
    adaptive.method[1L],criterion[1L],verbose,parallelize,ncores,inits,control,
    list(...),screening=screening[1L],report=FALSE)
}
