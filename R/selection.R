# Internal organization: see inst/DEVELOPMENT.md.
# Execute forward, backward or both staged constraint searches.
.pmbfa_adaptive <- function(fit_model, G, family, method, criterion,
                             screening="none", npar_model=NULL,
                             has_fit=function(name)FALSE) {
  normal <- family == "Normal"
  forward <- list(concentration=11L, variance.sharing=c(2L,6L),
                  isotropy=c(3L,7L), loadings=c(1L,5L), skewness=9L)
  backward <- list(skewness=9L, loadings=c(1L,5L), isotropy=c(3L,7L),
                   variance.sharing=c(2L,6L), concentration=11L)
  if (normal) {
    forward <- forward[c("variance.sharing", "isotropy", "loadings")]
    backward <- backward[c("loadings", "isotropy", "variance.sharing")]
  }
  directions <- if (method == "both") c("forward", "backward") else method
  paths <- vector("list", length(directions)); names(paths) <- directions
  results <- vector("list", length(directions)); records <- list()
  for (direction in directions) {
    current <- if (direction == "forward")
      if (normal) "CCC-CCC" else "CCC-CCC-C-C" else
        if (normal) "UUU-UUU" else "UUU-UUU-U-U"
    current <- .pmbfa_canonical(current, G)
    stages <- if (direction == "forward") forward else backward
    path <- vector("list", length(stages)); names(path) <- names(stages)
    for (stage in names(stages)) {
      positions <- stages[[stage]]
      candidates <- .pmbfa_stage_names(current,positions,G)
      pair <- direction == "backward" && screening == "nested" &&
        length(positions) == 2L && length(candidates) == 4L && G > 1L
      if (pair) {
        # Loose model and the two direct parents precede their common child.
        constraints <- vapply(candidates,function(name)
          sum(substring(name,positions,positions) == "C"),integer(1))
        candidates <- candidates[order(constraints)]
      }
      fits <- list(); rows <- vector("list",length(candidates))
      for (i in seq_along(candidates)) {
        name <- candidates[i]
        screen <- list(screened=FALSE,lower.bound=NA_real_,
                       parents=character(),reason="")
        is_joint <- pair && all(substring(name,positions,positions) == "C")
        if (is_joint && !has_fit(name)) {
          parent_names <- vapply(positions,function(position) {
            parent <- name; substr(parent,position,position) <- "U"; parent
          },character(1))
          if (all(parent_names %in% names(fits))) {
            incumbent <- .pmbfa_best(fits,criterion)$score
            screen <- .pmbfa_nested_screen(name,fits[parent_names],incumbent,
                                            criterion,npar_model)
          }
        }
        if (screen$screened) {
          record <- data.frame(direction=direction,stage=stage,model=name,
            lower.bound=screen$lower.bound,
            incumbent.score=incumbent,
            parents=paste(screen$parents,collapse=", "),
            reason=screen$reason,stringsAsFactors=FALSE)
          records[[length(records)+1L]] <- record
          rows[[i]] <- data.frame(model=name,score=NA_real_,converged=NA,
            selected=FALSE,screened=TRUE,lower.bound=screen$lower.bound,
            parents=record$parents,status="screened",stringsAsFactors=FALSE)
        } else {
          fit <- fit_model(name); fits[[name]] <- fit
          rows[[i]] <- data.frame(model=name,
            score=if(isTRUE(fit$converged))fit[[tolower(criterion)]]else Inf,
            converged=isTRUE(fit$converged),selected=FALSE,screened=FALSE,
            lower.bound=screen$lower.bound,
            parents=paste(screen$parents,collapse=", "),
            status=fit$status,stringsAsFactors=FALSE)
        }
      }
      best <- .pmbfa_best(fits, criterion)
      current <- best$params$moddat$name
      path[[stage]] <- do.call(rbind,rows)
      path[[stage]]$selected <- path[[stage]]$model == current
    }
    paths[[direction]] <- path
    results[[match(direction, directions)]] <- best
  }
  best <- .pmbfa_best(results, criterion)
  skipped <- if(length(records))do.call(rbind,records)else data.frame(
    direction=character(),stage=character(),model=character(),
    lower.bound=numeric(),incumbent.score=numeric(),parents=character(),
    reason=character(),stringsAsFactors=FALSE)
  best$selection <- list(method=method,paths=paths,screening=screening,
    screened.models=skipped,
    screening.caveat=if(screening == "nested")
      "Nested screening assumes parent optima; finite multistart likelihoods are not certified global upper bounds." else NULL)
  best
}

# Reuse one initialization bank while fitting and caching constraint models.
.pmbfa_run_tuple <- function(x, G, q, r, family, labels, modelSubset,
                             adaptive, adaptive.method, criterion, seed,
                             inits, control, screening="none") {
  tuple_seed <- .pmbfa_tuple_seed(seed, G, q, r, family)
  control$seed <- tuple_seed
  control$criterion <- criterion
  initial_name <- if (family == "Normal") "UUU-UUU" else "UUU-UUU-U-U"
  bank <- if (is.null(inits)) tryCatch(
    .pmbfa_with_seed(tuple_seed,
      pmbfa_initializations(x, G, q, r, family, aargs=control,
                            labels=as.numeric(labels))),
    error=function(error) error) else
      if (!is.null(inits$moddat)) list(inits) else inits
  if (inherits(bank, "error") || !is.list(bank) || !length(bank)) {
    message <- if (inherits(bank, "error")) conditionMessage(bank) else
      "Initialization produced no usable starting values."
    failed <- .pmbfa_failure(x, G, q, r, family, initial_name, message)
    failed$criterion <- criterion
    return(failed)
  }
  for (start in bank) {
    if (!is.list(start) || !all(c("moddat", "gmp", "args") %in% names(start)))
      stop("Each initialization must contain moddat, gmp, and args.", call.=FALSE)
    expected <- c(G=G, q=q, r=r)
    for (key in names(expected))
      if (!is.null(start$moddat[[key]]) &&
          !identical(as.integer(start$moddat[[key]]), as.integer(expected[[key]])))
        stop("Initialization dimensions do not match G, q, and r.", call.=FALSE)
  }
  cache <- new.env(parent=emptyenv())
  fit_model <- function(name) {
    name <- .pmbfa_canonical(name, G)
    if (exists(name, envir=cache, inherits=FALSE)) return(get(name, envir=cache))
    fits <- lapply(bank, function(start) {
      start$moddat$name <- name
      start$moddat$family <- family
      start$moddat$class <- as.numeric(labels)
      start$moddat$icl <- criterion == "ICL"
      fitter <- if (family == "Normal") pmbfaEM else pmsbfaEM
      tryCatch(fitter(x, start, control), error=function(error)
        .pmbfa_failure(x,G,q,r,family,name,conditionMessage(error), start,
                       "numerical_failure"))
    })
    converged <- vapply(fits, function(fit) isTRUE(fit$converged), logical(1))
    ll <- vapply(fits, .pmbfa_final_ll, numeric(1))
    index <- if (any(converged)) which.max(replace(ll, !converged, -Inf)) else
      which.max(ll)
    best <- fits[[index]]
    best$params$moddat$name <- name
    best$criterion <- criterion
    best$score <- if (isTRUE(best$converged)) best[[tolower(criterion)]] else Inf
    best$starts <- data.frame(start=seq_along(fits), loglik=ll,
                              converged=converged)
    assign(name, best, envir=cache)
    if (isTRUE(control$verbose))
      cat(sprintf("G=%d q=%d r=%d %s %s=%.2f (%s)\n", G,q,r,name,
                   criterion,best$score,best$status))
    best
  }
  if (!is.null(modelSubset) || !adaptive) {
    models <- if (is.null(modelSubset)) .pmbfa_models(family) else modelSubset
    models <- unique(vapply(models, .pmbfa_canonical, character(1), G=G))
    best <- .pmbfa_best(lapply(models, fit_model), criterion)
    best$selection <- list(method="exhaustive", models=models)
  } else best <- .pmbfa_adaptive(fit_model,G,family,adaptive.method,criterion,
    screening=screening,
    npar_model=function(name)pmbfa_model_npar(dim(x)[1],dim(x)[2],G,q,r,
                                             name,family),
    has_fit=function(name)exists(name,envir=cache,inherits=FALSE))
  names <- ls(cache, all.names=TRUE)
  best$model.criteria <- data.frame(model=names,
    bic=vapply(names, function(name) get(name,cache)$bic, numeric(1)),
    icl=vapply(names, function(name) get(name,cache)$icl, numeric(1)),
    converged=vapply(names, function(name)
      isTRUE(get(name,cache)$converged), logical(1)), stringsAsFactors=FALSE)
  best$model.criteria$score <- ifelse(best$model.criteria$converged,
    best$model.criteria[[tolower(criterion)]], Inf)
  best$initialization <- list(seed=tuple_seed, nstarts=length(bank),
                              shared.across.models=TRUE)
  map <- attr(labels, "label.map")
  if (!is.null(map)) best$label.map <- map[seq_len(G)]
  best
}

