# Internal organization: see inst/DEVELOPMENT.md.
# Validate and normalize positive or bounded integer requests.
.pmbfa_integer <- function(x, name, lower=1L, upper=.Machine$integer.max) {
  if (!is.numeric(x) || !length(x) || anyNA(x) ||
      any(!is.finite(x)) || any(x != floor(x)) ||
      any(x < lower | x > upper))
    stop(name, " must contain integer values in [", lower, ", ", upper,
         "].", call.=FALSE)
  unique(as.integer(x))
}

# Find the conservative identifiable diagonal-factor rank limit.
.pmbfa_max_factors <- function(d) {
  candidates <- seq_len(max(1L, d - 1L))
  candidates <- candidates[(d - candidates)^2 > d + candidates]
  if (length(candidates)) max(candidates) else 0L
}

# Map known labels and preserve unknown observations.
.pmbfa_labels <- function(labels, N, G) {
  if (is.null(labels)) return(rep(0, N))
  if (!is.atomic(labels) || length(labels) != N)
    stop("class must have one label per observation.", call.=FALSE)
  unknown <- is.na(labels) | as.character(labels) == "0"
  known <- unique(labels[!unknown])
  if (length(known) > min(G))
    stop("The number of known classes exceeds G.", call.=FALSE)
  out <- rep(0, N)
  if (is.numeric(labels) && all(labels[!unknown] == floor(labels[!unknown])) &&
      all(labels[!unknown] >= 1 & labels[!unknown] <= min(G))) {
    out[!unknown] <- labels[!unknown]
    map <- seq_len(max(G))
  } else {
    known <- sort(known)
    out[!unknown] <- match(labels[!unknown], known)
    map <- c(as.character(known),
             as.character(seq.int(length(known)+1L, max(G))))
    map <- map[seq_len(max(G))]
  }
  attr(out, "label.map") <- map
  out
}

# Merge, validate and normalize fitting controls and legacy aliases.
.pmbfa_controls <- function(control, dots, verbose, N, n, p) {
  defaults <- list(max.iter=10000L, min.iter=2L, tol=NULL,
                   mag=log10(sqrt(N)*n*p)+3,
                   stopping="increment", ncandid=200L, nstart=20L,
                   nrep=1L, init.iter=5L, verbose=verbose,
                   min.component=6, variance.floor=0.1,
                   nu.min=4.1, nu.max=1e12, shape.min=0.5,
                   backtracks=12L,pcg.max.iter=0L,pcg.tol=1e-9,
                   init_method="nearest",criterion="BIC",profile=FALSE,
                   pcg.preconditioner="diagonal")
  if (!is.list(control) || (length(control) &&
      (is.null(names(control)) || any(!nzchar(names(control))))))
    stop("control must be a named list.", call.=FALSE)
  # Legacy controls passed directly through ... remain usable.
  for (name in intersect(names(dots), names(defaults))) {
    if (is.null(control[[name]])) control[[name]] <- dots[[name]]
  }
  if (!is.null(control$aitken)) {
    control$stopping <- if (isTRUE(control$aitken)) "aitken" else "increment"
    control$aitken <- NULL
  }
  unknown <- setdiff(names(control), c(names(defaults), "seed"))
  if (length(unknown))
    stop("Unknown control argument(s): ", paste(unknown, collapse=", "),
         call.=FALSE)
  for (name in names(control)) defaults[name] <- control[name]
  for (name in c("max.iter", "min.iter", "ncandid", "nstart", "nrep",
                 "init.iter", "backtracks", "pcg.max.iter")) {
    defaults[[name]] <- .pmbfa_integer(defaults[[name]], name,
      if (name %in% c("max.iter","init.iter","backtracks","pcg.max.iter"))
        0L else 1L)
    if (length(defaults[[name]]) != 1L)
      stop(name, " must be scalar.", call.=FALSE)
  }
  if (defaults$nrep > defaults$nstart || defaults$nstart > defaults$ncandid)
    stop("Require nrep <= nstart <= ncandid.", call.=FALSE)
  if (!is.null(defaults$tol) &&
      (!is.numeric(defaults$tol) || length(defaults$tol) != 1L ||
       !is.finite(defaults$tol) || defaults$tol <= 0))
    stop("tol must be NULL or a positive finite scalar.", call.=FALSE)
  for (name in c("mag", "min.component", "variance.floor", "nu.min",
                 "nu.max", "shape.min", "pcg.tol")) {
    value <- defaults[[name]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < 0 || (name != "min.component" && value == 0))
      stop(name, " must be finite and ",
        if (name == "min.component") "nonnegative." else "positive.",call.=FALSE)
  }
  if (defaults$nu.max < defaults$nu.min)
    stop("nu.max must be at least nu.min.",call.=FALSE)
  defaults$stopping <- match.arg(defaults$stopping,
                                  c("increment", "aitken"))
  defaults$init_method <- match.arg(defaults$init_method,c("nearest","kmeans"))
  defaults$pcg.preconditioner <- match.arg(defaults$pcg.preconditioner,
                                          c("diagonal","component"))
  defaults$criterion <- match.arg(toupper(defaults$criterion),c("BIC","ICL"))
  defaults$verbose <- isTRUE(defaults$verbose)
  defaults$profile <- isTRUE(defaults$profile)
  defaults
}

