# Internal organization: see inst/DEVELOPMENT.md.
# Enumerate the normal or skewed constraint models.
.pmbfa_models <- function(family) {
  bits <- if (family == "Normal") 6L else 8L
  values <- expand.grid(rep(list(c("C", "U")), bits),
                        stringsAsFactors=FALSE)
  apply(values, 1L, function(z) {
    name <- paste0(paste0(z[1:3], collapse=""), "-",
                   paste0(z[4:6], collapse=""))
    if (bits == 8L) name <- paste0(name, "-", z[7], "-", z[8])
    name
  })
}

# Remove redundant cross-component flags when there is only one component.
.pmbfa_canonical <- function(name, G) {
  if (G != 1L) return(name)
  positions <- c(1L, 2L, 5L, 6L)
  if (nchar(name) == 11L) positions <- c(positions, 9L, 11L)
  for (position in positions) substr(name, position, position) <- "C"
  name
}

# Generate the constraint candidates for one adaptive stage.
.pmbfa_stage_names <- function(name, positions, G) {
  bits <- expand.grid(rep(list(c("C", "U")), length(positions)),
                      stringsAsFactors=FALSE)
  unique(apply(bits, 1L, function(values) {
    candidate <- name
    for (i in seq_along(positions))
      substr(candidate, positions[i], positions[i]) <- values[i]
    .pmbfa_canonical(candidate, G)
  }))
}

# Compute the optional backward nested-model screening decision and bound.
.pmbfa_nested_screen <- function(child, parents, incumbent, criterion,
                                  npar_model) {
  no <- list(screened=FALSE, lower.bound=NA_real_, parents=character(),
             reason="")
  if (length(parents) != 2L || !is.function(npar_model) ||
      !is.finite(incumbent) ||
      !all(vapply(parents,function(fit)isTRUE(fit$converged),logical(1))))
    return(no)
  metadata <- lapply(parents,function(fit)fit$params$moddat)
  keys <- c("n","p","N","G","q","r","family")
  if (!all(vapply(keys,function(key)
      identical(metadata[[1]][[key]],metadata[[2]][[key]]),logical(1))))
    return(no)
  names <- vapply(metadata,`[[`,character(1),"name")
  nested <- vapply(names,function(parent) {
    c <- strsplit(child, "",fixed=TRUE)[[1]]
    p <- strsplit(parent,"",fixed=TRUE)[[1]]
    if (length(c) != length(p)) return(FALSE)
    changed <- which(c != p)
    length(changed) == 1L && c[changed] == "C" && p[changed] == "U"
  },logical(1))
  if (!all(nested) || anyDuplicated(names)) return(no)
  ll <- vapply(parents,.pmbfa_final_ll,numeric(1))
  if (any(!is.finite(ll))) return(no)
  k <- npar_model(child)
  N <- metadata[[1]]$N
  if (length(k) != 1L || !is.finite(k) || !is.numeric(N) ||
      length(N) != 1L || !is.finite(N) || N <= 1) return(no)
  lower <- -2*min(ll) + k*log(N)
  # This is a bound only if parent likelihoods are globally maximized. With
  # finite multistart EM, it is an explicitly optional screening heuristic.
  # The same BIC bound can screen ICL since its classification penalty >= 0.
  margin <- 8*.Machine$double.eps*max(1,abs(lower),abs(incumbent))
  list(screened=is.finite(lower) && lower > incumbent+margin,
       lower.bound=lower, parents=names,
       reason=if (criterion == "ICL") "nested_BIC_bound_for_ICL" else
         "nested_BIC_bound")
}

