# Internal organization: see inst/DEVELOPMENT.md.
# Evaluate an expression with a temporary seed and restore the caller RNG state.
.pmbfa_with_seed <- function(seed, expression) {
  if (is.null(seed)) return(force(expression))
  had_seed <- exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir=.GlobalEnv)
  on.exit(if (had_seed) assign(".Random.seed", old_seed, envir=.GlobalEnv)
          else if (exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE))
            rm(".Random.seed", envir=.GlobalEnv), add=TRUE)
  set.seed(seed)
  force(expression)
}

# Derive a reproducible seed for a component/factor/family combination.
.pmbfa_tuple_seed <- function(seed, G, q, r, family) {
  if (is.null(seed)) return(NULL)
  as.integer((as.double(seed) + 104729*G + 13007*q + 433*r +
                97*match(family, c("Normal", "ST", "GH", "NIG", "VG"))) %%
               (.Machine$integer.max-1L) + 1L)
}

