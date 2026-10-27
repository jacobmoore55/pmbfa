# Internal organization: see inst/DEVELOPMENT.md.
# Fit skewed component mixtures for a chosen distribution family and specified dimensions.
pmsbfa <- function(x, G=2L, q=1:2, r=1:2, class=NULL, modelSubset=NULL,
                    family="ST", seed=123456L, adaptive=TRUE,
                    adaptive.method=c("both","forward","backward"),
                    criterion=c("BIC","ICL"), verbose=FALSE,
                    parallelize=FALSE, ncores=NULL, inits=NULL,
                    control=list(), screening=c("none","nested"), ...) {
  requested_call <- match.call()
  started <- proc.time()[["elapsed"]]
  if (missing(q) && length(dim(x)) == 3L)
    q <- seq_len(max(1L,min(2L,.pmbfa_max_factors(dim(x)[1]))))
  if (missing(r) && length(dim(x)) == 3L)
    r <- seq_len(max(1L,min(2L,.pmbfa_max_factors(dim(x)[2]))))
  result <- .pmbfa_fit(x,G,q,r,class,modelSubset,family,seed,adaptive,
    adaptive.method[1L],criterion[1L],verbose,parallelize,ncores,inits,control,
    list(...),screening=screening[1L])
  .pmbfa_record_result(result,requested_call,.pmbfa_elapsed(started))
}

# Fit normal component mixtures for specified component and factor counts.
pmbfa <- function(x, G=2L, q=1:2, r=1:2, class=NULL, modelSubset=NULL,
                   seed=123456L, adaptive=TRUE,
                   adaptive.method=c("both","forward","backward"),
                   criterion=c("BIC","ICL"), verbose=FALSE,
                   parallelize=FALSE, ncores=NULL, inits=NULL,
                   control=list(), screening=c("none","nested"), ...) {
  requested_call <- match.call()
  started <- proc.time()[["elapsed"]]
  if (missing(q) && length(dim(x)) == 3L)
    q <- seq_len(max(1L,min(2L,.pmbfa_max_factors(dim(x)[1]))))
  if (missing(r) && length(dim(x)) == 3L)
    r <- seq_len(max(1L,min(2L,.pmbfa_max_factors(dim(x)[2]))))
  result <- .pmbfa_fit(x,G,q,r,class,modelSubset,"Normal",seed,adaptive,
    adaptive.method[1L],criterion[1L],verbose,parallelize,ncores,inits,control,
    list(...),screening=screening[1L])
  .pmbfa_record_result(result,requested_call,.pmbfa_elapsed(started))
}

