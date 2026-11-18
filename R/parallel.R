# Socket workers load the installed namespace and its registered native symbols.
# Run grid tuples in socket workers using the installed package namespace.
.pmbfa_psock <- function(tuples, arguments, ncores) {
  workers <- .pmbfa_cluster(ncores)
  on.exit(parallel::stopCluster(workers), add=TRUE)
  parallel::parLapply(workers, seq_len(nrow(tuples)),
    function(i, tuples, arguments) {
      runner <- getFromNamespace(".pmbfa_timed_tuple", "pmbfa")
      do.call(runner, c(arguments, as.list(tuples[i,,drop=FALSE])))
    }, tuples=tuples, arguments=arguments)
}

.pmbfa_worker_count <- function(parallelize,ncores,tasks) {
  if (!is.logical(parallelize) || length(parallelize) != 1L ||
      is.na(parallelize)) stop("parallelize must be TRUE or FALSE.", call.=FALSE)
  if (!is.null(ncores)) {
    ncores <- .pmbfa_integer(ncores,"ncores")
    if (length(ncores) != 1L) stop("ncores must be scalar.",call.=FALSE)
  }
  if (!parallelize || tasks == 1L) return(1L)
  if (is.null(ncores))
    ncores <- max(1L,parallel::detectCores(logical=FALSE)-1L,na.rm=TRUE)
  min(ncores,tasks)
}

.pmbfa_cluster <- function(ncores) {
  workers <- parallel::makeCluster(ncores)
  on.exit(parallel::stopCluster(workers))
  parallel::clusterCall(workers,function(paths) {
    .libPaths(paths)
    loadNamespace("pmbfa")
    NULL
  },.libPaths())
  on.exit(NULL)
  workers
}
