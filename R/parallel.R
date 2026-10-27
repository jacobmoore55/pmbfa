# Socket workers load the installed namespace and its registered native symbols.
# Run grid tuples in socket workers using the installed package namespace.
.pmbfa_psock <- function(tuples, arguments, ncores) {
  workers <- parallel::makeCluster(ncores)
  on.exit(parallel::stopCluster(workers), add=TRUE)
  parallel::clusterCall(workers, function(paths) {
    .libPaths(paths)
    loadNamespace("pmbfa")
    NULL
  }, .libPaths())
  parallel::parLapply(workers, seq_len(nrow(tuples)),
    function(i, tuples, arguments) {
      runner <- getFromNamespace(".pmbfa_timed_tuple", "pmbfa")
      do.call(runner, c(arguments, as.list(tuples[i,,drop=FALSE])))
    }, tuples=tuples, arguments=arguments)
}
