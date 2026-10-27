# Wall-clock timing and the original user invocation are recorded at API boundaries.
# Calculate nonnegative elapsed wall-clock seconds from a recorded start.
.pmbfa_elapsed <- function(started) {
  unname(max(0,proc.time()[["elapsed"]]-started))
}

# Add per-tuple elapsed time to the fitting-engine result.
.pmbfa_timed_tuple <- function(...) {
  started <- proc.time()[["elapsed"]]
  fit <- .pmbfa_run_tuple(...)
  fit$elapsed <- .pmbfa_elapsed(started)
  fit
}

# Attach original-call and whole-request metadata without adding list entries.
.pmbfa_record_result <- function(result,requested_call,elapsed) {
  if(inherits(result,"pmbfa_fits")) {
    for(i in seq_along(result))result[[i]]$call <- requested_call
    attr(result,"call") <- requested_call
    attr(result,"elapsed") <- elapsed
  } else {
    result$call <- requested_call
    result$elapsed <- elapsed
  }
  result
}
