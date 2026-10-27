# Internal organization: see inst/DEVELOPMENT.md.
# Print one selected converged model at a public fitting/search boundary.
.pmbfa_report_fit <- function(fit,verbose=fit$control$verbose) {
  if (isTRUE(verbose) && isTRUE(fit$converged) &&
      length(fit$score) == 1L && is.finite(fit$score)) {
    m <- fit$params$moddat
    cat(sprintf("G=%d q=%d r=%d %s %s=%.2f (%s)\n",
      m$G,m$q,m$r,m$name,fit$criterion,fit$score,fit$status))
  }
  invisible(fit)
}

# Print a fit summary, original call and elapsed seconds.
print.pmbfa <- function(x,...) {
  m <- x$params$moddat
  cat(sprintf("%s bilinear factor mixture: G=%d, q=%d, r=%d\n",
               m$family,m$G,m$q,m$r))
  cat(sprintf("  Model: %s | %s | cycles: %d\n",m$name,x$status,x$iter))
  cat(sprintf("  logLik: %.2f | BIC: %.2f | ICL: %.2f | parameters: %s\n",
               .pmbfa_final_ll(x),x$bic,x$icl,as.character(x$npar)))
  if (!is.null(x$elapsed))cat(sprintf("  Elapsed: %.2f seconds\n",x$elapsed))
  if (!is.null(x$call))cat("  Call: ",paste(deparse(x$call,width.cutoff=80L),collapse="\n        "),"\n",sep="")
  if (length(x$message) && nzchar(x$message)) cat("  ",x$message,"\n",sep="")
  invisible(x)
}

# Print the ordered model list and whole-request metadata.
print.pmbfa_fits <- function(x,...) {
  table <- do.call(rbind,lapply(x,function(fit) {
    m <- fit$params$moddat
    data.frame(G=m$G,q=m$q,r=m$r,family=m$family,model=m$name,
               score=sprintf("%.2f",fit$score),seconds=if(is.null(fit$elapsed))NA_character_ else sprintf("%.2f",fit$elapsed),converged=fit$converged,status=fit$status)
  }))
  print(table,row.names=FALSE,...)
  if (!is.null(attr(x,"elapsed")))cat(sprintf("Total elapsed: %.2f seconds\n",attr(x,"elapsed")))
  if (!is.null(attr(x,"call")))cat("Call: ",paste(deparse(attr(x,"call"),width.cutoff=80L),collapse="\n      "),"\n",sep="")
  invisible(x)
}

# Return the accepted log likelihood with parameter and observation counts.
logLik.pmbfa <- function(object,...) {
  structure(if (length(object$loglik)) tail(object$loglik,1L) else NA_real_,
             df=object$npar,nobs=object$params$moddat$N,class="logLik")
}

# Return the stored raw BIC, including for diagnostic nonconverged fits.
BIC.pmbfa <- function(object,...) object$bic

# Return the number of matrix observations.
nobs.pmbfa <- function(object,...) object$params$moddat$N

