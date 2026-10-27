\name{pmbfa-control}
\alias{pmbfa-control}
\title{AECM and Initialization Controls}
\description{Named entries accepted in the control list of the fitting and search functions.}
\details{
\describe{
\item{max.iter, min.iter}{Maximum and minimum completed cycles; defaults 10000 and 2. Zero maximum evaluates only the supplied initial state.}
\item{tol, mag}{Positive absolute likelihood stopping threshold, or NULL. An explicit tol is never reset. If NULL, the threshold is calibrated once as \code{max(1,abs(initial_loglik))*10^(-mag)}, with default \code{mag=log10(sqrt(N)*n*p)+3}.}
\item{stopping}{\code{"increment"} or \code{"aitken"}. Aitken uses an admissible extrapolated remaining likelihood gain; it does not extrapolate parameters.}
\item{ncandid, nstart, nrep, init.iter}{Initialization candidate count, shortlisted count, retained starts, and short-fit cycles. Defaults 200, 20, 1 and 5. Require \code{nrep <= nstart <= ncandid}.}
\item{init_method}{\code{"nearest"} by default; \code{"kmeans"} calls Lloyd k-means during initialization. Failed k-means calls fall back to the nearest candidate.}
\item{variance.floor}{Lower bound on diagonal row and column residual terms; default 0.1.}
\item{nu.min, nu.max}{Skew-t degrees-of-freedom bounds; defaults 4.1 and 1e12.}
\item{shape.min}{Positive lower bound for NIG/VG shape updates; default 0.5.}
\item{min.component}{Every accepted state's effective component sizes must be strictly greater than this value; default 6. A proposed collapse is rejected before convergence testing.}
\item{pcg.tol, pcg.max.iter}{Shared-skewness residual tolerance (1e-9) and per-attempt inner iteration limit (0 means n*p). The true residual is verified.}
\item{pcg.preconditioner}{\code{"diagonal"} by default, or \code{"component"}. The latter can help aligned spatial covariance systems but may be slower with heterogeneous scales. Failed component attempts retry diagonal from the original starting matrix. The per-attempt budget applies to both attempts.}
\item{backtracks}{Maximum likelihood-monotonicity damping attempts; default 12.}
\item{verbose, profile}{Print progress, or record numerical phase timing; defaults FALSE. Profiling also returns \code{pcg_iterations_by_cycle} for accepted cycles. Stored values retain full precision even though fitted-object printing uses two decimals.}
\item{criterion}{Internal score label, overwritten by the public criterion argument.}
}
}
\section{Numerical interpretation}{Absolute variance floors depend on the observation units. They prevent collapse but do not impose an upper condition-number limit on a scale matrix. Numerical failures are returned with the last accepted state and an infinite selection score. No automatic diagonal jitter or undocumented change to the likelihood is applied.}
\seealso{\link{pmbfa}, \link{search_pmbfa}}
