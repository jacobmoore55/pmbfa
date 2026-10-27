\name{pmsbfa}
\alias{pmsbfa}
\title{Fit Skewed Bilinear Factor Mixtures}
\description{Fits parsimonious mixtures of matrix-variate skewed bilinear factor analyzers for specified component and factor counts. Supported families are skew-t, generalized hyperbolic, normal inverse Gaussian and variance-gamma, with a normal-family option. Skewness and concentration can be common or component-specific.}
\usage{
pmsbfa(x, G = 2L, q = 1:2, r = 1:2, class = NULL,
  modelSubset = NULL, family = "ST", seed = 123456L, adaptive = TRUE,
  adaptive.method = c("both", "forward", "backward"),
  criterion = c("BIC", "ICL"), verbose = FALSE, parallelize = FALSE,
  ncores = NULL, inits = NULL, control = list(),
  screening = c("none", "nested"), ...)
}
\arguments{
\item{x}{A finite numeric array with dimensions \code{c(n, p, N)}, with \code{N > 10}.}
\item{G}{One or more positive integer component counts.}
\item{q, r}{Row and column factor counts: Lambda has dimensions \code{n by q} and Delta \code{p by r}. Vectors request a grid.}
\item{class}{Optional known labels, one per observation. Zero and \code{NA} denote unknown labels.}
\item{modelSubset}{Optional character vector of constraint strings. Supplying it fits those models directly instead of an adaptive constraint search.}
\item{family}{For \code{pmsbfa}, one of \code{"ST"}, \code{"GH"}, \code{"NIG"}, \code{"VG"}, or \code{"Normal"}.}
\item{seed}{A scalar integer seed, or \code{NULL}. Fixed seeds are reproducible and restore the caller's RNG state.}
\item{adaptive}{Logical, or a search direction string. \code{FALSE} fits all allowed constraint models when \code{modelSubset} is absent.}
\item{adaptive.method}{Constraint search direction: \code{"both"} (the default), \code{"forward"}, or \code{"backward"}. Both paths run when adaptive selection is enabled and no modelSubset is supplied.}
\item{criterion}{\code{"BIC"} or \code{"ICL"}; smaller values are preferred.}
\item{verbose}{Print model and iteration progress.}
\item{parallelize}{Parallelize component/factor grid combinations. A single native fit is serial.}
\item{ncores}{Optional positive integer worker count.}
\item{inits}{Optional parameter list or bank of lists for one component/factor combination. See Details.}
\item{control}{Named algorithm controls; see \link{pmbfa-control}.}
\item{screening}{\code{"none"} by default, or \code{"nested"} for optional backward nested-model screening.}
\item{...}{Controls can also be passed directly. Legacy aliases \code{rG}, \code{rq}, \code{rr}, \code{icl}, \code{track}, \code{quiet}, and \code{aitken} are recognized. Unknown arguments error.}
}
\details{
Each native call executes a complete alternating expectation conditional maximization (AECM) fit. All constraint models for the same \code{(G,q,r,family)} reuse the same initialization bank. A previously fitted constraint model is not a warm start for another model.

A covariance block contains three C/U flags: common versus component-specific loadings, common versus component-specific diagonal residuals, and isotropic versus diagonal residuals. Normal model names have the form \code{"UUU-UUU"}. Skewed model names append common/component-specific skewness and concentration flags, for example \code{"UUU-CUU-C-U"}. For GH the final flag concerns omega; lambda remains component-specific.

Adaptive selection uses staged comparisons, and factor selection is also a heuristic unless its exhaustive option is requested. Optional nested screening uses converged local parent fits, which do not certify a global likelihood bound. It can therefore change a finite-start search result. Skipped models have records in \code{selection$screened.models}, rather than fabricated fits.

Supplied initializations contain \code{moddat}, \code{gmp}, and component \code{args}. The metadata identify G, q, r, family and model; \code{gmp} contains \code{Sigmak}, \code{Psik}, \code{Lambdak}, \code{Deltak}, and \code{pik}. Each component argument contains M and, for skewed families, A and its shape parameters.

The defaults preserve residual variance floors of 0.1, skew-t degrees of freedom of at least 4.1, and effective component sizes strictly greater than 6. Invalid input formats error. Numerical failures retain the last accepted full-cycle parameter state and its matching likelihood and posterior. Such a fit has \code{converged=FALSE} and selection score \code{Inf}, even when its raw criterion is finite. User interrupts propagate.
}
\value{
A \code{pmbfa} object for one component/factor combination. Skewed fits also inherit from \code{pmsbfa}. Grids return a \code{pmbfa_fits} list ordered by selection score. A fit contains parameters, initial and accepted-cycle likelihoods, iteration count, posterior probabilities, assignments, convergence status and message, parameter count, raw BIC/ICL, selection score, and numerical diagnostics. The call and elapsed fields are described below. Raw information criteria remain available after numerical stopping when a finite accepted likelihood exists. The criterion score is infinite for a fit that did not converge.
}

\section{Call and elapsed time}{
For a single returned fit, \code{call} is the original matched user call and \code{elapsed} is wall-clock time in seconds for the complete public call, including initialization and model selection.

A \code{pmbfa_fits} list stores the overall call and wall-clock time in \code{attr(result,"call")} and \code{attr(result,"elapsed")}. Each contained fit also has the original user call. Its \code{elapsed} covers its component/factor combination for direct fitting, or the complete factor search for its G in a search result. The original call may request a grid; the selected combination is recorded in \code{params$moddat}.

For parallel grids, overall elapsed time includes worker overhead and is not the sum of worker times. These timings measure complete user requests rather than isolated native iteration speed. The recorded call contains the user's expression; the data must still be available in the evaluation environment to rerun it.
}

\seealso{\link{search_pmbfa}, \link{pmbfa-control}}
\examples{
set.seed(1)
x <- array(rnorm(4 * 4 * 60), c(4, 4, 60))
fit <- pmsbfa(x, G = 1, q = 1, r = 1, family = "ST",
  modelSubset = "CCC-CCC-C-C",
  control = list(ncandid = 3, nstart = 1, init.iter = 1, max.iter = 20))
fit$call
fit$elapsed
}
