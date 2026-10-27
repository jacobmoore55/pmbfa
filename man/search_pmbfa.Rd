\name{search_pmbfa}
\alias{search_pmbfa}
\title{Search Normal Mixture Component and Factor Counts}
\description{Compares normal bilinear factor mixture models across candidate component counts and row/column factor counts. It returns the selected factor and covariance model for each requested component count, ordered by the information criterion.}
\usage{
search_pmbfa(x, G, max_q = NULL, max_r = NULL,
  search.method = c("adaptive", "grid"), ...)
}
\arguments{
\item{x}{A finite array with dimensions \code{c(n,p,N)}.}
\item{G}{One or more positive integer component counts.}
\item{max_q, max_r}{Maximum factor counts. Defaults use the identifiable diagonal factor-analysis rank limit.}
\item{search.method}{\code{"adaptive"} (the default) uses coarse-to-fine coordinate and neighbor comparisons; \code{"grid"} fits every factor pair up to the supplied maxima.}
\item{...}{Arguments forwarded to \code{pmbfa} or \code{pmsbfa}, including constraint selection, criterion, initialization and controls. Specify factor limits here instead of forwarding q or r.}
}
\details{By default, covariance constraints are searched adaptively in both forward and backward directions: \code{adaptive=TRUE} and \code{adaptive.method="both"}. These arguments are forwarded through \code{...}. Supplying \code{modelSubset} fits only the listed constraint models and bypasses that adaptive constraint search. The separate \code{search.method} argument controls factor-count search. Its adaptive option is local and does not guarantee the best factor combination on the entire grid. Constraint selection within each factor combination follows the arguments forwarded in \code{...}. Failed or unconverged fits cannot win a criterion comparison against a converged fit. If all candidates fail, the most informative available object is returned for diagnosis.}
\section{Progress output}{With \code{verbose=TRUE}, one summary is printed for the best converged factor/constraint model for each tested \code{G}, using the requested criterion (BIC by default). Intermediate constraint and factor candidates are not reported. Per-cycle likelihood output is enabled separately with \code{control=list(track_likelihood=TRUE)}.}
\value{For one G, the selected \code{pmbfa} object, including its \code{factor.search} table and search method. For a vector of G, a criterion-ordered \code{pmbfa_fits} list.}

\section{Call and elapsed time}{
For a single returned fit, \code{call} is the original matched user call and \code{elapsed} is wall-clock time in seconds for the complete public call, including initialization and model selection.

A \code{pmbfa_fits} list stores the overall call and wall-clock time in \code{attr(result,"call")} and \code{attr(result,"elapsed")}. Each contained fit also has the original user call. Its \code{elapsed} covers its component/factor combination for direct fitting, or the complete factor search for its G in a search result. The original call may request a grid; the selected combination is recorded in \code{params$moddat}.

For parallel grids, overall elapsed time includes worker overhead and is not the sum of worker times. These timings measure complete user requests rather than isolated native iteration speed. The recorded call contains the user's expression; the data must still be available in the evaluation environment to rerun it.
}

\section{Balanced normal-mixture demonstration}{
\code{demo("normal_mixture", package="pmbfa", ask=FALSE)} generates 200 normal PMBFA observations of 8 by 8 matrices in two balanced components, with distinct locations, true \code{q=r=2}, and covariance model \code{CCC-CUC}. Its only fitting call is \code{search_pmbfa(x, G = 1:4, verbose = TRUE)}. This exercises the default adaptive factor search, automatically determined factor limits, both forward/backward constraint searches, and default initialization and stopping controls. No \code{modelSubset} restricts the fitted covariance models. The demo prints criteria, convergence, timing and the recorded call.
}
\seealso{\link{pmbfa}, \link{pmbfa-control}}
\examples{
set.seed(2)
x <- array(rnorm(4 * 4 * 60), c(4, 4, 60))
fit <- search_pmbfa(x, G = 1, max_q = 1, max_r = 1,
  modelSubset = "CCC-CCC", search.method = "grid",
  control = list(ncandid = 3, nstart = 1, init.iter = 1,
                 max.iter = 50, tol = 1e-3))
fit$factor.search
}
