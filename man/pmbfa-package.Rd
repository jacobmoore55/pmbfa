\name{pmbfa-package}
\alias{pmbfa-package}
\docType{package}
\title{Parsimonious Mixtures of Bilinear Factor Analyzers}
\description{Fits mixtures of normal and skewed matrix-valued observations using low-rank row and column covariance factors. Component counts, factor dimensions, covariance sharing, skewness and concentration constraints can be compared using BIC or ICL. Numerical fitting and initialization run in compiled C++; the adaptive searches remain in R.}
\details{The four fitting entry points are \code{pmbfa}, \code{pmsbfa}, \code{search_pmbfa}, and \code{search_pmsbfa}. Native routines are registered and compiled during package installation.}
\section{Getting started}{
Use \code{pmbfa} and \code{pmsbfa} for specified component/factor combinations. Use \code{search_pmbfa} and \code{search_pmsbfa} to compare factor counts. All four record the original call and complete elapsed time. Each function has its own help page; \link{pmbfa-control} describes controls.

Run \code{demo("normal_mixture", package="pmbfa", ask=FALSE)} for normal PMBFA data with true \code{q=r=2}, constraints \code{CCC-CUC}, and 200 matrices of dimensions 8 by 8 split equally between two components. The sole fitting call \code{search_pmbfa(x,G=1:4)} uses the default adaptive factor search and both forward/backward constraint directions.
}
\author{Jacob Moore \email{jacob.c.moore55@gmail.com}}
\seealso{\link{pmbfa}, \link{search_pmbfa}, \link{pmbfa-control}}
\keyword{package}
