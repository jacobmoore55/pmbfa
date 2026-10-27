\name{pmbfa-methods}
\alias{print.pmbfa}
\alias{print.pmbfa_fits}
\alias{logLik.pmbfa}
\alias{BIC.pmbfa}
\alias{nobs.pmbfa}
\title{Methods for Bilinear Factor Mixture Fits}
\description{Print fits and retrieve their accepted likelihood, raw BIC, and observation count.}
\usage{
\method{print}{pmbfa}(x, ...)
\method{print}{pmbfa_fits}(x, ...)
\method{logLik}{pmbfa}(object, ...)
\method{BIC}{pmbfa}(object, ...)
\method{nobs}{pmbfa}(object, ...)
}
\arguments{
\item{x, object}{A fitted \code{pmbfa} object, or a \code{pmbfa_fits} list for its print method.}
\item{...}{Additional print arguments or reserved generic arguments.}
}
\value{Printing invisibly returns its input. \code{logLik} returns a \code{logLik} object with parameter and observation counts. \code{BIC} returns the raw stored BIC, including for a nonconverged fit when available. \code{nobs} returns N.}
\details{Skewed fits inherit from \code{pmsbfa} and \code{pmbfa}; the common methods apply to both families. Printed likelihoods, criteria and elapsed seconds use two decimals. Fits also print their stored call when present; grid printing reports overall elapsed time and call. A finite raw BIC does not make a nonconverged fit eligible for selection.}
\seealso{\link{pmbfa}, \link[stats]{BIC}}
