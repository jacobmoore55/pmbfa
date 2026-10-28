# pmbfa

[![R-CMD-check](https://github.com/jacobmoore55/pmbfa/actions/workflows/R-CMD-check.yaml/badge.svg?branch=main)](https://github.com/jacobmoore55/pmbfa/actions/workflows/R-CMD-check.yaml)

`pmbfa` fits parsimonious mixtures of normal and skewed matrix-variate bilinear
factor analyzers. Each observation is a matrix. Low-rank factors describe row
and column dependence; model constraints share loadings, residual variances,
skewness or concentration across components. Fitting and initialization run in
C++, with model searches in R.

## Functions

| Function | Use |
|---|---|
| `pmbfa()` | Fit normal mixtures for chosen component and factor counts. |
| `pmsbfa()` | Fit ST, GH, NIG or VG mixtures, with a normal-family option. |
| `search_pmbfa()` | Compare normal mixtures across component and factor counts. |
| `search_pmsbfa()` | Compare skewed mixtures across component and factor counts. |

Each function has its own help page. `help("pmbfa-package")` introduces the
package; `help("pmbfa-control")` describes algorithm controls.

## Install

Requires R >= 4.1, Rcpp >= 1.1.0 and a C++17 compiler for source installation.
Install from GitHub:

```r
install.packages("remotes")
remotes::install_github("jacobmoore55/pmbfa")
library(pmbfa)
```

Alternatively, from the package source directory (for example, a local
repository clone):

```r
install.packages(".", repos = NULL, type = "source")
library(pmbfa)
```

### Help pages after reinstalling

If R reports that `pmbfa/help/pmbfa.rdb` is corrupt immediately after updating
the package, restart R before loading the new version. In RStudio, use
**Session > Restart R**. A running session can retain references to help
databases from the previous installation.

Then check the installed help in the fresh session:

```r
library(pmbfa)
help("search_pmbfa", package = "pmbfa")
```

If the error persists in a fresh session, close sessions using the package and
reinstall it from source. To check every installed help topic without fitting a
model, run `tools::Rd_db("pmbfa")` in a fresh R session.

## A rank-two PMBFA demo

The data are generated from a **normal PMBFA** model, with 200 observations of
8 by 8 matrices split into two balanced components. The true factor dimensions
are q=r=2 and the covariance model is CCC-CUC: common row loadings and common
isotropic row residuals, common column loadings, and component-specific
isotropic column residuals. Locations differ between components.

```r
demo("normal_mixture", package = "pmbfa", ask = FALSE, echo = FALSE)
```

`echo=FALSE` turns off R's echo of the demo source code. The demo prints
exactly four selected-model summaries, one for each G in 1:4.

The complete data generation and search code is:

```r
set.seed(2026)
true_class <- rep(1:2,each=100L)
true_locations <- list(
  matrix(-.75,8L,8L),
  .75+outer(seq(-.5,.5,length.out=8L),seq(-.5,.5,length.out=8L)))

# Common rank-two row/column loadings. Their columns are independent.
true_Lambda <- .7*cbind(rep(c(1,-1),4L),rep(c(1,1,-1,-1),2L))
true_Delta <- .6*cbind(seq(-1,1,length.out=8L),rep(c(1,-1,-1,1),2L))
true_sigma <- .8
true_psi <- c(.5,1.2)
# CCC: common row loadings and one common isotropic residual variance.
true_row_scale <- tcrossprod(true_Lambda)+diag(true_sigma,8L)
# CUC: common column loadings, component-specific isotropic residuals.
true_column_scales <- lapply(true_psi,function(psi)
  tcrossprod(true_Delta)+diag(psi,8L))
true_model <- list(family="Normal",G=2L,q=2L,r=2L,model="CCC-CUC",
  Lambda=true_Lambda,Delta=true_Delta,sigma=true_sigma,psi=true_psi,
  row_scale=true_row_scale,column_scales=true_column_scales,M=true_locations)

# If Z has independent standard normal entries, L_row %*% Z %*% U_col
# has matrix-normal row scale Sigma and column scale Psi_g.
row_chol <- t(chol(true_row_scale))
column_chols <- lapply(true_column_scales,chol)
x <- array(0,c(8L,8L,200L))
for(i in seq_len(200L)) {
  g <- true_class[i]
  Z <- matrix(rnorm(64L),8L,8L)
  x[,,i] <- true_locations[[g]]+row_chol%*%Z%*%column_chols[[g]]
}
ordering <- sample.int(200L)
x <- x[,,ordering,drop=FALSE]
true_class <- true_class[ordering]

# Exercise the public defaults for factor search, constraints and initialization.
models <- search_pmbfa(x, G = 1:4, verbose = TRUE)

selected_model <- models[[1]]
```

### Two different search settings

`search.method` controls the factor-count search and defaults to `"adaptive"`.
The demo passes no search overrides. Factor limits are derived from the matrix
dimensions (up to four per side for these 8 by 8 observations); the true q=r=2
values specify the generating model, not a fitted search restriction.

Covariance-constraint search separately defaults to `adaptive=TRUE` and
`adaptive.method="both"`. No `modelSubset` is supplied, so both forward and
backward constraint paths run. Supplying `modelSubset` would bypass this search.

The returned selected fit records both histories:

```r
selected_model$selection$paths$forward
selected_model$selection$paths$backward
```

Initialization, seed, stopping tolerance and iteration budget also use their
public defaults. The demo prints only the selected-model summaries. Its
`models`, `comparison`, `selected_model` and `accuracy` objects can be inspected
afterward; fits also store their calls and elapsed times.

## Progress, calls and runtime

With `verbose=TRUE`, `pmbfa()` and `pmsbfa()` print the best constraint model
for each tested `(G,q,r)` combination. `search_pmbfa()` and `search_pmsbfa()`
print the best factor/constraint model for each tested `G`, after that G's
search finishes. Summaries describe selected converged fits and use the
requested criterion: BIC by default, or ICL when specified. For example:

```text
G=2 q=2 r=2 CCC-CUC BIC=40955.94 (converged)
```

The legacy argument `quiet=FALSE` enables the same winner summaries;
`quiet=TRUE` suppresses them. An explicit `quiet` setting overrides `verbose`
and `control$verbose`.

Set `control=list(track_likelihood=TRUE)` to also print the log likelihood and
its change after each accepted AECM cycle. This option defaults to `FALSE`
and works independently of `verbose`. The likelihood curve is always stored
in the fit. Initialization screening and short refinements remain quiet.

For a single returned fit:

```r
fit$call       # original matched user call
fit$elapsed    # wall-clock seconds for the complete public call
```

For a grid or vector of G:

```r
attr(models, "call")
attr(models, "elapsed")       # total public-request time
models[[1]]$call
models[[1]]$elapsed           # search time for this G
```

Each list entry remains a fit object. Direct-fit grids record per-(G,q,r)
times; searches record the complete factor search for each G. The list-level
time includes validation and orchestration/worker overhead. Stored calls retain
the original public expression instead of an internal `do.call` with literal
data. The selected parameters identify each entry's actual model dimensions.
These are end-to-end wall-clock timings; native phase timing remains available
with `control=list(profile=TRUE)`. Printing uses two decimal places while stored
values retain their precision. Diagnostic failure objects also carry metadata.

## Bounds and numerical behavior

Defaults remain `variance.floor=0.1`, `nu.min=4.1`, and strictly more than six
effective observations per component. Every proposed cycle is checked before
acceptance. Failed or unconverged fits retain diagnostic states/raw criteria
with infinite selection scores. Explicit `tol` is never overwritten. Optional
nested screening remains off by default because parent local optima do not
certify global likelihood bounds.

## Development and tests

The package source tree contains package metadata, R/C++ modules, help pages,
`demo/`, `inst/`, and `tests/`. Historical research scripts, saved fits and
benchmark outputs are maintained separately from the package repository.

`inst/DEVELOPMENT.md` maps the modules; `inst/FUNCTIONS.md` describes every R
function and the native implementation groups. Rcpp compiles at installation.
No external matrix library or runtime compilation is needed.

```sh
Rscript -e 'Rcpp::compileAttributes(".")'
R CMD build .
R CMD check --no-manual pmbfa_0.1.7.tar.gz
```

Package tests cover independent densities/moments, all normal/ST constraints,
numerical failure and component-size guards, adaptive searches, S3 methods,
RNG preservation, parallel workers, and call/timing metadata.

`inst/PERFORMANCE.md` records a measured skew-t search comparison with the
legacy implementation, including peak resident memory, search differences,
and a fixed-model timing comparison.

GitHub Actions runs package checks on Linux, macOS and Windows, plus Linux
checks with R-devel and the preceding R release. Checks run for pushes to
`main`, pull requests and manual workflow requests. Workflow action revisions
are pinned in `.github/workflows/R-CMD-check.yaml`.

## License

GPL (>= 2), as specified in `DESCRIPTION`. A copy of GPL version 2 is included
in `COPYING`; the package may be distributed under GPL version 2 or later.
