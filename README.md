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

The demo compares component counts 1 through 4 and displays the selected fit
and its classifications. Run the companion skew-t example with
`demo("skew_mixture", package = "pmbfa")`.

Both factor searches use adaptive factor selection and forward and backward
constraint selection by default. To run independent component counts concurrently,
use, for example, `search_pmsbfa(x, G = 1:4, parallelize = TRUE, ncores = 2)`.
Parallel searches report their winners in the supplied component-count order
after the workers finish. Returned fits are ordered by the selected criterion.

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
R CMD check --no-manual pmbfa_0.1.8.tar.gz
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
