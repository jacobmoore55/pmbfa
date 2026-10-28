# Skew-t search timing and memory

Measurements on 2026-10-27 used R 4.5.0, Rcpp 1.1.2, macOS arm64 and 8 GiB
physical RAM. These are single full-search runs on one machine, rather than
general performance guarantees. Compilation and package loading were excluded
from fit elapsed times. Each run used a fresh, separate R process. Resident
memory was sampled every 0.25 seconds; peaks include process dependencies.

The data are the normal-mixture demo's 200 observations of 8 by 8 matrices,
with two balanced components, true q=r=2 and covariance model CCC-CUC. Here the
fitted family is ST, the default for search_pmsbfa.

## Default search comparison

| Implementation | Fit elapsed, seconds | Sampled peak RSS, MiB | Factor pairs tested |
|---|---:|---:|---:|
| Local legacy MVBFA-6 | 79.768 | 256.8 | 20 |
| pmbfa 0.1.6 | 477.244 | 388.4 | 30 |
| pmbfa 0.1.7 | 406.948 | 94.5 | 30 |

The current call was `models <- search_pmsbfa(x, G=1:4, verbose=TRUE)`.
The legacy call was `models <- search_pmsbfa(x, G=1:4)` in a separate session
loaded from the local MVBFA-6 sources and freshly compiled helpers. The
legacy build used a diagnostic compatibility shim around external Rcpp
headers; the numerical source was retained.

The working-array change reduced observed elapsed time by 14.7% and sampled
peak resident memory by 75.7%. All four returned parameter sets, likelihood
curves, posteriors, factor-search tables and constraint-score tables were
bit-for-bit identical between 0.1.6 and 0.1.7. Both returned approximately
230 KiB of fit objects. Peak working memory should be distinguished from
the size of the returned result.

### Why the legacy default is faster here

These default calls perform different searches:

* The available legacy skewed selector searches backward. The package searches
  both forward and backward by default.
* The legacy factor search visits five pairs per G on this example. The
  package visits 7, 7, 9 and 7 pairs, with additional boundary and neighbor
  exploration.
* Legacy constraint stages reuse fitted parameters, with an initialization
  fallback. The package fits candidates from a fixed shared initialization
  bank for each (G,q,r), as specified for comparable constraint starts.
* Legacy intermediate adaptive fits can use a likelihood tolerance ten times
  larger, followed by tighter final stages. Its explicit tolerance is also
  reset at iteration 10. The package honors its stopping controls throughout.
* The legacy selector can skip jointly constrained candidates using fitted
  parent likelihoods. Nested screening is opt-in in the package because
  finite parent fits do not certify global likelihood bounds.
* The implementations differ in conditional-moment arithmetic. For example,
  the package uses a four-point Bessel-order derivative; legacy Elog mixing
  accumulation skips posterior weights at or below 1e-12.

These choices were retained in 0.1.7. The new search is still slower than the
local legacy default on this example. Reverting search coverage, starting
values or stopping rules would change the fitting procedure and can change
which model is selected.

## Fixed-model comparison

A separate comparison used the same data, G=2, q=r=2, model CCC-CUC-C-C and
the same perturbed starting parameters. It performed eight completed update
cycles, excluding initialization and search. Five timed repetitions followed
one warm-up. Native fits used max.iter=8 and min.iter=8; legacy max.iter=9
accounts for its initial likelihood evaluation before eight updates. Both
used tol=1e-15, and the legacy run ended before its iteration-10 tolerance reset.

| Implementation | Median elapsed, seconds |
|---|---:|
| Legacy update loop | 0.052 |
| Native 0.1.6 | 0.044 |
| Native 0.1.7 | 0.041 |

The native versions returned identical parameters and curves. Legacy and
native calculations are close but do not implement identical numerical
updates: their nine likelihood values differed by at most 0.001602 in this
short test. This is an illustrative fixed-model comparison, not a claim that
every native fit is faster.

## Allocation change

A CPU stack sample during the original search found 180 of 261 leaf samples
inside the Bessel kernel, with additional time in allocation and garbage
collection. R's ordinary bessel_k creates and releases an R work vector on
each evaluation. The work is released; the issue measured here is allocation
churn and a high working-memory peak.

Version 0.1.7 supplies a C++ buffer to bessel_k_ex, which calls the same
Rmath kernel. Each fit owns and releases its buffer. There is no persistent
global Bessel cache. This also applies to posterior-moment evaluation and GH
shape updates.

An Rprofmem check on one 120-observation ST density evaluation recorded:

| Version | R bytes allocated at nu=6 | R bytes allocated at nu=200 |
|---|---:|---:|
| 0.1.6 | 5,888 | 228,608 |
| 0.1.7 | 5,888 | 5,888 |

Rprofmem does not include the C++ buffer; that buffer scales with the largest
required Bessel order in the fit, rather than the number of evaluations.
RSS measurements above include native memory as well as R allocations.
The regression test checks that increasing the Bessel work length does not
reintroduce per-observation R working-array allocations.

The interface and lifetime can be inspected in the
[Rmath source](https://github.com/wch/r-source/blob/f8d4d7d48051860cc695b99db9be9cf439aee743/src/nmath/bessel_k.c)
and [R's memory-allocation documentation](https://stat.ethz.ch/R-manual/R-devel/doc/manual/R-exts.html#Transient-storage-allocation).

## Criterion comparison

Both default searches selected G=2, q=r=2, CCC-CUC-C-C. Legacy reported BIC
41210.43; scoring its fitted parameters with the package's density and
parameter count gives 41215.73. The package fit gives 41216.08. Raw legacy
and package criteria therefore should not be compared without reconciling
parameter counting. The near equality for G=2 also does not imply that the
two searches visited the same alternatives.

## Legacy source provenance

The comparison uses the local Desktop/MVBFA-6 copy, whose skewed forward
selector is absent. A newer legacy copy could behave differently. SHA-256:

* MVBFA.R: `85a2c39d12472f3ecda7f362dc6de09b8eac5edfe51d1369e5ac76203af7fa61`
* helpers.cpp: `000b533f9bfdc4d557b8bfaed88993250777c5d47fadd991c2d91f94f6df9b80`

To reproduce the package side, generate x using the demo's data-generation
section and time the ST search in a fresh session. The demo itself fits the
Normal family; timing that call measures a different workload.
