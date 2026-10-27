# pmbfa 0.1.5

* Prints only the selected constraint model for each (G,q,r) combination in pmbfa/pmsbfa, and the selected factor/constraint model for each G in search_pmbfa/search_pmsbfa, when verbose=TRUE.
* Emits parallel-grid model summaries in the parent process so they are visible with both fork and socket workers.

# pmbfa 0.1.4

* Makes verbose output one criterion summary per converged model fit.
* Adds control$track_likelihood (FALSE by default) for independent per-cycle likelihood output; initialization remains quiet.
* Updates the demo to use search_pmbfa(x, G = 1:4, verbose = TRUE) and keeps diagnostic assertions in tests.
* Verifies that both factor searches default to adaptive factor selection and both directions of adaptive constraint selection.

# pmbfa 0.1.3

* Makes the demo fit with exactly search_pmbfa(x, G=1:4), exercising all public defaults.
* Clarifies the default adaptive factor search, factor limits, initialization and both covariance-search directions.

# pmbfa 0.1.2

* Corrects the demo to generate normal PMBFA data with true q=r=2 and covariance model CCC-CUC.
* Runs both forward and backward adaptive covariance searches in the G=1:4 demo, without fixing modelSubset.
* Clarifies the distinct factor-search and covariance-search arguments and tests the existing both-direction defaults.

# pmbfa 0.1.1

* Adds original user calls and elapsed wall-clock seconds to all four entry points, including grid-level metadata and per-result timing.
* Provides separate function help pages, an expanded package overview, and an executable balanced 8 by 8 matrix-mixture demo searching G=1:4.
* Keeps the Desktop source tree focused on the package; historical standalone tests, fits, and research outputs are archived separately.

# pmbfa 0.1.0

* Replaces the legacy package scaffold with the current normal/skewed AECM implementation.
* Exports four fitting and search functions, with registered S3 and native methods.
* Separates fitting, initialization, distributions, linear algebra, and model selection into modules.
* Preserves variance and shape bounds; enforces effective component sizes strictly above the configured threshold at every accepted cycle.
* Includes package-level numerical, constraint, failure, selection, and parallel tests.
