# Function reference

Only the four fitting/search functions are exported. Registered S3 methods support fitted objects. Other R and native functions are implementation details; their signatures may change.

| R function | File | Purpose |
|---|---|---|
| `pmsbfa` | `R/api.R` | Fit skewed component mixtures for a chosen distribution family and specified dimensions. |
| `pmbfa` | `R/api.R` | Fit normal component mixtures for specified component and factor counts. |
| `.pmbfa_diag_logK` | `R/diagnostics.R` | Evaluate a stable scaled log Bessel K value for mean diagnostics. |
| `.pmbfa_prior_mixing_mean` | `R/diagnostics.R` | Compute the family-specific prior expectation of the mixing variable. |
| `pmbfa_component_mean` | `R/diagnostics.R` | Calculate a fitted population component mean matrix or its scalar pixel average. |
| `pmbfa_mean_diagnostics` | `R/diagnostics.R` | Compare population, posterior-based and empirical means, scores and conditioning. |
| `.pmbfa_failure` | `R/fit.R` | Construct a consistent diagnostic fit after an initialization or numerical failure. |
| `.pmbfa_final_ll` | `R/fit.R` | Retrieve the last finite accepted likelihood, or minus infinity if none exists. |
| `.pmbfa_best` | `R/fit.R` | Choose by the requested criterion while excluding nonconverged fits. |
| `.pmbfa_fit` | `R/fit.R` | Validate a public fitting request and orchestrate its component/factor grid. |
| `.pmbfa_with_seed` | `R/initialization.R` | Evaluate an expression with a temporary seed and restore the caller RNG state. |
| `.pmbfa_tuple_seed` | `R/initialization.R` | Derive a reproducible seed for a component/factor/family combination. |
| `.pmbfa_report_fit` | `R/methods.R` | Print one selected converged model at a public fitting/search boundary. |
| `print.pmbfa` | `R/methods.R` | Print a fit summary, original call and elapsed seconds. |
| `print.pmbfa_fits` | `R/methods.R` | Print the ordered model list and whole-request metadata. |
| `logLik.pmbfa` | `R/methods.R` | Return the accepted log likelihood with parameter and observation counts. |
| `BIC.pmbfa` | `R/methods.R` | Return the stored raw BIC, including for diagnostic nonconverged fits. |
| `nobs.pmbfa` | `R/methods.R` | Return the number of matrix observations. |
| `.pmbfa_models` | `R/models.R` | Enumerate the normal or skewed constraint models. |
| `.pmbfa_canonical` | `R/models.R` | Remove redundant cross-component flags when there is only one component. |
| `.pmbfa_stage_names` | `R/models.R` | Generate the constraint candidates for one adaptive stage. |
| `.pmbfa_nested_screen` | `R/models.R` | Compute the optional backward nested-model screening decision and bound. |
| `.pmbfa_psock` | `R/parallel.R` | Run grid tuples in socket workers using the installed package namespace. |
| `.pmbfa_factor_search` | `R/search.R` | Search factor dimensions, cache comparisons and rank results for each G. |
| `search_pmsbfa` | `R/search.R` | Search skewed mixture component counts and factor dimensions. |
| `search_pmbfa` | `R/search.R` | Search normal mixture component counts and factor dimensions. |
| `.pmbfa_search_fit` | `R/search.R` | Fit an internal search candidate without capturing expanded literal-data calls. |
| `.pmbfa_adaptive` | `R/selection.R` | Execute forward, backward or both staged constraint searches. |
| `.pmbfa_run_tuple` | `R/selection.R` | Reuse one initialization bank while fitting and caching constraint models. |
| `.pmbfa_elapsed` | `R/timing.R` | Calculate nonnegative elapsed wall-clock seconds from a recorded start. |
| `.pmbfa_timed_tuple` | `R/timing.R` | Add per-tuple elapsed time to the fitting-engine result. |
| `.pmbfa_record_result` | `R/timing.R` | Attach original-call and whole-request metadata without adding list entries. |
| `.pmbfa_integer` | `R/validation.R` | Validate and normalize positive or bounded integer requests. |
| `.pmbfa_max_factors` | `R/validation.R` | Find the conservative identifiable diagonal-factor rank limit. |
| `.pmbfa_labels` | `R/validation.R` | Map known labels and preserve unknown observations. |
| `.pmbfa_controls` | `R/validation.R` | Merge, validate and normalize fitting controls and legacy aliases. |

## Registered native entry points

| Internal wrapper | Purpose |
|---|---|
| `pmbfaEM` | Execute the complete normal AECM algorithm. |
| `pmsbfaEM` | Execute the complete skewed/general-family AECM algorithm. |
| `pmbfa_initializations` | Generate, screen and refine a reusable initialization bank. |
| `pmbfa_evaluate` | Evaluate a parameter state without fitting cycles. |
| `pmbfa_model_npar` | Count identifiable model parameters before fitting. |
| `pmbfa_posterior_moments` | Evaluate responsibilities and conditional mixing moments for diagnostics. |
| `pmbfa_core_version` | Identify the compiled numerical core. |

## C++ implementation functions

| Source module | Functions and responsibility |
|---|---|
| `matrix.cpp` | `transpose`, `multiply`, `dot`, `identity`: basic matrix operations; `chol`, `forward`, `backward`, `solve`: small positive-definite solves; `from_r`, `to_r`: data conversion; `Precision` constructor and `left`/`right`: low-rank inverse covariance operators. |
| `parameters.cpp` | `Controls` constructor: numerical controls; `check_name`, `parse`: input validation; `project`: impose model sharing; `count`: parameter count; `serialize`: return parameter lists; `blend`: damp a proposed cycle. |
| `distributions.cpp` | `logbk`, `dlogbk`: stable Bessel values/order derivatives; `gig`: conditional moments; `log_integral`: mixing integral; `shape_root`, `maximize`: safeguarded scalar updates; `update_family`: update concentration/shape parameters. |
| `expectation.cpp` | `quadratic_distance`: low-rank Mahalanobis calculation; `residual`: centered matrices; `expectation`: densities, posteriors and sufficient mixing moments. |
| `location.cpp` | `weighted_operator`: shared-A linear operator; `pcg_attempt`, `pcg`: warm-start solves with verified residuals and optional fallback; `update_location`: update locations and skewness. |
| `scale.cpp` | `sufficient`: stream factor sufficient statistics; `update_scale`: update row or column loadings and bounded diagonal residuals. |
| `em.cpp` | `fit`: complete accepted-cycle AECM driver, likelihood stopping, component-size guard and result construction. |
| `initialization.cpp` | `option_int`, `option_positive`: validate initialization controls; `random_index`: R-RNG index; `centers`, `nearest_assignment`, `repair_empty`: assignments; `make_params`: initial parameters; `fit_score`, `keep_candidate`: rank starts; `initialize`: candidate/refinement pipeline. |

The header contains shared structures and declarations; `api.cpp` contains thin Rcpp entry points. `RcppExports` files are generated, rather than maintained manually.
