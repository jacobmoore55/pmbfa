# Mathematical implementation notes

This note derives the parameter counts and streamed covariance updates from the
supplied `main.tex`, rather than treating the existing R implementation as a
numerical reference. Mathematical statements below assume generic, locally
identifiable factor-analyzer covariance parametrizations. Rank-deficient
loadings and overlarge factor ranks need separate treatment.

## Constraint interpretation and parameter counts

For a covariance block with observed dimension `d`, factor rank `k`, and three
letters `abc`:

- `a=C` shares its loading matrix across all components; `a=U` separates it.
- `b=C` shares its positive diagonal residual matrix across all components.
- `c=C` restricts each residual matrix to be isotropic; `c=U` permits its `d`
  diagonal entries to differ. This letter does not imply cross-group sharing.

Let `Lmult = 1` for `a=C`, otherwise `G`; `Dmult = 1` for `b=C`, otherwise `G`;
and `Dwidth = 1` for `c=C`, otherwise `d`. The block dimension is

    Lmult * (d*k - k*(k-1)/2) + Dmult * Dwidth.

The subtraction `k*(k-1)/2` removes continuous orthogonal rotations separately
for every independently fitted loading matrix. The manuscript tables at
`main.tex:231–264` give exactly these eight dimensions.

For the complete mixture, before the Kronecker scale adjustment:

    dimension = (G-1) + G*n*p + row_dimension + column_dimension
                + skew_dimension + shape_dimension.

`skew_dimension` is zero for Normal, `n*p` for shared skewness, and `G*n*p`
otherwise. Under the existing GH convention, the shared concentration letter
shares `omega` only, with `lambda_g` always fitted separately:

| Family | Final concentration letter C | Final concentration letter U |
| --- | ---: | ---: |
| Normal | 0 | 0 |
| ST | 1 | G |
| NIG | 1 | G |
| VG | 1 | G |
| GH | G+1 | 2G |

These dimensions presume these shape parameters are estimated, not fixed.
The existing `count_pars` omits the mixing-proportion contribution `G-1`.

### Kronecker scale dimension

Write `S_g = L_g L_g' + D_g` and `P_g = H_g H_g' + E_g`. The observed
component scale is `P_g tensor S_g`. Its componentwise scale invariance is

    L_g -> sqrt(c_g)*L_g,   D_g -> c_g*D_g,
    H_g -> H_g/sqrt(c_g),   E_g -> E_g/c_g,

for `c_g>0`. If a nonzero loading matrix is shared across groups, preservation
of that equality forces all `c_g` to agree. Sharing a positive residual matrix
likewise forces all `c_g` to agree. Isotropy by itself imposes no equality among
the `c_g`.

Consequently, for positive factor ranks and generic locally identifiable
factor parametrizations, subtract

    G, if both row and column blocks start with UU;
    1, otherwise.

Thus only the four covariance combinations from
`{UUC, UUU} x {UUC, UUU}` have `G` independent scale gauges. A shared loading
matrix alone is enough to reduce the gauge dimension to one, even when both
residual matrices are heterogeneous across components. For zero factors, a
shared empty loading matrix is vacuous and must not tie scale gauges; positive
shared residual matrices still do so.

Why this counts *all* generic scale redundancies: equality of positive definite
Kronecker products implies their factors differ by a positive scalar. Generic
local identifiability of each factor-analyzer representation then permits only
the stated scale change plus the loading rotations already removed above.

The normalization discussion at `main.tex:732` is too broad: independently
normalizing each component in a partially shared model can violate its shared
loading/residual constraints. Normalization is not needed to compute the
dimension or to run AECM.

If normalization is used, apply one common `c` in shared models, or one `c_g`
per component in fully independent models. With a residual variance lower
bound `f>0`, admissibility requires

    max_{g,j}(f / D_g[j]) <= c <= min_{g,j}(E_g[j] / f)

within the group of components being scaled. The current state makes `c=1`
feasible. A desired normalization outside this interval must be omitted or
clamped to the interval; normalizing and then flooring the entries changes the
fitted component density. The safest default is no automatic normalization.

### Numerical verification

`verify_mathematics.R` builds an independent map from raw loadings and
residuals to vectorized symmetric Kronecker scales. Central-difference
Jacobians were checked for all 64 covariance models at each of:

- `n=p=4`, `q=r=1`, `G=3`;
- `n=p=5`, `q=r=2`, `G=2`.

All 128 generic Jacobian ranks agree with the dimension formula above, including
rotation discounts in the second configuration. This is numerical evidence in
addition to the generic scale argument, not a proof for singular or
overparameterized cases. The Ledermann-style rank limit currently used in
`max_factors` is relevant for unrestricted diagonal factor analyzers; formulas
must not be represented as exact identifiable dimensions for arbitrary ranks.

## Streamed row-factor sufficient statistics

For one component, let:

    R_i = X_i - M,   z_i = responsibility,
    a_i = E[W_i | X_i,g],   b_i = E[1/W_i | X_i,g],
    J = (Psi + Delta Delta')^{-1},
    V = (I + Lambda' Sigma^{-1} Lambda)^{-1},
    B = V Lambda' Sigma^{-1},   N_g = sum_i z_i.

Here `Sigma` and `Psi` denote the diagonal residual matrices; `J` uses the full
opposite covariance. The only large scatter implicit in the row CM step is

    C = sum_i z_i * [ b_i R_i J R_i'
                       - R_i J A' - A J R_i' + a_i A J A' ].

One does not need to allocate this `n x n` matrix. The necessary statistics are

    cdiag = diag(C),
    T = C B',
    Q = N_g * p * V + B T.

Their sizes are `n`, `n x q`, and `q x q`. Symmetrize `Q` to remove roundoff
asymmetry before an SPD solve.

### Stable scatter representation

For `b_i>0`, define

    H_i = sqrt(b_i) R_i - A/sqrt(b_i),
    kappa = sum_i z_i * (a_i - 1/b_i).

Then the same scatter is

    C = sum_i z_i H_i J H_i' + kappa A J A'.

The conditional-moment inequality `a_i*b_i >= 1` makes `kappa>=0` in exact
arithmetic. This representation avoids subtracting several large scatter
matrices. If only `a_sum=sum_i z_i*a_i` was retained, then
`kappa=a_sum-sum_i z_i/b_i`; tiny negative roundoff can be clamped with a
relative tolerance, but a materially negative value indicates invalid moments.

For each observation, stream:

    HJ = H_i J,
    cdiag += z_i * rowSums(H_i * HJ),
    T += z_i * HJ * (B H_i)'.

Finally add the same quantities for `A`, weighted by `kappa`. Apply `J` through
its diagonal-plus-low-rank inverse representation, without allocating a dense
`p x p` inverse. Normal models simply have `A=0`, `a_i=b_i=1`.

### Loading and residual updates under all constraints

For component-specific loadings:

    Lambda_new,g = T_g Q_g^{-1}.

For shared loadings, the general rowwise update for row `j` is

    Lambda_new[j,] = (sum_g T_g[j,] / Sigma_old,g[j])
                      (sum_g Q_g / Sigma_old,g[j])^{-1}.

This simplifies to one unweighted aggregate solve when the residual diagonal
is shared (`CCC`, `CCU`), and one scalar-weighted aggregate solve for `CUC`.
Only `CUU` needs separate small solves for each observed row. Weighting uses
the old residual variances for this CM step.

For each component define the residual numerator vector

    u_g = cdiag_g - 2*rowSums(Lambda_new,g * T_g)
            + rowSums((Lambda_new,g Q_g) * Lambda_new,g).

This is the general valid expression for all loading constraints. Do not
universally replace it by `cdiag-diag(Lambda_new T')`: that simplification uses
the componentwise loading normal equation and fails for shared loadings with
heterogeneous residuals.

Residual CM updates are:

- Unique anisotropic: `Sigma_new,g = u_g / (p*N_g)`.
- Unique isotropic: repeat `sum(u_g)/(n*p*N_g)` down the diagonal.
- Shared anisotropic: every component gets `sum_g u_g/(p*N)`.
- Shared isotropic: every component gets `sum_g sum(u_g)/(n*p*N)`.

A fixed lower bound on each variance is imposed by taking the elementwise
maximum with the bound after the corresponding unconstrained variance
maximization. For isotropic/shared cases, aggregate first, then floor the
result to preserve the specified equalities.

The column CM step is exactly the transpose analogue: exchange
`n <-> p`, `q <-> r`, `Lambda <-> Delta`, `Sigma <-> Psi`, and use `X_i'`,
`M'`, and `A'`. Its opposite inverse must use the newly updated row covariance.

`verify_mathematics.R` also verifies 20 independent random skewed/Normal
examples against explicitly constructed dense scatters and residual matrices.
All streamed statistics and residual diagonals agreed within `1e-10`.

## Outer and inner expectations

The supplied manuscript and original implementation compute `z_i`, `a_i`,
`b_i`, and `c_i` once per outer iteration. The row-factor expectations then use
the newly updated location and skewness, and the column-factor expectations
use the newly updated row covariance (`main.tex:573` explicitly calls for this
refresh). The old R loop confirms that ordering at `MVBFA.R:2350–2365`.

There is no need to recompute the expensive GIG/Bessel posterior moments
between covariance CM stages. Holding the outer missing-data expectations
fixed defines an expected complete-data objective; the factor-analyzer
updates act as inner EM/conditional updates that increase this same objective.
The stable scatter representation above exposes that positive semidefinite
weighted Gaussian scatter explicitly. This preserves the manuscript's
algorithm while permitting streaming and avoids introducing an unrequested
different schedule of outer posterior updates.

## Other mathematical conventions

- Preserve the manuscript's lower-is-better BIC:
  `-2*loglik + dimension*log(N)`.
- Existing ICL is the MAP-classification approximation:
  `BIC - 2*sum_i log(max_g posterior[i,g])`, not the soft entropy variant.
- For known observation labels, the observed objective must use
  `log(pi_label)+log_density_label` for that observation, not the unlabelled
  mixture log-sum followed only by overwriting its responsibilities.
- The skew-t definition permits `nu>0`, whereas the old R routine uses a floor
  `4.1`. With `X=M+W*A+sqrt(W)*V`, a nonzero skewness matrix requires `nu>4`
  for finite second moments. A numerical/root-finding fallback must not be
  silently confused with this statistical restriction.

## Native-kernel audit and practical numerical limits

An independent dense likelihood calculation at moderate, nondegenerate
parameters was compared to the compiled kernel for Normal, ST, GH, NIG, and
VG, including two known class labels. Absolute differences were at most
`2.3e-13`. The repeatable check is `verify_density.R`.

The Bessel routine uses R's scaled `K_nu` when representable, with a leading
uniform asymptotic and first Debye correction as a fallback for large orders.
That fallback is an approximation; differences of extremely large log-Bessel
values lose precision when obtaining moment ratios and order derivatives.
On a 4-by-5 skew-t example, one complete CM cycle succeeded for starting
`nu=10, 1e3, 1e5, 1e7`, while starting `nu=1e9` and `1e12` returned a safe
`numerical_error` with zero accepted iterations because the resulting GIG
moments violated `E[W] E[1/W] >= 1`. The numerical status and last valid state
are preserved. Very large configured bounds are therefore not a promise of
accurate evaluation throughout that range; additional high-order asymptotic
moment methods would be needed to extend it reliably.

## Corrected Normal shared-loading residual update

The old Normal implementation's `norm_sig_expr` and `norm_psi_expr` simplify
the residual numerator to `diag(C)-diag(Lambda_new*T')` (and its column
analogue). This is valid when each component's loading update separately
satisfies `Lambda_new,g*Q_g=T_g`. A common loading matrix with heterogeneous
component residual variances instead satisfies a weighted sum of these normal
equations. The simplified identity need not hold for each component, so that
implementation can update Normal residual variances incorrectly in `CUU`
models. The new kernel uses the complete residual expression documented above
for both Normal and skewed cases.

`test_multifactor.R` independently constructs full covariance and scatter
matrices in R, performs two full Gaussian AECM cycles with per-row shared
loading solves and the complete residual expression, and compares the entire
likelihood curve, final posterior, and loading/residual parameter arrays with
the native kernel. For `n=6,p=7,q=r=2,N=90,G=2` and covariance models
`CUU-CUU`, `CUU-UUU`, and `UUU-CUU`, all likelihood curves agree within
`1.819e-12`, with parameter and posterior comparisons also passing their
strict tolerances. The three corresponding skew-t regressions still agree
with the old skewed implementation within `3.638e-12`. Differences from the
old Normal trajectory in these cases are therefore an intentional correction,
not accepted numerical drift.

## Follow-up: common covariance shortcut and mean diagnostics

Write `ng=sum(z)`, `sa=sum(z*a)`, `sb=sum(z*b)`, where
`a=E[W|X,g]`, `b=E[1/W|X,g]`, `T=sum(z*X)`, and `B=sum(z*b*X)`.
Eliminating the component location from the common-skewness normal equation
leaves

```
sum_g (sa_g - ng_g^2/sb_g) Sigma_g^-1 A Psi_g^-1
  = sum_g Sigma_g^-1 (T_g - ng_g B_g/sb_g) Psi_g^-1.
```

When both scale matrices are shared (row and column loading/sharing flags
both CC), the common precisions cancel. The native implementation uses the
exact direct solution

```
A = sum_g (T_g - ng_g B_g/sb_g) / sum_g (sa_g - ng_g^2/sb_g)
M_g = (B_g - ng_g A)/sb_g.
```

This removes an iterative linear solve without changing the CM objective. The
four combinations of common isotropic/diagonal residuals were checked for all
four skew families against independently computed conditional moments.

For component mean diagnostics, set `abar_g=sa_g/ng_g` and
`score_A_g=T_g-ng_g*M_g-sa_g*A_g`. The exact algebraic identity is

```
(M_g + E_prior[W_g]*A_g) - T_g/ng_g
  = (E_prior[W_g] - abar_g)*A_g - score_A_g/ng_g.
```

It applies even before convergence. At a stationary solution with separate A,
`score_A_g=0`; it does not require `E_prior[W_g]=abar_g`. With shared A, only
`sum_g Sigma_g^-1 score_A_g Psi_g^-1=0` is required. Hard-cluster empirical
means add the effect of replacing soft responsibilities with assignments.
The diagnostic helper reports each contribution and verifies the identity.

The density and conditional-moment calculations now reuse their identical
log-Bessel normalizer within an observation/component. This removes duplicate
work and adds no approximation. The isolated reuse change gave bitwise equal
parameters, likelihood curves, and posteriors in checks of all five families
and a full-size Olivetti cycle.
