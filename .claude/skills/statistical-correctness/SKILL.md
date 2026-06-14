---
name: statistical-correctness
description: Verify statistical computation correctness for the TMLE estimator. Use when reviewing or changing getEIC.R, doTmleUpdate.R, getPropScore.R, getHazEstimate.R, getInitialEstimate.R, or rcpp_getCleverCovariate.cpp.
user-invocable: false
---

# Statistical Correctness Verification for `concrete`

## Required Invariants

After any computation change, verify all of the following hold:

| Invariant | How to Check |
|---|---|
| `h_l(t,i) ≥ 0` | `stopifnot(all(unlist(Estimates$Hazards) >= 0))` |
| `S(t,i) ∈ (1e-12, 1]` | `stopifnot(all(Estimates$EvntFreeSurv >= 1e-12), all(Estimates$EvntFreeSurv <= 1))` |
| `F_j(τ,i) ∈ [0,1]` | Check cumulative incidence columns in ConcreteEst |
| `Σ_j F_j(τ,i) ≤ 1` | Sum across cause-specific CIs at each time, check ≤ 1 |
| `mean(IC[j,τ]) ≈ 0` | After convergence: `abs(SummEIC$PnEIC) <= SummEIC$seEIC / (sqrt(n) * log(n))` |
| Nuisance weights ≥ `MinNuisance` | `stopifnot(all(NuisanceWeight >= MinNuisance, na.rm = TRUE))` |

## Convergence Criterion

```r
# Converged when TRUE for all (j, tau):
abs(SummEIC$PnEIC) <= SummEIC$seEIC / (sqrt(n) * log(n))
```

## C++ vs R Reference Comparison

When modifying `computeIC()` or `updateHazardsCpp()`, verify against the R reference:

```r
library(concrete)
set.seed(42)
n <- 100
dt <- data.table::data.table(
  id = 1:n,
  time = rexp(n, 0.5),
  event = sample(0:2, n, replace = TRUE, prob = c(0.3, 0.4, 0.3)),
  trt = rbinom(n, 1, 0.5),
  x1 = rnorm(n),
  x2 = rnorm(n)
)
args <- formatArguments(dt, "time", "event", "trt",
                        TargetTime = quantile(dt$time, 0.6),
                        TargetEvent = 1:2)
est  <- getInitialEstimate(args$Data, args$Model, args$CVFolds,
                            args$MinNuisance, args$TargetEvent,
                            args$TargetTime, args$Regime,
                            ReturnModels = FALSE)

# C++ path (via getEIC)
cpp_result <- concrete:::getEIC(est, args$Data, args$Regime,
                                 args$TargetEvent, args$TargetTime,
                                 args$MinNuisance, GComp = FALSE)

# Compare PnEIC values — should be identical to machine precision
print(cpp_result[[1]]$SummEIC)
```

## Numerical Red Flags

Watch for these in code review:

- **Division by near-zero survival**: any `F / S` or `h / S` without a floor guard → should be `max(S, 1e-12)`
- **Log of near-zero hazard**: `log(h)` where h can be 0 → guard or use `log1p`
- **Nuisance weight explosion**: if `MinNuisance` is too small relative to n, weights blow up; default is `5 / (sqrt(n) * log(n))`
- **Non-monotone cumulative incidence**: `F_j(t)` must be non-decreasing in t; a decreasing value signals a bug in the hazard-to-CI integration
- **`hTotal` collapse error**: the C++ optimization in `computeIC()` that collapses the L-loop via `hTotal` is only valid because the integral is linear in the cause-specific hazards — do not change this structure without re-deriving

## Tests to Run After Any Computation Change

```bash
# C++ reference comparison (most targeted)
Rscript -e 'devtools::test(filter = "getIC", reporter = "minimal")'

# Numerical stability edge cases
Rscript -e 'devtools::test(filter = "numerical-stability", reporter = "minimal")'

# End-to-end TMLE
Rscript -e 'devtools::test(filter = "doConcrete", reporter = "minimal")'
```

## Reference

Methodology: Rytgaard et al. (2023) <doi:10.1111/biom.13856> and Rytgaard & van der Laan (2023) <doi:10.1007/s10985-022-09576-2>

Preprint: arXiv:2310.19197
