# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

`concrete` is an R package implementing continuous-time Targeted Minimum Loss-Based Estimation (TMLE) for survival analysis with competing risks. It computes covariate-adjusted marginal cumulative incidence estimates in right-censored survival settings.

## Build and Development Commands

```bash
# Install dependencies and build
R CMD INSTALL .

# Run all tests
R -e "devtools::test()"

# Run a single test file
R -e "testthat::test_file('tests/testthat/test-doConcrete.R')"

# Check package (full R CMD check)
R -e "devtools::check()"

# Generate documentation from roxygen comments
R -e "devtools::document()"

# Build and install from source
R -e "devtools::install()"
```

## Architecture

### Core Workflow (Main Pipeline)

The package follows a three-step estimation pipeline:

1. **`formatArguments()`** (`R/formatArguments.R`) - Validates and formats input data, creates a `ConcreteArgs` object containing:
   - Formatted data.table with event times, types, treatment, covariates
   - Intervention specifications via `makeITT()` for intent-to-treat or custom functions
   - Cross-validation fold setup via `origami` package
   - Model specifications for propensity scores and hazards

2. **`doConcrete()`** (`R/doConcrete.R`) - Runs the TMLE estimation, returns `ConcreteEst` object:
   - Calls `getInitialEstimate()` for nuisance parameter estimation
   - Calls `getEIC()` to compute efficient influence curves
   - Calls `doTmleUpdate()` for one-step TMLE updates until convergence

3. **`getOutput()`** (`R/getOutput.R`) - Extracts final estimates, returns `ConcreteOut` object:
   - Computes absolute risks, risk differences (RD), relative risks (RR)
   - Provides both pointwise and simultaneous confidence intervals

### Internal Estimation Components

- **`getInitialEstimate.R`** - Orchestrates initial nuisance parameter estimation
- **`getPropScore.R`** - Treatment propensity score estimation using SuperLearner
- **`getHazEstimate.R`** - Cause-specific hazard estimation using Cox regression ensembles
- **`getEIC.R`** - Efficient influence curve computation for all target estimands
- **`doTmleUpdate.R`** - One-step TMLE fluctuation updates

### C++ Performance Code

`src/rcpp_getCleverCovariate.cpp` contains RcppArmadillo implementations for:
- `getCleverCovariate()` - Clever covariate calculation for EIC
- `getHazLS()` - Hazard integration helper
- `updateHazardsCpp()` - TMLE hazard update step

## Key Data Structures

- **`ConcreteArgs`** - Environment object holding validated inputs (class assigned in `formatArguments()`)
- **`ConcreteEst`** - List of estimates per intervention with hazards, survival curves, EICs (class assigned in `doConcrete()`)
- **`ConcreteOut`** - data.table of final estimates with CIs (class assigned in `getOutput()`)

## Dependencies

Core: `data.table`, `survival`, `SuperLearner`, `origami`, `Rcpp`, `RcppArmadillo`

Suggested for extended functionality: `glmnet`, `xgboost`, `ranger`, `ggplot2`

## Testing Approach

Tests use `testthat` and run on subsets of `survival::pbc` data. Test files cover the main pipeline functions. Tests can take several minutes due to TMLE iterations.
