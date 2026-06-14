---
name: r-package-conventions
description: R package development patterns and conventions. Use when writing R functions, editing DESCRIPTION/NAMESPACE, managing imports, or dealing with R CMD check issues.
---

# R Package Development Conventions for `concrete`

## Function Structure
```r
#' Compute Cause-Specific Absolute Risk
#'
#' @description Computes TMLE-targeted cumulative incidence estimates
#'   for specified cause(s) at given time point(s).
#'
#' @param ConcreteArgs A `ConcreteArgs` object from [formatArguments()].
#' @param Interventions Named list of intervention functions.
#' @param Times Numeric vector of evaluation times.
#'
#' @return A `ConcreteEst` object (S3) containing risk estimates, influence
#'   curve values, and inference results.
#'
#' @examples
#' \donttest{
#' data(pbc)
#' # ... example code ...
#' }
#'
#' @seealso [formatArguments()], [plot.ConcreteEst()]
#' @family estimation
#' @export
computeRisk <- function(ConcreteArgs, Interventions, Times) {
  # Input validation first
  checkmate::assert_class(ConcreteArgs, "ConcreteArgs")
  checkmate::assert_numeric(Times, lower = 0, any.missing = FALSE)

  # Implementation...
}
```

## DESCRIPTION Management
- `Imports`: packages whose functions you call (data.table, Rcpp, survival, etc.)
- `Suggests`: packages used only in tests, vignettes, or examples (testthat, knitr, covr)
- `LinkingTo`: Rcpp (for C++ compilation)
- Version bounds: add only when a specific feature is needed: `data.table (>= 1.14.0)`
- Alphabetize entries within each field

## Common R CMD check Fixes
- "no visible binding for global variable 'X'" → Add `X <- NULL` at top of function (data.table columns) or use `.data$X`
- "Namespace in Imports field not imported from" → Actually use the package or move to Suggests
- "S3 method documentation" → Ensure `@method` and `@export` are both present for S3 methods
- "'library' or 'require' call in package code" → Remove; use NAMESPACE imports
