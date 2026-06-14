---
name: rcpp-patterns
description: C++ and Rcpp conventions for R packages. Use when writing or editing C++ code in src/, debugging compilation errors, or optimizing R code with C++.
---

# Rcpp Patterns for `concrete`

## File Structure
- Source files: `src/*.cpp`
- Auto-generated: `src/RcppExports.cpp` and `R/RcppExports.R` — never edit these
- Regenerate after any change to `// [[Rcpp::export]]` signatures:
  ```bash
  Rscript -e 'Rcpp::compileAttributes()'
  ```

## Template
```cpp
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector myFunction(NumericMatrix x, NumericVector weights) {
  int n = x.nrow();
  int p = x.ncol();
  NumericVector result(n);

  for (int i = 0; i < n; i++) {
    // Check for user interrupt in long loops
    if (i % 1000 == 0) Rcpp::checkUserInterrupt();

    double sum = 0.0;
    for (int j = 0; j < p; j++) {
      // Guard against NA
      if (NumericVector::is_na(x(i, j)) || NumericVector::is_na(weights[j])) {
        result[i] = NA_REAL;
        break;
      }
      sum += x(i, j) * weights[j];
    }
    result[i] = sum;
  }
  return result;
}
```

## Key Rules
- Always guard against `NA`/`NaN`/`Inf` — R passes them silently into C++
- Use `Rcpp::checkUserInterrupt()` in long loops so users can Ctrl-C
- Use `Rcpp::stop("message")` for errors — it throws an R-level error, not a C++ exception
- Prefer `NumericVector`/`NumericMatrix` over raw `SEXP`
- Never use `std::cout` — use `Rcpp::Rcout` if you must print (but prefer R-level messaging)
- Memory: Rcpp objects are garbage-collected by R; don't `delete` them
- CRAN portability: no platform-specific headers, no C++17+ features unless guarded in SystemRequirements
