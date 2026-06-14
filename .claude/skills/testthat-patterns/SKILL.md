---
name: testthat-patterns
description: Testing patterns for statistical R packages using testthat. Use when writing, running, or debugging tests, or when improving test coverage.
---

# testthat Patterns for `concrete`

## File Setup
```r
# tests/testthat/test-estimation.R
test_that("TMLE estimate matches known value on simple case", {
  withr::local_seed(42)

  # Create minimal reproducible data
  n <- 200
  dt <- data.table::data.table(
    id = 1:n,
    time = rexp(n),
    event = sample(0:2, n, replace = TRUE),
    trt = rbinom(n, 1, 0.5),
    x1 = rnorm(n), 
    x2 = rbinom(n, 1, 0.5) # glmnet needs at least 2 covariates
  )

  # Run estimation
  result <- someFunctionFromConcrete(dt)

  # Test structure

  expect_s3_class(result, "ExpectedClass")
  expect_named(result, c("estimates", "influence_curves", "inference"))

  # Test numerical correctness with tolerance
  expect_equal(result$estimates[1], 0.312, tolerance = 0.05)

  # Test dimensions
  expect_length(result$estimates, expected_length)
})
```

## Testing Statistical Methods
```r
test_that("TMLE is consistent (converges to truth with large n)", {
  skip_on_cran()  # slow test
  withr::local_seed(123)

  # Generate data from known DGP where truth is analytically computable
  n <- 5000
  # ... simulate data where true risk = 0.3 ...

  result <- estimateRisk(...)
  expect_equal(result$risk, 0.3, tolerance = 0.05)
})

test_that("confidence intervals have correct coverage", {
  skip_on_cran()
  withr::local_seed(456)

  covers <- replicate(200, {
    # ... simulate and estimate ...
    truth >= ci_lower && truth <= ci_upper
  })

  # 95% CI should cover ~95% of the time (allow sampling variability)
  expect_gt(mean(covers), 0.88)
  expect_lt(mean(covers), 0.99)
})
```

## Edge Case Patterns
```r
test_that("handles single-event data (no competing risks)", {
  # Only event type 1 and censoring (0), no type 2
  # ...
})

test_that("errors informatively on empty data", {
  expect_error(someFunction(data.table()), regexp = "at least")
})

test_that("data.table input is not modified by reference", {
  dt <- data.table::data.table(x = 1:5, y = rnorm(5))
  dt_copy <- data.table::copy(dt)
  someFunction(dt)
  expect_identical(dt, dt_copy)
})
```

## Key Expectations
- `expect_equal(x, y, tolerance = tol)` — floating point (use tol = 1e-6 to 1e-2 depending on context)
- `expect_identical(x, y)` — exact match (integers, strings, structure)
- `expect_error(expr, regexp = "pattern")` — verify error message
- `expect_warning(expr, regexp = "pattern")` — verify warning
- `expect_message(expr, regexp = "pattern")` — verify message()
- `expect_s3_class(obj, "ClassName")` — check S3 class
- `expect_true(all(result >= 0 & result <= 1))` — range check for probabilities
