Write testthat tests for the specified function or file: $ARGUMENTS

## Steps
1. Read the source code for the target function(s)
2. Read any existing tests in `tests/testthat/` that cover this code
3. Identify untested behaviors, edge cases, and code paths
4. Write new `test_that()` blocks covering:
   - **Basic correctness**: known-input → known-output with hand-verifiable values
   - **Edge cases**: empty input, single observation, all-censored data, missing values
   - **Type checking**: wrong input types produce informative errors
   - **Dimensions**: output has expected length/shape for varying input sizes
   - **Statistical properties**: for stochastic methods, test on simulated data with large n and generous tolerances
   - **Regression**: if fixing a bug, add a test that would have caught it
5. Run the new tests: `Rscript -e 'devtools::test(filter = "<test_file>")'`
6. Fix any issues until all new tests pass
7. Report what was tested and any gaps that remain

## Conventions
- File naming: `tests/testthat/test-<topic>.R` matching the R source file
- Use `withr::local_seed()` for reproducible randomness, not bare `set.seed()`
- Use `expect_equal(..., tolerance = ...)` for floating-point comparisons
- Use `expect_error(..., regexp = "...")` to verify error messages
- Keep individual tests fast (< 2 seconds each); mark slow tests with `skip_on_cran()`
- For functions that use `data.table`, test that inputs are not modified by reference
