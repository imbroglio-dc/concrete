Analyze test coverage and fill gaps for: $ARGUMENTS (or the whole package if blank).

## Steps
1. Run `Rscript -e 'cov <- covr::package_coverage(); print(cov); covr::report(cov, file = "/tmp/covr-report.html")'`
2. Identify the lowest-coverage files and functions
3. For each significant gap (< 80% line coverage):
   - Read the untested code paths
   - Write targeted testthat tests that exercise them
   - Place tests in the appropriate `tests/testthat/test-*.R` file
4. Re-run coverage to confirm improvement
5. Report before/after coverage percentages for each file touched

## Priorities
- Exported functions should have ≥ 90% coverage
- Error-handling branches (input validation, edge cases) are high-value test targets
- C++ code coverage via R-level tests — ensure the R wrappers exercise all C++ branches
