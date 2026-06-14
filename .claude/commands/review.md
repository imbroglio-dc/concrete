Review the current code changes for quality, correctness, and CRAN compliance.

## Steps
1. Run `git diff main` (or `git diff HEAD~1` if on main) to see what changed
2. For each changed file, review for:
   - **R code**: correctness, roxygen2 completeness, style, CRAN compliance, input validation
   - **C++ code**: memory safety, NA handling, Rcpp conventions, portability
   - **Tests**: adequate coverage of new/changed behavior, meaningful assertions
   - **Documentation**: accuracy, completeness, runnable examples
3. Run `Rscript -e 'devtools::check()'` to verify no new issues
4. Run `Rscript -e 'devtools::test()'` to verify tests pass
5. Provide a structured review:
   - **Must fix**: Bugs, CRAN violations, incorrect statistics
   - **Should fix**: Missing tests, incomplete docs, performance issues
   - **Suggestions**: Style improvements, refactoring opportunities

Use the r-reviewer agent for the code review analysis.
