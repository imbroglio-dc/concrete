Plan and implement a new feature: $ARGUMENTS

## Workflow

### Phase 1: Plan
1. Understand the feature request and its statistical/methodological basis
2. Survey the existing codebase to identify where the feature fits
3. Outline the implementation plan:
   - Which files need to change
   - New functions needed (name, signature, purpose)
   - C++ additions if performance-critical
   - Changes to the user-facing API
   - Backward compatibility considerations
4. **Present the plan and wait for approval before proceeding**

### Phase 2: Implement
1. Create a feature branch: `feature/<description>`
2. Write the implementation following CLAUDE.md coding standards
3. Add complete roxygen2 documentation for any new/changed exports
4. Run `Rscript -e 'Rcpp::compileAttributes()'` if C++ was touched
5. Run `Rscript -e 'devtools::document()'`

### Phase 3: Test & Validate
1. Write testthat tests covering the new functionality
2. Run `Rscript -e 'devtools::test()'` — all tests must pass
3. Run `Rscript -e 'devtools::check()'` — no new issues
4. Summarize what was built and suggest vignette updates if appropriate
