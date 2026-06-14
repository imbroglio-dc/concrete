Prepare this package for CRAN submission. Run through the full pre-submission checklist.

## Checklist (execute each step)

### 1. Version & Metadata
- Read `DESCRIPTION` and verify: version bump is appropriate (major.minor.patch), all authors listed, License is valid, URL and BugReports point to GitHub
- Verify `cran-comments.md` exists and is current (create/update if not)
- Check that `NEWS.md` has an entry for this version

### 2. Documentation
- Run `Rscript -e 'devtools::document()'`
- Verify all exported functions have complete roxygen2 docs
- Verify all `@examples` run without error
- Build vignettes: `Rscript -e 'devtools::build_vignettes()'`

### 3. Tests
- Run `Rscript -e 'devtools::test()'` — all must pass
- Check coverage: `Rscript -e 'covr::package_coverage()'` — flag any critical gaps

### 4. R CMD check
- Run `Rscript -e 'devtools::check(args = c("--as-cran"))'`
- Must produce: 0 errors, 0 warnings
- Document any NOTEs in `cran-comments.md`

### 5. Multi-platform (guidance)
- Remind to check via `rhub::check_for_cran()` or `devtools::check_win_devel()` before actual submission
- Check that compiled code (src/) is portable (no platform-specific headers)

### 6. Final Cleanup
- No leftover `browser()`, `debug()`, `cat()` calls in R code
- No `Depends` that should be `Imports`
- `.Rbuildignore` excludes dev files (`.github/`, `paper/`, etc.)
- `inst/CITATION` is up to date if present

### 7. Submission
- Generate tarball: `R CMD build .`
- Provide the `devtools::submit_cran()` command but DO NOT run it — leave that to the human

Report a summary of pass/fail for each section and any issues found.
