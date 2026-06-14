---
name: r-github-actions
description: GitHub Actions CI/CD for R packages. Use when setting up, debugging, or updating CI workflows, or when preparing for CRAN submission.
---

# GitHub Actions for R Packages

## Standard R CMD check Workflow
Use the r-lib/actions collection. The canonical setup:

```yaml
# .github/workflows/R-CMD-check.yaml
name: R-CMD-check
on:
  push:
    branches: [main, develop]
  pull_request:
    branches: [main]

jobs:
  R-CMD-check:
    runs-on: ${{ matrix.config.os }}
    name: ${{ matrix.config.os }} (R ${{ matrix.config.r }})
    strategy:
      fail-fast: false
      matrix:
        config:
          - {os: ubuntu-latest, r: 'release'}
          - {os: ubuntu-latest, r: 'devel'}
          - {os: macOS-latest, r: 'release'}
          - {os: windows-latest, r: 'release'}

    steps:
      - uses: actions/checkout@v4
      - uses: r-lib/actions/setup-r@v2
        with:
          r-version: ${{ matrix.config.r }}
      - uses: r-lib/actions/setup-r-dependencies@v2
        with:
          extra-packages: any::rcmdcheck
      - uses: r-lib/actions/check-r-package@v2
```

## Additional Useful Workflows
- **test-coverage.yaml**: Run covr and upload to Codecov
- **pkgdown.yaml**: Build and deploy package documentation site
- **lint.yaml**: Run lintr on PRs

## Debugging CI Failures
- Read the full log; R CMD check errors are often near the bottom
- "package not available" → add to `extra-packages` or check CRAN availability
- Compilation failures on Windows → check for POSIX-only headers in src/
- Timeout → mark slow tests with `skip_on_cran()` and check example runtimes
