Lint the R codebase and fix style issues.

## Steps
1. Run `Rscript -e 'lintr::lint_package()'` and capture all results
2. Group issues by severity and type
3. Auto-fix straightforward issues:
   - Trailing whitespace, missing terminal newlines
   - Inconsistent spacing around operators
   - `T`/`F` → `TRUE`/`FALSE`
   - Trailing commas, semicolons
4. For issues requiring judgment (line length, variable naming, complexity):
   - Fix where the change is unambiguous and safe
   - Flag others for human review
5. Re-run the linter to confirm improvements
6. Report what was fixed and what remains

If `$ARGUMENTS` is provided, lint only that file.
