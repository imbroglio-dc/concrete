Run `R CMD check` on this package and fix any problems found.

## Steps
1. First run `Rscript -e 'devtools::document()'` to ensure documentation is current
2. Run `Rscript -e 'devtools::check()'` and capture all output
3. Parse the output for ERRORs, WARNINGs, and NOTEs
4. For each issue found:
   - Identify the root cause in the source code
   - Apply the fix
   - Re-run the check to confirm it's resolved
5. Summarize what was found and fixed

## CRAN Standards
- Zero ERRORs and zero WARNINGs are mandatory
- NOTEs should be explained or resolved where possible
- Common issues: missing documentation, undeclared imports, examples exceeding time limits, non-portable code in src/

If the check passes cleanly, report the result and suggest any improvements that wouldn't show up as check failures (e.g., missing `@examples`, weak test coverage on new code).
