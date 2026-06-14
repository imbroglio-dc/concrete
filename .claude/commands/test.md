Run the testthat suite and address any failures.

## Steps
1. Run `Rscript -e 'devtools::test()'` and capture output
2. If all tests pass, report the summary (number of tests, skips, time)
3. If any tests fail:
   - Read the failing test file to understand what's being asserted
   - Read the relevant source code
   - Determine if the bug is in the source code or the test expectation
   - Fix the root cause (prefer fixing source code unless the test expectation is clearly wrong)
   - Re-run to confirm the fix
4. Report what was found and fixed

Optional argument: `$ARGUMENTS` — if provided, run only matching test files:
```
Rscript -e 'devtools::test(filter = "$ARGUMENTS")'
```
