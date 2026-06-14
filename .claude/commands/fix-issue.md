Fix GitHub issue: $ARGUMENTS

## Steps
1. Read the issue description from the provided issue number or description
2. Reproduce the problem:
   - Write a minimal script that triggers the bug or demonstrates the missing feature
   - Run it to confirm the current behavior
3. Identify the root cause in the source code
4. Implement the fix or feature:
   - Make minimal, focused changes
   - Follow the coding standards in CLAUDE.md
5. Write a regression test that would have caught this issue
6. Run `Rscript -e 'devtools::test()'` to ensure nothing else broke
7. Run `Rscript -e 'devtools::check()'` to ensure no new check issues
8. Update documentation if the fix changes user-facing behavior
9. Summarize: what the issue was, what caused it, what was changed, and what tests were added

## Git workflow
- Create a descriptive branch name: `fix/<issue-number>-<brief-description>`
- Make atomic commits with clear messages
- Note the issue number in commit messages: `Fix #<n>: <description>`
