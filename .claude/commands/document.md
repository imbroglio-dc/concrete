Update roxygen2 documentation for this package.

## Steps
1. Scan all R files in `R/` for exported functions (those with `@export` tags)
2. Check that every exported function has complete roxygen2 documentation:
   - `@title` — concise one-line title
   - `@description` — what it does and when to use it
   - `@param` for every parameter — type, meaning, default behavior
   - `@return` — what the function returns, including class and structure
   - `@examples` — runnable example (wrap slow ones in `\donttest{}`)
   - `@export`
3. For any function missing documentation or with incomplete docs, write or complete the roxygen block
4. Run `Rscript -e 'devtools::document()'` to regenerate `man/` and `NAMESPACE`
5. Run `Rscript -e 'devtools::check_man()'` or inspect for documentation warnings
6. Report what was added or changed

If `$ARGUMENTS` is provided, focus only on the specified file or function name.

## Style
- Use present tense ("Computes...", "Extracts...", not "This function computes...")
- Parameter descriptions: start with the expected type/class, then explain meaning
- Link related functions with `@seealso`
- Use `@family` tags to group related functions
