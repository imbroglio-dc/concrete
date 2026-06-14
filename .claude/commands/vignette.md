Write or update a package vignette: $ARGUMENTS

## Steps
1. If updating an existing vignette, read it first and identify what needs to change
2. If creating a new vignette:
   - Create `vignettes/<topic>.Rmd` with proper YAML header:
     ```yaml
     ---
     title: "<Title>"
     output: rmarkdown::html_vignette
     vignette: >
       %\VignetteIndexEntry{<Title>}
       %\VignetteEngine{knitr::rmarkdown}
       %\VignetteEncoding{UTF-8}
     ---
     ```
   - Ensure `knitr` and `rmarkdown` are in `Suggests:` in DESCRIPTION
3. Write content that:
   - Starts with the problem the user is trying to solve (not the API)
   - Uses a realistic worked example (simulated data or a built-in dataset)
   - Shows complete, runnable code chunks
   - Explains the statistical interpretation of results
   - Follows a logical progression: setup → analysis → results → interpretation
4. Build the vignette: `Rscript -e 'devtools::build_vignettes()'`
5. Fix any errors in code chunks
6. Report what was created or changed

## Style
- Write for a biostatistician audience who knows survival analysis but not necessarily TMLE
- Use `concrete`'s actual API — don't invent functions that don't exist
- Keep code chunks under 20 lines; explain between chunks
- Use `eval=FALSE` sparingly — runnable code builds trust
