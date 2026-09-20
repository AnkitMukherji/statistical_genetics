# Implementation plan

1. Inspect the existing repository, R environment, and all articles linked from
   <https://chr1swallace.github.io/coloc/index.html>.
2. Write six independent, annotated R lessons: data structures, ABF enumeration,
   sensitivity, SuSiE, variant-specific priors, and deprecated conditioning.
   Cover the introductory article in the reading guide.
3. Use the package's bundled examples; save figures and tables, and provide a
   runner that records console output and package versions.
4. Populate the existing project `renv` library and snapshot the dependencies.
   Version 6 is needed for variant-specific priors; pin the source revision if
   the configured CRAN mirror does not provide it.
5. Execute every lesson, check posterior normalization and the examples' main
   scientific conclusions, inspect plots, and document observed results and
   differences from the website.
6. Commit the section and dependency lockfile and push to the existing GitHub
   repository. Keep generated bulky R objects and installed libraries local.
