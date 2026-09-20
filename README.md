# Statistical Genetics

Reproducible R analyses and learning notes, using the project-level `renv` environment.

## Tutorials

- [Colocalisation with coloc](coloc/README.md): six annotated lessons covering data
  structures, single-variant colocalisation, prior sensitivity, SuSiE, variant-specific
  priors, and deprecated conditioning. Includes plots, result tables, and checks.

From this repository's root:

```sh
Rscript -e 'renv::restore()'
Rscript coloc/run_all.R
```
