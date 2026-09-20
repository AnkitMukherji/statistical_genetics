# Sources and adaptation notes

Consulted 20 September 2026. The starting point is the
[coloc website](https://chr1swallace.github.io/coloc/index.html), displaying version
6.0.1. The scripts adapt the examples in its linked articles by Chris Wallace and
Jeffrey Pullin. Article links are included in every lesson and in the reading
guide. The package and bundled datasets are installed through `renv`, rather than
vendored into this repository.

## Source provenance

- [Source repository](https://github.com/chr1swallace/coloc)
- [Pinned installed source, version 6.0.3](https://github.com/chr1swallace/coloc/tree/17e0b168face59bff1df1d59d927bfa74072b578)
- [Original sensitivity Rmd, including hidden data construction](https://github.com/chr1swallace/coloc/blob/17e0b168face59bff1df1d59d927bfa74072b578/vignettes/a04_sensitivity.Rmd)
- [Original SuSiE Rmd, including weak-signal examples](https://github.com/chr1swallace/coloc/blob/17e0b168face59bff1df1d59d927bfa74072b578/vignettes/a06_SuSiE.Rmd)

The original package declares GPL. The adapted teaching code and documentation
in this `coloc/` section are distributed under GPL version 3 or later; see
[LICENSE](LICENSE). This does not change the license of unrelated repository
content or third-party packages. Credit for the underlying methods, example data,
and tutorial designs remains with the coloc authors.

## Deliberate adaptations

1. Explicit `D1`–`D4` objects replace `attach()`, preventing hidden name lookup.
2. Each lesson can run alone from the repository root. A fixed seed, saved tables,
   headless PNG plots, and result checks make execution inspectable.
3. The article numbering is reordered to teach SuSiE before legacy conditioning.
4. Input-format case-control examples are explicitly labeled as such. No real
   binary-trait data are inferred from the simulated quantitative effects.
5. Sensitivity's hidden variance/N changes are reproduced literally. They are
   illustrative perturbations, not a consistent sample-size simulation.
6. A numerical prior grid and H4 credible-set checks extend the tutorials. These
   are labeled as teaching additions, and the grid is not a prescribed prior range.
7. The weak SuSiE examples are executed, including the deliberately low coverage.
   Their inference is described separately from the main 95% analysis.
8. TSS positions use seed 20260920. The upstream draw has no local seed, so its
   printed weighted posterior is not an exact numeric target for this adaptation.
   Nearest endpoint matching is reported if simulated positions slightly exceed
   the empirical density grid; this follows the upstream nearest-bin calculation.
9. Full model objects are reproducible local files. Small tables, plots, and logs
   are versioned to allow reading the course on GitHub without running R.
10. The alignment histogram uses the same Z-product/LD calculation as
    `check_alignment()`, with shorter labels to prevent overlapping annotations.

## Papers cited by the tutorials

- Giambartolomei et al. (2014), [Bayesian test for colocalisation between pairs of
  genetic association studies using summary statistics](https://doi.org/10.1371/journal.pgen.1004383).
- Wallace (2020), [Eliciting priors and relaxing the single causal variant
  assumption in colocalisation analyses](https://doi.org/10.1371/journal.pgen.1008720).
- Wang et al. (2020), [A simple new approach to variable selection in regression,
  with application to genetic fine mapping](https://doi.org/10.1111/rssb.12388).
- Wallace (2021), [A more accurate method for colocalisation analysis allowing
  for multiple causal variants](https://doi.org/10.1371/journal.pgen.1009440).
- Pullin and Wallace (2025), [Variant-specific priors clarify colocalisation
  analysis](https://doi.org/10.1371/journal.pgen.1011697).

Use `citation("coloc")` and `citation("susieR")` for citations associated with the
installed package versions when reporting your own analyses.
