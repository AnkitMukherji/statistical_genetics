# Colocalisation worked examples

Restore the root renv project. Start R from the repository root and run an individual lesson with Rscript coloc/scripts/FILE.R, or use the parent runner. See the coloc overview for sources and output locations.

## Script and notebook guide

| File | Purpose |
| --- | --- |
| [`01_data_structures.R`](01_data_structures.R) | Builds and validates example coloc datasets and checks variant and allele alignment. |
| [`02_single_causal_variant.R`](02_single_causal_variant.R) | Works through single-trait fine-mapping and two-trait colocalisation under a single causal variant assumption. |
| [`03_prior_sensitivity.R`](03_prior_sensitivity.R) | Checks how posterior support changes when the prior for a shared signal changes. |
| [`04_multiple_signals_susie.R`](04_multiple_signals_susie.R) | Fits multi-signal SuSiE models and compares signal pairs using appropriate LD. |
| [`05_variant_specific_priors.R`](05_variant_specific_priors.R) | Demonstrates how per-variant prior weights change posterior support. |
| [`06_deprecated_conditioning.R`](06_deprecated_conditioning.R) | Reproduces the historical conditioning approach and illustrates sensitivity to secondary-signal selection. |

[Repository overview](../../README.md)
