# Statistical Genetics: Evidence, Assumptions, and Reproducibility

R analyses and annotated lessons focused on interpreting genetic association evidence. The current module explores **Bayesian colocalisation**, from input alignment and model assumptions to prior sensitivity and multiple causal signals.

## Colocalisation module

[Start with the complete tutorial](coloc/README.md), then compare the [recorded results](coloc/RESULTS.md) with the [source and adaptation notes](coloc/SOURCES.md).

| Lesson | Focus |
| --- | --- |
| [Data structures](coloc/scripts/01_data_structures.R) | Input fields, shared variants, and beta/LD alignment |
| [Single causal variant](coloc/scripts/02_single_causal_variant.R) | Regional hypotheses, fine mapping, and conditional credible sets |
| [Prior sensitivity](coloc/scripts/03_prior_sensitivity.R) | How sharing priors change posterior support |
| [Multiple signals](coloc/scripts/04_multiple_signals_susie.R) | SuSiE and signal-pair colocalisation |
| [Variant-specific priors](coloc/scripts/05_variant_specific_priors.R) | Annotation-informed prior weights |
| [Historical conditioning](coloc/scripts/06_deprecated_conditioning.R) | Understanding a deprecated workflow |

## Run and inspect

From the repository root:

```bash
Rscript -e 'renv::restore()'
Rscript coloc/run_all.R
Rscript coloc/validate_results.R
```

The root `.Rprofile` activates the recorded environment. Do not add `--vanilla`, which bypasses that activation. See the module guide for the pinned package source, R-version requirements, individual lesson commands, and output behavior.

The committed run records **R 4.6.1, coloc 6.0.3, and susieR 0.14.2**. The examples use simulated datasets bundled with coloc; they are methodological demonstrations rather than findings about real traits.

## Interpretation principles

Read all regional hypotheses, separate regional sharing support from variant rankings conditional on sharing, and report sensitivity to priors. For multiple signals, interpret each signal pair with appropriately matched LD and allele coding.

Colocalisation supports shared genetic architecture under specified assumptions. It does not by itself prove trait causality, a molecular mechanism, or a mediating gene. The recorded analysis also examines a discrepancy between upstream prose and computed prior-sensitivity results, making interpretation checkable against exported tables.

Related: [population pharmacogenomics](https://github.com/AnkitMukherji/Medical_Relevance_GI) · [phasing and imputation](https://github.com/AnkitMukherji/Phasing_Imputation).
