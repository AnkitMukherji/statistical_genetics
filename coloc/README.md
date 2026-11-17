# Learning colocalisation with coloc

Six runnable lessons adapted from the [official coloc tutorial collection](https://chr1swallace.github.io/coloc/index.html).
Start with the reading order below, open a script, run its sections, and compare
your output with the committed tables and figures. The examples use **simulated
data bundled with coloc**, so no GWAS downloads or private data are needed.

## Run the lessons

Start a terminal or R session in the **statistical_genetics repository root**, one
directory above this README. Its `.Rprofile` activates the existing `renv` project.

```sh
# On another machine, restore the versions recorded in the root renv.lock first:
Rscript -e 'renv::restore()'

# Run all six lessons, including the historical conditioning example:
Rscript coloc/run_all.R

# Or run one lesson independently:
Rscript coloc/scripts/02_single_causal_variant.R

# Recheck already-generated results without refitting:
Rscript coloc/validate_results.R
```

In RStudio, open the repository as your working directory, restart R so that
`.Rprofile` runs, and use `source("coloc/scripts/02_single_causal_variant.R")` or
execute its lines in order. Each lesson loads fresh data and can run independently.
Do not use `Rscript --vanilla` for these lessons: that bypasses project activation.
The scripts do not install packages or modify the lockfile while running.

The tested environment uses **R 4.6.1, coloc 6.0.3, and susieR 0.14.2**. The website
displays coloc 6.0.1. The configured CRAN mirror supplied only 5.2.3, so this project
pins coloc's upstream commit
[`17e0b168face59bff1df1d59d927bfa74072b578`](https://github.com/chr1swallace/coloc/tree/17e0b168face59bff1df1d59d927bfa74072b578)
in `renv.lock`. Restoring that GitHub source needs network access; running the
lessons after restoration does not. `renv` restores packages, not the R executable.

## Reading order and tutorial coverage

| Lesson | Script | Upstream article | Main output to inspect |
| --- | --- | --- | --- |
| Orientation | This guide | [Introduction](https://chr1swallace.github.io/coloc/articles/a01_intro.html) | Hypothesis table below |
| 1. Data structures | [01_data_structures.R](scripts/01_data_structures.R) | [Data](https://chr1swallace.github.io/coloc/articles/a02_data.html) | `alignment.csv`, overlap/alignment plots |
| 2. Single causal variant | [02_single_causal_variant.R](scripts/02_single_causal_variant.R) | [Enumeration](https://chr1swallace.github.io/coloc/articles/a03_enumeration.html) | `hypothesis_posteriors.csv`, `shared_variant_credible_set_95.csv` |
| 3. Prior sensitivity | [03_prior_sensitivity.R](scripts/03_prior_sensitivity.R) | [Sensitivity](https://chr1swallace.github.io/coloc/articles/a04_sensitivity.html) | `p12_grid.csv`, `H4_vs_p12.png` |
| 4. Multiple signals | [04_multiple_signals_susie.R](scripts/04_multiple_signals_susie.R) | [SuSiE](https://chr1swallace.github.io/coloc/articles/a06_SuSiE.html) | `signal_pair_posteriors.csv`, credible-set tables |
| 5. Annotation priors | [05_variant_specific_priors.R](scripts/05_variant_specific_priors.R) | [Variant-specific priors](https://chr1swallace.github.io/coloc/articles/a07_variant-specific-priors.html) | `variant_weights.csv`, `prior_comparison.csv` |
| 6. Historical method | [06_deprecated_conditioning.R](scripts/06_deprecated_conditioning.R) | [Deprecated conditioning](https://chr1swallace.github.io/coloc/articles/a05_conditioning.html) | `selected_signals.csv`, `conditional_pair_posteriors.csv` |

The introduction is conceptual rather than a separate runnable analysis. Lesson 6
is retained to cover the linked articles; use lesson 4 for new multiple-signal work.
The deprecated article's executed examples use conditioning, so this section does
not invent an additional masking workflow.

## The question coloc answers

Nearby associated SNPs can reflect one causal variant because of linkage
disequilibrium (LD). Two traits with nearby association peaks therefore need not
share a causal variant. Coloc compares explanations of the regional data:

| Hypothesis | Meaning under the single-causal-variant model |
| --- | --- |
| H0 | Neither trait has a causal variant in this region. |
| H1 | Only trait 1 has a causal variant. |
| H2 | Only trait 2 has a causal variant. |
| H3 | Both traits have causal variants, but they are different. |
| H4 | Both traits share one causal variant. |

`PP.H4.abf` is posterior support for H4 under the fitted model and priors. Inspect
**all five probabilities**: weak support for H4 can reflect weak association
evidence, rather than positive evidence for two different variants. Strong H3
supports distinct signals. Colocalisation alone does not establish that one trait
causes the other, identify a mechanism, or prove which gene mediates an effect.

### 1. Understand the inputs

Read the `str(minimum)` output alongside the script. The vectors describe SNPs;
`type`, `N`, and `sdY` describe the study or phenotype. A list can contain both,
whereas a data frame would repeat scalar fields down its rows.

| Field | Interpretation |
| --- | --- |
| `snp` | Unique variant identifiers used to match studies. |
| `position` | Coordinates for plotting; the toy data initially use indices 1–500. |
| `beta` | Estimated effect; use log odds ratios for a binary phenotype. |
| `varbeta` | Sampling variance of the effect estimate, equal to `SE^2`. |
| `type` | `"quant"` or `"cc"`. |
| `sdY` | Phenotype standard deviation for a quantitative trait. |
| `N`, `MAF` | Sample size and allele frequency information used by relevant input routes. |
| `pvalues`, `s` | Two-sided association p-values; `s` is the sample case fraction for a binary p-value input. |
| `LD` | Signed SNP correlation matrix needed for SuSiE/conditioning. |

The case-control examples demonstrate the interface only; their effects are
copied from simulated quantitative data. They are not binary phenotype analyses.
The deliberately broken input should print a caught validation error and then
continue. The truncated datasets should have 201 shared SNPs.

**Try:** inspect the good and deliberately misaligned LD plots. Why must an
effect-allele flip change both beta and a row/column of LD? Alignment diagnostics
can expose inconsistencies, but they cannot certify harmonization from SNP names.

### 2. Separate regional evidence from candidate-variant ranking

`finemap.abf()` analyzes one trait and includes a `null` row. Its `SNP.PP`
probabilities include that no-association model. `coloc.abf()` then compares the two
traits. Its `SNP.PP.H4` is a different quantity: the probability distribution over
shared SNPs **conditional on H4**. A variant can rank first even in a region with
little evidence for H4.

The script sorts `SNP.PP.H4` and keeps SNPs through the first cumulative probability
of at least 0.95. This creates a conditional credible set, not a list of proven
causal variants. The numerical checks include the boundary-crossing SNP.

**Try:** compare the region's H4 with its top SNP's `SNP.PP.H4`. Explain why both
quantities are needed when reporting a candidate shared variant.

### 3. Treat priors as part of the analysis

`p1`, `p2`, and `p12` encode per-SNP prior association/sharing probabilities. They
are not association p-values or the final H1/H2/H4 probabilities. Regional priors
also depend on the number of SNPs and possible configurations.

This lesson holds `p1` and `p2` fixed and varies `p12`. The `H4 > 0.5` and compound
rules are teaching examples, not universal thresholds. Increasing support for
sharing a priori can change the conclusion, especially with less precise data.
Choose a defensible range before looking for a preferred result.

**Try:** find the first passing `p12` in `p12_grid.csv` for each dataset. The source
article's weaker-data construction is hidden in its Rmd; it is made explicit here.
Its prose and plotted/numerical example are not perfectly consistent, so use the
computed grid rather than copying a claimed threshold.

### 4. Read SuSiE results one signal pair at a time

`runsusie()` fits each trait using summary statistics and LD. The saved objects
can be reused in `coloc.susie()` without refitting. In this example D3 has two
credible sets and D4 has one, giving two comparisons. Read `hit1`, `hit2`, and
the posterior columns together. Do not assume that a particular row number
always corresponds to the same signal across versions.

The weakened example produces no credible set at 95% coverage and one at 10% in
the tested environment. Lowering coverage changes the inferential target. It
does not create strong evidence or justify reporting a 10% set as a 95% set.
An empty set does not demonstrate that two traits have distinct causal variants.

**Try:** compare the baseline ABF result with the two SuSiE rows. Which additional
signal was obscured by the baseline's single-variant assumption?

### 5. Use external annotation information carefully

Weights redistribute prior support across variants. The lesson first doubles
weights for SNPs 1–100, then uses the bundled eQTLGen distance-to-transcription-
start-site (TSS) density. Weights must follow each dataset's SNP order; when
reordering or subsetting SNPs, update weights too.

The TSS example assigns random coordinates near a PTPN22 coordinate to simulated
SNPs. It does not analyze real PTPN22 associations. A fixed seed makes this
adaptation reproducible; the website's unseeded draw need not give the same H4.
Nearest-density-bin matching follows the source; endpoint matches are reported.

**Try:** compare uniform and weighted posteriors. Could a prior made from the
same association statistics double-count evidence? Prefer independent annotation
and document why it is relevant to the trait.

### 6. Understand why conditioning is historical

The old workflow selects signals and conditions on other selected variants.
Changing `pthr` can omit a secondary signal or add weak ones. Errors in selection
can affect the later colocalisation comparisons. The upstream documentation
recommends SuSiE instead; this lesson explains older analyses rather than setting
the default for new work.

## Moving from toy examples to your own region

Prepare one region for one trait in each study. Use compatible genome builds,
unambiguous variant IDs, harmonized alleles, valid effect variances, and dense
coverage after appropriate study QC. Do not LD-clump or retain only significant
SNPs for coloc. Inspect which SNPs are shared and whether missing variants leave
holes around an association peak. A whole chromosome or multiple genes' eQTL
records pooled together are not a single coloc dataset.

For SuSiE, match each study to suitable LD, including ancestry, variant order,
and effect-allele coding. `r^2` loses the signs needed for the model. The helpers
check matrix dimensions, names, symmetry, diagonal, and finite values; they do
not resolve ambiguous alleles, validate ancestry, or replace scientific QC.

## Files and reproducibility

- `scripts/`: annotated independent lessons; start here to learn the R workflow.
- `R/helpers.R`: shared output and validation utilities.
- `run_all.R`: executes lessons in separate processes and records `run.log` files.
- `validate_results.R`: checks scientific outcomes and posterior normalization.
- `outputs/<lesson>/`: committed CSV tables, PNG figures, and execution logs.
- `outputs/validation.txt` and `outputs/session_info.txt`: verification and versions.
- `outputs/<lesson>/objects.rds`: full objects created locally, ignored by Git.
  For example, `readRDS("coloc/outputs/04_multiple_signals_susie/objects.rds")$S3`
  retrieves the fitted D3 SuSiE object. Rerunning rebuilds these files.

Results are overwritten when their lesson is rerun. See [RESULTS.md](RESULTS.md)
for the tested numerical outcomes and [SOURCES.md](SOURCES.md) for attribution,
upstream references, and deliberate adaptations. The root lockfile contains the
dependencies; no second `renv` project is created inside this section.

## Run and inspect the code

- [Lesson script guide](scripts/README.md) explains each worked example.
- [Shared utilities](R/README.md) documents the functions used across lessons.
