# Lesson 4: decompose multiple signals with SuSiE before comparing traits.
# Source: https://chr1swallace.github.io/coloc/articles/a06_SuSiE.html
# Run: Rscript coloc/scripts/04_multiple_signals_susie.R
source("coloc/R/helpers.R")
out <- start_lesson("04_multiple_signals_susie")
d <- example_data()
D3 <- d$D3
D4 <- d$D4
check_ld(D3)
check_ld(D4)

# A high regional H4 from ABF can hide an additional signal in one trait.
baseline <- coloc::coloc.abf(D3, D4)
save_table(posterior_row(baseline, "single_variant_baseline"), out, "abf_baseline")
save_plot(out, "regional_associations", {
  par(mfrow = c(2, 1))
  coloc::plot_dataset(D3, main = "D3: multiple signals")
  coloc::plot_dataset(D4, main = "D4")
})
save_plot(out, "baseline_sensitivity", coloc::sensitivity(baseline, "H4 > 0.9"))

# Fit once per trait and reuse. An LD matrix must reflect the ancestry and effect
# allele coding of THAT study; separate studies may need different matrices.
S3 <- coloc::runsusie(D3)
S4 <- coloc::runsusie(D4)
stopifnot(isTRUE(S3$converged), isTRUE(S4$converged))
print(summary(S3))
print(summary(S4))
save_table(summary(S3)$cs, out, "D3_credible_sets")
save_table(summary(S4)$cs, out, "D4_credible_sets")

# Each row compares ONE signal from D3 with ONE signal from D4. It is not a new
# H4 for the whole locus, and row numbers do not imply stable biological labels.
fit <- coloc::coloc.susie(S3, S4)
check_posteriors(fit$summary)
print(fit$summary)
save_table(fit$summary, out, "signal_pair_posteriors")
for (i in seq_len(nrow(fit$summary))) {
  save_plot(out, paste0("pair_", i, "_sensitivity"), {
    coloc::sensitivity(fit, "H4 > 0.9", row = i, dataset1 = D3, dataset2 = D4)
  })
}

# Optional example in the upstream text: weaken the data and reduce coverage.
# Reproduce it, but 10% coverage is an illustration, not a discovery standard.
D5 <- D3
D5$varbeta <- D5$varbeta * 4
D5$N <- D5$N / 2
S5_default <- coloc::runsusie(D5)
S5_low <- coloc::runsusie(D5, coverage = 0.1)
stopifnot(isTRUE(S5_default$converged), isTRUE(S5_low$converged))
print(summary(S5_default))
print(summary(S5_low))
save_table(data.frame(coverage = c(0.95, 0.1),
                      credible_sets = c(length(S5_default$sets$cs), length(S5_low$sets$cs))),
           out, "weak_signal_coverage")
save_plot(out, "weak_signal_associations", {
  par(mfrow = c(2, 1))
  coloc::plot_dataset(D3, main = "Original D3")
  coloc::plot_dataset(D5, main = "Artificially weakened D5")
})
weak_fit <- NULL
if (length(S5_low$sets$cs) > 0 && length(S4$sets$cs) > 0) {
  weak_fit <- coloc::coloc.susie(S5_low, S4)
  check_posteriors(weak_fit$summary)
  save_table(weak_fit$summary, out, "weak_signal_pair_posteriors")
} else {
  cat("No credible set available for the weak comparison; this is not evidence for H3.\n")
}
save_objects(list(S3 = S3, S4 = S4, coloc = fit, S5_default = S5_default,
                  S5_low = S5_low, weak_coloc = weak_fit), out)
