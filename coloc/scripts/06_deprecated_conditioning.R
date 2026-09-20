# Historical lesson: reproduce the DEPRECATED conditioning tutorial.
# Source: https://chr1swallace.github.io/coloc/articles/a05_conditioning.html
# Run: Rscript coloc/scripts/06_deprecated_conditioning.R
# Use lesson 4 (SuSiE) for new multiple-signal analyses.
source("coloc/R/helpers.R")
out <- start_lesson("06_deprecated_conditioning")
cat("HISTORICAL EXAMPLE: conditioning is deprecated upstream; prefer SuSiE.\n")
d <- example_data()
D3 <- d$D3
D4 <- d$D4
check_ld(D3)
check_ld(D4)
baseline <- coloc::coloc.abf(D3, D4)
save_table(posterior_row(baseline, "single_variant_baseline"), out, "abf_baseline")
save_plot(out, "baseline_sensitivity", coloc::sensitivity(baseline, "H4 > 0.9"))

# pthr controls which secondary signals get selected. This demonstration shows
# both losing a weaker signal and admitting additional weak signals.
scenarios <- list(D3_default = list(data = D3, pthr = 1e-6),
                  D4_default = list(data = D4, pthr = 1e-6),
                  D3_stringent = list(data = D3, pthr = 1e-20),
                  D4_permissive = list(data = D4, pthr = 0.1))
selected <- lapply(names(scenarios), function(label) {
  x <- scenarios[[label]]
  signals <- coloc::finemap.signals(x$data, method = "cond", pthr = x$pthr)
  print(signals)
  data.frame(scenario = label, pthr = x$pthr, snp = names(signals), z = unname(signals))
})
save_table(do.call(rbind, selected), out, "selected_signals")

# Each comparison conditions out the other selected signals before assessing
# whether the remaining pair shares a causal variant. Selection errors propagate.
fit <- coloc::coloc.signals(D3, D4, method = "cond", p12 = 1e-6, pthr = 1e-6)
print(fit)
check_posteriors(fit$summary)
save_table(fit$summary, out, "conditional_pair_posteriors")
for (i in seq_len(nrow(fit$summary))) {
  save_plot(out, paste0("pair_", i, "_sensitivity"), coloc::sensitivity(fit, "H4 > 0.9", row = i))
}
save_objects(list(baseline = baseline, conditioning = fit), out)
