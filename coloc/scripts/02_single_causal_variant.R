# Lesson 2: fine mapping and colocalisation under at most one causal SNP per trait.
# Source: https://chr1swallace.github.io/coloc/articles/a03_enumeration.html
# Run: Rscript coloc/scripts/02_single_causal_variant.R
source("coloc/R/helpers.R")
out <- start_lesson("02_single_causal_variant")
d <- example_data()
D1 <- d$D1
D2 <- d$D2
coloc::check_dataset(D1)
coloc::check_dataset(D2)

# 1. Single-trait fine mapping includes a null row (no causal variant).
fm <- coloc::finemap.abf(dataset = D1)
print(fm[21:30, ])
print(tail(fm, 3))
stopifnot(abs(sum(fm$SNP.PP) - 1) < 1e-8, "null" %in% fm$snp)
save_table(fm, out, "finemap_D1")

# 2. Compare two traits. p1/p2/p12 are SNP-level association priors, not p-values
# and not the final region-level probabilities H0-H4.
fit <- coloc::coloc.abf(D1, D2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
print(fit)
post <- posterior_row(fit, "D1_vs_D2")
check_posteriors(post)
save_table(post, out, "hypothesis_posteriors")
save_table(fit$results, out, "all_snps")

# 3. SNP.PP.H4 ranks candidate shared variants CONDITIONAL on H4 being true.
# A high SNP.PP.H4 alone cannot establish that the region colocalises.
ranked <- fit$results[order(fit$results$SNP.PP.H4, decreasing = TRUE), ]
ranked$cumulative_PP_H4 <- cumsum(ranked$SNP.PP.H4)
k <- which(ranked$cumulative_PP_H4 >= 0.95)[1]
credible_set <- ranked[seq_len(k), ]
stopifnot(abs(sum(ranked$SNP.PP.H4) - 1) < 1e-8,
          tail(credible_set$cumulative_PP_H4, 1) >= 0.95)
print(subset(fit$results, SNP.PP.H4 > 0.01))
print(credible_set[, c("snp", "SNP.PP.H4", "cumulative_PP_H4")])
save_table(credible_set, out, "shared_variant_credible_set_95")
save_plot(out, "regional_associations", {
  par(mfrow = c(2, 1))
  coloc::plot_dataset(D1, main = "D1")
  coloc::plot_dataset(D2, main = "D2")
})
save_plot(out, "hypothesis_posteriors", {
  barplot(unlist(post[paste0("H", 0:4)]), ylim = c(0, 1),
          ylab = "Posterior probability", main = "D1 versus D2: single causal variant")
})
save_objects(list(finemap = fm, coloc = fit, credible_set = credible_set), out)
cat("P(H4) =", post$H4, "; conditional 95% set size =", k, "\n")
