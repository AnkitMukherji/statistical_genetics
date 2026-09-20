# Lesson 5: distribute prior support unequally across variants.
# Source: https://chr1swallace.github.io/coloc/articles/a07_variant-specific-priors.html
# Run: Rscript coloc/scripts/05_variant_specific_priors.R
source("coloc/R/helpers.R")
out <- start_lesson("05_variant_specific_priors")
d <- example_data()
D1 <- d$D1
D2 <- d$D2

# 1. Reproduce the example giving the first 100 of 500 SNPs double the weight.
# Weight vectors follow the SNP ORDER IN THE CORRESPONDING INPUT DATASET.
# Supply independent annotation information; choosing weights from these same
# association results would reuse the evidence and bias the comparison.
weights <- c(rep(2, 100), rep(1, 400))
stopifnot(length(weights) == length(D1$snp))
uniform <- coloc::coloc.abf(D1, D2)
weighted <- coloc::coloc.abf(D1, D2, prior_weights1 = weights)
equal <- coloc::coloc.abf(D1, D2, prior_weights1 = rep(1, length(D1$snp)))
stopifnot(isTRUE(all.equal(uniform$summary, equal$summary, tolerance = 1e-10)))

# 2. Use the package's empirical eQTL-to-TSS distance density as annotation.
env <- new.env()
utils::data("eqtlgen_density_data", package = "coloc", envir = env)
density <- env$eqtlgen_density_data
save_plot(out, "eqtl_tss_density", {
  plot(density$rel_dist, density$dens_value, type = "l",
       xlab = "Distance from TSS (bp)", ylab = "Density", main = "eQTLGen distance annotation")
})

# These coordinates are RANDOMLY ASSIGNED to simulated SNPs, as in the vignette.
# This is not a real PTPN22 locus analysis. The fixed seed makes our run repeatable.
positions <- sample(113371759:114371759, length(D1$snp), replace = FALSE)
D1$position <- D2$position <- positions
tss <- 113871759
# Match each relative distance to the nearest sampled density location, as in
# the tutorial. The empirical grid stops slightly short of the simulated window;
# report endpoint matching explicitly instead of silently extrapolating.
tss_weights <- function(position, tss, density) {
  distance <- position - tss
  stopifnot(all(is.finite(distance)))
  outside <- distance < min(density$rel_dist) | distance > max(density$rel_dist)
  if (any(outside)) {
    cat(sum(outside), "distances lie outside the density grid; using the nearest endpoint,",
        "as in the tutorial. Reconsider this policy for real data.\n")
  }
  index <- vapply(distance, function(x) which.min(abs(density$rel_dist - x)), integer(1))
  w <- density$dens_value[index]
  stopifnot(all(is.finite(w)), all(w >= 0), sum(w) > 0)
  w / sum(w)
}
w1 <- tss_weights(D1$position, tss, density)
tss_fit <- coloc::coloc.abf(D1, D2, prior_weights1 = w1)
comparison <- rbind(posterior_row(uniform, "uniform"),
                    posterior_row(weighted, "first_100_double_weight"),
                    posterior_row(tss_fit, "TSS_distance_trait1"))
check_posteriors(comparison)
print(comparison)
save_table(comparison, out, "prior_comparison")
save_table(data.frame(snp = D1$snp, simulated_position = positions,
                      distance_to_TSS = positions - tss, block_weight = weights,
                      normalized_TSS_weight = w1), out, "variant_weights")
save_table(tss_fit$results, out, "TSS_weighted_snps")
save_plot(out, "TSS_weights", {
  plot(positions - tss, w1, pch = 16, cex = 0.65,
       xlab = "Simulated distance from TSS (bp)", ylab = "Normalized prior weight",
       main = "Trait 1 annotation weights")
})
save_objects(list(uniform = uniform, block_weighted = weighted, TSS_weighted = tss_fit), out)
