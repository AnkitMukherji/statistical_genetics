# Check scientific outcomes from a completed run, without refitting models.
# Run: Rscript coloc/validate_results.R
source("coloc/R/helpers.R")
read_result <- function(lesson, filename) {
  utils::read.csv(file.path("coloc/outputs", lesson, paste0(filename, ".csv")), check.names = FALSE)
}
abf <- read_result("02_single_causal_variant", "hypothesis_posteriors")
check_posteriors(abf)
# Reference D1/D2 H4 is approximately 0.999166, not an exact platform-specific bit pattern.
stopifnot(abf$nsnps == 500, abs(abf$H4 - 0.9991661) < 1e-5)
fm <- read_result("02_single_causal_variant", "finemap_D1")
stopifnot(abs(sum(fm$SNP.PP) - 1) < 1e-8, "null" %in% fm$snp)
cs <- read_result("02_single_causal_variant", "shared_variant_credible_set_95")
stopifnot(sum(cs$SNP.PP.H4) >= 0.95,
          sum(head(cs$SNP.PP.H4, -1)) < 0.95)
scan <- read_result("03_prior_sensitivity", "p12_grid")
check_posteriors(scan)
for (label in unique(scan$analysis)) {
  x <- subset(scan, analysis == label)
  stopifnot(all(diff(x$H4[order(x$p12)]) >= -1e-10))
}
su <- read_result("04_multiple_signals_susie", "signal_pair_posteriors")
check_posteriors(su)
stopifnot(nrow(su) == 2L, any(su$PP.H4.abf > 0.9), any(su$PP.H3.abf > 0.9))
weights <- read_result("05_variant_specific_priors", "variant_weights")
stopifnot(nrow(weights) == 500L, all(weights$normalized_TSS_weight >= 0),
          abs(sum(weights$normalized_TSS_weight) - 1) < 1e-8)
check_posteriors(read_result("05_variant_specific_priors", "prior_comparison"))
check_posteriors(read_result("06_deprecated_conditioning", "conditional_pair_posteriors"))
writeLines(c("PASS: ABF reproduces the reference D1/D2 posterior.",
             "PASS: fine-mapping and H0-H4 probabilities are normalized.",
             "PASS: the H4 credible set reaches 95% with the crossing SNP included.",
             "PASS: H4 increases with p12 in both sensitivity examples.",
             "PASS: SuSiE identifies a shared pair and a distinct pair of signals.",
             "PASS: variant weights and deprecated conditioning results are valid."),
           "coloc/outputs/validation.txt")
# Load the actual analytical dependencies before recording session information.
invisible(lapply(c("coloc", "susieR"), requireNamespace, quietly = TRUE))
writeLines(sub("[[:blank:]]+$", "", capture.output(sessionInfo())),
           "coloc/outputs/session_info.txt")
cat("Scientific result checks passed.\n")
