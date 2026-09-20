# Lesson 1: make valid inputs before interpreting a coloc result.
# Source: https://chr1swallace.github.io/coloc/articles/a02_data.html
# Run: Rscript coloc/scripts/01_data_structures.R
source("coloc/R/helpers.R")
out <- start_lesson("01_data_structures")
d <- example_data()
D1 <- d$D1
D2 <- d$D2

# 1. A dataset is a LIST: variant-level vectors plus trait-level scalar values.
# varbeta is SE^2, not SE. sdY is the phenotype SD, not the SE of beta.
minimum <- D1[c("beta", "varbeta", "snp", "position", "type", "sdY")]
str(minimum)
coloc::check_dataset(minimum, warn.minp = 1e-10)
save_table(data.frame(snp = D1$snp, position = D1$position, beta = D1$beta,
                      SE = sqrt(D1$varbeta), varbeta = D1$varbeta, MAF = D1$MAF),
           out, "example_summary_statistics")
save_plot(out, "association", coloc::plot_dataset(minimum))

# This is only a FORMAT example. Relabeling simulated quantitative effects does
# not turn them into real case-control data. Real binary-trait beta = log(OR).
minimum_cc <- D1[c("beta", "varbeta", "snp", "position")]
minimum_cc$type <- "cc"
coloc::check_dataset(minimum_cc)

# 2. Alternative inputs: if sdY is unknown, provide N and MAF for estimation.
# Only use sdY = 1 if the phenotype really was standardized to variance one.
no_sdY <- D1[c("beta", "varbeta", "snp", "position", "type", "N", "MAF")]
coloc::check_dataset(no_sdY)

# Without effect estimates, a two-sided p-value, MAF, N, and type are needed.
# Prefer beta/varbeta when available, since their precision reflects imputation.
pvalue_data <- D1[c("MAF", "snp", "position", "type", "sdY", "N")]
pvalue_data$pvalues <- 2 * pnorm(-abs(D1$beta / sqrt(D1$varbeta)))
coloc::check_dataset(pvalue_data)
pvalue_cc <- pvalue_data
pvalue_cc$type <- "cc"
pvalue_cc$sdY <- NULL
pvalue_cc$s <- 0.5 # fraction of the study sample that are cases, not prevalence
coloc::check_dataset(pvalue_cc)

# 3. Only shared SNPs contribute. Truncate here ONLY to teach overlap diagnostics;
# real analyses should retain dense, quality-controlled regional coverage.
short1 <- lapply(D1[c("beta", "varbeta", "snp", "position")], `[`, 200:500)
short2 <- lapply(D1[c("beta", "varbeta", "snp", "position")], `[`, 1:400)
short1$type <- short2$type <- "cc"
coloc::check_dataset(short1)
coloc::check_dataset(short2)
stopifnot(length(intersect(short1$snp, short2$snp)) == 201L)
save_plot(out, "shared_coverage", coloc::plot_datasets(minimum, minimum_cc))
save_plot(out, "incomplete_overlap", coloc::plot_datasets(short1, short2))

# 4. SuSiE needs signed correlations (r, not r^2), with SNP IDs and order matching
# beta. Flipping an effect allele requires flipping beta AND its LD row/column.
check_ld(D1)
signs <- sample(c(-1, 1), length(D1$beta), replace = TRUE)
good <- D1
good$beta <- good$beta * signs
good$LD <- good$LD * outer(signs, signs)
bad <- good
bad$beta <- bad$beta * sample(c(-1, 1), length(bad$beta), replace = TRUE)
good_score <- coloc::check_alignment(good, do_plot = FALSE)
bad_score <- coloc::check_alignment(bad, do_plot = FALSE)
# Redraw the same diagnostic with shorter labels: the package's two-line
# subtitle overlaps its x-axis label on this graphics device.
draw_alignment <- function(dataset, label, score) {
  z <- dataset$beta / sqrt(dataset$varbeta)
  ratio <- (outer(z, z) / dataset$LD)[abs(dataset$LD) > 0.2]
  hist(ratio, main = label, xlab = "Z-score product / LD", col = "grey85")
  abline(v = 0, col = "red")
  legend("topright", legend = sprintf("Positive ratios: %.1f%%", 100 * score), bty = "n")
}
save_plot(out, "allele_alignment", {
  par(mfrow = c(1, 2), mar = c(5, 4, 4, 1), oma = c(2, 0, 0, 0), cex = 0.85)
  draw_alignment(good, "Consistent alleles and LD", good_score)
  draw_alignment(bad, "Effects flipped without updating LD", bad_score)
  mtext("Consistent coding tends to give more positive ratios; symmetry can flag poor alignment.",
        side = 1, outer = TRUE, line = 0.3)
}, width = 1800)
save_table(data.frame(example = c("consistent", "deliberately_misaligned"),
                      alignment_score = c(good_score, bad_score)), out, "alignment")

# A malformed input must stop, rather than produce a plausible-looking result.
broken <- minimum
broken$varbeta <- broken$varbeta[-1]
err <- tryCatch(coloc::check_dataset(broken), error = identity)
stopifnot(inherits(err, "error"))
cat("Expected validation error:", conditionMessage(err), "\n")
save_objects(list(minimum = minimum, no_sdY = no_sdY, pvalue_data = pvalue_data), out)
cat("All input checks passed; overlap example has 201 shared SNPs.\n")
