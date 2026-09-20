# Lesson 3: does the conclusion depend on the prior probability of sharing?
# Source: https://chr1swallace.github.io/coloc/articles/a04_sensitivity.html
# Run: Rscript coloc/scripts/03_prior_sensitivity.R
source("coloc/R/helpers.R")
out <- start_lesson("03_prior_sensitivity")
d <- example_data()
D1 <- d$D1
D2 <- d$D2
strong <- coloc::coloc.abf(D1, D2, p12 = 1e-6)
save_plot(out, "original_sensitivity", coloc::sensitivity(strong, rule = "H4 > 0.5"))

# Reproduce the source Rmd's HIDDEN construction of its weaker examples exactly.
# These artificial changes illustrate uncertainty, not a coherent sample-size
# simulation: multiplying variance by 2 does not correspond to quartering N.
D1a <- D1
D1a$varbeta <- D1$varbeta * sqrt(4)
D1a$N <- D1$N / 4
D2a <- D2
D2a$varbeta <- D2$varbeta * sqrt(2)
D2a$N <- D2$N / 2
weak <- coloc::coloc.abf(D1a, D2a, p12 = 1e-6)
save_plot(out, "weaker_sensitivity", coloc::sensitivity(weak, rule = "H4 > 0.5"))
save_plot(out, "weaker_compound_rule", {
  coloc::sensitivity(weak, rule = "H4 > 3*H3 & H0 < 0.1")
})

# Extension: export numerical sensitivity so the lesson is readable without R.
# Hold p1 and p2 fixed. This grid is illustrative, not a recommended prior range.
grid <- 10^seq(-8, -4, length.out = 41)
scan_priors <- function(a, b, label) {
  rows <- lapply(grid, function(p12) {
    # Capture routine coloc printing during the sweep; warnings remain visible.
    invisible(capture.output(fit <- coloc::coloc.abf(a, b, p12 = p12)))
    row <- posterior_row(fit, label)
    row$p12 <- p12
    row$passes_H4_0.5 <- row$H4 > 0.5
    row$passes_compound <- row$H4 > 3 * row$H3 & row$H0 < 0.1
    row
  })
  do.call(rbind, rows)
}
scan <- rbind(scan_priors(D1, D2, "original"), scan_priors(D1a, D2a, "weaker"))
check_posteriors(scan)
save_table(scan, out, "p12_grid")
save_table(rbind(posterior_row(strong, "original"), posterior_row(weak, "weaker")),
           out, "posteriors_at_p12_1e-6")
save_plot(out, "H4_vs_p12", {
  a <- subset(scan, analysis == "original")
  b <- subset(scan, analysis == "weaker")
  plot(a$p12, a$H4, log = "x", type = "l", lwd = 2, ylim = c(0, 1),
       xlab = "Per-SNP shared prior p12", ylab = "P(H4 | data)", col = "navy")
  lines(b$p12, b$H4, col = "tomato3", lwd = 2)
  abline(h = 0.5, lty = 2)
  legend("bottomright", c("Original", "Weaker", "Illustrative 0.5 rule"),
         col = c("navy", "tomato3", "black"), lty = c(1, 1, 2), bty = "n")
})
save_objects(list(original = strong, weaker = weak), out)
print(subset(scan, passes_H4_0.5)[c(1, sum(scan$analysis == "original" & scan$passes_H4_0.5) + 1), ])
cat("Use the exported grid to read thresholds; do not assume the prose's threshold is exact.\n")
