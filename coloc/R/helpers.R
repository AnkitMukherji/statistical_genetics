# Shared file handling and checks. The statistical steps remain in each lesson.
# Run from the statistical_genetics repository root so .Rprofile activates renv.

start_lesson <- function(name) {
  if (!file.exists("renv.lock") || !dir.exists("coloc/scripts")) {
    stop("Start R in the statistical_genetics repository root; see coloc/README.md.")
  }
  if (!requireNamespace("coloc", quietly = TRUE) ||
      utils::packageVersion("coloc") < "6.0.1") {
    stop("This course needs coloc >= 6.0.1. Run renv::restore() from the repository root.")
  }
  set.seed(20260920)
  output <- file.path("coloc", "outputs", name)
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  cat("\nLesson:", name, "\ncoloc:", as.character(utils::packageVersion("coloc")), "\n")
  output
}

example_data <- function() {
  env <- new.env()
  utils::data("coloc_test_data", package = "coloc", envir = env)
  env$coloc_test_data
}

save_table <- function(x, output, name) {
  utils::write.csv(as.data.frame(x), file.path(output, paste0(name, ".csv")),
                   row.names = FALSE, na = "NA")
}

save_plot <- function(output, name, code, width = 1600, height = 1100) {
  grDevices::png(file.path(output, paste0(name, ".png")),
                width = width, height = height, res = 150)
  on.exit(grDevices::dev.off())
  eval(substitute(code), envir = parent.frame())
  invisible(NULL)
}

save_objects <- function(x, output) {
  # RDS objects are useful for exploration but excluded from Git; rerun to rebuild.
  saveRDS(x, file.path(output, "objects.rds"))
}

posterior_row <- function(result, label) {
  values <- result$summary
  data.frame(analysis = label, nsnps = unname(values["nsnps"]),
             H0 = unname(values["PP.H0.abf"]), H1 = unname(values["PP.H1.abf"]),
             H2 = unname(values["PP.H2.abf"]), H3 = unname(values["PP.H3.abf"]),
             H4 = unname(values["PP.H4.abf"]))
}

check_posteriors <- function(x) {
  # Accept either coloc's PP.H*.abf columns or the exported short column names.
  x <- as.data.frame(x)
  columns <- paste0("PP.H", 0:4, ".abf")
  if (!all(columns %in% names(x))) columns <- paste0("H", 0:4)
  p <- as.matrix(x[, columns, drop = FALSE])
  stopifnot(nrow(p) > 0, all(is.finite(p)), all(p >= 0), all(p <= 1),
            all(abs(rowSums(p) - 1) < 1e-8))
  invisible(TRUE)
}

check_ld <- function(dataset) {
  coloc::check_dataset(dataset, req = c("beta", "varbeta", "snp", "type", "N", "LD"))
  r <- dataset$LD
  stopifnot(!anyDuplicated(dataset$snp), all(is.finite(r)),
            identical(rownames(r), dataset$snp), identical(colnames(r), dataset$snp),
            isSymmetric(r, tol = 1e-8), all(abs(diag(r) - 1) < 1e-8),
            all(abs(r) <= 1 + 1e-8))
  # These checks cannot establish that the alleles were harmonized correctly.
  invisible(TRUE)
}
