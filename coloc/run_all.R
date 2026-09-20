#!/usr/bin/env Rscript
# Run from the repository root: Rscript coloc/run_all.R
# Each lesson runs in a separate R process, avoiding accidental shared objects.
if (!file.exists("renv.lock") || !dir.exists("coloc/scripts")) {
  stop("Run from the statistical_genetics repository root: Rscript coloc/run_all.R")
}
dir.create("coloc/outputs", recursive = TRUE, showWarnings = FALSE)
scripts <- sort(list.files("coloc/scripts", pattern = "^[0-9].*\\.R$", full.names = TRUE))
stopifnot(length(scripts) == 6L)
for (script in scripts) {
  name <- tools::file_path_sans_ext(basename(script))
  output <- file.path("coloc/outputs", name)
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  log <- file.path(output, "run.log")
  cat("Running", script, "...\n")
  status <- system2(file.path(R.home("bin"), "Rscript"), shQuote(script), stdout = log, stderr = log)
  # R's console tables contain trailing spaces; normalize only line endings for Git.
  writeLines(sub("[[:blank:]]+$", "", readLines(log, warn = FALSE)), log)
  if (status != 0L) {
    cat(tail(readLines(log, warn = FALSE), 30), sep = "\n")
    stop("Lesson failed; inspect ", log, call. = FALSE)
  }
  cat("  Complete:", output, "\n")
}
source("coloc/validate_results.R")
cat("All six lessons and scientific checks passed.\n")
