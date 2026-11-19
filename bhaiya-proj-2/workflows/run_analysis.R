#!/usr/bin/env Rscript
# workflows/run_analysis.R
# Blueprint Priority 5: R analysis driver. Runs every figure script in order.
# Each script is self-guarding: missing inputs are reported and the script
# exits 0 without fabricating a figure. This driver never patches data.

repo_root <- if (length(commandArgs(trailingOnly = FALSE)) == 0) {
  normalizePath(getwd(), mustWork = TRUE)
} else {
  normalizePath(".", mustWork = TRUE)
}

scripts <- c(
  "src/R/summarize_targets.R",
  "src/R/plot_protein_characteristics.R",
  "src/R/plot_epitope_predictions.R",
  "src/R/plot_bcell_predictions.R",
  "src/R/plot_conservation.R",
  "src/R/plot_population_coverage.R",
  "src/R/plot_candidate_priorities.R",
  "src/R/plot_ranking_sensitivity.R",
  "src/R/plot_workflow.R"
)

for (s in scripts) {
  path <- file.path(repo_root, s)
  cat(sprintf("\n=== %s ===\n", s))
  if (!file.exists(path)) {
    cat(sprintf("SKIPPED: %s not present\n", path))
    next
  }
  sys <- paste("Rscript", shQuote(path))
  rc <- system(sys)
  if (rc != 0) cat(sprintf("WARNING: %s exited with status %d\n", s, rc))
}

cat("\nR analysis driver finished. See results/figures/ for outputs.\n")
