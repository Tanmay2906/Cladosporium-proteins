#!/usr/bin/env Rscript
# plot_population_coverage.R
# Blueprint Part 14 / Figure 6: population coverage.
#   Input:  results/tables/population_coverage_results.csv
#   Output: results/figures/population_coverage_summary.png (NOT produced)
#
# STATUS: not_evaluated. The IEDB population-coverage service is not scriptable
# in this environment (legacy tools_api/population/ -> HTTP 403; nextgen
# POST /pipeline -> HTTP 405; probed 2026-11-14). Per the blueprint rule, a
# missing predictor output is never treated as a negative result and is never
# fabricated: this script therefore writes NO figure and exits 0.
#
# Uses ONLY ggplot2 + readr. No package installation.

suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
})

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
repo_root <- if (length(file_arg) == 1) {
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), "..", ".."), mustWork = FALSE)
} else {
  normalizePath(getwd(), mustWork = TRUE)
}

csv_path <- file.path(repo_root, "results", "tables", "population_coverage_results.csv")
fig_path <- file.path(repo_root, "results", "figures", "population_coverage_summary.png")

if (!file.exists(csv_path)) {
  cat(sprintf("SKIPPED: %s not present (population coverage not_evaluated)\n", csv_path))
  quit(status = 0)
}

res <- read_csv(csv_path, show_col_types = FALSE, progress = FALSE)
names(res) <- tolower(gsub("[^A-Za-z0-9]+", "_", names(res)))

status_col <- intersect(c("status", "analysis_status", "final_status"), names(res))[[1]]
statuses <- unique(res[[status_col]])

if (all(statuses == "not_evaluated")) {
  cat("SKIPPED: population coverage is not_evaluated (service not scriptable; HTTP 403/405 documented)\n")
  cat(sprintf("         No figure written to %s — fabricating coverage values is prohibited\n", fig_path))
  quit(status = 0)
}

# If a future run supplies real coverage data, render it then.
cat("NOTE: population_coverage_results.csv contains evaluated data; render per blueprint Fig 6.\n")
quit(status = 0)
