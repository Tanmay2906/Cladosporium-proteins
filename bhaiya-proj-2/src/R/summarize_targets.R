#!/usr/bin/env Rscript
# summarize_targets.R
# Blueprint Priority 5 / R analysis: target summary table.
#   Input:  results/tables/target_protein_characteristics.csv
#           data/processed/target_metadata.csv
#   Output: results/tables/target_summary.csv
#
# Uses ONLY readr + dplyr (+ base R). No package installation.
# Never invents data: only columns that exist are read; missing inputs are
# reported and the script exits 0.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
repo_root <- if (length(file_arg) == 1) {
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), "..", ".."), mustWork = FALSE)
} else {
  normalizePath(getwd(), mustWork = TRUE)
}

char_csv <- file.path(repo_root, "results", "tables", "target_protein_characteristics.csv")
meta_csv <- file.path(repo_root, "data", "processed", "target_metadata.csv")
out_csv <- file.path(repo_root, "results", "tables", "target_summary.csv")

if (!file.exists(char_csv)) {
  cat(sprintf("SKIPPED: %s not present\n", char_csv))
  quit(status = 0)
}

char_df <- read_csv(char_csv, show_col_types = FALSE, progress = FALSE) %>%
  select(record_id, length, molecular_weight, isoelectric_point,
         instability_index, aliphatic_index, gravy)

merged <- char_df
if (file.exists(meta_csv)) {
  meta_df <- read_csv(meta_csv, show_col_types = FALSE, progress = FALSE) %>%
    select(target_id, record_id, is_valid, sequence_sha256, reference_status)
  merged <- meta_df %>% left_join(char_df, by = "record_id")
}

write_csv(merged, out_csv)
cat(sprintf("WROTE: %s (%d targets)\n", out_csv, nrow(merged)))
