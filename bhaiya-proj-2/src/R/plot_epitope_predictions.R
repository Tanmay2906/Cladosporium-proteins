#!/usr/bin/env Rscript
# plot_epitope_predictions.R
# Resilient summary plots for IEDB MHC prediction tables.
#   MHC-I: results/tables/mhci_all_predictions.csv
#   MHC-II: results/tables/mhcii_all_predictions.csv
# Outputs: results/figures/mhci_prediction_overview.png / mhcii_prediction_overview.png
#
# If an input file is missing, prints "SKIPPED: <file> not present" and continues.
# Exits 0 unless an internal error occurs. NEVER invents data: only columns that
# exist in the file are plotted; unrecognisable files are reported, not patched.
#
# Uses ONLY ggplot2 + dplyr + readr (+ base R). No package installation.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
})

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
repo_root <- if (length(file_arg) == 1) {
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), "..", ".."), mustWork = FALSE)
} else {
  normalizePath(getwd(), mustWork = TRUE)
}

fig_dir <- file.path(repo_root, "results", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

pick <- function(df, cands) {
  hit <- intersect(cands, names(df))
  if (length(hit)) hit[[1]] else NA_character_
}

process_one <- function(input_csv, label, out_png) {
  if (!file.exists(input_csv)) {
    cat(sprintf("SKIPPED: %s not present\n", input_csv))
    return(invisible(FALSE))
  }
  df <- read_csv(input_csv, show_col_types = FALSE, progress = FALSE)
  if (nrow(df) == 0) { cat(sprintf("SKIPPED: %s present but contains 0 rows\n", input_csv)); return(invisible(FALSE)) }

  names(df) <- tolower(gsub("[^A-Za-z0-9]+", "_", names(df)))
  allele_col <- pick(df, c("allele", "mhc_allele", "hla", "allele_name", "mhci_allele", "mhcii_allele", "mhc", "locus_allele"))
  pep_col    <- pick(df, c("peptide", "peptide_sequence", "epitope", "pept", "sequence"))
  rank_col   <- pick(df, c("percentile_rank", "percentile", "rank", "percentile_rank_icore", "icore_percentile_rank", "pct_rank"))

  made <- character(0)

  if (!is.na(allele_col)) {
    cnt <- df |> count(.data[[allele_col]], name = "n") |> arrange(desc(n))
    cnt[[allele_col]] <- factor(cnt[[allele_col]], levels = rev(cnt[[allele_col]]))
    p1 <- ggplot(cnt, aes(x = n, y = .data[[allele_col]])) +
      geom_col(fill = "steelblue") +
      geom_text(aes(label = n), hjust = -0.15, size = 3) +
      labs(x = "Number of predictions", y = NULL,
           title = sprintf("%s — predictions per allele", label),
           subtitle = basename(input_csv)) +
      theme_bw(base_size = 11) +
      theme(panel.grid.major.y = element_blank()) +
      expand_limits(x = max(cnt$n) * 1.1)
    fig_counts <- sub("_overview\\.png$", "_counts_per_allele.png", out_png)
    ggsave(fig_counts, p1, width = 9, height = max(5, 0.35 * nrow(cnt) + 1.5), dpi = 300)
    made <- c(made, fig_counts)
  } else {
    cat(sprintf("NOTE: %s — no allele-like column found; skipping per-allele bar chart.\n", basename(input_csv)))
  }

  if (!is.na(rank_col)) {
    df$.rank <- suppressWarnings(as.numeric(df[[rank_col]]))
    df_r <- df[!is.na(df$.rank), , drop = FALSE]
    if (nrow(df_r) > 0) {
      if (!is.na(allele_col)) {
        p2 <- ggplot(df_r, aes(x = .rank)) +
          geom_histogram(bins = 40, fill = "darkorange", colour = "white") +
          facet_wrap(as.formula(paste("~", allele_col)), scales = "free_y") +
          labs(x = "Percentile rank (lower = stronger)", y = "Count",
               title = sprintf("%s — percentile-rank distribution by allele", label),
               subtitle = basename(input_csv)) +
          theme_bw(base_size = 11)
      } else {
        p2 <- ggplot(df_r, aes(x = .rank)) +
          geom_histogram(bins = 40, fill = "darkorange", colour = "white") +
          labs(x = "Percentile rank (lower = stronger)", y = "Count",
               title = sprintf("%s — percentile-rank distribution", label),
               subtitle = basename(input_csv)) +
          theme_bw(base_size = 11)
      }
      fig_rank <- sub("_overview\\.png$", "_rank_distribution.png", out_png)
      ggsave(fig_rank, p2, width = 11, height = 7, dpi = 300)
      made <- c(made, fig_rank)
    } else {
      cat(sprintf("NOTE: %s — rank column '%s' had no numeric values; skipping histogram.\n", basename(input_csv), rank_col))
    }
  } else {
    cat(sprintf("NOTE: %s — no percentile-rank-like column found; skipping histogram.\n", basename(input_csv)))
  }

  if (length(made) == 0) {
    cat(sprintf("SKIPPED: %s — no recognisable allele/percentile columns to plot (never-invent-data policy).\n", input_csv))
    return(invisible(FALSE))
  }
  for (f in made) cat(sprintf("WROTE: %s\n", f))
  invisible(TRUE)
}

process_one(file.path(repo_root, "results", "tables", "mhci_all_predictions.csv"),
            "MHC-I", file.path(fig_dir, "mhci_prediction_overview.png"))
process_one(file.path(repo_root, "results", "tables", "mhcii_all_predictions.csv"),
            "MHC-II", file.path(fig_dir, "mhcii_prediction_overview.png"))
quit(save = "no", status = 0)
