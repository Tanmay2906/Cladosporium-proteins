#!/usr/bin/env Rscript
# plot_bcell_predictions.R
# Blueprint Part 14 / Figure 4: linear B-cell prediction tracks.
#   Input:  results/tables/bcell_residue_scores.csv (BepiPred-2.0 residue scores)
#           results/tables/bcell_predicted_regions.csv (threshold >= 0.5, len >= 7)
#   Output: results/figures/bcell_prediction_tracks.png
#
# Uses ONLY ggplot2 + dplyr + readr (+ base R). No package installation.
# Never invents data: only columns that exist are plotted; missing input
# files are reported and the script exits 0 (no fabricated figure).

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

scores_csv <- file.path(repo_root, "results", "tables", "bcell_residue_scores.csv")
regions_csv <- file.path(repo_root, "results", "tables", "bcell_predicted_regions.csv")
fig_path <- file.path(repo_root, "results", "figures", "bcell_prediction_tracks.png")

if (!file.exists(scores_csv)) {
  cat(sprintf("SKIPPED: %s not present\n", scores_csv))
  quit(status = 0)
}

scores <- read_csv(scores_csv, show_col_types = FALSE, progress = FALSE)
names(scores) <- tolower(gsub("[^A-Za-z0-9]+", "_", names(scores)))

# Required columns, with fallback picks
pick <- function(df, cands) { hit <- intersect(cands, names(df)); if (length(hit)) hit[[1]] else NA_character_ }
pos_col <- pick(scores, c("position_0based", "position", "residue_position"))
score_col <- pick(scores, c("score", "prediction_score", "value"))
prot_col <- pick(scores, c("source_protein", "protein_id", "target"))
if (is.na(pos_col) || is.na(score_col) || is.na(prot_col)) {
  cat("SKIPPED: bcell_residue_scores.csv lacks required columns (position/score/source_protein)\n")
  quit(status = 0)
}

scores$.pos <- as.numeric(scores[[pos_col]]) + 1   # 0-based -> 1-based
scores$.score <- as.numeric(scores[[score_col]])
scores$.protein <- scores[[prot_col]]

# Regions (optional overlay)
regions <- NULL
if (file.exists(regions_csv)) {
  rr <- read_csv(regions_csv, show_col_types = FALSE, progress = FALSE)
  names(rr) <- tolower(gsub("[^A-Za-z0-9]+", "_", names(rr)))
  st <- pick(rr, c("start_aa", "start", "region_start"))
  en <- pick(rr, c("end_aa", "end", "region_end"))
  pr <- pick(rr, c("source_protein", "protein_id", "target"))
  if (!is.na(st) && !is.na(en) && !is.na(pr)) {
    regions <- data.frame(
      protein = rr[[pr]],
      start = as.numeric(rr[[st]]),
      end = as.numeric(rr[[en]]),
      stringsAsFactors = FALSE
    )
  }
}

p <- ggplot(scores, aes(x = .pos, y = .score)) +
  geom_line(colour = "grey55", linewidth = 0.4) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "firebrick", linewidth = 0.5) +
  annotate("text", x = Inf, y = 0.5, label = "threshold 0.5", vjust = -0.6, hjust = 1.1,
           colour = "firebrick", size = 3) +
  labs(
    title = "BepiPred-2.0 linear B-cell prediction tracks",
    subtitle = "Per-residue scores (IEDB tools-API, 2026-11-14); dashed line = 0.5 threshold",
    x = "Residue position (1-based)",
    y = "BepiPred-2.0 score",
    caption = "Computational hypothesis only; not evidence of antibody recognition."
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())

if (!is.null(regions) && nrow(regions) > 0) {
  p <- p + geom_rect(
    data = regions,
    mapping = aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
    fill = "steelblue", alpha = 0.15, inherit.aes = FALSE
  )
}

p <- p + facet_wrap(~ .protein, scales = "free_x", ncol = 1)

ggsave(fig_path, p, width = 10, height = 6, dpi = 300)
cat(sprintf("WROTE: %s\n", fig_path))
