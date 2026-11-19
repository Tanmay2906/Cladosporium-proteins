#!/usr/bin/env Rscript
# plot_ranking_sensitivity.R
# Blueprint Part 14 / Figure 8: ranking sensitivity across prespecified scenarios.
#   Input:  results/tables/ranking_sensitivity.csv (sensitivity_analysis.py)
#   Output: results/figures/ranking_sensitivity.png
#
# Uses ONLY ggplot2 + dplyr + readr (+ base R). No package installation.
# Never invents data: bars are drawn from the scenario table; if the file is
# missing the script reports it and exits 0 (no fabricated figure).

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

csv_path <- file.path(repo_root, "results", "tables", "ranking_sensitivity.csv")
fig_path <- file.path(repo_root, "results", "figures", "ranking_sensitivity.png")

if (!file.exists(csv_path)) {
  cat(sprintf("SKIPPED: %s not present\n", csv_path))
  quit(status = 0)
}

df <- read_csv(csv_path, show_col_types = FALSE, progress = FALSE)
names(df) <- tolower(gsub("[^A-Za-z0-9]+", "_", names(df)))

stopifnot(all(c("scenario", "retained_candidates", "jaccard_vs_primary_retained",
                "jaccard_top25_vs_primary") %in% names(df)))

df$scenario <- factor(df$scenario, levels = df$scenario[order(df$jaccard_vs_primary_retained)])

p <- ggplot(df, aes(x = scenario, y = jaccard_vs_primary_retained, fill = retained_candidates)) +
  geom_col(width = 0.7) +
  geom_point(aes(y = jaccard_top25_vs_primary), shape = 21, size = 3.4,
             fill = "white", colour = "grey25") +
  coord_flip() +
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.25)) +
  scale_fill_gradient(low = "#d7301f", high = "#1b7837", name = "Retained\ncandidates") +
  labs(
    title = "Ranking sensitivity across prespecified decision rules",
    subtitle = "Bars: Jaccard overlap of retained sets vs PRIMARY; circles: top-25 Jaccard",
    x = NULL,
    y = "Jaccard overlap vs PRIMARY (<500 nM full panel)",
    caption = "PRIMARY = IC50 < 500 nM over the full 10-allele panel (MHC-I + MHC-II). Missing predictors never excluded candidates."
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 9))

ggsave(fig_path, p, width = 10, height = 6.2, dpi = 300)
cat(sprintf("WROTE: %s\n", fig_path))
