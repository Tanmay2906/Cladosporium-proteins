#!/usr/bin/env Rscript
# plot_conservation.R
# Blueprint Part 14 / Figure 5: candidate conservation.
#   Input:  results/tables/epitope_conservation.csv (exact-containment stats)
#   Output: results/figures/epitope_conservation.png
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

csv_path <- file.path(repo_root, "results", "tables", "epitope_conservation.csv")
fig_path <- file.path(repo_root, "results", "figures", "epitope_conservation.png")

if (!file.exists(csv_path)) {
  cat(sprintf("SKIPPED: %s not present\n", csv_path))
  quit(status = 0)
}

df <- read_csv(csv_path, show_col_types = FALSE, progress = FALSE)
names(df) <- tolower(gsub("[^A-Za-z0-9]+", "_", names(df)))

pick <- function(d, cands) { hit <- intersect(cands, names(d)); if (length(hit)) hit[[1]] else NA_character_ }

pep_col   <- pick(df, c("peptide", "peptide_sequence", "epitope"))
cnt_col   <- pick(df, c("homolog_count", "n_homologs", "count"))
tot_col   <- pick(df, c("total_homologs", "n_total", "total"))
tgt_col   <- pick(df, c("target_id", "target", "protein_id", "source_protein"))
frac_col  <- pick(df, c("conservation_fraction", "fraction", "conserved_fraction"))
ref_col   <- pick(df, c("reference_only", "is_reference"))

if (is.na(pep_col) || is.na(cnt_col) || is.na(tot_col)) {
  cat("SKIPPED: epitope_conservation.csv lacks required columns (peptide/homolog_count/total_homologs)\n")
  quit(status = 0)
}

df$.cnt <- as.numeric(df[[cnt_col]])
df$.tot <- as.numeric(df[[tot_col]])
if (!is.na(frac_col)) {
  df$.frac <- as.numeric(df[[frac_col]])
} else {
  df$.frac <- ifelse(df$.tot > 0, df$.cnt / df$.tot, NA_real_)
}
if (!is.na(tgt_col)) df$.target <- df[[tgt_col]] else df$.target <- "all"
if (!is.na(ref_col)) df$.ref <- as.logical(df[[ref_col]]) else df$.ref <- FALSE

# PKS1 verified target view (C. sphaerospermum); reference-only rows excluded
df_plot <- df[!df$.ref & df$.target == "PKS1", ]
n_pep <- nrow(df_plot)
n_conserved <- sum(df_plot$.cnt > 0)
n_homologs <- max(df_plot$.tot, na.rm = TRUE)
sub_title <- sprintf("PKS1 (C. sphaerospermum verified target) — %d/%d peptides exactly conserved in %d Cladosporium PKS homologues",
                     n_conserved, n_pep, n_homologs)

p <- ggplot(df_plot, aes(x = .frac)) +
  geom_histogram(bins = 30, fill = "darkgreen", colour = "white", alpha = 0.85) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  labs(
    title = "Exact-sequence peptide conservation across fungal homologues",
    subtitle = sub_title,
    x = "Exact-containment conservation fraction (0 = absent, 1 = present in all homologues)",
    y = "Number of peptides",
    caption = "Computational hypothesis only; not evidence of immunogenicity. SOD1REF rows excluded (reference-only)."
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave(fig_path, p, width = 10, height = 6, dpi = 300)
cat(sprintf("WROTE: %s\n", fig_path))
