#!/usr/bin/env Rscript
# plot_candidate_priorities.R
# Blueprint Part 14 / Figure 7: candidate evidence matrix (top-25 triage).
#   Input:  results/tables/candidate_master_table.csv (Part 10 output)
#   Output: results/figures/candidate_evidence_matrix.png
#
# Uses ONLY ggplot2 + dplyr + readr (+ base R + tidyr). No package installation.
# Never invents data: every cell is read from the master table; predictors that
# were never executed are shown as "not_evaluated", never as negative results.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
})

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
repo_root <- if (length(file_arg) == 1) {
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), "..", ".."), mustWork = FALSE)
} else {
  normalizePath(getwd(), mustWork = TRUE)
}

csv_path <- file.path(repo_root, "results", "tables", "candidate_master_table.csv")
fig_path <- file.path(repo_root, "results", "figures", "candidate_evidence_matrix.png")

if (!file.exists(csv_path)) {
  cat(sprintf("SKIPPED: %s not present\n", csv_path))
  quit(status = 0)
}

m <- read_csv(csv_path, show_col_types = FALSE, progress = FALSE)
names(m) <- tolower(gsub("[^A-Za-z0-9]+", "_", names(m)))

pick <- function(d, cands) { hit <- intersect(cands, names(d)); if (length(hit)) hit[[1]] else NA_character_ }
cid  <- pick(m, c("candidate_id"))
tgt  <- pick(m, c("target_protein", "source_protein"))
trk  <- pick(m, c("triage_rank"))
ni   <- pick(m, c("n_alleles_mhc_i_below_500"))
nii  <- pick(m, c("n_alleles_mhc_ii_below_500"))
bov  <- pick(m, c("bcell_predicted_region_overlap"))
hex  <- pick(m, c("human_exact_match_n"))
cfr  <- pick(m, c("conservation_fraction"))
scr  <- pick(m, c("screening_status"))
pcv  <- pick(m, c("population_coverage_context"))

top <- m %>%
  filter(.data[[tgt]] == "PKS1") %>%
  arrange(.data[[trk]]) %>%
  slice_head(n = 25)

# Class-specific binder columns are NA for wrong-class candidates (structural,
# not missing predictor output): shown as "n/a", never as a negative result.
bind_flag <- function(x) ifelse(is.na(x), "n/a", ifelse(x > 0, "yes", "no"))

ev <- data.frame(
  candidate_id = top[[cid]],
  triage_rank  = top[[trk]],
  `MHC-I binder <500 nM`      = bind_flag(top[[ni]]),
  `MHC-II binder <500 nM`     = bind_flag(top[[nii]]),
  `B-cell region overlap`     = ifelse(as.logical(top[[bov]]), "yes", "no"),
  `Conserved in homologs`     = ifelse(top[[cfr]] > 0, "yes", "no"),
  `Human exact match`         = ifelse(top[[hex]] > 0, "yes", "no"),
  `Antigenicity/allergenicity/toxicity` = "not_evaluated",
  `Population coverage`       = "not_evaluated",
  check.names = FALSE
)

ev_long <- ev %>%
  pivot_longer(cols = -c("candidate_id", "triage_rank"),
               names_to = "evidence", values_to = "status")

ev_long$evidence <- factor(ev_long$evidence, levels = colnames(ev)[-(1:2)])
ev_long$status <- factor(ev_long$status, levels = c("yes", "no", "not_evaluated", "n/a"))

p <- ggplot(ev_long, aes(x = evidence, y = reorder(candidate_id, triage_rank))) +
  geom_tile(aes(fill = status), colour = "white", linewidth = 0.6) +
  geom_text(aes(label = ifelse(status == "not_evaluated", "n/e", "")),
            colour = "grey20", size = 2.6) +
  scale_fill_manual(values = c("yes" = "#1b7837", "no" = "#d9d9d9",
                               "not_evaluated" = "#92c5de", "n/a" = "#f0f0f0"),
                    labels = c("yes" = "evidence present", "no" = "absent",
                               "not_evaluated" = "not evaluated (no predictor)",
                               "n/a" = "not applicable (wrong MHC class)")) +
  labs(
    title = "Candidate evidence matrix — top 25 triage candidates (PKS1)",
    subtitle = "Evidence categories are independent; no weighted composite score (blueprint Part 10 rule)",
    x = NULL, y = "Candidate (triage order)",
    caption = "n/e = not_evaluated: predictor not executable in this environment; missing output is never treated as a negative result."
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1, size = 8.5),
    axis.text.y = element_text(size = 7.5),
    legend.position = "bottom",
    legend.title = element_blank()
  )

ggsave(fig_path, p, width = 11, height = 9, dpi = 300)
cat(sprintf("WROTE: %s\n", fig_path))
