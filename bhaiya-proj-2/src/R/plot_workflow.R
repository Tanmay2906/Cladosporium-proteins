#!/usr/bin/env Rscript
# plot_workflow.R
# Blueprint Part 14 / Figure 1: project workflow and data provenance diagram.
#   Output: results/figures/workflow_diagram.png
#
# Static diagram drawn with ggplot2; content describes the pipeline as actually
# executed (statuses reflect verified on-disk state, not intent). No data beyond
# what is documented in docs/data_provenance.md is drawn. No package installation.

suppressPackageStartupMessages({
  library(ggplot2)
})

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
repo_root <- if (length(file_arg) == 1) {
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), "..", ".."), mustWork = FALSE)
} else {
  normalizePath(getwd(), mustWork = TRUE)
}

fig_path <- file.path(repo_root, "results", "figures", "workflow_diagram.png")

W <- 1.52    # box half-width
H <- 1.42    # box half-height
XPOS <- c(1.72, 5.16, 8.6, 12.04)   # column centres

nodes <- data.frame(
  x = c( XPOS[1], XPOS[2], XPOS[3], XPOS[4],
         XPOS[4], XPOS[3], XPOS[2], XPOS[1] ),
  y = c( 4.55, 4.55, 4.55, 4.55,
         1.35, 1.35, 1.35, 1.35 ),
  label = c(
    "INPUTS\nPKS1 QBZ73626.1 (verified,\n235 aa partial type III PKS)\nSOD1REF A0AB34KRD1\n(REFERENCE_ONLY)\nHLA panel: 10 alleles",
    "Part 3 - CANDIDATES\n1,849 candidates\n(8-11 mer MHC-I,\n15 mer MHC-II)\npeptide_candidates.csv",
    "Part 4 - HUMAN SIMILARITY\nvs Swiss-Prot 2026_04\n(20,429 reviewed records)\n25 exact matches;\n0 binders < 500 nM",
    "Part 5 + 6 - IEDB (executed\n2026-11-14)\nnetmhcpan_ba-4.1 /\nnetmhciipan_ba-4.1 /\nBepiPred-2.0",

    "Part 7 + 8 - GAPS\n(documented)\nPart 7: no documented API\n-> screening = not_evaluated\nPart 8: 0/1,127 PKS1 peptides\nconserved (50 homologs)",
    "Part 9 - NOT EVALUABLE\nPopulation coverage returns\n403 (legacy) / 405 (nextgen)\ndocumented, never\nfabricated",
    "Part 10 - PRIORITISATION\nExclusion log: 1,849 decisions\n378 prioritized: 176 retain +\n202 retain_with_caveat\n(RULE-01..RULE-05)",
    "OUTPUTS\ncandidate_master_table.csv\ncandidate_exclusion_log.csv\ncandidate_allele_associations.csv\nprioritized_candidates.csv\nranking_sensitivity.csv"
  ),
  fill = c("#eef3fa", "#dcebdc", "#dcebdc", "#c7e9c0",
           "#f6e8d8", "#f6d8d8", "#dcebdc", "#eef3fa"),
  colour = c("#8aa6c8", "#4d9b4d", "#4d9b4d", "#1e7a1e",
             "#c88a3a", "#c0392b", "#4d9b4d", "#8aa6c8"),
  stringsAsFactors = FALSE
)

# Snake arrows: left -> right on top row, down at right edge, right -> left on bottom row
arrows <- data.frame(
  x    = c(XPOS[1] + W, XPOS[2] + W, XPOS[3] + W,   XPOS[4],  XPOS[4] - W,  XPOS[3] - W,  XPOS[2] - W),
  y    = c(4.55,        4.55,        4.55,          4.55,      1.35,         1.35,         1.35),
  xend = c(XPOS[2] - W, XPOS[3] - W, XPOS[4] - W,   XPOS[4],   XPOS[3] + W,  XPOS[2] + W,  XPOS[1] + W),
  yend = c(4.55,        4.55,        4.55,          1.35,      1.35,         1.35,         1.35)
)

p <- ggplot() +
  geom_segment(data = arrows,
               aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
               colour = "#444444", linewidth = 0.85) +
  geom_rect(data = nodes,
            aes(xmin = x - W, xmax = x + W, ymin = y - H, ymax = y + H),
            fill = nodes$fill, colour = nodes$colour, linewidth = 0.9) +
  geom_text(data = nodes, aes(x = x, y = y, label = label),
            size = 2.05, lineheight = 1.05, colour = "#1a1a1a") +
  coord_cartesian(xlim = c(0, 13.76), ylim = c(-0.75, 6.85), expand = FALSE) +
  theme_void() +
  labs(
    title = "C. sphaerospermum immunoinformatics workflow — executed state",
    subtitle = "Green = executed & verified  ·  amber = framework / not_evaluated  ·  red = not evaluable (documented)  ·  blue = inputs / outputs",
    caption = "All IEDB calls executed 2026-11-14 with byte-identical raw responses + sha256 logs. SOD1REF is REFERENCE_ONLY_NOT_C_sphaerospermum throughout. Computational hypotheses only; no wet-lab validation. See docs/data_provenance.md."
  ) +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5, margin = margin(b = 2)),
    plot.subtitle = element_text(size = 8, hjust = 0.5, colour = "#444444", margin = margin(b = 5)),
    plot.caption = element_text(size = 7, hjust = 0.5, colour = "#555555", margin = margin(t = 5))
  )

ggsave(fig_path, p, width = 13.76, height = 6.6, dpi = 300)
cat(sprintf("WROTE: %s\n", fig_path))
