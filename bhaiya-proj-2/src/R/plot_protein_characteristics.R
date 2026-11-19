#!/usr/bin/env Rscript
# plot_protein_characteristics.R
# Reads results/tables/target_protein_characteristics.csv and facets bar charts
# by parameter: length, molecular weight, theoretical pI, instability index, GRAVY.
# Saves results/figures/protein_characteristics.png
#
# Uses ONLY ggplot2 + dplyr + readr (+ base R). Exits with status 1 and a clear
# error message if the input CSV is missing or has no recognisable columns.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
})

# --- Locate repo root from this script's path (src/R/<this>), fallback to cwd ---
args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
repo_root <- if (length(file_arg) == 1) {
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), "..", ".."), mustWork = FALSE)
} else {
  normalizePath(getwd(), mustWork = TRUE)
}

input_csv <- file.path(repo_root, "results", "tables", "target_protein_characteristics.csv")
output_png <- file.path(repo_root, "results", "figures", "protein_characteristics.png")

if (!file.exists(input_csv)) {
  stop(sprintf(
    "ERROR: expected input table not found: %s\nRun src/python/characterize_proteins.py first (see README.md Quickstart).",
    input_csv
  ), call. = FALSE)
}

df <- read_csv(input_csv, show_col_types = FALSE, progress = FALSE)
if (nrow(df) == 0) stop(sprintf("ERROR: %s contains 0 rows.", input_csv), call. = FALSE)

# Normalise column names for tolerant matching across producer versions.
names(df) <- tolower(gsub("[^A-Za-z0-9]+", "_", names(df)))

pick <- function(cands) {
  hit <- intersect(cands, names(df))
  if (length(hit)) hit[[1]] else NA_character_
}

record_col <- pick(c("record_id", "record", "id", "name", "protein", "target", "protein_name"))
param_map <- c(
  "Length (aa)"          = pick(c("length", "n_aa", "sequence_length", "len")),
  "Molecular weight (Da)" = pick(c("molecular_weight", "mw", "molecular_weight_da", "mass", "mol_wt")),
  "Theoretical pI"       = pick(c("isoelectric_point", "pi", "pI", "theoretical_pi")),
  "Instability index"    = pick(c("instability_index", "instability", "instability_idx")),
  "GRAVY"                = pick(c("gravy", "grand_average_hydropathy", "grand_average_of_hydropathy"))
)

if (is.na(record_col)) {
  df$.record <- paste0("record_", seq_len(nrow(df)))
  record_col <- ".record"
  message("NOTE: no recognisable record-id column; using positional labels.")
}

missing <- names(param_map)[is.na(param_map)]
if (length(missing) > 0) {
  message(sprintf("WARNING: no column recognised for: %s — those panels will be absent.",
                  paste(missing, collapse = ", ")))
}
present <- param_map[!is.na(param_map)]
if (length(present) == 0) {
  stop(sprintf("ERROR: none of the expected parameter columns found in %s.\nColumns present: %s",
               input_csv, paste(names(df), collapse = ", ")), call. = FALSE)
}

long <- do.call(rbind, lapply(names(present), function(lbl) {
  data.frame(
    record_id   = as.character(df[[record_col]]),
    parameter   = lbl,
    value       = suppressWarnings(as.numeric(df[[present[[lbl]]]])),
    stringsAsFactors = FALSE
  )
}))
long <- long[!is.na(long$value), , drop = FALSE]
long$parameter <- factor(long$parameter, levels = names(present))
long$record_id <- factor(long$record_id, levels = unique(as.character(df[[record_col]])))

p <- ggplot(long, aes(x = value, y = record_id, fill = record_id)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  facet_wrap(~ parameter, scales = "free", ncol = 3) +
  labs(x = NULL, y = NULL, title = "Protein characteristics by parameter",
       subtitle = basename(input_csv)) +
  theme_bw(base_size = 11) +
  theme(strip.text = element_text(face = "bold"),
        panel.grid.major.y = element_blank())

dir.create(dirname(output_png), recursive = TRUE, showWarnings = FALSE)
ggsave(output_png, plot = p, width = 11, height = 7, dpi = 300)
cat(sprintf("WROTE: %s\n", output_png))
quit(save = "no", status = 0)
