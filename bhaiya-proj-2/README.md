# *Cladosporium sphaerospermum* Immunoinformatics — Computational Epitope & Target Characterisation

A computational-only project to identify and triage candidate vaccine/therapeutic epitopes and immunologically relevant protein targets from *Cladosporium sphaerospermum*, starting from verified sequence evidence. Every claim in this repository is a computational hypothesis; **no wet-lab validation exists or is implied**. The only verified *C. sphaerospermum* PKS-class record is QBZ73626.1, a partial type III polyketide synthase (isolate DWEL1; unpublished direct submission MK213995.1); no SOD1 could be verified for the species (0 hits across UniProt/NCBI/local screen), so genus-level references (*C. halotolerans* A0AB34KRD1 / A0AB34KVX2) are used strictly as references — never as targets. Targets, methods (NetMHCpan/NetMHCIIpan BA 4.1 via IEDB, BepiPred 2.0), panel, exclusion policy, sensitivity scenarios and honest gaps are documented in `docs/methodology.md`, `docs/data_provenance.md` and `docs/limitations.md`.

> **No GitHub configured yet — repository not initialised.** No `.git` directory is present (by instruction). Provenance rests on `data/manifests/retrieval_manifest.csv` + file checksums recorded in `docs/data_provenance.md`, not on version control.

## Repository structure

```
bhaiya-proj-2/
├── README.md                        # this file
├── target_selection_rationale.md    # Pin 1 verification report (blueprint Part 1) — read first
├── target_metadata.csv              # target metadata table (4 rows incl. excluded SOD1_Csphaero)
├── config/                          # PROTECTED (other pins): hla_panel.csv, prediction_methods.yaml,
│                                      #   human_proteome.txt (provenance manifest, NOT sequences)
├── data/
│   ├── manifests/retrieval_manifest.csv   # R001–R020 ledger of every external retrieval (2026-11-13)
│   ├── processed/target_evidence_matrix.csv # per-target evidence (domain/localisation/expression/literature)
│   ├── raw/                          # PROTECTED: PKS1.fasta (QBZ73626.1), SOD1.fasta (A0AB34KRD1, reference)
│   └── reference/c_sphaerospermum_ncbi_proteome.fasta  # 324 records, taxid 92950 (R017)
├── docs/                            # blueprint deliverables (Pin 3)
│   ├── methodology.md               # Part 17 Methodological Decisions (D-01…D-11)
│   ├── data_provenance.md           # every dataset/tool/API: dates, URLs, checksums, R-refs, PENDING owners
│   └── limitations.md               # L-01…L-24, stated plainly
├── src/
│   ├── python/                      # PROTECTED (other pins): validate_sequences.py, characterize_proteins.py,
│   │                                #   generate_peptides.py, verify_target_metadata.py
│   └── R/                           # Pin 3 plotting (ggplot2 + dplyr + readr only; no installs)
│       ├── plot_protein_characteristics.R   # needs results/tables/target_protein_characteristics.csv (errors if absent)
│       ├── plot_epitope_predictions.R      # needs mhci/mhcii_all_predictions.csv (SKIPPED if absent)
│       ├── plot_conservation.R             # needs epitope_conservation.csv (SKIPPED if absent)
│       └── plot_candidate_priorities.R     # needs candidate_master_table.csv (SKIPPED if absent)
├── results/
│   ├── tables/putative_sod_homology_screen.csv  # "no hits above threshold" (R020)
│   ├── figures/                     # empty until R scripts run with inputs present
│   └── logs/                        # empty until exclusion log is created at first exclusion
├── notebooks/  tests/  workflows/   # currently empty scaffold dirs
```

## Current status note (2026-11-14)

Two prerequisites have advanced since this README was first written:
- `data/reference/human_proteome.fasta` (UniProtKB/Swiss-Prot 2026_04, 20429 reviewed records, sha256 recorded in `data/reference/human_proteome_manifest.csv`) has been **downloaded** and used for whole-protein and peptide-level similarity screening (see `results/tables/human_protein_similarity.csv` and `candidate_human_similarity.csv`).
- `results/tables/target_protein_characteristics.csv` has been **produced** and `results/figures/protein_characteristics.png` exists.

**IEDB MHC-I/MHC-II predictions were executed 2026-11-14** (user-approved): `src/python/run_mhci_predictions.py` + `run_mhcii_predictions.py` POST the two validated target sequences (the only public data sent) to the documented IEDB tools-API. Raw responses are preserved byte-identically (`data/raw/iedb/`, per-call sha256 in `results/logs/`), and results are tabulated in `results/tables/mhci_all_predictions.csv` (7,440 rows) and `mhcii_all_predictions.csv` (1,805 rows). Row-level cross-checks reproduced every candidate exactly. These are **computational hypotheses only**, not evidence of presentation or immunogenicity. **B-cell (D-07), conservation (D-08), prioritisation (D-09) and sensitivity were also executed 2026-11-14**: BepiPred-2.0 via IEDB (R024), NCBI-eutils exact-containment conservation (R025), RULE-01..05 prioritisation + 8 sensitivity scenarios (R026). Population coverage remains not_evaluated (D-13).

## Quickstart

Prereqs (verified 2026-11-14): Python 3.13.5, Biopython 1.88, pandas 3.1.1; R 4.4.2 with ggplot2, dplyr, readr. **Do not install packages ad hoc; do not run IEDB/NCBI calls not yet approved.**

```bash
cd /Users/aaashish/Projects/bhaiya-proj-2

# 1) Validate sequences (Pin-1 files; read-only for other pins)
python3 src/python/validate_sequences.py --input data/raw/PKS1.fasta --output results/tables/pks1_validation.csv
python3 src/python/validate_sequences.py --input data/raw/SOD1.fasta --output results/tables/sod1_validation.csv

# 2) Characterise proteins (produces results/tables/target_protein_characteristics.csv)
python3 src/python/characterize_proteins.py --input data/raw/PKS1.fasta --output results/tables/target_protein_characteristics.csv

# 3) Generate peptide candidates (no predictions are run by this step)
python3 src/python/generate_peptides.py --input data/raw/PKS1.fasta --output results/tables/pks1_mhci_candidates.csv --target_kind MHCI --lengths 8 9 10 11
python3 src/python/generate_peptides.py --input data/raw/PKS1.fasta --output results/tables/pks1_mhcii_candidates.csv --target_kind MHCII --lengths 15

# 4) Figures (R)
Rscript src/R/plot_protein_characteristics.R    # protein characteristics figure
Rscript src/R/plot_epitope_predictions.R        # MHC-I/MHC-II summary figures
Rscript src/R/plot_bcell_predictions.R          # Fig 4 — BepiPred-2.0 tracks
Rscript src/R/plot_conservation.R               # Fig 5 — conservation histogram (PKS1)
Rscript src/R/plot_population_coverage.R        # prints SKIPPED; no figure (not_evaluated)
Rscript src/R/plot_candidate_priorities.R       # Fig 7 — evidence matrix (top 25)
Rscript src/R/plot_ranking_sensitivity.R        # Fig 8 — ranking sensitivity scenarios
Rscript src/R/plot_workflow.R                   # Fig 1 — workflow diagram
Rscript workflows/run_analysis.R                # driver: runs all of the above
```

Similarity screening, human-proteome retrieval, the IEDB MHC-I/MHC-II prediction runs, B-cell, conservation, prioritisation, sensitivity and the full figure set were executed and are documented/verified (R021–R027). Still pending/not evaluable: NCBI BLAST REST cross-check (D-04, manual step), population coverage (D-13, service 403/405), six-genome mining (PENDING), structural/advanced analyses (not scheduled). See `docs/methodology.md`, `docs/data_provenance.md`, `docs/limitations.md` and `final_compliance_report.md`.

## Blueprint-part completion checklist

Part numbering follows the orchestrator brief; the blueprint document itself is **not stored in this repository** (see `docs/limitations.md` L-21). Checkmarks reflect this repository only.

| Blueprint area (as named in brief/docs) | Status | Evidence |
|---|---|---|
| Part 1 — Target discovery & verification | ✅ done (with documented caveats) | `target_selection_rationale.md`; R001–R020; `target_evidence_matrix.csv` |
| Part 17 — Methodological decisions | ✅ written (this docs set) | `docs/methodology.md`, `docs/data_provenance.md`, `docs/limitations.md` |
| Retrieval ledger / provenance | ✅ R001–R027 (2026-11-13/14) | `data/manifests/retrieval_manifest.csv`; human proteome `HP001` in `data/reference/human_proteome_manifest.csv`; checksums `data/manifests/data_checksums.txt` |
| Sequence validation & characterisation | ✅ run via `workflows/run_core_pipeline.py` (stages 1–4, exit state as logged) | `results/tables/*_validation_report.csv`, `target_protein_characteristics.csv`; 3,702 checks PASS / 0 FAIL |
| Homology screen (local) | ✅ run; 0 hits above threshold | R020; `putative_sod_homology_screen.csv` |
| Human-proteome similarity screen | ✅ run; whole-protein (73 retained) + peptide exact-match (25 hits) | `human_protein_similarity.csv`, `candidate_human_similarity.csv`; R021 |
| NCBI BLAST REST cross-check | ⬜ PENDING (manual browser step; exact-containment conservation used instead) | D-04; D-11; L-28 |
| MHC-I predictions (NetMHCpan BA 4.1) | ✅ executed 2026-11-14 (7,440 rows; all candidates reproduced exactly) | `results/tables/mhci_all_predictions.csv`; R022; D-05/D-12 |
| MHC-II predictions (NetMHCIIpan BA 4.1) | ✅ executed 2026-11-14 (1,805 rows; all candidates reproduced exactly) | `results/tables/mhcii_all_predictions.csv`; R023; D-06/D-12 |
| B-cell (BepiPred 2.0 via IEDB) | ✅ executed 2026-11-14 (389 residue scores; 8 regions) | `results/tables/bcell_residue_scores.csv`, `bcell_predicted_regions.csv`; R024; D-07 |
| Part 7 screening (antigenicity/allergenicity/toxicity) | ⚠️ framework done; all `not_evaluated` (no documented API) | `immunological_screening.csv` (1,857 rows, none negative) |
| Conservation (exact containment) | ✅ executed 2026-11-14 (PKS1 0/1,127 conserved; SOD1REF 350/722) | `epitope_conservation.csv`; `data/reference/homolog_metadata.csv`; R025; D-08 |
| Population coverage (Part 9) | ❌ not_evaluated (service 403/405; nothing fabricated; Fig 6 not produced) | `config/populations.yaml`; `population_coverage_results.csv`; D-13 |
| Prioritisation + exclusion log (Part 10) | ✅ executed 2026-11-14 (1,849 decisions; 378 prioritized; RULE-01..05) | `candidate_master_table.csv`, `candidate_exclusion_log.csv`, `prioritized_candidates.csv`, `candidate_allele_associations.csv`; R026; D-09 |
| Sensitivity analysis (Part 10) | ✅ executed 2026-11-14 (8 scenarios) | `ranking_sensitivity.csv`; D-09; EL leg documented, not run (L-25/L-28) |
| Figures (Part 14) | ✅ 10 real PNGs (Figs 1–5, 7, 8 + MHC-I/II counts/rank); Fig 6 absent by design | `results/figures/`; `src/R/plot_*.R`; R024–R027 |
| Structural analysis | ⬜ not scheduled (unlock = complete cores only) | D-10 |
| Six genome assemblies mining | ⬜ PENDING retrieval/annotation | L-05; D-01 |

Citation placeholder: **not created** (by instruction). No licence file added yet. No external sharing, no publication claim.
