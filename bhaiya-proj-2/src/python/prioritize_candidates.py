#!/usr/bin/env python3
"""Blueprint Part 10: Integrated candidate prioritisation.

Builds, from the verified stage tables already on disk:

  results/tables/candidate_master_table.csv      (one row per candidate)
  results/tables/candidate_allele_associations.csv  (candidate x allele, Part 10)
  results/tables/candidate_exclusion_log.csv     (every decision + reason)
  results/tables/prioritized_candidates.csv      (deterministic triage order)

Hierarchical decision process (prespecified, NOT a weighted composite score;
blueprint Part 10 'Avoid arbitrary composite scores'):

  RULE-01  source_protein == SOD1REF
           -> exclude_under_prespecified_rule (reference-only C. halotolerans
           homolog; never C. sphaerospermum evidence)
  RULE-02  human_exact_match_count >= 1 (exact match under the D-04 exact-match
           screen)
           -> exclude_under_prespecified_rule (potential autoimmunity screen)
  RULE-03  No allele below 500 nM (MHC-I or MHC-II)
           -> review_required (insufficient evidence; NOT 'failed')
  RULE-04  Best allele below 5000 nM but none below 500 nM
           -> retain_with_caveat
  DEFAULT  At least one allele below 500 nM
           -> retain

Prioritisation is a deterministic triage order (NOT a validated score):
  1. number of panel alleles with IC50 < 500 nM (desc)
  2. best MHC-I percentile_rank / MHC-II rank (asc)
  3. epitope_conservation fraction (desc)
  4. candidate_id (asc; stable tie-break)

Population coverage and antigenicity/allergenicity/toxicity are
`not_evaluated` (documented; blueprint rule: missing output never becomes a
negative result) and therefore cannot exclude any candidate.
No prediction is recomputed here; all evidence comes from the tables listed
in docs/data_provenance.md (R022/R023, HP001, conservation run 2026-11-14).
"""

import argparse
import sys
import os

import pandas as pd

MHCI = "results/tables/mhci_all_predictions.csv"
MHCII = "results/tables/mhcii_all_predictions.csv"
CANDIDATES = "data/processed/peptide_candidates.csv"
SIM = "results/tables/candidate_human_similarity.csv"
CONSERV = "results/tables/epitope_conservation.csv"
BCELL = "results/tables/bcell_predicted_regions.csv"
SCREEN = "results/tables/immunological_screening.csv"

BIND_THRESH_NM = 500
WEAK_THRESH_NM = 5000


def load_tables():
    m1 = pd.read_csv(MHCI)
    m2 = pd.read_csv(MHCII)
    cands = pd.read_csv(CANDIDATES)
    sim = pd.read_csv(SIM, comment="#")
    con = pd.read_csv(CONSERV)
    bcell = pd.read_csv(BCELL)
    screen = pd.read_csv(SCREEN)
    return m1, m2, cands, sim, con, bcell, screen


def build_association_table(m1, m2):
    """One row per (candidate, allele); both classes; score/rank/type preserved."""
    a1 = m1.rename(columns={
        "ic50_nM": "score", "percentile_rank": "rank", "method": "method_name",
    })
    a1 = a1.assign(score_type="ic50_nM", method="netmhcpan_ba", method_version="4.1",
                   score_type_rank="percentile_rank")
    a2 = m2.rename(columns={
        "ic50_nM": "score", "rank": "rank", "method": "method_name",
    })
    a2 = a2.assign(score_type="ic50_nM", method="netmhciipan_ba", method_version="4.1",
                   score_type_rank="rank")
    keep = ["candidate_id", "source_protein", "epitope_class", "peptide", "start_aa",
            "end_aa", "peptide_length", "allele", "method", "method_version",
            "score_type", "score", "rank", "reference_only", "reference_status"]
    assoc = pd.concat([a1[keep], a2[keep]], ignore_index=True)
    return assoc


def bcell_overlap(cands, bcell):
    """Boolean: candidate peptide overlaps any B-cell predicted region of the SAME protein."""
    overlap = {}
    regions = {p: [] for p in bcell.source_protein.unique()}
    for r in bcell.itertuples():
        regions[r.source_protein].append((r.start_aa, r.end_aa))
    for row in cands.itertuples():
        flags = []
        for start, end in regions.get(row.source_protein, []):
            if row.start_aa <= end and row.end_aa >= start:
                flags.append((start, end))
        overlap[row.candidate_id] = bool(flags)
    return overlap


def main():
    m1, m2, cands, sim, con, bcell, screen = load_tables()

    # --- Aggregate best evidence per candidate ---------------------------
    m1_agg = m1.groupby("candidate_id").agg(
        best_mhc_i_ic50=("ic50_nM", "min"),
        best_mhc_i_rank=("percentile_rank", "min"),
        n_alleles_mhc_i_below_500=("ic50_nM", lambda s: (s < BIND_THRESH_NM).sum()),
        n_alleles_mhc_i_below_5000=("ic50_nM", lambda s: (s < WEAK_THRESH_NM).sum()),
        best_allele_mhc_i=("allele", lambda s: s.loc[m1.loc[s.index, "ic50_nM"].idxmin()]),
    )
    m2_agg = m2.groupby("candidate_id").agg(
        best_mhc_ii_ic50=("ic50_nM", "min"),
        best_mhc_ii_rank=("rank", "min"),
        n_alleles_mhc_ii_below_500=("ic50_nM", lambda s: (s < BIND_THRESH_NM).sum()),
        n_alleles_mhc_ii_below_5000=("ic50_nM", lambda s: (s < WEAK_THRESH_NM).sum()),
        best_allele_mhc_ii=("allele", lambda s: s.loc[m2.loc[s.index, "ic50_nM"].idxmin()]),
    )

    # --- Base candidate table -------------------------------------------
    base = cands.rename(columns={"peptide_sequence": "peptide"}).copy()
    base = base.merge(m1_agg, on="candidate_id", how="left").merge(m2_agg, on="candidate_id", how="left")
    base["n_alleles_i_or_ii_below_500"] = (
        base["n_alleles_mhc_i_below_500"].fillna(0) + base["n_alleles_mhc_ii_below_500"].fillna(0)
    ).astype(int)
    base["n_alleles_i_or_ii_below_5000"] = (
        base["n_alleles_mhc_i_below_5000"].fillna(0) + base["n_alleles_mhc_ii_below_5000"].fillna(0)
    ).astype(int)
    base["best_any_ic50"] = base[["best_mhc_i_ic50", "best_mhc_ii_ic50"]].min(axis=1)
    base["best_any_rank"] = base[["best_mhc_i_rank", "best_mhc_ii_rank"]].min(axis=1)

    # --- Join human similarity, conservation, screening, B-cell ---------
    sim_j = sim[["candidate_id", "human_exact_match_count"]].rename(columns={"human_exact_match_count": "human_exact_match_n"})
    base = base.merge(sim_j, on="candidate_id", how="left").fillna({"human_exact_match_n": 0})
    con_j = con[["peptide", "conservation_fraction", "total_homologs"]].rename(
        columns={"total_homologs": "conservation_total_homologs",
                 "conservation_fraction": "conservation_fraction"})
    base = base.merge(con_j, on="peptide", how="left")
    base["conservation_fraction"] = base["conservation_fraction"].fillna(0.0)
    scr = screen.set_index("candidate_id")["final_screening_status"].to_dict()
    base["screening_status"] = base["candidate_id"].map(scr).fillna("not_evaluated")
    overlap = bcell_overlap(cands, bcell)
    base["bcell_predicted_region_overlap"] = base["candidate_id"].map(overlap).fillna(False)
    base["population_coverage_context"] = "not_evaluated"

    # --- Hierarchical decisions (prespecified) ---------------------------
    def decide(row):
        if row["source_protein"] == "SOD1REF":
            return ("exclude_under_prespecified_rule",
                    "RULE-01: C. halotolerans A0AB34KRD1 reference homologue; "
                    "REFERENCE_ONLY_NOT_C_sphaerospermum — not C. sphaerospermum evidence")
        if row["human_exact_match_n"] >= 1:
            return ("exclude_under_prespecified_rule",
                    f"RULE-02: exact match in human reviewed proteome (n={row['human_exact_match_n']}); "
                    "potential autoimmunity screen — excluded under prespecified rule")
        if row["n_alleles_i_or_ii_below_500"] == 0:
            if row["n_alleles_i_or_ii_below_5000"] > 0:
                return ("retain_with_caveat",
                        "RULE-04: no allele below 500 nM; weakest-tier binding below 5000 nM only; "
                        "insufficient evidence for vaccine prioritisation")
            return ("review_required",
                    "RULE-03: no panel allele below 500 nM (MHC-I or MHC-II); insufficient evidence — "
                    "NOT a negative result, review required")
        return ("retain", "RULE-05: at least one panel allele with IC50 < 500 nM; no human exact match; "
                          "reference-authentic C. sphaerospermum candidate")

    dec = base.apply(decide, axis=1, result_type="expand")
    base["decision"] = dec[0]
    base["decision_reason"] = dec[1]

    # --- Deterministic triage order -------------------------------------
    # Sort first, then map the 1-based rank back onto the original frame so the
    # index alignment can never scramble ranks across candidates.
    order = base.sort_values(
        ["n_alleles_i_or_ii_below_500", "best_any_rank", "conservation_fraction", "candidate_id"],
        ascending=[False, True, False, True]
    )["candidate_id"].tolist()
    rank_of = {cid: r + 1 for r, cid in enumerate(order)}
    base["triage_rank"] = base["candidate_id"].map(rank_of).astype(int)

    # --- Outputs ---------------------------------------------------------
    out_dir = "results/tables"

    master_cols = ["candidate_id", "peptide", "source_protein", "start_aa", "end_aa",
                   "peptide_length", "epitope_class",
                   "best_mhc_i_ic50", "best_mhc_i_rank", "best_allele_mhc_i",
                   "best_mhc_ii_ic50", "best_mhc_ii_rank", "best_allele_mhc_ii",
                   "n_alleles_mhc_i_below_500", "n_alleles_mhc_ii_below_500",
                   "n_alleles_i_or_ii_below_500", "n_alleles_i_or_ii_below_5000",
                   "bcell_predicted_region_overlap", "human_exact_match_n",
                   "conservation_fraction", "conservation_total_homologs",
                   "screening_status", "population_coverage_context",
                   "decision", "decision_reason", "triage_rank",
                   "reference_status"]
    master = base[master_cols].rename(columns={"source_protein": "target_protein"})
    master.to_csv(f"{out_dir}/candidate_master_table.csv", index=False)

    # Exclusion log: every candidate, every decision, with rule and reason
    log = base[["candidate_id", "source_protein", "peptide", "epitope_class",
                "decision", "decision_reason", "reference_status"]].copy()
    log["exclusion_rule"] = log["decision_reason"].str.split(":", n=1).str[0]
    log.to_csv(f"{out_dir}/candidate_exclusion_log.csv", index=False)

    # Candidate x allele association table (Part 10: do not discard associations)
    assoc = build_association_table(m1, m2)
    assoc.to_csv(f"{out_dir}/candidate_allele_associations.csv", index=False)

    # Prioritised list: retain/retain_with_caveat, ordered by triage_rank
    prio = base[base["decision"].isin(["retain", "retain_with_caveat"])].sort_values(
        ["triage_rank", "candidate_id"]).reset_index(drop=True)
    prio["triage_rank"] = range(1, len(prio) + 1)
    prio_cols = ["triage_rank",
                 "candidate_id", "peptide", "source_protein", "epitope_class",
                 "start_aa", "end_aa", "peptide_length",
                 "n_alleles_i_or_ii_below_500", "best_any_ic50", "best_any_rank",
                 "conservation_fraction", "human_exact_match_n",
                 "bcell_predicted_region_overlap",
                 "decision", "decision_reason", "reference_status"]
    prio[prio_cols].to_csv(f"{out_dir}/prioritized_candidates.csv", index=False)

    # Summary
    print(f"Master table: {len(master)} rows")
    print(f"Association table: {len(assoc)} rows")
    print(f"Exclusion log: {len(log)} rows")
    print(f"Prioritised (retain/retain_with_caveat): {len(prio)} rows")
    print()
    print("Decisions:")
    print(base["decision"].value_counts().to_string())
    print()
    print("SOD1REF rows must all be excluded:",
          (base[base.source_protein == "SOD1REF"]["decision"] == "exclude_under_prespecified_rule").all())
    return 0


if __name__ == "__main__":
    sys.exit(main())
