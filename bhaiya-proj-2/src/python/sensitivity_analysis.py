#!/usr/bin/env python3
"""Blueprint Part 10 (sensitivity leg): ranking robustness under prespecified
alternative decision rules.

PRIMARY rule (used for the master triage):
    PKS1 candidates, no human exact match, >=1 panel allele with IC50 < 500 nM
    (MHC-I or MHC-II).

Scenarios (all prespecified, no data-driven tuning):
    STRICT        threshold 50 nM instead of 500 nM
    LENIENT       threshold 5000 nM instead of 500 nM
    RANK05        percentile/rank <= 0.5 instead of IC50 threshold
    NO_DQ         drop DQ alleles from the panel
    NO_DRB101     drop DRB1*01:01 from the panel
    CLASS_I_ONLY  MHC-I evidence only
    CLASS_II_ONLY MHC-II evidence only

For each scenario: retained candidate count and Jaccard overlap of the retained
set with the PRIMARY retained set (plus top-25 Jaccard for ranking stability).

Uses only verified tables on disk; no network calls; no recomputation of
predictions. Missing evidence is never treated as a negative result: candidates
with no prediction for a class are simply not retained via that class.
"""

import sys

import pandas as pd

MHCI = "results/tables/mhci_all_predictions.csv"
MHCII = "results/tables/mhcii_all_predictions.csv"
SIM = "results/tables/candidate_human_similarity.csv"
MASTER = "results/tables/candidate_master_table.csv"
OUT = "results/tables/ranking_sensitivity.csv"

BIND = 500
WEAK = 5000
STRICT = 50
RANK_CUT = 0.5


def binders(threshold, rank_cut=None):
    """Return {candidate_id: set of supporting (class, allele)} per scenario rule."""
    sup = {}
    for path, cls, rankcol in ((MHCI, "MHC-I", "percentile_rank"), (MHCII, "MHC-II", "rank")):
        df = pd.read_csv(path)
        if rank_cut is not None:
            hit = df[df[rankcol] <= rank_cut]
        else:
            hit = df[df["ic50_nM"] < threshold]
        for cid, grp in hit.groupby("candidate_id"):
            sup.setdefault(cid, set()).update((cls, a) for a in grp["allele"].unique())
    return sup


def retained(support, excluded_human, panel=None, classes=None):
    out = []
    for cid, pairs in support.items():
        if cid in excluded_human:
            continue
        if panel is not None and not any(a in panel for _, a in pairs):
            continue
        if classes is not None and not any(c in classes for c, _ in pairs):
            continue
        out.append(cid)
    return set(out)


def jaccard(a, b):
    if not a and not b:
        return 1.0
    return len(a & b) / len(a | b)


def main():
    m1 = pd.read_csv(MHCI)
    m2 = pd.read_csv(MHCII)
    sim = pd.read_csv(SIM, comment="#")
    master = pd.read_csv(MASTER)
    master["best_any_ic50"] = master[["best_mhc_i_ic50", "best_mhc_ii_ic50"]].min(axis=1)

    human_exact = set(sim.loc[sim["human_exact_match_count"] > 0, "candidate_id"])
    pks1_ids = set(master.loc[master["target_protein"] == "PKS1", "candidate_id"])

    alleles_all = sorted(set(m1["allele"]) | set(m2["allele"]))
    dq = {a for a in alleles_all if a.startswith("HLA-DQ")}
    drb101 = {a for a in alleles_all if a.startswith("DRB1*01:01") or a == "DRB1*01:01"}
    class_i = {a for a in m1["allele"].unique()}
    class_ii = {a for a in m2["allele"].unique()}

    scenarios = []
    # PRIMARY: IC50 < 500 nM, full panel, both classes
    s = binders(BIND)
    primary = retained(s, human_exact) & pks1_ids
    scenarios.append(("PRIMARY", "IC50 < 500 nM, full 10-allele panel, MHC-I+MHC-II", primary))

    s = binders(STRICT)
    scenarios.append(("STRICT", "IC50 < 50 nM, full panel, MHC-I+MHC-II",
                      retained(s, human_exact) & pks1_ids))

    s = binders(WEAK)
    scenarios.append(("LENIENT", "IC50 < 5000 nM, full panel, MHC-I+MHC-II",
                      retained(s, human_exact) & pks1_ids))

    s = binders(None, rank_cut=RANK_CUT)
    scenarios.append(("RANK05", "percentile/rank <= 0.5, full panel, MHC-I+MHC-II",
                      retained(s, human_exact) & pks1_ids))

    s = binders(BIND)
    scenarios.append(("NO_DQ", "IC50 < 500 nM, DQ alleles removed",
                      retained(s, human_exact, panel=set(alleles_all) - dq) & pks1_ids))

    s = binders(BIND)
    scenarios.append(("NO_DRB101", "IC50 < 500 nM, DRB1*01:01 removed",
                      retained(s, human_exact, panel=set(alleles_all) - drb101) & pks1_ids))

    s = binders(BIND)
    scenarios.append(("CLASS_I_ONLY", "IC50 < 500 nM, MHC-I evidence only",
                      retained(s, human_exact, classes={"MHC-I"}) & pks1_ids))

    s = binders(BIND)
    scenarios.append(("CLASS_II_ONLY", "IC50 < 500 nM, MHC-II evidence only",
                      retained(s, human_exact, classes={"MHC-II"}) & pks1_ids))

    rows = []
    # Top-25 for every scenario (including PRIMARY): strongest binding evidence
    # among the retained set, deterministic tie-break by candidate_id.
    primary_sorted = master[master["candidate_id"].isin(primary)].sort_values(
        ["best_any_ic50", "candidate_id"])["candidate_id"].tolist()
    top25_primary = set(primary_sorted[:25])
    for name, desc, rset in scenarios:
        # Scenario top-25: strongest binding evidence among that scenario's retained set.
        scen_sorted = master[master["candidate_id"].isin(rset)].sort_values(
            ["best_any_ic50", "candidate_id"])["candidate_id"].tolist()
        top25_scen = set(scen_sorted[:25])
        rows.append({
            "scenario": name,
            "rule": desc,
            "retained_candidates": len(rset),
            "jaccard_vs_PRIMARY_retained": round(jaccard(rset, primary), 4),
            "overlap_with_PRIMARY": len(rset & primary),
            "jaccard_top25_vs_PRIMARY": round(jaccard(top25_scen, top25_primary), 4),
        })

    out = pd.DataFrame(rows)
    out.to_csv(OUT, index=False)
    print(out.to_string(index=False))
    print(f"\nPRIMARY retained: {len(primary)} | wrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
