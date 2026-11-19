#!/usr/bin/env python3
"""Blueprint Part 7: Antigenicity, Allergenicity, Toxicity screening.

Probes for documented APIs (VaxiJen, AlgPred, ToxinPred) — none found
(IEDB tools-API does not expose these predictors; direct web servers have no
documented REST API). Per the blueprint rule, missing predictor output is
recorded as `not_evaluated`, never as a negative result, and no candidate is
excluded on the basis of a missing predictor.

Outputs: results/tables/immunological_screening.csv

Columns (blueprint Part 7 recommended fields):
  candidate_id, antigenicity_method, antigenicity_score, antigenicity_status,
  allergenicity_method, allergenicity_status, toxicity_method, toxicity_status,
  human_similarity_flag, final_screening_status, decision_reason
"""

import argparse
import sys
import os
import pandas as pd

PART = "Part 7"

NOT_EVALUATED = "not_evaluated"
SCREENING_NOTE = (
    "not_evaluated: no documented API for VaxiJen/AlgPred/ToxinPred at IEDB "
    "tools-API or direct endpoints as of 2026-11-14; blueprint rule applied "
    "(missing predictor output is never treated as a negative result)"
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mhci", default="results/tables/mhci_all_predictions.csv")
    parser.add_argument("--mhcii", default="results/tables/mhcii_all_predictions.csv")
    parser.add_argument("--bcell", default="results/tables/bcell_predicted_regions.csv")
    parser.add_argument("--similarity", default="results/tables/candidate_human_similarity.csv")
    parser.add_argument("--output", default="results/tables/immunological_screening.csv")
    args = parser.parse_args()

    # Load candidate IDs from each prediction source
    ids = set()
    for path in [args.mhci, args.mhcii, args.bcell]:
        if os.path.exists(path):
            df = pd.read_csv(path)
            if "candidate_id" in df.columns:
                ids.update(df["candidate_id"].unique())
    if not ids:
        print("No prediction candidates found in any source table.", file=sys.stderr)
        return 2

    # Load human-similarity flags (exact-match screen)
    sim_flags = {}
    if os.path.exists(args.similarity):
        sim = pd.read_csv(args.similarity, comment="#")
        for row in sim.itertuples():
            sim_flags[row.candidate_id] = "exact_match_in_human_proteome" if row.human_exact_match_count > 0 else "no_exact_match"

    rows = []
    for cid in sorted(ids):
        rows.append({
            "candidate_id": cid,
            "antigenicity_method": NOT_EVALUATED,
            "antigenicity_score": None,
            "antigenicity_status": NOT_EVALUATED,
            "allergenicity_method": NOT_EVALUATED,
            "allergenicity_status": NOT_EVALUATED,
            "toxicity_method": NOT_EVALUATED,
            "toxicity_status": NOT_EVALUATED,
            "human_similarity_flag": sim_flags.get(cid, "not_evaluated"),
            "final_screening_status": NOT_EVALUATED,
            "decision_reason": SCREENING_NOTE,
        })

    out = pd.DataFrame(rows)
    out_path = os.path.abspath(args.output)
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    out.to_csv(out_path, index=False)
    print(f"Wrote {len(out)} screening rows to {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
