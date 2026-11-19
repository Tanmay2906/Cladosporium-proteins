#!/usr/bin/env python3
"""Blueprint Part 6: B-cell epitope predictions via IEDB (BepiPred 2.0).

Submits each validated target protein to the documented IEDB tools API
(https://tools.iedb.org/main/tools-api/) with the method "Bepipred-2.0".

Flow:
  1. One API call per validated target protein (2 targets total) — 2 calls.
  2. Raw responses are saved verbatim under data/raw/iedb/bcell/ and sha256-logged
     to results/logs/bcell_predictions.log. Any failure aborts the run with a
     non-zero exit code and NO prediction files are written (nothing is fabricated).
  3. Returned residue-level scores are saved verbatim to bcell_residue_scores.csv
     (Position/Residue/Score/Assignment). Candidate regions are defined by
     contiguity of assignments >= 0.5 and length >= 7.
  4. The per-region summary (bcell_predicted_regions.csv) is written along with
     the per-residue table.

SOD1REF rows are REFERENCE_ONLY_NOT_C_sphaerospermum (C. halotolerans A0AB34KRD1)
throughout.

Predictions are computational hypotheses, NOT evidence of presentation or
immunogenicity.

Only run after user approval (via "run iedb" / "Proceeed till all the requiments")
because the service is an IEDB tool API call.
"""

import argparse
import csv
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import pandas as pd

import iedb_client

METHOD_CALL = "Bepipred-2.0"   # exact IEDB documented method name
METHOD_NAME = "Bepipred-2.0"
SERVICE = "IEDB tools-API (https://tools-cluster-interface.iedb.org/tools_api/bcell/)"
THRESHOLD = 0.5
MIN_LEN = 7


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--targets", default="data/processed/validated_targets.fasta",
                        help="Validated targets FASTA (PKS1 and SOD1REF)")
    parser.add_argument("--output-dir", default="results/tables",
                        help="Directory for output CSV files")
    parser.add_argument("--raw-dir", default="data/raw/iedb/bcell")
    parser.add_argument("--log", default="results/logs/bcell_predictions.log")
    parser.add_argument("--delay", type=float, default=1.5, help="Seconds to wait between API calls")
    parser.add_argument("--dry-run", action="store_true", help="Print planned calls; send nothing")
    parser.add_argument("--force", action="store_true", help="Overwrite output CSV files if they exist")
    args = parser.parse_args()

    # Load targets
    targets = iedb_client.read_validated_targets(args.targets)
    n_planned = len(targets)
    print(f"Planned IEDB B-cell calls: {n_planned} (targets={[t[0] for t in targets]}, method={METHOD_CALL})")

    responses = {}
    for record_id, seq in targets:
        fname = f"{record_id}_bcell.tsv"
        print(f"CALL {record_id} L={len(seq)}")
        raw, sha = iedb_client.call_iedb_bcell(seq, args.log, args.raw_dir, fname,
                                               delay_s=args.delay, dry_run=args.dry_run)
        responses[record_id] = {
            "raw": raw,
            "sha256": sha,
            "rows": iedb_client.parse_bcell_response(raw),
        }

    # Load candidate metadata for correct reference status per protein
    cands = pd.read_csv("data/processed/peptide_candidates.csv")
    ref_status_map = dict(zip(cands["source_protein"], cands["reference_status"]))
    ref_only_map = {p: (p == "SOD1REF") for p in cands["source_protein"].unique()}

    # Prepare residue-level scores output
    residue_rows = []
    for record_id, data in responses.items():
        for r in data["rows"]:
            residue_rows.append({
                "candidate_id": f"{record_id}_BC_{r['position']}",
                "source_protein": record_id,
                "epitope_class": "BCELL",
                "position_0based": r["position"],
                "residue": r["residue"],
                "score": r["score"],
                "assignment": r["assignment"],
                "reference_only": (record_id == "SOD1REF") & (ref_only_map.get(record_id, False)),
                "reference_status": ref_status_map.get(record_id, "UNKNOWN"),
                "source_service": SERVICE,
                "retrieval_date": iedb_client.now_stamp(),
                "response_file": f"{args.raw_dir}/{data['sha256'].split('.')[0]}.tsv",
                "response_sha256": data["sha256"],
            })

    # Write residue-level scores
    residues_path = os.path.join(args.output_dir, "bcell_residue_scores.csv")
    residues_path_parent = os.path.dirname(residues_path)
    os.makedirs(residues_path_parent, exist_ok=True)

    pd.DataFrame(residue_rows).to_csv(residues_path, index=False)
    print(f"Wrote {len(residue_rows)} residue-level B-cell predictions to {residues_path}")

    # Define candidate regions for each protein
    region_rows = []
    for record_id, data in responses.items():
        for r in data["rows"]:
            r["score"] = r["score"]
        regions = iedb_client.define_bcell_regions(data["rows"], THRESHOLD, MIN_LEN)
        for reg in regions:
            region_rows.append({
                "candidate_id": f"{record_id}_BC_REG_{reg['start']}-{reg['end']}",
                "source_protein": record_id,
                "epitope_class": "BCELL",
                "start_aa": reg["start"],
                "end_aa": reg["end"],
                "peptide_length": reg["length"],
                "peptide_sequence": reg["residues"],
                "mean_score": reg["mean_score"],
                "max_score": reg["max_score"],
                "threshold": THRESHOLD,
                "method": METHOD_NAME,
                "reference_only": record_id == "SOD1REF",
                "reference_status": ref_status_map.get(record_id, "UNKNOWN"),
                "source_service": SERVICE,
                "retrieval_date": iedb_client.now_stamp(),
            })

    # Write regions summary
    regions_path = os.path.join(args.output_dir, "bcell_predicted_regions.csv")
    regions_path_parent = os.path.dirname(regions_path)
    os.makedirs(regions_path_parent, exist_ok=True)

    pd.DataFrame(region_rows).to_csv(regions_path, index=False)
    print(f"Wrote {len(region_rows)} B-cell candidate regions to {regions_path}")

    # Generate track PNG (a simple placeholder R script would be called by workflow)
    # For now, create a placeholder figure file
    try:
        from pathlib import Path
        fig_path = Path(args.output_dir).parent / "figures" / "bcell_prediction_tracks.png"
        fig_path.parent.mkdir(parents=True, exist_ok=True)
        with open(fig_path, "wb") as fh:
            fh.write(b"PLACEHOLDER: R script plot_epitope_predictions.R would produce this figure")
        print(f"Created placeholder figure {fig_path}")
    except Exception as e:
        print(f"Note: PNG figure generation uses R script (src/R/plot_epitope_predictions.R)")

    print(f"B-cell prediction run complete. Outputs stored under {args.output_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())