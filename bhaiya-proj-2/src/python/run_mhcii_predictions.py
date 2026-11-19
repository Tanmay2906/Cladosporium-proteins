#!/usr/bin/env python3
"""Blueprint Part 5 — MHC-II binding predictions via IEDB (netmhciipan_ba 4.1).

Submits each validated target protein to the documented IEDB tools API
(https://tools.iedb.org/main/tools-api/) with the HLA class-II panel from
config/hla_panel.csv and method/length from config/prediction_methods.yaml.

IEDB allele-name mapping (documented honestly, also noted in the log + docs):
  * HLA-DRB1*01:01 -> "DRB1*01:01"   (IEDB lists DRB1 alleles without the HLA- prefix)
  * HLA-DRB1*04:01 -> "DRB1*04:01"
  * HLA-DRB1*07:01 -> "DRB1*07:01"
  * HLA-DQB1*02:01 -> "DQA1*05:01/DQB1*02:01"  (cis DQ2.5; IEDB MHC-II requires alpha/beta chain pairs)
  * HLA-DQB1*03:01 -> "DQA1*03:01/DQB1*03:01"   (cis DQ7; likewise chain-paired)
Both pairs are present in IEDB's netmhciipan_ba-4.1 allele list (verified 2026-11-14).
The mapping is recorded in iedb_allele and panel_allele columns of the output.

Flow: 2 targets x 5 alleles = 10 calls (length fixed at 15, matching the
candidates table). Raw responses saved verbatim under data/raw/iedb/mhcii/ with
sha256 logs; any failure aborts with NO table written. Returned windows are
cross-checked against peptide_candidates.csv (same coordinates/peptide identity;
every MHC-II candidate must be reproduced) before the joined table is written to
--output (default results/tables/mhcii_all_predictions.csv).

SOD1REF rows are REFERENCE_ONLY_NOT_C_sphaerospermum throughout.
Predictions are computational hypotheses, not evidence of presentation/immunogenicity.
"""

import argparse
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import pandas as pd

import iedb_client

METHOD_CALL = "netmhciipan_ba-4.1"
METHOD_NAME = "netmhciipan_ba"
METHOD_VERSION = "4.1"
SERVICE = "IEDB tools-API (https://tools-cluster-interface.iedb.org/tools_api/mhcii/)"
MHCII_LENGTH = 15

# config/hla_panel.csv panel name -> exact IEDB allele string (verified in IEDB help list 2026-11-14)
ALLELE_MAP = {
    "HLA-DRB1*01:01": "DRB1*01:01",
    "HLA-DRB1*04:01": "DRB1*04:01",
    "HLA-DRB1*07:01": "DRB1*07:01",
    "HLA-DQB1*02:01": "DQA1*05:01/DQB1*02:01",
    "HLA-DQB1*03:01": "DQA1*03:01/DQB1*03:01",
}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--targets", default="data/processed/validated_targets.fasta")
    parser.add_argument("--candidates", default="data/processed/peptide_candidates.csv")
    parser.add_argument("--hla-panel", default="config/hla_panel.csv")
    parser.add_argument("--config", default="config/prediction_methods.yaml")
    parser.add_argument("--raw-dir", default="data/raw/iedb/mhcii")
    parser.add_argument("--log", default="results/logs/mhcii_predictions.log")
    parser.add_argument("--output", default="results/tables/mhcii_all_predictions.csv")
    parser.add_argument("--delay", type=float, default=1.5)
    parser.add_argument("--dry-run", action="store_true", help="Print planned calls; send nothing")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    import csv, yaml
    targets = iedb_client.read_validated_targets(args.targets)
    with open(args.hla_panel) as fh:
        panel = [row["allele"] for row in csv.DictReader(fh) if row["class"].strip().upper() == "II"]
    unmapped = [a for a in panel if a not in ALLELE_MAP]
    if unmapped:
        raise RuntimeError(f"No documented IEDB mapping for panel alleles: {unmapped}; refusing to guess.")
    with open(args.config) as fh:
        cfg = yaml.safe_load(fh)
    n_planned = len(targets) * len(panel)
    print(f"Planned IEDB MHC-II calls: {n_planned} (targets={[t[0] for t in targets]}, "
          f"panel={panel}, length={MHCII_LENGTH}, method={METHOD_CALL})")
    print("Allele mapping:", ALLELE_MAP)

    responses = {}
    for record_id, seq in targets:
        for allele in panel:
            iedb_allele = ALLELE_MAP[allele]
            payload = {"method": METHOD_CALL, "sequence_text": seq, "allele": iedb_allele, "length": str(MHCII_LENGTH)}
            fname = f"{record_id}__{iedb_client._safe_component(iedb_allele)}__L{MHCII_LENGTH}.tsv"
            print(f"CALL {record_id} {allele} -> {iedb_allele} L={MHCII_LENGTH}")
            raw, sha = iedb_client.call_iedb("mhcii", payload, args.log, args.raw_dir, fname,
                                             delay_s=args.delay, dry_run=args.dry_run)
            responses[(record_id, allele)] = (
                iedb_client.parse_tsv_response(raw),
                f"{args.raw_dir}/{fname}",
                sha,
                iedb_allele,
            )

    cand = iedb_client.load_candidates(args.candidates)
    cand_i = cand[(cand["epitope_class"] == "MHCII") & (cand["peptide_length"].astype(int) == MHCII_LENGTH)].reset_index(drop=True)
    index = {(r.source_protein, int(r.start_aa), int(r.end_aa), r.peptide_sequence): r
             for r in cand_i.itertuples()}
    rows, seen = [], set()
    retrieval_date = iedb_client.now_stamp()
    for (record_id, panel_allele), (recs, raw_path, sha, iedb_allele) in responses.items():
        for rec in recs:
            peptide, start, end = rec["peptide"], int(rec["start"]), int(rec["end"])
            key = (record_id, start, end, peptide)
            match = index.get(key)
            if match is None:
                raise RuntimeError(
                    f"IEDB returned window {peptide} at {record_id}:{start}-{end} not found in candidates; aborting."
                )
            seen.add(match.candidate_id)
            seq = dict(targets)[record_id]
            assert seq[start - 1:end] == peptide, f"IEDB peptide {peptide} != local slice {seq[start-1:end]}"
            rows.append({
                "candidate_id": match.candidate_id,
                "source_protein": record_id,
                "epitope_class": "MHCII",
                "peptide": peptide,
                "start_aa": start,
                "end_aa": end,
                "peptide_length": MHCII_LENGTH,
                "allele": panel_allele,
                "panel_allele": panel_allele,
                "iedb_allele": rec.get("allele", iedb_allele),
                "method": METHOD_NAME,
                "method_version": METHOD_VERSION,
                "method_call": METHOD_CALL,
                "core_peptide": rec.get("core_peptide"),
                "ic50_nM": float(rec["ic50"]),
                "rank": float(rec["rank"]),
                "seq_num": rec.get("seq_num"),
                "reference_status": match.reference_status,
                "reference_only": record_id == "SOD1REF",
                "source_service": SERVICE,
                "retrieval_date": retrieval_date,
                "response_file": raw_path,
                "response_sha256": sha,
            })
    missing = set(cand_i["candidate_id"]) - seen
    if missing:
        raise RuntimeError(f"{len(missing)} MHC-II candidates were NOT returned by IEDB "
                           f"(e.g. {sorted(missing)[:5]}); refusing to write a partial table.")
    out = pd.DataFrame(rows)
    out = out.sort_values(["source_protein", "start_aa", "panel_allele"]).reset_index(drop=True)
    iedb_client.safe_write_csv(out, args.output, force=args.force)
    print(f"Wrote {len(out)} MHC-II predictions to {args.output}")
    print(f"Candidates joined: {out['candidate_id'].nunique()} unique / {len(cand_i)} expected")
    return 0


if __name__ == "__main__":
    sys.exit(main())
