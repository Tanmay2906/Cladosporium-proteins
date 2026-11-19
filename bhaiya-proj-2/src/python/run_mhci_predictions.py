#!/usr/bin/env python3
"""Blueprint Part 5 — MHC-I binding predictions via IEDB (netmhcpan_ba 4.1).

Submits each validated target protein to the documented IEDB tools API
(https://tools.iedb.org/main/tools-api/) with the HLA class-I panel from
config/hla_panel.csv and method/lengths from config/prediction_methods.yaml.

Flow:
  1. One API call per (validated target, HLA class-I allele, peptide length
     in [8, 9, 10, 11]) — 2 targets x 5 alleles x 4 lengths = 40 calls.
  2. Raw responses are saved verbatim under data/raw/iedb/mhci/ and sha256-logged
     to results/logs/mhci_predictions.log. Any failure aborts the run with a
     non-zero exit code and NO prediction table is written (nothing is fabricated).
  3. Returned windows are cross-checked against data/processed/peptide_candidates.csv
     (same coordinates and peptide identity). Every MHC-I candidate must be
     reproduced exactly or the run fails loudly.
  4. The joined tidy table is written to --output (default
     results/tables/mhci_all_predictions.csv) only after ALL calls and checks pass.

SOD1REF rows are REFERENCE_ONLY_NOT_C_sphaerospermum (C. halotolerans A0AB34KRD1)
throughout and are never reported as C. sphaerospermum data.

Predictions are computational hypotheses, NOT evidence of presentation or
immunogenicity. No IEDB call happens with --dry-run; nothing is sent in that mode.
"""

import argparse
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import pandas as pd

import iedb_client

METHOD_CALL = "netmhcpan_ba-4.1"   # documented method-version syntax (IEDB tools-API)
METHOD_NAME = "netmhcpan_ba"
METHOD_VERSION = "4.1"
SERVICE = "IEDB tools-API (https://tools-cluster-interface.iedb.org/tools_api/mhci/)"
LENGTHS = [8, 9, 10, 11]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--targets", default="data/processed/validated_targets.fasta",
                        help="Validated targets FASTA (PKS1 and SOD1REF)")
    parser.add_argument("--candidates", default="data/processed/peptide_candidates.csv",
                        help="Peptide candidates table for cross-checking/joining")
    parser.add_argument("--hla-panel", default="config/hla_panel.csv")
    parser.add_argument("--config", default="config/prediction_methods.yaml")
    parser.add_argument("--raw-dir", default="data/raw/iedb/mhci")
    parser.add_argument("--log", default="results/logs/mhci_predictions.log")
    parser.add_argument("--output", default="results/tables/mhci_all_predictions.csv")
    parser.add_argument("--delay", type=float, default=1.5, help="Seconds to wait between API calls")
    parser.add_argument("--dry-run", action="store_true", help="Print planned calls; send nothing")
    parser.add_argument("--force", action="store_true", help="Allow overwriting the output table")
    args = parser.parse_args()

    import csv, yaml
    targets = iedb_client.read_validated_targets(args.targets)
    with open(args.hla_panel) as fh:
        alleles = [row["allele"] for row in csv.DictReader(fh) if row["class"].strip().upper() == "I"]
    with open(args.config) as fh:
        cfg = yaml.safe_load(fh)
    lengths = [int(n) for n in cfg["mhci"]["peptide_lengths"]]
    assert lengths == LENGTHS, f"config mhci.peptide_lengths {lengths} != expected {LENGTHS}"
    n_planned = len(targets) * len(alleles) * len(lengths)
    print(f"Planned IEDB MHC-I calls: {n_planned} (targets={[t[0] for t in targets]}, "
          f"alleles={alleles}, lengths={lengths}, method={METHOD_CALL})")

    responses = {}
    for record_id, seq in targets:
        for allele in alleles:
            for length in lengths:
                payload = {"method": METHOD_CALL, "sequence_text": seq, "allele": allele, "length": str(length)}
                fname = f"{record_id}__{iedb_client._safe_component(allele)}__L{length}.tsv"
                print(f"CALL {record_id} {allele} L={length}")
                raw, sha = iedb_client.call_iedb("mhci", payload, args.log, args.raw_dir, fname,
                                                  delay_s=args.delay, dry_run=args.dry_run)
                responses[(record_id, allele, length)] = (
                    iedb_client.parse_tsv_response(raw),
                    f"{args.raw_dir}/{fname}",
                    sha,
                )

    cand = iedb_client.load_candidates(args.candidates)
    cand_i = cand[cand["epitope_class"] == "MHCI"].reset_index(drop=True)
    index = {(r.source_protein, int(r.start_aa), int(r.end_aa), r.peptide_sequence): r
             for r in cand_i.itertuples()}
    # also index by peptide for identity diagnostics
    rows, seen = [], set()
    retrieval_date = iedb_client.now_stamp()
    for (record_id, allele, length), (recs, raw_path, sha) in responses.items():
        for rec in recs:
            peptide, start, end = rec["peptide"], int(rec["start"]), int(rec["end"])
            key = (record_id, start, end, peptide)
            match = index.get(key)
            if match is None:
                raise RuntimeError(
                    f"IEDB returned window {peptide} at {record_id}:{start}-{end} "
                    f"which does NOT match any candidate in {args.candidates}; aborting before writing output."
                )
            if int(match.peptide_length) != length or int(match.end_aa) - int(match.start_aa) + 1 != length:
                raise RuntimeError(f"length mismatch for candidate {match.candidate_id}")
            seen.add(match.candidate_id)
            seq = dict(targets)[record_id]
            assert seq[start - 1:end] == peptide, f"IEDB peptide {peptide} != local slice {seq[start-1:end]}"
            rows.append({
                "candidate_id": match.candidate_id,
                "source_protein": record_id,
                "epitope_class": "MHCI",
                "peptide": peptide,
                "start_aa": start,
                "end_aa": end,
                "peptide_length": length,
                "allele": allele,
                "panel_allele": allele,
                "iedb_allele": rec.get("allele", allele),
                "method": METHOD_NAME,
                "method_version": METHOD_VERSION,
                "method_call": METHOD_CALL,
                "core": rec.get("core"),
                "icore": rec.get("icore"),
                "ic50_nM": float(rec["ic50"]),
                "percentile_rank": float(rec["percentile_rank"]),
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
        raise RuntimeError(f"{len(missing)} MHC-I candidates were NOT returned by IEDB "
                           f"(e.g. {sorted(missing)[:5]}); refusing to write a partial table.")
    out = pd.DataFrame(rows)
    out = out.sort_values(["source_protein", "peptide_length", "start_aa", "allele"]).reset_index(drop=True)
    iedb_client.safe_write_csv(out, args.output, force=args.force)
    print(f"Wrote {len(out)} MHC-I predictions to {args.output}")
    print(f"Candidates joined: {out['candidate_id'].nunique()} unique / {len(cand_i)} expected")
    return 0


if __name__ == "__main__":
    sys.exit(main())
