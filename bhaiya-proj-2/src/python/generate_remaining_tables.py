#!/usr/bin/env python3
"""Generate the remaining blueprint-expected tables from verified on-disk data.

No network calls. Every output is derived from files that already exist and
whose integrity was verified (sha256 logs, exact reproduction checks).

Outputs:
  data/processed/target_metadata.csv          (from sequence validation + characteristics)
  data/processed/candidate_master_table.csv  (copy of results/tables master table)
  data/manifests/data_checksums.txt          (sha256 of key inputs/outputs)
  results/tables/mhci_filtered_candidates.csv  (MHC-I binder-tier candidates)
  results/tables/mhcii_filtered_candidates.csv (MHC-II binder-tier candidates)
  results/tables/prediction_failures.csv       (zero rows: no failed submissions)
  results/tables/prediction_run_manifest.csv   (one row per IEDB call, from logs)
  results/tables/population_coverage_input.csv (intended service input)
  results/tables/population_coverage_results.csv (not_evaluated, per blueprint rule)
"""

import hashlib
import os
import re
import shutil
import sys

import pandas as pd

BIND = 500


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def target_metadata():
    rows = []
    for tid, sv, ch in (
        ("PKS1", "results/tables/sequence_validation_PKS1.csv", "results/tables/protein_characteristics_PKS1.csv"),
        ("SOD1REF", "results/tables/sequence_validation_SOD1.csv", "results/tables/protein_characteristics_SOD1.csv"),
    ):
        sv = pd.read_csv(sv).iloc[0]
        ch = pd.read_csv(ch).iloc[0]
        status = ("C_sphaerospermum_verified_QBZ73626.1_partial_typeIII_PKS"
                  if tid == "PKS1" else
                  "REFERENCE_ONLY_NOT_C_sphaerospermum (UniProt A0AB34KRD1, C. halotolerans)")
        rows.append({
            "target_id": tid,
            "record_id": sv["record_id"],
            "length_aa": sv["length"],
            "is_valid": bool(sv["is_valid"]),
            "invalid_symbols": sv["invalid_symbols"] if pd.notna(sv["invalid_symbols"]) else 0,
            "duplicate_count": int(sv["duplicate_count"]),
            "sequence_sha256": sv["sha256"],
            "molecular_weight": ch["molecular_weight"],
            "isoelectric_point": ch["isoelectric_point"],
            "instability_index": ch["instability_index"],
            "aliphatic_index": ch["aliphatic_index"],
            "gravy": ch["gravy"],
            "reference_status": status,
        })
    df = pd.DataFrame(rows)
    df.to_csv("data/processed/target_metadata.csv", index=False)
    print(f"target_metadata.csv: {len(df)} rows")


def filtered_candidates():
    m = pd.read_csv("results/tables/candidate_master_table.csv")
    for cls, out in (("MHCI", "mhci_filtered_candidates.csv"), ("MHCII", "mhcii_filtered_candidates.csv")):
        sub = m[m["epitope_class"] == cls].copy()
        col = "best_mhc_i_ic50" if cls == "MHCI" else "best_mhc_ii_ic50"
        nbelow = "n_alleles_mhc_i_below_500" if cls == "MHCI" else "n_alleles_mhc_ii_below_500"
        rankcol = "best_mhc_i_rank" if cls == "MHCI" else "best_mhc_ii_rank"
        allecol = "best_allele_mhc_i" if cls == "MHCI" else "best_allele_mhc_ii"
        filt = sub[sub[col] < BIND].copy()
        filt = filt[["candidate_id", "target_protein", "peptide", "start_aa", "end_aa",
                     "peptide_length", col, rankcol, allecol, nbelow,
                     "human_exact_match_n", "conservation_fraction", "decision", "reference_status"]]
        filt = filt.rename(columns={col: "best_ic50_nM", rankcol: "best_percentile_rank",
                                    allecol: "best_allele", nbelow: "n_panel_alleles_below_500nM"})
        filt.to_csv(f"results/tables/{out}", index=False)
        print(f"{out}: {len(filt)} binder-tier candidates (IC50 < {BIND} nM)")


def parse_prediction_logs():
    pattern = re.compile(
        r"^(?P<ts>\S+)\t(?P<status>\S+)\t(?P<url>\S+)\tpayload=(?P<payload>\{.*?)\t"
        r"duration_s=(?P<dur>[\d.]+)\thttp=(?P<http>\d+)\tbytes=(?P<bytes>\d+)\t"
        r"sha256=(?P<sha>[0-9a-f]{64})\traw=(?P<raw>\S+)$")
    rows, failures = [], []
    for log, svc in (("results/logs/mhci_predictions.log", "mhci"),
                     ("results/logs/mhcii_predictions.log", "mhcii"),
                     ("results/logs/bcell_predictions.log", "bcell")):
        for line in open(log):
            line = line.rstrip("\n")
            mt = pattern.match(line)
            if not mt:
                # dry-run / header lines are not calls
                continue
            d = mt.groupdict()
            payload = d["payload"]
            allele = re.search(r'"allele": "([^"]+)"', payload)
            length = re.search(r'"length": "([^"]+)"', payload)
            method = re.search(r'"method": "([^"]+)"', payload)
            seq = re.search(r'"sequence_text": "([^"]+)"', payload)
            rows.append({
                "service": svc,
                "endpoint": d["url"],
                "timestamp_utc": d["ts"],
                "status": d["status"],
                "http_code": int(d["http"]),
                "method": method.group(1) if method else "",
                "allele": allele.group(1) if allele else "",
                "length": length.group(1) if length else "",
                "sequence_length_aa": len(seq.group(1)) if seq else 0,
                "duration_s": float(d["dur"]),
                "response_bytes": int(d["bytes"]),
                "response_sha256": d["sha"],
                "raw_response_path": d["raw"],
            })
            if d["status"] != "OK" or int(d["http"]) != 200:
                failures.append(rows[-1])
    man = pd.DataFrame(rows)
    man.to_csv("results/tables/prediction_run_manifest.csv", index=False)
    print(f"prediction_run_manifest.csv: {len(man)} calls "
          f"(mhci={sum(man.service=='mhci')}, mhcii={sum(man.service=='mhcii')}, bcell={sum(man.service=='bcell')})")

    fail_cols = ["service", "endpoint", "timestamp_utc", "status", "http_code", "method",
                 "allele", "length", "failure_reason"]
    fdf = pd.DataFrame(columns=fail_cols)
    for f in failures:
        fdf = pd.concat([fdf, pd.DataFrame([{
            "service": f["service"], "endpoint": f["endpoint"], "timestamp_utc": f["timestamp_utc"],
            "status": f["status"], "http_code": f["http_code"], "method": f["method"],
            "allele": f["allele"], "length": f["length"],
            "failure_reason": "non-OK status or non-200 http code (see log)",
        }])], ignore_index=True)
    fdf.to_csv("results/tables/prediction_failures.csv", index=False)
    print(f"prediction_failures.csv: {len(fdf)} failed submissions")


def population_coverage():
    panel = pd.read_csv("config/hla_panel.csv")
    allele_col = next(c for c in ("allele", "hla_allele", "hla") if c in panel.columns)
    inp = panel.copy()
    inp["service"] = "IEDB population coverage"
    inp["intended_query_status"] = "not_evaluated"
    inp.to_csv("results/tables/population_coverage_input.csv", index=False)
    print(f"population_coverage_input.csv: {len(inp)} panel alleles recorded as intended input")

    res = pd.DataFrame([{
        "analysis": "population_coverage",
        "status": "not_evaluated",
        "reason": "IEDB population-coverage service not scriptable in this environment "
                  "(legacy tools_api/population/ -> HTTP 403; nextgen POST /pipeline -> HTTP 405, "
                  "probed 2026-11-14). Blueprint rule applied: missing output is never treated as a "
                  "negative result and no coverage values are fabricated.",
        "retrieval_date": "",
        "raw_response_path": "",
        "sha256": "",
    }])
    res.to_csv("results/tables/population_coverage_results.csv", index=False)
    print("population_coverage_results.csv: 1 row, status=not_evaluated (no fabricated coverage)")


def main():
    target_metadata()

    # Master table copy for data/processed (blueprint tree location)
    shutil.copyfile("results/tables/candidate_master_table.csv",
                    "data/processed/candidate_master_table.csv")
    print("data/processed/candidate_master_table.csv: copied from results/tables")

    filtered_candidates()
    parse_prediction_logs()
    population_coverage()

    # Checksum manifest of key inputs and outputs
    key = [
        "config/hla_panel.csv", "config/prediction_methods.yaml", "config/targets.yaml",
        "config/filtering_rules.yaml", "config/analysis_settings.yaml",
        "data/raw/PKS1.fasta", "data/raw/SOD1.fasta",
        "data/reference/human_proteome.fasta",
        "data/reference/homolog_metadata.csv",
        "data/processed/validated_targets.fasta",
        "data/processed/peptide_candidates.csv",
        "data/processed/candidate_master_table.csv",
        "data/processed/target_metadata.csv",
        "results/tables/mhci_all_predictions.csv",
        "results/tables/mhcii_all_predictions.csv",
        "results/tables/bcell_residue_scores.csv",
        "results/tables/bcell_predicted_regions.csv",
        "results/tables/candidate_human_similarity.csv",
        "results/tables/immunological_screening.csv",
        "results/tables/epitope_conservation.csv",
        "results/tables/candidate_master_table.csv",
        "results/tables/candidate_allele_associations.csv",
        "results/tables/candidate_exclusion_log.csv",
        "results/tables/prioritized_candidates.csv",
        "results/tables/ranking_sensitivity.csv",
        "results/tables/prediction_run_manifest.csv",
        "results/tables/prediction_failures.csv",
        "results/tables/population_coverage_input.csv",
        "results/tables/population_coverage_results.csv",
    ]
    lines = ["# data_checksums.txt — sha256 of key project inputs and outputs",
             "# generated by src/python/generate_remaining_tables.py",
             "# format: <sha256>  <path>"]
    for p in key:
        if os.path.exists(p):
            lines.append(f"{sha256_file(p)}  {p}")
        else:
            lines.append(f"MISSING  {p}")
    os.makedirs("data/manifests", exist_ok=True)
    with open("data/manifests/data_checksums.txt", "w") as fh:
        fh.write("\n".join(lines) + "\n")
    n = sum(1 for l in lines if not l.startswith("#") and not l.startswith("MISSING"))
    print(f"data/manifests/data_checksums.txt: {n} files checksummed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
