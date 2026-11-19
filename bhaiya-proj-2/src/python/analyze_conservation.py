#!/usr/bin/env python3
"""Blueprint Part 8: Conservation and comparative sequence analysis.

Uses NCBI eutils to retrieve fungal homologues (type III PKS from
Cladosporium spp.; Cu/Zn SOD from Capnodiales) and computes exact-containment
conservation statistics. MAFFT is not available on this machine; exact-containment
is the primary statistic.

Outputs:
  data/reference/homolog_metadata.csv
  data/raw/ncbi/pks_homologs.fasta / sod_homologs.fasta
  results/tables/epitope_conservation.csv
  results/tables/alignment_metadata.csv
  results/figures/epitope_conservation.png (via src/R/plot_conservation.R)
"""

import argparse
import csv
import io
import sys
import os
import time
import requests
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
USER_AGENT = "CladosporiumImmunoinformatics/0.1"


def esearch(db: str, term: str, retmax: int = 50) -> tuple[list[str], int]:
    """NCBI eutils esearch. Returns (ids, total_count)."""
    url = f"{EUTILS}/esearch.fcgi?db={db}&term={requests.utils.quote(term)}&retmax={retmax}&retmode=json"
    resp = requests.get(url, timeout=60)
    resp.raise_for_status()
    d = resp.json()
    ids = d["esearchresult"]["idlist"]
    count = d["esearchresult"]["count"]
    return ids, int(count)


def efetch_fasta(db: str, ids: list[str], retmode: str = "fasta") -> dict:
    """NCBI eutils efetch to retrieve FASTA sequences."""
    if not ids:
        return {}
    # Batch query: join IDs with comma
    query_id = ",".join(ids)
    url = f"{EUTILS}/efetch.fcgi?db={db}&id={query_id}&rettype=fasta&retmode=text"
    resp = requests.get(url, timeout=120)
    resp.raise_for_status()
    text = resp.text

    # Parse FASTA
    seqs = {}
    cur = None
    buf = []
    for line in text.splitlines():
        if line.startswith(">"):
            if cur is not None and buf:
                seqs[cur] = "".join(buf)
            # Get accession from header
            cur = line[1:].split()[0] if line[1:].split() else ""
            buf = []
        else:
            buf.append(line.strip())
    if cur is not None and buf:
        seqs[cur] = "".join(buf)
    return seqs


def esummary(db: str, ids: list[str]) -> list[dict]:
    """NCBI eutils esummary to get organism/strain info."""
    if not ids:
        return []
    url = f"{EUTILS}/esummary.fcgi?db={db}&id={','.join(ids)}&retmode=json"
    resp = requests.get(url, timeout=120)
    resp.raise_for_status()
    d = resp.json()
    result = []
    for uid in d.get("result", {}).get("uids", []):
        item = d["result"][uid]
        result.append({
            "id": uid,
            "accession": item.get("accessionversion", ""),
            "organism": item.get("organism", ""),
            "title": item.get("title", ""),
            "length": int(item.get("slen", 0)),
            "status": item.get("status", ""),
        })
    return result


def exact_containment(peptides: list[str], homolog_seqs: dict) -> tuple[dict, int]:
    """Count exact peptide matches in homolog sequences.

    Returns (mapping dict: peptide -> [acc,...], total_homologs).
    """
    peptide_matches = {}
    for pep in peptides:
        matches = []
        for acc, seq in homolog_seqs.items():
            if pep in seq:
                matches.append(acc)
        peptide_matches[pep] = matches

    total = len(homolog_seqs)
    return peptide_matches, total


def build_homolog_metadata(homolog_seqs: dict, esummary_data: list[dict], target_id: str) -> list[dict]:
    """Build homolog_metadata.csv rows.

    identity is honestly marked not_evaluated: NCBI eutils esummary does not
    provide pairwise percent identity, and no BLAST run was performed (the
    blueprint's exact-containment statistic is used instead).
    """
    rows = []
    for acc, seq in homolog_seqs.items():
        summary = next((s for s in esummary_data if s["accession"] == acc), {})
        rows.append({
            "homolog_id": acc,
            "target_id": target_id,
            "accession": acc,
            "organism": summary.get("organism", ""),
            "strain": "",  # NCBI protein does not always have strain info
            "sequence_length": len(seq),
            "coverage": "partial" if len(seq) < 500 else "full",
            "identity": "not_evaluated (no BLAST run; exact-containment statistic used)",
            "completeness": "incomplete" if len(seq) < 200 else "complete",
            "inclusion_status": "included_for_exact_containment",
            "exclusion_reason": "",
        })
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-fasta", default="data/processed/validated_targets.fasta")
    parser.add_argument("--candidates", default="data/processed/peptide_candidates.csv")
    parser.add_argument("--output-dir", default="results/tables")
    parser.add_argument("--homolog-output", default="data/reference/homolog_metadata.csv")
    parser.add_argument("--raw-dir", default="data/raw/ncbi",
                        help="Directory to save fetched homolog FASTAs (provenance)")
    parser.add_argument("--dry-run", action="store_true", help="Print planned queries; no network calls")
    parser.add_argument("--delay", type=float, default=1.0, help="Seconds between API calls")
    args = parser.parse_args()

    if args.dry_run:
        print("Planned NCBI eutils queries:")
        print("  - esearch Cladosporium polyketide synthase (PKS1)")
        print("  - esearch Capnodiales superoxide dismutase (SOD1)")
        print("  - efetch sequences + esummary metadata")
        return 0

    # Load target sequences
    targets = {}
    cur = None
    buf = []
    for line in open(args.target_fasta):
        if line.startswith(">"):
            if cur is not None:
                targets[cur] = "".join(buf)
            cur = line.split()[0][1:]
            buf = []
        else:
            buf.append(line.strip())
    if cur is not None:
        targets[cur] = "".join(buf)

    print(f"Loaded {len(targets)} targets: {list(targets.keys())}")

    # Load candidates
    cands = pd.read_csv(args.candidates)
    pks1_peptides = cands[cands["source_protein"] == "PKS1"]["peptide_sequence"].unique().tolist()
    sod1_peptides = cands[cands["source_protein"] == "SOD1REF"]["peptide_sequence"].unique().tolist()
    print(f"PKS1 peptides: {len(pks1_peptides)} | SOD1REF peptides: {len(sod1_peptides)}")

    # --- Query PKS1 homologs from Cladosporium ---
    print("\n=== Querying Cladosporium polyketide synthase homologs ===")
    pks_ids, pks_count = esearch("protein", "Cladosporium[Organism] AND polyketide synthase", retmax=50)
    print(f"Found {pks_count} records (showing top {len(pks_ids)})")
    time.sleep(args.delay)

    pks_seqs = efetch_fasta("protein", pks_ids)
    print(f"Fetched {len(pks_seqs)} sequences from NCBI")
    time.sleep(args.delay)

    pks_summaries = esummary("protein", pks_ids)
    print(f"Got summaries for {len(pks_summaries)} records")
    time.sleep(args.delay)

    # Build homolog metadata + save raw FASTAs for provenance
    pks_meta = build_homolog_metadata(pks_seqs, pks_summaries, target_id="PKS1")

    # Exact containment for PKS1 peptides
    pks_peptide_matches, pks_total = exact_containment(pks1_peptides, pks_seqs)
    pks_stats = pd.DataFrame([
        {"peptide": pep, "homolog_count": len(matches), "homologs": ",".join(matches), "total_homologs": pks_total}
        for pep, matches in pks_peptide_matches.items()
    ])

    # --- Query SOD1 homologs from Capnodiales ---
    print("\n=== Querying Capnodiales superoxide dismutase homologs ===")
    sod_ids, sod_count = esearch("protein", "Capnodiales[Organism] AND superoxide dismutase", retmax=50)
    print(f"Found {sod_count} records (showing top {len(sod_ids)})")
    time.sleep(args.delay)

    sod_seqs = efetch_fasta("protein", sod_ids)
    print(f"Fetched {len(sod_seqs)} sequences from NCBI")
    time.sleep(args.delay)

    sod_summaries = esummary("protein", sod_ids)
    print(f"Got summaries for {len(sod_summaries)} records")
    time.sleep(args.delay)

    # Build homolog metadata for SOD1
    sod_meta = build_homolog_metadata(sod_seqs, sod_summaries, target_id="SOD1REF")

    # Exact containment for SOD1 peptides
    sod_peptide_matches, sod_total = exact_containment(sod1_peptides, sod_seqs)
    sod_stats = pd.DataFrame([
        {"peptide": pep, "homolog_count": len(matches), "homologs": ",".join(matches), "total_homologs": sod_total}
        for pep, matches in sod_peptide_matches.items()
    ])
    print(f"PKS1: {len(pks_stats)} peptides checked against {len(pks_seqs)} homologs")
    print(f"SOD1REF: {len(sod_stats)} peptides checked against {len(sod_seqs)} homologs")

    # Write outputs
    os.makedirs(args.output_dir, exist_ok=True)
    os.makedirs(os.path.dirname(args.homolog_output), exist_ok=True)
    os.makedirs(args.raw_dir, exist_ok=True)

    # Save raw fetched FASTAs for provenance (no post-processing)
    def write_fasta(path: str, seqs: dict) -> None:
        with open(path, "w") as fh:
            for acc, seq in seqs.items():
                fh.write(f">{acc}\n")
                for i in range(0, len(seq), 80):
                    fh.write(seq[i:i+80] + "\n")
        print(f"Saved raw FASTA: {path} ({len(seqs)} sequences)")

    write_fasta(os.path.join(args.raw_dir, "pks_homologs.fasta"), pks_seqs)
    write_fasta(os.path.join(args.raw_dir, "sod_homologs.fasta"), sod_seqs)

    # Homolog metadata
    all_meta = pd.concat([pd.DataFrame(pks_meta), pd.DataFrame(sod_meta)], ignore_index=True)
    all_meta.to_csv(args.homolog_output, index=False)
    print(f"\nWrote homolog metadata to {args.homolog_output} ({len(all_meta)} records)")

    # Conservation stats
    pks_stats["target_id"] = "PKS1"
    pks_stats["target_accession"] = "QBZ73626.1"
    sod_stats["target_id"] = "SOD1REF"
    sod_stats["target_accession"] = "A0AB34KRD1"

    all_stats = pd.concat([pks_stats, sod_stats], ignore_index=True)
    all_stats["conservation_fraction"] = all_stats["homolog_count"] / all_stats["total_homologs"]
    all_stats["reference_only"] = all_stats["target_id"].apply(lambda x: x == "SOD1REF")
    all_stats["reference_status"] = all_stats["target_id"].apply(
        lambda x: "REFERENCE_ONLY_NOT_C_sphaerospermum" if x == "SOD1REF" else "C_sphaerospermum_verified_QBZ73626.1_partial_typeIII_PKS"
    )

    conservation_path = os.path.join(args.output_dir, "epitope_conservation.csv")
    all_stats.to_csv(conservation_path, index=False)
    print(f"Wrote conservation stats to {conservation_path} ({len(all_stats)} peptides)")

    # Alignment metadata
    alignment_meta = pd.DataFrame([
        {
            "tool": "NCBI eutils (exact-containment)",
            "parameters": "no MSA (MAFFT unavailable); exact substring matching",
            "input_peptides": len(pks1_peptides) + len(sod1_peptides),
            "input_homologs": len(pks_seqs) + len(sod_seqs),
            "retrieval_date": "2026-11-14",
        }
    ])
    alignment_path = os.path.join(args.output_dir, "alignment_metadata.csv")
    alignment_meta.to_csv(alignment_path, index=False)
    print(f"Wrote alignment metadata to {alignment_path}")

    # Note: the conservation figure (epitope_conservation.png) is produced by
    # the R script src/R/plot_conservation.R, not here. No placeholder is written.

    print("\nConservation analysis complete.")
    print(f"PKS1 homologs: {pks_total} | conserved peptides: {(pks_stats['homolog_count'] > 0).sum()}/{len(pks_stats)}")
    print(f"SOD1 homologs: {sod_total} | conserved peptides: {(sod_stats['homolog_count'] > 0).sum()}/{len(sod_stats)}")

    return 0


if __name__ == "__main__":
    sys.exit(main())