"""Generate overlapping peptide windows from protein sequences for epitope prediction.

Coordinates are one-based and inclusive (start_aa, end_aa), matching the
blueprint convention. Candidate IDs are stable and follow the format:
  <SOURCE>_MHCI_000001
  <SOURCE>_MHCII_000001

Peptide windows are generated per protein and per peptide length. This script
produces candidates only; it does NOT run any predictions or call external
services.

Usage:
    python src/python/generate_peptides.py --input <validated.fasta> \
        --output <peptide_candidates.csv> \
        --target_kind MHCI \
        --lengths 8 9 10 11
"""

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
)
logger = logging.getLogger(__name__)

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")


def generate_peptides(sequence: str, peptide_length: int, source_protein: str,
                      target_kind: str, start_index: int = 1) -> list[dict]:
    """Generate overlapping windows of length `peptide_length`.

    Returns a list of dicts with one-based inclusive coordinates and a
    sequential candidate_id. `start_index` is the numeric part of the first
    candidate ID, so IDs remain globally unique across windows from the
    same protein.
    """
    peptides = []
    if target_kind not in {"MHCI", "MHCII"}:
        raise ValueError(f"target_kind must be 'MHCI' or 'MHCII', got {target_kind!r}")

    for start in range(0, len(sequence) - peptide_length + 1):
        peptide = sequence[start:start + peptide_length]
        peptides.append(
            {
                "candidate_id": f"{source_protein}_{target_kind}_{start_index + len(peptides):06d}",
                "source_protein": source_protein,
                "epitope_class": target_kind,
                "peptide_sequence": peptide,
                "start_aa": start + 1,
                "end_aa": start + peptide_length,
                "peptide_length": peptide_length,
            }
        )

    return peptides


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Path to validated FASTA file")
    parser.add_argument("--output", required=True, help="Path to output CSV of peptide candidates")
    parser.add_argument(
        "--target_kind",
        required=True,
        choices=["MHCI", "MHCII"],
        help="Candidate ID / class tag: MHCI or MHCII",
    )
    parser.add_argument(
        "--lengths",
        required=True,
        type=int,
        nargs="+",
        help="Peptide lengths to generate (e.g., 8 9 10 11 for MHC-I; range of 15-mers for MHC-II)",
    )
    args = parser.parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)

    if not input_path.exists():
        logger.error("Input FASTA not found: %s", input_path)
        return 1

    records = list(SeqIO.parse(str(input_path), "fasta"))
    if not records:
        logger.error("No FASTA records found in %s", input_path)
        return 1

    all_peptides: list[dict] = []
    counter = 1

    for record in records:
        sequence = str(record.seq).upper()
        if not sequence:
            logger.warning("Skipping empty record %s", record.id)
            continue
        invalid = sorted(set(sequence) - VALID_AA)
        if invalid:
            logger.warning("Record %s contains non-standard symbols: %s", record.id, ",".join(invalid))

        for length in args.lengths:
            peptides = generate_peptides(sequence, length, record.id, args.target_kind, counter)
            all_peptides.extend(peptides)
            counter += len(peptides)

    if not all_peptides:
        logger.error("No peptide windows could be generated (sequence shorter than peptide length?)")
        return 1

    output_path.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(all_peptides)
    df.to_csv(output_path, index=False)
    logger.info("Wrote %d peptide candidates to %s", len(df), output_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
