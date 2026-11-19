"""Validate FASTA files for protein sequences.

Reads a FASTA file, checks each record, and writes a per-record report with:
  record_id, length, invalid_symbols, is_valid, sha256.

Ambiguous residues (X, B, Z, U, etc.) are reported in `invalid_symbols`
and cause `is_valid` to be False, but they do NOT cause the record to be
silently discarded or modified. Sequences are never trimmed or rewritten.

Usage:
    python src/python/validate_sequences.py --input <fasta_path> --output <report.csv>
"""

import argparse
import hashlib
import logging
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
)
logger = logging.getLogger(__name__)


def calculate_sha256(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("utf-8")).hexdigest()


def validate_fasta(fasta_path: str | Path) -> pd.DataFrame:
    """Return a per-record validation report as a DataFrame."""
    records = list(SeqIO.parse(str(fasta_path), "fasta"))

    if not records:
        raise ValueError(f"No FASTA records found in {fasta_path}")

    results = []
    seen_sequences = {}

    for record in records:
        sequence = str(record.seq).upper()
        invalid = sorted(set(sequence) - VALID_AA)
        duplicates = seen_sequences.get(sequence, [])

        results.append(
            {
                "record_id": record.id,
                "length": len(sequence),
                "invalid_symbols": ",".join(invalid),
                "duplicate_count": len(duplicates),
                "is_valid": len(sequence) > 0 and len(invalid) == 0,
                "sha256": calculate_sha256(sequence),
            }
        )
        seen_sequences.setdefault(sequence, []).append(record.id)

    report = pd.DataFrame(results)
    logger.info("Validated %d records from %s", len(report), fasta_path)
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Path to input FASTA file")
    parser.add_argument("--output", required=True, help="Path to output CSV report")
    args = parser.parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)

    if not input_path.exists():
        logger.error("Input FASTA not found: %s", input_path)
        return 1

    output_path.parent.mkdir(parents=True, exist_ok=True)

    report = validate_fasta(input_path)
    report.to_csv(output_path, index=False)
    logger.info("Wrote validation report to %s", output_path)
    print(report.to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
