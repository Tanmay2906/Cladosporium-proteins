"""Characterise protein sequences with Bio.SeqUtils.ProtParam (skeleton).

Calculates: sequence length, molecular weight, theoretical pI, instability index,
aliphatic index, GRAVY, and amino-acid composition.

This is a SKELETON. It has NOT been run against the target proteins, and no
target sequence has been confirmed yet. Do not treat outputs as final until
Pin 1 reports verified target sequences.

Usage:
    python src/python/characterize_proteins.py --input <validated.fasta> --output <characteristics.csv>
"""

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
)
logger = logging.getLogger(__name__)


def characterize_record(record):
    sequence = str(record.seq)
    analysis = ProteinAnalysis(sequence)

    return {
        "record_id": record.id,
        "length": len(sequence),
        "molecular_weight": analysis.molecular_weight(),
        "isoelectric_point": analysis.isoelectric_point(),
        "instability_index": analysis.instability_index(),
        "aliphatic_index": analysis.aliphatic_index(),
        "gravy": analysis.gravy(),
        "amino_acid_composition": analysis.get_amino_acids_percent(),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Path to validated FASTA file")
    parser.add_argument("--output", required=True, help="Path to output CSV characteristics table")
    args = parser.parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)

    if not input_path.exists():
        logger.error("Input FASTA not found: %s", input_path)
        return 1

    records = list(SeqIO.parse(str(input_path), "fasta"))
    if not records:
        logger.error("No FASTA records in %s", input_path)
        return 1

    results = [characterize_record(r) for r in records]
    df = pd.DataFrame(results)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False)
    logger.info("Wrote characteristics for %d records to %s", len(df), output_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
