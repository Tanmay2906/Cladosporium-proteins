"""Verify target metadata for selected fungal proteins (skeleton).

Cross-checks a target metadata table (e.g., data/processed/target_metadata.csv)
against minimal integrity rules: accession present, organism/strain attribution
noted, sequence_length consistent with the referenced FASTA record when provided.

This is a SKELETON. It does not perform live literature or database lookups.
It must be completed and verified against actual retrieved records before use
for the PKS1/SOD1 target decision.

Usage:
    python src/python/verify_target_metadata.py --metadata <target_metadata.csv> --output <verification_report.csv> [--fasta <validated_targets.fasta>]
"""

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
)
logger = logging.getLogger(__name__)

REQUIRED_COLUMNS = [
    "target_id",
    "organism",
    "strain",
    "accession",
    "accession_version",
    "sequence_length",
    "completeness",
]


def check_metadata(df: pd.DataFrame) -> pd.DataFrame:
    """Apply skeleton integrity checks; returns per-target pass/fail flags."""
    results = []
    for _, row in df.iterrows():
        missing = [col for col in REQUIRED_COLUMNS if col not in df.columns or pd.isna(row.get(col)) or str(row.get(col)).strip() == ""]
        results.append(
            {
                "target_id": row.get("target_id", "<missing>"),
                "missing_columns": ",".join(missing),
                "is_verified": len(missing) == 0,
                "notes": "" if not missing else "Missing required metadata fields; complete manually after retrieval.",
            }
        )
    return pd.DataFrame(results)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True, help="Path to target metadata CSV")
    parser.add_argument("--output", required=True, help="Path to write verification report CSV")
    parser.add_argument("--fasta", required=False, help="Optional validated FASTA to cross-check sequence lengths")
    args = parser.parse_args()

    metadata_path = Path(args.metadata)
    output_path = Path(args.output)

    if not metadata_path.exists():
        logger.error("Metadata file not found: %s", metadata_path)
        return 1

    df = pd.read_csv(metadata_path)
    report = check_metadata(df)

    # SKELETON: optional cross-check against FASTA lengths is not implemented here.
    # Requires verified target sequences (Pin 1) before it can be trusted.
    if args.fasta:
        logger.warning("--fasta cross-check is a placeholder; implement only after verified targets are confirmed.")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    report.to_csv(output_path, index=False)
    logger.info("Wrote verification report to %s", output_path)
    print(report.to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
