#!/usr/bin/env python3
"""Blueprint Part 5: parse and normalise raw IEDB epitope prediction outputs.

SCAFFOLD ONLY. This script is not implemented: the raw IEDB MHC-I/MHC-II
output tables do not exist yet, so there is nothing to parse, and no parsed
table may be fabricated. When implemented, it should read the raw IEDB exports
produced by run_mhci_predictions.py and run_mhcii_predictions.py, normalise
IC50/percentile columns to a common schema, join back to the candidate
metadata (candidate_id, source_protein, start_aa, end_aa), and write the
parsed table to --output for Part 7 screening.

Exit code 2 (NOT_IMPLEMENTED) until the raw IEDB outputs exist.
"""

import argparse
import sys

PART = "Part 5"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mhci", default="results/tables/mhci_predictions.csv",
                        help="Raw IEDB MHC-I output table (from run_mhci_predictions.py)")
    parser.add_argument("--mhcii", default="results/tables/mhcii_predictions.csv",
                        help="Raw IEDB MHC-II output table (from run_mhcii_predictions.py)")
    parser.add_argument("--candidates", default="data/processed/peptide_candidates.csv",
                        help="Candidate metadata to join back onto predictions")
    parser.add_argument("--output", default="data/processed/epitope_predictions_parsed.csv",
                        help="Where to write the normalised prediction table")
    args = parser.parse_args()
    print("NOT_IMPLEMENTED: Part 5 awaits targets/config/data: "
          "raw IEDB MHC-I/MHC-II output tables do not exist "
          "(run_mhci_predictions.py / run_mhcii_predictions.py have not produced them); "
          "nothing to parse.")
    return 2


if __name__ == "__main__":
    sys.exit(main())
