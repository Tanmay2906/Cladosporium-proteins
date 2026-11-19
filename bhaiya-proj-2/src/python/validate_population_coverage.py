#!/usr/bin/env python3
"""Blueprint Part 9: population coverage validation for selected epitopes.

SCAFFOLD ONLY. This script is not implemented: the screened allele-epitope
table (Part 7) and population allele-frequency data do not exist, so coverage
cannot be computed and none may be fabricated. When implemented, it should
compute population coverage for the selected allele set using the frequency
source designated by the blueprint Part 9 (e.g. AFND), report per-population
coverage to --output, and flag panels falling below the configured threshold.

Exit code 2 (NOT_IMPLEMENTED) until the screened candidates and frequency
data exist.
"""

import argparse
import sys

PART = "Part 9"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--screened", default="data/processed/screened_candidates.csv",
                        help="Screened candidate table (from screen_candidates.py)")
    parser.add_argument("--hla-panel", default="config/hla_panel.csv",
                        help="HLA allele panel used for prediction")
    parser.add_argument("--output", default="results/tables/population_coverage.csv",
                        help="Where to write per-population coverage results")
    args = parser.parse_args()
    print("NOT_IMPLEMENTED: Part 9 awaits targets/config/data: "
          "screened allele-epitope table (Part 7) and population allele-frequency "
          "data are not available; population coverage cannot be computed.")
    return 2


if __name__ == "__main__":
    sys.exit(main())
