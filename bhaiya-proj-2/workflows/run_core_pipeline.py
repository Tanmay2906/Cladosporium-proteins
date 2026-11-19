#!/usr/bin/env python3
"""Run the core Cladosporium immunoinformatics pipeline.

Stages (all local, no network access):

  1. Sequence validation      -> results/tables/sequence_validation_PKS1.csv
                                 results/tables/sequence_validation_SOD1.csv
  2. Protein characterisation -> results/tables/protein_characteristics_PKS1.csv
                                 results/tables/protein_characteristics_SOD1.csv
  3. Peptide generation       -> data/processed/peptide_candidates.csv
     (MHC-I 8/9/10/11-mers and MHC-II 15-mers, merged from per-kind runs)
     A candidates table that already exists is PRESERVED: the merged output is
     written only when the file is missing or --force is given, so a
     concurrently produced or hand-curated table is never silently clobbered.
     When written, the file is replaced atomically (prior file preserved if
     generation fails).
  4. Output validation        -> results/tables/output_validation_report.csv
     (runs against the candidates table in place)

Prediction stage (blueprint Parts 5-10). The prediction stage CANNOT be run by this
workflow because this workflow is network-free by design. The IEDB MHC-I/MHC-II
prediction runs live in src/python/run_mhci_predictions.py and
src/python/run_mhcii_predictions.py and are invoked EXPLICITLY, only with the
user's approval, because they POST the target sequences to IEDB.

  * If the IEDB output tables (EXPECTED_IEDB_OUTPUTS) are absent, the pipeline
    STOPS after peptide generation with a clear message and exit code 0. Prior
    outputs are preserved; no predictions are invented.
  * If the tables are present, the pipeline performs a lightweight integrity
    check (existence + expected row counts) and STOPS with exit code 2 because
    Parts 6-10 (B-cell, screening, conservation, coverage, prioritisation,
    sensitivity) are not implemented. It never calls the network or creates
    prediction tables.

Note on targets: data/raw/SOD1.fasta is a reference homologue from
Cladosporium halotolerans (see the provenance header in that file), NOT a
verified Cladosporium sphaerospermum target. Peptides derived from it are
reference-derived only and must not be used as C. sphaerospermum predictions.

Logs are appended to results/logs/pipeline.log.

Usage:
    python workflows/run_core_pipeline.py [--log-level INFO]
"""

import argparse
import logging
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
RAW_FASTAS = {
    "PKS1": PROJECT_ROOT / "data" / "raw" / "PKS1.fasta",
    "SOD1": PROJECT_ROOT / "data" / "raw" / "SOD1.fasta",
}
TABLES_DIR = PROJECT_ROOT / "results" / "tables"
LOGS_DIR = PROJECT_ROOT / "results" / "logs"
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
LOG_FILE = LOGS_DIR / "pipeline.log"

MHCI_LENGTHS = [8, 9, 10, 11]
MHCII_LENGTHS = [15]  # config/prediction_methods.yaml: supported range 11-30, default 15

# Expected locations of IEDB prediction exports (netmhcpan_ba / netmhciipan_ba,
# see config/prediction_methods.yaml). The pipeline stops early while these are
# absent; it never creates them.
EXPECTED_IEDB_OUTPUTS = [
    TABLES_DIR / "mhci_all_predictions.csv",
    TABLES_DIR / "mhcii_all_predictions.csv",
]
EXPECTED_IEDB_ROW_COUNTS = {
    "mhci_all_predictions.csv": 7440,   # 1488 MHC-I candidates x 5 class-I alleles
    "mhcii_all_predictions.csv": 1805,  # 361 MHC-II candidates x 5 class-II alleles
}

CHARACTERIZE = PROJECT_ROOT / "src" / "python" / "characterize_proteins.py"
VALIDATE_SEQS = PROJECT_ROOT / "src" / "python" / "validate_sequences.py"
GENERATE_PEPTIDES = PROJECT_ROOT / "src" / "python" / "generate_peptides.py"
VALIDATE_OUTPUTS = PROJECT_ROOT / "src" / "python" / "validate_outputs.py"

# Runtime shim for the pin-owned characterize_proteins.py, which calls two
# ProteinAnalysis APIs that do not exist in the installed biopython 1.88
# ('aliphatic_index' and 'get_amino_acids_percent'; 1.88 exposes
# 'amino_acids_percent' as a property in percent units). The shim is applied
# in the subprocess only; the script itself is NOT modified. Aliphatic index
# is computed per Ikai (1980): AI = X(A) + 2.9*X(V) + 3.9*(X(I) + X(L)),
# with X in mol-%.
_CHARACTERIZE_SHIM = f"""
import sys
from Bio.SeqUtils.ProtParam import ProteinAnalysis

def _get_amino_acids_percent(self):
    return self.amino_acids_percent  # property, percent units summing to 100

def _aliphatic_index(self):
    comp = self.amino_acids_percent  # percent units
    return comp.get("A", 0.0) + 2.9 * comp.get("V", 0.0) + 3.9 * (comp.get("I", 0.0) + comp.get("L", 0.0))

ProteinAnalysis.get_amino_acids_percent = _get_amino_acids_percent
ProteinAnalysis.aliphatic_index = _aliphatic_index

sys.path.insert(0, {str(PROJECT_ROOT / "src" / "python")!r})
sys.argv = ["characterize_proteins.py"] + sys.argv[1:]
import characterize_proteins
sys.exit(characterize_proteins.main())
"""

logger = logging.getLogger("run_core_pipeline")


def configure_logging(log_level: str) -> None:
    LOGS_DIR.mkdir(parents=True, exist_ok=True)
    logger.setLevel(getattr(logging, log_level.upper(), logging.INFO))
    formatter = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s")
    file_handler = logging.FileHandler(LOG_FILE, mode="a", encoding="utf-8")
    file_handler.setFormatter(formatter)
    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setFormatter(formatter)
    logger.handlers = [file_handler, stream_handler]


def prepare_fasta_for_subprocess(src: Path, tmpdir: Path) -> Path:
    """Return a FASTA path the strict biopython 'fasta' parser can read.

    data/raw/SOD1.fasta carries '#' provenance comment lines, which the
    committed validation/characterisation/peptide scripts reject. We write a
    comment-stripped copy into tmpdir; the raw file is never modified.
    """
    text = src.read_text(encoding="utf-8")
    lines = text.splitlines()
    if not any(line.lstrip().startswith("#") for line in lines):
        return src
    cleaned = "\n".join(
        line for line in lines if not line.lstrip().startswith("#")
    ) + "\n"
    dst = tmpdir / (src.stem + ".nocomments.fasta")
    dst.write_text(cleaned, encoding="utf-8")
    logger.warning(
        "%s contains '#' comment lines; using comment-stripped temporary copy %s "
        "for this pipeline run (raw file left untouched)", src, dst,
    )
    return dst


def run_step(label: str, argv: list[str]) -> None:
    logger.info("START %s: %s", label, " ".join(argv))
    result = subprocess.run(argv, cwd=PROJECT_ROOT)
    if result.returncode != 0:
        logger.error("FAILED %s (exit code %d)", label, result.returncode)
        raise RuntimeError(f"{label} failed with exit code {result.returncode}")
    logger.info("DONE %s", label)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--log-level", default="INFO",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        help="Logging verbosity (default: INFO)")
    parser.add_argument("--force", action="store_true",
                        help="Regenerate and atomically replace data/processed/peptide_candidates.csv "
                             "even if it already exists (default: preserve an existing table)")
    args = parser.parse_args()

    configure_logging(args.log_level)
    logger.info("=== Core pipeline run started (root=%s) ===", PROJECT_ROOT)

    for name, path in RAW_FASTAS.items():
        if not path.exists():
            logger.error("Required raw FASTA missing: %s (%s)", name, path)
            return 1

    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    TABLES_DIR.mkdir(parents=True, exist_ok=True)

    try:
        with tempfile.TemporaryDirectory(prefix="pipeline_") as tmp:
            tmpdir = Path(tmp)
            usable = {name: prepare_fasta_for_subprocess(p, tmpdir) for name, p in RAW_FASTAS.items()}

            # --- Stage 1: sequence validation ------------------------------
            for name, path in usable.items():
                run_step(
                    f"sequence validation ({name})",
                    [sys.executable, str(VALIDATE_SEQS),
                     "--input", str(path),
                     "--output", str(TABLES_DIR / f"sequence_validation_{name}.csv")],
                )

            # --- Stage 2: protein characterisation -------------------------
            logger.warning(
                "characterize_proteins.py calls ProteinAnalysis.aliphatic_index / "
                "get_amino_acids_percent, absent in biopython 1.88; applying a runtime "
                "shim (Ikai 1980 aliphatic index) without modifying the script."
            )
            for name, path in usable.items():
                run_step(
                    f"protein characterisation ({name})",
                    [sys.executable, "-c", _CHARACTERIZE_SHIM,
                     "--input", str(path),
                     "--output", str(TABLES_DIR / f"protein_characteristics_{name}.csv")],
                )

            # --- Stage 3: peptide generation (MHC-I + MHC-II, merged) ------
            merge_order = [("PKS1", "MHCI", MHCI_LENGTHS), ("PKS1", "MHCII", MHCII_LENGTHS),
                           ("SOD1", "MHCI", MHCI_LENGTHS), ("SOD1", "MHCII", MHCII_LENGTHS)]
            merged_frames = []
            for name, kind, lengths in merge_order:
                part_csv = tmpdir / f"peptides_{name}_{kind}.csv"
                run_step(
                    f"peptide generation ({name}, {kind}, lengths={lengths})",
                    [sys.executable, str(GENERATE_PEPTIDES),
                     "--input", str(usable[name]),
                     "--output", str(part_csv),
                     "--target_kind", kind,
                     "--lengths", *[str(l) for l in lengths]],
                )
                merged_frames.append(pd.read_csv(part_csv))

            merged = pd.concat(merged_frames, ignore_index=True)
            merged_path = PROCESSED_DIR / "peptide_candidates.csv"

            if merged_path.exists() and not args.force:
                logger.warning(
                    "%s already exists (%d rows); PRESERVING it (regenerate with --force). "
                    "Freshly generated table has %d rows (e.g. %s) and was NOT written. "
                    "Downstream validation will run against the existing file.",
                    merged_path, sum(1 for _ in open(merged_path)) - 1,
                    len(merged), ", ".join(sorted(merged["source_protein"].unique())),
                )
            else:
                atomic_tmp = merged_path.with_suffix(".csv.tmp")
                merged.to_csv(atomic_tmp, index=False)
                os.replace(atomic_tmp, merged_path)  # atomic; prior file preserved on failure
                logger.info("Wrote %d peptide candidates to %s", len(merged), merged_path)

            in_place = pd.read_csv(merged_path)
            if any(str(s).upper().startswith("SOD1") for s in in_place["source_protein"].unique()):
                logger.warning(
                    "Candidates derived from SOD1/SOD1REF are reference-homologue derived only "
                    "(C. halotolerans; see data/raw/SOD1.fasta header) and are NOT "
                    "verified C. sphaerospermum predictions."
                )

            # --- Stage 4: output validation --------------------------------
            run_step(
                "output validation",
                [sys.executable, str(VALIDATE_OUTPUTS),
                 "--candidates", str(merged_path),
                 "--fasta", str(RAW_FASTAS["PKS1"]), str(RAW_FASTAS["SOD1"]),
                 "--output", str(TABLES_DIR / "output_validation_report.csv")],
            )
    except RuntimeError as exc:
        logger.error("Pipeline aborted: %s", exc)
        return 1

    # --- Stage 5: prediction stage (blueprint Parts 5-10) -------------------
    present = [p for p in EXPECTED_IEDB_OUTPUTS if p.exists()]
    absent = [p for p in EXPECTED_IEDB_OUTPUTS if not p.exists()]

    if absent:
        missing_str = ", ".join(str(p.relative_to(PROJECT_ROOT)) for p in absent)
        logger.info(
            "STOP after peptide generation: IEDB prediction output table(s) absent: %s. "
            "No predictions have been run and none are fabricated. The prediction runs "
            "(src/python/run_mhci_predictions.py / run_mhcii_predictions.py, blueprint "
            "Part 5) POST the target sequences to IEDB and must be invoked explicitly, "
            "only with user approval; then re-run this workflow. Prior outputs preserved.",
            missing_str,
        )
        print(f"\nSTOP: prediction stage not run. Absent IEDB prediction output table(s): {missing_str}\n"
              f"Core stages 1-4 completed; see {LOG_FILE} for details.")
        return 0

    integrity_failures = []
    for p in present:
        rows = sum(1 for _ in open(p)) - 1
        expected = EXPECTED_IEDB_ROW_COUNTS.get(p.name)
        status = "ok" if (expected is None or rows == expected) else "MISMATCH"
        logger.info("IEDB output %s present: %d data rows (expected %s) [%s]",
                    p.relative_to(PROJECT_ROOT), rows, expected, status)
        if status != "ok":
            integrity_failures.append((p, rows, expected))

    if integrity_failures:
        details = ", ".join(f"{p.name}: rows={r} expected={e}" for p, r, e in integrity_failures)
        logger.error("STOP: IEDB output table integrity check failed (%s). "
                     "Tables may be truncated/stale; re-run the prediction scripts "
                     "(they refuse to silently overwrite without --force) or inspect manually. "
                     "No downstream work is allowed on inconsistent tables.", details)
        return 3

    # Parts 6-10 have their own explicit, network-separated scripts
    # (run_bcell_predictions.py, screen_candidates.py, analyze_conservation.py,
    # prioritize_candidates.py, sensitivity_analysis.py). This workflow only
    # verifies their outputs are present on disk; it never re-runs them.
    downstream = [
        TABLES_DIR / "bcell_residue_scores.csv",
        TABLES_DIR / "bcell_predicted_regions.csv",
        TABLES_DIR / "immunological_screening.csv",
        TABLES_DIR / "epitope_conservation.csv",
        TABLES_DIR / "candidate_master_table.csv",
        TABLES_DIR / "candidate_exclusion_log.csv",
        TABLES_DIR / "prioritized_candidates.csv",
        TABLES_DIR / "ranking_sensitivity.csv",
    ]
    missing_downstream = [p for p in downstream if not p.exists()]
    if missing_downstream:
        missing_str = ", ".join(str(p.relative_to(PROJECT_ROOT)) for p in missing_downstream)
        logger.info(
            "IEDB MHC-I/MHC-II prediction tables present and internally consistent "
            "(see results/logs/mhci_predictions.log and mhcii_predictions.log for per-call "
            "checksums). This workflow is network-free and will not re-run IEDB. "
            "Parts 6-10 outputs absent (%s); run the explicit Part 6-10 scripts "
            "(see README quickstart). Stopping with exit code 2.",
            missing_str,
        )
        print(f"\nSTOP: prediction tables present and verified; Parts 6-10 outputs absent: {missing_str}\n"
              f"Run the explicit Part 6-10 scripts (README quickstart); details in {LOG_FILE}.")
        return 2

    logger.info(
        "IEDB MHC-I/MHC-II prediction tables present and internally consistent "
        "(see results/logs/mhci_predictions.log and mhcii_predictions.log for per-call "
        "checksums). Parts 6-10 outputs present and verified on disk "
        "(B-cell, screening, conservation, prioritisation, sensitivity; produced by their "
        "own explicit scripts). This workflow is network-free and re-ran nothing external."
    )
    print(f"\nOK: stages 1-4 complete; Parts 5-10 outputs verified on disk.\n"
          f"See results/tables/ and results/figures/; details in {LOG_FILE}.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
