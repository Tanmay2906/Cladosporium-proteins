"""Cross-check validated pipeline outputs against the raw FASTA sources.

Verifies that every peptide in the peptide-candidates table is internally
consistent and reconstructs *exactly* from the raw FASTA records:

  * required columns are present and no required values are missing;
  * candidate_id is unique;
  * peptide_length equals end_aa - start_aa + 1 (one-based inclusive);
  * start_aa/end_aa are in range for the source protein;
  * the peptide sequence equals the source-protein slice at
    [start_aa..end_aa] taken from the raw FASTA, one-based inclusive.

Source resolution for the reconstruction check is deterministic: first an
exact record-id match across the raw FASTAs, then a filename-stem hint
(e.g. 'PKS1' -> PKS1.fasta, 'SOD1REF' -> SOD1.fasta), which is needed when a
candidates table labels sources differently from the raw record ids. The
report notes which resolution was used for every peptide.

No predictions are generated or read here; this script only checks that the
candidate table faithfully reproduces the raw input sequences. It exits 0 when
every check passes and non-zero (1) on any failure, including a missing or
unreadable candidates file.

Usage:
    python src/python/validate_outputs.py \
        --candidates data/processed/peptide_candidates.csv \
        --fasta data/raw/PKS1.fasta data/raw/SOD1.fasta \
        --output results/tables/output_validation_report.csv
"""

import argparse
import logging
import sys
from io import StringIO
from pathlib import Path

import pandas as pd
from Bio import SeqIO

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
)
logger = logging.getLogger(__name__)

REQUIRED_COLUMNS = [
    "candidate_id",
    "source_protein",
    "epitope_class",
    "peptide_sequence",
    "start_aa",
    "end_aa",
    "peptide_length",
]

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")

REPORT_COLUMNS = ["check", "candidate_id", "source_protein", "status", "detail"]


def parse_fasta_tolerant(path: Path) -> list:
    """Parse a FASTA file, tolerating leading '#' comment lines.

    data/raw/SOD1.fasta carries provenance notes as '#' comment lines, which
    the strict biopython 'fasta' parser rejects. Comment lines are stripped in
    an in-memory copy; the original file is never modified.
    """
    try:
        return list(SeqIO.parse(str(path), "fasta"))
    except ValueError as exc:
        text = path.read_text(encoding="utf-8")
        if not any(line.lstrip().startswith("#") for line in text.splitlines()):
            raise
        logger.warning(
            "%s: %s -- retrying with '#' comment lines stripped (source file untouched)",
            path, exc,
        )
        cleaned = "\n".join(
            line for line in text.splitlines() if not line.lstrip().startswith("#")
        )
        return list(SeqIO.parse(StringIO(cleaned), "fasta"))


def load_source_records(
    fasta_paths: list[Path],
) -> tuple[dict[str, list[tuple[str, str]]], dict[str, list[tuple[str, str]]]]:
    """Parse the raw FASTAs.

    Returns (by_id, by_file):
      by_id   record_id -> list of (sequence, source_file)
      by_file str(path) -> list of (sequence, record_id)
    """
    by_id: dict[str, list[tuple[str, str]]] = {}
    by_file: dict[str, list[tuple[str, str]]] = {}
    for path in fasta_paths:
        if not path.exists():
            logger.error("Source FASTA not found: %s", path)
            continue
        recs = parse_fasta_tolerant(path)
        by_file[str(path)] = [(str(r.seq).upper(), r.id) for r in recs]
        for record in recs:
            by_id.setdefault(record.id, []).append((str(record.seq).upper(), str(path)))
    return by_id, by_file


def resolve_source(
    source_protein: str,
    by_id: dict[str, list[tuple[str, str]]],
    by_file: dict[str, list[tuple[str, str]]],
    fasta_paths: list[Path],
) -> tuple[list[tuple[str, str]], str] | tuple[None, str]:
    """Resolve a candidate's source_protein label to raw-FASTA record(s).

    Resolution order (deterministic, and the report notes which was used):
      1. exact record-id match across the raw FASTAs;
      2. filename-stem hint: the upper-cased source_protein equals or starts
         with a FASTA filename stem (e.g. 'PKS1' -> PKS1.fasta, 'SOD1REF' ->
         SOD1.fasta). The hint is needed because the candidates table may label
         sources differently from raw record ids; no sequence content is
         assumed, the mapping is purely by filename and is logged.

    Returns (sources, resolution_note) where sources is a list of
    (sequence, source_file); or (None, reason) if unresolvable.
    """
    if source_protein in by_id:
        return by_id[source_protein], "exact record id"

    sp_upper = source_protein.strip().upper()
    for path in fasta_paths:
        key = str(path)
        if key not in by_file:
            continue
        stem = path.stem.upper()
        if sp_upper == stem or sp_upper.startswith(stem):
            sources = [(seq, key) for seq, _rid in by_file[key]]
            return sources, f"filename-stem hint {stem!r} -> {key} (label not a raw-fasta record id)"

    return None, f"source_protein {source_protein!r} not resolvable to any provided raw FASTA"


def check_frame(
    df: pd.DataFrame,
    by_id: dict[str, list[tuple[str, str]]],
    by_file: dict[str, list[tuple[str, str]]],
    fasta_paths: list[Path],
) -> list[dict]:
    """Run all checks; return report rows (dicts with REPORT_COLUMNS keys)."""
    rows: list[dict] = []

    def row(check, candidate_id, source_protein, status, detail):
        rows.append(
            {
                "check": check,
                "candidate_id": candidate_id if candidate_id is not None else "GLOBAL",
                "source_protein": source_protein if source_protein is not None else "",
                "status": status,
                "detail": detail,
            }
        )

    # --- Check 1: required columns present ---------------------------------
    missing_cols = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing_cols:
        row(
            "required_columns",
            None,
            "",
            "FAIL",
            f"missing required column(s): {','.join(missing_cols)}",
        )
        # Per-peptide reconstruction is impossible without the required
        # columns; report and stop here.
        return rows
    row("required_columns", None, "", "PASS", f"all required columns present: {','.join(REQUIRED_COLUMNS)}")

    # --- Check 2: no missing required values --------------------------------
    req_df = df[REQUIRED_COLUMNS]
    missing_mask = req_df.isna().any(axis=1)
    empty_mask = req_df.apply(
        lambda col: col.astype(str).str.strip() == ""
    ).any(axis=1)
    missing_mask = missing_mask | empty_mask
    if missing_mask.any():
        bad_ids = df.loc[missing_mask, "candidate_id"].astype(str).tolist()
        row(
            "missing_values",
            None,
            "",
            "FAIL",
            f"{len(bad_ids)} row(s) with missing/empty required values: {','.join(bad_ids[:20])}"
            + ("..." if len(bad_ids) > 20 else ""),
        )
    else:
        row("missing_values", None, "", "PASS", "no missing values in required columns")

    # --- Check 3: candidate_id uniqueness -----------------------------------
    ids = df["candidate_id"].astype(str)
    dupes = sorted(ids[ids.duplicated()].unique())
    if dupes:
        row(
            "unique_candidate_id",
            None,
            "",
            "FAIL",
            f"{len(dupes)} duplicated candidate_id value(s): {','.join(dupes[:20])}"
            + ("..." if len(dupes) > 20 else ""),
        )
    else:
        row("unique_candidate_id", None, "", "PASS", f"all {len(df)} candidate_id values unique")

    # --- Checks 4-7: per-peptide consistency and reconstruction -------------
    for idx, peptide_row in df.iterrows():
        cid = str(peptide_row["candidate_id"])
        sp = str(peptide_row["source_protein"])
        # Rows with missing required values are already reported by the
        # missing_values check; skip them here to avoid cascade noise.
        if any(pd.isna(peptide_row[c]) or str(peptide_row[c]).strip() == "" for c in REQUIRED_COLUMNS):
            continue
        pep = str(peptide_row["peptide_sequence"]).strip().upper()

        try:
            start_aa = int(peptide_row["start_aa"])
            end_aa = int(peptide_row["end_aa"])
            plen = int(peptide_row["peptide_length"])
        except (TypeError, ValueError):
            row("coordinate_types", cid, sp, "FAIL", "start_aa/end_aa/peptide_length are not integers")
            continue

        # length consistency (one-based inclusive)
        if plen == end_aa - start_aa + 1 and len(pep) == plen:
            row("peptide_length_consistency", cid, sp, "PASS",
                f"peptide_length={plen} == end_aa-start_aa+1 == len(peptide_sequence)")
        else:
            row("peptide_length_consistency", cid, sp, "FAIL",
                f"peptide_length={plen}, end_aa-start_aa+1={end_aa - start_aa + 1}, "
                f"len(peptide_sequence)={len(pep)}")

        # invalid symbols in the peptide itself
        bad_symbols = sorted(set(pep) - VALID_AA)
        if bad_symbols:
            row("peptide_alphabet", cid, sp, "FAIL",
                f"non-standard residue(s) in peptide_sequence: {','.join(bad_symbols)}")

        # resolve source protein
        sources, resolution = resolve_source(sp, by_id, by_file, fasta_paths)
        if sources is None:
            row("reconstruct_from_source", cid, sp, "FAIL", resolution)
            continue

        matched = False
        best_detail = ""
        for seq, src_file in sources:
            if start_aa < 1:
                best_detail = f"start_aa={start_aa} < 1 (coordinates are one-based inclusive)"
                continue
            if end_aa > len(seq):
                best_detail = f"end_aa={end_aa} exceeds source length of {len(seq)} aa ({src_file})"
                continue
            fragment = seq[start_aa - 1 : end_aa]
            if fragment == pep:
                matched = True
                best_detail = f"exact reconstruction from {src_file} ({resolution}; {len(seq)} aa)"
                break
            best_detail = (
                f"sequence mismatch at {src_file} ({resolution}): "
                f"expected {fragment!r}, got {pep!r}"
            )

        row("reconstruct_from_source", cid, sp, "PASS" if matched else "FAIL", best_detail)

    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--candidates",
        default="data/processed/peptide_candidates.csv",
        help="Peptide candidates table to validate",
    )
    parser.add_argument(
        "--fasta",
        nargs="+",
        default=["data/raw/PKS1.fasta", "data/raw/SOD1.fasta"],
        help="Raw FASTA source(s) used to reconstruct peptides (one-based inclusive coordinates)",
    )
    parser.add_argument(
        "--output",
        default="results/tables/output_validation_report.csv",
        help="Path to write the validation report CSV",
    )
    args = parser.parse_args()

    candidates_path = Path(args.candidates)
    fasta_paths = [Path(p) for p in args.fasta]
    output_path = Path(args.output)

    by_id, by_file = load_source_records(fasta_paths)
    if not by_file:
        logger.error("No FASTA records could be loaded from: %s", ", ".join(str(p) for p in fasta_paths))

    rows: list[dict]
    if not candidates_path.exists():
        logger.error("Candidates file not found: %s", candidates_path)
        rows = [
            {
                "check": "file_exists",
                "candidate_id": "GLOBAL",
                "source_protein": "",
                "status": "FAIL",
                "detail": f"candidates file not found: {candidates_path}",
            }
        ]
    else:
        try:
            df = pd.read_csv(candidates_path)
        except Exception as exc:  # unreadable/corrupt CSV is a failure, not an exception traceback
            logger.error("Could not read candidates file %s: %s", candidates_path, exc)
            rows = [
                {
                    "check": "file_readable",
                    "candidate_id": "GLOBAL",
                    "source_protein": "",
                    "status": "FAIL",
                    "detail": f"unreadable candidates file: {exc}",
                }
            ]
        else:
            logger.info("Loaded %d candidate rows from %s", len(df), candidates_path)
            rows = [{
                "check": "file_exists",
                "candidate_id": "GLOBAL",
                "source_protein": "",
                "status": "PASS",
                "detail": f"{candidates_path} ({len(df)} rows)",
            }] + check_frame(df, by_id, by_file, fasta_paths)

    report = pd.DataFrame(rows, columns=REPORT_COLUMNS)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    report.to_csv(output_path, index=False)

    n_fail = int((report["status"] == "FAIL").sum())
    n_pass = int((report["status"] == "PASS").sum())
    logger.info("Wrote validation report to %s (%d PASS, %d FAIL)", output_path, n_pass, n_fail)
    print(report["status"].value_counts().to_string())
    print(f"report: {output_path}")

    if n_fail:
        logger.error("Output validation FAILED: %d failing check(s); see %s", n_fail, output_path)
        return 1
    logger.info("Output validation passed: all %d checks OK", n_pass)
    return 0


if __name__ == "__main__":
    sys.exit(main())
