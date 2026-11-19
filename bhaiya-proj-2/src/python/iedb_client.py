#!/usr/bin/env python3
"""Shared IEDB tools-API client for the Cladosporium immunoinformatics project.

Implements the documented IEDB Analysis Resource RESTful interface
(https://tools.iedb.org/main/tools-api/, accessed 2026-11-14) used by
src/python/run_mhci_predictions.py and src/python/run_mhcii_predictions.py.

Non-negotiable rules for this module:
  * Raw responses are saved byte-for-byte and never edited.
  * Any failure (non-200 HTTP, transport error, garbled body) fails the run
    loudly; NO prediction values are ever synthesized or interpolated.
  * The final joined table is written only after every call in the panel
    succeeded AND cross-checks against data/processed/peptide_candidates.csv
    pass exactly.
"""

import hashlib
import json
import os
import sys
import time
from datetime import datetime, timezone

import requests

API_BASE = "https://tools-cluster-interface.iedb.org/tools_api"
USER_AGENT = "CladosporiumImmunoinformatics/0.1 (local computational research)"
DEFAULT_TIMEOUT = 120
MAX_RETRIES = 2  # retries on transport/5xx errors only; 4xx fails immediately
RETRY_BACKOFF_S = 5.0


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _safe_component(text: str) -> str:
    return (
        text.replace("/", "_").replace("*", "").replace(":", "")
        .replace(" ", "_").replace("-", "")
    )


def call_iedb(endpoint: str, payload: dict, log_path: str,
              raw_dir: str, raw_filename: str,
              delay_s: float = 1.5, dry_run: bool = False) -> tuple[bytes, str]:
    """POST to the IEDB tools API with logging; persist raw response bytes.

    Returns (response_bytes, sha256_hex). Raises RuntimeError on any failure;
    callers must abort the run and must not fabricate missing values.
    """
    os.makedirs(raw_dir, exist_ok=True)
    os.makedirs(os.path.dirname(log_path), exist_ok=True)
    target = f"{API_BASE}/{endpoint.lstrip('/')}/"
    if dry_run:
        msg = f"[DRY-RUN] would POST {target} with payload={payload}"
        print(msg)
        with open(log_path, "a") as fh:
            fh.write(f"{utc_now_iso()}\tDRY_RUN\t{json.dumps(payload, sort_keys=True)}\n")
        raise SystemExit(0)

    last_err = None
    for attempt in range(1, MAX_RETRIES + 2):
        started = time.time()
        try:
            resp = requests.post(
                target,
                data=payload,
                headers={"User-Agent": USER_AGENT},
                timeout=DEFAULT_TIMEOUT,
            )
            duration = time.time() - started
            body = resp.content
            sha = hashlib.sha256(body).hexdigest()
            if resp.status_code == 200 and body.strip():
                raw_path = os.path.join(raw_dir, raw_filename)
                with open(raw_path, "wb") as fh:
                    fh.write(body)
                with open(log_path, "a") as fh:
                    fh.write(
                        f"{utc_now_iso()}\tOK\t{target}\tpayload={json.dumps(payload, sort_keys=True)}"
                        f"\tduration_s={duration:.1f}\thttp={resp.status_code}\tbytes={len(body)}\tsha256={sha}\traw={raw_path}\n"
                    )
                if delay_s > 0:
                    time.sleep(delay_s)
                return body, sha
            last_err = RuntimeError(
                f"HTTP {resp.status_code} for {target} (attempt {attempt}): {body[:500]!r}"
            )
        except requests.RequestException as exc:  # transport-level failure
            duration = time.time() - started
            last_err = RuntimeError(f"transport error for {target} (attempt {attempt}): {exc}")
            with open(log_path, "a") as fh:
                fh.write(f"{utc_now_iso()}\tERROR\t{target}\tpayload={json.dumps(payload, sort_keys=True)}\terror={exc}\n")
        if attempt <= MAX_RETRIES:
            time.sleep(RETRY_BACKOFF_S * attempt)
    raise last_err


def parse_tsv_response(raw: bytes) -> list[dict]:
    """Parse an IEDB TSV response. Raises ValueError if the shape is unexpected."""
    text = raw.decode("utf-8", errors="strict")
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if len(lines) < 2:
        raise ValueError(f"IEDB response has no data rows (got {len(lines)} lines): {text[:300]!r}")
    header = lines[0].split("\t")
    rows = []
    for ln in lines[1:]:
        parts = ln.split("\t")
        if len(parts) != len(header):
            raise ValueError(f"row/column mismatch: expected {len(header)} columns, got {len(parts)}: {ln[:300]!r}")
        rows.append(dict(zip(header, parts)))
    return rows


# --------------------------------------------------------------------------- #
# B-cell (BepiPred 2.0) — same IEDB tools-API contract, different endpoint     #
# --------------------------------------------------------------------------- #
BCALL_METHOD = "Bepipred-2.0"
BCALL_THRESHOLD = 0.5  # BepiPred 2.0 recommended threshold (0 = non-epitope, 1 = epitope)
BCALL_MIN_REGION_LEN = 7  # BepiPred 2.0 recommended minimum contiguous epitope length


def call_iedb_bcell(sequence: str, log_path: str, raw_dir: str,
                    raw_filename: str, delay_s: float = 1.5,
                    dry_run: bool = False) -> tuple[bytes, str]:
    """POST a single protein to IEDB's B-cell endpoint (BepiPred 2.0).

    Returns (raw_response_bytes, sha256_hex). Same non-negotiable rules as
    call_iedb(): raw bytes are saved verbatim; any failure aborts the run;
    nothing is synthesized. The IEDB B-cell endpoint returns a TSV with
    Position(0-based) / Residue / Score / Assignment columns.
    """
    payload = {"method": BCALL_METHOD, "sequence_text": sequence}
    return call_iedb("bcell", payload, log_path, raw_dir, raw_filename,
                     delay_s=delay_s, dry_run=dry_run)


def parse_bcell_response(raw: bytes) -> list[dict]:
    """Parse a BepiPred-2.0 TSV response.

    Columns (verified 2026-11-14): Position (0-based), Residue, Score,
    Assignment ('E' = predicted epitope, '.' = non-epitope).
    """
    text = raw.decode("utf-8", errors="strict")
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if len(lines) < 2:
        raise ValueError(f"B-cell response has no data rows: {text[:300]!r}")
    header = lines[0].split("\t")
    if header != ["Position", "Residue", "Score", "Assignment"]:
        raise ValueError(f"unexpected B-cell response header: {header!r}; {text[:300]!r}")
    rows = []
    for ln in lines[1:]:
        parts = ln.split("\t")
        if len(parts) != 4:
            raise ValueError(f"B-cell row/column mismatch: {ln[:300]!r}")
        rows.append({
            "position": int(parts[0]),
            "residue": parts[1],
            "score": float(parts[2]),
            "assignment": parts[3],
        })
    return rows


def define_bcell_regions(rows: list[dict], threshold: float,
                         min_len: int) -> list[dict]:
    """Combine contiguous residues at/above threshold into candidate regions.

    Returns a list of dicts with start (1-based), end (1-based), length,
    mean_score, max_score, residues, and the original residue-level rows.
    Residues below threshold break contiguity. Residues exactly at
    threshold are included (>= threshold).
    """
    regions = []
    current = []
    for r in rows:
        if r["score"] >= threshold:
            current.append(r)
        else:
            if len(current) >= min_len:
                regions.append(current)
            current = []
    if len(current) >= min_len:
        regions.append(current)
    out = []
    for cur in regions:
        scores = [r["score"] for r in cur]
        out.append({
            "start": cur[0]["position"] + 1,  # 0-based -> 1-based
            "end": cur[-1]["position"] + 1,
            "length": len(cur),
            "mean_score": round(sum(scores) / len(scores), 4),
            "max_score": round(max(scores), 4),
            "residues": "".join(r["residue"] for r in cur),
            "n_residues": len(cur),
        })
    return out


def read_validated_targets(path: str) -> list[tuple[str, str]]:
    from Bio import SeqIO
    records = []
    for rec in SeqIO.parse(path, "fasta"):
        records.append((rec.id, str(rec.seq)))
    if not records:
        raise ValueError(f"No FASTA records in {path}")
    return records


def safe_write_csv(df, path: str, force: bool = False) -> None:
    if os.path.exists(path) and not force:
        raise FileExistsError(f"{path} exists; pass --force to overwrite. Prior results are preserved by default.")
    tmp = f"{path}.tmp.{os.getpid()}"
    df.to_csv(tmp, index=False)
    os.replace(tmp, path)


def load_candidates(path: str):
    import pandas as pd
    df = pd.read_csv(path)
    required = {"candidate_id", "source_protein", "epitope_class", "peptide_sequence",
                "start_aa", "end_aa", "peptide_length", "reference_status"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path} is missing required columns: {sorted(missing)}")
    return df


def now_stamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%d")
