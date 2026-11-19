"""Tests for the sequence-validation stage (src/python/validate_sequences.py)
and the metadata-integrity check (src/python/verify_target_metadata.py).

All sequences here are ARTIFICIAL test strings, not real biological sequences,
per the project blueprint.
"""

import hashlib

import pandas as pd
import pytest

from validate_sequences import validate_fasta
from verify_target_metadata import check_metadata

# Artificial 20-residue test sequence (standard AA alphabet, not a real protein).
ARTIFICIAL_SEQ = "ACDEFGHIKLMNPQRSTVWY"


def write_fasta(tmp_path, content: str, name: str = "test.fasta"):
    path = tmp_path / name
    path.write_text(content, encoding="utf-8")
    return path


def test_valid_fasta(tmp_path):
    fasta = write_fasta(tmp_path, f">artificial_1 some description\n{ARTIFICIAL_SEQ}\n")
    report = validate_fasta(fasta)

    assert len(report) == 1
    row = report.iloc[0]
    assert row["record_id"] == "artificial_1"
    assert row["length"] == len(ARTIFICIAL_SEQ)
    assert row["is_valid"] == True  # noqa: E712 - pandas bool, keep explicit
    assert row["invalid_symbols"] == ""
    assert row["sha256"] == hashlib.sha256(ARTIFICIAL_SEQ.encode("utf-8")).hexdigest()


def test_empty_fasta_raises(tmp_path):
    fasta = write_fasta(tmp_path, "")
    with pytest.raises(ValueError, match="No FASTA records"):
        validate_fasta(fasta)


def test_header_only_record_is_invalid(tmp_path):
    fasta = write_fasta(tmp_path, ">empty_record\n")
    report = validate_fasta(fasta)
    assert len(report) == 1
    assert report.iloc[0]["length"] == 0
    assert report.iloc[0]["is_valid"] == False  # noqa: E712


def test_invalid_aa_symbols(tmp_path):
    fasta = write_fasta(tmp_path, ">bad_seq\nACDEFJ*X\n")  # J, *, X are not standard AAs
    report = validate_fasta(fasta)
    row = report.iloc[0]
    assert row["is_valid"] == False  # noqa: E712
    assert set(row["invalid_symbols"].split(",")) == {"*", "J", "X"}
    # Record is reported, not silently dropped or rewritten.
    assert row["length"] == 8


def test_multiple_records(tmp_path):
    content = (
        ">rec_a\nACDEFGHIK\n"
        ">rec_b\nLMNPQRSTV\n"
        ">rec_c\nWYACDEFGH\n"
    )
    fasta = write_fasta(tmp_path, content)
    report = validate_fasta(fasta)
    assert list(report["record_id"]) == ["rec_a", "rec_b", "rec_c"]
    assert list(report["length"]) == [9, 9, 9]
    assert report["is_valid"].all()


def test_duplicate_sequence_flagged(tmp_path):
    content = (
        ">dup_1\nACDEFGHIK\n"
        ">dup_2\nACDEFGHIK\n"
        ">dup_3\nACDEFGHIK\n"
    )
    fasta = write_fasta(tmp_path, content)
    report = validate_fasta(fasta)
    assert list(report["duplicate_count"]) == [0, 1, 2]


def test_missing_metadata_flagged():
    # A target-metadata row missing required fields must fail verification.
    df = pd.DataFrame(
        [
            {
                "target_id": "ARTIFICIAL_T1",
                "organism": "Cladosporium sphaerospermum",
                "strain": "",                      # missing/empty strain
                "accession": "",                   # missing accession
                "accession_version": "ART.1",
                "sequence_length": 123,
                "completeness": "partial",
            }
        ]
    )
    report = check_metadata(df)
    row = report.iloc[0]
    assert row["is_verified"] == False  # noqa: E712
    missing = set(row["missing_columns"].split(","))
    assert {"strain", "accession"} <= missing


def test_wrong_length_is_reported_not_hidden(tmp_path):
    # validate_fasta does not judge against an expected length; it reports the
    # true length so a mismatch against metadata is visible downstream.
    fasta = write_fasta(tmp_path, f">artificial_wrong_len\n{ARTIFICIAL_SEQ}A\n")
    report = validate_fasta(fasta)
    row = report.iloc[0]
    assert row["length"] == 21
    assert row["length"] != 20  # mismatch vs. an arbitrary expected 20-aa length
    assert row["is_valid"] == True  # noqa: E712 - symbols valid; length is surfaced, not enforced here


def test_malformed_retrieval_response_rejected(tmp_path):
    # An HTTP error body / XML error page saved where a FASTA was expected
    # (e.g. a failed NCBI efetch) contains no '>' records.
    fasta = write_fasta(tmp_path, "HTTP/1.1 429 Too Many Requests\n<Error>API rate limit exceeded</Error>\n")
    with pytest.raises(ValueError):
        validate_fasta(fasta)
