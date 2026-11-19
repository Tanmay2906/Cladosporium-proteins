"""Tests for the peptide-generation stage (src/python/generate_peptides.py).

All sequences are ARTIFICIAL test strings, per the project blueprint. The
20-residue ARTIFICIAL_SEQ below is the standard AA alphabet as a literal test
string; it is not a real biological sequence.
"""

import pytest

from generate_peptides import generate_peptides

# Artificial 20-residue test sequence.
ARTIFICIAL_SEQ = "ACDEFGHIKLMNPQRSTVWY"


def test_artificial_sequence_is_20_residues():
    assert len(ARTIFICIAL_SEQ) == 20


def test_window_counts_mhc_lengths():
    expected = {8: 13, 9: 12, 10: 11, 11: 10}  # 20 - L + 1
    for length, count in expected.items():
        peptides = generate_peptides(ARTIFICIAL_SEQ, length, "ARTIFICIAL", "MHCI")
        assert len(peptides) == count, f"length {length}"


def test_window_counts_mhc_ii_15mer():
    peptides = generate_peptides(ARTIFICIAL_SEQ, 15, "ARTIFICIAL", "MHCII")
    assert len(peptides) == 6  # 20 - 15 + 1


def test_window_counts_total_overlapping():
    total = 0
    for length in (8, 9, 10, 11):
        total += len(generate_peptides(ARTIFICIAL_SEQ, length, "ARTIFICIAL", "MHCI"))
    assert total == 46


def test_coordinate_correctness_one_based_inclusive():
    for length in (8, 9, 10, 11, 15):
        peptides = generate_peptides(ARTIFICIAL_SEQ, length, "ARTIFICIAL", "MHCI")
        for i, pep in enumerate(peptides):
            assert pep["start_aa"] == i + 1
            assert pep["end_aa"] == i + length
            assert pep["peptide_length"] == length
            # One-based inclusive slice of the source sequence.
            assert pep["peptide_sequence"] == ARTIFICIAL_SEQ[pep["start_aa"] - 1 : pep["end_aa"]]


def test_coordinate_bounds():
    peptides = generate_peptides(ARTIFICIAL_SEQ, 8, "ARTIFICIAL", "MHCI")
    assert peptides[0]["start_aa"] == 1
    assert peptides[0]["end_aa"] == 8
    assert peptides[-1]["start_aa"] == 13
    assert peptides[-1]["end_aa"] == 20


def test_candidate_ids_unique_and_sequential():
    counter_start = 100
    peptides = generate_peptides(ARTIFICIAL_SEQ, 8, "ARTIFICIAL", "MHCI",
                                 start_index=counter_start)
    ids = [p["candidate_id"] for p in peptides]
    assert len(set(ids)) == len(ids)
    assert ids[0] == f"ARTIFICIAL_MHCI_{counter_start:06d}"
    assert ids[-1] == f"ARTIFICIAL_MHCI_{counter_start + len(ids) - 1:06d}"


def test_ids_unique_across_lengths_with_running_counter():
    counter = 1
    all_ids = []
    for length in (8, 9, 10, 11):
        peptides = generate_peptides(ARTIFICIAL_SEQ, length, "ARTIFICIAL", "MHCI",
                                     start_index=counter)
        all_ids.extend(p["candidate_id"] for p in peptides)
        counter += len(peptides)
    assert len(set(all_ids)) == len(all_ids) == 46


def test_invalid_target_kind_rejected():
    with pytest.raises(ValueError, match="target_kind"):
        generate_peptides(ARTIFICIAL_SEQ, 9, "ARTIFICIAL", "MHCIII")


def test_peptide_longer_than_sequence_yields_no_windows():
    assert generate_peptides("ACDEFG", 20, "ARTIFICIAL", "MHCII") == []
