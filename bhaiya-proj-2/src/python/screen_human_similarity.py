"""Whole-proteome similarity screen of target sequences vs the human proteome.

Method (documented in results/logs/pin2_seq_similarity.log):
  - Bio.Align.PairwiseAligner, mode='local' (Smith-Waterman-like), BLOSUM62,
    open_gap_score=-11, extend_gap_score=-1.
  - EVERY human reviewed (Swiss-Prot) sequence from data/reference/human_proteome.fasta
    (UniProtKB 2026_04, 20429 records) is aligned against each query.
  - 25 human records contain selenocysteine ('U'), which is absent from BLOSUM62's
    alphabet; for alignment scoring only, 'U' is substituted in-memory with 'X'
    (files on disk are untouched). Documented in the log.
  - Metrics per alignment (best local alignment returned by the aligner):
      alignment_length      = number of alignment columns (including gap columns)
      identity_pct          = identical-residue columns / alignment_length * 100  (BLAST-style)
      target_coverage_pct   = aligned target residues (non-gap) / target length * 100
      human_coverage_pct    = aligned human residues (non-gap) / human subject length * 100
      bitscore_proxy        = raw BLOSUM62 local alignment score. NOT a BLAST bitscore;
                              a proxy for ranking only.
  - Retention filter: identity_pct >= 30.0 AND target_coverage_pct >= 40.0
  - Output: results/tables/human_protein_similarity.csv (all retained hits)
            results/tables/human_protein_similarity_top50.csv (top 50 by bitscore_proxy per query)
"""

import csv
import sys
import time
from pathlib import Path

import pandas as pd
from Bio import SeqIO, Align
from Bio.Align import substitution_matrices

TARGETS = [('PKS1', 'PKS1'), ('SOD1REF', 'SOD1REF')]
HUMAN = 'data/reference/human_proteome.fasta'
OUT_ALL = 'results/tables/human_protein_similarity.csv'
OUT_TOP50 = 'results/tables/human_protein_similarity_top50.csv'

IDENTITY_MIN = 30.0
COVERAGE_MIN = 40.0


def metrics(aln, target_len, human_len):
    a = aln[0]
    b = aln[1]
    cols = len(a)
    matches = sum(1 for x, y in zip(a, b) if x == y and x != '-')
    t_non_gap = cols - a.count('-')
    h_non_gap = cols - b.count('-')
    return {
        'identity_pct': round(100.0 * matches / cols, 4) if cols else 0.0,
        'target_coverage_pct': round(100.0 * t_non_gap / target_len, 4),
        'human_coverage_pct': round(100.0 * h_non_gap / human_len, 4),
        'alignment_length': cols,
        'bitscore_proxy': round(float(aln.score), 4),
    }


def subj_id(record):
    # e.g. "sp|P69905|HBA_HUMAN ..." -> "P69905"
    parts = record.id.split('|')
    return parts[1] if len(parts) >= 2 else record.id


def alignable(seq):
    # BLOSUM62 has no 'U' (selenocysteine); substitute in-memory with 'X' for scoring.
    return seq.upper().replace('U', 'X')


def main():
    t_start = time.time()
    matrix = substitution_matrices.load('BLOSUM62')
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.substitution_matrix = matrix
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1

    targets = {}
    for r in SeqIO.parse('data/processed/validated_targets.fasta', 'fasta'):
        targets[r.id] = str(r.seq)

    human_records = list(SeqIO.parse(HUMAN, 'fasta'))
    print(f'human records: {len(human_records)}', flush=True)

    all_rows = []
    scored = {q: 0 for q, _ in TARGETS}
    for qid, tname in TARGETS:
        tseq = targets[tname]
        for rec in human_records:
            hseq = str(rec.seq)
            aln = aligner.align(alignable(tseq), alignable(hseq))[0]
            m = metrics(aln, len(tseq), len(hseq))
            scored[qid] += 1
            if m['identity_pct'] >= IDENTITY_MIN and m['target_coverage_pct'] >= COVERAGE_MIN:
                all_rows.append({'query_id': qid, 'human_subject_id': subj_id(rec),
                                 'target_name': tname, **m})
        print(f'query {qid}: scored {scored[qid]} human subjects, '
              f'elapsed {time.time()-t_start:.1f}s', flush=True)

    df = pd.DataFrame(all_rows, columns=['query_id', 'human_subject_id', 'target_name',
                                         'identity_pct', 'target_coverage_pct',
                                         'human_coverage_pct', 'alignment_length',
                                         'bitscore_proxy'])
    df = df.sort_values(['query_id', 'bitscore_proxy'], ascending=[True, False]).reset_index(drop=True)
    df.to_csv(OUT_ALL, index=False)

    top50 = df.groupby('query_id', group_keys=False).head(50).reset_index(drop=True)
    top50.to_csv(OUT_TOP50, index=False)

    print(f'retained hits: {len(df)} -> {OUT_ALL}')
    print(f'top50 rows: {len(top50)} -> {OUT_TOP50}')
    print(df.groupby('query_id').size())
    print(f'total elapsed {time.time()-t_start:.1f}s')
    return 0


if __name__ == '__main__':
    sys.exit(main())
