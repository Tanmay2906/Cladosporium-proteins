# pin2_seq — Final report (Cladosporium immunoinformatics)

Date: 2026-11-14 | Worker: pin2_seq (pane w5:p6) | CWD: /Users/aaashish/Projects/bhaiya-proj-2
Scope honoured: no git, no IEDB remote calls, no network except UniProtKB REST (task 4). No sequences fabricated; raw files under data/raw/ untouched.
SOD1 labelling: every SOD1-derived row is `REFERENCE_ONLY_NOT_C_sphaerospermum` (UniProt A0AB34KRD1, *C. halotolerans*, not a *C. sphaerospermum* target).

## Stage counts

1. **Sequence validation** — 2 records processed, 2 valid.
   - PKS1 (QBZ73626.1): 235 aa, sha256 f70f60d0cbb25ccc6980a2689271fc88f532d9705dfaf2b490d47014ce42f43c, valid.
   - SOD1REF (A0AB34KRD1): 154 aa, sha256 ca2ff2b0dba341bac4ee0ef9dd26f85998b61ea28dc21f4fe4f9e177ae1dc97b, valid.
   - **Failure recorded:** `src/python/validate_sequences.py` on the raw `data/raw/SOD1.fasta` exits 1 — Biopython's `fasta` parser rejects the leading `#` comment lines in that file. Raw file unchanged. The sequence was validated from a comment-stripped derivative (`data/processed/SOD1_comments_stripped.fasta`; only the leading `#` lines removed, header and sequence byte-identical). The two-run requirement was therefore executed as: PKS1 on raw input, SOD1REF on the documented derivative.
2. **Protein characteristics** — 2 records.
   - PKS1: 24,974.22 Da, pI 5.69, instability 39.49 (unstable by the >40 rule? No: 39.49 < 40, borderline), aliphatic index 92.21, GRAVY +0.0072.
   - SOD1REF: 15,771.04 Da, pI 6.21, instability 32.74, aliphatic index 61.43, GRAVY −0.553.
   - **Failure recorded:** the script is incompatible with installed Biopython 1.88 (`ProteinAnalysis.aliphatic_index` does not exist; `get_amino_acids_percent` was renamed to cached_property `amino_acids_percent`). Script NOT modified; a logged shim injects (a) aliphatic index per Ikai 1980 (AI = Ala% + 2.9·Val% + 3.9·(Ile%+Leu%), mole %) and (b) `get_amino_acids_percent = amino_acids_percent` alias. First shim run still failed (second rename); resolved in second attempt. No values other than these were substituted.
   - amino_acid_composition.csv: 40 rows (2 records × 20 amino acids), exact counts + percents, with REFERENCE_ONLY flag on SOD1REF rows.
3. **Peptide candidates** — 1,849 total.
   - MHC-I: 1,488 (PKS1 906 = 228+227+226+225 for 8/9/10/11-mers; SOD1REF 582 = 147/146/145/144).
   - MHC-II 15-mers: 361 (PKS1 221; SOD1REF 140).
   - IDs: `PKS1_MHCI_000001..`, `PKS1_MHCII_000001..`, `SOD1REF_MHCI_000907..`, `SOD1REF_MHCII_000222..` (SOD1REF numbering continues the per-kind counter after PKS1's; prefixes are as required). One-based inclusive coordinates, `peptide_length` retained.
4. **Human proteome download** — SUCCESS.
   - UniProtKB/Swiss-Prot release **2026_04** (release date 2026-11-04; API deployment 2026-11-10), retrieved 2026-11-14, HTTP 200.
   - URL: https://rest.uniprot.org/uniprotkb/stream?query=organism_id:9606+AND+reviewed:true&format=fasta
   - 20,429 canonical reviewed sequences (NOT "~54k"), 13,670,961 bytes, uncompressed sha256 `13852e174364bf1c705206c38eacae57210f5bca9955f05be6b4ff48f727d37b`, md5 `3064f6cbf0e881570a2663b5a838f2ff`, gzip checksum (`gzip -n -c | shasum -a 256`, deterministic) `e80a300628162d7296a594afc95c32541be8c9e1ba788c1229ccbda17746a433` (7,507,552 gzip bytes). Recorded in `data/reference/human_proteome_manifest.csv`.
   - The stream returned plain FASTA (no pre-gzipped file exists from this endpoint), hence the gzip checksum is computed as documented in the manifest.
   - Reported record count discrepancy noted: task referenced "54k human proteins"; reviewed Swiss-Prot 2026_04 contains 20,429. Percentages downstream use 20,429 as denominator (column name `match_percent_of_54k_human_proteins` retained for spec compliance; denominator documented in file header).
5. **Whole-proteome similarity screen** — complete, all 20,429 human records scored against both queries (40,858 local alignments), none skipped.
   - Method: Bio.Align.PairwiseAligner, mode='local', BLOSUM62, open_gap_score=−11, extend_gap_score=−1; identity_pct = identical columns / alignment columns (incl. gap columns); coverage = aligned residues / full length of that sequence; bitscore_proxy = raw BLOSUM62 local score (NOT a BLAST bitscore). 25 human records contain selenocysteine `U` (absent from BLOSUM62) → substituted in-memory with `X` for scoring only; no files modified.
   - Retained (identity ≥ 30% AND target coverage ≥ 40%): **73 hits — PKS1 9, SOD1REF 64.**
   - Notable: SOD1REF aligns to human superoxide dismutase [Cu-Zn] P00441 at 55.3% identity, 98.7% target coverage (score proxy 435) — confirms A0AB34KRD1 is close to the human enzyme; a further reason the SOD1REF branch is reference-only. PKS1 retained hits are all borderline (~30–32% identity, 40–51% coverage) — no strong human PKS-class counterpart.
   - top50 file holds 59 rows (PKS1 contributes only 9). Runtime ≈ 55 s; full proteome covered, no truncation.
   - Helper script added: `src/python/screen_human_similarity.py` (new file, unrelated scripts untouched).
6. **Peptide-level exact-match screen** — all 1,849 candidates screened; exact substring only (concatenated proteome, `#` separator, `str.count`/`str.find`). Approximate matches out of scope (stated in file header comments).
   - 25 peptides have ≥1 exact match in the human proteome: **1 PKS1** (`PKS1_MHCI_000074` AVVLGNGL → Q8NH49) and **24 SOD1REF** (MHC-I only; 23 of the 24 match human SOD1 P00441; `SOD1REF_MHCI_000922` GTVTFEQA → Q9P2K5). Zero MHC-I 10/11-mer PKS1 other hits; **0 MHC-II 15-mer exact matches**; 0 matches for 9/11-mer PKS1 peptides beyond the one 8-mer.
   - Percent denominator = 20,429 (see note 4); each match = 0.004895%.
7. **This report.**

## Absolute paths of outputs

- /Users/aaashish/Projects/bhaiya-proj-2/results/tables/sequence_validation_report.csv
- /Users/aaashish/Projects/bhaiya-proj-2/data/processed/validated_targets.fasta (PKS1, SOD1REF; case-normalised to upper; descriptions carry REFERENCE_ONLY flag)
- /Users/aaashish/Projects/bhaiya-proj-2/data/processed/SOD1_comments_stripped.fasta (task-1 derivative)
- /Users/aaashish/Projects/bhaiya-proj-2/results/tables/target_protein_characteristics.csv
- /Users/aaashish/Projects/bhaiya-proj-2/results/tables/amino_acid_composition.csv
- /Users/aaashish/Projects/bhaiya-proj-2/data/processed/peptide_candidates.csv
- /Users/aaashish/Projects/bhaiya-proj-2/data/reference/human_proteome.fasta
- /Users/aaashish/Projects/bhaiya-proj-2/data/reference/human_proteome_manifest.csv
- /Users/aaashish/Projects/bhaiya-proj-2/results/tables/human_protein_similarity.csv (73 retained hits)
- /Users/aaashish/Projects/bhaiya-proj-2/results/tables/human_protein_similarity_top50.csv (59 rows; ≤50 per query)
- /Users/aaashish/Projects/bhaiya-proj-2/results/tables/candidate_human_similarity.csv (1,849 rows + comment header)
- /Users/aaashish/Projects/bhaiya-proj-2/src/python/screen_human_similarity.py (new helper)
- Logs: /Users/aaashish/Projects/bhaiya-proj-2/results/logs/pin2_seq_{validate,characterize,peptides,humanproteome,similarity,pepsim}.log

## Failures / caveats (summary)

1. Stock validator cannot read raw SOD1.fasta (leading `#` comments) → validated via documented comment-stripped derivative; raw untouched. (Task 1)
2. `characterize_proteins.py` is broken under Biopython 1.88 → logged-shim run, script unmodified; aliphatic index per Ikai 1980. (Task 2)
3. "54k human proteins" denominator does not exist in reviewed Swiss-Prot 2026_04 (20,429 canonical sequences) → documented; percentages computed against 20,429. (Tasks 4 and 6)
4. 25 human sequences contain `U`; aligned with in-memory `X` substitution for scoring only. (Task 5)
5. bitscore_proxy is a raw BLOSUM62 local score, not a BLAST bitscore; no E-values were computed. (Task 5)
6. `config/human_proteome.txt` still says "NOT YET RETRIEVED" — left unmodified per minimal-touch; the manifest CSV (data/reference/human_proteome_manifest.csv) is the authoritative provenance record.
7. No git operations, no IEDB/NCBI calls; only UniProtKB REST used.
8. All SOD1REF outputs are REFERENCE_ONLY_NOT_C_sphaerospermum; PKS1 is a partial (235 aa), unpublished direct submission (QBZ73626.1) — downstream claims must stay hypothesis-generating.
