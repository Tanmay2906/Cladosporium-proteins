# Methodological Decisions — *Cladosporium sphaerospermum* Immunoinformatics Blueprint

Blueprint Part 17 deliverable. Status: **decision register + rationale**; numerical results are marked PENDING_RESULTS until the corresponding pipeline steps run (as of 2026-11-14: sequences, peptides, human-proteome download, similarity screens, and **IEDB MHC-I/MHC-II prediction runs, B-cell (D-07), conservation (D-08), prioritisation (D-09) and sensitivity** were executed 2026-11-14; population coverage (D-13) remains not_evaluated by design).

Scope note (honesty): the blueprint document itself is **not stored in this repository**. References to blueprint Parts 1, 10 and 17 (and the review policy categories) derive from the orchestrator brief (pane w5:p6) and from project documents on disk (`target_selection_rationale.md`, `target_metadata.csv`, `data/manifests/retrieval_manifest.csv`, `data/processed/target_evidence_matrix.csv`, `config/`, `src/python/`). Where a blueprint Part number is cited without an on-disk counterpart, that is stated. Verified facts (accessions, counts, checksums) were re-checked on 2026-11-13/14; anything not yet executed is marked **PENDING_RESULTS** or **PENDING**.

How to read each decision: Decision / Rationale / Alternatives considered / Revisit condition / Single source of truth (the on-disk path where the decision is revisited and updated).

---

## 0. Decision register

| ID | Decision (one line) | Source of truth (revisit here) |
|---|---|---|
| D-01 | Target identity fixed per verified evidence; references are reference-only | `target_selection_rationale.md`; `data/processed/target_evidence_matrix.csv` |
| D-02 | Ambiguous residues retained with flags, never deleted | `src/python/validate_sequences.py`; `docs/methodology.md` (this file) |
| D-03 | Human reference = UniProtKB/Swiss-Prot UP000005640, reviewed only (retrieved: 2026_04, 2026-11-14) | `config/human_proteome.txt`; R021; `data/reference/human_proteome_manifest.csv` |
| D-04 | Similarity search: documented local method (+optional NCBI BLAST REST cross-check); prespecified review policy | `data/manifests/retrieval_manifest.csv` (R020); `docs/data_provenance.md` |
| D-05 | MHC-I: 5 alleles, netmhcpan_ba 4.1, lengths 8–11 — executed 2026-11-14 | `config/hla_panel.csv`; `config/prediction_methods.yaml`; R022; `results/tables/mhci_all_predictions.csv` |
| D-06 | MHC-II: 5 alleles, netmhciipan_ba 4.1, 15-mer flanks, preserve core — executed 2026-11-14 | `config/hla_panel.csv`; `config/prediction_methods.yaml`; R023; `results/tables/mhcii_all_predictions.csv` |
| D-07 | B-cell: BepiPred 2.0 via IEDB antibody-epitope prediction | `config/prediction_methods.yaml`; `docs/data_provenance.md` |
| D-08 | Conservation = exact-sequence fraction; within- vs across-species separate | `docs/methodology.md` D-08; `results/tables/epitope_conservation.csv` (DONE 2026-11-14) |
| D-09 | Prioritisation by evidence categories, not weighted sum; exclusion log mandatory; sensitivity scenarios | `docs/methodology.md` D-09; `results/tables/candidate_master_table.csv` (DONE 2026-11-14) |
| D-10 | Structural analysis only if core complete; not scheduled yet | `docs/methodology.md` D-10 |
| D-11 | Source-of-truth path register for every decision | Section D-11 below |
| D-12 | IEDB allele naming / DQ α/β cis-chain-pair mapping / byte-identical raw-response retention | `docs/methodology.md` D-12; `config/prediction_methods.yaml` (`iedb_allele_mapping`); R023 |

---

## D-01 Target identity

**Decision.**
- **PKS1-class anchor (include, with caveat):** `QBZ73626.1` — "type III polyketide synthase, **partial**", isolate **DWEL1**, 235 aa, `no-ends` (5′ ends missing), codon 1..235 (`cond_enzymes`, InterPro cl09938 / CDD:447866), GenBank direct submission **MK213995.1** (Sanger dideoxy; Manoharan et al., *unpublished*, 22-Nov-2018), taxid 92950, GI 1613841415. Sequence saved byte-exact at `data/raw/PKS1.fasta` (sha256 `6f2ef2c1cb5f9e088d2de68fe6067b855279297c052fc93bd6cc2eaa59cb1b4e`). **This is a type III PKS; it is NOT the canonical iterative type I "PKS1".** The name "PKS1" is project nomenclature, not annotation.
- **SOD1 in *C. sphaerospermum*: NOT VERIFIED — excluded.** UniProt search (R001/R002) returns no *sphaerospermum* SOD1; NCBI esearch with organism filter returns **count 0** for SOD1/superoxide/dismutase (R010), catalase (R011), glutathione (R012), melanin (R013), and `PKS1[All Fields] OR SOD1[All Fields]` (R005). A local BLOSUM62 screen of the *C. sphaerospermum* NCBI proteome (324 records) vs a Cu-Zn SOD1 reference returned **0 hits above threshold** (R020; `results/tables/putative_sod_homology_screen.csv`).
- **Reference-only (NOT targets):** `A0AB34KRD1` (*C. halotolerans* TM138-S3, SOD1 Cu-Zn, 154 aa, EC 1.15.1.1, complete; `data/raw/SOD1.fasta`) and `A0AB34KVX2` (*C. halotolerans* TM138-S3, non-reducing iterative type I PKS1, 2175 aa, WGS-derived, CAUTION "WGS preliminary data"). Both are genus-level references for domain architecture / comparative epitope context. Every downstream result derived from them must be labelled reference-based.
- **Six unannotated genome assemblies (brief-supplied; not yet integrated, not a basis for claims):** `GCA_001931905.1` (strain **IMV00045**, Chernobyl Unit-4 turbine hall wall, collection-date **1998-06-16**); `GCA_000261425.2` (UM 843); `GCA_019191205.1`; `GCA_023621355.1`; `GCA_052323905.1`; `GCA_000261425.1`. Note: `GCA_000261425.1`/`.2` are versioned releases of one assembly (UM 843) — deduplicate at analysis time. These carry **no retrievable accession-level annotations locally**; they are listed as *candidate in-silico mining sources*, not as targets.
- **Never** synthesise a consensus or hybrid PKS or SOD sequence, and never back-substitute the reference into the target set (blueprint Part 1 guidance: do not manufacture sequences; document the limitation and reassess — see `target_selection_rationale.md` §2).

**Rationale.** The only verified *C. sphaerospermum* PKS-class record is QBZ73626.1; treating the *halotolerans* references as targets would fabricate species-level identity. The six assemblies are unannotated; annotating them is real retrieval+analysis work, not a documentation assumption. Domain evidence (InterPro/CDD/PROSITE/Pfam, R015) supports class calls for the three raw sequences only.

**Alternatives considered.** (a) Use A0AB34KVX2 as the PKS1 target — rejected: species mismatch; explicitly reference-only (blueprint: reference ≠ target). (b) Back-fill SOD1 from the genus reference — rejected: same; count 0 verified locally. (c) Drop the PKS anchor — rejected: QBZ73626.1 is genuine (NCBI-verified, though unpublished); caveat labelling is the correct handling. (d) Wait for genome mining before any downstream work — rejected as a blocker; adopted as a parallel PENDING workstream (see D-01 revisit). (e) Download the unquoted 641-record NCBI set — rejected: ~317 records are other organisms (see `retrieval_manifest.csv` R017 premise note); the strict taxid-92950 set (324 records) was downloaded instead.

**Revisit condition.** Any of: complete *C. sphaerospermum* DWEL1 or equivalent proteome/assembly with an iterative type I PKS1 or Cu-Zn SOD1; annotation of any of the six GCAs revealing either target class; new UniProt/NCBI releases; FungiDB access becoming available (R016 limitation cleared); any retraction/update of MK213995.1.

**Single source of truth:** `target_selection_rationale.md` (reassess and rewrite in place) cross-checked against `data/processed/target_evidence_matrix.csv`, `target_metadata.csv`, `data/manifests/retrieval_manifest.csv`.

---

## D-02 Ambiguous residues

**Decision.** Ambiguous / non-standard residues (X, B, Z, U, J, O, lowercase artefacts, etc.) are **retained with flags**, never deleted, replaced, or trimmed. `src/python/validate_sequences.py` reports `invalid_symbols`, sets `is_valid=False`, and keeps the record and its `sha256`; downstream steps may segregate flagged records but may not edit them. Verified locally 2026-11-14: both `data/raw/PKS1.fasta` (235 aa) and `data/raw/SOD1.fasta` (154 aa) contain **only the 20 standard residues** (nonstandard = {}), so no current flag exists — the policy still binds all future records (genome-mined fragments, new accessions).

**Rationale.** Deleting residues shifts coordinates, fabricates a database record that does not exist, and silently breaks alignment/conservation bookkeeping (blueprint Part 17: sequence-integrity principle). Retaining with flags preserves provenance and keeps every exclusion auditable.

**Alternatives considered.** (a) Drop invalid records — rejected: lossy and unauditable. (b) Replace ambiguous residues with Ala/gap — rejected: fabrication; corrupts coordinates and conservation denominators. (c) Mask silently in downstream tools — rejected: flags must propagate to tables (`invalid_symbols`, `is_valid`). (d) BLAST/masked-residue tolerant mode — noted as acceptable *annotation* of a flagged record, never as silent correction.

**Revisit condition.** If any future record actually carries an ambiguous residue, document the flag in the validation report and carry `is_valid` into the candidate tables; policy itself changes only via blueprint amendment.

**Single source of truth:** `src/python/validate_sequences.py` (implementation) and `docs/methodology.md` (this policy text, this file).

---

## D-03 Human reference proteome

**Decision.** All human-similarity screening uses **UniProtKB/Swiss-Prot reviewed entries only**, reference proteome **UP000005640** (*Homo sapiens*). **No TrEMBL/predicted human entries.** The FASTA **was retrieved 2026-11-14** via the UniProtKB REST stream (`https://rest.uniprot.org/uniprotkb/stream?query=organism_id:9606+AND+reviewed:true&format=fasta`): release **2026_04**, **20,429** canonical reviewed records, 13,670,961 bytes, sha256 `13852e174364bf1c705206c38eacae57210f5bca9955f05be6b4ff48f727d37b` (local copy `data/reference/human_proteome.fasta`; manifest `data/reference/human_proteome_manifest.csv` entry HP001; ledger R021). 25 records contain U and are handled as X. The reviewed denominator is 20,429 (NOT "~54k"). Screening (D-04) executed 2026-11-14: whole-protein + peptide exact-match screens (`results/tables/human_protein_similarity.csv`, `candidate_human_similarity.csv`).

**Rationale.** Reviewed Swiss-Prot entries are the lowest-error human background; including predicted entries would inflate "human similarity" noise and could falsely exclude fungal peptides (blueprint: a similarity screen against a noisy reference is not a safety screen).

**Alternatives considered.** (a) Full UniProtKB human (Swiss-Prot + TrEMBL) — rejected: unreviewed entries; false-positive off-target flags. (b) RefSeq human proteome (GCF_000001405) — viable cross-check, secondary only, to be added if the screen is borderline; not the primary reference. (c) NCBI nr / human subset via BLAST REST — acceptable as the D-04 cross-check, but reference set identity must remain Swiss-Prot-reviewed for the denominator. (d) TrEMBL-inclusive download — rejected (same as (a)).

**Revisit condition.** If UniProt changes UP000005640 composition at a release boundary, re-screen with a new HP manifest row; keep release/size/checksum pinned per download.

**Single source of truth:** `data/reference/human_proteome_manifest.csv` (HP001; never stores sequences) + `data/manifests/retrieval_manifest.csv` (R021) + `docs/data_provenance.md` §A #17.

---

## D-04 Similarity search & review policy

**Decision.** Two complementary layers, both prespecified:
1. **Documented local method (primary):** local alignment of the reference (A0AB34KRD1) against all 324 records of `data/reference/c_sphaerospermum_ncbi_proteome.fasta` (sha256 `1aaf2b5c753d29c6fec6bc432905dc32a2c37115fc3307b5d11b50c8d2c67222`) using `Bio.Align.PairwiseAligner` (Biopython 1.88, local mode, BLOSUM62, gap open −10, gap extend −0.5), thresholds **identity ≥ 30% AND reference coverage ≥ 40%**. Implemented and run 2026-11-13; result: **0 hits above threshold**; 30 records reach coverage ≥ 40% (best 64.29%, YFQ81086.1 / YFQ81084.1, identity 23.31%); 42 records coverage ≥ 20%; output `results/tables/putative_sod_homology_screen.csv` (explicit "no hits above threshold" row). No record may be labelled SOD-like.
2. **NCBI BLAST REST (cross-check, PENDING):** `blastp` via `https://blast.ncbi.nlm.nih.gov/Blast.cgi` against taxid 92950 (*C. sphaerospermum*) and against UniProtKB/Swiss-Prot for the human screen (D-03). RID/RTOE, parameters, database, and date recorded at run time; results appended to `docs/data_provenance.md` and the manifest. **Not run yet.**

**Review policy categories (from blueprint Part 4/17 brief; applied at screening time, thresholded, never hand-waved):**
- `no_flag_under_defined_screen` — below prespecified thresholds under both layers; logged; no exclusion (exclusion would be unjustified).
- `review_required` — meets ≥1 prespecified partial criterion (e.g., identity ≥ 30% OR coverage ≥ 40% but not both; or within 5 percentage points of a threshold; or any human-screen hit at any identity) → routed to human review **before** any inclusion/exclusion decision; outcome recorded in the exclusion log.
- `excluded_under_prespecified_rule` — fails the prespecified rule with evidence logged (e.g., human-screen exact match ≥ 8–11 aa MHC-I / ≥ 15-mer MHC-II core at ≥ the prespecified identity; or non-target contamination) → excluded with reason, date, and rule ID.

**Rationale.** The local method is cheap, fully reproducible, and already executed; BLAST REST is the community-standard cross-check against a much larger corpus (catches what 324 short Sanger fragments cannot: genomes, longer ORFs). Prespecified categories prevent post-hoc threshold-picking (blueprint: thresholds arbitrary but prespecified — see D-09 and `docs/limitations.md`).

**Alternatives considered.** (a) DIAMOND — faster; unnecessary at 324 records; retained option for scaled genome mining of the six GCAs. (b) HMMER/Pfam profile search — good for domain calling, not for peptide-exact conservation; reserved for annotation support. (c) PSI-BLAST — unnecessary sensitivity for this corpus. (d) BLAST REST as primary — rejected: NCBI-side parameters drift and are harder to pin for audit; kept as cross-check. (e) Local BLAST+ database build — equivalent to (d) offline; acceptable substitute if NCBI REST is unavailable; document either way.

**Revisit condition.** Human proteome downloaded (D-03); BLAST REST cross-check executed; six GCAs annotated and screened; threshold sensitivity results (D-09 scenarios) show rank instability near the 30/40 line → prespecify a review band, do not move thresholds post hoc.

**Single source of truth:** `data/manifests/retrieval_manifest.csv` (R017–R020) + `results/tables/putative_sod_homology_screen.csv`; cross-check results → `docs/data_provenance.md`; exclusion decisions → `results/tables/candidate_exclusion_log.csv` (DONE 2026-11-14).

---

## D-05 MHC-I prediction

**Decision.** Panel (from `config/hla_panel.csv`, AFND-frequency-rationale for European and South Asian populations, well-supported by IEDB):
- `HLA-A*02:01`, `HLA-B*07:02`, `HLA-C*07:01`, `HLA-A*01:01`, `HLA-B*08:01`.

Method: **netmhcpan_ba (NetMHCpan 4.1, BA mode)** via IEDB, **peptide lengths 8, 9, 10, 11** (`config/prediction_methods.yaml` → `mhci`). BA (binding affinity) is the primary screen. **NetMHCpan EL (elution) is a secondary sensitivity check only**, run after BA results are reviewed, never in place of BA for primary ranking. IEDB recommended methods (2023.09 onwards; BA + EL pair) are noted in the config comment. Candidate IDs per `src/python/generate_peptides.py` format (`<SOURCE>_MHCI_000001`), 1-based inclusive coordinates. **Executed 2026-11-14** — user-approved IEDB run; 2 targets × 5 alleles × 4 lengths = 40 API calls; raw responses preserved in `data/raw/iedb/mhci/`; joined table `results/tables/mhci_all_predictions.csv` (7,440 rows: 1,488 candidates × 5 alleles; all candidates reproduced exactly with matching coordinates and peptide identity).

**Rationale.** BA estimates intrinsic binding, the stablest quantity across IEDB releases; EL additionally weighs processing/presentation likelihood but is more model-dependent. The panel covers frequent European/South Asian alleles (per AFND rationale strings in `config/hla_panel.csv`); HLA-C*07:01 is included as one of the most frequent HLA-C alleles. 8–11-mers are the canonical MHC-I length band.

**Alternatives considered.** (a) EL as primary — rejected: config explicitly retains BA primary; EL reserved for sensitivity. (b) MHCflurry 2.0 — viable, but blueprint/config standard is IEDB NetMHCpan; using two primary tools would fork the audit trail. (c) Wider panel (A*03, B*35, …) — rejected for primary; reserved as D-09 sensitivity scenario "panel trim/extend". (d) Lengths 9 only — rejected; 8–11 retained per config. (e) Mixing DQ/DP or murine alleles — out of scope; human vaccine-facing goal.

**Revisit condition.** Any IEDB method/version change (config header warns it can change); results unstable under BA vs EL sensitivity (D-09); panel trim scenario materially changes candidate ranks → re-document in `config/hla_panel.csv` and this section.

**Single source of truth:** `config/hla_panel.csv` (panel + rationale) + `config/prediction_methods.yaml` (method/version) ; IEDB retrieval rows appended to `data/manifests/retrieval_manifest.csv` at run time (R022; supersedes "PENDING").

---

## D-06 MHC-II prediction

**Decision.** Panel (`config/hla_panel.csv`, same population rationale):
- `HLA-DRB1*01:01`, `HLA-DRB1*04:01`, `HLA-DRB1*07:01`, `HLA-DQB1*02:01`, `HLA-DQB1*03:01`.

Method: **netmhciipan_ba (NetMHCIIpan 4.1, BA mode)** via IEDB; 15-mer peptides in flanked context (default IEDB flank setting; as-submitted). Only length 15 was executed (matches the candidates table exactly; IEDB-supported range 11–30). **The predicted 9-mer binding core reported by IEDB/NetMHCIIpan is preserved verbatim in all outputs** (`core_peptide` column) — cores are never re-derived, re-windowed, or truncated for convenience. Candidate IDs per `generate_peptides.py` (`<SOURCE>_MHCII_000001`). **Executed 2026-11-14** — user-approved IEDB run; 2 targets × 5 alleles × 1 length = 10 API calls; raw responses in `data/raw/iedb/mhcii/`; joined table `results/tables/mhcii_all_predictions.csv` (1,805 rows: 361 candidates × 5 alleles; all candidates reproduced exactly with matching coordinates and peptide identity).

**Rationale.** Flanked 15-mers are the IEDB default and capture ~all DR/DQ binders' register contexts; preserving the predicted core keeps the register auditable against IEDB output and avoids off-by-frame errors downstream (conservation of the core is the biologically meaningful unit for D-08).

**Alternatives considered.** (a) 13-mer or 17-mer windows — retained as D-09 sensitivity options, not primary. (b) NetMHCIIpan EL primary — rejected (same as D-05). (c) Core-only re-prediction — rejected: loses register; corrupts conservation mapping. (d) DP alleles — out of scope for primary panel.

**Revisit condition.** IEDB method change; core often predicted at window edges (possible signal that flank length is wrong for short proteins like QBZ73626.1's partial sequence — note the 235 aa sequence with 15-mer windows still yields valid cores, but edge cores deserve `review_required` flags); sensitivity scenario "window length" executed.

**Single source of truth:** `config/hla_panel.csv` + `config/prediction_methods.yaml` (mhcii block) + manifest rows at run time (R023; supersedes "PENDING") + allele mapping register (D-12).

---

## D-07 B-cell epitope prediction

**Decision.** Method: **BepiPred 2.0 (random-forest), accessed through the IEDB antibody-epitope prediction tool** (antibody epitope prediction endpoint), linear (sequential) B-cell epitopes only. Exact tool name/version, IEDB endpoint URL, parameters, download/execution date, and checksums of the output are recorded in `docs/data_provenance.md` (R024; executed 2026-11-14). Scores are used as an evidence category (D-09), never as a standalone gate. Conformational epitopes (DiscoTope etc.) are out of scope until a target structure exists (D-10).

**Rationale.** BepiPred 2.0 is the blueprint-named method and is available via IEDB without local installation, keeping the audit trail inside one toolchain. Linear epitopes are well-defined without 3-D structure.

**Alternatives considered.** (a) ABCpred — older RNN method; not blueprint-named. (b) BepiPred 1.0 — superseded by 2.0. (c) DiscoTope / conformational methods — require structure; excluded until D-10 is unlocked. (d) ElliPro — same structural dependency. (e) Local install — unnecessary via IEDB; version pinning weaker, so run-time recording is mandatory.

**Revisit condition.** Tool/version recorded at run time drifts from this section (IEDB retirement/change) → update `docs/data_provenance.md` + this section; D-09 candidate table build shows B-cell score non-independence from hydrophilicity → re-examine before ranking.

**Single source of truth:** `docs/data_provenance.md` (run-time record) + `config/prediction_methods.yaml` (to be extended with a `bcell:` block by the executing pin at run time — NOT done yet) + `docs/methodology.md` D-07.

---

## D-08 Conservation

**Decision.** Conservation statistic (per candidate peptide): **exact-sequence conservation fraction** =
> eligible homologues containing the exact peptide substring (case-insensitive on the 20-letter alphabet) ÷ total eligible homologues.

Reported **separately within-species** (*C. sphaerospermum* homologue set) and **across-species** (genus-level Cladosporium homologue set), never pooled silently. "Eligible homologue" is prespecified at run time and recorded in `docs/data_provenance.md` (planned: full-length or near-full-length homologues passing D-04 screening against UniProtKB/NCBI nr-taxid restricted sets; fragments retained only if they fully contain the peptide window). No partial-credit matching; no Blosum-weighted conservation in the primary metric. `results/tables/epitope_conservation.csv` is the output (1,849 rows; DONE 2026-11-14); `src/R/plot_conservation.R` produces `results/figures/epitope_conservation.png`.

**Rationale.** Exact-fraction is transparent, reproducible without further modelling, and directly interpretable for breadth of protection (a peptide conserved across congeners is a broader vaccine hypothesis). Splitting within/across species separates the DWEL1-specific signal from genus-wide signal — critical because the anchor itself is a partial, unpublished sequence.

**Alternatives considered.** (a) Position weight matrix / BLOSUM62 column scores — reserved as sensitivity, not primary. (b) Percent identity of the peptide window — near-equivalent but admits partial matches; rejected for primary. (c) PSSM (PSI-BLAST) — needs iteration; heavy for 3 proteins; rejected for primary. (d) Presence/absence of protein (orthogroup) as proxy — too coarse for 9–15-mer cores. (e) Pooling species — rejected: hides the within- vs across-species asymmetry that is itself a finding.

**Revisit condition.** Homologue eligibility set changes (new releases, six-GCA annotation); within-species n too small to report (n < ~5 → report with explicit low-n caveat, never silently pool); exact-fraction shown to be flat (peptide in 100% of homologues) → add identity-fraction secondary metric, documented.

**Single source of truth:** `docs/methodology.md` D-08 (definition) + `docs/data_provenance.md` (eligibility + retrieval) + `results/tables/epitope_conservation.csv` (data).

---

## D-09 Prioritisation, exclusion log, sensitivity — EXECUTED 2026-11-14

**Decision.** Candidates are assigned **evidence categories**, not a weighted sum. No numeric weights are introduced anywhere in this project (any such weights would be fabricated precision). Execution (2026-11-14, `src/python/prioritize_candidates.py`, no network calls) applied this prespecified hierarchy top-down; first match wins:

1. **RULE-01 — reference-only invariant first:** `source_protein == SOD1REF` → `exclude_under_prespecified_rule` (C. halotolerans A0AB34KRD1; REFERENCE_ONLY_NOT_C_sphaerospermum). All 722 SOD1REF candidates excluded; logged.
2. **RULE-02 — human exact-match invariant:** `human_exact_match_count >= 1` in reviewed Swiss-Prot → `exclude_under_prespecified_rule` (autoimmunity screen). Exactly 1 PKS1 candidate and 24 SOD1REF candidates carry exact matches; PKS1 → RULE-02, SOD1REF → RULE-01 (invariant order).
3. **RULE-03 — no binding evidence:** no panel allele below 500 nM (MHC-I or MHC-II) → `review_required` (748 candidates). Never reported as a negative result.
4. **RULE-04 — weak tier only:** best IC50 < 5000 nM but none < 500 nM → `retain_with_caveat` (202 candidates).
5. **RULE-05 — binding evidence:** ≥1 panel allele < 500 nM, no human exact match → `retain` (176 candidates).

Prioritised list (`results/tables/prioritized_candidates.csv`, 378 rows = 176 retain + 202 retain_with_caveat) is a deterministic **triage order**, not a score: (i) number of panel alleles with IC50 < 500 nM (desc); (ii) best MHC-I/MHC-II percentile rank (asc); (iii) conservation fraction (desc); (iv) candidate_id (asc, stable). The exclusion log (`results/tables/candidate_exclusion_log.csv`, 1,849 rows) records one row per candidate with rule prefix, reason and reference status; the candidate × allele association table (`results/tables/candidate_allele_associations.csv`, 9,245 rows = 7,440 MHC-I + 1,805 MHC-II) keeps every allele association (blueprint Part 10: do not discard associations). A copy of the master table is at `data/processed/candidate_master_table.csv`.

**Sensitivity scenarios (prespecified; EXECUTED 2026-11-14 as `src/python/sensitivity_analysis.py` → `results/tables/ranking_sensitivity.csv`):**
- PRIMARY: IC50 < 500 nM, full 10-allele panel, MHC-I+MHC-II → 176 retained (Jaccard 1.0 by construction).
- STRICT (< 50 nM) → 33 retained; Jaccard vs PRIMARY 0.1875; top-25 Jaccard 1.0000.
- LENIENT (< 5000 nM) → 378 retained; Jaccard 0.4656; top-25 Jaccard 1.0000.
- RANK05 (rank/percentile ≤ 0.5) → 35 retained; Jaccard 0.1530.
- NO_DQ (DQ alleles removed) → 159 retained; Jaccard 0.9034.
- NO_DRB101 (DRB1*01:01 removed) → 176 retained; Jaccard 1.0000.
- CLASS_I_ONLY → 42 retained; Jaccard 0.2386; CLASS_II_ONLY → 134; Jaccard 0.7614.
Interpretation: the PRIMARY top-25 is stable under affinity-threshold relaxation and DQ/DRB1 removal, but fragile to class restriction (top-25 Jaccard 0.25 for CLASS_I_ONLY) and to rank-vs-affinity rules. BA-vs-EL concordance remains NOT run (EL not executed; documented, not fabricated). Population coverage remains not_evaluated (service HTTP 403/405; blueprint rule applied).

**Rationale.** Categories survive heterogeneous tool outputs and keep every decision auditable back to a rule with a name; weighted sums manufacture false comparability between BA units, EL units, BepiPred outputs and conservation fractions. Blueprint Part 10 specifies the hierarchy; the exclusion log makes it reviewable.

**Alternatives considered.** (a) Weighted scoring / desirability functions — rejected: fabricated weights. (b) Combined p-value (Fisher) across tools — rejected: tools not statistically independent and not all p-valued. (c) ML classifier on published epitopes — rejected: no labelled set for fungal antigens; overfitting certain. (d) Simple allele-count rank as overall score — used only within triage tie-breaks, not as a validated score.

**Revisit condition.** Any sensitivity scenario flips ≥10% of retained/excluded assignments → categories re-examined and re-documented (not patched); exclusion log audit by orchestrator pane w5:p6; new blueprint Part 10 amendment.

**Single source of truth:** `docs/methodology.md` D-09 + `results/tables/candidate_master_table.csv` + `results/tables/candidate_exclusion_log.csv` + `src/R/plot_candidate_priorities.R` + `src/R/plot_ranking_sensitivity.R` (all present, 2026-11-14).

---

## D-10 Structural analysis

**Decision.** Structure-based assessment is performed **only if the predicted peptide core is complete** (MHC-I 9-mer core fully inside the sequence; MHC-II 9-mer core fully inside the 15-mer with flanks on both sides, and the 15-mer itself not crossing a sequence end). Even then, outputs are **hypothesis criteria**, documented as such: e.g., stable binding groove compatibility (NetMHCpan already proxies this), candidate core plausibility vs known HLA groove chemistry, and (optionally) AlphaFold/template checks for β-sheet/β-hairpin context in fungal enzymes. **Structural analysis is NOT scheduled yet** and must not gate or delay D-05/D-06/D-07/D-08. No structure-derived immunogenicity claim may be made from this computational project alone.

**Rationale.** QBZ73626.1 is a partial (235 aa) type III PKS with no solved structure; a partial sequence plus predicted cores is insufficient evidence for structural claims. Blueprint: structural work requires completeness and explicitly-stated hypotheses.

**Alternatives considered.** (a) Model full type III PKS from p-hydroxybenzoate synthase templates (PDB family) — deferred; only under "core complete" unlock. (b) Run it now "for completeness" — rejected: wastes effort on candidates that will fail D-05/D-08; violates ordering. (c) Dimer/homodimer modelling for SOD1 reference — reference-only anyway; no decision-relevant output.

**Revisit condition.** First candidate-master table built AND at least one candidate passes D-09 category A with complete cores → structure work may be proposed (with criteria written before looking at any structure); never before.

**Single source of truth:** `docs/methodology.md` D-10 (policy); any future structural run documented in `docs/data_provenance.md` + this section.

---

## D-12 IEDB API allele naming, DQ chain pairing, and raw-response retention (added 2026-11-14)

**Decision.** At run time, the following were established and recorded (config `prediction_methods.yaml` `iedb_allele_mapping`; `src/python/run_mhcii_predictions.py` docstring; R023):
- MHC-I allele names in `config/hla_panel.csv` are exactly those listed by IEDB (verified against the `netmhcpan_ba-4.1` help allele list, 2026-11-14) and are submitted verbatim.
- MHC-II: IEDB lists DRB1 alleles **without** the `HLA-` prefix (`DRB1*01:01` etc. — verified against the `netmhciipan_ba-4.1` allele list) and accepts DQ **only as α/β chain pairs**. The panel's DQ entries were therefore mapped to documented cis pairs: `HLA-DQB1*02:01 → DQA1*05:01/DQB1*02:01` (DQ2.5) and `HLA-DQB1*03:01 → DQA1*03:01/DQB1*03:01` (DQ7). The untagged panel alleles themselves were not altered; both the panel allele and the submitted chain pair are carried (`allele`/`panel_allele` vs `iedb_allele`) in `results/tables/mhcii_all_predictions.csv`.
- Every raw IEDB response is saved byte-identically under `data/raw/iedb*/…` with per-call sha256 in `results/logs/mh*_predictions.log`; joined tables are written only after (a) all panel calls succeed and (b) every candidate is reproduced exactly (IEDB window coordinates == locally generated 1-based inclusive coordinates; peptide slices identical). No response is imputed, interpolated, or synthesized; any failure aborts without writing a table.

**Rationale.** Chain pairing is required by the IEDB MHC-II model; inventing single-chain DQ typings would silently predict the wrong molecule. Byte-identical retention makes every published number auditable against the service's actual bytes.

**Alternatives considered.** (a) Submit DQ as bare `HLA-DQB1*02:01` — rejected (not in IEDB allele list; would error). (b) Drop DQ alleles from the panel — rejected: panel unchanged, mapping documented instead. (c) Predict with both common DQ α partners for each DQB1 — rejected for primary panel; logged as a sensitivity option. (d) Trust parsed values without checksum retention — rejected: fails the provenance ledger rule (§D of `docs/data_provenance.md`).

**Revisit condition.** IEDB changes the allele list or chain-pair rules → re-map and bump D-12; sensitivity scenario adds alternate DQ2/DQ7 partners → D-09 revisited.

---

## D-13 Population coverage — NOT EVALUABLE (documented, 2026-11-14)

**Decision.** Population coverage (blueprint Part 9) requires the IEDB population-coverage service. Both documented interfaces were probed on 2026-11-14 and are not scriptable from this environment:
- legacy `https://tools-cluster-interface.iedb.org/tools_api/population/` → HTTP 403 Forbidden;
- next-gen `https://nextgen-tools.iedb.org/pipeline` → HTTP 405 Not Allowed.

Per the blueprint rule, a missing predictor output is never treated as a negative result: every candidate keeps `not_evaluated` for this evidence category, no candidate is excluded for it, and **no coverage figure is produced** (fabricating population frequencies or coverage fractions is prohibited). The intended input is preserved (`results/tables/population_coverage_input.csv` = the 10-allele panel), the result table records the not_evaluated status with reason (`results/tables/population_coverage_results.csv`), and `src/R/plot_population_coverage.R` deliberately writes no PNG (it prints SKIPPED and exits 0).

**Re-run condition.** When a documented, authenticated IEDB population-coverage API (or an approved local frequency dataset, e.g., allele-frequency nets with cited populations) is available, execute with the panel in `config/hla_panel.csv`, record endpoint/version/dates/checksums here, then produce Fig 6 from the real output. Until then Fig 6 remains absent by design, not by omission.

**Single source of truth:** `config/populations.yaml` + `results/tables/population_coverage_results.csv` + this section.

---

## D-11 Source-of-truth register (where each decision is revisited)

| Decision | Primary revisit file(s) | Secondary / evidence | Owner until re-assigned |
|---|---|---|---|
| D-01 target identity | `target_selection_rationale.md` | `target_metadata.csv`; `data/processed/target_evidence_matrix.csv`; `data/manifests/retrieval_manifest.csv` (R001–R020) | Pin 1 (target discovery, 2026-11-13) |
| D-02 ambiguous residues | `docs/methodology.md` (this file) | `src/python/validate_sequences.py`; `results/tables/sequence_validation_*.csv`; 3,702 checks PASS / 0 FAIL | DONE 2026-11-13 |
| D-03 human reference | `config/human_proteome.txt` | retrieval ledger row R021 + `HP001` in `data/reference/human_proteome_manifest.csv`; `docs/data_provenance.md` §A #17 | VERIFIED 2026-11-14 (Swiss-Prot 2026_04, 20,429 records, sha256 13852e17…) |
| D-04 similarity + review policy | `data/manifests/retrieval_manifest.csv` | `results/tables/putative_sod_homology_screen.csv`; `results/tables/human_protein_similarity.csv`; `results/tables/candidate_human_similarity.csv`; `results/tables/candidate_exclusion_log.csv` (DONE 2026-11-14); `docs/data_provenance.md` | Pin 2 / orchestrator |
| D-05 MHC-I | `config/hla_panel.csv`; `config/prediction_methods.yaml` | manifest rows R022; `results/tables/mhci_all_predictions.csv`; `results/logs/mhci_predictions.log` | orchestrator (run 2026-11-14) |
| D-06 MHC-II | `config/hla_panel.csv`; `config/prediction_methods.yaml` | manifest rows R023; `results/tables/mhcii_all_predictions.csv`; `results/logs/mhcii_predictions.log`; `docs/methodology.md` D-12 | orchestrator (run 2026-11-14) |
| D-07 B-cell / BepiPred 2.0 | `docs/data_provenance.md` | `results/tables/bcell_residue_scores.csv` (389 rows); `results/tables/bcell_predicted_regions.csv` (8 regions); `data/raw/iedb/bcell/` | DONE 2026-11-14 |
| D-08 conservation | `docs/methodology.md` D-08 | `results/tables/epitope_conservation.csv` (DONE 2026-11-14); `data/reference/homolog_metadata.csv`; `src/R/plot_conservation.R` | DONE 2026-11-14 |
| D-09 prioritisation & exclusion | `docs/methodology.md` D-09 | `results/tables/candidate_master_table.csv` (DONE 2026-11-14); `results/tables/candidate_exclusion_log.csv`; `results/tables/prioritized_candidates.csv`; `results/tables/ranking_sensitivity.csv`; `src/R/plot_candidate_priorities.R`; `src/R/plot_ranking_sensitivity.R` | DONE 2026-11-14 |
| D-10 structure | `docs/methodology.md` D-10 | any future structural outputs; `docs/data_provenance.md` | PENDING — not scheduled |
| D-11 this register | `docs/methodology.md` | whole repo | Pin 3 (docs/R) |
| D-12 IEDB naming / DQ mapping / raw retention | `docs/methodology.md` D-12 | `config/prediction_methods.yaml` `iedb_allele_mapping`; `src/python/run_mhcii_predictions.py`; R023 | orchestrator (2026-11-14) |

## Cross-cutting premises (read before any claim)

- Blueprint document not stored in this repository; Part numbering cited from orchestrator brief (pane w5:p6).
- `data/raw/*.fasta` and `src/python/*.py` belong to other pins; this document cites but never modifies them.
- Every accession, count, and checksum quoted above was verified on disk or via live APIs on 2026-11-13 (manifest) / re-checked or executed 2026-11-14 (this file; IEDB runs R022/R023); anything not yet run is PENDING with owner in D-11.
- No claim of immunogenicity, protection, safety, or clinical relevance may be made: all outputs are computational hypotheses.
