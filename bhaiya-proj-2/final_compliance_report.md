# Cladosporium sphaerospermum Immunoinformatics Project — Final Compliance Report

**Report date:** 2026-11-14 · **Orchestration:** Herdr panes w5:p2/p3/p4 (workers), w5:p6 (orchestrator/verifier) · **Scope:** local only; no GitHub/VCS, no publication artifacts (user instruction).

## 1. Overall status

**The project is COMPLETE to the fullest extent achievable in this environment**, with the following honest exceptions, none of which are silent:

| # | Exception | Status | Why (documented) |
|---|---|---|---|
| E1 | Population coverage (Part 9 / Fig 6) | `not_evaluated` — no figure | IEDB service not scriptable: legacy `tools_api/population/` → HTTP 403; nextgen `POST /pipeline` → HTTP 405 (probed 2026-11-14). `config/populations.yaml`, `docs/methodology.md` D-13 |
| E2 | Antigenicity / allergenicity / toxicity (Part 7 screens, 1,857 rows) | `not_evaluated` — never a negative result | No documented API for VaxiJen/AlgPred/ToxinPred at IEDB or otherwise as of 2026-11-14. `results/tables/immunological_screening.csv` |
| E3 | BA vs EL sensitivity concordance | Not run | EL methods were not executed (user-approved scope = BA primary). Recorded in L-25/L-28 |
| E4 | NCBI BLAST REST cross-check (Part 4) | Not run | Manual browser/RID step; exact-containment (NCBI eutils) used for Part 8 instead. D-04 |
| E5 | Six genome assemblies mining | PENDING (no verified *C. sphaerospermum* SOD1 anywhere; 0 verified homologs) | All 6 assemblies un-annotated. L-05, D-01 |
| E6 | Structural / docking / MD / immune simulation (Parts 16, 19 "Priority 6") | Not scheduled | Insufficient target structural evidence; D-10 unlock criteria unmet |
| E7 | Manuscript / GitHub / publication | Prohibited by user instruction | Local outputs only |

Everything else in the blueprint file set exists and was verified on disk.

## 2. Verified results (all on disk, checksums in `data/manifests/data_checksums.txt` + provenance ledger R001–R027)

### 2.1 Targets & sequences
- **PKS1** (verified *C. sphaerospermum* partial type III PKS, QBZ73626.1, isolate DWEL1, MK213995.1, 235 aa) — the only verified *C. sphaerospermum* protein found anywhere.
- **SOD1REF** (A0AB34KRD1, *C. halotolerans*) — `REFERENCE_ONLY_NOT_C_sphaerospermum` in every derived file, enforcement verified by audit (all 722 rows excluded under RULE-01).
- No *C. sphaerospermum* SOD1 exists in UniProt/NCBI/local screens (0 hits).
- Sequence validation 1,849-candidate framing; 3,702 checks PASS / 0 FAIL; 19/19 pytest.

### 2.2 Candidates
- **1,849 peptide candidates**: 1,488 MHC-I (8–11 mers) + 361 MHC-II (15 mers); coordinates validated against source sequences.

### 2.3 Human-proteome screen (D-04)
- Swiss-Prot 2026_04 (20,429 reviewed records; sha256 in manifest; HP001). 25 exact-match candidates (1 PKS1 → excluded RULE-02; 24 SOD1REF → excluded RULE-01). 0 binders <500 nM are human-exact-matches.

### 2.4 IEDB predictions (executed 2026-11-14, user-approved)
- MHC-I `netmhcpan_ba-4.1`: 40 calls; 7,440 prediction rows; every candidate reproduced exactly.
- MHC-II `netmhciipan_ba-4.1`: 10 calls; 1,805 rows; every candidate reproduced exactly.
- BepiPred-2.0 (bcell endpoint): 6 raw calls (2 targets, byte-identical re-runs during debugging); 389 residue scores; 8 predicted regions (PKS1 5; SOD1REF 3, reference-only).
- 56 HTTP-200 calls logged (40+10+6); `prediction_failures.csv` = 0 rows; raw responses byte-identical, sha256-verified; manifest `results/tables/prediction_run_manifest.csv`.

### 2.5 Conservation (D-08, executed; MAFFT unavailable → exact-containment)
- PKS1: **0/1,127** peptides exactly conserved in 50 Cladosporium PKS homologs (109 found, top 50 fetched).
- SOD1REF: **350/722** conserved in 9 Capnodiales SOD homologs (reference context only).
- Raw FASTAs saved (`data/raw/ncbi/`); homolog metadata at `data/reference/homolog_metadata.csv`; identity column honestly `not_evaluated (no BLAST)`.

### 2.6 Prioritisation (D-09, executed 2026-11-14)
| Decision | Count | Rule |
|---|---|---|
| `exclude_under_prespecified_rule` | 723 (722 SOD1REF RULE-01 + 1 PKS1 RULE-02) | reference-only / human exact match |
| `review_required` | 748 | RULE-03 (no allele <500 nM) — **not a negative result** |
| `retain_with_caveat` | 202 | RULE-04 (<5000 nM only) |
| `retain` | 176 | RULE-05 |

- Prioritized list: 378 rows (`prioritized_candidates.csv`); deterministic triage (allele count → best rank → conservation → candidate_id); no weighted composite scores.
- Exclusion log: 1,849 rows (`candidate_exclusion_log.csv`); candidate × allele associations: 9,245 rows (`candidate_allele_associations.csv`).

### 2.7 Sensitivity (Part 10; `ranking_sensitivity.csv`)
- PRIMARY 176 → STRICT(50 nM) 33 (top-25 stable, Jaccard 1.0); LENIENT(5000 nM) 378; RANK05 35; NO_DQ 159 (0.9034); NO_DRB101 176 (1.0000); CLASS_I_ONLY 42 (top-25 Jaccard 0.25); CLASS_II_ONLY 134 (0.7614). Full interpretation in D-09.

### 2.8 Figures (10 real PNGs; Fig 6 absent by design)
Fig 1 workflow diagram · Fig 2 protein characteristics · Fig 3 MHC-I/II counts + rank distributions (2×2 PNGs) · Fig 4 BepiPred tracks · Fig 5 conservation histogram · Fig 7 evidence matrix (top 25) · Fig 8 ranking sensitivity. All R-script-generated, checksummed, verified as real PNGs.

## 3. Provenance & integrity
- R001–R027 ledger; `data/manifests/data_checksums.txt` (29 files); every IEDB/NCBI call logged with sha256; raw responses byte-identical.
- SOD1REF reference-only discipline verified in every table audited (all 722 excluded by invariant rule).
- No fabricated data; no invented predictor outputs; `not_evaluated` ≠ negative anywhere.

## 4. Blueprint Parts 11–19
- Part 15 (3D construct/vaccine design), Part 16 (structural/docking/MD), Part 19 ("what would make it impressive"): **not executed** — require validated candidates and structural evidence; honest scope exclusion (user instruction prohibits publication artifacts; scientific scope requires target completeness we do not have).
- Part 18 (resources/docs): satisfied by `docs/` set + this report.

## 5. Input → output, one line each
**Inputs:** 2 target FASTAs (PKS1 verified partial; SOD1REF reference-only) + Swiss-Prot 2026_04 + 10-allele HLA panel + IEDB (BA 4.1, BepiPred-2.0) + NCBI eutils.
**Outputs:** 1,849-candidate table with IEDB binding predictions (7,440 + 1,805 rows) → BepiPred-2.0 regions (8) → conservation stats → RULE-01..05 prioritisation (176 retain + 202 caveat) → 10 real figures → full audit trail. All computational hypotheses; **nothing is a validated epitope; nothing may be cited as vaccine/drug evidence without wet-lab work.**
