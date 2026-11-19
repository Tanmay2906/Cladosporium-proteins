# Scientific Limitations — *Cladosporium sphaerospermum* Immunoinformatics Blueprint

Compiled 2026-11-14 by Pin 3 (docs/R) from verified project evidence. Every limitation is stated plainly; none is papered over by wording. Where a limitation is already recorded in `target_selection_rationale.md`, `data/manifests/retrieval_manifest.csv` or `config/human_proteome.txt`, the path is cited. Blueprint Part numbering is per the orchestrator brief (the blueprint document itself is not stored in this repository).

## L-01 — The PKS anchor is a *partial type III* PKS, not the canonical PKS1
QBZ73626.1 (235 aa; `cond_enzymes`, codon 1..235, 5′ ends missing, `no-ends`) is the only verified *C. sphaerospermum* PKS-class record — **not** the iterative type I "PKS1" the project name anticipates. Any domain/active-site or mechanism language must be framed on type III chemistry. (Evidence: R008/R009/R014; `data/raw/PKS1.fasta`; `target_selection_rationale.md` §1.2, §2.)

## L-02 — The PKS record is an unpublished direct submission
MK213995.1 (Sanger dideoxy; Manoharan G, Ramakrishnan D, Marimuthu J — *unpublished*, 22-Nov-2018; National Institute of Ocean Technology, Chennai). Peer review absent; sequence could be withdrawn/revised. All claims derived from it are hypothesis-generating only. (Evidence: R009; `data/processed/target_evidence_matrix.csv` row 1.)

## L-03 — No functional or melanin-pathway assay exists
No melanin production assay, no enzyme kinetics, no expression data, no localisation experiment (predicated cytoplasmic, no targeting signal) for QBZ73626.1 or for *C. sphaerospermum* targets of any kind. The record was found via a `polyketide` query bucket — **not** via a melanin bucket (R013 count 0). The inference "QBZ73626.1 is a melanin-pathway enzyme" is **unsupported**; "virulence factor" is unsupported. (Evidence: R009, R013, R019; evidence matrix; rationale §4.)

## L-04 — SOD1 was not found in *C. sphaerospermum* — zero hits anywhere
UniProt: no *sphaerospermum* SOD1 (R001/R002; genus only: A0AB34KRD1 *C. halotolerans*). NCBI with organism filter: count 0 for SOD1/superoxide/dismutase (R010), catalase (R011), glutathione (R012), melanin (R013), oxidative-stress terms. Local BLOSUM62 screen vs the 324-record NCBI proteome: **0 hits above threshold** (identity ≥30% AND reference coverage ≥40%); best coverage 64.29% at 23.31% identity (YFQ81086.1/YFQ81084.1) — RNA-polymerase-II/beta-tubulin fragments, not SODs. **Absence of evidence is not evidence of absence**: the proteome set is dominated by short Sanger/WGS fragments (mean 76.1 aa), lacks isolate/strain tokens, and excludes any undeposited *T*-lineage genomes. (Evidence: R005, R010–R013, R017–R020; `results/tables/putative_sod_homology_screen.csv`; rationale §1.2, §4.)

## L-05 — The six unannotated genome assemblies are unverified and unused
GCA_001931905.1 (IMV00045, Chernobyl Unit-4 turbine-hall wall, collected 1998-06-16), GCA_000261425.2 (UM 843), GCA_019191205.1, GCA_023621355.1, GCA_052323905.1, GCA_000261425.1 (earlier version of the UM 843 assembly — duplicated, deduplicate) are **not present locally**, were **not queried this session**, and carry no local annotation. Any statement about their PKS/SOD content today would be fabricated. They are listed only as candidate mining sources. (Brief-supplied; pending retrieval/annotation.)

## L-06 — Halotolerans references are NOT *sphaerospermum* targets
A0AB34KRD1 (SOD1, 154 aa, complete, PE=3 inferred-from-homology) and A0AB34KVX2 (non-reducing type I PKS1, 2175 aa, **WGS-derived, "CAUTION: WGS preliminary data"**) are *C. halotolerans* TM138-S3 — a different species. They may inform domain models and genus-level comparison only. Neither may be substituted for a verified target (blueprint Part 1 / design constraint; rationale §2–§3). Any metric computed on them must be labelled reference-based.

## L-07 — The reference screen’s scope is the named proteome only
`data/reference/c_sphaerospermum_ncbi_proteome.fasta` (324 records, taxid 92950, all `.1` versions) covers named NCBI records; it does **not** include: undeposited genomes, the six GCAs (L-05), UniProt *sphaerospermum* (23 records, **enumeration not completed** — stream endpoint returned HTTP 400, documented), or any SRA-only assemblies (DWEL1 SRA status unknown). A complete proteome could change every downstream result. (Evidence: R017 premise note — quoted 641-record set was an unquoted OR-match including ~317 other-organism records and was NOT downloaded; R018 isolate tokens absent except DWEL1; rationale §1.2, §5.)

## L-08 — Predictions are not immunogenicity
NetMHCpan/NetMHCIIpan BA (and EL) outputs are binding/presentation likelihoods for HLA molecules. They do **not** measure T-cell activation, cytotoxicity, memory, protection in an animal model, or clinical efficacy. BepiPred 2.0 linear-epitope scores are sequence propensities, not antibody binding evidence. No wet-lab work exists or is implied by this repository. Ranking = computational hypothesis only.

## L-09 — A human-similarity screen is not a safety proof
D-03/D-04 human screening (UniProtKB/Swiss-Prot UP000005640, reviewed only) can flag obvious homologies to human proteins/peptides. It cannot detect: cross-reactivity via non-homologous mimicry, post-translational epitope chemistry, unlisted human proteoforms, tissue-specific expression cross-reactivity, or off-target T-cell epitopes without homology. Passing the screen ⇒ "no flagged homology", never "safe". The screen has **not been run** (human proteome NOT retrieved; `config/human_proteome.txt` status NOT YET RETRIEVED 2026-11-13).

## L-10 — Homology-screen thresholds are arbitrary but prespecified
Identity ≥30% AND reference coverage ≥40% (R020) were chosen as prespecified operational cut-offs, not from a validated ROC. Near-boundary behaviour is unquantified; different defensible thresholds (e.g., 25/30, 35/50) would admit different record sets (see D-09 sensitivity scenarios; 30 records ≥40% coverage but ≤~25% identity are currently excluded). Prespecification is the only protection against post-hoc tuning; it does not make the threshold scientifically validated.

## L-11 — FungiDB inaccessible — cross-evidence base is incomplete
The FungiDB portal at `https://fungidb.org/fungidb/app/` is an interactive JS application; no unauthenticated programmatic protein-record retrieval succeeded (2026-11-13, R016). A third, independent evidence source was therefore unavailable; UniProtKB + NCBI eutils + InterPro-via-xrefs only. Genus annotation completeness in FungiDB may differ materially from both.

## L-12 — Human proteome retrieval is deliberately held and currently absent
See D-03/L-09. Release, size, checksums: unknown (record at download). Screen, EL-vs-BA concordance, and every "human-similarity" flag shown in any downstream table would inherit this absence — until retrieval, the column does not exist.

## L-13 — Conservation metrics are exact-match and low-n
D-08 counts exact substring matches within eligible homologues, split within- vs across-species. It ignores conservative substitutions (a Thr→Ser peptide scores 0 outside exact copies) and is subject to tiny denominators: within-species *sphaerospermum* homologues may number 0–1 today (there is no complete proteome); across-species genus-wide congeners dominate. A "100% conserved" call with n=1 is meaningless; the D-08 rule mandates explicit low-n caveats and never pooling across species buckets.

## L-14 — MHC panel is European/South-Asian-biased; coverage incomplete
Panel alleles (A*02:01, A*01:01, B*07:02, B*08:01, C*07:01, DRB1*01:01/*04:01/*07:01, DQB1*02:01/*03:01) reflect European/South Asian frequencies (AFND rationale strings, `config/hla_panel.csv`). Frequencies elsewhere (e.g., East/Southeast Asian, African, Native American) differ; no inference may be made for those populations. HLA-C is underpowered generally in epitope literature. Minor alleles are excluded by design, with a D-09 sensitivity trim scenario as the only concession.

## L-15 — No BA/EL cross-method concordance yet; BA primary by fiat
Per `config/prediction_methods.yaml`, BA is retained for primary screening and EL only as a secondary sensitivity check. This is a defensible but conventional choice; concordance between BA and EL rankings (and hence the robustness of any category assignment) is unmeasured until D-09 scenarios run.

## L-16 — BepiPred 2.0 version pinned only at run time; linear epitopes only
Exact tool/version/endpoint are **not yet recorded** (PENDING, D-07). IEDB endpoints change without notice (config header warns). Linear epitopes ignore conformational/discontinuous epitopes entirely; conformational prediction requires a structure (not available; D-10 locked). Antibody-epitope ≠ T-cell epitope — these outputs must not be merged.

## L-17 — This project is computational-only
No wet-lab, animal, ex-vivo, or clinical component exists or is planned here. Nothing in this repository constitutes medical, diagnostic, or veterinary guidance; do not present it as such.

## L-18 — Evidence hierarchy weakness: predicted-only annotations
A0AB34KRD1 protein existence = 3 (inferred from homology), no experimental characterisation; A0AB34KVX2 existence = 4 (predicted), WGS preliminary; both UniProt TrEMBL unreviewed. Only domain-level cross-checks (InterPro/CDD/PROSITE/Pfam, R015) support class calls. There is no peptide-coincubation, no antifungal assay, no IgG titre, no T-cell assay anywhere.

## L-19 — UniProt stream endpoint failure and 23-record gap
The *sphaerospermum* UniProtKB record enumeration could not be produced (`stream` HTTP 400, documented in rationale §1.3). The claim "UniProt has 23 records, zero named PKS1/SOD1" rests on search-API evidence, not exhaustive enumeration — adequate, but weaker than a complete dump.

## L-20 — No version control / reproducibility rests on manual ledgers
No git (by instruction), no CI, no container pinning. Reproducibility rests solely on `data/manifests/retrieval_manifest.csv`, `config/human_proteome.txt`, `docs/data_provenance.md` and file checksums. Bit-rot of the manifest (hand edits) would be undetectable except by checksum diff.

## L-21 — Blueprint not in-repo; methodological citation risk
Blueprint Parts referenced here (1, 10, 17) and review-policy category names are cited from the orchestrator brief. If the blueprint text differs, these docs must be re-aligned; no blueprint PDF/DOCX is stored in this repository to check against.

## L-22 — Local environment drift risk
Checksums and versions (Python 3.13.5, Biopython 1.88, pandas 3.1.1, R 4.4.2, ggplot2/dplyr/readr pre-installed) are recorded 2026-11-14. Package upgrades (notably pandas ≥3 and Bio.Align API drift) may break `src/python/*.py`; scripts are skeletons/pending runs and must be re-verified before any production run.

## L-23 — Type III PKS chemistry may exclude melanin-like products
Type III (chalcone/stilbene-like, "cond_enzymes") enzymes typically act on small-molecule acyl-CoA substrates; whether the DWEL1 enzyme touches a melanin pathway is unknown (L-02, L-03). Any vaccine-antigen narrative built on "melanin enzyme" language would exceed the evidence; permissible language is "partial type III polyketide synthase, unpublished, functional assignment probable by domain family only".

## L-24 — Low-complexity/composition artifacts
Short fragments (min 4 aa, R017) include likely spurious/ORF-noise records; inclusion of such fragments in any denominator (e.g., D-08 eligibility) would distort fractions. Eligibility rules must exclude sub-window fragments and non-protein artifacts before denominators are fixed — prespecify at run time; do not decide post hoc.

## L-25 — IEDB BA predictions executed 2026-11-14; what they can and cannot support
The IEDB calls were executed against NetMHCpan 4.1 BA and NetMHCIIpan 4.1 BA. **Binding-affinity predictions are not presentation, processing, or immunogenicity measurements.** High affinity does not imply in-vivo generation of the peptide by proteasomal processing, TAP transport, HLA stabilisation, TCR recognition, or any cellular/physiological outcome. Percentile-rank and IC50 cut-points (e.g. 50/500 nM; NetMHCpan-style %Rank ≤0.5/≤2) are tooling conventions for triage, not validated thresholds for this organism or this protein pair. EL-mode (elution/epitope) methods weight processing likelihood and would re-rank candidates; BA-vs-EL is prespecified as a sensitivity scenario (D-09) but is NOT run. Any candidate "binder" highlighted in figures is a computational hypothesis only. SOD1REF columns everywhere remain REFERENCE_ONLY_NOT_C_sphaerospermum (*C. halotolerans* A0AB34KRD1) — no C. sphaerospermum SOD1 predictions may be cited. The 235-aa PKS1 is a partial type III PKS fragment — every prediction on it inherits partial-sequence risk, and candidate coordinates are reported against QBZ73626.1 (NCBI residue numbering), not against any full-length orthologue.

## L-26 — DQ allele mapping is a documented assumption, not a panel change
**IEDB MHC-II** accepts DQ alleles only as α/β chain pairs, so the approved panel entries `HLA-DQB1*02:01` and `HLA-DQB1*03:01` were submitted as `DQA1*05:01/DQB1*02:01` (cis DQ2.5) and `DQA1*03:01/DQB1*03:01` (cis DQ7) respectively, and DRB1 alleles were submitted without the `HLA-` prefix, exactly as IEDB's allele list spells them (verified 2026-11-14). Alternate defensible partners (e.g., DQA1*02:01/DQB1*02:01 for DQ4; DQA1*06:01/DQB1*03:01) exist and would score differently; the panel was NOT changed, only its submission names — but downstream claims must repeat this caveat verbatim. Any other DQA1/DQB1 partner is an outstanding sensitivity question, deliberately logged here and in D-12.

## L-27 — IEDB service availability and silent drift
IEDB tools-API responses were captured byte-identically with per-call sha256 (50 raw files, 0 mismatches on re-check) — but the exact model parameters inside IEDB's 2023.09/4.1 releases are opaque to this project and may change without notice; a re-run at a later date may produce different values for identical inputs. Any future comparison across time must record IEDB's documented release at that moment or it is invalid. Raw responses are the ground truth for this run; do not re-derive them from the figures.

## L-28 — Prioritisation rests on BA-only binding evidence; EL, coverage, antigenicity and structure are open
The Part 10 hierarchy (D-09) ranks candidates on NetMHCpan/NetMHCIIpan **BA** affinity, the BepiPred-2.0 region overlap, the human exact-match screen and exact-containment conservation. **EL (eluted-ligand) methods were not executed**; the BA-vs-EL sensitivity leg is therefore documented but absent, and BA-only triage can favour peptides that bind yet are not naturally presented. Population coverage is not_evaluated (D-13, service 403/405); antigenicity/allergenicity/toxicity predictors have no documented API in this environment (Part 7, all not_evaluated); no within-*C. sphaerospermum* SOD1 homolog exists (zero verified records), so SOD1REF conservation is reference-context only; no structural, docking, MD or immune-simulation evidence exists (Part 16 advanced analyses not scheduled). The 176 `retain` / 202 `retain_with_caveat` candidates are triage ordering for review, not validated epitopes; any downstream use must re-state every caveat in L-25/L-26/L-27 verbatim.
