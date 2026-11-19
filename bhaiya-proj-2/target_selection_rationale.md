# Target Selection Rationale — *Cladosporium sphaerospermum* Immunoinformatics Blueprint (Part 1: Target Discovery)

Date: 2026-11-13  |  Agent: Pin 1 — Target Discovery  |  Scope: biological target verification for PKS1 / SOD1 (primary) and close-alternatives (melanin / oxidative-stress defence).

## 1. What was verified, where, and how

All verification was done against **live databases** on 2026-11-13, using the UniProtKB REST API and NCBI eutils (esearch/efetch/esummary). InterPro evidence was cross-checked through UniProtKB xrefs (InterPro, Pfam, PROSITE, CDD, FunFam, TIGR04532, NCBI GenBank "cond_enzymes"). FungiDB was attempted but could not be queried as an unauthenticated machine endpoint (portal is a JS frontend); its absence is documented in the manifest (R016) rather than worked around.

### 1.1 UniProt — PKS1 / SOD1 homologues in the genus *Cladosporium* (any species)

| Target | Organism | Strain | Accession (entry/seq version) | Length (aa) | Reviewed | Protein existence | Gene | Orf |
|---|---|---|---|---|---|---|---|---|
| PKS1 homologue | *C. halotolerans* (not sphaerospermum) | TM138-S3 | A0AB34KVX2 (entry v7, seq v1) | 2175 | TrEMBL | 4: Predicted | PKS1 | NWR41_04342 |
| SOD1 homologue | *C. halotolerans* (not sphaerospermum) | TM138-S3 | A0AB34KRD1 (entry v7, seq v1) | 154 | TrEMBL | 3: Inferred from homology | sod1 | WHR41_05553 |
| — (out-of-genus best negative control) | *Dothistroma septosporum* | NZE10 / CBS 128990 | N1PC72_DOTSN / sod1 analogue | 2214 | TrEMBL | 3: Inferred from homology | pks1 | DOTSEDRAFT_47338 |
| — (out-of-genus) | *Fulvia fulva* | (tomato leaf mold) | A0A9Q8UVR0 | 2215 | TrEMBL | 4: Predicted | — | CLAFUR5_12824 |

Key details for the genus-relevant entries:

- **A0AB34KVX2 (*C. halotolerans* PKS1)** — annotated "Non-reducing polyketide synthase PKS1"; catalytic reaction Rhea RHEA:22836 (a fatty acyl-[ACP] + malonyl-[ACP] → 3-oxoacyl-[ACP]); cofactor pantetheine 4′-phosphate; domains KS3 (PROSITE PS52004, 370-802), PKS/mFAS dehydratase (PS52019, 1292-1602), 2× Carrier (PS50075, 1666-1743 / 1797-1874), N- and C-terminal hotdog folds, active sites 542/1324/1515, phosphopantetheine 1703/1834; FunFam 1.10.1200.10/3.10.129.110 (sterigmatocystin-PKS-like); TIGR04532 PT_fungal_PKS; CDD cd00833; IPR030918 PT_fungal_PKS; NCBI RefSeq XP_069230389.1. CAUTION flag: WGS-derived preliminary sequence; **verified domain architecture only, no experimental localisation or expression data for this strain**.
- **A0AB34KRD1 (*C. halotolerans* SOD1)** — recommended name "**Superoxide dismutase [Cu-Zn]**", EC 1.15.1.1; function annotation ("destroys radicals..."), homodimer (ARBA00011738); cofactors Cu 1+ ×1 (ChEBI CHEBI:23378) and Zn2+ ×1 (CHEBI:29105); domain Superoxide dismutase copper/zinc binding (Pfam PF00080; 14-150); InterPro IPR001424/IPR036423/IPR018152/IPR024134; CDD cd00305; C-terminal disordered (MobiDB 125-154); decision logic inferred-from-homology (evidence ECO:0000256/RU000393, experimental: NONE). **Complete sequence (154 aa, MD5 8E003FE5E9060056B5542B38A427F62B, crc64 8E2025D407DCB70D)**.

**No *C. sphaerospermum* record in UniProt carries PKS1 or SOD1.** This matches the project premise exactly (23 UniProtKB records for *C. sphaerospermum*; 0 named PKS1/SOD1).

### 1.2 NCBI Protein — PKS1 / SOD1 in *C. sphaerospermum* specifically

- **PKS1** — esearch `"Cladosporium sphaerospermum"[Organism] AND (PKS1[All Fields] OR SOD1[All Fields])` → **count 0**. Broadening to `"Cladosporium sphaerospermum"[Organism] AND (melanin OR polyketide)` → **count 1**: **QBZ73626.1**, "type III polyketide synthase, **partial** [Cladosporium sphaerospermum]", isolate **DWEL1**, length **235 aa** (5′ ends missing; `no-ends` per esummary), code MK213995.1, Sanger dideoxy (Manoharan et al., *unpublished* / direct submission, 22-Nov-2018), feature `cond_enzymes` (CDD:447866; InterPro cl09938), taxon 92950. **This is a type III PKS, not the canonical iterative type I "PKS1"** that the project's expectation appears to assume; it is, however, the only genuine verified PKS-class sequence from *C. sphaerospermum*.
- **SOD1** — esearch with organism filter, none of the buckets return a hit:
  - `"Cladosporium sphaerospermum"[Organism] AND (SOD1 OR superoxide OR dismutase)` → **count 0**
  - `"Cladosporium sphaerospermum"[Organism] AND catalase` → **count 0**
  - `"Cladosporium sphaerospermum"[Organism] AND glutathione` → **count 0**
  - `"Cladosporium sphaerospermum"[Organism] AND melanin` → **count 0**
  - `"Cladosporium sphaerospermum"[Organism] AND oxidative` → **count 0**
- **Nucleotide caution**: `MK213995.1` is unpublished Sanger-direct-submission sequence — flagged with that evidence status, not hidden.

### 1.3 Alternative stress-associated proteins in *C. sphaerospermum* (step 3 actions)

Searches were run for melanin-biosynthesis and oxidative-stress defence functions.
- Polyketide / PKS bucket → QBZ73626.1 (type III PKS) — retained.
- Melanin, catalase, glutathione, oxidative, superoxide, dismutase buckets → **0 records**.
- UniProt stream `organism:"Cladosporium sphaerospermum"` could not be produced programmatically (REST `stream` returned HTTP 400 for the organism-quoted query and enumeration of the 23 records could not be retrieved within this session); NCBI record set was used instead. Documented as a limitation.
- **FungiDB** — the public portal is an interactive JS app; programmatic/REST access to protein records without an authenticated Bioproject/session was not possible. Not used as evidence; noted in manifest (R016). InterPro evidence (via UniProt and CDD xrefs) serves as the structural/domain cross-check.

## 2. Verification outcomes vs. blueprint points 186-188

The blueprint states: *"Do not manufacture a sequence or combine different homologues into a synthetic target. Document the limitation and reassess the target set."* Both rules apply:

1. **PKS1 (as the project names it) does NOT verify for *C. sphaerospermum*.** The only verified PKS-class sequence from the species is the **partial type III PKS QBZ73626.1 (235 aa)**. The complete, domain-validated iterative type I PKS1 (A0AB34KVX2) comes from **C. halotolerans**, a different species. It is therefore **reference-only** and may not be presented as a verified *C. sphaerospermum* PKS1.
2. **SOD1 does NOT verify for *C. sphaerospermum*.** Exactly 0 records across all superoxide/dismutase/oxidative/catalase/glutathione buckets; the only Cladosporium-genus SOD1 (A0AB34KRD1, 154 aa, Cu-Zn) is from **C. halotolerans**. It is likewise **reference-only**.
3. **No synthetic/homologue-stitched sequence was produced.** `data/raw/PKS1.fasta` holds QBZ73626.1 exactly as submitted; `data/raw/SOD1.fasta` contains the A0AB34KRD1 sequence explicitly labelled as a reference from *C. halotolerans*, **not** as a *C. sphaerospermum* target.
4. **Limitation is documented here and in the evidence matrix**; target set is reassessed below.

## 3. Reassessed target set (final)

| Rank | Target | Species | Status | Role in the blueprint | Decision |
|---|---|---|---|---|---|
| 1 | **QBZ73626.1** (type III PKS, partial, 235 aa) | *C. sphaerospermum* | Verified, partial, unpublished | PKS class target; anchor for downstream epitope / immunoinformatics pipeline | **include, with explicit "partial type III PKS" caveat** |
| 2 | **A0AB34KRD1** (Cu-Zn SOD1, 154 aa, complete) | *C. halotolerans* (reference) | Verified, complete, genus-level | Reference for domain/active-site modelling and for genus-level SOD epitope comparison; NOT the target species | **include as reference, not as *C. sphaerospermum* target** |
| 3 | **A0AB34KVX2** (non-reducing PKS1, 2175 aa, complete, WGS-derived) | *C. halotolerans* (reference) | Verified, complete, genus-level, WGS-CACUTION | Reference for domain architecture of classic fungal PKS1 class | **include as reference, not as *C. sphaerospermum* target** |
| 4 | SOD1 in *C. sphaerospermum* | *C. sphaerospermum* | ✗ not found | — | **exclude; documented limitation** |

### Decisions (final target set)

1. **Retain PKS1-class work** using QBZ73626.1 as the verified *C. sphaerospermum* anchor, augmented by A0AB34KVX2 as reference domain architecture. **Do not** combine them into one synthetic sequence.
2. **Do not proceed with an unverified SOD1** for *C. sphaerospermum*. Proceed only with the genus-level reference A0AB34KRD1 for comparative purposes; flag every downstream result (epitope priority score, conservation, hydrophobicity) as reference-based.
3. **Before downstream blast/epitope selection** for *C. sphaerospermum*, obtain the complete *C. sphaerospermum* proteome (e.g., from the DWEL1 isolate or a related Sanger/WGS assembly deposited in the NCBI SRA/GenBank) and re-run the same three verification rounds; if a complete iterative PKS1 and/or a canonical Cu-Zn SOD1 is found, replace the anchor in this file and in the FASTA, matrix, and manifest (all via the retrieval pipeline, never by patching).
4. **Literature caution**: the only *C. sphaerospermum* PKS record (MK213995.1) is an unpublished direct submission. Any biological claims derived from it must be framed as hypothesis-generating.

## 4. Confidence summary

- **PKS1**: moderate for the class call (type III PKS, partial); high for the reference homologue's domain model; **low for any claim that QBZ73626.1 is a melanin-pathway enzyme or virulence factor** (no functional assay, no expression data, unpublished source).
- **SOD1**: **non-verified for *C. sphaerospermum*** (count = 0 in all relevant buckets); high confidence in the genus reference (A0AB34KRD1) from annotation perspective, low for species-specific conclusions.
- **Localisation evidence** for every entry above is predicted (cytoplasm / no targeting signal); **no experimental localisation, expression, or immunological assay exists** for this organism.
- **Literature evidence** for *C. sphaerospermum* targets is **effectively absent**; strongest literature support comes from out-of-genus *Verticillium dahliae* VdPKS1 studies and *Gioti et al. 2020* (genus draft genome).

## 5. Next steps required from the project

1. Obtain/verify *C. sphaerospermum* DWEL1 complete proteome (NCBI SRA / GenBank accession MK213995.1, BioProject implied by direct submission).
2. Re-run verification with exact isolation source; if complete PKS1 / SOD1 found, promote them over the reference entries and update all manifests via this retrieval pipeline.
3. Cross-source InterPro domain calls for both classes with experimental localisation assays (e.g., GFP-tagging) before epitope prioritisation.
4. Do **not** proceed with immunogenicity scoring on a hybrid/synthetic PKS or on unverified SOD1; holds blueprint integrity.

---
**Bottom line**: The only verified PKS-class target from *C. sphaerospermum* is QBZ73626.1, a partial type III PKS. No SOD1 (and no catalase or glutathione) is verifiable for this species today. Genus-level references (A0AB34KVX2 / A0AB34KRD1, both *C. halotolerans*) are documented, domain-validated, and explicitly labelled as references — not as *C. sphaerospermum* targets. The target set is reduced to one class anchor plus two references, and the gap is documented rather than papered over, in line with blueprint points 186-188.
