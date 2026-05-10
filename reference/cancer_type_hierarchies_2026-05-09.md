# 29 BIDIFAC+ TCGA Cancer Types — Candidate-Hierarchy Reference Table

**Compiled**: 2026-05-09
**Cohort**: 29 TCGA cancer types as used in Lock/Park/Hoadley BIDIFAC+ (AOAS 2022) and Samorodnitsky/Hoadley/Lock (BMC Bioinformatics 2022). BIDIFAC+ merges COAD + READ into a single "CORE" unit and excludes GBM, LAML, UVM relative to the broader 33-type TCGA Pan-Cancer Atlas.

**Purpose**: Reference for structural-hierarchy-candidate identification against the C[j,k] inclusion vectors emerging from BIDIFAC+. Compile-only deliverable; contested assignments flagged.

---

## Master reference table

| TCGA | Full name (WHO-aligned) | Tissue: Coarsest (anatomical system, 5) | Tissue: Mid (organ, 11) | Tissue: Finest (specific tissue, 29) | Epi: Coarsest (3) | Epi: Mid (6) | Epi: Finest (10) | Germ: Coarsest (3) | Germ: Mid (sub-divisions) | Sex composition | Hoadley 2018 supercluster | Hoadley 2018 iCluster (28) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **ACC** | Adrenocortical carcinoma | Endocrine | Adrenal | Adrenal cortex | Epithelial | Adenocarcinoma (steroid-producing) | Adrenocortical | Mesoderm | Intermediate mesoderm (urogenital ridge) | Both | (mixed CNS/endocrine) | C9 (with KICH) [2A] |
| **BLCA** | Bladder urothelial carcinoma | Urogenital | Bladder | Bladder urothelium | Epithelial | Transitional/urothelial | Urothelial carcinoma | Endoderm — *contested with mesoderm* † | Cloacal endoderm + mesoderm | Both | Pan-squamous (subset) | Distributed: C25 / C27 / C20 [2B] |
| **BRCA** | Breast invasive carcinoma | Reproductive (female-dominant) | Breast | Breast (mammary gland) | Epithelial | Adenocarcinoma | Ductal/lobular (mostly luminal) | Ectoderm — *with mesodermal stromal contribution* † | Surface ectoderm | Female (>95%) ‡ | Pan-gyn (in part) | C19 (luminal) / C2 (HER2-amp) / split [2B, Fig 2C] |
| **CESC** | Cervical squamous + endocervical adenoCa | Reproductive (female) | Gynecological | Cervix | Epithelial | Squamous + Adenocarcinoma (mixed) | Cervical SCC + endocervical adenoCa | Endoderm + mesoderm — *contested* † | Cloacal endoderm + Müllerian (lateral plate mesoderm) | Female | Pan-squamous (HPV) | C27 (Pan-SCC HPV) [2B] |
| **CHOL** | Cholangiocarcinoma | Digestive | Hepatobiliary + pancreas | Bile duct | Epithelial | Adenocarcinoma | Cholangiocarcinoma | Endoderm | Foregut endoderm | Both | (Pan-GI fringe) | Diverse — <50% in one cluster [2B, text] |
| **CORE** | Colorectal (COAD + READ merged) | Digestive | GI tract (lower) | Colorectum | Epithelial | Adenocarcinoma | Colorectal adenoCa | Endoderm | Hindgut endoderm | Both | Pan-GI | C4 (Pan-GI CRC) / C18 (Pan-GI MSI) [2B] |
| **DLBC** | Diffuse large B-cell lymphoma | Hematological | Lymphoid | Lymph node B-cell | Hematological | Hematopoietic | B-cell lymphoma | Mesoderm | Hematopoietic mesoderm | Both | (Mesenchymal/immune) | C3 / C20 [2B] |
| **ESCA** | Esophageal carcinoma | Digestive | GI tract (upper) | Esophagus | Epithelial | Squamous + Adenocarcinoma (split) ‡ | ESCC + Esophageal adenoCa | Endoderm | Foregut endoderm | Both | Pan-GI / Pan-squamous (mixed) | C1 / C4 / C18 / squamous clusters [2B] |
| **HNSC** | Head and neck squamous cell carcinoma | Head & neck | Head/neck | Oropharyngeal/laryngeal mucosa | Epithelial | Squamous | HNSCC (HPV+/HPV−) | Ectoderm + endoderm — *contested* † | Surface ectoderm + foregut endoderm (depending on subsite) | Both | Pan-squamous | C10 / C25 / C27 (HPV+) [2B] |
| **KICH** | Kidney chromophobe | Urogenital | Kidney | Distal nephron / collecting duct | Epithelial | Chromophobe carcinoma | Chromophobe RCC | Mesoderm | Intermediate mesoderm (metanephric) | Both | (clusters with ACC) | C9 [2B] |
| **KIRC** | Kidney renal clear cell carcinoma | Urogenital | Kidney | Proximal tubule | Epithelial | Clear-cell carcinoma | Clear-cell RCC | Mesoderm | Intermediate mesoderm (metanephric) | Both | Pan-kidney | C28 [2B] |
| **KIRP** | Kidney renal papillary | Urogenital | Kidney | Proximal tubule (papillary) | Epithelial | Papillary carcinoma | Papillary RCC | Mesoderm | Intermediate mesoderm (metanephric) | Both | Pan-kidney | C28 [2B] |
| **LGG** | Lower-grade glioma (WHO grades 2–3) | CNS | Brain/CNS | Glia (astrocyte/oligodendrocyte) | Non-epithelial | Glial / neuroepithelial | Astrocytoma + oligodendroglioma | Ectoderm | Neural ectoderm (neural tube) | Both | (CNS/endocrine) | C11 (IDH1-mut) / C23 (IDH1-wt; with GBM) [2B] |
| **LIHC** | Liver hepatocellular carcinoma | Digestive | Hepatobiliary + pancreas | Liver hepatocyte | Epithelial | Hepatocellular | HCC | Endoderm | Foregut endoderm | Both | Single-type dominant | C26 [2B] |
| **LUAD** | Lung adenocarcinoma | Thoracic | Lung/pleura | Lung peripheral epithelium | Epithelial | Adenocarcinoma | Lung adenoCa | Endoderm | Foregut endoderm (respiratory diverticulum) | Both | Single-type dominant | C14 [2B] |
| **LUSC** | Lung squamous | Thoracic | Lung/pleura | Lung central airway epithelium | Epithelial | Squamous | Lung SCC | Endoderm | Foregut endoderm | Both | Pan-squamous | C10 (Pan-SCC, predominantly LUSC) [2B] |
| **MESO** | Mesothelioma | Thoracic | Lung/pleura | Mesothelium (pleura/peritoneum) | Epithelial — *contested* † | Mesenchymal-mesothelial (boundary) | Pleural/peritoneal mesothelioma | Mesoderm | Lateral plate mesoderm (somatopleure) | Both | (uncommon, fringe) | Diverse [2B] |
| **OV** | Ovarian serous cystadenocarcinoma | Reproductive (female) | Gynecological | Ovary (or fallopian tube fimbria) ‡ | Epithelial | Serous papillary adenocarcinoma (mostly) | High-grade serous ovarian Ca | Mesoderm | Intermediate mesoderm (gonadal) + Müllerian (lateral plate) | Female | Pan-gyn | C6 (single-type dominant) [2B] |
| **PAAD** | Pancreatic ductal adenocarcinoma | Digestive | Hepatobiliary + pancreas | Pancreatic duct | Epithelial | Adenocarcinoma | PDAC | Endoderm | Foregut endoderm | Both | (Pan-GI fringe) | Diverse [2B] |
| **PCPG** | Pheochromocytoma + paraganglioma | Endocrine | Adrenal | Adrenal medulla / paraganglion | Non-epithelial | Neuroendocrine | Pheo / paraganglioma | Ectoderm | Neural crest | Both | (CNS/endocrine) | C5 [2B] |
| **PRAD** | Prostate adenocarcinoma | Urogenital | Male reproductive | Prostate | Epithelial | Adenocarcinoma | Prostate adenoCa | Endoderm + mesoderm — *contested* † | Urogenital sinus endoderm + mesoderm | Male | Single-type dominant | C16 [2B] |
| **SARC** | Sarcoma (multi-subtype) | Soft tissue | Soft tissue / mesenchymal | Multiple subtypes (LMS, DDLPS, MFS, UPS, MPNST, SS) | Non-epithelial | Mesenchymal | Soft tissue sarcoma | Mesoderm | Paraxial + lateral plate mesoderm (subtype-dependent) | Both | (Mesenchymal/immune) | C20 / C3 [2B] |
| **SKCM** | Skin cutaneous melanoma | Dermatological | Skin | Skin (melanocyte) | Non-epithelial | Melanocytic (neural-crest derivative) | Cutaneous melanoma | Ectoderm | Neural crest | Both | (Pan-skin/UV) | C15 (with UVM) [2B] |
| **STAD** | Stomach adenocarcinoma | Digestive | GI tract (upper) | Stomach (gastric epithelium) | Epithelial | Adenocarcinoma | Gastric adenoCa | Endoderm | Foregut endoderm | Both | Pan-GI | C1 (EBV-CIMP) / C18 (MSI) [2B] |
| **TGCT** | Testicular germ cell tumor | Reproductive (male) | Male reproductive | Testis (germ cell) | Non-epithelial | Germ cell | Seminoma + non-seminoma | Mesoderm — *germ cells originate from primordial germ cells* † | Intermediate mesoderm (gonad) + PGCs (extra-embryonic origin) | Male | (germ cell, fringe) | Diverse [2B] |
| **THCA** | Thyroid carcinoma | Endocrine | Thyroid | Thyroid follicular cell | Epithelial | Adenocarcinoma (follicular/papillary) | PTC + FTC | Endoderm | Foregut endoderm | Both | Single-type dominant | C12 [2B] |
| **THYM** | Thymoma + thymic carcinoma | Endocrine + lymphoid (mixed) | Thymus | Thymic epithelium | Epithelial | Thymic epithelial (sui generis) | Thymoma A/AB/B/C | Endoderm + mesoderm — *contested* † | Pharyngeal pouch endoderm (thymic epithelium) + mesoderm (stroma) | Both | (CNS/endocrine fringe) | Distinct cluster [2B] |
| **UCEC** | Uterine corpus endometrial carcinoma | Reproductive (female) | Gynecological | Endometrium | Epithelial | Endometrioid + Serous (split) | Endometrioid + Serous endometrial Ca | Mesoderm | Müllerian (lateral plate mesoderm) | Female | Pan-gyn | C8 (single-type dominant) [2B] |
| **UCS** | Uterine carcinosarcoma | Reproductive (female) | Gynecological | Uterine corpus (mixed carcinoma + sarcoma) | Mixed epithelial + non-epithelial — *biphasic* † | Carcinosarcoma (mixed) | Müllerian carcinosarcoma | Mesoderm | Müllerian (lateral plate mesoderm) | Female | (diverse, fringe) | Diverse [2B] |

**Legend**:
† = contested or ambiguous assignment (see "Contested or ambiguous assignments" section).
‡ = explanatory note in cohort/cell.
[2A], [2B], [2C] = Figure references in Hoadley et al. 2018, Cell 173:291-304.

---

## Hierarchy 1 — Tissue of origin (3 granularities)

### Granularity 1a — Coarsest (5 anatomical-system bins)

Five-bin top-level WHO-aligned groupings (the user's brief asks for 5, with the caveat that "if a different 5-category top-level is more standard"). The 5 used here are an honest collapse of WHO Classification of Tumors 5th-edition volume structure:

1. **Digestive** (foregut to hindgut, plus accessory organs): ESCA, STAD, CORE, LIHC, CHOL, PAAD
2. **Urogenital + Reproductive** (kidneys, bladder, both reproductive tracts): BLCA, KICH, KIRC, KIRP, PRAD, TGCT, BRCA, CESC, OV, UCEC, UCS
3. **Endocrine** (hormone-producing glands): ACC, PCPG, THCA, THYM
4. **Thoracic** (lung, pleura): LUAD, LUSC, MESO
5. **Other** (CNS, hematologic, skin, soft tissue, head & neck): LGG, DLBC, SKCM, SARC, HNSC

*Imbalance note*: This 5-bin schema produces cell counts of 6/11/4/3/5 — the "Other" bin and the urogenital+reproductive bin are the largest. There is no fully balanced 5-bin partition over 29 cancer types because TCGA's pan-cancer composition is itself biased toward urogenital and digestive tumors. WHO does not publish a 5-bin top-level; this is a methodologically imposed coarsening.

### Granularity 1b — Mid (organ-level, 11 bins)

| # | Bin | Members |
|---|---|---|
| 1 | Lung + pleura | LUAD, LUSC, MESO |
| 2 | Kidney | KICH, KIRC, KIRP |
| 3 | Hepatobiliary + pancreas | LIHC, CHOL, PAAD |
| 4 | GI tract | ESCA, STAD, CORE |
| 5 | Breast | BRCA |
| 6 | Gynecological (uterus + cervix + ovary) | OV, UCEC, UCS, CESC |
| 7 | Male reproductive + bladder | PRAD, TGCT, BLCA |
| 8 | Endocrine | ACC, PCPG, THCA, THYM |
| 9 | CNS | LGG |
| 10 | Skin | SKCM |
| 11 | Other (soft tissue, hematologic, head & neck) | SARC, DLBC, HNSC |

### Granularity 1c — Finest (specific tissue, 29 bins)

One bin per cancer type at this resolution; see master table column "Tissue: Finest". This is the "trivial" granularity — it provides no compression beyond the cancer-type label itself, but is included for completeness so that "tissue of origin at finest resolution" is well-defined as a hierarchy level.

---

## Hierarchy 2 — Epithelial classification (3 granularities)

### Granularity 2a — Coarsest (3 bins)

1. **Epithelial** (24 of 29): ACC, BLCA, BRCA, CESC, CHOL, CORE, ESCA, HNSC, KICH, KIRC, KIRP, LIHC, LUAD, LUSC, OV, PAAD, PRAD, STAD, THCA, THYM, UCEC, MESO ‡, UCS ‡
   - ‡ MESO and UCS are boundary cases — see contested section
2. **Non-epithelial** (4 of 29): LGG (glial/neuroepithelial), SARC (mesenchymal), SKCM (melanocytic — neural crest), TGCT (germ cell), PCPG (neuroendocrine — neural crest)
3. **Hematological** (1 of 29): DLBC

### Granularity 2b — Mid (6 bins by broad histology)

1. **Squamous**: HNSC, LUSC, CESC (mostly), ESCA (~50%) — *split for ESCA, CESC*
2. **Adenocarcinoma (glandular)**: BRCA, CORE, LUAD, PRAD, PAAD, STAD, UCEC, OV, ACC, THCA, LIHC ‡, CHOL ‡, ESCA (subset)
   - ‡ HCC and cholangiocarcinoma have hepatocellular and biliary subtypes that are technically adenocarcinoma-like
3. **Transitional / urothelial**: BLCA
4. **Glial / neuroepithelial**: LGG
5. **Mesenchymal / mesothelial / germ cell** (non-epithelial cancers): SARC, MESO, TGCT, UCS (mixed epi + mesench)
6. **Hematopoietic / neural-crest / sui generis**: DLBC, SKCM, PCPG, KICH ‡, KIRC ‡, KIRP ‡, THYM
   - ‡ Renal carcinomas don't fit cleanly into squamous/adeno/transitional; "RCC" is its own histology class

### Granularity 2c — Finest (10 histological-subtype bins)

1. Squamous cell carcinoma: HNSC, LUSC, CESC (subset), ESCA (subset)
2. Glandular adenocarcinoma: BRCA, CORE, LUAD, PRAD, PAAD, STAD, ACC, THCA
3. Hepatocellular carcinoma: LIHC
4. Cholangiocarcinoma: CHOL
5. Renal cell carcinoma (clear, papillary, chromophobe): KIRC, KIRP, KICH
6. Urothelial carcinoma: BLCA
7. Endometrioid + serous (gynecological adenoCa): UCEC, OV, CESC (subset)
8. Glial: LGG
9. Mesenchymal / mesothelial / sarcoma: SARC, MESO, UCS
10. Other (germ cell, melanocytic, neuroendocrine, hematopoietic, thymic): TGCT, SKCM, PCPG, DLBC, THYM

---

## Hierarchy 3 — Embryonic germ layer (2 granularities)

This is the **most contested** hierarchy. Many TCGA tumor types arise from anatomical sites with mixed embryonic origin (e.g., the bladder body is mesodermal but the trigone/urothelial lining has endodermal contribution from the cloaca). Compile-only assignments below cite developmental biology references where available.

### Granularity 3a — Coarsest (3 bins)

1. **Endoderm**: ESCA, STAD, CORE, LIHC, CHOL, PAAD, LUAD, LUSC, THCA (foregut/midgut/hindgut)
2. **Mesoderm**: KICH, KIRC, KIRP, OV, ACC, SARC, DLBC, MESO, UCEC, UCS, TGCT
3. **Ectoderm**: BRCA, SKCM, LGG, PCPG (surface + neural + neural crest)

**Contested 4 cells** (assigned with contested marker in master table):
- **BLCA**: trigone (endoderm-derived) + body (mesoderm-derived) — *contested*
- **CESC**: cervix has cloacal endoderm origin for ectocervical glands + Müllerian (lateral plate mesoderm) origin for endocervix — *contested*
- **HNSC**: oropharyngeal subsites range from ectoderm-derived (oral cavity) to endoderm-derived (deep pharynx) — *contested*
- **THYM**: thymic epithelium derives from pharyngeal pouch endoderm; thymic stroma derives from mesoderm — *contested*
- **PRAD**: prostate epithelium has endoderm + mesoderm origins — *contested*

### Granularity 3b — Mid (sub-divisions — 7 bins)

| # | Bin | Members |
|---|---|---|
| 1 | Foregut endoderm | ESCA, STAD, LIHC, CHOL, PAAD, LUAD, LUSC, THCA, THYM (epithelium) |
| 2 | Hindgut endoderm | CORE |
| 3 | Surface ectoderm | BRCA, SKCM (epidermal context) |
| 4 | Neural ectoderm (neural tube) | LGG |
| 5 | Neural crest | PCPG, SKCM (melanocytes specifically) |
| 6 | Intermediate mesoderm (urogenital ridge) | ACC, KICH, KIRC, KIRP, OV, TGCT |
| 7 | Lateral plate mesoderm (somatopleure / Müllerian) | MESO, UCEC, UCS |
| – | Hematopoietic mesoderm | DLBC |
| – | Paraxial / lateral plate mesoderm (subtype-dependent) | SARC |
| – | Mixed / contested | BLCA, CESC, HNSC, PRAD, THYM (see 3a) |

*Note*: SKCM has dual representation under both surface ectoderm (epidermal context) and neural crest (melanocytic lineage). Most authorities classify SKCM under neural crest because the cell-of-origin (melanocyte) is neural-crest-derived.

---

## Hierarchy 4 — Sex composition of cohort (1 granularity, ternary)

| Category | TCGA cancer types |
|---|---|
| **Both sexes** | ACC, BLCA, CHOL, CORE, DLBC, ESCA, HNSC, KICH, KIRC, KIRP, LGG, LIHC, LUAD, LUSC, MESO, PAAD, PCPG, SARC, SKCM, STAD, THCA, THYM |
| **Female-only** | CESC, OV, UCEC, UCS |
| **Male-only** | PRAD, TGCT |
| **Mixed but female-dominant (>95%)** | BRCA |

*Note on BRCA*: TCGA includes ~12 male BRCA samples out of ~1,100 total (~1%). BIDIFAC+ does not stratify; for the binary sex-restriction granularity, BRCA can be treated as either "both" or "female" depending on convention. Most analyses treat BRCA as functionally female-only.

---

## Hierarchy 5 — Hoadley 2018 cell-of-origin (mixed granularity)

Source: Hoadley et al. 2018, "Cell-of-Origin Patterns Dominate the Molecular Classification of 10,000 Tumors from 33 Types of Cancer," Cell 173(2):291-304. Figures 2A-D and Table S6.

**⚠ Circularity warning**: The Hoadley 2018 cluster assignments are **derived from TCGA itself** via integrated multi-omic clustering (iCluster on the same 4 omics platforms BIDIFAC+ uses). Using these assignments as a structural prior for a TCGA-based replication is methodologically circular — the prior contains exactly the data signal the model is trying to recover. This hierarchy is included for completeness because it's the most-cited cell-of-origin classification in pan-cancer methodology, but it is the LEAST defensible candidate as an external structural prior. Note this limitation prominently in any methodology paper.

### Granularity 5a — Coarsest (4 cell-of-origin patterns reported in Hoadley 2018)

The paper reports **4 major cell-of-origin patterns**, not 5 as suggested by the user's brief:

1. **Pan-gastrointestinal** (Pan-GI): CORE, STAD, ESCA (subset), CHOL (subset), PAAD (subset)
2. **Pan-gynecological** (Pan-gyn): BRCA, OV, UCEC, CESC (in part), UCS
3. **Pan-kidney**: KIRC, KIRP, KICH (with ACC nearby)
4. **Pan-squamous**: HNSC, LUSC, CESC (HPV+), BLCA (subset), ESCA (subset)

*Cancer types not assigned to a pan-cluster*: ACC, DLBC, LGG, LIHC, LUAD, MESO, PCPG, PRAD, SARC, SKCM, TGCT, THCA, THYM — these form single-type-dominant or mixed clusters.

### Granularity 5b — Mid (no canonical 14-15 cluster intermediate)

The user's brief mentioned ~14-15 clusters at mid-level, but Hoadley 2018 does not report such a level. The paper goes from 4 cell-of-origin patterns directly to 28 iClusters. **There is no canonical mid-level granularity in Hoadley 2018.** Any 14-15 grouping would be a post-hoc collapse and would need its own methodology citation. Flagging as "not in Hoadley 2018."

### Granularity 5c — Finest (28 integrated clusters)

Per Figure 2B / Table S6:

| iCluster | Composition |
|---|---|
| C1 | STAD (EBV-CIMP) |
| C2 | BRCA (HER2-amp / ERBB2-amplified) |
| C3 | Mesenchymal/immune (multi-type) |
| C4 | Pan-GI / CRC (predominantly CORE/COAD/READ) |
| C5 | CNS / endocrine (PCPG-related) |
| C6 | OV (single-type dominant) |
| C7 | (mixed) |
| C8 | UCEC (single-type dominant) |
| C9 | KICH + ACC |
| C10 | Pan-SCC (predominantly LUSC) |
| C11 | LGG (IDH1-mutant) |
| C12 | THCA (single-type dominant) |
| C13 | Mixed (chr8 del); BRCA subset |
| C14 | LUAD (single-type dominant) |
| C15 | SKCM + UVM (melanomas) |
| C16 | PRAD (single-type dominant) |
| C17 | (smaller cluster) |
| C18 | Pan-GI MSI (STAD/COAD MSI-high) |
| C19 | BRCA (luminal) |
| C20 | Mixed (stromal/immune; 25 tumor types) |
| C21 | (smaller cluster) |
| C22 | (smaller cluster) |
| C23 | GBM + LGG (IDH1-wt) ‡ |
| C24 | LAML ‡ |
| C25 | Pan-SCC (Chr11 amplification) |
| C26 | LIHC (single-type dominant) |
| C27 | Pan-SCC HPV (CESC, HNSC HPV+, BLCA subset) |
| C28 | Pan-kidney (KIRC + KIRP) |

‡ C23 contains GBM (not in BIDIFAC+'s 29) and LGG. C24 contains LAML (not in BIDIFAC+'s 29). Six cancer types reported in Hoadley as having <50% representation in any single cluster (i.e., genuinely diverse): **BLCA, UCS, HNSC, ESCA, STAD, CHOL** (Figure 2B + accompanying text).

---

## Sources and citations

| Hierarchy / granularity | Canonical reference |
|---|---|
| Tissue of origin (all 3 granularities) | WHO Classification of Tumours, 5th edition (2019–ongoing); the multi-volume series organizes tumors by anatomical region and provides histological subtype taxonomy. The 5-bin coarsest is a methodological collapse — WHO does not publish a 5-bin top level. |
| Epithelial classification (all 3 granularities) | WHO Classification of Tumours 5e + Robbins & Cotran Pathologic Basis of Disease, 10th edition (Kumar, Abbas, Aster, 2020). |
| Embryonic germ layer | Gilbert SF, Barresi MJF. *Developmental Biology*, 12th edition (Sinauer/OUP 2019). Sadler TW. *Langman's Medical Embryology*, 14th edition (Wolters Kluwer 2019). |
| Sex composition | TCGA cohort metadata (cBioPortal study summaries; TCGA Network publications per cancer type). |
| Hoadley 2018 cell-of-origin (4 patterns + 28 iClusters) | Hoadley KA, Yau C, Hinoue T, et al. *Cell* 2018;173(2):291-304. Figures 2A–D, Table S6. PMC: [PMC5957518](https://pmc.ncbi.nlm.nih.gov/articles/PMC5957518/). |

---

## Contested or ambiguous assignments

1. **BLCA — germ layer (endoderm vs mesoderm)**. Bladder trigone derives from cloacal endoderm; bladder body derives from intraembryonic mesoderm. Most modern developmental references (Sadler 2019) assign the bulk of bladder urothelium to endoderm, with the muscularis propria mesodermal. *Clinically*, bladder cancer is treated as a urothelial / transitional malignancy with mixed origins. **Recommendation: mark contested; do not assign single germ-layer column without footnote.**

2. **BRCA — germ layer (ectoderm with mesodermal stroma)**. Mammary gland epithelium is surface-ectoderm-derived; mammary stroma is mesodermal. The cancer arises in the epithelial compartment, so single-germ-layer assignment to ectoderm is defensible but loses the stromal mesodermal contribution that's known to influence tumor biology. **Recommendation: ectoderm at coarsest granularity; flag stromal contribution.**

3. **CESC — germ layer (endoderm + mesoderm split)**. Cervix has cloacal endoderm (ectocervix) + Müllerian/lateral plate mesoderm (endocervix) origins. CESC tumors include both squamous (ectocervical) and adenocarcinoma (endocervical) histologies. **Recommendation: mark contested; the ECCO/CESC adenocarcinoma split itself is a sub-classification.**

4. **HNSC — germ layer (ectoderm + endoderm by anatomical subsite)**. Oral cavity (anterior 2/3 of tongue, floor of mouth) is surface-ectoderm-derived; oropharynx + posterior tongue + larynx is foregut-endoderm-derived. **Recommendation: mark contested; HNSC is genuinely heterogeneous by subsite at the germ-layer level.**

5. **KIRC, KIRP, KICH — RCC histological classification at finest level**. These are usually grouped under "renal cell carcinoma" but are histologically distinct — clear cell vs papillary vs chromophobe arise from different parts of the nephron. The "Pan-kidney" Hoadley cluster groups KIRC and KIRP but assigns KICH to C9 (with ACC), not the Pan-kidney cluster. **Recommendation: at finest histology granularity, distinct subtypes; at organ granularity, all three are kidney.**

6. **LGG — cell-of-origin (multiple competing classifications)**. Lower-grade gliomas have IDH1-mutant and IDH1-wildtype subtypes that have very different biology. WHO 2021 CNS tumor classification (5th ed) groups them under "diffuse glioma" but separates by molecular subtype. Hoadley 2018 places IDH1-mut LGG in C11 and IDH1-wt LGG in C23 (with GBM). **Recommendation: mark contested at coarsest histology; IDH1 status should be a separate column if relevant.**

7. **MESO — epithelial status (mesothelial — boundary case)**. Mesothelial cells are epithelial-like (express keratin, form sheets) but derived from mesoderm. WHO classifies mesothelioma under "mesothelial tumors" in the volume on lung/pleura/thymus/heart, and the histology has both epithelioid and sarcomatoid variants. **Recommendation: mark contested at epithelial-vs-non-epithelial granularity; at finest histology, "mesothelial" is its own category.**

8. **PRAD — germ layer (endoderm + mesoderm)**. Prostate develops from urogenital sinus endoderm (epithelium) + mesenchyme (stroma). Most authorities assign prostatic epithelium to endoderm. **Recommendation: endoderm at coarsest with mesodermal stromal note.**

9. **TGCT — germ cell origin**. Testicular germ cell tumors arise from primordial germ cells (PGCs), which are extra-embryonic in origin (yolk sac → genital ridge migration). At a strict embryological-germ-layer level, PGCs precede the tripartite germ-layer establishment. *Practically*, TGCT is most often assigned to mesoderm because PGCs colonize the gonad (intermediate mesoderm) and acquire mesodermal context. **Recommendation: mesoderm with footnote; PGC origin is technically pre-germ-layer.**

10. **THYM — germ layer (endoderm + mesoderm)**. Thymic epithelium derives from 3rd pharyngeal pouch endoderm; thymic stroma is mesoderm-derived. Thymoma is a tumor of the epithelial component, but THYM tumors include subtypes with varying epithelial vs stromal involvement. **Recommendation: endoderm for thymic epithelial origin; flag dual contribution.**

11. **UCS — biphasic (carcinoma + sarcoma)**. Uterine carcinosarcoma is a biphasic tumor with both carcinomatous (epithelial) and sarcomatous (mesenchymal) components. The current consensus is that UCS is a metaplastic carcinoma (clonally epithelial, with sarcoma-like differentiation) rather than a true mixed tumor — but classification has been debated. WHO 5e classifies UCS under "uterine sarcomas" but acknowledges the metaplastic origin debate. **Recommendation: mark contested at epithelial-vs-non-epithelial; the dual nature is itself the data.**

12. **THYM — Hoadley iCluster assignment**. Hoadley 2018 places THYM in a fringe cluster (not one of the named pan-clusters). Specific iCluster assignment varies between paper figures — Table S6 should be authoritative. **Recommendation: cite Table S6 directly for THYM.**

---

## TCGA cohort notes

1. **CORE = COAD + READ merged**. BIDIFAC+ treats colon adenocarcinoma (COAD) and rectum adenocarcinoma (READ) as a single unit "CORE" because they are biologically near-indistinguishable at the molecular level (Hoadley 2018 confirms; both fall in the same Pan-GI iClusters C4 / C18). TCGA elsewhere keeps them separate, so cross-paper comparisons need this re-mapping.

2. **GBM, LAML, UVM excluded from BIDIFAC+'s 29**. These three appear in Hoadley 2018's 33-type cohort but not in BIDIFAC+. Likely reasons:
   - **GBM** (glioblastoma): may have been merged into "diffuse glioma" with LGG, or excluded because methylation panel coverage differed.
   - **LAML** (acute myeloid leukemia): hematologic — many TCGA solid-tumor analyses exclude leukemia for its different sample-collection pipeline (peripheral blood) and non-comparable mutation calling.
   - **UVM** (uveal melanoma): small cohort (~80 samples), often pooled with cutaneous melanoma in pan-cancer studies but pulled out by Hoadley.

3. **BIDIFAC+'s 29-type list is documented in Lock/Park/Hoadley AOAS 2022**. The exact list with sample sizes per type is in their supplementary materials (Lock 2022 Annals of Applied Statistics, Vol. 16 Issue 1).

4. **Sample-size imbalance across the 29 types**. From TCGA (post-QC sample counts available in cBioPortal):
   - Largest: BRCA (~1,100), KIRC (~600), LUAD (~570), HNSC (~520), LGG (~520), UCEC (~550)
   - Smallest: CHOL (~36), DLBC (~48), MESO (~87), UCS (~57), KICH (~66)
   - The 10× sample-size range across types is itself a hierarchy candidate (large vs small cohorts may behave differently in BIDIFAC+'s shared/specific decomposition).

---

## Additional candidate hierarchies (flagged but NOT in main table)

Per the user's brief: don't pad, but flag if obviously relevant.

1. **Mutation burden tier** (low/medium/high TMB). Three rough bins from the TCGA Pan-Cancer Atlas Mutational Heterogeneity paper (Lawrence et al. 2013, Nature; Bailey et al. 2018, Cell). Low: THCA, PCPG, TGCT, ACC, LGG, KICH; High: SKCM, LUAD, LUSC, BLCA, CESC (HPV+), HNSC, STAD, CORE; Medium: rest. **Relevant because** TMB clustering tracks closely with environmental-mutagen exposure (UV, smoking, MSI), which is itself a candidate structural prior.

2. **Differentiation state / stem-vs-differentiated**. Hoadley 2018 implicitly captures this in their "stemness signature" analysis but doesn't bin it as a separate hierarchy. Gentles et al. 2015 (Nat Med) provides a defensible pan-cancer stemness classification. **Relevant because** stemness signatures are themselves NMF latent factors with parameters scaling per cancer type.

3. **Microsatellite stability tier (MSI-H / MSS)**. Binary or three-bin (MSI-H / MSS / MSI-L). Most relevant for STAD and CORE; smaller minorities in UCEC, others. Tracks with mutation burden but is methodologically distinct.

4. **Viral etiology**. Three bins: HPV-driven (CESC, HNSC subset), HBV/HCV-driven (LIHC subset), EBV-driven (STAD subset). Methodologically interesting because viral integration is a discrete latent state.

5. **Anatomical laterality** (paired vs unpaired organs). KIRC/KIRP/KICH/OV/TGCT/PCPG are bilateral; others are not. Probably methodologically irrelevant but flagged for completeness.

These five additional candidates are not in the main reference table per the user's "compile, don't synthesize, don't pad" directive, but should be considered if the main 5 hierarchies don't yield a clean structural prior recovery.

---

*Compiled 2026-05-09 as input to the methodology contribution to cancer genomics; pairs with the BIDIFAC+ paper for structural-candidate identification work. No analysis performed; reference document only.*
