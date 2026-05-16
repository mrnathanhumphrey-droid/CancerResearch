# Paper 5 — Data Access

## Primary substrate: Cavalli et al. 2017

| accession | platform | n | content | URL |
|---|---|---|---|---|
| **GSE85217** | Affymetrix Human Gene 1.1 ST | 763 | gene expression, RMA-normalized | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE85217 |
| **GSE85212** | Illumina HumanMethylation450 | 763 | DNA methylation β-values | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE85212 |

GEO samples are joined to clinical data via `Sample_title` (GEO pData `title`) ↔ `Study_ID` (mmc2 column 0), both formatted as `MB_SubtypeStudy_NNNNN`. Verified 2026-05-15 against mmc2 row inspection.

## Supplementary tables (placed in this folder 2026-05-15)

Downloaded from Elsevier CDN, PII `S1535610817302015` (Cavalli 2017 *Cancer Cell* DOI 10.1016/j.ccell.2017.05.005):

| file | source URL | size | content |
|---|---|---|---|
| `cavalli_supplementary_S1.xlsx` | `https://ars.els-cdn.com/content/image/1-s2.0-S1535610817302015-mmc2.xlsx` | 217 KB | **Table S1** (clinical + subgroup + 12-subtype + OS): 1 sheet "Sheet1", 765 rows × 57 cols, header on row 2, 763 sample rows. Columns: Study_ID, Age, AgeGroup, Gender, histology, Met status, Dead, OS (years), Subgroup, Subtype, c_2..c_8, SNP6_data, plus 40 per-chromosome-arm CNA event columns. |
| `cavalli_supplementary_S2_subgroup_specific.xlsx` | `https://ars.els-cdn.com/content/image/1-s2.0-S1535610817302015-mmc3.xlsx` | 35 KB | Subgroup-specific results: 5 sheets (Subgroup, WNT, SHH, Group 3, Group 4) |
| `cavalli_supplementary_S3_topgenes.xlsx` | `https://ars.els-cdn.com/content/image/1-s2.0-S1535610817302015-mmc4.xlsx` | 1.5 MB | Top genes (expression + methylation) per subgroup |
| `cavalli_supplementary_S4_subtype_specific.xlsx` | `https://ars.els-cdn.com/content/image/1-s2.0-S1535610817302015-mmc5.xlsx` | 656 KB | Subtype-specific results: 12 sheets, one per 12-subtype (WNT alpha/beta, SHH alpha/beta/gamma/delta, Group 3 alpha/beta/gamma, Group 4 alpha/beta/gamma) |

The `fetch_cavalli_data.R` script reads only `cavalli_supplementary_S1.xlsx` for the clinical pipeline; the other three are shipped here for completeness and downstream sensitivity analyses.

## Expected on-disk artifacts after `fetch_cavalli_data.R`

```
data/
├── cavalli_supplementary_S1.xlsx                     ← Table S1 (placed 2026-05-15)
├── cavalli_supplementary_S2_subgroup_specific.xlsx   ← Table S2
├── cavalli_supplementary_S3_topgenes.xlsx            ← Table S3
├── cavalli_supplementary_S4_subtype_specific.xlsx    ← Table S4
├── cavalli_expr.rds                                  ← from GSE85217 (probes × 763)
├── cavalli_expr_pdata.rds                            ← GEO pData (763 × ~35 cols)
├── cavalli_meth.rds                                  ← from GSE85212 (CpG × 763)
├── cavalli_meth_pdata.rds                            ← GEO pData (763 × ~35 cols)
├── cavalli_clinical.rds                              ← parsed mmc2 (763 × 10 cols)
└── cavalli_supp_sha256.txt                           ← provenance hash
```

After `run_bidifac_plus.R`, additional artifacts:

```
data/
├── cavalli_clinical_aligned.rds        ← clinical rows reordered to match GEO sample order
├── bidifac_components.rds              ← per-sample scores (samples × n_components)
├── bidifac_components_shared.rds       ← shared components only
├── bidifac_components_per_block.rds    ← per-block components only
├── bidifac_loadings.rds                ← block-level loadings
├── bidifac_diagnostics.rds             ← rank-selection trace
├── train_data.rda                      ← stratified-80 Lock-style train
├── test_data.rda                       ← stratified-20 Lock-style test
├── target_covs_primary.rda             ← K=3 target covariates, L1-pooled at L2 granularity
└── target_covs_sensitivity.rda         ← K=3 target covariates, L1 4-subgroup
```

## Deferred (external validation phase)

- **Sharma et al. 2019** (n = 852 with matched transcriptome of the 1501 methylomes) — data deposit not fully resolved at pre-reg time; either confirmed via the paper's "Data availability" section or contacted via corresponding author.
- **MAGIC / St. Jude Integrative Portal** (n = 898 across ACNS0331 + SJMB03 + ACNS0332) — dbGaP-controlled raw data; portal access at https://viz.stjude.cloud/st-jude-childrens-research-hospital/visualization/medulloblastoma-integrative-analysis-portal~2882 for clinical + methylation-derived subgroup labels.

Neither is fetched by Paper 5 primary.
