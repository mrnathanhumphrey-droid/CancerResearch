# Finer-Resolution Prediction Test for Cov 8.2
## 2026-05-10 — discovery work, ANOVA-only, no Gibbs

This is a **prediction test**, not a validation. The methodology validated in
Paper 1 (Δ=+6.26 nats, CI excludes zero) used a coarse 3-group structural
prior on cov 8.2 (epithelial classification at 2a_3bin). The methodology
*predicts* that finer biological resolution should organize residual
variation that the coarse prior didn't capture. This file tests that
prediction at two resolutions:

- **PRIMARY (Option A)** — tissue-within-group: ANOVA cov-8.2 residuals
  (per-cancer beta minus epi/2a group mean) against tissue_origin/1b_11bin
  labels, computed within each of the 3 epi/2a groups.
- **SENSITIVITY (Option C)** — hybrid: tissues with >1 TCGA cancer use
  modal molecular subtype as label; single-cancer tissues use the tissue
  label.

**This file does not change Paper 1's validated narrow result.** The
test informs whether Paper 1's discussion (or a future Paper 2 framing)
should claim methodology-resolution beyond the validated 3-group level.

---

## 1. Pre-committed assignments

Subtype assignments locked **before computing**, sourced from canonical
TCGA flagship papers per cancer. When modal subtype is unclear or no
canonical TCGA subtype exists at the cohort size, the cancer name itself
is used as the subtype label (effectively makes that cancer its own
subtype-bin). The locked list:

| cancer | tissue (1b_11bin) | modal subtype | source |
|---|---|---|---|
| ACC | Endocrine | ACC (cancer name; CIMP cluster modal uncertain) | TCGA ACC |
| BLCA | Male repro+bladder | Luminal | TCGA BLCA |
| BRCA | Breast | LumA | PAM50 |
| CESC | Gynecological | Squamous_HPV+ | TCGA CESC |
| CHOL | Hepatobiliary+pancreas | CHOL (cancer name) | no canonical subtype |
| CORE | GI tract | CMS2 | Consensus Molecular Subtypes |
| DLBC | Soft+heme+H&N | GCB | TCGA DLBC |
| ESCA | GI tract | ESCC1 | TCGA ESCA squamous-majority |
| HNSC | Soft+heme+H&N | Classical | TCGA HNSC HPV-neg modal |
| KICH | Kidney | Chromophobe | only one chromophobe subtype |
| KIRC | Kidney | ccA | ClearCode34 |
| KIRP | Kidney | Type1 | TCGA KIRP papillary |
| LGG | CNS | IDH1mut | TCGA LGG |
| LIHC | Hepatobiliary+pancreas | LIHC (cancer name) | iCluster heterogeneous |
| LUAD | Lung+pleura | TRU | TCGA LUAD |
| LUSC | Lung+pleura | Classical | Wilkerson 2010 |
| OV | Gynecological | Mesenchymal | Verhaak/Tothill |
| PAAD | Hepatobiliary+pancreas | Classical | Bailey/Moffitt |
| PCPG | Endocrine | Pseudohypoxia | TCGA PCPG cluster 1 |
| PRAD | Male repro+bladder | ERG_fusion+ | TCGA PRAD |
| SARC | Soft+heme+H&N | LMS | TCGA SARC |
| SKCM | Skin | BRAF_mut | TCGA SKCM |
| STAD | GI tract | CIN | TCGA STAD |
| TGCT | Male repro+bladder | Seminoma | TCGA TGCT |
| THCA | Endocrine | BRAF_papillary | TCGA THCA |
| THYM | Endocrine | B2 | WHO type B2 modal |
| UCEC | Gynecological | Endometrioid | TCGA UCEC |
| MESO, UCS | — | (contested for epi/2a, dropped) | — |

The hybrid label (Option C) is built per cancer:
- If tissue contains >1 TCGA cancer (Lung+pleura, Kidney, Hepatobiliary+pancreas,
  GI tract, Gynecological, Male repro+bladder, Endocrine, Soft+heme+H&N) →
  use modal subtype.
- If tissue contains 1 TCGA cancer (Breast, CNS, Skin) → use tissue label.

---

## 2. Residual setup

26 cancers retained (29 minus contested MESO/UCS at epi/2a; minus 1
because cov 8.2 isn't measured for one of the cancers — DLBC is in
Hematological group as the only member).

Per-group means at epi/2a (cov 8.2 baseline beta):
- **Epithelial** (n=21): mean β = +0.0020
- **Hematological** (n=1): β = −0.128 (DLBC alone; residual = 0)
- **Non-epithelial** (n=4): mean β = +0.0037

Residuals for each cancer (= β − group mean):

| cancer | group | residual | tissue | hybrid |
|---|---|---|---|---|
| ACC | Epithelial | +0.0213 | Endocrine | ACC |
| THYM | Epithelial | +0.0108 | Endocrine | B2 |
| CHOL | Epithelial | +0.0139 | Hepatobiliary+pancreas | CHOL |
| PRAD | Epithelial | +0.0039 | Male repro+bladder | ERG_fusion+ |
| THCA | Epithelial | +0.0011 | Endocrine | BRAF_papillary |
| ESCA | Epithelial | +0.0007 | GI tract | ESCC1 |
| KIRC | Epithelial | +0.0000 | Kidney | ccA |
| PAAD | Epithelial | +0.0001 | Hepatobiliary+pancreas | Classical |
| LIHC | Epithelial | −0.0010 | Hepatobiliary+pancreas | LIHC |
| LUAD | Epithelial | −0.0007 | Lung+pleura | TRU |
| STAD | Epithelial | −0.0013 | GI tract | CIN |
| CORE | Epithelial | −0.0015 | GI tract | CMS2 |
| HNSC | Epithelial | −0.0021 | Soft+heme+H&N | Classical |
| UCEC | Epithelial | −0.0031 | Gynecological | Endometrioid |
| CESC | Epithelial | −0.0036 | Gynecological | Squamous_HPV+ |
| BLCA | Epithelial | −0.0037 | Male repro+bladder | Luminal |
| KIRP | Epithelial | −0.0043 | Kidney | Type1 |
| LUSC | Epithelial | −0.0050 | Lung+pleura | Classical |
| OV | Epithelial | −0.0066 | Gynecological | Mesenchymal |
| KICH | Epithelial | −0.0087 | Kidney | Chromophobe |
| BRCA | Epithelial | −0.0101 | Breast | Breast |
| DLBC | Hematological | 0.000 | Soft+heme+H&N | GCB |
| SARC | Non-epithelial | +0.0030 | Soft+heme+H&N | LMS |
| SKCM | Non-epithelial | +0.0032 | Skin | Skin |
| PCPG | Non-epithelial | −0.0023 | Endocrine | Pseudohypoxia |
| TGCT | Non-epithelial | −0.0039 | Male repro+bladder | Seminoma |

---

## 3. PRIMARY (Option A): tissue within group

| epi/2a group | n cancers | n tissues | n singleton | n multi | η² obs | perm null median | perm p95 | p_emp | underpowered? |
|---|---|---|---|---|---|---|---|---|---|
| **Epithelial** | **21** | **9** | **2** | **7** | **0.608** | **0.381** | **0.685** | **0.117** | **no** |
| Hematological | 1 | 1 | 1 | 0 | — | — | — | — | YES (n<4) |
| Non-epithelial | 4 | 4 | 4 | 0 | 1.000 | 1.000 | 1.000 | 1.000 | YES (all singleton) |

Only the Epithelial group is testable. Hematological is degenerate (n=1).
Non-epithelial is degenerate (all 4 cancers in different tissues —
SARC/Soft+heme, SKCM/Skin, PCPG/Endocrine, TGCT/Male repro+bladder — so η²=1
trivially with no within-group variance to estimate).

**Epithelial result:** observed η² = 0.608. Permutation null median = 0.381,
95th pct = 0.685. Observed η² is between the null median and 95th pct.
**p_emp = 0.117** (one-sided permutation p). Not significant at α=0.05.

Standard ANOVA F-test: F = 2.328, p = 0.091 (similar conclusion).

The bin-imbalance among the 21 Epithelial cancers across 9 tissues
(2 singletons: Breast/BRCA only; Skin not in epi/2a Epithelial group) is
producing 38% of the observed η² under permutation alone. The actual
PIP-vs-tissue signal in the residuals adds another ≈ 23 percentage points
of η² over chance, but this excess is not statistically certified at n=21.

### Tissue-coherent pattern in Epithelial residuals (qualitative)

Residual sums per tissue within Epithelial group:

| tissue | n | sum residuals | mean residual | direction |
|---|---|---|---|---|
| Endocrine | 3 (ACC, THCA, THYM) | +0.0332 | +0.0111 | **all positive** |
| Hepatobiliary+pancreas | 3 (CHOL, LIHC, PAAD) | +0.0130 | +0.0043 | mixed; CHOL drives |
| Male repro+bladder | 2 (BLCA, PRAD) | +0.0002 | +0.0001 | mixed |
| GI tract | 3 (CORE, ESCA, STAD) | −0.0021 | −0.0007 | nearly zero |
| Soft+heme+H&N | 1 (HNSC) | −0.0021 | — | — |
| Lung+pleura | 2 (LUAD, LUSC) | −0.0057 | −0.0029 | both negative |
| Gynecological | 3 (CESC, OV, UCEC) | −0.0133 | −0.0044 | **all negative** |
| Kidney | 3 (KICH, KIRC, KIRP) | −0.0130 | −0.0043 | KICH most negative; mixed |
| Breast | 1 (BRCA) | −0.0101 | — | — |

There is a directional pattern: Endocrine cancers cluster positive,
Gynecological cluster negative, Lung negative, Kidney mostly negative,
Breast (singleton) negative, Hepatobiliary mixed driven by CHOL. The
pattern is qualitatively consistent with tissue-organized residual
variation — but at n=21 with effect sizes < 0.022 in absolute terms,
the permutation test cannot distinguish the pattern from chance at
α=0.05 (p_emp = 0.117).

---

## 4. SENSITIVITY (Option C): hybrid label

| epi/2a group | n cancers | n labels | n singleton | n multi | η² obs | perm null median | perm p95 | p_emp | underpowered? |
|---|---|---|---|---|---|---|---|---|---|
| Epithelial | 21 | 19 | **18** | **1** | 0.988 | 0.953 | 0.998 | 0.232 | **YES** (>50% singleton) |
| Hematological | 1 | 1 | 1 | 0 | — | — | — | — | YES |
| Non-epithelial | 4 | 4 | 4 | 0 | 1.000 | 1.000 | 1.000 | 1.000 | YES |

Within Epithelial group, the only non-singleton hybrid label is "Classical"
(LUSC + HNSC + PAAD). All 18 other cancers are in singleton hybrid bins.
Per the pre-registered degeneracy rule, this is **underpowered** (>50%
singletons).

η² = 0.988 is mathematically dominated by the singleton structure —
the inflated observed η² is mostly just "every cancer in own bin."
Permutation null median = 0.953 confirms this: even under PIP-randomization,
the bin structure forces η² near 1. The observed-vs-null gap is 0.04,
within the random fluctuation band (p_emp = 0.232).

**Option C does not add interpretable resolution at this n.** The
finer-than-tissue partition produces too many singleton bins for the
ANOVA to test against. This isn't evidence against the prediction;
it's evidence that the dataset is too small to certify
finer-than-tissue resolution given the constraint of subtype-modal labels.

---

## 5. Option A vs Option C

| epi/2a group | η² Option A | η² Option C | Δ (C − A) | underpowered? |
|---|---|---|---|---|
| Epithelial | 0.608 | 0.988 | +0.380 | **C only** |
| Hematological | — | — | — | both |
| Non-epithelial | 1.000 | 1.000 | 0.000 | both |

The +0.38 η² gap between Option C and Option A in the Epithelial group is
**not** evidence for added explanatory power from molecular subtype — it
is the direct artifact of going from 9 tissue bins (with 7 multi-cancer
bins, providing within-bin variance) to 19 hybrid labels (18 singletons,
providing no within-bin variance). The η² inflation from singleton-bin
saturation is large in absolute terms but not informative.

**At n=21 in the Epithelial group, the test cannot distinguish Option C from
Option A at the methodology-resolution level.** The pre-registered
disposition holds: methodology certified at the resolution Option A tests
(coarse epi/2a × within-group tissue, with p_emp=0.117 in Epithelial);
finer resolutions remain predictions awaiting larger n.

---

## 6. Disposition

Per the pre-registered interpretation guidance:

> "If primary η² > 0.30 AND p_emp < 0.05 in at least 2 of 3 groups: the
> methodology predicts at tissue-within-group resolution."

**Not met.** Only 1 of 3 groups testable; that one has η²=0.608 with
p_emp=0.117 (not below 0.05).

> "If primary η² is at chance (within permutation null) across all 3
> groups: the methodology operates at the coarse resolution only."

**Not met either.** The Epithelial group's observed η²=0.608 is above the
null median (0.381) — directionally consistent with the prediction, but
not statistically certified.

> "If primary η² is high in some groups but at chance in others: the
> methodology predicts at finer resolution within some epithelial groups
> but not others."

**Closest fit, but qualified by underpowered cells.** The Epithelial
group shows η²=0.608 with directional consistency (Endocrine cluster
positive, Gynecological cluster negative, Lung negative). The other two
groups are degenerate at this n.

### Honest reading

At n=21–26 (the testable Epithelial group), the finer-resolution
prediction is **directionally consistent** with the data (η²=0.608,
above null median 0.381) but **not certified at conventional rigor**
(p_emp=0.117). The non-Epithelial / Hematological groups are degenerate
at this n. The hybrid (subtype-when-possible) resolution is past the
data's certifiable resolution at this sample size.

This result:

- Does **not** support claiming methodology operates at tissue-within-group
  resolution as a finding. The data doesn't certify it at α=0.05.
- Does **not** support claiming a resolution ceiling at the coarse epi/2a
  level. The data is directionally consistent with finer resolution.
- **Does** support flagging "tissue-within-group resolution shows
  directional consistency in the testable group but is not statistically
  certified at this n" as a candidate Paper 2 (or Paper 1 discussion-section)
  follow-up. Larger pan-cancer cohorts (e.g., HTAN, AACR GENIE) would
  resolve this question.

The validated narrow result from Paper 1 stands either way: Paper 1's
3-covariate (or 4-covariate) structural prior at coarse epi/2a granularity
outperforms the unconstrained baseline on held-out predictive
performance.

---

## 7. Files

| file | contents |
|---|---|
| `scripts/extension/02_finer_resolution_cov8_2.py` | reproducer (ANOVA + permutation null) |
| `results/extension/finer_resolution_residuals.csv` | 26 cancer × {epi_2a, tissue, hybrid, beta_mean, group_mean, residual} |
| `results/extension/option_a_anova.csv` | per-group ANOVA + permutation null for Option A |
| `results/extension/option_c_anova.csv` | per-group ANOVA + permutation null for Option C |
| `results/extension/option_a_vs_c_comparison.csv` | side-by-side η² comparison |
