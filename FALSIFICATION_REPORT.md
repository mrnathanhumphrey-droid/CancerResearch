# Falsification Report — Structural-Hierarchy Screening
## 2026-05-09 — adversarial stress-tests of yesterday's screening verdict

**Yesterday's headline (now retracted):** "30 covariates have at least one
non-Hoadley structural hit at η² > 0.30; tissue_origin at organ-level shows
granularity-robust ρ=0.72; the methodology has empirical traction."

**This report's verdict:** Under multiple-testing correction, block-presence
circularity check, granularity-nesting null, and cohort-size partialing,
**no candidate survives the adversarial battery at conventional rigor.** The
screening produces η² values at the chance rate. The hypothesis cannot be
rejected outright (large effect sizes on a few candidates remain), but it
also cannot be supported at the level required for publication.

---

## Test A. Permutation null on every (cov, hier, gran) screening row

For each test, hold the cancer-set and bin-assignment fixed; shuffle which
PIP belongs to which cancer (B=1000 per test). Recompute η². Empirical
p = P(η²_perm ≥ η²_obs).

**Result (319 valid tests):**

| group | n | p<0.05 | expected | p<0.01 | expected |
|---|---|---|---|---|---|
| ALL | 319 | **18** | 16.0 | 3 | 3.2 |
| Non-Hoadley | 259 | 15 | 13.0 | 2 | 2.6 |
| Hoadley | 60 | 3 | 3.0 | 1 | 0.6 |

Mean observed η² = **0.388**. Mean permutation-null η² (median over B perms) = **0.363**.

**Interpretation:** The bin-imbalance among cancers in each test produces 93%
of the observed η² on average. Survivor count is exactly at the chance rate.
The screening's apparent signal is mostly artifact of which cancers happen
to have which covariate present, intersected with which cancers happen to
fall in which bins.

**File:** `falsification_perm_null.csv`

---

## Test B. BIDIFAC+ block-presence confound

Most BIDIFAC+ blocks are tissue-restricted (e.g., block 6 is in 1 cancer,
block 11 is in 4, block 13 is in 8). If block presence is itself
tissue-clustered, the per-cancer PIP profiles inherit tissue structure
through BIDIFAC+ — circular.

For each (block, hierarchy, granularity), I computed:
- block_presence_η²: ANOVA of the binary block-presence vector vs hierarchy bins
- avg_pip_screening_η²: mean of screening η² values for covariates from this block at this (hier, gran)
- Spearman ρ between the two over all blocks

**Result:**

| hierarchy | ρ(block_presence_η², pip_screening_η²) | n cells with both >0.30 |
|---|---|---|
| **hoadley** | **+0.597** | 23 |
| epithelial_class | +0.349 | 20 |
| germ_layer | +0.231 | 11 |
| sex_composition | +0.183 | 1 |
| **tissue_origin** | **−0.075** | 21 |

**Interpretation:**
- Hoadley's signal is **substantially confounded** by block-presence circularity (ρ=+0.60, ~36% shared variance with what BIDIFAC+ encoded a priori).
- Epithelial_class is partially confounded (ρ=+0.35).
- **Tissue_origin is CLEAN of this confound** (ρ=−0.075). If there's tissue-origin signal in PIPs, it isn't BIDIFAC+'s own clustering bleeding through. (But see Test D below — most "tissue-origin signal" turns out to fail other checks.)

Many high block-presence-η² cells are blocks present in only 1–4 cancers,
which trivially dominate any hierarchy partition. This inflates the Hoadley
correlation in particular because Hoadley's iCluster granularity has many
single-cancer bins.

**File:** `falsification_block_presence.csv`

---

## Test C. Granularity-nesting null on Spearman ρ

Yesterday's claim: tissue_origin shows ρ=+0.72 between 1a (5 bins) and 1b
(11 bins) → "evidence of real structure." But 1a is a NESTED aggregation of
1b (every 1b bin is a subset of one 1a bin), so the rankings are
mathematically correlated by construction.

Permutation: shuffle PIPs within each covariate, recompute the per-(hier, gran)
η², then ρ between granularities. Repeat B=1000.

**Result:**

| hierarchy | granularities | ρ_obs | perm_median | perm_p95 | p_emp |
|---|---|---|---|---|---|
| **tissue_origin** | (1a, 1b) | **+0.872** | +0.769 | +0.868 | **0.042** |
| epithelial_class | (2a, 2b) | +0.322 | +0.466 | +0.718 | 0.799 |
| epithelial_class | (2a, 2c) | +0.008 | +0.274 | +0.582 | 0.918 |
| epithelial_class | (2b, 2c) | +0.475 | +0.634 | +0.838 | 0.861 |
| germ_layer | (3a, 3b) | +0.256 | +0.564 | +0.770 | 0.968 |

**Interpretation:**
- The nesting artifact alone produces null-median ρ values of 0.27–0.77.
  These are NOT zero — granularity nesting structurally inflates ρ.
- **Tissue_origin's observed ρ=0.87 is barely above its null median 0.77** (excess = 0.10, p=0.042). Only marginally beats chance.
- All other granularity pairs have observed ρ **below the null median** —
  the apparent rankings are LESS correlated than nesting alone would predict.
  Granularity is load-bearing for those hierarchies (no robust signal).

Yesterday's "tissue_origin granularity-robust" claim is essentially refuted.
The 0.72 was mostly nesting, with a small marginal excess at p=0.042 (which
itself doesn't survive multiple-testing).

**File:** `falsification_granularity.csv`

---

## Test D. High-resolution permutation null on the 16 survivors (B=10000)

Tighten p-values for the 16 nominal survivors of Test A. Bonferroni threshold
at α=0.05 with 319 tests is **0.000157**.

**Result (sorted by p_emp):**

| cov | hier | gran | n | η_obs | p_emp (B=10000) | Bonferroni? |
|---|---|---|---|---|---|---|
| 3.2 | sex_composition | 4_4bin | 19 | 0.631 | 0.0026 | no |
| 8.2 | epithelial_class | 2a_3bin | 26 | 0.814 | 0.0057 | no |
| 20.1 | germ_layer | 3b_mid | 16 | 0.813 | 0.0072 | no |
| 10.1 | epithelial_class | 2c_10bin | 25 | 0.900 | 0.0185 | no |
| 3.1 | hoadley | 5a_4plus_unassigned | 16 | 0.579 | 0.0201 | no |
| 2.1 | hoadley | 5c_iclusters | 16 | 1.000 | 0.0225 | no |
| 36.1 | epithelial_class | 2a_3bin | 6 | 0.927 | 0.0253 | no |
| 31.1 | tissue_origin | 1b_11bin | 8 | 1.000 | 0.0378 | no |
| 1.1 | epithelial_class | 2c_10bin | 27 | 0.643 | 0.0385 | no |
| 26.1 | epithelial_class | 2c_10bin | 6 | 0.957 | 0.0419 | no |
| 26.1 | epithelial_class | 2b_6bin | 6 | 0.957 | 0.0431 | no |
| 22.1 | tissue_origin | 1b_11bin | 12 | 0.845 | 0.0461 | no |
| 8.2 | germ_layer | 3b_mid | 21 | 0.884 | 0.0487 | no |
| 8.3 | sex_composition | 4_4bin | 28 | 0.482 | 0.0502 | no |
| 26.1 | epithelial_class | 2a_3bin | 6 | 0.873 | 0.165 | no |
| 26.1 | hoadley | 5a_4plus_unassigned | 3 | 0.974 | 0.346 | no |

**Interpretation:**
- **None of the 16 survive Bonferroni correction.** Lowest p=0.0026, threshold 0.000157.
- The single most-significant hit at B=1000 (cov 26.1/hoadley, p=0.000) is now p=0.346 — was 0/1000 luck on 3 cancers.
- Even Benjamini-Hochberg FDR at q=0.05 fails: lowest p=0.0026 needs to be below 1/319 × 0.05 = 0.000157 for k=1.
- Effect sizes (η²) remain large for cov 8.2, 10.1, 1.1, 20.1 — these COULD be real but cannot be certified at n=29 cancers across 319 simultaneous tests.

**File:** `falsification_high_res_perm.csv`

---

## Test E. Cohort-size confound

TCGA cohort sizes range 30 (CHOL) to 837 (BRCA), a 28× spread. Bayesian
posterior precision varies with n; if small-n cancers cluster within
particular hierarchy bins, η² may be detecting cohort size, not biology.

### E.1 Per-cancer mean(PIP) vs log(n_samples)
- Pearson r = −0.151 (p=0.435), Spearman ρ = −0.057 (p=0.769) — **no global confound** at the cancer level.

### E.2 Cohort-size segregation across tissue_origin/1b_11bin
- ANOVA log(n) by tissue bin: F=0.93, p=0.53, **η²=0.340**
- Bin means range from log_n=4.43 (Hepatobiliary+pancreas: CHOL/PAAD/LIHC, all small) to 6.73 (Breast: BRCA alone)
- **Tissue bins DO segregate by cohort size at the screening's threshold (η²>0.30).** Not significant by F-test (small n_groups), but the effect size is at the level we'd flag.

### E.3 Partial-correlation η² for top survivors
After regressing PIP on log(n_samples) per-cancer, recomputed η² on residuals:

| cov | hier/gran | n | r(PIP, log_n) | η_orig | η_partial | Δ |
|---|---|---|---|---|---|---|
| 8.2 | epi 2a | 26 | **−0.61** (p=0.0009) | 0.814 | 0.664 | **−0.150** |
| 8.2 | germ 3b | 21 | −0.60 (p=0.004) | 0.884 | 0.829 | −0.055 |
| 10.1 | epi 2c | 25 | −0.21 (p=0.32) | 0.900 | 0.943 | +0.044 |
| 1.1 | epi 2c | 27 | −0.03 (p=0.89) | 0.643 | 0.637 | −0.007 |
| 20.1 | germ 3b | 16 | −0.26 (p=0.33) | 0.813 | 0.835 | +0.022 |
| 22.1 | tissue 1b | 12 | −0.03 (p=0.93) | 0.845 | 0.852 | +0.007 |
| 3.2 | sex 4 | 19 | −0.06 (p=0.82) | 0.631 | 0.652 | +0.021 |
| 3.1 | hoadley 5a | 16 | **−0.69** (p=0.003) | 0.579 | 0.637 | +0.058 |

**Interpretation:**
- Cohort confound is real for **cov 8.2 / epithelial 2a** (Δ=−0.15) and **cov 3.1 / hoadley 5a** (PIP-log_n r=−0.69 strong, but partialing INCREASES η² → cohort isn't driving the apparent Hoadley bin signal).
- For most other survivors, partialing log_n shifts η² by < 0.06 — cohort isn't the driver.
- Even cov 8.2 / epi 2a survives partialing at η²=0.66 (still very large).

**File:** `falsification_cohort.csv`

---

## Test F. Survivor independence

The 16 survivors aren't 16 independent confirmations.

### F.1 Cancer-set Jaccard overlap

| metric | value |
|---|---|
| Pairs with Jaccard > 0.5 | **32 / 120** |
| Pairs with Jaccard = 1.00 (identical cancer sets) | 3 |
| Highest pair-overlap clusters | cov 26.1 across epi 2a/2b/2c (all Jaccard=1.0); cov 8.2/epi 2a vs cov 10.1/epi 2c (Jaccard=0.96) vs cov 1.1/epi 2c (Jaccard=0.93); cov 8.3/sex 4 vs cov 8.2/epi 2a (Jaccard=0.93) |

### F.2 BIDIFAC+ block diversity
- 12 distinct covariate IDs in 16 survivors → **10 distinct BIDIFAC+ blocks**
- Block 26 alone provides 4 of 16 survivors (n=6 cancers in each, all Jaccard=1.0)
- Block 8 provides 3 (cov 8.2 × 2 hierarchies + cov 8.3)
- Block 3 provides 2

### F.3 PCA on survivor-covariate PIP vectors (29 cancers × 12 distinct covs)
- Singular values: 1.359, 1.069, 0.839, 0.526, 0.472, 0.377, 0.280, 0.213, 0.174, ...
- Variance explained: 41%, 25%, 16%, 6%, 5%, 3%, ...
- **5 components for 90% variance** (out of 12 input covariates)

**Interpretation:** The 16 "survivors" are testing roughly **5 independent
dimensions** of the per-cancer PIP space. The "16 hits" framing
overcounts by a factor of ~3.

**File:** `falsification_independence.csv`

---

## Synthesis — what survived

| Check | Result |
|---|---|
| Permutation null (B=1000) | At chance rate (16/319 vs expected 16) |
| Bonferroni correction (B=10000) | 0 of 16 survive |
| Block-presence confound | Hoadley confounded (ρ=+0.60), epithelial partial (ρ=+0.35), tissue clean (ρ=−0.075) |
| Granularity-nesting null | tissue_origin marginal (p=0.042), all others fail |
| Cohort-size partialing | Most survivors lose < 0.06 η², cov 8.2/epi 2a loses 0.15 |
| Independence | 16 survivors → ~5 effective independent dimensions |

**No candidate passes all five filters at conventional rigor.**

The largest, cleanest effect-size candidates that survive multiple checks
(but fail multiple-testing correction) are:

| cov | hier/gran | n | η_partial | comment |
|---|---|---|---|---|
| 10.1 | epi/2c | 25 | 0.94 | from near-universal block 10; clean of cohort, block-presence, granularity-nesting |
| 1.1 | epi/2c | 27 | 0.64 | from universal block 1; clean of all confounds |
| 8.2 | epi/2a + germ/3b | 26/21 | 0.66/0.83 | multi-hierarchy confirmation; cohort partialing reduces but doesn't kill |
| 20.1 | germ/3b | 16 | 0.84 | mid-coverage block 20 |

These four candidates have effect sizes in the η²=0.6–0.9 range across
n=16–27 cancers — far above the screening's η²>0.30 threshold — but at
permutation p-values 0.006–0.034 (none below the 0.000157 Bonferroni line).

**The honest reading:** with n=29 cancer types and 319 simultaneous
hierarchy-bin tests, the screening is fundamentally underpowered to certify
structural decomposition with publication-grade rigor. Effect-size signals
exist on a small handful of covariate-by-hierarchy cells, but they cannot
be distinguished from chance under the proper controls.

---

## What this means for the methodology paper

1. **Yesterday's headline ("30 candidates, methodology has traction") is retracted.** The screening's apparent signal is mostly bin-imbalance + block-presence + granularity-nesting artifacts.
2. **Hoadley dominance is real but circular** — the apparent +0.38 raw η² gap shrinks to +0.08 under proper filters, and even that is largely accounted for by block-presence inheritance.
3. **The methodology framing needs to change**:
   - "Closed-form structural decomposition" of the spike-and-slab requires a far more powerful test than ANOVA on per-cancer PIPs.
   - The right validation is the **reduced-parameter Gibbs replacement** the user originally proposed — fit the spike-and-slab with structural priors at tissue/epithelial granularity instead of free per-cancer PIPs, and compare held-out predictive log-likelihood. That tests whether the structural prior captures enough of the cross-cancer borrowing.
   - The validation requires either (a) a power calculation indicating n=29 is sufficient for the structural-prior contrast, or (b) extending to a larger pan-cancer cohort (HTAN, AACR GENIE).
4. **What CAN be claimed** (at the abstract level): the per-cancer PIPs do
   show large effect sizes against tissue-class hierarchies for a handful
   of BIDIFAC+ blocks, but these are not statistically certified. They are
   **leads worth following**, not findings.

---

## Files in this run

- `falsification_perm_null.csv` — 636 rows, B=1000 permutation null per test
- `falsification_block_presence.csv` — 510 rows, block-presence vs PIP-screening η² coupling
- `falsification_granularity.csv` — 5 hierarchy-pair granularity-nesting nulls
- `falsification_cohort.csv` — partial-correlation η² for top 8 survivors
- `falsification_high_res_perm.csv` — 16 survivors at B=10000
- `falsification_independence.csv` — Jaccard / block diversity / PCA summary
- `falsify_0[1-6]_*.py` — pipeline scripts
- `falsify_0[1-6]_*.log` — script logs

**Replication baseline (Lock paper A) is unchanged and still locked. The
screening is the honest answer; the methodology contribution requires a
different validation approach.**
