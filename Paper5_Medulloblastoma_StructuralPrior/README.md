# Paper 5 — Within-Cancer Structural-Prior Methodology Transfer to Medulloblastoma

**Status (2026-05-18): PRELIMINARY RESULT, NEEDS FURTHER VALIDATION.**

A full leakage-clean rerun completed (24-chain Gibbs, 6 specs, all R-hat <
1.01). Primary spec exceeds baseline by **+6.7 nats on held-out n=118**
(+0.058 per patient, Lock-comparable). But a projection-scale defect in
the BIDIFAC+ test-fold scores (53 of 55 components show test_sd / train_sd
> 1.5, median 5.1×, max 714×) makes the held-out LPPD partially
artifact-laden. Under four reasonable post-hoc rescalings, primary Δ
ranges from **+0.7 to +9.2 nats**; sensitivity Δ ranges from **−24.2 to
+13.0 nats**. The direction of the primary-vs-sensitivity comparison
flips between rescalings. **No defensible verdict yet.**

The full clean rerun + diagnostic trail is committed for reproducibility;
the paper-disposition decision is deferred until a Tikhonov-regularized
re-projection + Gibbs refit lands. This is cancer research — better to
publish nothing than to publish a number we cannot stand behind.

See `DEVIATIONS.md` for the full bug-discovery trail and
`results/lppd_rescaling_diagnostic.csv` for the rescaling sensitivity
table.

---

## What this paper tests

The Lock 2022 hierarchical-spike-and-slab structural-prior methodology,
validated on Paper 1 against a pan-cancer (29-cancer) cohort with
external epithelial-classification and germ-layer hierarchies, is here
applied to a **single-cancer** corpus with the structural-prior hierarchy
being the **within-cancer molecular subtype taxonomy**.

Substrate: medulloblastoma. Cohort: Cavalli et al. 2017, *Cancer Cell*,
n = 763 primary tumors with matched DNA methylation (Illumina 450k) and
gene expression (Affymetrix Hu Gene 1.1 ST) — both fully public via GEO
accessions GSE85212 and GSE85217. Hierarchy: Cavalli's 12-subtype
molecular classification (primary granularity), with the 4-subgroup WHO
consensus as a coarser sensitivity granularity.

This is a **within-cancer test** complementary to Paper 4 Phase 1's
**between-hierarchy test** on the same 29-cancer pan-cancer substrate.
Paper 4 found the methodology is hierarchy-specific (0/3 outperform
under matched specs). Paper 5 asks whether the methodology transfers to
a substrate where the hierarchy is naturally within-cancer rather than
imposed externally.

---

## Result snapshot (NOT FINAL)

Leakage-clean Gibbs run, 24 chains × 6 specs, all R-hat < 1.01:

| spec              | LPPD     | Δ-vs-baseline | per-patient |
|-------------------|----------|---------------|-------------|
| baseline_train    | -145.08  | 0.00          | —           |
| primary (12 L2, L1-pooled) | -138.40  | **+6.68**     | +0.057      |
| sensitivity (4 L1, identity) | -132.08  | **+13.00**    | +0.110      |
| primary_drop_t1   | -143.23  | +1.85         |             |
| primary_drop_t2   | -138.28  | +6.80         |             |
| primary_drop_t3   | -143.72  | +1.36         |             |

The numbers above are **as-computed**. They are not the final verdict.
See the sensitivity-to-rescaling table below.

### Why the numbers above are not the verdict

The BIDIFAC+ projection of held-out (n=118) test scores onto
training-derived module loadings inflates the test-fold standard
deviation:

| metric | value |
|---|---|
| Components with test_sd / train_sd in [0.5, 1.5] (stable) | 2 of 55 |
| Median test_sd / train_sd ratio | 5.13 |
| Max ratio | 713.8 (column M02_joint_r19, d_min/d_max = 0.0002) |

Root cause: per-rank projection V_test = (X_test^T u_j) / d_j is
sensitive to small singular values inside each BIDIFAC+ module. Modules
with wide d-spectra produce stable scores at high ranks and noise-
amplified scores at low ranks. The locked target covariates (cov 43, 25,
23) sit at intermediate ranks (ratios 8.93 / 5.74 / 3.96), neither
stable nor catastrophic.

A post-hoc diagnostic re-evaluated LPPD across four reasonable rescalings
of the test fold:

| spec              | orig | own | clip | drop_worst |
|-------------------|------|-----|------|------------|
| primary           | +6.68 | +6.79 | +9.24 | +0.69 |
| sensitivity       | +13.00 | +3.84 | +0.71 | −24.19 |
| primary_drop_t1   | +1.85 | +1.54 | +1.23 | +2.31 |
| primary_drop_t2   | +6.80 | +5.43 | +9.83 | +0.73 |
| primary_drop_t3   | +1.36 | +3.81 | +6.04 | −1.87 |

- `orig`: current pipeline output (train-fold standardization applied to test)
- `own`: test fold z-scored by its own fold stats (folds symmetrized)
- `clip`: test scores clipped to |z| < 5 per column
- `drop_worst`: columns with test_sd / train_sd > 5 zeroed in test (β unchanged)

Primary is not robust (range +0.7 to +9.2). Sensitivity is even less
robust (range −24 to +13). Primary vs sensitivity flips direction
between rescalings.

### What this means

The leakage-clean rerun got past three pre-reg-fatal bugs (DEVIATIONS.md
Entry 001): the LPPD units mismatch, the spec collapse, and the
BIDIFAC+ + screening leakage. But a fourth issue surfaced in
post-fit diagnostics — the projection of test-fold scores is not
numerically stable for low-d module components. The Gibbs spike-and-slab
shrinks β toward zero for the noisiest components (PIPs 0.07–0.21 for
the worst 10), but the locked target covariates are force-included by
the structural prior at their as-projected scale, and that scale is 4–9×
the training scale.

The cleanest fix: re-project test scores with Tikhonov regularization
(λI added to the per-module Gram matrix), refit Gibbs end-to-end, and
recompute LPPD on coordinate-aligned train and test. This is roughly
another 4–8 hours of compute. The current `_clean` results are the
intermediate state.

---

## Why medulloblastoma + Cavalli 2017

1. **Three nested granularity levels are native to the data.** Lock 2022
   needed to scan 5 external hierarchies × 11 granularities to find the
   epi-class / germ-layer combination. Medulloblastoma comes pre-
   annotated with a 4-subgroup / 12-subtype / 8-subtype-within-G3-G4
   (Sharma 2019) hierarchy — the granularity sweep is built in, not
   imposed.
2. **WNT singleton problem has a clean fix at L2.** The 4-subgroup level
   leaves WNT (n = 71, ~4 events) as a structural-prior singleton with
   no neighbors to pool with. The 12-subtype level splits WNT into
   WNT-α and WNT-β, which share a parent group.
3. **Data is fully public and matched 1:1.** GSE85212 (methylation) and
   GSE85217 (expression) cover the same 763 samples; no dbGaP gating.
4. **Two-block BIDIFAC+ scales down from Lock's four-block setup with
   no algorithm change.**

---

## Folder contents

```
Paper5_Medulloblastoma_StructuralPrior/
├── README.md                ← this file
├── PRE_REGISTRATION.md      ← locked pre-reg (methodology + hierarchy + decision rules)
├── DEVIATIONS.md            ← Entry 001: three bugs discovered post-fit + fixes
├── data/                    ← Cavalli expression + methylation + clinical
├── reference/               ← subtype distribution CSV, locked target-covariate IDs,
│                              split indices. paper5_target_covariates_clean.csv is
│                              the SHA-anchored clean-screening lock.
├── scripts/                 ← runners, v1 + v2 (leakage-clean), diagnostics
├── results/                 ← v1 LPPD (contaminated, retained for audit), logfix
│                              LPPD (units-fix on v1 chains), _clean LPPD (leakage
│                              + units fixed), rescaling-sensitivity diagnostic
└── logs/                    ← orchestrator + per-step run logs from overnight pipeline
```

---

## Pre-reg highlights (load-bearing decisions; original pre-reg unchanged)

- **Primary granularity: L2 12-subtype Cavalli.** L1 4-subgroup is the
  sensitivity row.
- **Target covariates: 3.** Original lock at `paper5_target_covariates.csv`
  (covs 20/15/57) was screened on contaminated full-data BIDIFAC+
  components. Re-lock at `paper5_target_covariates_clean.csv` (covs
  43/25/23) was screened on training-fold-only BIDIFAC+. Both lock files
  are preserved for reproducibility; the clean lock supersedes the
  contaminated one per DEVIATIONS.md Entry 001.
- **Per-covariate CHECK 1 + CHECK 2 marginal-contribution cascade
  committed at pre-reg time.** Under the orig rescaling, CHECK 1 passes
  2/3 (t2 ≈ 0); under the own rescaling, CHECK 1 passes 3/3.
- **Held-out split: 80/20 stratified at L2, seed 20260515.** Unchanged.
- **Gibbs parameters: 4 chains × 100k iters.** Unchanged. Clean-rerun
  seeds 20260518–20260521.
- **Decision rule: matches if |Δ| ≤ 2 nats OR CI brackets zero;
  outperforms if Δ > +2 AND CI excludes zero.** Pending refit before any
  decision rule is applied.

---

## What's not in Paper 5 (deferred)

- Sharma 2019 G3/G4 8-subtype L3 granularity → external-validation phase.
- MAGIC n = 898 → MAGIC requires dbGaP access; deferred.
- Other within-cancer subtype taxonomies → separate pre-registrations.
- PFS endpoint → OS only.

---

## Position in the corpus

| paper | substrate | hierarchy | verdict |
|---|---|---|---|
| Paper 1 | 29-cancer TCGA pan-cancer | epi-class + germ-layer (external) | **3-cov MATCHES with parameter reduction** (post-disposition); 4-cov sensitivity +6.26 nats |
| Paper 2 | same as Paper 1 | operator-composition extension | FALSIFIES; L1 retro = cancer-identity-dominant |
| Paper 3 | breast Xenium + Chromium (MOFA-FLEX) | spatial-niche | FALSIFIED across 4 operator families; substrate boundary |
| Paper 4 Phase 1 | same as Paper 1 | 3 alternative pan-cancer hierarchies | PAPER1_HIERARCHY_SPECIFIC; 0/3 outperform |
| **Paper 5** | **Cavalli n = 763 medulloblastoma** | **within-cancer 12-subtype** | **PRELIMINARY +6.7 nat primary Δ; needs Tikhonov re-projection + refit before any verdict** |
