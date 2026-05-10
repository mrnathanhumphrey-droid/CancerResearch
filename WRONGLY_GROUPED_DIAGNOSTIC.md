# Wrongly-Grouped Diagnostic
## Paper 2 prerequisite (b) — 2026-05-10

This diagnostic was filed before any Paper 2 extension Gibbs run. It examines
whether the structural-prior groupings in Paper 1 force together cancers
whose baseline-fitted (per-cancer-free) betas point in different directions
— a "wrong grouping" pattern that would predict poorer group-shared posterior
performance.

**Important constraint:** per `PRE_REGISTRATION_PAPER2.md`, the diagnostic
informs the *discussion* of extension findings but does **not** modify the
operator-composition specification list. The 18 specs are locked
independently of these results.

---

## 1. Method

For each target covariate `k ∈ {8.2, 10.1, 1.1, 20.1}` and each candidate
hierarchy/granularity `(h, g) ∈ {epi/2a_3bin, epi/2c_10bin, germ/3b_mid,
tissue/1b_11bin}`:

1. From the Paper 1 baseline-on-train chains (4 × 100k iters, post-burn 50k,
   per-cancer free fit), extract the per-cancer posterior mean of `β_c,k`.
2. Within each hierarchy group `g`, compute:
   - `n_cancers` (in group with covariate `k` present)
   - `group_mean` of per-cancer posterior means
   - `within_group_sd` of per-cancer posterior means (between-cancer-within-group dispersion)
   - `within_group_range` (max − min)
   - `sign_agreement` (fraction of cancers in group with same sign as group mean)
3. Aggregate per `(cov, hierarchy)`:
   - `eta_squared_baseline_betas` = SS_between / SS_total over per-cancer means
   - `pooled_within_sd`
   - `mean_sign_agreement` weighted by group size

Low sign agreement and high within-group SD relative to between-group SD →
wrongly grouped.

---

## 2. Per-(covariate, hierarchy) aggregate

Sorted ascending by mean sign agreement:

| cov_id | hierarchy | n_cancers | n_groups | n_singleton | η² (baseline) | pooled within SD | mean sign agreement |
|---|---|---|---|---|---|---|---|
| **10.1** | **epi/2a_3bin** | 25 | 3 | 1 | 0.122 | 0.0206 | **0.480** ⚠ |
| **10.1** | **germ/3b_mid** | 21 | 7 | 3 | 0.135 | 0.0251 | **0.571** ⚠ |
| 8.2 | epi/2a_3bin | 26 | 3 | 1 | 0.134 | 0.0139 | 0.769 |
| 10.1 | tissue/1b_11bin | 27 | 10 | 2 | 0.165 | 0.0230 | 0.778 |
| 20.1 | epi/2a_3bin | 18 | 3 | 2 | 0.883 | 0.0422 | 0.778 |
| 8.2 | tissue/1b_11bin | 28 | 10 | 2 | 0.377 | 0.0133 | 0.786 |
| 10.1 | epi/2c_10bin | 25 | 9 | 4 | 0.759 | 0.0127 | 0.800 |
| 20.1 | tissue/1b_11bin | 20 | 8 | 1 | 0.486 | 0.0991 | 0.800 |
| 8.2 | epi/2c_10bin | 26 | 9 | 4 | 0.149 | 0.0160 | 0.808 |
| 20.1 | germ/3b_mid | 16 | 5 | 2 | 0.947 | 0.0325 | 0.812 |
| 1.1 | tissue/1b_11bin | 29 | 11 | 3 | 0.516 | 0.0513 | 0.828 |
| 20.1 | epi/2c_10bin | 18 | 8 | 3 | 0.497 | 0.1072 | 0.833 |
| 1.1 | epi/2a_3bin | 27 | 3 | 1 | 0.169 | 0.0581 | 0.852 |
| 1.1 | epi/2c_10bin | 27 | 10 | 5 | 0.555 | 0.0505 | 0.852 |
| 8.2 | germ/3b_mid | 21 | 7 | 3 | 0.217 | 0.0160 | 0.857 |
| 1.1 | germ/3b_mid | 22 | 8 | 4 | 0.758 | 0.0363 | 0.909 |

### Cells flagged (sign agreement < 0.70 with ≥ 2 non-singleton groups)

Two cells flagged:

| cov | hier | n_cancers | mean sign agreement |
|---|---|---|---|
| **10.1** | **epi/2a_3bin** | 25 | **0.480** |
| **10.1** | **germ/3b_mid** | 21 | **0.571** |

---

## 3. Worst within-group disagreements

### cov 10.1 @ epi/2a_3bin

| group | n | group_mean | within_group_sd | sign_agreement |
|---|---|---|---|---|
| **Epithelial** | **21** | **−0.0019** | 0.0207 | **0.38** ⚠ |
| Non-epithelial | 3 | −0.0226 | 0.0193 | 1.00 |

Epithelial group cancer-level posterior means (cov 10.1):
- ACC = −0.023, BLCA = +0.001, BRCA = +0.002, CESC = +0.006, CHOL = **+0.036**
- CORE = −0.000, ESCA = −0.031, HNSC = −0.000, KICH = −0.006, KIRC = −0.001
- KIRP = +0.001, **LIHC = −0.073**, LUAD = +0.001, LUSC = +0.001, OV = +0.016
- PAAD = +0.010, PRAD = +0.012, STAD = +0.004, THCA = +0.005, THYM = −0.002, UCEC = +0.000

Effects are tiny (group mean ≈ 0). Cancer-level betas span [−0.073, +0.036].
Eight cancers have negative sign, twelve have positive sign, one is zero —
the "Epithelial" group is essentially a sign-mixed bin with no consistent
direction at the level cov 10.1's per-cancer fit can resolve.

### cov 10.1 @ germ/3b_mid

| group | n | group_mean | within_group_sd | sign_agreement |
|---|---|---|---|---|
| **Foregut endoderm** | **8** | **−0.0057** | 0.0325 | **0.25** ⚠ |
| Intermediate mesoderm | 5 | −0.0024 | 0.0140 | 0.60 |
| Lateral plate mesoderm | 3 | +0.0017 | 0.0039 | 0.67 |

Foregut endoderm group (cov 10.1):
- CHOL = **+0.036**, ESCA = −0.031, **LIHC = −0.073**
- LUAD = +0.001, LUSC = +0.001
- PAAD = +0.010, STAD = +0.004, THCA = +0.005

The example flagged in the Paper 2 brief — "LUSC and LUAD pulling opposite
directions in foregut endoderm group" — does not appear at cov 10.1: LUAD
and LUSC both have posterior mean ≈ +0.001 (same direction, near zero).
The actual within-group disagreement is **CHOL = +0.036 vs LIHC = −0.073**
— a 0.11-unit spread between two hepatobiliary cancers within the same
foregut-endoderm group. CHOL and LIHC are cholangiocarcinoma (bile duct) and
hepatocellular carcinoma (liver) — both arise from hepatobiliary tissue
but with substantively different histology and behavior, and apparently
substantively different cov-10.1 effect at the per-cancer level.

The "Foregut endoderm" group at germ_layer/3b_mid is the largest mixing bin
in the diagnostic: it pools eight cancers spanning esophagus, stomach,
liver, bile duct, pancreas, lung (×2), and thyroid. These cancers are not
biologically equivalent at the gene-expression-module level for cov 10.1,
and the group-shared posterior mean wipes out their distinct
cancer-specific contributions.

---

## 4. Cells **not** flagged (sign agreement ≥ 0.70)

Most cells pass the sign-agreement threshold:

- **cov 8.2 at all four hierarchies** (sign agreement 0.77 to 0.86). Paper 1's
  primary spec used cov 8.2 @ epi/2a (sign agreement 0.77) and the held-out
  Δ was +6.26 nats, supporting the spec choice.
- **cov 8.2 @ germ/3b_mid (Paper 1 secondary, sign agreement 0.86)** — actually
  *more* consistent than the primary anchor. Paper 1's secondary held-out Δ
  was inconclusive (+4.33, CI [−1.86, +10.97]) but the issue is **not**
  wrong-grouping. Likely a smaller effect size given the broader cancer set
  (22 cancers vs 26 at epi/2a).
- **cov 1.1 @ all four hierarchies** (sign agreement 0.83 to 0.91). The
  highest sign-agreement cell is cov 1.1 @ germ/3b_mid (0.91).
- **cov 20.1 @ all four** (sign agreement 0.78 to 0.81). Effects are
  larger (group_mean range −0.50 to +0.20) so signs are more stable.

---

## 5. Implications for the Paper 2 extension

The pre-registration locks the 18 specifications regardless of diagnostic
findings. The diagnostic informs the **interpretation** of post-compute
results:

1. **Cov 10.1 results may be inflated by regularization rather than recovered
   structure.** With baseline per-cancer betas tiny and sign-mixed within
   epi/2a_3bin and germ/3b_mid groups, group-shared priors will shrink
   per-cancer betas toward their group means (≈ 0). LPPD may improve
   relative to baseline simply because the group prior is a near-zero
   regularizer for tiny effects, not because the group recovers a real
   compositional pattern. This is a known failure mode for the extension
   claim and should be flagged in the Paper 2 discussion if cov-10.1
   compositions show the largest improvements.

2. **Cov 8.2 results are the cleanest test.** Sign agreement is moderately
   high at all four hierarchies (0.77–0.86); within-group SD is small
   (0.013–0.025). Operator compositions on cov 8.2 will be the most
   informative for the extension claim because the per-cancer signal is
   resolvable enough to detect compositional structure if it exists.

3. **Cov 1.1 results may be over-regularized in the opposite direction.**
   At germ/3b_mid, sign agreement is 0.91 — extremely consistent — but
   pooled within-group SD is small (0.036). The group prior may recover
   a near-baseline fit (improvement marginal). Cov 1.1 compositions may
   "match" Paper 1 rather than improve it.

4. **The Paper 2 brief's LUSC/LUAD-at-foregut-endoderm hypothesis is not
   supported by the diagnostic** for cov 10.1 (both have posterior mean
   ≈ +0.001). The actual largest within-group disagreement at germ/3b_mid
   foregut endoderm is **CHOL vs LIHC** (a 0.11-unit posterior-mean spread).
   This refines the wrongly-grouped hypothesis: it's not lung-vs-lung
   disagreement, it's hepatobiliary-vs-other-foregut disagreement.

5. **Tissue origin (1b_11bin) tends to have higher sign agreement than germ
   layer or coarse epithelial classification** for the targets considered
   (e.g., cov 8.2: 0.79 at tissue vs 0.77 at epi/2a; cov 10.1: 0.78 at
   tissue vs 0.48 at epi/2a). Pair B (epi/2c × tissue/1b_11bin) is
   plausibly cleaner than Pair A (epi/2a × germ/3b) for the operator-
   composition tests on cov 8.2 and cov 10.1.

---

## 6. Decision unmodified

The Paper 2 specifications are unchanged. All 18 fire as pre-registered.
This document is committed before any extension Gibbs run as required by
the pre-registration's prerequisite (b).

---

## 7. Files

| file | contents |
|---|---|
| `runs/run_extension_2026-05-10/per_cancer_baseline_betas.csv` | per-cancer × per-cov posterior summary from Paper 1 baseline (104 rows) |
| `runs/run_extension_2026-05-10/wrongly_grouped_diagnostic.csv` | per-(cov, hier, group) within-group statistics (114 rows) |
| `runs/run_extension_2026-05-10/diag_within_group_summary.csv` | per-(cov, hier) aggregate (16 rows) |
| `runs/run_extension_2026-05-10/00_wrongly_grouped_diagnostic.R` | reproducer script |

(Paths are in the working tree — public-repo copies under `results/diagnostic/`.)
