# Lock-Disposition Diagnostic Report — Cancer Paper 1

**Generated:** 2026-05-11 00:25 -04:00
**Pre-registration:** `LOCK_DIAGNOSTIC_PRE_REGISTRATION.md` (committed 2026-05-11 00:20:17 -04:00 before any diagnostic output viewed)
**Working data:** validation per-patient LPPDs at `results/validation/loglik_per_patient_{baseline_train,primary,primary_3cov}.csv`; finer-resolution PIP test outputs at `results/extension/{option_a_anova,option_c_anova,finer_resolution_residuals}.csv`.

---

## Headline

**Final disposition: 3-cov MATCHES baseline with parameter reduction, coarse-resolution structural signal (D1 = WEAKENED, D2 = Reading A — both checks agree).**

Cancer Paper 1 ships as: *structural priors derived from epithelial classification and embryonic germ layer match the per-cancer fitted-hierarchical baseline on held-out predictive log-likelihood, with parameter reduction; cov 8.2's structural signal is directionally consistent at subtype-within-tissue granularity but does not reach conventional significance at the finer scale.*

The 4-cov "+6.26 nats outperforms" verdict is **retracted from the headline claim** and documented as a sensitivity: the +6.26-nat margin is partially carried by cov 10.1, whose individual contribution to the LPPD is statistically indistinguishable from zero by both per-cancer (CHECK 1) and bootstrap (CHECK 2) diagnostics.

---

## Diagnostic 1 — cov 8.2 multi-resolution claim

### Input

Finer-resolution test executed at `scripts/extension/02_finer_resolution_cov8_2.py`. **The test as run is an ANOVA on cov 8.2's per-cancer PIP residuals (after epi/2a group-mean subtraction), not a held-out Gibbs/LPPD test.** The locked decision rules were written assuming LPPD Δ-with-CI; the actual diagnostic is η²-with-permutation-p. This is documented here as a pre-registered edge case and the closest locked category is applied.

### Convergence / interpretability check

For the ANOVA-form diagnostic, "interpretable" means the partition has at least one bin with ≥ 2 cancers (otherwise η² saturates mechanically and the permutation null cannot bound the observation).

| Group | n | Option A (tissue/1b) bins | Option C (hybrid) bins | Interpretable |
|---|---|---|---|---|
| Epithelial | 21 | 9 (2 singleton, 7 multi) | 19 (18 singleton, 1 multi) | A: yes / C: degenerate |
| Hematological | 1 | 1 (singleton) | 1 (singleton) | no |
| Non-epithelial | 4 | 4 (all singletons) | 4 (all singletons) | no |

The only group with statistical power to test the multi-resolution claim is **Epithelial** under **Option A** (tissue-within-epithelial-class, 21 cancers, 9 tissue bins with 7 containing 2+ cancers each). Option C's hybrid label (subtype where tissue is multi-cancer, tissue otherwise) is statistically degenerate in this group: 18 of 19 labels are singletons, so η² mechanically saturates near 1.

### Results

| Test | η²_obs | perm_median | perm_p95 | p_emp | F | p_F |
|---|---|---|---|---|---|---|
| Option A: residuals ~ tissue (Epithelial group, n=21) | 0.608 | 0.381 | 0.685 | **0.117** | 2.328 | 0.091 |
| Option C: residuals ~ hybrid (Epithelial group, n=21, degenerate) | 0.988 | 0.953 | 0.998 | **0.232** | 9.079 | 0.104 |

Both tests have η² above the permutation null median (directional consistency with the coarse-resolution screening signal). Neither test reaches p < 0.05 under any of the three statistics (permutation, F-test, or permutation 95th percentile).

### Classification per locked rules

- **VALIDATED** (rule: same direction AND CI excludes zero) → **NO**. p_emp = 0.117 (Option A) and 0.232 (Option C); neither excludes zero.
- **WEAKENED** (rule: same direction AND CI overlaps zero) → **YES**. η²_obs > perm_median for both tests; both p_emp > 0.05.
- **FALSIFIED** (rule: effect disappears or reverses) → **NO**. η²_obs is well above perm_median, not at it.

**D1 disposition: WEAKENED.**

The multi-resolution claim narrows from "biological structure persists at subtype-within-tissue granularity" to **"coarse-resolution structural signal in cov 8.2 with directional consistency at finer resolution that does not reach conventional statistical significance at the data available."**

Pre-registered caveat documented: Option C's ANOVA is statistically degenerate (18/19 singleton bins) and cannot itself validate or falsify the multi-resolution claim; Option A on tissue-within-epi-class is the clean sensitivity check and gives the same WEAKENED disposition.

Source CSVs: `results/extension/option_a_anova.csv`, `results/extension/option_c_anova.csv`, `results/extension/finer_resolution_residuals.csv`.

---

## Diagnostic 2 — cov 10.1 disposition

### CHECK 1 — Cross-group effect analysis

Per-cancer LPPD contribution of cov 10.1, computed as (4-cov per-patient LPPD aggregated by cancer) − (3-cov per-patient LPPD aggregated by cancer). Both LPPDs derived from existing Gibbs chains for the primary (4-cov) and primary_3cov specs.

**Distribution summary (n = 29 cancers):**
- Positive contribution (cov 10.1 helps): 15 cancers. Largest: LIHC +0.902, LUSC +0.626, SARC +0.543.
- Negative contribution (cov 10.1 hurts): 14 cancers. Largest: STAD −0.439, KIRP −0.274, CHOL −0.257.
- |Contribution| > 1 nat: **0 cancers** (largest is LIHC at +0.902).
- |Contribution| > 0.5 nat: 3 cancers, all positive (LIHC, LUSC, SARC).
- Binomial test against H₀ = random ± (15 of 29 positive): two-sided p = 1.000.

**Decision rule:** 3+ cancers with |contrib| > 1 nat in same direction AND binomial p < 0.05 → Reading B.

Result: 0 cancers cross |1 nat|. Binomial test on full distribution is at H₀ exactly. At the relaxed |0.5 nat| threshold (3 positive of 3 with |contrib| > 0.5), binomial p = 0.250 — also fails p < 0.05.

**CHECK 1 verdict: Reading A.**

Source: `results/diagnostic/cov_10_1_cross_group.csv`.

### CHECK 2 — Bootstrap distribution of cov 10.1 marginal

Patient-level bootstrap, B = 1000, seed 20260510 (matches existing bootstrap conventions). For each resample: marginal = (Σ 4-cov per-patient LPPD over resample) − (Σ 3-cov per-patient LPPD over resample).

| Metric | Value |
|---|---|
| Observed marginal | **+1.420 nats** |
| Bootstrap 95% CI | **[−1.655, +4.400]** |
| Bootstrap median | +1.459 |
| Excludes zero (positive direction) | **FALSE** |
| Excludes zero (negative direction) | FALSE |

**Decision rule:** Bootstrap CI excludes zero positive → Reading B. CI brackets zero → Reading A.

Result: CI brackets zero. Observed marginal is positive but cannot be statistically distinguished from zero.

**CHECK 2 verdict: Reading A.**

Sources: `results/diagnostic/cov_10_1_bootstrap.csv` (full B=1000 distribution), `results/diagnostic/cov_10_1_bootstrap_summary.csv` (CI + verdict), `scripts/diagnostic/check2_bootstrap.R` (reproducibility).

### Combined Diagnostic 2 classification

CHECK 1 = Reading A. CHECK 2 = Reading A. **No conflict; conservative-default rule not invoked.**

**D2 disposition: Reading A.** Cov 10.1's contribution to the 4-cov spec is neither concentrated in specific cancers (no per-cancer crosses |1 nat|; 15/14 split is chance-tier) nor statistically separable from zero at the spec level (bootstrap CI brackets zero). Final spec drops cov 10.1.

---

## Combined disposition per locked decision tree

D1 = WEAKENED × D2 = Reading A → **3-cov MATCHES with coarse-resolution structural signal**.

### What the paper claims now

| Element | Status |
|---|---|
| Headline | "Structural priors derived from epithelial classification and embryonic germ layer **match** the per-cancer fitted-hierarchical baseline on held-out predictive log-likelihood, with parameter reduction." |
| Primary spec | **3-cov** (1.1, 8.2 @ epi/2a; 20.1 @ germ/3b) |
| Primary Δ | **+4.84 nats vs baseline-on-train**, bootstrap 95% CI **[−0.82, +10.05]** (CI brackets zero → "matches" disposition; point estimate consistent with structural prior performing at least as well as baseline) |
| Free-parameter reduction | 3 × 29 = 87 per-cancer βs → 22 group-shared βs after dropping contested cells: **75% reduction** on the three structural-prior covariates |
| Cov 8.2 multi-resolution claim | "Coarse-resolution structural signal with directional consistency at finer resolution that does not reach conventional statistical significance" |
| 4-cov spec | Reported as sensitivity; explicitly noted that the +6.26-nat margin is partially carried by cov 10.1, whose individual contribution is statistically indistinguishable from zero |

### What changes from the existing README

The README headline currently says "Validation outperforms baseline by **+6.26 nats**" (4-cov). The Lock disposition changes this to:

- Headline becomes "matches with parameter reduction" using the 3-cov spec.
- The 4-cov result is preserved as a sensitivity analysis with the cov-10.1-noise caveat documented in the body of the paper.
- Section 5 (validation methodology) keeps both primary and sensitivity specs; section 1 (headline result) and section 6 (reconciliation with screening) get rewritten for the matches framing.
- Section 2 ("what was tested and what wasn't") gets a new bullet documenting that cov 10.1's marginal contribution failed both per-cancer and bootstrap diagnostics, motivating the drop.

The README rewrite itself is a separate edit; this document only records the disposition.

---

## Deviations from pre-registration

One pre-registered edge case:

- **Diagnostic 1 was an ANOVA-on-PIP-residuals test, not a Gibbs/LPPD test.** The locked decision rules described the substantive question (does cov 8.2 effect persist at subtype-within-tissue granularity?) and provided VALIDATED / WEAKENED / FALSIFIED categories framed in LPPD Δ-with-CI language. The actual test available is η²-with-permutation-p. The pre-registration's EDGE-CASE clause ("If diagnostic outputs are inconsistent with pre-registered decision categories ... document the issue and report which combination is closest. Do not silently classify edge cases.") was invoked. The closest category was applied (WEAKENED) and the test-form mismatch is documented explicitly here.

No other deviations.

---

## Compute usage

| Step | Wall time |
|---|---|
| Diagnostic 1 (Option C ANOVA result already complete) | n/a — pre-existing |
| CHECK 1 (per-cancer LPPD diff from existing per-patient files) | seconds |
| CHECK 2 (patient-level bootstrap B=1000 on existing per-patient LPPDs) | ~2 seconds |
| Total new compute | < 1 minute |
