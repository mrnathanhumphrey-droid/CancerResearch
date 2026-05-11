# Lock-Disposition Pre-Registration — Cancer Paper 1

**Committed:** 2026-05-11 00:20:17 -04:00
**Status when committed:** Diagnostic outputs NOT YET viewed. Decision rules locked here are applied after this file is written.

This file is the locked pre-registration. Any deviation from these rules must be documented explicitly in `LOCK_DIAGNOSTIC_DISPOSITION.md` as a deviation, not silently applied.

---

## Background

Cancer Paper 1 validation Gibbs completed on 4 pre-registered covariates with 4 hierarchies × granularities. Primary 4-cov spec (10.1, 1.1, 8.2 → epi/2a, 20.1 → germ/3b) produced Δ = +6.26 nats vs baseline-on-train, bootstrap 95% CI [+0.34, +12.29], "OUTPERFORMS" per the original validation pre-registration.

Two open diagnostics block final disposition:

1. **Diagnostic 1 — Option C finer-resolution test on cov 8.2.** Already launched as an ANOVA-on-PIP-residuals test (`scripts/extension/02_finer_resolution_cov8_2.py`). Confirms or falsifies the multi-resolution biological signal claim.

2. **Diagnostic 2 — Cov 10.1 disposition.** The 3-cov vs 4-cov bootstrap (`scripts/validation/02b_compute_loglik_3cov.R`) showed dropping cov 10.1 brings the result to Δ = +4.84 nats; the CI overlapped zero. Cov 10.1's individual within-group sign agreement is at chance (0.48). The 4-cov "outperforms" verdict may be partially carried by a noise-tier covariate. Background: `COV_10_1_INVESTIGATION.md`.

---

## Diagnostic 1 decision rules (cov 8.2 multi-resolution claim)

The Option C result will be pulled and classified per the locked categories below. Convergence is checked first (R-hat ≤ 1.05, ESS adequate) — for the ANOVA-form Option C, "convergence" means: ANOVA test is not degenerate (i.e., at least 2 non-singleton bins exist in the partition).

- **VALIDATED:** cov 8.2 effect persists at subtype-within-tissue granularity. Operationally: finer-resolution test point estimate is in the same direction as the primary screening effect AND the test's permutation/bootstrap CI excludes zero (or p < 0.05 against the structural null). Multi-resolution biological structure claim stands.

- **WEAKENED:** effect attenuates at finer resolution but maintains directionality. Operationally: point estimate same direction as primary AND CI overlaps zero (or 0.05 ≤ p < permutation-null upper bound). Multi-resolution claim narrows to "coarse-resolution structure with directional consistency at finer resolution."

- **FALSIFIED:** effect disappears or reverses at finer resolution. Operationally: point estimate opposite direction OR effectively at permutation null (eta² within permutation median ± mad, p approaches null mean). Multi-resolution claim retracts. Cov 8.2 stays in 4-cov spec but methodology framing in paper narrows.

- **EDGE CASE (locked):** If the Option C result is uninterpretable due to convergence failure, near-degenerate bin structure (singletons dominating to the point the ANOVA is statistically meaningless), or any other reason that prevents clean classification — document the issue and report which combination is closest. Do not silently classify edge cases.

## Diagnostic 2 decision rules (Reading A vs Reading B on cov 10.1)

Two checks run independently. Decisions combined per the conservative-default rule below.

**CHECK 1 — Cross-group effect analysis on cov 10.1.**
- Compute per-cancer LPPD contribution: (4-cov per-patient LPPD aggregated by cancer) − (3-cov per-patient LPPD aggregated by cancer).
- Output to `results/diagnostic/cov_10_1_cross_group.csv`.
- Decision: If **3 or more cancers** show |individual contribution| > 1 nat **in the same direction** AND a binomial test against null "random ±" yields p < 0.05 — **Reading B holds**. Otherwise **Reading A holds**.

**CHECK 2 — Bootstrap distribution of cov 10.1 marginal at the spec level.**
- For B = 1000 patient-level resamples of the held-out test set, compute (sum 4-cov per-patient LPPD over resampled patients) − (sum 3-cov per-patient LPPD over resampled patients).
- Output to `results/diagnostic/cov_10_1_bootstrap.csv`.
- Decision: If bootstrap 95% CI on the marginal **excludes zero in the positive direction** (i.e., 2.5th percentile > 0) — **Reading B holds**. If CI brackets zero — **Reading A holds**.

**Conflict resolution (locked):**
If CHECK 1 and CHECK 2 disagree, **default to Reading A** on conservative grounds: noise-tier individual contribution under multiple diagnostics is more parsimoniously explained as fit artifact than as real signal with subtle structure.

**Reading A (cov 10.1 dropped):** Final spec is 3-cov, "matches" disposition. Paper claim is "structural priors match the fitted-hierarchical baseline with parameter reduction." Quantify reduction: cov 1.1 + cov 8.2 + cov 20.1 free `β_jk` drop from 3 × 29 = 87 (per-cancer) to 22 (group-shared, after dropping contested cells) — 75% reduction on those three.

**Reading B (cov 10.1 kept):** Final spec stays 4-cov, "outperforms" disposition. Original Δ = +6.26 nats, CI [+0.34, +12.29], with the principled rationale that cov 10.1 was pre-registered, passed Bonferroni at screening, and contributes via cross-group pooling even if its within-group sign agreement is at chance.

## Combined disposition decision tree (LOCKED)

| Diagnostic 1 | Diagnostic 2 | Disposition |
|---|---|---|
| VALIDATED | Reading B (keep) | **4-cov OUTPERFORMS with multi-resolution structural signal** |
| VALIDATED | Reading A (drop) | **3-cov MATCHES with multi-resolution signal in cov 8.2** |
| WEAKENED | Reading B | **4-cov OUTPERFORMS with coarse-resolution signal, narrowed framing** |
| WEAKENED | Reading A | **3-cov MATCHES with coarse-resolution signal** |
| FALSIFIED | Reading B | **4-cov OUTPERFORMS at coarse resolution only; multi-resolution claim retracted** |
| FALSIFIED | Reading A | **3-cov MATCHES at coarse resolution only** |

## Constraints

- Pre-registration committed (this file) before any diagnostic output is viewed.
- Diagnostic 1 result incorporation happens first.
- Diagnostic 2 CHECK 1 and CHECK 2 both run before Reading A vs B resolution.
- No spec is re-fit, replaced, or augmented after this file is written. Existing Gibbs chains and per-patient LPPDs are the working data.

## Deliverables

- `LOCK_DIAGNOSTIC_DISPOSITION.md` — top-level report with (a) D1 classification, (b) D2 classification, (c) combined disposition, (d) all check outputs documented
- `results/diagnostic/cov_10_1_cross_group.csv` — per-cancer LPPD contribution from cov 10.1 (CHECK 1)
- `results/diagnostic/cov_10_1_bootstrap.csv` — bootstrap distribution of cov 10.1 marginal (CHECK 2)
- `results/diagnostic/cov_10_1_bootstrap_summary.csv` — CHECK 2 CI + verdict
- `scripts/diagnostic/check2_bootstrap.R` — CHECK 2 reproducibility

(The original brief also requested a Lock email draft; per agent-role boundary the email is user-drafted, not in this deliverable set.)

---

*End of pre-registration. Diagnostic outputs read starting now.*
