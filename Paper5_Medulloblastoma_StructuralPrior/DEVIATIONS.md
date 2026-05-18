# Paper 5 — Pre-registration Deviations

Living log of methodology deviations from the locked PRE_REGISTRATION.md
(committed at SHA `83ea4e6`, covariate lock at SHA `34db626`). Each entry
documents the deviation, when it was discovered, why it was necessary,
and how it was handled.

---

## Entry 001 — Two pre-reg-fatal bugs discovered post-fit, full pipeline rebuild

**Discovered:** 2026-05-17, immediately after the 24-chain Gibbs landed
and held-out LPPD was first computed.

**Symptom:** v1 reported `primary` outperforming `baseline_train` by
+349.115 nats on n=118 test patients (per-patient Δ = +2.96 nats), versus
Lock 2022 Paper 1's per-patient Δ ≈ +0.04 nats — a 60× anomaly. Two
sensitivity-equivalent specs (`primary` and `sensitivity`) reported
bit-identical LPPDs to 4 decimals.

**Root causes identified:**

### Bug 1 — Held-out LPPD units mismatch

`scripts/compute_held_out_loglik_paper5.R` line 47 evaluated:
```r
ll[uncens] <- dnorm(y_test[uncens], mean = mu[uncens], sd = sigma, log = TRUE)
```
where `y_test` is raw survival in years (from
`pmax(clin_s$os_time_years[keep], 1/365)`) and `mu`/`sigma` come from
`out$betas` and `out$sigma2`. The Lock sampler is log-normal AFT; it
applies `log(Y)` internally during fitting, so the stored β are in
log-time. The v1 LPPD therefore evaluated `dnorm(y_years, mu_log_time,
sigma_log_time)` — incoherent units.

**Magnitude:** recomputing LPPD on the v1 (contaminated) chains with
`dnorm(log(y_test), mu, sigma)`:

| spec              | v1 Δ (buggy) | logfix Δ | shift   |
|-------------------|--------------|----------|---------|
| baseline_train    | 0.00         | 0.00     | 0.00    |
| primary           | +349.12      | +10.87   | -338.24 |
| sensitivity       | +349.12      | +10.87   | -338.24 |
| primary_drop_t1   | +285.45      | +7.43    | -278.02 |
| primary_drop_t2   | +171.93      | +6.71    | -165.23 |
| primary_drop_t3   | +318.48      | +9.78    | -308.70 |

The units bug accounted for **97%** of the v1 headline. Per-patient Δ on
the units-corrected v1 chains: +0.092 nats (≈ 2× Lock), order of
magnitude consistent with literature, not anomalous.

**Fix:** `compute_held_out_loglik_paper5_v2.R` applies `log()` before
`dnorm`/`pnorm`. Stored as
`results/loglik_summary_logfix.csv` for the v1-chains-units-fixed run.

### Bug 2 — Sensitivity spec collapse

`scripts/build_paper5_split.R` lines 121–127 built `target_covs_primary`
and `target_covs_sensitivity` using the same `Covariates/Survival/Censored`
lists (12 L2 subtypes) and the same `group_per_cancer = l1_per_subtype`
vector. As a result, the two specs are bit-identical — `primary` and
`sensitivity` reported the same LPPD to all digits in both the buggy and
logfix recomputes.

**Intended sensitivity** (per PRE_REGISTRATION.md §3): cancer-type axis
= 4 L1 subgroups (WNT/SHH/Group3/Group4) with identity hierarchy (no
further pooling), versus primary's 12 L2 subtypes pooled by L1.

**Fix:** `build_paper5_split_v2.R` builds:
- `train_data_clean.rda` / `test_data_clean.rda` — L2 12-subtype lists for
  primary + ablation specs
- `train_data_L1_clean.rda` / `test_data_L1_clean.rda` — L1 4-subgroup
  lists for sensitivity
- `target_covs_primary_clean.rda` — L2 with L1-pooling structural prior
- `target_covs_sensitivity_clean.rda` — L1 with identity hierarchy

`run_gibbs_paper5_v2.R` loads the appropriate train_data per spec.

### Bug 3 — BIDIFAC+ + screening leakage (the originally-flagged issue)

`run_bidifac_plus.R` fit BIDIFAC+ on the full Cavalli cohort (n=763)
including the held-out 20% test fold. Components therefore encode
test-set covariance. `screening_paper5.R` then ran ANOVA on PIPs from a
baseline Gibbs fit on full data (n=763), so the locked target covariates
{20, 15, 57} were selected with test-set information. Both contaminate
held-out LPPD.

**Decision:** the locked covariate IDs cannot be preserved across a
clean rerun (BIDIFAC+ training-only produces a different module
decomposition; identities don't transfer). The cleanest path is full
pipeline rebuild with re-screening on training data only, accepting that
this deviates from the §10 constraint 3 ("no covariate substitution
after lock"). The deviation is principled — the original lock was made
on a methodologically invalid screening — not cherry-picking.

**Fix:**
1. `run_bidifac_plus_leakage_clean.R` — fits BIDIFAC+ on n=494
   training-fold samples only, projects n=118 test samples onto the
   training-derived module loadings via least-squares
   (V_test = X_test^T L (L^T L)^-1). Variance prefilter and row-mean
   imputation use training-fold statistics only.
2. `screening_paper5_v2.R` — ANOVA on PIPs from a Gibbs fit using
   training rows of the leakage-clean components. Top-3 covariates by
   η² re-locked at `reference/paper5_target_covariates_clean.csv`.
3. Both v2 build script and v2 Gibbs runner load the clean artifacts.

**Open verdict:** the units-fixed +10.87 nat Δ on contaminated chains
mixes leakage-driven gain with any real structural-prior effect. The
leakage-clean rerun (in progress as of this entry) will resolve this.

### Pre-reg constraints touched

- **§10 constraint 3** ("no covariate substitution"): deviated from. The
  original lock at SHA `34db626` was on a contaminated screening; treating
  it as binding would have locked in the leakage. New lock at
  `paper5_target_covariates_clean.csv`, committed before clean Gibbs
  fires.
- **§3 sensitivity definition**: clarified — sensitivity = L1 4-subgroup
  identity hierarchy, materially distinct from primary, per the original
  pre-reg intent. v1 implementation accidentally collapsed the two.
- **§4 decision rule** (±2 nats, CI excludes zero): unchanged. The new
  scale is +0.09/pt vs Lock's +0.04/pt, no longer anomalous; the ±2 nat
  threshold is now a reasonable check rather than trivially exceeded.

### Audit trail

Files written before clean rerun (preserved for review, not used in v2):
- `results/loglik_summary.csv` — v1 buggy (units + leakage)
- `results/loglik_summary_logfix.csv` — v1 chains, units fix only
- `results/gibbs_<spec>/chain_[1-4].rda` — v1 contaminated chains

Clean rerun outputs (new):
- `data/bidifac_components_clean.rds`
- `reference/paper5_target_covariates_clean.csv`
- `results/gibbs_<spec>_clean/chain_[1-4].rda`
- `results/loglik_summary_clean.csv`
