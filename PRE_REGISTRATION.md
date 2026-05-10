# Pre-Registration — Structural-Prior Validation

**Locked:** 2026-05-09, before any validation Gibbs run.
**Filed alongside:** the screening + falsification cascade summarized in
`FALSIFICATION_REPORT.md`. The decision rules below were committed in
writing prior to fitting any modified-Gibbs model.

---

## 1. Hypothesis under test

> Structural priors derived from epithelial classification and embryonic
> germ layer can replace the per-cancer fitted PIPs for four pre-specified
> covariates without loss of held-out predictive performance.

Operationalized as: a Gibbs sampler in which the per-cancer `β_jk` for the
four target covariates is replaced by a group-shared `β_gk` (group defined
by external hierarchy assignment, fitted from pooled within-group residual
data), with PIP forced to 1 within group and 0 outside. All other covariates
retain the original Lock spike-and-slab structure. The model is fit on an
80% training set and evaluated on a 20% held-out test set.

---

## 2. Pre-specified target covariates

Four covariates only; **no expansion based on validation results**.

| cov_id | primary specification | secondary specification |
|---|---|---|
| 10.1 | epithelial_class / 2c_10bin | (same) |
| 1.1  | epithelial_class / 2c_10bin | (same) |
| 8.2  | epithelial_class / 2a_3bin  | germ_layer / 3b_mid     |
| 20.1 | germ_layer / 3b_mid         | (same) |

Primary and secondary specifications differ only in the cov 8.2 hierarchy
(chosen because it had nominal hits in both epithelial_class/2a and
germ_layer/3b in the screening). The primary specification is the headline;
the secondary is sensitivity.

---

## 3. Hierarchy assignments and contested-cell handling

Group assignments per cancer follow `reference/cancer_type_hierarchies_2026-05-09.md`.
Cancers flagged contested in that table are **excluded per-hierarchy** for
the relevant target covariate: their PIP is forced to 0 (not in group).
Cancers in non-contested groups have PIP forced to 1, with their per-cancer
beta replaced by the group-shared sample.

Specifically:
- epithelial_class contested: MESO, UCS (2 of 29 cancers)
- germ_layer contested: BLCA, BRCA, CESC, HNSC, PRAD, TGCT, THYM (7 of 29)

These are pre-specified from the cancer-type hierarchy reference file and
not adjusted based on validation outcomes.

---

## 4. Held-out evaluation protocol

- **Split:** stratified random 20% of patients per cancer cohort, with seed
  `42`. Each cancer cohort retains both training and test representation.
- **Training set:** 5,386 patients (80% of all 6,748).
- **Test set:** 1,362 patients.
- **Test-set composition:** preserved across all three model fits — baseline
  on 80%, primary, secondary all evaluate on the same 1,362 patients.

---

## 5. Gibbs run parameters

Identical to the published Lock 2022 baseline except for the four target
covariates' structural-prior treatment.

- Chains: 4, seeds 20260509, 20260510, 20260511, 20260512
- Iterations: 100,000 per chain
- Burn-in: 50,000 (post-burn samples = 50,001:100,000 per chain)
- Priors on retained covariates: identical to baseline (`spike_priorvar = 1/10000`,
  `betatilde_priorvar_intercept = 100`, `betatilde_priorvar_coefficient = 1`,
  `lambda2_priorshape_intercept = 1`, `lambda2_priorrate_intercept = 1`,
  `lambda2_priorshape_coefficient = 5`, `lambda2_priorrate_coefficient = 1`,
  `sigma2_priorshape = sigma2_priorrate = 0.01`)
- `pi_generation`: `"shared_across_cancers"` (matches baseline)
- Posterior thinning for held-out evaluation: every 400 post-burn iterations,
  4 chains × 125 thinned samples = 500 posterior draws per model.

---

## 6. Predictive metric

**Held-out log-likelihood (LPPD):** for each test patient `i` in cancer `c`
with covariates `X_i`, observed survival `Y_i` (or right-censored at `C_i`),
and posterior sample `θ_s = (β_c^s, σ^{2,s})`:

$$
\log p(Y_i \mid \theta_s) = \begin{cases}
\log \mathcal{N}(\log Y_i; \mu_i^s, \sigma^s) & \text{uncensored} \\
\log \Phi\!\left(\frac{\mu_i^s - \log C_i}{\sigma^s}\right) & \text{censored}
\end{cases}
$$

with $\mu_i^s = X_i^\top \beta_c^s$. Pointwise predictive density:

$$
\text{LPPD}_i = \log \frac{1}{S} \sum_{s=1}^{S} p(Y_i \mid \theta_s)
$$

Total = `sum_i LPPD_i` over all 1,362 test patients.

---

## 7. Decision rule

Bootstrap 95% CI on the held-out LPPD difference between specifications,
B=1000 patient-level resamples (uniform with replacement from the test set).

A specification:

| outcome | conditions | paper disposition |
|---|---|---|
| **matches** baseline | observed Δ ∈ ±2 nats AND bootstrap 95% CI includes zero | narrow methodology contribution paper exists |
| **outperforms** baseline | observed Δ > +2 nats AND bootstrap 95% CI excludes zero (positive side) | stronger paper than expected |
| **underperforms** baseline | observed Δ < −2 nats OR bootstrap 95% CI excludes zero (negative side) | no cancer paper from this work; methodology didn't translate; document negative result |

The threshold is **±2 nats** in absolute log-likelihood and CI exclusion of
zero — both required for a non-null verdict.

The primary specification's verdict determines the paper disposition. The
secondary specification is reported as sensitivity.

---

## 8. Convergence halt rule

If any monitored parameter (β̃, σ², π — total 137 monitored params) has
**R-hat > 1.05** in any of the three runs (baseline-on-train, primary, or
secondary), the analysis halts and reports the convergence failure. No
interpretation is pushed through inadequate sampling.

R-hat is computed via the `posterior` R package using the standard
split-R̂ convention on 50,000 post-burn samples × 4 chains.

---

## 9. Pre-registered constraints

Locked before the validation Gibbs runs:

1. **No covariate added beyond the four pre-specified.** The temptation to
   "expand the structural-prior treatment to other covariates that look
   promising in the screening" is exactly the data-mining failure mode the
   falsification cascade was designed to prevent.

2. **No granularity tuning based on intermediate results.** The
   (covariate, hierarchy, granularity) assignments above are fixed.

3. **No decision-rule adjustment.** The ±2 nat threshold and CI-excludes-zero
   requirement are committed before any held-out evaluation.

4. **One held-out split.** A single 80/20 stratified split with seed 42; not
   k-fold. (k-fold acknowledged as a stronger design but out of scope for
   this pre-registration.)

5. **No re-running with adjusted priors or alternative starting values.**
   Convergence diagnostics are checked once per spec; if R-hat > 1.05,
   halt rather than re-fit.

6. **Both specifications run regardless of primary outcome.** Secondary is
   not contingent on primary's verdict.

---

## 10. What is reported regardless of outcome

- Held-out total LPPD for baseline-on-train, primary, secondary
- Per-cancer LPPD breakdown (29 rows)
- Bootstrap CI on three pairwise differences: primary vs baseline, secondary
  vs baseline, primary vs secondary
- Convergence diagnostics for all three fits (R-hat, bulk-ESS, tail-ESS,
  per-chain mean cross-check)
- Wall-clock per chain
- Verdict per the decision rule, applied to the primary specification

---

## 11. What was reported (post-hoc, for reference)

This pre-registration was filed before the validation Gibbs runs. The
outcomes are reported in `VALIDATION_RESULTS.md`. The decision rule fired:

> **OUTPERFORMS** for the primary specification.
> Δ = +6.26 nats vs baseline-on-train, bootstrap 95% CI [+0.34, +12.29].

Disposition: stronger paper than expected.

The secondary specification: Δ = +4.33 nats, 95% CI [−1.86, +10.97],
verdict **inconclusive** (same direction, CI overlaps zero).

These post-hoc outcomes are reported here for completeness and transparency;
the pre-registration itself was filed before the runs and the decision rule
above was applied unmodified.
