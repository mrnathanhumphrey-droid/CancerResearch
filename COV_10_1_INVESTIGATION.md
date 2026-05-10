# Cov 10.1 Investigation — Reading A vs Reading B Test
## 2026-05-10 — follow-up to the wrongly-grouped diagnostic

After `WRONGLY_GROUPED_DIAGNOSTIC.md` flagged cov 10.1 as having near-random
sign agreement (0.48 at epi/2a; 0.57 at germ/3b), an additional Gibbs run was
fired to test whether cov 10.1 was contributing meaningful signal to Paper 1's
primary specification or fitting noise.

---

## 1. Test design

Rerun the primary validation with **three** of the four pre-specified target
covariates — drop cov 10.1 — holding everything else identical:

- 4 chains × 100,000 iters
- Same training / test split (seed 42)
- Same priors, same hyperparameters
- Three remaining target covariates: 8.2 / epi/2a, 1.1 / epi/2c, 20.1 / germ/3b
- Same modified Gibbs sampler (`scripts/validation/sampler_structural_prior.R`)

Compare held-out LPPD vs the 4-covariate primary spec and the unconstrained
baseline.

**Reading A:** if Δ (3-cov − 4-cov) ≈ 0 within bootstrap noise → cov 10.1
wasn't contributing meaningful signal; the +6.26-nat improvement came from
the other three. Cov 10.1 should be reconsidered in Paper 1's headline.

**Reading B:** if Δ (3-cov − 4-cov) is meaningfully negative → cov 10.1 was
contributing real signal despite the near-random within-group sign agreement.
The 4-covariate spec stands as-is.

---

## 2. Result

| spec | held-out total LPPD | Δ vs baseline | bootstrap 95% CI | CI excludes 0? |
|---|---|---|---|---|
| baseline (per-cancer SS, 80% train) | −1038.69 | — | — | — |
| 4-cov primary | −1032.43 | **+6.26** | [+0.34, +12.29] | **YES** |
| 3-cov primary (cov 10.1 dropped) | −1033.85 | +4.84 | [−0.82, +10.05] | NO |

**Δ (3-cov − 4-cov) = −1.42 nats**, bootstrap 95% CI **[−4.40, +1.65]**.
The CI on the difference includes zero.

### Convergence (3-cov primary)

R-hat max 1.0047 (cleaner than Paper 1's primary 1.0110). Bulk-ESS min 520,
tail-ESS min 1846. Halt rule not triggered.

---

## 3. Result — Reading A holds, with structural detail

The bootstrap CI on (3-cov − 4-cov) includes zero (point estimate −1.42 nats,
CI [−4.40, +1.65]). Per the test design, this is **Reading A**: cov 10.1's
contribution is statistically indistinguishable from zero at this n.

The structural detail beyond the binary verdict:

- Cov 10.1's bootstrap one-sided p (proportion of bootstrap samples where
  3-cov ≥ 4-cov) ≈ 0.16. Not significant; not vanishingly small either.
- Without cov 10.1, the 3-cov spec's improvement over baseline (+4.84 nats)
  has 95% CI [−0.82, +10.05] that **includes zero**. The 3-cov spec does
  **not** independently meet the original Paper 1 decision rule (Δ > +2
  *and* CI excludes zero).
- So the 4-cov spec is passing the decision rule with cov 10.1 contributing
  enough point-estimate (+1.42) and tightening (CI lower bound shifts from
  −0.82 to +0.34) to clear the threshold. Remove cov 10.1 and the rest of
  the model passes the +2-nat point-estimate criterion but not the
  CI-excludes-zero criterion.

---

## 4. Per-cancer breakdown (top movers)

The cancers most affected by removing cov 10.1 (positive = 3-cov outperforms
4-cov for that cancer; negative = 3-cov worse):

| cancer | Δ_4cov | Δ_3cov | Δ_3cov − Δ_4cov |
|---|---|---|---|
| LIHC | +1.47 | +0.57 | **−0.90** |
| LUSC | +1.43 | +0.81 | **−0.63** |
| SARC | +1.19 | +0.65 | **−0.54** |
| OV | +0.75 | +0.39 | −0.36 |
| BLCA | +0.23 | −0.09 | −0.33 |
| HNSC | +0.53 | +0.37 | −0.17 |
| ... | ... | ... | ... |
| KIRC | +0.42 | +0.56 | +0.14 |
| BRCA | +0.21 | +0.34 | +0.13 |
| SKCM | −0.02 | +0.20 | +0.22 |
| CHOL | +0.08 | +0.34 | +0.27 |
| KIRP | −1.14 | −0.87 | +0.27 |
| STAD | −0.84 | −0.40 | +0.44 |

Cov 10.1 was contributing signal in LIHC / LUSC / SARC / OV / BLCA / HNSC —
six cancers concentrated in lung-pleura, hepatobiliary, gynecological, soft
tissue, and head-neck. It was actively hurting fits in CHOL / KIRP / STAD /
SKCM. The aggregate balance is roughly zero across the 29 cancers, but the
per-cancer redistribution is real.

---

## 5. Statistical summary

The data without disposition framing:

- **The 4-cov primary spec's improvement over baseline** is +6.26 nats with
  bootstrap 95% CI [+0.34, +12.29] — point estimate well above the +2-nat
  threshold; CI excludes zero only by ≈0.3 nats at the lower bound.
- **Removing cov 10.1** drops the point estimate by 1.42 nats. The CI on
  this drop ([−4.40, +1.65]) includes zero — the per-test difference is
  not statistically distinguishable from zero at B=1000.
- **The 3-cov spec's improvement over baseline** is +4.84 nats with CI
  [−0.82, +10.05] — point estimate above +2 but CI now overlaps zero at
  the lower bound by ≈0.8 nats.

Three numerical relationships matter:

1. The 4-cov spec is at the lower edge of the CI-excludes-zero criterion
   (lower bound +0.34).
2. Cov 10.1's individual contribution is in the noise band (CI includes
   zero, point estimate +1.42).
3. Removing cov 10.1 moves the spec from "CI excludes zero" to "CI includes
   zero" — i.e., cov 10.1 is the marginal covariate by which the 4-cov
   spec passes its decision rule.

Per-cancer redistribution: cov 10.1 contributes positive LPPD in LIHC, LUSC,
SARC, OV, BLCA, HNSC; negative LPPD in CHOL, KIRP, STAD. Aggregate is the
+1.42-nat point estimate; redistribution within is real.

The 3-covariate Gibbs run, its convergence diagnostics, and its bootstrap
CI are deposited in this repository for downstream reviewers and
disposition decisions.

---

## 6. Files

| file | contents |
|---|---|
| `scripts/validation/01b_build_primary_3cov.R` | constructs `target_covs_primary_3cov.rda` (drops cov 10.1) |
| `scripts/validation/02b_compute_loglik_3cov.R` | LPPD + bootstrap CI for the 3-cov spec |
| `scripts/validation/05b_convergence_3cov.R` | convergence diagnostics for the 3-cov chains |
| `results/validation/loglik_3cov_summary.csv` | total LPPD per spec |
| `results/validation/per_cancer_loglik_3cov_comparison.csv` | per-cancer breakdown |
| `results/validation/loglik_per_patient_primary_3cov.csv` | per-test-patient LPPD for 3-cov spec |

(Chain RDAs at `gibbs_primary_3cov/chain_[1-4].rda` are not in the public
repo; they will be deposited at Zenodo alongside the other Paper 1 chain
RDAs.)
