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

## 3. Verdict — Reading A

The bootstrap CI on (3-cov − 4-cov) includes zero, so dropping cov 10.1 does
**not** significantly hurt held-out performance. Reading A holds.

But the verdict has a caveat with paper-defensibility consequences:

- Cov 10.1's contribution to LPPD is statistically indistinguishable from
  zero (CI on its contribution: roughly [−4.40, +1.65] nats, point estimate
  +1.42 nats, p ≈ 0.16 in the bootstrap one-sided test).
- However, **without** cov 10.1, the 3-cov spec's improvement over baseline
  (+4.84 nats) has 95% CI [−0.82, +10.05] that **includes zero**. The
  3-cov spec does **not** meet the original Paper 1 decision rule
  (Δ > +2 *and* CI excludes zero).

So: cov 10.1 was nudging the 4-cov result across the pre-registered
"outperforms" threshold. Without it, the rest of the structural-prior
treatment passes the +2-nat point-estimate threshold but not the
CI-excludes-zero threshold.

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

## 5. Two honest options for Paper 1

### Option (a): Keep the 4-covariate spec, document the caveat

Paper 1's "outperforms" claim stands at Δ=+6.26 nats with CI excluding zero,
but the methods explicitly note:

- Cov 10.1's individual contribution is borderline (point estimate +1.42 nats,
  bootstrap CI on its removal [−1.65, +4.40], includes zero).
- Without cov 10.1, the spec is "matches" not "outperforms" under the
  original decision rule.
- The 4-covariate result is therefore at the edge of the decision rule's
  resolution rather than well above it.

This is the more conservative paper claim with full transparency about the
limit-of-the-test situation.

### Option (b): Drop cov 10.1, paper becomes "matches"

The headline becomes: 3-covariate structural priors (cov 8.2 / 1.1 / 20.1)
match the unconstrained baseline on held-out predictive performance with a
73%-equivalent parameter reduction (3 × 29 = 87 free betas → 3 × {3, 10, 8} =
21 group-shared betas).

This is a narrower methodology contribution: the structural priors don't
*hurt* held-out performance, even with the parameter reduction. The
"outperforms" claim is dropped; the paper makes the "matches" claim instead.

---

## 6. Recommendation

The data supports either option. The choice is a paper-defensibility call,
not a statistical one.

- **Option (a)** is defensible if the paper foregrounds the bootstrap CI's
  limit-of-resolution character honestly. Risks: a careful reviewer may
  zero in on the cov-10.1 caveat and require the 3-cov sensitivity test be
  the headline.
- **Option (b)** is more conservative and sidesteps that risk entirely. The
  scientific contribution is narrower but the claim is rock-solid.

The 3-covariate Gibbs run, its convergence diagnostics, and its bootstrap CI
are deposited in this repository so reviewers can evaluate either framing
directly.

---

## 7. Files

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
