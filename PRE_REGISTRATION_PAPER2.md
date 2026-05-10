# Pre-Registration — Operator-Composition Extension (Paper 2)

**Locked:** 2026-05-10. **Compute deferred** pending three prerequisites
(see §12). This document is committed to the public repository before any
extension Gibbs run fires. The defensibility of the extension claim
depends on this commit predating compute.

**Companion:** `PRE_REGISTRATION.md` (Paper 1, locked 2026-05-09); the
Paper 1 result is the comparison anchor for this extension.

---

## 1. Hypothesis under test

> Structural priors at single (hierarchy, granularity) cells operate as
> block-projection operators on the 29-cancer state space. The methodology
> may extend to operator **compositions** — combinations of structural
> priors that, applied jointly, recover patterns single-operator priors
> miss.

Paper 2 question: does composing structural-prior operators across
hierarchies improve held-out predictive performance beyond what the best
single-operator prior achieves?

Operationalized as: a Gibbs sampler in which the per-cancer `β_jk` and
`γ_jk` for a target covariate are determined jointly by *two* hierarchy
operators (P_h1, P_h2) under one of three pre-specified composition rules,
rather than by a single hierarchy operator as in Paper 1. The model is fit
on the same 80% training set used in Paper 1 and evaluated on the same
20% held-out test set (seed 42).

---

## 2. Pre-specified specifications (18 total, locked)

### 2.1 Three operator-composition rules

For target covariate `c`, let:
- P_h1(c) → group label per cancer under hierarchy/granularity `h1`
- P_h2(c) → group label per cancer under hierarchy/granularity `h2`
- β̂_g_k(h) → group-shared posterior mean for group `g` of hierarchy `h`,
  fit from pooled within-group residual data (the Paper 1 mechanism)

The three rules:

**RULE 1 — sum (equal-vote averaging).**
For each cancer `c`, the prior beta is:
$$
\beta_{c,k}^{\text{rule1}} = \tfrac{1}{2}\bigl(\hat\beta_{P_{h_1}(c), k} + \hat\beta_{P_{h_2}(c), k}\bigr)
$$
Each hierarchy contributes a group-mean estimate; the final per-cancer prior
is the equal-weight average. PIP is forced to 1 if the cancer is in
non-contested groups under *both* hierarchies; 0 otherwise.

**RULE 2 — sequential projection (hierarchical refinement).**
Apply h1 first to determine primary group. Within each h1-group, apply h2
as a sub-grouping. The per-cancer prior is the per-(h1, h2)-joint group
mean:
$$
\beta_{c,k}^{\text{rule2}} = \hat\beta_{(P_{h_1}(c),\, P_{h_2}(c)),\, k}
$$
fit by pooling residual data only across cancers sharing both labels. PIP
is forced to 1 if the joint group has at least 2 cancers; 0 otherwise.

**RULE 3 — intersection-only with fallback.**
Same as RULE 2 when the joint (h1, h2) group has at least 2 cancers; for
singleton joint groups (cancer alone in its (h1, h2) cell), fall back to
the h1 group mean:
$$
\beta_{c,k}^{\text{rule3}} =
\begin{cases}
\hat\beta_{(P_{h_1}(c),\, P_{h_2}(c)),\, k} & \text{if joint group } |g| \geq 2 \\
\hat\beta_{P_{h_1}(c),\, k} & \text{otherwise}
\end{cases}
$$
PIP is forced to 1 if the cancer is in a non-contested h1-group; 0 otherwise.
Avoids overfitting on tiny joint cells while exploiting joint structure
where data supports it.

### 2.2 Three target covariates

Fixed from Paper 1's surviving four (no expansion):

| cov_id | rationale |
|---|---|
| **8.2**  | strongest Paper 1 result; unambiguous single-hierarchy winner (epi/2a; Δ=+6.26 nats) |
| **10.1** | Paper 1 partial survivor; single-hierarchy result at epi/2c |
| **1.1**  | Paper 1 partial survivor; single-hierarchy result at epi/2c |

Cov 20.1 from Paper 1 is **not included** as a target covariate in this
pre-registration. It can be added in a future pre-registration only.

### 2.3 Two hierarchy pairs

| pair | h1 | h2 | rationale |
|---|---|---|---|
| **A** | epithelial_class / 2a_3bin | germ_layer / 3b_mid | the two hierarchies where Paper 1 tested cov 8.2 (primary vs. secondary specifications) |
| **B** | epithelial_class / 2c_10bin | tissue_origin / 1b_11bin | epi at fine granularity composed with anatomical tissue origin (the only hierarchy that escaped the BIDIFAC+ block-presence confound in the Paper 1 falsification cascade) |

Hierarchy assignments per cancer follow `reference/cancer_type_hierarchies_2026-05-09.md`.

### 2.4 The 18 specifications

```
3 covariates × 3 rules × 2 hierarchy pairs = 18 specifications.
```

| spec | covariate | rule | pair |
|---|---|---|---|
| 1 | 8.2 | 1 | A |
| 2 | 8.2 | 1 | B |
| 3 | 8.2 | 2 | A |
| 4 | 8.2 | 2 | B |
| 5 | 8.2 | 3 | A |
| 6 | 8.2 | 3 | B |
| 7 | 10.1 | 1 | A |
| 8 | 10.1 | 1 | B |
| 9 | 10.1 | 2 | A |
| 10 | 10.1 | 2 | B |
| 11 | 10.1 | 3 | A |
| 12 | 10.1 | 3 | B |
| 13 | 1.1 | 1 | A |
| 14 | 1.1 | 1 | B |
| 15 | 1.1 | 2 | A |
| 16 | 1.1 | 2 | B |
| 17 | 1.1 | 3 | A |
| 18 | 1.1 | 3 | B |

**All 18 specifications are reported regardless of result.** No additions,
no substitutions, no granularity tuning during compute.

---

## 3. Held-out evaluation protocol

Same as Paper 1: stratified random 20% per cancer, seed 42, 1,362 test
patients, 5,386 training patients. The exact same patient-level split is
used so the LPPD comparisons against Paper 1 are direct.

---

## 4. Gibbs run parameters

Identical to Paper 1 baseline except for the target covariate's
operator-composition treatment.

- Chains: 4 per specification, seeds 20260509–20260512
- Iterations: 100,000 per chain
- Burn-in: 50,000
- Priors on retained covariates: identical to baseline
- `pi_generation`: `"shared_across_cancers"`
- Posterior thinning for held-out evaluation: every 400 post-burn iters,
  4 chains × 125 thinned samples = 500 posterior draws per spec

---

## 5. Predictive metric

Held-out log-likelihood (LPPD), identical to Paper 1 §6:

$$
\text{LPPD} = \sum_{i=1}^{n_{\text{test}}} \log \frac{1}{S}\sum_{s=1}^{S} p(Y_i \mid \theta_s)
$$

uncensored: log-normal density on log(Y); censored: log survival function.

---

## 6. Decision rule

Each of the 18 specifications is compared against:
1. The **Paper 1 single-operator result** for the same target covariate
   (the relevant primary or secondary specification of Paper 1 — see §6.1)
2. The **baseline replication** (`gibbs_baseline_train` from Paper 1)

### 6.1 Paper 1 anchors per target covariate

| covariate | Paper 1 anchor | held-out total LPPD |
|---|---|---|
| 8.2 | primary spec (epi/2a) | −1032.43 |
| 10.1 | primary spec (epi/2c, in primary multi-cov spec) | (component of −1032.43) |
| 1.1 | primary spec (epi/2c, in primary multi-cov spec) | (component of −1032.43) |

**Note on the anchor for 10.1 and 1.1:** Paper 1's primary specification
includes all four covariates jointly. The "Paper 1 single-operator anchor"
for cov 10.1 and cov 1.1 in this extension test is therefore the primary-spec
total LPPD (−1032.43). The extension question for these two covariates is
whether composing operators on the same single covariate yields LPPD better
than what the primary spec achieves with single-operator structural priors
on all four covariates. (Equivalently: does composing on cov 10.1 alone
recover predictive performance comparable to or better than composing
single-operator priors on the full Paper 1 four-covariate set?) This anchor
choice is conservative: it sets a higher bar for the extension claim than
isolating each covariate's individual single-operator contribution would.

### 6.2 Outcomes

For each specification, against the Paper 1 anchor:

| outcome | conditions |
|---|---|
| **improves on Paper 1** | observed Δ > +2 nats AND bootstrap p < 0.00278 (Bonferroni-corrected, see §7) |
| **matches Paper 1** | \|Δ\| ≤ 2 nats AND nominal bootstrap 95% CI includes zero |
| **underperforms Paper 1** | observed Δ < −2 nats OR nominal bootstrap CI excludes zero on the negative side |

Specifications that have observed Δ > +2 nats and nominal bootstrap CI
excluding zero but **fail** the Bonferroni threshold (0.00278) are reported
as **nominal-only** improvements and are NOT counted toward the
"validation" condition (§8).

---

## 7. Bonferroni multiple-testing correction

The 18 specifications are jointly tested. The Bonferroni-corrected
significance threshold is α = 0.05 / 18 = **0.00278**.

The bootstrap p-value for "improvement over Paper 1" is computed as:
$$
p_{\text{boot}} = \frac{|\{b : \Delta_b \leq 0\}|}{B},\quad B = 1000
$$
(the proportion of bootstrap samples in which the spec did *not* improve
over Paper 1). For an "improves on Paper 1" verdict, we require
`p_boot < 0.00278`.

The pre-registered nominal CI-excludes-zero criterion (Paper 1's rule) is
**also** computed and reported, but the validation/falsification
condition (§8) uses Bonferroni.

---

## 8. Validation / falsification criteria for the extension

For the operator-composition extension to be **validated**:

1. **At least 2 of the 18 specifications** must "improve on Paper 1"
   under the Bonferroni-corrected criterion.
2. **At least 1 of those improvements** must come from RULE 2 or RULE 3
   (sequential projection or intersection-only with fallback).

Condition 2 prevents validation from depending solely on RULE 1 (sum),
which could be artifactual regularization (averaging two priors → smaller
prior variance → tighter regularization → marginal LPPD improvement
without recovering distinct compositional structure).

If **either condition fails**, the extension is **falsified**: operator
composition does not generalize beyond single-operator structural priors
in this domain at this sample size. The result is documented as a
null-extension finding and the operator-composition framing is shelved at
the cancer-domain level.

---

## 9. Convergence halt rule

If R-hat > 1.05 on any monitored parameter for any specification: report
the convergence failure for that specification but **DO NOT replace** it
with a re-fit using different priors or starting values. The
pre-registered specification list is locked; the convergence failure is
part of the result.

A specification with R-hat > 1.05 is reported as **convergence-failed**
and excluded from the validation-condition counting in §8 (it cannot
"improve on Paper 1" with bad sampling). It is reported in the umbrella
table as part of the 18 results.

---

## 10. Pre-registered constraints (LOCKED)

1. **No specification added** beyond the pre-specified 18.
2. **No granularity tuning** based on intermediate or final results.
3. **No covariate substitution.** The three target covariates are fixed.
4. **No hierarchy substitution.** The two hierarchy pairs are fixed.
5. **No decision-rule adjustment.** The +2-nat / Bonferroni-0.00278
   thresholds and the validation conditions in §8 are committed before
   any extension compute.
6. **All 18 specifications are reported equally.** No "promising-looking"
   results get extra documentation; no "ugly-looking" results get hidden
   in appendices.
7. **No follow-up tests pursued based on intermediate results.**
   If validation falsifies, the falsification is the result; it does not
   open the door to "let's try other compositions."
8. **The Paper 1 single-operator anchor is the comparison anchor.**
   Comparisons against the baseline replication are also reported but the
   validation condition uses the Paper 1 anchor.

---

## 11. What is reported regardless of outcome

For each of the 18 specifications:

- Held-out total LPPD
- Δ vs Paper 1 anchor (single-operator result for the same target covariate)
- Δ vs baseline replication (`gibbs_baseline_train`)
- Bootstrap 95% CI on each Δ (B=1000, patient-level resample with replacement)
- Bootstrap p-value for "Δ ≤ 0"
- Per-cancer LPPD breakdown (29 rows per spec)
- Convergence diagnostics: R-hat max, R-hat > 1.01 count, R-hat > 1.05
  count, bulk-ESS min, tail-ESS min, divergent-transition count
- Wall-clock per chain
- Verdict per the decision rule (improves / matches / underperforms /
  nominal-only / convergence-failed)
- Per-cancer LPPD difference vs Paper 1 anchor (which cancers improve, which worsen)

The umbrella deliverables:

| file | contents |
|---|---|
| `EXTENSION_VALIDATION_RESULTS.md` | full report with all 18 specifications + verdict |
| `extension_decision_log.md` | explicit decision-rule application per spec |
| `operator_composition_summary.csv` | tabular: spec, rule, pair, cov, Δ_paper1, CI, Bonferroni p, verdict |
| `gibbs_extension_<spec>/` | chain RDAs (Zenodo, not in git) |
| `loglik_per_patient_<spec>.csv` | per-test-patient LPPD |
| `bootstrap_extension_<spec>.csv` | B=1000 bootstrap diffs |
| `per_cancer_loglik_<spec>.csv` | per-cancer LPPD breakdown |
| `convergence_<spec>.csv` | per-spec convergence diagnostics |

---

## 12. Prerequisites for compute

**Compute does not begin until ALL three prerequisites are satisfied:**

(a) **Paper 1 has shipped.** Either:
   - Cancer Paper 1 manuscript sent to the original authors, OR
   - Public release (this GitHub repository pushed to remote + Zenodo deposit
     for the Paper 1 chain RDAs created)

(b) **Wrongly-grouped diagnostic is complete and documented.** This is
   a separate diagnostic — upstream prerequisite work, not part of this
   pre-registration. The diagnostic identifies whether the Paper 1 secondary
   specification's germ-layer grouping has systematic group-mixing (e.g.,
   LUSC/LUAD pulling opposite directions in the foregut endoderm group).
   The diagnostic results inform the *discussion* of this extension's
   findings but **do not modify the operator-composition specifications**.
   The pre-registration locks the test before the diagnostic informs the
   design.

(c) **This pre-registration document is committed** to the public
   repository with a timestamped commit, prior to any extension compute.

The pre-registration commit is **load-bearing**. The extension claim's
defensibility depends on the pre-registration being verifiable as
predating the compute. Without the timestamped commit, the result is
fishing regardless of how clean the validation looks.

---

## 13. Compute budget

18 Gibbs runs at Paper 1 cost level (~35 minutes wall per spec with
4 parallel chains). Sequential execution (no parallelism across specs to
avoid resource contention).

- **Sequential: ~10–11 hours wall** (18 × 35 min)
- **Per-chain parallel within spec: ~35 min wall per spec** (already true
  in Paper 1)

If compute hardware allows, parallelizing 2-3 specs concurrently is
acceptable provided memory and CPU contention does not degrade per-chain
performance below the Paper 1 baseline. The pre-registration does not
require strict sequential execution; it requires that all 18 specifications
complete before any analysis or interpretation.

---

## 14. Disposition rules (locked)

If the extension is **validated** (per §8):
- Cancer Paper 2 writes up the operator-composition framework + the
  validating specifications, with the falsified specifications honestly
  reported as null results.
- A cross-domain test (e.g., SP500) follows as a separate phase,
  pre-registered separately.

If the extension is **falsified** (per §8):
- Cancer Paper 2 does not exist as a methodology-extension paper.
- The result is documented as a null-extension finding.
- The operator-composition framing is shelved at the cancer-domain level.
- Any cross-domain test does not proceed under this framework.

There is no third option. The validation condition is binary and dispositive.

---

*Pre-registered 2026-05-10. Commit hash to be assigned by the next git
commit including this file. Compute deferred pending §12 prerequisites.*
