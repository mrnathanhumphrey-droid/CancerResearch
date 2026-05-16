# Pre-Registration — Paper 5: Within-Cancer Structural-Prior Methodology Transfer to Medulloblastoma

**Date drafted:** 2026-05-15
**To be locked:** before any Paper 5 Gibbs run (target: GitHub commit, public, timestamped).
**Filed alongside:** [../PRE_REGISTRATION.md](../PRE_REGISTRATION.md) (Paper 1), [../PRE_REGISTRATION_PAPER2.md](../PRE_REGISTRATION_PAPER2.md) (Paper 2 — falsifies), [../PRE_REGISTRATION_PAPER3_MOFAFLEX_NICHE.md](../PRE_REGISTRATION_PAPER3_MOFAFLEX_NICHE.md) (Paper 3 — falsifies on substrate), [../PRE_REGISTRATION_PAPER4_PHASE1.md](../PRE_REGISTRATION_PAPER4_PHASE1.md) (Paper 4 Phase 1 — hierarchy-specific).

---

## 1. Hypothesis under test

> The Lock 2022 hierarchical-spike-and-slab structural-prior methodology that delivered Paper 1's "matches with parameter reduction" verdict on a 29-cancer pan-cancer cohort (with epithelial-classification and embryonic-germ-layer as the structural-prior hierarchies) transfers, with comparable or improved effect, to a single-cancer corpus where the structural-prior hierarchy is the **within-cancer molecular subtype taxonomy**.

The substrate is medulloblastoma. The taxonomy is the Cavalli 2017 12-subtype molecular classification (with the WHO-consensus 4-subgroup as a coarser granularity for the sweep). The methodology is otherwise identical to Paper 1.

This is a within-cancer hierarchy test, **complementary to Paper 4's between-hierarchy generalization test** (Paper 4 held the 29-cancer substrate fixed and varied the hierarchy across three alternatives; Paper 5 holds the structural-prior mechanism fixed and varies the substrate to a single cancer with a hierarchically-organized molecular taxonomy).

---

## 2. Pre-specified substrate

**Primary substrate:** Cavalli et al. 2017, *Cancer Cell* 31:737–754. n = 763 primary medulloblastomas with matched DNA methylation (Illumina Infinium HumanMethylation450) and gene expression (Affymetrix Human Gene 1.1 ST Array).

| GEO accession | platform | n | content |
|---|---|---|---|
| **GSE85212** | Illumina Infinium HumanMethylation450 | 763 | DNA methylation β-values |
| **GSE85217** | Affymetrix Human Gene 1.1 ST Array | 763 | gene expression (RMA-normalized) |

Sample-level matching is 1:1 between the two accessions. Clinical and survival data are reported in Cavalli 2017 supplementary table S1 with completeness as follows (per the paper): age 95.7%, histology 76.9%, metastatic status 75.2%, **survival 82%** (≈626 patients with survival info; the back-of-envelope event total ≈ 170–210 deaths assuming pan-MB 5-yr OS ≈ 70–75%).

Confirmed before Paper 5 Gibbs fire: exact per-subgroup and per-12-subtype with-survival n and event counts go into `reference/cavalli_subtype_distribution.csv` produced by [scripts/build_cavalli_metadata.R](scripts/build_cavalli_metadata.R) on first data pull. The pre-registration is locked on the methodology and hierarchy; the data-pull step is a deterministic resolution of TBD numerical values, not a free parameter.

**External validation cohort (deferred):** Sharma et al. 2019, *Acta Neuropathologica* 138:309–326. n = 852 with matched DNA methylation + transcriptome (subset of 1,501 methylomes). Data deposit not fully resolved at pre-reg time; either confirmed via the paper's "Data availability" section before Paper 5 ships, contacted via the corresponding author, or deferred to a follow-up phase. Sharma's 8-subtype refinement within Group 3 / Group 4 is the level-3 granularity in §3 below.

**MAGIC / St. Jude Integrative Portal (n = 898, ACNS0331 / SJMB03 / ACNS0332 pooled):** not the primary substrate because raw genomic data is dbGaP-controlled and methylation + RNA-seq overlap is not the full 898. If used in a later phase, **trials must be pooled with `trial ∈ {ACNS0331, SJMB03, ACNS0332}` as a covariate**, and the trial × subgroup interaction must be reported as a sensitivity check (the 2022 MAGIC integrated analysis found no cross-trial EFS/PFS difference within matched dose strata — strong evidence the interaction is small, but not zero).

---

## 3. Pre-specified hierarchy granularities (sweep)

Three nested granularities, all derivable from Cavalli 2017's published per-sample labels:

| level | name | n bins | source | bin labels (TBD names confirmed from supp. table S1) |
|---|---|---|---|---|
| **L1** | 4-subgroup (WHO consensus) | 4 | Taylor 2012 consensus | {WNT, SHH, Group3, Group4} |
| **L2** | 12-subtype Cavalli | 12 | Cavalli 2017 SNF on methylation + expression | {WNT-α, WNT-β, SHH-α, SHH-β, SHH-γ, SHH-δ, G3-α, G3-β, G3-γ, G4-α, G4-β, G4-γ} (Greek-letter labels per Cavalli; exact label set confirmed from GSE85217 series_matrix metadata + Cavalli supplementary table S1 before Gibbs fire) |
| **L3** | 8-subtype within G3/G4 (Sharma 2019) | 8 within G3/G4 (still 4 in WNT/SHH) | Sharma 2019 meta-analysis | {I, II, III, IV, V, VI, VII, VIII} for G3/G4; WNT/SHH retained at L1 |

**Primary granularity for Paper 5: L2 (12-subtype Cavalli).** Locked.

Rationale: L1 (4-subgroup) leaves WNT as a singleton — n = 71 with ~4 events. The structural-prior mechanism requires the constrained subgroup to share information with biologically-related neighbors via a group-shared β. At L1, WNT has no neighbors → structural prior degenerates for WNT. At L2, WNT-α and WNT-β are siblings sharing a parental WNT label → structural prior pools them. The Lock-2022 analog: L1 maps to a hierarchy that puts each cancer in its own group (degenerate, like 1c), while L2 maps to the granularity sweep that gave Paper 1 its match-with-parameter-reduction verdict.

**Secondary granularity for sensitivity: L1 (4-subgroup).** Reported as a labelled sensitivity row, the same way Paper 1's 4-cov spec is preserved as a sensitivity to the 3-cov primary post-disposition.

**L3 (8-subtype Sharma G3/G4 refinement):** deferred to external-validation phase. Mentioned for corpus continuity.

---

## 4. Pre-specified target covariates

**Covariates are derived from a two-block BIDIFAC+ decomposition** of (methylation 450k β-values, expression Affy Gene 1.1 ST) on Cavalli's 763 samples, run before any structural-prior Gibbs fit. This mirrors Lock 2022's setup, scaled down from 4 omic blocks (Park-Lock-Hoadley 2022 BIDIFAC+ on TCGA) to 2 blocks (the only continuous-valued omics universally available in Cavalli).

**Screening pass:** ANOVA of per-subtype PIPs against the L2 (12-subtype) hierarchy on all BIDIFAC+-derived covariates, run before target-covariate selection. Number of BIDIFAC+ covariates expected: 30–50 (smaller than Lock's 68 because 2 blocks not 4). Number of survivors at η² > 0.30 and p < 0.05: TBD before target-covariate selection.

**Pre-specified target covariate count for the joint structural-prior validation: 3.** Locked.

Rationale for K = 3: Paper 1 fixed 4 covariates pre-reg and retracted cov 10.1 post-disposition to 3. Paper 5 starts at 3 to bake in the lesson. **All 3 target covariates are selected by largest-effect-size from the screening pass, with no covariate selection contingent on intermediate results.** The 3 covariate IDs (BIDIFAC+ component indices) are saved in `reference/paper5_target_covariates.csv` at the moment of screening completion, before any structural-prior Gibbs fit. This file's commit timestamp is the pre-registration of the covariate selection.

**Per-covariate marginal-contribution diagnostic (locked alongside the joint test, not post-hoc — Paper 1's hard-won lesson):**

For each of the 3 target covariates, the following two checks are reported alongside the joint LPPD verdict:

- **CHECK 1 (per-subgroup LPPD contribution):** for each of the 4 subgroups (WNT, SHH, G3, G4), compute (4-cov LPPD if K = 4 had been chosen) − (3-cov LPPD with that covariate ablated). Verdict A if no subgroup crosses |1 nat| AND binomial test on the 4 directional contributions is at chance; Verdict B if ≥ 1 subgroup crosses |1 nat| AND binomial p < 0.05.
- **CHECK 2 (patient-level bootstrap of the per-covariate marginal):** B = 1000 resamples, paired on the same per-patient LPPDs. Verdict A if 95% CI brackets zero; Verdict B if CI excludes zero in the direction of contribution.

If both CHECK 1 and CHECK 2 return Verdict A for a target covariate, that covariate is flagged as a noise-tier carrier and the corresponding ablated spec (K = 2) is reported as the headline. If both return Verdict B, the covariate is confirmed as load-bearing. Mixed verdicts trigger a conservative fallback to the ablated spec.

This is the same CHECK 1 / CHECK 2 cascade that drove Paper 1's 4-cov → 3-cov disposition, **committed at pre-reg time rather than after the joint result lands.**

---

## 5. Held-out evaluation protocol

- **Split:** stratified random 20% per L2 12-subtype, seed `20260515`. Ensures every subtype is represented in both training and held-out. Stratification at L2 (not L1) because Paper 5's primary granularity is L2; per-subtype balance is the relevant constraint.
- **Train / held-out counts:** ≈ 80% / 20% of with-survival n. Exact counts written to `reference/paper5_split_indices.rda` at split time; reported in §10 deliverables.
- **Baseline-on-train:** Cavalli's 763 samples fit with Lock-2022's unmodified `HierarchicalLogNormalSpikeSlab` on the 80% training subset, treating each L2 12-subtype as a "cancer" level (analog to Lock's 29 TCGA cancers).
- **Structural-prior fits:** the same sampler with the L2 (primary) or L1 (sensitivity) hierarchy injected via the modified Gibbs override step from Paper 1 (`scripts/validation/sampler_structural_prior.R` reused).

---

## 6. Gibbs run parameters

**Identical to Paper 1 with seed offsets to avoid collision with Papers 1–4.**

- Chains: 4 per spec, seeds 20260515, 20260516, 20260517, 20260518.
- Iterations: 100,000 per chain.
- Burn-in: 50,000.
- Priors: identical to Paper 1 (Inv-Gamma on σ², Normal on β̃, Bernoulli on γ).
- `pi_generation`: `"shared_across_cancers"` (where "cancer" here means L2 subtype).
- Posterior thinning for held-out evaluation: every 400 post-burn iters → 500 posterior draws per model.
- Specifications fit: (a) baseline-on-train (no structural prior); (b) 3-cov primary (L2 hierarchy); (c) 3-cov sensitivity (L1 hierarchy).

---

## 7. Predictive metric

**Identical to Paper 1.** Held-out log-likelihood (LPPD), sum over held-out patients, under the log-normal accelerated-failure-time spike-and-slab likelihood (right-censored).

---

## 8. Decision rule per specification

Bootstrap 95% CI on held-out LPPD difference per spec vs baseline-on-train, B = 1000 patient-level resamples, paired indices.

| outcome | conditions | per-spec disposition |
|---|---|---|
| **matches** baseline | observed Δ ∈ ±2 nats AND CI brackets zero | structural prior is methodologically inert at the joint level; per-covariate CHECK 1/CHECK 2 still reported |
| **outperforms** baseline | observed Δ > +2 nats AND CI excludes zero (positive) | structural prior generalizes to within-cancer subtype taxonomy on this substrate |
| **underperforms** baseline | observed Δ < −2 nats OR CI excludes zero (negative) | structural prior degrades performance; within-cancer transfer fails on this substrate |

**Joint Paper 5 disposition:**

| joint pattern | overall reading |
|---|---|
| primary (L2) outperforms | strong: within-cancer hierarchical structure carries the Paper-1 signal at 12-subtype granularity |
| primary (L2) matches AND parameter reduction ≥ 60% | matches-with-parameter-reduction (analog to Paper 1's post-disposition headline) |
| primary (L2) matches AND parameter reduction < 60% | no clean methodology lift; report as inconclusive |
| primary (L2) underperforms | within-cancer transfer fails; substrate-boundary finding |
| any sensitivity (L1) row contradicts primary | flag interaction between hierarchy granularity and structural-prior efficacy |

**Decision rule for L2 primary's parameter reduction:** free per-subtype βs for the 3 target covariates drop from 3 × 12 = 36 to 3 × (number of bins at L1) = 3 × 4 = 12, a 67% reduction. (Contested-cell handling: at L2, no contested cells expected if Cavalli's labels are complete; if any L2 subtype is unrepresented in the training fold after splitting, that subtype's PIP is forced to 0 for the structural-prior fit and noted.)

---

## 9. Convergence halt rule

**Identical to Paper 1.** R-hat > 1.05 on any monitored parameter halts that specification and reports the convergence failure. The other specifications proceed independently. No re-running with adjusted priors or alternative starts.

---

## 10. Pre-registered constraints

Locked before any Paper 5 Gibbs run:

1. **Substrate is Cavalli n = 763 only.** MAGIC, Sharma, and any other corpus are deferred to follow-up phases with their own pre-registrations.
2. **Hierarchy granularities are L1 (4-subgroup) and L2 (12-subtype Cavalli) only.** L3 (Sharma) is deferred.
3. **Three target covariates, selected by largest-effect-size from a BIDIFAC+ screening pass against L2.** Once committed to `reference/paper5_target_covariates.csv`, no covariate substitution.
4. **One held-out split (seed 20260515, 80/20 stratified at L2).** No k-fold expansion at this phase.
5. **No granularity tuning based on intermediate results.** L2 primary / L1 sensitivity is the spec list.
6. **No decision-rule adjustment.** ±2 nat threshold and CI-excludes-zero requirement committed pre-fit.
7. **All specifications run regardless of any one's outcome.**
8. **Per-covariate CHECK 1 + CHECK 2 cascade committed at pre-reg time**, reported alongside the joint test, not post-hoc.
9. **No re-running with adjusted priors or alternative starts.** Convergence checked once per spec; halt rather than re-fit.

---

## 11. What is reported regardless of outcome

- Held-out total LPPD for each of: baseline-on-train, L2 primary, L1 sensitivity.
- Per-subtype LPPD breakdown (12 rows at L2; 4 rows at L1).
- Bootstrap CI on three pairwise differences: each spec vs baseline-on-train.
- Bootstrap CI on the primary-vs-sensitivity comparison.
- **Per-covariate CHECK 1 + CHECK 2 cascade results** for each of the 3 target covariates.
- Convergence diagnostics for all fits (R-hat, bulk-ESS, tail-ESS, per-chain mean cross-check, σ² posterior agreement).
- Wall-clock per chain.
- Per-spec verdict per the §8 decision rule.
- Joint Paper 5 disposition.
- Free-parameter count before and after structural prior (analog to Paper 1's 87 → 22).
- BIDIFAC+ screening output (all candidate covariates × per-subtype PIP) — full table at `results/paper5_screening_results.csv`.

---

## 12. Provenance

Data sources, all public:

| source | platform | accession / DOI | role |
|---|---|---|---|
| Cavalli et al. 2017 — expression | Affy Hu Gene 1.1 ST | GEO **GSE85217** | primary expression block |
| Cavalli et al. 2017 — methylation | Illumina HumanMethylation450 | GEO **GSE85212** | primary methylation block |
| Cavalli et al. 2017 — clinical + 12-subtype labels + OS | supplementary table S1 (mmc2.xlsx) | DOI 10.1016/j.ccell.2017.05.005; Elsevier CDN `https://ars.els-cdn.com/content/image/1-s2.0-S1535610817302015-mmc2.xlsx` (fetched 2026-05-15; SHA-256 recorded by fetch_cavalli_data.R) | per-sample subtype assignment, survival |

Verified Cavalli mmc2 column mapping (header on row 2; 763 sample rows starting row 3):

| mmc2 column | normalized name | type | role |
|---|---|---|---|
| `Study_ID` | `sample_id` | char | join key to GEO `Sample_title` |
| `Age` | `age_at_dx` | numeric | age at diagnosis (years) |
| `Gender` | `sex` | char | M / F |
| `histology` | `histology` | char | Classic / Desmoplastic / Large-cell-anaplastic / etc. |
| `Met status (1 Met, 0 M0)` | `metastatic` | int | 1 = metastatic, 0 = M0, NA |
| `Dead` | `os_event` | int | 1 = death, 0 = censored |
| `OS (years)` | `os_time_years` | numeric | overall survival in years |
| `Subgroup` | `subgroup` | char | L1 4-bin: WNT / SHH / Group3 / Group4 |
| `Subtype` | `subtype` | char | L2 12-bin: {WNT,SHH,Group3,Group4}_{alpha,beta,gamma,delta} |
| Sharma et al. 2019 — G3/G4 8-subtype labels (deferred) | meta-analysis | DOI 10.1007/s00401-019-02020-0 | L3 granularity, external-validation phase |
| Taylor et al. 2012 — 4-subgroup consensus (reference only) | consensus position paper | DOI 10.1007/s00401-011-0922-z | L1 granularity definition |

Build scripts (to be added under `scripts/`):

- [scripts/fetch_cavalli_data.R](scripts/fetch_cavalli_data.R) — downloads GSE85217 + GSE85212 normalized matrices + supplementary clinical table; writes `data/cavalli_expr.rds`, `data/cavalli_meth.rds`, `data/cavalli_clinical.rds`.
- [scripts/build_cavalli_metadata.R](scripts/build_cavalli_metadata.R) — emits `reference/cavalli_subtype_distribution.csv` (per-subgroup + per-subtype n, with-survival n, event counts, median follow-up).
- [scripts/run_bidifac_plus.R](scripts/run_bidifac_plus.R) — two-block BIDIFAC+ on methylation + expression; writes BIDIFAC+ components.
- [scripts/screening_paper5.R](scripts/screening_paper5.R) — ANOVA of per-subtype PIPs against L2; emits `results/paper5_screening_results.csv` and `reference/paper5_target_covariates.csv` (top 3 by effect size).
- [scripts/build_paper5_split.R](scripts/build_paper5_split.R) — stratified 80/20 split at L2, seed 20260515; writes `reference/paper5_split_indices.rda`.
- [scripts/run_gibbs_paper5.R](scripts/run_gibbs_paper5.R) — Gibbs runner per spec × chain (reuses Paper 1's `sampler_structural_prior.R`).
- [scripts/compute_held_out_loglik_paper5.R](scripts/compute_held_out_loglik_paper5.R) — held-out LPPD per spec.
- [scripts/check1_paper5.R](scripts/check1_paper5.R) — per-subgroup LPPD contribution per target covariate.
- [scripts/check2_paper5.R](scripts/check2_paper5.R) — patient-level bootstrap of per-covariate marginal.
- [scripts/bootstrap_paper5.R](scripts/bootstrap_paper5.R) — B = 1000 paired bootstrap on joint LPPD differences.
- [scripts/convergence_paper5.R](scripts/convergence_paper5.R) — R-hat + ESS per spec; halt-or-pass.

This pre-registration is locked at commit SHA: **TBD on GitHub push**. The GitHub commit timestamp is the methods-section citation.

---

## 13. Event-count power assumption (locked at pre-reg time)

Back-of-envelope (to be replaced by `reference/cavalli_subtype_distribution.csv` exact numbers before Gibbs fire):

| subgroup | Cavalli n | 5-yr OS (literature) | est. events (Cavalli with-survival cohort) |
|---|---|---|---|
| WNT | 71 | 91–95% | ~4 |
| SHH | 233 | 65–85% | ~45–55 |
| Group 3 | 144 | 41–61% | ~55–70 |
| Group 4 | 326 | 70–80% | ~70–80 |
| **total (Cavalli with-survival ≈ 626)** | 763 | — | **~170–210** |

The structural-prior mechanism is expected to be most beneficial for the low-event subgroups (WNT specifically), via parameter-pooling through the L2 12-subtype hierarchy. This is the analog to Paper 1's pattern, where low-n cancers like KIRP and THCA gained most from group-shared βs. WNT-α and WNT-β as L2 siblings sharing a β is the load-bearing structural element that makes the methodology applicable here.

If the actual event count from `cavalli_subtype_distribution.csv` is ≥ 50% below the estimate range above (i.e., < 85 total events), the methodology power is too low to expect a clean OUTPERFORMS verdict; the spec proceeds anyway but the §8 disposition pre-commits to "underpowered" framing if Δ falls within ±2 nats. **This power-aware framing clause is the pre-registration's only event-count-contingent rule.**

---

## 14. What is explicitly not in this pre-registration

- Generalization to other within-cancer subtype taxonomies (DLBCL LymphGen, glioma WHO 2021 IDH/1p19q hierarchy, breast PAM50, prostate grade group, etc.) — separate pre-registrations per substrate.
- Generalization to all BIDIFAC+ components beyond the 3 target covariates — not claimed.
- Survival endpoint beyond overall survival (OS). PFS is not the primary endpoint; if used, separate pre-reg.
- K-fold cross-validation — single train/held-out split only at this phase.
- MAGIC (n = 898) or Sharma (n = 852) as primary substrates — deferred.
- Subgroup-specific structural-prior performance reporting (per-subgroup OUTPERFORMS vs MATCHES verdicts) — these are diagnostic, not headline. Joint Paper 5 disposition at §8 is the load-bearing claim.
