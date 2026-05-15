# Pre-Registration — Paper 4 / Phase 1: Structural-Prior Validation Across Three New External Hierarchies

**Date drafted:** 2026-05-15
**To be locked:** before any Phase 1 Gibbs run (target: GitHub commit, public, timestamped)
**Filed alongside:** [PRE_REGISTRATION.md](PRE_REGISTRATION.md) (Paper 1), [PRE_REGISTRATION_PAPER2.md](PRE_REGISTRATION_PAPER2.md) (Paper 2), [PRE_REGISTRATION_PAPER3_MOFAFLEX_NICHE.md](PRE_REGISTRATION_PAPER3_MOFAFLEX_NICHE.md) (Paper 3, falsified)

---

## 1. Hypothesis under test

> The Paper 1 +6.26 nat structural-prior result, obtained on the epithelial-class and germ-layer hierarchies, generalizes to other externally-published pan-cancer hierarchies of comparable cohort scale that are *not* derived from the BIDIFAC+ blocks (and therefore not subject to the circularity flagged for Hoadley supercluster).

Three independently-published TCGA pan-cancer hierarchies are tested as structural-prior replacements for the same four target covariates as Paper 1.

---

## 2. Pre-specified target covariates

**Same four covariates as Paper 1.** No covariate added, removed, or substituted.

| cov_id | Paper 1 hierarchy (for reference) |
|---|---|
| 10.1 | epithelial_class / 2c_10bin |
| 1.1 | epithelial_class / 2c_10bin |
| 8.2 | epithelial_class / 2a_3bin (primary) or germ_layer / 3b_mid (secondary) |
| 20.1 | germ_layer / 3b_mid |

Phase 1 replaces the Paper 1 hierarchies with each of three NEW hierarchies (separately, one per specification).

---

## 3. Three new hierarchies + cancer assignments

### Hierarchy A — Thorsson 2018 dominant immune subtype (C1-C6)
**Source:** Thorsson V. et al., "The Immune Landscape of Cancer," *Immunity* 48(4):812-830 (2018). DOI: 10.1016/j.immuni.2018.03.023. Per-tumor immune subtype assignments from supplementary table mmc2 (PanImmune_MS sheet, "Immune Subtype" column). Fetched directly from Elsevier CDN (`https://ars.els-cdn.com/content/image/1-s2.0-S1074761318301213-mmc2.xlsx`).
**Method of cancer-level assignment:** for each BIDIFAC+ cancer, the most-common immune subtype across primary tumors.
**Group structure:** {C1, C2, C3, C4, C5}, plus contested = {DLBC, THYM} (no Thorsson data).

### Hierarchy B — Malta 2018 stemness tertile (low / med / high)
**Source:** Malta T.M. et al., "Machine Learning Identifies Stemness Features Associated with Oncogenic Dedifferentiation," *Cell* 173(2):338-354 (2018). DOI: 10.1016/j.cell.2018.03.034. Per-tumor mDNAsi (DNA-methylation stemness index) from supplementary table mmc1 (StemnessScores_DNAmeth sheet, "mDNAsi" column). Fetched from Elsevier CDN (`https://ars.els-cdn.com/content/image/1-s2.0-S0092867418303581-mmc1.xlsx`).
**Method of cancer-level assignment:** per-cancer mean mDNAsi across primary tumors, then tertile-split across all 29 BIDIFAC+ cancers. Cutoffs (computed before Gibbs runs, on the 2026-05-15 build): **low ≤ 0.142, med ≤ 0.196, high > 0.196**.
**Group structure:** {low, med, high}. No contested.

### Hierarchy C — Sanchez-Vega 2018 dominant oncogenic pathway
**Source:** Sanchez-Vega F. et al., "Oncogenic Signaling Pathways in The Cancer Genome Atlas," *Cell* 173(2):321-337 (2018). DOI: 10.1016/j.cell.2018.03.035. Per-tumor binary pathway alteration matrix (10 pathways: Cell Cycle, HIPPO, MYC, NOTCH, NRF2, PI3K, RTK-RAS, TP53, TGF-Beta, WNT) from supplementary table S4 (Pathway level sheet). Fetched from Washington University Digital Commons (`https://digitalcommons.wustl.edu/cgi/viewcontent.cgi?filename=3&article=8671&context=open_access_pubs&type=additional`).
**Method of cancer-level assignment:** per-cancer alteration rate per pathway, then dominant = pathway with highest rate.
**Group structure:** {Cell Cycle, HIPPO, PI3K, RTK-RAS, TP53, WNT} (6 groups; MYC, NOTCH, NRF2, TGF-Beta not dominant in any BIDIFAC+ cancer).

### Combined cancer × hierarchy assignment table

| cancer | Thorsson | Stemness | Sanchez-Vega |
|---|---|---|---|
| ACC | C4 | low | WNT |
| BLCA | C1 | high | Cell Cycle |
| BRCA | C2 | med | PI3K |
| CESC | C2 | med | PI3K |
| CHOL | C3 | low | RTK-RAS |
| CORE | C1 | med | WNT |
| DLBC | _contested_ | high | Cell Cycle |
| ESCA | C2 | high | Cell Cycle |
| HNSC | C2 | high | Cell Cycle |
| KICH | C3 | med | TP53 |
| KIRC | C3 | low | Cell Cycle |
| KIRP | C3 | low | RTK-RAS |
| LGG | C5 | low | HIPPO |
| LIHC | C4 | high | Cell Cycle |
| LUAD | C3 | high | RTK-RAS |
| LUSC | C1 | high | Cell Cycle |
| MESO | C1 | low | Cell Cycle |
| OV | C2 | high | TP53 |
| PAAD | C1 | med | RTK-RAS |
| PCPG | C3 | high | RTK-RAS |
| PRAD | C3 | low | WNT |
| SARC | C1 | med | TP53 |
| SKCM | C1 | med | RTK-RAS |
| STAD | C2 | med | Cell Cycle |
| TGCT | C2 | high | RTK-RAS |
| THCA | C3 | low | RTK-RAS |
| THYM | _contested_ | low | RTK-RAS |
| UCEC | C1 | med | PI3K |
| UCS | C1 | med | TP53 |

Saved in canonical form at [`reference/phase1_hierarchy_assignments.csv`](reference/phase1_hierarchy_assignments.csv). The CSV is the authoritative source; the table above is for reading.

### Pre-specified contested-cell exclusions

Hierarchy A (Thorsson): **DLBC, THYM** excluded (no Thorsson data — these are not in the panimmune cohort). Their per-covariate PIP is forced to 0 (not in group) for the Phase 1 Hierarchy A specification only; they retain Paper 1-style spike-and-slab structure in the other two specifications.

Hierarchy B (Malta): no contested.

Hierarchy C (Sanchez-Vega): no contested.

---

## 4. Held-out evaluation protocol

**Identical to Paper 1** — the held-out test set used for the +6.26 nat result is reused.

- Split seed `42`, stratified random 20% per cancer cohort. 5,386 training patients, 1,362 test patients.
- All Phase 1 specs evaluate on the same 1,362-patient test set.
- The baseline-on-train fit is reused from Paper 1 (no refit needed). If Paper 1 chain RDAs are not available on disk (TBD), a baseline-on-train refit is added to Phase 1.

---

## 5. Gibbs run parameters

**Identical to Paper 1.**

- Chains: 4 per spec, seeds 20260515, 20260516, 20260517, 20260518 (one offset per spec to avoid seed collision with Paper 1's 20260509-20260512).
- Iterations: 100,000 per chain.
- Burn-in: 50,000.
- Priors: identical to Paper 1.
- `pi_generation`: `"shared_across_cancers"`.
- Posterior thinning for held-out evaluation: every 400 post-burn iters → 500 posterior draws per model.

---

## 6. Predictive metric

**Identical to Paper 1.** Held-out log-likelihood (LPPD), sum over 1,362 test patients. Per the log-normal AFT spike-and-slab likelihood (right-censored).

---

## 7. Decision rule per specification

Bootstrap 95% CI on held-out LPPD difference per spec vs Paper 1 baseline-on-train, B=1000 patient-level resamples.

| outcome | conditions | per-spec disposition |
|---|---|---|
| **matches** baseline | observed Δ ∈ ±2 nats AND bootstrap 95% CI includes zero | this hierarchy is methodologically inert (does not improve, does not degrade) |
| **outperforms** baseline | observed Δ > +2 nats AND bootstrap 95% CI excludes zero (positive) | this hierarchy generalizes the +6.26 nat result; methodology is hierarchy-portable |
| **underperforms** baseline | observed Δ < −2 nats OR bootstrap 95% CI excludes zero (negative) | this hierarchy degrades performance; methodology is hierarchy-specific |

**Joint disposition across the three hierarchies:**

| joint pattern | overall reading |
|---|---|
| **3/3 outperform** | strong corpus result: methodology is general across external hierarchies of this class |
| **2/3 outperform, 1 inert** | corpus partial-positive: most external hierarchies generalize, one does not |
| **1/3 outperform, 2 inert** | weak corpus positive: only specific hierarchies work |
| **0/3 outperform** | Paper 1 was hierarchy-specific (specific to epithelial-class and germ-layer) |
| **any spec underperforms** | reported as evidence of method-hierarchy mismatch |

---

## 8. Convergence halt rule

**Identical to Paper 1.** R-hat > 1.05 in any of the three Phase 1 specifications halts that specification and reports the convergence failure. The other specifications proceed independently.

---

## 9. Pre-registered constraints

Locked before the Phase 1 Gibbs runs:

1. **No covariate added beyond the four pre-specified.** Same as Paper 1.
2. **No hierarchy tuning based on intermediate results.** The cancer-level assignments above are fixed (saved as `reference/phase1_hierarchy_assignments.csv`).
3. **No decision-rule adjustment.** ±2 nat threshold and CI-excludes-zero requirement are committed before any Phase 1 evaluation.
4. **One held-out split (the Paper 1 split).** Reused for direct comparability to the +6.26 nat result.
5. **All three specifications run regardless of any one's outcome.** No spec is contingent on another.
6. **No re-running with adjusted priors or alternative starting values.** Convergence is checked once per spec; if R-hat > 1.05, halt rather than re-fit.

---

## 10. What is reported regardless of outcome

- Held-out total LPPD for each of: baseline-on-train (from Paper 1 or re-fit), Phase 1 spec A (Thorsson), Phase 1 spec B (Malta), Phase 1 spec C (Sanchez-Vega).
- Per-cancer LPPD breakdown (29 rows) per spec.
- Bootstrap CI on three pairwise differences: each Phase 1 spec vs baseline-on-train.
- Bootstrap CI on three pairwise comparisons among Phase 1 specs.
- Convergence diagnostics for all Phase 1 fits (R-hat, bulk-ESS, tail-ESS, per-chain mean cross-check).
- Wall-clock per chain.
- Per-spec verdict per the decision rule.
- Joint disposition across the three specs.

---

## 11. Provenance

Hierarchies built from external published sources, downloaded from authoritative depositories:

- Thorsson 2018 immune subtypes: Elsevier CDN, `1-s2.0-S1074761318301213-mmc2.xlsx`, MD5 `3f1705d92a3b6a82fdc9eae5450ebd03`, fetched 2026-05-15.
- Malta 2018 stemness scores: Elsevier CDN, `1-s2.0-S0092867418303581-mmc1.xlsx`, MD5 `dc8fe343da78689110176aa3f94008d3`, fetched 2026-05-15.
- Sanchez-Vega 2018 pathway alterations: Washington University Digital Commons mirror (Cell mmc4 equivalent), `viewcontent.cgi?filename=3&article=8671&context=open_access_pubs&type=additional`, MD5 `d4080d1148ee89f6bb041c244e4ac7a5`, fetched 2026-05-15.

Build script: [scripts/validation/build_phase1_hierarchies.py](scripts/validation/build_phase1_hierarchies.py). Output: [reference/phase1_hierarchy_assignments.csv](reference/phase1_hierarchy_assignments.csv).

This pre-registration is locked at commit SHA: **TBD on GitHub push**. The GitHub commit timestamp is the methods-section citation.
