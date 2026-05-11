# Paper 2 — Operator-composition extension results

Generated: 2026-05-10 22:21:59

Anchor: Paper 1 primary 4-cov spec, total LPPD = -1032.43.
Baseline-on-train (diagnostic anchor) total LPPD = -1038.69.
Bootstrap B = 1000, seed 20260510, Bonferroni alpha = 0.002778 (= 0.05 / 18).

## Pre-registered decision rule

Per `PRE_REGISTRATION_PAPER2.md` (committed before compute fired):

- **IMPROVES**: Delta > +2 nats AND CI excludes 0 below AND Bonferroni p < 0.00278
- **NOMINAL-ONLY**: passes nominal 95% CI but fails Bonferroni (NOT counted toward validation)
- **MATCHES**: |Delta| <= 2 AND CI includes 0
- **UNDERPERFORMS**: Delta < -2 OR CI excludes 0 below
- **CONVERGENCE-FAILED**: R-hat > 1.05 anywhere monitored (halt-don't-replace, excluded)

**Validation iff**: >=2 specs IMPROVE AND >=1 of those is RULE 2 or RULE 3.

## Verdict

**FALSIFIES**

- Total specs: 18
- CONVERGENCE-FAILED (excluded per pre-reg): 0
- IMPROVES (Bonferroni p<0.00278): 0
- of those, RULE 2 or RULE 3: 0
- NOMINAL-ONLY: 0
- MATCHES: 0
- UNDERPERFORMS: 12
- INCONCLUSIVE: 6

## Per-spec results table

| spec | cov | rule | pair | R-hat max | conv | bonf | Delta vs P1 | 95% CI | bonf p | total LPPD | verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|
| spec_13_cov1_1_rule1_pairA | 1.1 | R1 | A | 1.01 | Y | N | -7.57 | [-13.55, -1.85] | 0.0160 | -1039.998 | UNDERPERFORMS |
| spec_14_cov1_1_rule1_pairB | 1.1 | R1 | B | 1.01 | Y | N | -4.20 | [-9.14, +1.21] | 0.1359 | -1036.626 | INCONCLUSIVE |
| spec_15_cov1_1_rule2_pairA | 1.1 | R2 | A | 1.01 | Y | N | -8.99 | [-15.25, -3.00] | 0.0040 | -1041.421 | UNDERPERFORMS |
| spec_16_cov1_1_rule2_pairB | 1.1 | R2 | B | 1.01 | Y | N | -7.45 | [-14.33, -1.00] | 0.0260 | -1039.878 | UNDERPERFORMS |
| spec_17_cov1_1_rule3_pairA | 1.1 | R3 | A | 1.01 | Y | N | -6.96 | [-12.52, -1.08] | 0.0180 | -1039.389 | UNDERPERFORMS |
| spec_18_cov1_1_rule3_pairB | 1.1 | R3 | B | 1.00 | Y | N | -3.58 | [-8.48, +1.57] | 0.1918 | -1036.006 | INCONCLUSIVE |
| spec_01_cov8_2_rule1_pairA | 8.2 | R1 | A | 1.02 | Y | N | -5.44 | [-10.95, +0.31] | 0.0659 | -1037.872 | INCONCLUSIVE |
| spec_02_cov8_2_rule1_pairB | 8.2 | R1 | B | 1.01 | Y | N | -6.77 | [-12.39, -1.09] | 0.0260 | -1039.197 | UNDERPERFORMS |
| spec_03_cov8_2_rule2_pairA | 8.2 | R2 | A | 1.01 | Y | N | -7.15 | [-13.22, -1.31] | 0.0160 | -1039.585 | UNDERPERFORMS |
| spec_04_cov8_2_rule2_pairB | 8.2 | R2 | B | 1.01 | Y | N | -6.88 | [-12.87, -1.10] | 0.0200 | -1039.313 | UNDERPERFORMS |
| spec_05_cov8_2_rule3_pairA | 8.2 | R3 | A | 1.01 | Y | N | -6.13 | [-11.93, -0.29] | 0.0400 | -1038.565 | UNDERPERFORMS |
| spec_06_cov8_2_rule3_pairB | 8.2 | R3 | B | 1.01 | Y | N | -7.85 | [-13.89, -1.88] | 0.0120 | -1040.281 | UNDERPERFORMS |
| spec_07_cov10_1_rule1_pairA | 10.1 | R1 | A | 1.01 | Y | N | -6.56 | [-12.73, -0.20] | 0.0460 | -1038.991 | UNDERPERFORMS |
| spec_08_cov10_1_rule1_pairB | 10.1 | R1 | B | 1.01 | Y | N | -5.45 | [-11.23, +0.34] | 0.0619 | -1037.878 | INCONCLUSIVE |
| spec_09_cov10_1_rule2_pairA | 10.1 | R2 | A | 1.00 | Y | N | -6.82 | [-13.32, -0.21] | 0.0400 | -1039.254 | UNDERPERFORMS |
| spec_10_cov10_1_rule2_pairB | 10.1 | R2 | B | 1.02 | Y | N | -7.48 | [-13.90, -1.27] | 0.0160 | -1039.91 | UNDERPERFORMS |
| spec_11_cov10_1_rule3_pairA | 10.1 | R3 | A | 1.01 | Y | N | -6.01 | [-11.94, +0.17] | 0.0559 | -1038.441 | INCONCLUSIVE |
| spec_12_cov10_1_rule3_pairB | 10.1 | R3 | B | 1.01 | Y | N | -5.37 | [-10.92, +0.62] | 0.0819 | -1037.805 | INCONCLUSIVE |

## Anchors used

- Paper 1 primary 4-cov LPPD: -1032.43 (per-patient at run_validation_2026-05-09/loglik_per_patient_primary.csv)
- baseline-on-train LPPD: -1038.69 (per-patient at run_validation_2026-05-09/loglik_per_patient_baseline_train.csv)

## Caveats

- Bootstrap p-values use mid-p adjustment ((B_tail + 1) / (B + 1)). Smallest achievable p with B=1000 is ~0.002.
- Anchor for all 18 specs is the Paper 1 4-covariate joint-fit total. Per-cov anchors would be different but no individual cov was fit alone in Paper 1.
- Per the pre-reg, no spec is replaced or re-fit if convergence fails.
