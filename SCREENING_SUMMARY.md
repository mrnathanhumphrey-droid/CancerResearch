# Structural-Hierarchy Screening — SUMMARY

**Run date**: 2026-05-09

**Source PIPs**: `runs/run_paper_a_gibbs_2026-05-09/per_cancer_beta_pip.csv` (551 rows, 29 cancers × variable covariate sets, 68 unique covariate IDs).

**Hierarchies tested**: 5 (tissue_origin, epithelial_class, germ_layer, sex_composition, hoadley) × 10 granularities. Granularity 1c (29 singleton bins) excluded as degenerate.

**Tests run**: 636 (cov × hierarchy × granularity), 233 non-degenerate after the n_singleton_groups > n_groups/2 filter.

**Statistical metric**: η² = SS_between / SS_total per one-way ANOVA. F and p reported but η² is the primary screening metric (n=29 with mostly small group counts).


---

## INTERPRETATION FRAMING (READ FIRST)

This is **screening only**, not validation. Triples flagged here are **candidates** for validation; validation requires a reduced-parameter Gibbs replacement using structural priors at the chosen (hierarchy, granularity), held-out predictive performance vs the baseline replication. That is a separate downstream phase.

**Do not** read these results as 'the closed-form decomposition is valid.' They identify the (hierarchy, granularity, covariate) cells with sufficient between-group variance to be **worth** validating.

---

## 1. Top-20 (covariate × hierarchy × granularity) triples by η²

All valid rows (degenerate ones flagged with † in the deg column). Per the brief, degenerate rows are not dropped — they still contribute to between-group variance.


### Non-Hoadley (external biological hierarchies)

| cov_id | hierarchy | granularity | n_groups | n_cancers | n_single | η² | F | p | deg |
|---|---|---|---|---|---|---|---|---|---|
| 26.1 | tissue_origin | 1b_11bin | 6 | 6 | 6 | 1.000 | — | — | † |
| 35.1 | epithelial_class | 2c_10bin | 4 | 4 | 4 | 1.000 | — | — | † |
| 37.1 | tissue_origin | 1b_11bin | 4 | 4 | 4 | 1.000 | — | — | † |
| 43.1 | tissue_origin | 1a_5bin | 3 | 3 | 3 | 1.000 | — | — | † |
| 43.1 | tissue_origin | 1b_11bin | 3 | 3 | 3 | 1.000 | — | — | † |
| 43.2 | tissue_origin | 1a_5bin | 3 | 3 | 3 | 1.000 | — | — | † |
| 43.2 | tissue_origin | 1b_11bin | 3 | 3 | 3 | 1.000 | — | — | † |
| 44.1 | tissue_origin | 1a_5bin | 5 | 5 | 5 | 1.000 | — | — | † |
| 44.1 | tissue_origin | 1b_11bin | 5 | 5 | 5 | 1.000 | — | — | † |
| 48.1 | epithelial_class | 2b_6bin | 3 | 3 | 3 | 1.000 | — | — | † |
| 48.1 | epithelial_class | 2c_10bin | 3 | 3 | 3 | 1.000 | — | — | † |
| 48.1 | germ_layer | 3b_mid | 3 | 3 | 3 | 1.000 | — | — | † |
| 31.1 | tissue_origin | 1b_11bin | 7 | 8 | 6 | 1.000 | 2830.28 | 0.0144 | † |
| 36.1 | germ_layer | 3b_mid | 4 | 5 | 3 | 0.998 | 140.69 | 0.0619 | † |
| 31.1 | germ_layer | 3b_mid | 3 | 5 | 2 | 0.994 | 164.82 | 0.00603 | † |
| 26.1 | germ_layer | 3a_3bin | 2 | 3 | 1 | 0.987 | 77.79 | 0.0719 |  |
| 26.1 | germ_layer | 3b_mid | 2 | 3 | 1 | 0.987 | 77.79 | 0.0719 |  |
| 50.1 | epithelial_class | 2c_10bin | 5 | 6 | 4 | 0.985 | 16.29 | 0.184 | † |
| 26.1 | epithelial_class | 2b_6bin | 3 | 6 | 2 | 0.957 | 33.11 | 0.00902 | † |
| 26.1 | epithelial_class | 2c_10bin | 3 | 6 | 2 | 0.957 | 33.11 | 0.00902 | † |


### Hoadley (circularity-flagged ⚠)

| cov_id | hierarchy | granularity | n_groups | n_cancers | n_single | η² | F | p | deg |
|---|---|---|---|---|---|---|---|---|---|
| 3.1 | hoadley | 5c_iclusters | 9 | 9 | 9 | 1.000 | — | — | † |
| 3.2 | hoadley | 5c_iclusters | 9 | 9 | 9 | 1.000 | — | — | † |
| 15.1 | hoadley | 5c_iclusters | 5 | 5 | 5 | 1.000 | — | — | † |
| 31.1 | hoadley | 5c_iclusters | 5 | 5 | 5 | 1.000 | — | — | † |
| 34.1 | hoadley | 5c_iclusters | 8 | 8 | 8 | 1.000 | — | — | † |
| 35.1 | hoadley | 5c_iclusters | 3 | 3 | 3 | 1.000 | — | — | † |
| 36.1 | hoadley | 5c_iclusters | 4 | 4 | 4 | 1.000 | — | — | † |
| 37.1 | hoadley | 5c_iclusters | 3 | 3 | 3 | 1.000 | — | — | † |
| 39.1 | hoadley | 5c_iclusters | 5 | 5 | 5 | 1.000 | — | — | † |
| 46.1 | hoadley | 5c_iclusters | 3 | 3 | 3 | 1.000 | — | — | † |
| 49.1 | hoadley | 5c_iclusters | 3 | 3 | 3 | 1.000 | — | — | † |
| 50.1 | hoadley | 5c_iclusters | 4 | 4 | 4 | 1.000 | — | — | † |
| 43.2 | hoadley | 5a_4plus_unassigned | 2 | 3 | 1 | 1.000 | 108080.11 | 0.00194 |  |
| 2.1 | hoadley | 5c_iclusters | 14 | 16 | 12 | 1.000 | 498.07 | 0.00201 | † |
| 4.1 | hoadley | 5c_iclusters | 10 | 12 | 8 | 0.998 | 107.29 | 0.00927 | † |
| 32.1 | hoadley | 5a_4plus_unassigned | 2 | 3 | 1 | 0.989 | 88.64 | 0.0674 |  |
| 1.1 | hoadley | 5c_iclusters | 14 | 16 | 12 | 0.987 | 11.86 | 0.0804 | † |
| 5.1 | hoadley | 5c_iclusters | 14 | 16 | 12 | 0.984 | 9.31 | 0.101 | † |
| 8.2 | hoadley | 5c_iclusters | 14 | 16 | 12 | 0.983 | 8.98 | 0.104 | † |
| 10.1 | hoadley | 5c_iclusters | 14 | 16 | 12 | 0.974 | 5.81 | 0.156 | † |


## 2. η² distribution: Hoadley vs non-Hoadley

Reported on (a) all valid rows and (b) non-degenerate rows.


### ALL valid rows

| stat | non-Hoadley | Hoadley | gap |
|---|---|---|---|
| 0th pct | 0.000 | 0.032 | +0.032 |
| 10th pct | 0.028 | 0.118 | +0.090 |
| 25th pct | 0.096 | 0.221 | +0.125 |
| 50th pct | 0.242 | 0.621 | +0.379 |
| 75th pct | 0.482 | 0.991 | +0.509 |
| 90th pct | 0.903 | 1.000 | +0.097 |
| 100th pct | 1.000 | 1.000 | +0.000 |
| mean | 0.338 | 0.606 | +0.268 |
| count | 259 | 60 | — |

**Median η² gap (Hoadley − non-Hoadley) = +0.379** (⚠ FLAG per user's >0.15 suspicion threshold)


### Non-degenerate only

| stat | non-Hoadley | Hoadley | gap |
|---|---|---|---|
| 0th pct | 0.000 | 0.032 | +0.032 |
| 10th pct | 0.021 | 0.114 | +0.093 |
| 25th pct | 0.069 | 0.127 | +0.058 |
| 50th pct | 0.199 | 0.280 | +0.080 |
| 75th pct | 0.380 | 0.533 | +0.154 |
| 90th pct | 0.681 | 0.831 | +0.150 |
| 100th pct | 0.987 | 1.000 | +0.013 |
| mean | 0.269 | 0.369 | +0.100 |
| count | 198 | 35 | — |

**Median η² gap (Hoadley − non-Hoadley) = +0.080** (OK per user's >0.15 suspicion threshold)


## 3. Broadly-informative hierarchies

Definition: at least one granularity producing η² > 0.30 on at least 5 covariates (non-degenerate rows only).

| hierarchy | granularity | n_covs with η²>0.30 | broadly_informative |
|---|---|---|---|
| epithelial_class | 2c_10bin | 17 | ✓ |
| tissue_origin | 1b_11bin | 16 | ✓ |
| hoadley | 5a_4plus_unassigned | 14 | ✓ |
| epithelial_class | 2b_6bin | 8 | ✓ |
| germ_layer | 3b_mid | 6 | ✓ |
| sex_composition | 4_4bin | 6 | ✓ |
| tissue_origin | 1a_5bin | 6 | ✓ |
| epithelial_class | 2a_3bin | 4 |  |
| germ_layer | 3a_3bin | 3 |  |
| hoadley | 5c_iclusters | 1 |  |

**Broadly-informative hierarchies**: epithelial_class, germ_layer, hoadley, sex_composition, tissue_origin


## 4. Structurally-explainable covariates

Definition: at least one (hierarchy, granularity) producing η² > 0.30 (non-degenerate, non-Hoadley).

**30 covariates** have at least one non-Hoadley (hierarchy, granularity) hit at η² > 0.30.

| cov_id | n_hits | best (hier, gran) | best η² |
|---|---|---|---|
| 11.1 | 5 | tissue_origin / 1b_11bin | 0.704 |
| 20.1 | 4 | germ_layer / 3b_mid | 0.813 |
| 26.1 | 4 | germ_layer / 3a_3bin | 0.987 |
| 8.2 | 4 | germ_layer / 3b_mid | 0.884 |
| 1.1 | 3 | epithelial_class / 2c_10bin | 0.643 |
| 22.1 | 3 | tissue_origin / 1b_11bin | 0.845 |
| 8.3 | 3 | sex_composition / 4_4bin | 0.482 |
| 48.1 | 3 | epithelial_class / 2a_3bin | 0.575 |
| 31.1 | 3 | epithelial_class / 2c_10bin | 0.833 |
| 2.1 | 2 | germ_layer / 3b_mid | 0.524 |
| 24.1 | 2 | epithelial_class / 2c_10bin | 0.443 |
| 15.1 | 2 | tissue_origin / 1b_11bin | 0.407 |
| 13.1 | 2 | epithelial_class / 2b_6bin | 0.336 |
| 0.5 | 2 | epithelial_class / 2c_10bin | 0.390 |
| 8.1 | 2 | tissue_origin / 1b_11bin | 0.371 |
| 44.1 | 2 | epithelial_class / 2b_6bin | 0.915 |
| 23.1 | 2 | tissue_origin / 1b_11bin | 0.324 |
| 25.1 | 2 | epithelial_class / 2b_6bin | 0.857 |
| 42.1 | 2 | tissue_origin / 1b_11bin | 0.540 |
| 3.2 | 2 | sex_composition / 4_4bin | 0.631 |
| 39.1 | 2 | tissue_origin / 1a_5bin | 0.482 |
| 36.1 | 2 | epithelial_class / 2b_6bin | 0.671 |
| 17.1 | 1 | sex_composition / 4_4bin | 0.371 |
| 10.1 | 1 | epithelial_class / 2c_10bin | 0.900 |
| 3.1 | 1 | epithelial_class / 2c_10bin | 0.440 |
| 28.1 | 1 | tissue_origin / 1b_11bin | 0.390 |
| 34.1 | 1 | tissue_origin / 1b_11bin | 0.442 |
| 43.1 | 1 | sex_composition / 4_4bin | 0.581 |
| 46.1 | 1 | sex_composition / 4_4bin | 0.311 |
| 5.1 | 1 | epithelial_class / 2c_10bin | 0.552 |


## 5. Hoadley vs best non-Hoadley η² per covariate

For each covariate where BOTH Hoadley and non-Hoadley produce η² > 0.20 (using ALL valid rows). Median gap > 0 = Hoadley wins; Hoadley dominance suggests TCGA-circularity rather than external biology.

| cov_id | best non-Hoadley | nh η² | Hoadley η² | gap |
|---|---|---|---|---|
| 4.1 | germ_layer/3b_mid | 0.278 | 0.998 | +0.720 |
| 49.1 | germ_layer/3b_mid | 0.369 | 1.000 | +0.631 |
| 15.1 | tissue_origin/1b_11bin | 0.407 | 1.000 | +0.593 |
| 3.1 | epithelial_class/2c_10bin | 0.440 | 1.000 | +0.560 |
| 34.1 | tissue_origin/1b_11bin | 0.442 | 1.000 | +0.558 |
| 0.5 | epithelial_class/2c_10bin | 0.390 | 0.909 | +0.519 |
| 2.1 | germ_layer/3b_mid | 0.524 | 1.000 | +0.476 |
| 8.1 | tissue_origin/1b_11bin | 0.371 | 0.824 | +0.453 |
| 5.1 | epithelial_class/2c_10bin | 0.552 | 0.984 | +0.432 |
| 3.2 | sex_composition/4_4bin | 0.631 | 1.000 | +0.369 |
| 8.3 | sex_composition/4_4bin | 0.482 | 0.849 | +0.367 |
| 1.1 | epithelial_class/2c_10bin | 0.643 | 0.987 | +0.344 |
| 39.1 | germ_layer/3b_mid | 0.735 | 1.000 | +0.265 |
| 20.1 | germ_layer/3b_mid | 0.813 | 0.973 | +0.159 |
| 8.2 | germ_layer/3b_mid | 0.884 | 0.983 | +0.099 |
| 10.1 | epithelial_class/2c_10bin | 0.900 | 0.974 | +0.075 |
| 46.1 | tissue_origin/1a_5bin | 0.952 | 1.000 | +0.048 |
| 24.1 | epithelial_class/2c_10bin | 0.443 | 0.462 | +0.019 |
| 50.1 | epithelial_class/2c_10bin | 0.985 | 1.000 | +0.015 |
| 36.1 | germ_layer/3b_mid | 0.998 | 1.000 | +0.002 |
| 31.1 | tissue_origin/1b_11bin | 1.000 | 1.000 | +0.000 |
| 11.1 | epithelial_class/2c_10bin | 0.813 | 0.813 | +0.000 |
| 35.1 | epithelial_class/2c_10bin | 1.000 | 1.000 | +0.000 |
| 37.1 | tissue_origin/1b_11bin | 1.000 | 1.000 | +0.000 |
| 43.2 | tissue_origin/1a_5bin | 1.000 | 1.000 | -0.000 |
| 26.1 | tissue_origin/1b_11bin | 1.000 | 0.974 | -0.026 |
| 42.1 | tissue_origin/1b_11bin | 0.540 | 0.501 | -0.039 |
| 22.1 | tissue_origin/1b_11bin | 0.845 | 0.686 | -0.159 |
| 13.1 | germ_layer/3b_mid | 0.702 | 0.299 | -0.403 |
| 48.1 | epithelial_class/2b_6bin | 1.000 | 0.565 | -0.435 |

**Median gap across 30 covariates: +0.087** (positive = Hoadley dominates; negative = external hierarchy outperforms).


## 6. Granularity robustness

Spearman ρ between covariate η² rankings at different granularities of the same hierarchy. High ρ = same covariates land on top regardless of bin choice (real structure). Low ρ = granularity is load-bearing (k=4 vs k=6 patterns).

| hierarchy | (g1, g2) | Spearman ρ | n_covs_compared |
|---|---|---|---|
| epithelial_class | (2a_3bin, 2b_6bin) | +0.253 | 18 |
| epithelial_class | (2a_3bin, 2c_10bin) | -0.238 | 16 |
| epithelial_class | (2b_6bin, 2c_10bin) | +0.598 | 22 |
| germ_layer | (3a_3bin, 3b_mid) | +0.253 | 16 |
| tissue_origin | (1a_5bin, 1b_11bin) | +0.716 | 19 |


---
## Files

- `screening_results.csv` — long-format results (every (cov, hier, gran) test, including degenerate rows)
- `hierarchies_long.csv` — cancer × hierarchy × granularity × group (with contested flag)
- `data_handling_log.md` — design choices: contested-cell drops, canonical group placements for split histologies, Hoadley iCluster primary-assignment choices
- `run_screening.py` — ANOVA pipeline; `hierarchies_long.py` — hierarchy CSV builder; `synthesize_summary.py` — this report.
