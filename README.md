# Pan-Cancer Hierarchical Spike-and-Slab — Replication & Structural-Prior Validation

A pre-registered held-out validation in which structural priors derived from
external biological hierarchies replace the per-cancer fitted PIPs for four
covariates of the Lock-lab pan-cancer survival spike-and-slab model. The
pre-registered decision rule fires **OUTPERFORMS** for the primary
specification.

---

## TL;DR

| model | held-out total LPPD | Δ vs baseline | bootstrap 95% CI | verdict |
|---|---|---|---|---|
| baseline (Samorodnitsky/Hoadley/Lock 2022, fit on 80% training) | −1038.69 | — | — | — |
| **primary** structural-prior | **−1032.43** | **+6.26 nats** | **[+0.34, +12.29]** | **OUTPERFORMS** |
| secondary structural-prior | −1034.36 | +4.33 nats | [−1.86, +10.97] | inconclusive |

Held-out test set: 1,362 patients across 29 cancer types, stratified random
20% per cohort, seed 42. Four chains × 100,000 Gibbs iters per model,
50k burn, 500 thinned posterior samples for LPPD evaluation, B=1000
patient-level bootstrap.

Decision rule (locked before any validation Gibbs run, see
`PRE_REGISTRATION.md`): a model "outperforms" the baseline if held-out LPPD
exceeds the baseline by > +2 nats *and* the bootstrap 95% CI on the
difference excludes zero. Primary specification meets both criteria.

The structural priors replace 116 free per-cancer beta parameters with 31
group-shared betas across the four target covariates (a 73% reduction).
Held-out predictive performance improves nonetheless — consistent with
productive regularization, not flexibility loss.

---

## Provenance

The replicated method:

> Samorodnitsky, S., Hoadley, K. A., & Lock, E. F. (2022).
> *A Hierarchical Spike-and-Slab Model for Pan-Cancer Survival Using Pan-Omic Data.*
> BMC Bioinformatics 23(1):235.
> Source: <https://github.com/sarahsamorodnitsky/HierarchicalSS_PanCanPanOmics>

The BIDIFAC+ inputs:

> Park, J., Lock, E. F., & Hoadley, K. A. (2022).
> *Bidimensional Linked Matrix Factorization for Pan-Omics Pan-Cancer Analysis.*
> Annals of Applied Statistics 16(1):484-501.

The 29 TCGA cancer types and pre-computed `XYC_V2_WithAge_StandardizedPredictors.rda`
inputs come directly from the upstream repository. No re-fit of BIDIFAC+ was
performed.

---

## What's in this repo

```
.
├── README.md                                 ← this file
├── README.pdf                                ← printable version
├── LICENSE                                   ← MIT (this work) + upstream attribution
├── PRE_REGISTRATION.md                       ← locked decision rules
├── REPLICATION_RESULTS.md                    ← baseline replication numbers
├── SCREENING_SUMMARY.md                      ← screening pass output
├── FALSIFICATION_REPORT.md                   ← 6-check adversarial audit of screening
├── VALIDATION_RESULTS.md                     ← pre-registered held-out validation
├── reference/
│   └── cancer_type_hierarchies_2026-05-09.md ← 29 cancers × 5 hierarchies × granularities
├── scripts/
│   ├── replication/                          ← R Gibbs runner + diagnostics
│   ├── screening/                            ← Python ANOVA + falsification cascade
│   └── validation/                           ← R modified-Gibbs + held-out evaluation
├── results/
│   ├── replication/                          ← per-cancer × per-cov posterior table + diagnostics
│   ├── screening/                            ← screening + 6 falsification CSVs
│   └── validation/                           ← LPPD, bootstrap, per-cancer, convergence
└── data/
    └── README.md                             ← pointers to XYC, TCGA, Zenodo (chain RDAs)
```

The chain RDA files (12 chains × ~580 MB = ~7 GB total) are **not in this
repository**. They will be deposited at Zenodo with a citable DOI; the
deposit URL goes here:

> **Zenodo DOI: TBD** (will be filled in after the public release of this repo)

Without the chain RDAs, downstream analyses (convergence, held-out LPPD,
bootstrap) cannot be recomputed from saved posteriors directly. The full
pipeline can still be re-run from scratch using the scripts here in
~2 hours wall on a modern workstation. See `data/README.md`.

---

## What was done, briefly

### 1. Replication baseline
Four-chain replication of `HierarchicalLogNormalSpikeSlab` (the upstream
sampler, unmodified) on the full pre-computed XYC inputs. R-hat max 1.0152
across 137 monitored params; per-chain σ² posteriors agree to 3 decimals.
Cancer-specific intercepts span THCA 10.28 (≈80 yr at standardized
predictors=0) to MESO 6.29 (≈1.5 yr) — clinically expected order. Pan-cancer
age effect β̃₀.₅ = −0.316 (PIP 0.91). See `REPLICATION_RESULTS.md`.

### 2. Structural-hierarchy screening
ANOVA of per-cancer PIPs against five candidate hierarchies (tissue origin,
epithelial classification, embryonic germ layer, sex composition, Hoadley
2018 supercluster/iCluster) at multiple granularities each. 636 tests total
across 68 BIDIFAC+ covariates. 16 nominal hits at η² > 0.30 with p < 0.05.
See `SCREENING_SUMMARY.md`.

### 3. Falsification cascade (six adversarial checks)
1. Permutation null on η² (B=1000): observed at chance rate
2. High-resolution perm on survivors (B=10000): 0/16 pass Bonferroni
3. BIDIFAC+ block-presence confound: Hoadley signal ρ=+0.60 with block presence (circular)
4. Granularity-nesting null on Spearman ρ: tissue origin marginal (p=0.042), all others fail
5. Cohort-size partialing: most survivors lose < 0.06 η², one loses 0.15
6. Independence: 16 survivors collapse to 5 effective dimensions; 32/120 pairs Jaccard > 0.5

No single (cov, hierarchy, granularity) cell is publishable as standalone
structural evidence at conventional rigor. See `FALSIFICATION_REPORT.md`.

### 4. Pre-registered validation
The screening result motivated a different question: do the four
largest-effect-size candidates *jointly* improve held-out predictive
performance when used as structural priors, even though no single candidate
is individually Bonferroni-significant? Pre-registered in
`PRE_REGISTRATION.md`, ran the modified Gibbs sampler at
`scripts/validation/sampler_structural_prior.R`, evaluated held-out LPPD on
the 20% test set with B=1000 bootstrap CI. Result: primary specification
outperforms baseline (Δ=+6.26 nats, CI excludes zero). See
`VALIDATION_RESULTS.md`.

---

## How to reproduce

Two paths:

### Fast path — verify reported numbers from CSV outputs
The small CSV outputs in `results/` are sufficient to regenerate every
table in the four reports without rerunning any Gibbs sampler. Inspect
`results/validation/loglik_summary.csv`, `bootstrap_ci_summary.csv`,
`per_cancer_loglik_wide.csv`, etc. directly.

### Full path — re-fit everything from raw data
~2 hours wall on a modern workstation, four parallel chains per Gibbs run,
three Gibbs runs sequentially. Step-by-step instructions in
`data/README.md`. Briefly:

```bash
# 1. Pre-computed XYC inputs (clone upstream repo)
git clone https://github.com/sarahsamorodnitsky/HierarchicalSS_PanCanPanOmics.git

# 2. Replication baseline (4 chains × 100k iters; ~36 min wall)
for c in 1 2 3 4; do
  Rscript scripts/replication/run_gibbs_chain.R $c &
done
wait
Rscript scripts/replication/diagnostics.R
Rscript scripts/replication/per_cancer_summary.R

# 3. Screening + 6 falsification checks (~30 min)
python scripts/screening/hierarchies_long.py
python scripts/screening/run_screening.py
python scripts/screening/falsify_0{1,2,3,4,5,6}_*.py
python scripts/screening/synthesize_summary.py

# 4. Pre-registered validation (3 Gibbs runs sequentially × 4 parallel chains; ~2 h)
Rscript scripts/validation/01_build_split_and_targets.R
for spec in baseline_train primary secondary; do
  for c in 1 2 3 4; do
    Rscript scripts/validation/run_chain.R $spec $c &
  done
  wait
done
Rscript scripts/validation/05_convergence.R       # halt if R-hat > 1.05
Rscript scripts/validation/02_compute_held_out_loglik.R
Rscript scripts/validation/04_per_cancer_breakdown.R
Rscript scripts/validation/03_bootstrap_ci.R
```

Software: R 4.6.0 with `posterior` 1.7.0 + the upstream BSFP dependencies;
Python 3.12 with pandas 2.3.3, scipy 1.17.1, numpy 2.4.4. The
`infinitefactor` package was archived from CRAN on 2025-12-25 and the
upstream code uses it indirectly through BSFP; a pure-R port of the
relevant functions is documented in the validation pipeline if needed.

Random seeds are fixed throughout: chains 1–4 use seeds 20260509–20260512;
the held-out split uses seed 42; bootstrap uses 20260510; permutation
nulls use 20260509 and 20260511.

---

## Caveats

- **Validation tests four pre-specified covariates only**, not all 68. The
  result does not generalize to the full covariate set.
- **Single 80/20 split**, not k-fold cross-validation. Tighter CIs and
  per-cancer subgroup tests would require k-fold.
- **Tissue origin** is not used in any validation specification. It escapes
  the BIDIFAC+ block-presence confound but only marginally beats the
  granularity-nesting null. A tissue-origin spec would need a separate
  pre-registration and held-out test.
- **KIRP and STAD** lose under both validation specifications. The gain is
  heterogeneous; understanding why these cancers worsen is a candidate
  follow-up.
- **n=29** is small. The screening was demonstrably underpowered for
  319 simultaneous tests at Bonferroni rigor. The validation succeeded by
  testing joint predictive contribution rather than per-cell significance.

---

## Citing

This repository (provisional):

> Humphrey, N. (2026). *Pan-Cancer Hierarchical Spike-and-Slab — Replication
> and Structural-Prior Validation.* GitHub: [URL TBD]. Zenodo: [DOI TBD].

The replicated method:

> Samorodnitsky, S., Hoadley, K. A., & Lock, E. F. (2022).
> *A Hierarchical Spike-and-Slab Model for Pan-Cancer Survival Using Pan-Omic Data.*
> BMC Bioinformatics 23(1):235.

The BIDIFAC+ inputs:

> Park, J., Lock, E. F., & Hoadley, K. A. (2022).
> *Bidimensional Linked Matrix Factorization for Pan-Omics Pan-Cancer Analysis.*
> Annals of Applied Statistics 16(1):484-501.

---

## License

MIT. The replication-baseline scripts in `scripts/replication/` are
patched verbatim from the upstream repository. Patches to upstream files
are documented inline in source. See `LICENSE`.
