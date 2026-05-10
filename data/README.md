# `data/` — Pointers, not files

This directory does not contain the raw TCGA omics, the pre-computed BIDIFAC+
inputs, or the chain RDA files. Those are large and live elsewhere.

---

## What's needed and where to get it

### 1. Pre-computed BIDIFAC+ inputs (1.5 MB)

`XYC_V2_WithAge_StandardizedPredictors.rda` ships with the upstream
replication repository:

> https://github.com/sarahsamorodnitsky/HierarchicalSS_PanCanPanOmics

Clone that repository to access the file directly. This is the only data
file required to reproduce the replication baseline (§3 of the top-level
README) and the validation Gibbs runs (§5). The pre-computed inputs contain
the post-BIDIFAC+, post-clinical-matching `Covariates` / `Survival` /
`Censored` lists for the 29-type pan-cancer cohort.

```bash
git clone https://github.com/sarahsamorodnitsky/HierarchicalSS_PanCanPanOmics.git
# the .rda files are at the repo root
ls HierarchicalSS_PanCanPanOmics/XYC_V2_WithAge_StandardizedPredictors.rda
```

### 2. Raw TCGA pan-cancer omics (4.2 GB) — only if rerunning BIDIFAC+

The four GDC PanCan-CellOfOrigin files, downloaded from
<https://gdc.cancer.gov/about-data/publications/PanCan-CellOfOrigin>:

| filename | UUID | size |
|---|---|---|
| `EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv` | `9a4679c3-855d-4055-8be9-3577ce10f66e` | 1.8 GB |
| `TCGA-RPPA-pancan-clean.txt` | `fcbb373e-28d4-4818-92f3-601ede3da5e1` | 19 MB |
| `pancanMiRs_EBadjOnProtocolPlatformWithoutRepsWithUnCorrectMiRs_08_04_16.csv` | `1c6174d9-8ffb-466e-b5ee-07b204c15cf8` | 65 MB |
| `jhu-usc.edu_PANCAN_merged_HumanMethylation27_HumanMethylation450.betaValue_whitelisted.tsv` | `d82e2c44-89eb-43d9-b6d3-712732bf6a53` | 2.3 GB |

URL pattern: `https://api.gdc.cancer.gov/data/<UUID>` — public, no auth required.

These files are **only needed if BIDIFAC+ is rerun from raw**. The
replication baseline and the validation in this repository use the
pre-computed `XYC_V2_*.rda` directly.

If rerunning BIDIFAC+:

```bash
mkdir -p data/PanTCGA
cd data/PanTCGA
curl -L -o PanCanExp.tsv \
  https://api.gdc.cancer.gov/data/9a4679c3-855d-4055-8be9-3577ce10f66e
curl -L -o TCGA-RPPA-pancan-clean.txt \
  https://api.gdc.cancer.gov/data/fcbb373e-28d4-4818-92f3-601ede3da5e1
curl -L -o PanCan_miRNA.csv \
  https://api.gdc.cancer.gov/data/1c6174d9-8ffb-466e-b5ee-07b204c15cf8
curl -L -o PANCAN_meth_merged.tsv \
  https://api.gdc.cancer.gov/data/d82e2c44-89eb-43d9-b6d3-712732bf6a53
```

The BIDIFAC+ pipeline is in the upstream replication repo at
`TCGA_BIDIFAC_plus.R`.

### 3. Chain RDA files (~7 GB) — Zenodo

Posterior chain RDAs for the three model fits are not in this repository.
They are too large (12 chains × ~580 MB = ~7 GB total) and are deposited at
Zenodo with a citable DOI:

> Zenodo deposit DOI: **TBD** — uploaded after public release of this repo.
> Will contain:
> - `gibbs_baseline_train_chain_[1-4].rda` (4 × ~588 MB) — replication on 80% training set
> - `gibbs_primary_chain_[1-4].rda` (4 × ~582 MB) — primary structural-prior spec
> - `gibbs_secondary_chain_[1-4].rda` (4 × ~580 MB) — secondary structural-prior spec
> - The original full-data baseline (`runs/run_paper_a_gibbs_2026-05-09/`),
>   4 × 588 MB, also included for replication comparison.

Without these RDAs, the convergence diagnostics, held-out log-likelihood,
bootstrap, and per-cancer breakdown cannot be recomputed from the saved
posteriors directly — but the full pipeline can be re-run from scratch
using the scripts in `scripts/` (~2 hours wall on a modern workstation).

### 4. TCGA-CDR clinical resource (1.5 MB) — only for downstream figures

`TCGA-CDR.csv` (Liu et al. 2018, Cell) is referenced in the upstream
`ScoresVsSubtypes.R` for figures comparing posterior scores to published
subtypes. Available from the supplementary materials of:

> Liu et al. *An Integrated TCGA Pan-Cancer Clinical Data Resource to Drive
> High-Quality Survival Outcome Analytics.* Cell. 2018.
> PMID: 29625055 — PMC6066282.

Not required for any analysis in this repository.

---

## Reproducibility chain

To reproduce all results in this repository starting from publicly available
data:

1. **Replication baseline** (no Zenodo needed): clone the upstream
   `HierarchicalSS_PanCanPanOmics` repo, install R packages (see §9 of the
   top-level README), and run:
   ```bash
   Rscript scripts/replication/run_gibbs_chain.R 1   # ... and 2, 3, 4
   Rscript scripts/replication/diagnostics.R
   Rscript scripts/replication/per_cancer_summary.R
   ```
   ~36 min wall total with 4 parallel chains.

2. **Screening + falsification** (only needs the per-cancer table from step 1
   and the hierarchy reference at `reference/cancer_type_hierarchies_2026-05-09.md`):
   ```bash
   python scripts/screening/hierarchies_long.py
   python scripts/screening/run_screening.py
   python scripts/screening/falsify_01_perm_null.py
   python scripts/screening/falsify_02_block_presence.py
   python scripts/screening/falsify_03_granularity.py
   python scripts/screening/falsify_04_cohort.py
   python scripts/screening/falsify_05_high_res_perm.py
   python scripts/screening/falsify_06_independence.py
   python scripts/screening/synthesize_summary.py
   ```
   ~30 min wall.

3. **Validation** (needs only the XYC pre-computed inputs):
   ```bash
   Rscript scripts/validation/01_build_split_and_targets.R
   for spec in baseline_train primary secondary; do
     for chain in 1 2 3 4; do
       Rscript scripts/validation/run_chain.R $spec $chain &
     done
     wait
   done
   Rscript scripts/validation/05_convergence.R         # halt if R-hat > 1.05
   Rscript scripts/validation/02_compute_held_out_loglik.R
   Rscript scripts/validation/04_per_cancer_breakdown.R
   Rscript scripts/validation/03_bootstrap_ci.R
   ```
   ~2 hours wall total with 4 parallel chains per spec, three specs run
   sequentially.

The pre-computed scripts assume the working directory is the repository
root. Paths in the scripts may need to be adjusted to point to your local
location of `XYC_V2_WithAge_StandardizedPredictors.rda` (search for
`HierarchicalSS_PanCanPanOmics` in each script).
