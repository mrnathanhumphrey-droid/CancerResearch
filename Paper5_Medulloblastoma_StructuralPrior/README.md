# Paper 5 — Within-Cancer Structural-Prior Methodology Transfer to Medulloblastoma

**Status (2026-05-15):** pre-registration locked at [PRE_REGISTRATION.md](PRE_REGISTRATION.md). No compute fired. Data not yet pulled. Scripts not yet implemented.

---

## What this paper tests

The Lock 2022 hierarchical-spike-and-slab structural-prior methodology, validated on Paper 1 against a pan-cancer (29-cancer) cohort with external epithelial-classification and germ-layer hierarchies, is here applied to a **single-cancer** corpus with the structural-prior hierarchy being the **within-cancer molecular subtype taxonomy**.

Substrate: medulloblastoma. Cohort: Cavalli et al. 2017, *Cancer Cell*, n = 763 primary tumors with matched DNA methylation (Illumina 450k) and gene expression (Affymetrix Hu Gene 1.1 ST) — both fully public via GEO accessions GSE85212 and GSE85217. Hierarchy: Cavalli's 12-subtype molecular classification (primary granularity), with the 4-subgroup WHO consensus as a coarser sensitivity granularity.

This is a **within-cancer test** complementary to Paper 4 Phase 1's **between-hierarchy test** on the same 29-cancer pan-cancer substrate. Paper 4 found the methodology is hierarchy-specific (0/3 outperform under matched specs). Paper 5 asks whether the methodology transfers to a substrate where the hierarchy is naturally within-cancer rather than imposed externally.

---

## Why medulloblastoma + Cavalli 2017

1. **Three nested granularity levels are native to the data.** Lock 2022 needed to scan 5 external hierarchies × 11 granularities to find the epi-class / germ-layer combination. Medulloblastoma comes pre-annotated with a 4-subgroup / 12-subtype / 8-subtype-within-G3-G4 (Sharma 2019) hierarchy — the granularity sweep is built in, not imposed.
2. **WNT singleton problem has a clean fix at L2.** The 4-subgroup level leaves WNT (n = 71, ~4 events) as a structural-prior singleton with no neighbors to pool with. The 12-subtype level splits WNT into WNT-α and WNT-β, which share a parent group — the structural-prior mechanism has a non-degenerate target.
3. **Data is fully public and matched 1:1.** GSE85212 (methylation) and GSE85217 (expression) cover the same 763 samples; no dbGaP gating. CRAN package `MBMethPred` already wraps the Cavalli data for classification, lowering data-engineering cost.
4. **Two-block BIDIFAC+ scales down from Lock's four-block setup with no algorithm change.** Continuous-valued omics on both blocks; algebraically clean methodology mirror.

---

## Folder contents

```
Paper5_Medulloblastoma_StructuralPrior/
├── README.md                ← this file
├── PRE_REGISTRATION.md      ← locked pre-reg (methodology + hierarchy + decision rules)
├── data/                    ← Cavalli expression + methylation + clinical (after fetch)
├── reference/               ← subtype distribution CSV, target-covariate IDs, split indices
├── scripts/                 ← runners (fetch, BIDIFAC+, screening, Gibbs, diagnostics)
└── results/                 ← LPPD, bootstrap, per-subtype breakdown, CHECK 1/CHECK 2
```

Subfolders are scaffolded but empty at pre-reg time. Population happens in the order: data fetch → BIDIFAC+ → screening → covariate selection → split → Gibbs fires → diagnostics.

---

## Pre-reg highlights (load-bearing decisions)

- **Primary granularity: L2 12-subtype Cavalli.** L1 4-subgroup is the sensitivity row.
- **Target covariates: 3.** Selected by largest-effect-size from a BIDIFAC+ screening pass against L2. Identity committed to `reference/paper5_target_covariates.csv` at screening completion. Bakes in Paper 1's hard-won lesson that 4-cov reduced to 3-cov post-disposition.
- **Per-covariate CHECK 1 + CHECK 2 marginal-contribution cascade committed at pre-reg time**, not post-hoc. CHECK 1 = per-subgroup LPPD contribution; CHECK 2 = patient-level bootstrap of per-covariate marginal. Both reported alongside the joint LPPD verdict for each target covariate.
- **Held-out split: 80/20 stratified at L2, seed 20260515.** No k-fold expansion.
- **Gibbs parameters: 4 chains × 100k iters, mirror Paper 1.** Seeds 20260515–20260518 to avoid collision with Papers 1–4.
- **Decision rule: matches if |Δ| ≤ 2 nats OR CI brackets zero; outperforms if Δ > +2 AND CI excludes zero.** Mirror Paper 1.
- **Power-aware framing clause:** if actual event count from data pull is < 85 (≥ 50% below the estimated 170–210 range), the spec proceeds but the disposition pre-commits to "underpowered" framing for any matches verdict.

---

## What's not in Paper 5 (deferred)

- Sharma 2019 G3/G4 8-subtype L3 granularity → external-validation phase.
- MAGIC n = 898 → MAGIC requires dbGaP access + trial-pooling-with-covariate setup; deferred.
- Other within-cancer subtype taxonomies (DLBCL LymphGen, glioma WHO IDH/1p19q, breast PAM50) → separate pre-registrations per substrate.
- PFS endpoint → OS only.

---

## Position in the corpus

| paper | substrate | hierarchy | verdict |
|---|---|---|---|
| Paper 1 | 29-cancer TCGA pan-cancer | epi-class + germ-layer (external) | **3-cov MATCHES with parameter reduction** (post-disposition); 4-cov sensitivity +6.26 nats |
| Paper 2 | same as Paper 1 | operator-composition extension | FALSIFIES; L1 retro = cancer-identity-dominant |
| Paper 3 | breast Xenium + Chromium (MOFA-FLEX) | spatial-niche | FALSIFIED across 4 operator families; substrate boundary |
| Paper 4 Phase 1 | same as Paper 1 | 3 alternative pan-cancer hierarchies | PAPER1_HIERARCHY_SPECIFIC; 0/3 outperform |
| **Paper 5** | **Cavalli n = 763 medulloblastoma** | **within-cancer 12-subtype** | **TBD (this pre-reg)** |
