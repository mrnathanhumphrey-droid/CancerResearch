"""Candidate-independence check.

Yesterday's headline: "30 covariates have at least one non-Hoadley structural
hit at eta² > 0.30." But many of these come from related covariates (same
BIDIFAC+ block) testing the same cancer subset. They aren't 30 independent
confirmations.

Tests:
1. Cancer-set Jaccard overlap matrix among the 16 survivor (cov, hier, gran)
   tests. If most pairs share >50% of their cancers, the "16 survivors" are
   really only a few independent dimensions.
2. PCA on the survivor-covariates' PIP vectors (29 cancers x N_survivor_covs).
   How many components capture 90% of variance?
3. Block-membership of the survivors: how many distinct BIDIFAC+ blocks?
"""
from pathlib import Path
from itertools import combinations
import numpy as np
import pandas as pd

PIP_CSV = Path(r"C:/FkCancer/runs/run_paper_a_gibbs_2026-05-09/per_cancer_beta_pip.csv")
HIER_CSV = Path(r"C:/FkCancer/runs/run_screening_2026-05-09/hierarchies_long.csv")
OUT_CSV  = Path(r"C:/FkCancer/runs/run_screening_2026-05-09/falsification_independence.csv")

# 16 survivors (cov, hier, gran)
SURVIVORS = [
    ("26.1", "hoadley", "5a_4plus_unassigned"),
    ("3.2",  "sex_composition", "4_4bin"),
    ("8.2",  "epithelial_class", "2a_3bin"),
    ("20.1", "germ_layer", "3b_mid"),
    ("10.1", "epithelial_class", "2c_10bin"),
    ("3.1",  "hoadley", "5a_4plus_unassigned"),
    ("2.1",  "hoadley", "5c_iclusters"),
    ("31.1", "tissue_origin", "1b_11bin"),
    ("8.2",  "germ_layer", "3b_mid"),
    ("1.1",  "epithelial_class", "2c_10bin"),
    ("26.1", "epithelial_class", "2b_6bin"),
    ("36.1", "epithelial_class", "2a_3bin"),
    ("26.1", "epithelial_class", "2a_3bin"),
    ("22.1", "tissue_origin", "1b_11bin"),
    ("8.3",  "sex_composition", "4_4bin"),
    ("26.1", "epithelial_class", "2c_10bin"),
]


def main():
    pip_df = pd.read_csv(PIP_CSV)
    pip_df["cov_id"] = pip_df["cov_id"].astype(str)
    hier_df = pd.read_csv(HIER_CSV)
    if hier_df["contested"].dtype != bool:
        hier_df["contested"] = (hier_df["contested"].astype(str).str.upper()
                                == "TRUE")

    # 1. Jaccard overlap matrix
    print("=== 1. Cancer-set Jaccard overlap among survivors ===")
    cancer_sets = []
    labels = []
    for cov, hier, gran in SURVIVORS:
        sub_pip = pip_df[pip_df["cov_id"] == cov][["cancer"]]
        sub_h = hier_df[(hier_df["hierarchy"] == hier) &
                        (hier_df["granularity"] == gran) &
                        (~hier_df["contested"])]["cancer"]
        cancers = set(sub_pip["cancer"]) & set(sub_h)
        cancer_sets.append(cancers)
        labels.append(f"{cov}/{hier[:6]}/{gran[:6]}")

    n = len(SURVIVORS)
    jacc = np.zeros((n, n))
    for i, j in combinations(range(n), 2):
        a, b = cancer_sets[i], cancer_sets[j]
        if len(a | b) == 0:
            continue
        jacc[i, j] = jacc[j, i] = len(a & b) / len(a | b)
    np.fill_diagonal(jacc, 1.0)

    print("Jaccard overlap matrix (rows = survivors):")
    print("Idx  " + "  ".join(f"{i:5d}" for i in range(n)))
    for i in range(n):
        print(f"{i:3d}  " + "  ".join(f"{jacc[i,j]:.2f}" for j in range(n))
              + f"  ({labels[i]})")

    # Pairs with > 50% Jaccard overlap (NOT independent)
    high_overlap = [(i, j, jacc[i, j])
                    for i in range(n) for j in range(i+1, n)
                    if jacc[i, j] > 0.5]
    print(f"\n  Pairs with Jaccard > 0.5 (not independent): "
          f"{len(high_overlap)} / {n*(n-1)//2}")
    print(f"  Pairs with Jaccard = 1.0 (identical cancer sets): "
          f"{sum(1 for _, _, j in high_overlap if j == 1.0)}")

    # 2. Block diversity
    print("\n=== 2. BIDIFAC+ block diversity of survivors ===")
    blocks = sorted({int(float(cov)) for cov, _, _ in SURVIVORS})
    print(f"  Survivors span {len(blocks)} distinct blocks: {blocks}")
    # cov-id breakdown
    block_counts = {}
    for cov, _, _ in SURVIVORS:
        b = int(float(cov))
        block_counts.setdefault(b, []).append(cov)
    for b in sorted(block_counts):
        print(f"  Block {b}: {len(block_counts[b])} survivors "
              f"({set(block_counts[b])})")

    # 3. PCA on the survivor-covariates' PIP vectors (29 cancers x distinct covs)
    print("\n=== 3. PCA on survivor-covariate PIP vectors ===")
    distinct_covs = sorted({cov for cov, _, _ in SURVIVORS},
                           key=lambda x: float(x))
    X = pip_df[pip_df["cov_id"].isin(distinct_covs)]
    pivot = X.pivot(index="cancer", columns="cov_id", values="pip").fillna(0)
    print(f"  Pivot table shape: {pivot.shape} (cancers x covariates)")
    # Center, then SVD
    centered = pivot - pivot.mean(axis=0)
    U, S, Vt = np.linalg.svd(centered, full_matrices=False)
    var_explained = (S**2) / (S**2).sum()
    cumvar = np.cumsum(var_explained)
    n_for_90 = int(np.argmax(cumvar >= 0.9)) + 1
    print(f"  Singular values: {S.round(3)}")
    print(f"  Variance explained: {var_explained.round(3)}")
    print(f"  Cumulative: {cumvar.round(3)}")
    print(f"  Components for 90% variance: {n_for_90}")

    # Save summary
    rows = []
    rows.append({
        "metric": "n_survivors",
        "value": len(SURVIVORS),
    })
    rows.append({
        "metric": "n_distinct_blocks",
        "value": len(blocks),
    })
    rows.append({
        "metric": "n_distinct_cov_ids",
        "value": len(distinct_covs),
    })
    rows.append({
        "metric": "n_cov_x_hier_x_gran_with_jacc_eq_1",
        "value": sum(1 for _, _, j in high_overlap if j == 1.0),
    })
    rows.append({
        "metric": "pca_n_components_90pct",
        "value": n_for_90,
    })
    rows.append({
        "metric": "pca_pct_variance_first_component",
        "value": round(float(var_explained[0]), 3),
    })
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(OUT_CSV, index=False)
    print(f"\nWrote {OUT_CSV}")


if __name__ == "__main__":
    main()
