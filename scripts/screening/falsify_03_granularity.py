"""Permutation null for granularity-robustness rho.

Yesterday I reported tissue_origin (1a, 1b) Spearman rho = +0.72 between
granularities and called it evidence of "real structure" because both
granularities ranked the same covariates as most-explainable.

PROBLEM: granularity 1a is a NESTED aggregation of 1b. Every 1b bin is a
subset of one 1a bin. This means rankings induced by 1a and 1b are
mathematically correlated by construction, even on noise. The observed
rho=+0.72 may be entirely the nesting artifact.

TEST: For each covariate, permute the cancer-to-PIP assignment, recompute
the per-(hier, gran) eta², take the per-covariate max-eta-or-eta-at-each-gran,
compute Spearman rho between the two granularities. Repeat 1000x. If observed
rho is in the null distribution, the apparent robustness is just nesting.

We do this for tissue_origin (1a, 1b), epithelial_class (2a, 2b), (2b, 2c),
(2a, 2c), and germ_layer (3a, 3b).
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

PIP_CSV    = Path(r"C:/FkCancer/runs/run_paper_a_gibbs_2026-05-09/per_cancer_beta_pip.csv")
HIER_CSV   = Path(r"C:/FkCancer/runs/run_screening_2026-05-09/hierarchies_long.csv")
SCREEN_CSV = Path(r"C:/FkCancer/runs/run_screening_2026-05-09/screening_results.csv")
OUT_CSV    = Path(r"C:/FkCancer/runs/run_screening_2026-05-09/falsification_granularity.csv")

B = 1000
SEED = 20260510


def eta_sq(values, group_labels):
    grand = values.mean()
    ss_total = np.sum((values - grand) ** 2)
    if ss_total <= 0:
        return np.nan
    ss_between = 0.0
    for g in np.unique(group_labels):
        mask = group_labels == g
        if mask.any():
            ss_between += mask.sum() * (values[mask].mean() - grand) ** 2
    return ss_between / ss_total


def main():
    rng = np.random.default_rng(SEED)
    pip_df = pd.read_csv(PIP_CSV)
    pip_df["cov_id"] = pip_df["cov_id"].astype(str)
    hier_df = pd.read_csv(HIER_CSV)
    if hier_df["contested"].dtype != bool:
        hier_df["contested"] = (hier_df["contested"].astype(str).str.upper()
                                == "TRUE")
    screen = pd.read_csv(SCREEN_CSV)
    screen["cov_id"] = screen["cov_id"].astype(str)
    screen["eta_squared"] = pd.to_numeric(screen["eta_squared"], errors="coerce")

    pairs = [
        ("tissue_origin", "1a_5bin", "1b_11bin"),
        ("epithelial_class", "2a_3bin", "2b_6bin"),
        ("epithelial_class", "2a_3bin", "2c_10bin"),
        ("epithelial_class", "2b_6bin", "2c_10bin"),
        ("germ_layer", "3a_3bin", "3b_mid"),
    ]

    results = []
    for hier, g1, g2 in pairs:
        # Observed rho from the screening
        a = screen[(screen["hierarchy"] == hier) &
                   (screen["granularity"] == g1)].set_index("cov_id")["eta_squared"]
        b = screen[(screen["hierarchy"] == hier) &
                   (screen["granularity"] == g2)].set_index("cov_id")["eta_squared"]
        joint = pd.concat([a, b], axis=1, join="inner").dropna()
        joint.columns = ["g1", "g2"]
        rho_obs = spearmanr(joint["g1"], joint["g2"]).statistic
        n_covs = len(joint)
        common_covs = list(joint.index)

        # Build hierarchy assignments
        sub_h1 = hier_df[(hier_df["hierarchy"] == hier)
                         & (hier_df["granularity"] == g1)
                         & (~hier_df["contested"])][["cancer", "group"]]
        sub_h2 = hier_df[(hier_df["hierarchy"] == hier)
                         & (hier_df["granularity"] == g2)
                         & (~hier_df["contested"])][["cancer", "group"]]
        # Permutation: for each covariate, shuffle PIPs across cancers, recompute
        # eta at each granularity, then rho.
        perm_rhos = np.empty(B)
        for bi in range(B):
            etas_g1, etas_g2 = [], []
            covs_kept = []
            for cov in common_covs:
                cov_pips = (pip_df[pip_df["cov_id"] == cov]
                            [["cancer", "pip"]]
                            .set_index("cancer")["pip"])
                cancers_with_cov = list(cov_pips.index)
                shuffled = rng.permutation(cov_pips.values)
                shuffled_pips = pd.Series(shuffled, index=cancers_with_cov)
                # eta at g1
                merged1 = sub_h1[sub_h1["cancer"].isin(cancers_with_cov)]
                if len(merged1) >= 3 and merged1["group"].nunique() >= 2:
                    e1 = eta_sq(
                        np.array([shuffled_pips[c] for c in merged1["cancer"]]),
                        merged1["group"].to_numpy())
                else:
                    e1 = np.nan
                # eta at g2
                merged2 = sub_h2[sub_h2["cancer"].isin(cancers_with_cov)]
                if len(merged2) >= 3 and merged2["group"].nunique() >= 2:
                    e2 = eta_sq(
                        np.array([shuffled_pips[c] for c in merged2["cancer"]]),
                        merged2["group"].to_numpy())
                else:
                    e2 = np.nan
                if not (np.isnan(e1) or np.isnan(e2)):
                    etas_g1.append(e1)
                    etas_g2.append(e2)
                    covs_kept.append(cov)
            if len(etas_g1) >= 5:
                perm_rhos[bi] = spearmanr(etas_g1, etas_g2).statistic
            else:
                perm_rhos[bi] = np.nan

        perm_rhos = perm_rhos[~np.isnan(perm_rhos)]
        if len(perm_rhos) == 0:
            continue
        p_emp = float(np.mean(perm_rhos >= rho_obs))
        results.append({
            "hierarchy": hier,
            "granularities": f"({g1}, {g2})",
            "n_covs_compared": n_covs,
            "rho_obs": rho_obs,
            "rho_perm_median": float(np.median(perm_rhos)),
            "rho_perm_p95": float(np.quantile(perm_rhos, 0.95)),
            "p_empirical_one_sided": p_emp,
            "B_used": int(len(perm_rhos)),
            "interpretation": "REAL" if p_emp < 0.05 and rho_obs > np.median(perm_rhos)
                              else "ARTIFACT (nesting / chance)",
        })
        print(f"  {hier:20s} ({g1}, {g2})  rho_obs={rho_obs:+.3f}  "
              f"perm_med={np.median(perm_rhos):+.3f}  "
              f"perm_p95={np.quantile(perm_rhos, 0.95):+.3f}  "
              f"p_emp={p_emp:.3f}", flush=True)

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(results).to_csv(OUT_CSV, index=False)
    print(f"\nWrote {OUT_CSV}")


if __name__ == "__main__":
    main()
