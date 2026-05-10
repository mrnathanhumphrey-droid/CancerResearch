"""Cohort-size confound check.

Cancers have wildly different sample sizes (BRCA n=837 vs DLBC n=32). Bayesian
posteriors at small n inflate uncertainty, which can systematically affect PIPs
in the spike-and-slab. If small-n cancers cluster within particular tissue
bins (e.g., endocrine cancers are mostly small cohorts), the per-cancer PIP
profile may be detecting cohort size, not biology.

Tests:
1. Per-cancer mean(PIP across covariates) vs log(n_samples): do small-n cancers
   have systematically different mean PIPs?
2. Per-tissue-bin (1b_11bin) mean log(n_samples): do tissue bins differ
   significantly in their cohort-size composition? (ANOVA on log(n) by bin)
3. Re-run the top survivor tests using PARTIAL CORRELATION: regress PIP on
   log(n_samples) per cancer, then ANOVA on residuals against tissue. If the
   eta_squared survives, the signal is not cohort-confounded.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

PIP_CSV    = Path(r"C:/FkCancer/runs/run_paper_a_gibbs_2026-05-09/per_cancer_beta_pip.csv")
HIER_CSV   = Path(r"C:/FkCancer/runs/run_screening_2026-05-09/hierarchies_long.csv")
N_CSV      = Path(r"C:/FkCancer/runs/run_screening_2026-05-09/n_per_cancer.csv")
OUT_CSV    = Path(r"C:/FkCancer/runs/run_screening_2026-05-09/falsification_cohort.csv")


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
    pip_df = pd.read_csv(PIP_CSV)
    pip_df["cov_id"] = pip_df["cov_id"].astype(str)
    hier_df = pd.read_csv(HIER_CSV)
    if hier_df["contested"].dtype != bool:
        hier_df["contested"] = (hier_df["contested"].astype(str).str.upper()
                                == "TRUE")
    n_df = pd.read_csv(N_CSV).set_index("cancer")
    n_df["log_n"] = np.log(n_df["n_samples"])

    print("=== Test 1: Per-cancer mean(PIP) vs log(n_samples) ===")
    mean_pip = pip_df.groupby("cancer")["pip"].mean()
    df1 = pd.concat([mean_pip, n_df["log_n"]], axis=1).dropna()
    df1.columns = ["mean_pip", "log_n"]
    r, p = stats.pearsonr(df1["log_n"], df1["mean_pip"])
    rho, p_rho = stats.spearmanr(df1["log_n"], df1["mean_pip"])
    print(f"  Pearson r = {r:+.3f} (p={p:.3f})")
    print(f"  Spearman rho = {rho:+.3f} (p={p_rho:.3f})")
    print(f"  Range: log_n {df1['log_n'].min():.2f}-{df1['log_n'].max():.2f}, "
          f"mean_pip {df1['mean_pip'].min():.3f}-{df1['mean_pip'].max():.3f}")

    print("\n=== Test 2: Per-tissue-bin (1b_11bin) cohort-size composition ===")
    tissue_1b = hier_df[(hier_df["hierarchy"] == "tissue_origin") &
                       (hier_df["granularity"] == "1b_11bin") &
                       (~hier_df["contested"])]
    bin_n = tissue_1b.merge(n_df.reset_index(), on="cancer")
    print("\n  Bin / mean log(n) / cancers:")
    for bin_name, sub in bin_n.groupby("group"):
        print(f"    {bin_name:25s}  mean_log_n={sub['log_n'].mean():.2f}  "
              f"({', '.join(f'{c}({n})' for c, n in zip(sub['cancer'], sub['n_samples']))})")
    eta_bin = eta_sq(bin_n["log_n"].to_numpy(), bin_n["group"].to_numpy())
    f_bin, p_bin = stats.f_oneway(*[g["log_n"].values for _, g in
                                    bin_n.groupby("group")])
    print(f"\n  ANOVA log(n) by tissue_origin/1b_11bin: F={f_bin:.2f} "
          f"p={p_bin:.3g}  eta²={eta_bin:.3f}")

    if eta_bin > 0.30:
        print(f"  WARNING: tissue bins segregate by cohort size at eta²={eta_bin:.2f}")
    else:
        print(f"  OK: tissue bins do NOT systematically segregate by cohort size")

    print("\n=== Test 3: Partial-correlation eta² for top survivors ===")
    # For each survivor (cov, hier, gran), regress PIP on log(n) per-cancer,
    # then run ANOVA on residuals vs hierarchy bin.
    survivors = [
        ("8.2", "epithelial_class", "2a_3bin"),
        ("8.2", "germ_layer", "3b_mid"),
        ("10.1", "epithelial_class", "2c_10bin"),
        ("1.1", "epithelial_class", "2c_10bin"),
        ("20.1", "germ_layer", "3b_mid"),
        ("22.1", "tissue_origin", "1b_11bin"),
        ("3.2", "sex_composition", "4_4bin"),
        ("3.1", "hoadley", "5a_4plus_unassigned"),
    ]
    rows = []
    for cov_id, hier, gran in survivors:
        sub_pip = pip_df[pip_df["cov_id"] == cov_id][["cancer", "pip"]]
        sub_h = hier_df[(hier_df["hierarchy"] == hier) &
                        (hier_df["granularity"] == gran) &
                        (~hier_df["contested"])][["cancer", "group"]]
        merged = sub_pip.merge(sub_h, on="cancer").merge(
            n_df.reset_index(), on="cancer")
        if len(merged) < 4:
            continue
        # Original eta²
        eta_orig = eta_sq(merged["pip"].to_numpy(), merged["group"].to_numpy())
        # Residual after regressing PIP on log_n
        slope, intercept, r_val, p_val, _ = stats.linregress(
            merged["log_n"], merged["pip"])
        residuals = merged["pip"] - (slope * merged["log_n"] + intercept)
        eta_resid = eta_sq(residuals.to_numpy(), merged["group"].to_numpy())
        rows.append({
            "cov_id": cov_id,
            "hierarchy": hier,
            "granularity": gran,
            "n": len(merged),
            "PIP_log_n_pearson_r": r_val,
            "p_PIP_log_n": p_val,
            "eta_orig": eta_orig,
            "eta_after_log_n_partialing": eta_resid,
            "delta_eta": eta_resid - eta_orig,
        })
    df3 = pd.DataFrame(rows)
    print(df3.to_string(index=False))

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    df3.to_csv(OUT_CSV, index=False)
    print(f"\nWrote {OUT_CSV}")
    return df1, df3


if __name__ == "__main__":
    main()
