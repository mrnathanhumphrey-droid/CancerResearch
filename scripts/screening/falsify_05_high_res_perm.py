"""High-resolution permutation null on the 16 survivor cells (B=10000)
to tighten p-values to ~0.0001 resolution.

Yesterday's 1000-perm pass left ambiguity for the lowest-p hits (some came back
as 0.000 = "less than 1 in 1000" but Bonferroni-correction at alpha=0.05 with
319 tests requires p < 0.000157, which we can only verify with B>=10000).
"""
from pathlib import Path
import numpy as np
import pandas as pd

PIP_CSV    = Path(r"C:/FkCancer/runs/run_paper_a_gibbs_2026-05-09/per_cancer_beta_pip.csv")
HIER_CSV   = Path(r"C:/FkCancer/runs/run_screening_2026-05-09/hierarchies_long.csv")
OUT_CSV    = Path(r"C:/FkCancer/runs/run_screening_2026-05-09/falsification_high_res_perm.csv")

B = 10000
SEED = 20260511

# 16 survivors from B=1000 pass (cov, hier, gran)
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

    rows = []
    for cov_id, hier, gran in SURVIVORS:
        sub_pip = pip_df[pip_df["cov_id"] == cov_id][["cancer", "pip"]]
        sub_h = hier_df[(hier_df["hierarchy"] == hier) &
                        (hier_df["granularity"] == gran) &
                        (~hier_df["contested"])][["cancer", "group"]]
        merged = sub_pip.merge(sub_h, on="cancer")
        if len(merged) < 3 or merged["group"].nunique() < 2:
            continue
        values = merged["pip"].to_numpy()
        labels = merged["group"].to_numpy()
        eta_obs = eta_sq(values, labels)
        eta_perm = np.empty(B)
        for b in range(B):
            eta_perm[b] = eta_sq(rng.permutation(values), labels)
        eta_perm = eta_perm[~np.isnan(eta_perm)]
        p_emp = float(np.mean(eta_perm >= eta_obs))
        # Bonferroni threshold = 0.05/319 = 0.000157
        bonf_pass = p_emp < 0.000157
        rows.append({
            "cov_id": cov_id,
            "hierarchy": hier,
            "granularity": gran,
            "n_cancers": int(len(merged)),
            "eta_obs": float(eta_obs),
            "perm_median": float(np.median(eta_perm)),
            "perm_p99": float(np.quantile(eta_perm, 0.99)),
            "perm_p999": float(np.quantile(eta_perm, 0.999)),
            "p_empirical": p_emp,
            "bonf_threshold_0.05_319": 0.05 / 319,
            "bonf_pass": "YES" if bonf_pass else "no",
            "B": int(len(eta_perm)),
        })
        print(f"  {cov_id} {hier:20s} {gran:25s}  n={len(merged):3d}  "
              f"eta_obs={eta_obs:.3f}  p_emp={p_emp:.5f}  "
              f"bonf_pass={'YES' if bonf_pass else 'no'}", flush=True)

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(OUT_CSV, index=False)
    print(f"\nWrote {OUT_CSV} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
