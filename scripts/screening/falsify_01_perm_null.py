"""Permutation null for every (cov, hier, gran) screening row.

For each test row:
  - Holds fixed: which cancers are tested, their bin assignments (group), n_groups
  - Permutes: which PIP value belongs to which cancer (shuffle PIPs across cancers
    that have this covariate present)
  - Recomputes eta_squared B times
  - Empirical p = (#perms with eta_perm >= eta_obs) / B

This directly tests: "does the cancer-to-PIP assignment matter, holding fixed
the structural-bin distribution of the cancers in the test?" If observed eta²
is high purely because of the bin-imbalance of the cancer set (e.g., 6 cancers
all from Pan-GI), permutation will reproduce the same eta² under shuffled PIPs.
If observed eta² is driven by REAL bin-vs-PIP signal, permutation eta² will be
much lower.

Output: falsification_perm_null.csv columns:
  cov_id, hierarchy, granularity, n_cancers_used, eta_obs, eta_perm_median,
  eta_perm_p95, p_empirical (one-sided), B (n perms run)
"""
import csv
from pathlib import Path
import numpy as np
import pandas as pd

PIP_CSV    = Path(r"C:/FkCancer/runs/run_paper_a_gibbs_2026-05-09/per_cancer_beta_pip.csv")
HIER_CSV   = Path(r"C:/FkCancer/runs/run_screening_2026-05-09/hierarchies_long.csv")
SCREEN_CSV = Path(r"C:/FkCancer/runs/run_screening_2026-05-09/screening_results.csv")
OUT_CSV    = Path(r"C:/FkCancer/runs/run_screening_2026-05-09/falsification_perm_null.csv")

B = 1000
SEED = 20260509


def eta_sq(values, group_labels):
    """Compute eta-squared given a values vector and an aligned group label vector."""
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
    # pandas may auto-cast contested -> bool
    if hier_df["contested"].dtype != bool:
        hier_df["contested"] = (hier_df["contested"].astype(str).str.upper()
                                == "TRUE")
    screen = pd.read_csv(SCREEN_CSV)
    screen["cov_id"] = screen["cov_id"].astype(str)

    rows = []
    n_total = len(screen)
    for i, row in screen.iterrows():
        cov_id  = row["cov_id"]
        hier    = row["hierarchy"]
        gran    = row["granularity"]
        eta_obs = row["eta_squared"]
        if pd.isna(eta_obs):
            rows.append({
                "cov_id": cov_id, "hierarchy": hier, "granularity": gran,
                "n_cancers_used": int(row["n_cancers_used"]),
                "eta_obs": "", "eta_perm_median": "", "eta_perm_p95": "",
                "p_empirical": "", "B": 0,
            })
            continue
        # Get the cancers and their PIPs for this covariate
        sub_pip = pip_df[pip_df["cov_id"] == cov_id][["cancer", "pip"]]
        # Group assignment for these cancers (drop contested)
        sub_h = hier_df[(hier_df["hierarchy"] == hier)
                        & (hier_df["granularity"] == gran)
                        & (~hier_df["contested"])]
        merged = sub_pip.merge(sub_h, on="cancer", how="inner")
        if len(merged) < 3:
            rows.append({
                "cov_id": cov_id, "hierarchy": hier, "granularity": gran,
                "n_cancers_used": len(merged),
                "eta_obs": float(eta_obs) if pd.notna(eta_obs) else "",
                "eta_perm_median": "", "eta_perm_p95": "", "p_empirical": "",
                "B": 0,
            })
            continue
        values = merged["pip"].to_numpy()
        labels = merged["group"].to_numpy()
        # Permutation null: shuffle PIPs among cancers, hold labels fixed
        eta_perm = np.empty(B)
        for b in range(B):
            shuffled = rng.permutation(values)
            eta_perm[b] = eta_sq(shuffled, labels)
        eta_perm = eta_perm[~np.isnan(eta_perm)]
        if len(eta_perm) == 0:
            rows.append({
                "cov_id": cov_id, "hierarchy": hier, "granularity": gran,
                "n_cancers_used": len(merged),
                "eta_obs": float(eta_obs), "eta_perm_median": "",
                "eta_perm_p95": "", "p_empirical": "", "B": 0,
            })
            continue
        p_emp = float(np.mean(eta_perm >= eta_obs))
        rows.append({
            "cov_id": cov_id, "hierarchy": hier, "granularity": gran,
            "n_cancers_used": len(merged),
            "eta_obs": float(eta_obs),
            "eta_perm_median": float(np.median(eta_perm)),
            "eta_perm_p95": float(np.quantile(eta_perm, 0.95)),
            "p_empirical": p_emp,
            "B": int(len(eta_perm)),
        })
        if (i + 1) % 100 == 0:
            print(f"  {i+1}/{n_total}", flush=True)

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    cols = ["cov_id", "hierarchy", "granularity", "n_cancers_used",
            "eta_obs", "eta_perm_median", "eta_perm_p95", "p_empirical", "B"]
    with OUT_CSV.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    print(f"\nWrote {OUT_CSV} ({len(rows)} rows)")

    # Quick diagnostic
    df = pd.DataFrame(rows)
    df = df[df["p_empirical"] != ""].copy()
    df["p_empirical"] = pd.to_numeric(df["p_empirical"])
    df["eta_obs"] = pd.to_numeric(df["eta_obs"])
    df["hoadley"] = df["hierarchy"] == "hoadley"

    print("\n=== Empirical p-value distribution ===")
    for label, sub in [("ALL", df), ("Non-Hoadley", df[~df["hoadley"]]),
                       ("Hoadley", df[df["hoadley"]])]:
        n = len(sub)
        if n == 0:
            continue
        n_sig05 = int((sub["p_empirical"] < 0.05).sum())
        n_sig01 = int((sub["p_empirical"] < 0.01).sum())
        expected_05 = n * 0.05
        expected_01 = n * 0.01
        print(f"  {label:15s}  n={n:3d}   p<0.05: {n_sig05}/{n} "
              f"(expect by chance: {expected_05:.1f})   "
              f"p<0.01: {n_sig01}/{n} (expect: {expected_01:.2f})")

    # Hits at η²>0.30 with empirical p<0.05
    print("\n=== High-η² hits that survive permutation null (p<0.05) ===")
    high = df[(df["eta_obs"] > 0.30)]
    high_sig = high[high["p_empirical"] < 0.05]
    print(f"  Of {len(high)} rows with η²>0.30, {len(high_sig)} have p<0.05")
    print(f"  Of those: {(~high_sig['hoadley']).sum()} non-Hoadley, "
          f"{high_sig['hoadley'].sum()} Hoadley")


if __name__ == "__main__":
    main()
