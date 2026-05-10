"""ANOVA screening of per-cancer PIPs against candidate structural hierarchies.

For each covariate × hierarchy × granularity:
  - Drop contested cancers (per the reference table)
  - Group remaining cancers by their structural-feature assignment
  - Run one-way ANOVA on PIP across groups
  - Report F, p, eta_squared, n_groups, n_singleton_groups, n_cancers_used,
    hoadley_flag, n_groups_used (groups with >=1 cancer after intersection)

Output:
  screening_results.csv   — long format, every (cov, hier, gran) test
  data_handling_log.md    — audit of design choices

Statistical notes:
  - eta_squared = SS_between / SS_total. With small group counts and small n,
    eta_squared is biased upward; we report it as the screening metric per the
    user brief, supplemented by F and p but those should not be the basis for
    interpretation.
  - Singleton groups (n=1 in a bin) contribute to between-group variance via
    their (group_mean - grand_mean)^2 term but provide no within-group estimate.
    We report n_singleton_groups; a row with n_singleton_groups > n_groups/2 is
    flagged degenerate.
  - When n_groups < 2 after the contested-filter intersection, we skip the test
    (not enough variance partitions to test).
"""
import csv
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

PIP_CSV  = Path(r"C:/FkCancer/runs/run_paper_a_gibbs_2026-05-09/per_cancer_beta_pip.csv")
HIER_CSV = Path(r"C:/FkCancer/runs/run_screening_2026-05-09/hierarchies_long.csv")
OUT_CSV  = Path(r"C:/FkCancer/runs/run_screening_2026-05-09/screening_results.csv")

HIERARCHIES_HOADLEY = {"hoadley"}


def eta_squared(values_by_group):
    """Compute eta-squared = SS_between / SS_total for one-way layout."""
    flat = np.concatenate(values_by_group)
    grand_mean = flat.mean()
    ss_total = np.sum((flat - grand_mean) ** 2)
    ss_between = sum(
        len(g) * (g.mean() - grand_mean) ** 2 for g in values_by_group)
    if ss_total <= 0:
        return np.nan
    return ss_between / ss_total


def main():
    pip_df = pd.read_csv(PIP_CSV)
    hier_df = pd.read_csv(HIER_CSV)
    print(f"PIP rows: {len(pip_df)}  unique covs: {pip_df['cov_id'].nunique()}"
          f"  unique cancers: {pip_df['cancer'].nunique()}")
    print(f"Hierarchy rows: {len(hier_df)}  hier×gran cells: "
          f"{hier_df.groupby(['hierarchy','granularity']).ngroups}")

    cov_ids = sorted(pip_df["cov_id"].unique(), key=lambda x: float(x))
    cells   = list(hier_df.groupby(["hierarchy", "granularity"]).groups.keys())

    results = []
    for cov in cov_ids:
        cov_pips = pip_df[pip_df["cov_id"] == cov][["cancer", "pip"]]
        cov_pips = cov_pips.set_index("cancer")["pip"]
        cancers_with_cov = set(cov_pips.index)
        for hier, gran in cells:
            sub = hier_df[(hier_df["hierarchy"] == hier)
                          & (hier_df["granularity"] == gran)
                          & (hier_df["contested"] == False)]  # drop contested
            sub = sub[sub["cancer"].isin(cancers_with_cov)]
            if len(sub) == 0:
                continue
            # Group PIPs
            groups = {}
            for _, r in sub.iterrows():
                groups.setdefault(r["group"], []).append(cov_pips[r["cancer"]])
            n_groups = len(groups)
            n_singleton = sum(1 for v in groups.values() if len(v) == 1)
            n_cancers_used = len(sub)
            n_groups_total = hier_df[(hier_df["hierarchy"] == hier)
                                     & (hier_df["granularity"] == gran)
                                     & (hier_df["contested"] == False)
                                     ]["group"].nunique()
            arrays = [np.array(v) for v in groups.values()]
            if n_groups < 2 or n_cancers_used < 3:
                F, pval, eta2 = (np.nan, np.nan, np.nan)
                degenerate = True
            else:
                try:
                    F, pval = stats.f_oneway(*arrays)
                except Exception:
                    F, pval = (np.nan, np.nan)
                eta2 = eta_squared(arrays)
                degenerate = (n_singleton > n_groups / 2.0)
            results.append({
                "cov_id": cov,
                "hierarchy": hier,
                "granularity": gran,
                "n_groups": n_groups,
                "n_groups_total": n_groups_total,
                "n_singleton_groups": n_singleton,
                "n_cancers_used": n_cancers_used,
                "F": F if F == F else "",
                "p_value": pval if pval == pval else "",
                "eta_squared": eta2 if eta2 == eta2 else "",
                "hoadley_flag": "TRUE" if hier in HIERARCHIES_HOADLEY else "FALSE",
                "degenerate_majority_singletons": "TRUE" if degenerate else "FALSE",
            })

    cols = ["cov_id", "hierarchy", "granularity", "n_groups", "n_groups_total",
            "n_singleton_groups", "n_cancers_used", "F", "p_value",
            "eta_squared", "hoadley_flag", "degenerate_majority_singletons"]
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with OUT_CSV.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(results)
    print(f"Wrote {OUT_CSV} ({len(results)} rows)")

    # Quick diagnostic
    df = pd.DataFrame(results)
    df_eta = df[df["eta_squared"] != ""].copy()
    df_eta["eta_squared"] = pd.to_numeric(df_eta["eta_squared"])
    print("\n=== Top 10 by eta² (any hierarchy) ===")
    top = df_eta.nlargest(10, "eta_squared")[
        ["cov_id", "hierarchy", "granularity", "n_groups", "n_cancers_used",
         "eta_squared", "hoadley_flag", "degenerate_majority_singletons"]]
    print(top.to_string(index=False))
    print("\n=== eta² distribution by hierarchy (median, mean, max) ===")
    by_h = df_eta.groupby("hierarchy")["eta_squared"].agg(
        ["count", "median", "mean", "max"])
    print(by_h.round(4))


if __name__ == "__main__":
    main()
