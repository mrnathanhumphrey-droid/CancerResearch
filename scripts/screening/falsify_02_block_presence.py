"""BIDIFAC+ block-presence confound check.

Hypothesis we're trying to falsify: the per-cancer PIP profiles reveal real
external structure (tissue origin) in the cross-cancer borrowing pattern.

ALTERNATIVE: the per-cancer PIP profiles are just re-detecting BIDIFAC+'s own
tissue-stratified module structure. BIDIFAC+ produces blocks of latent factors
on the 4 omics. Some blocks are tissue-specific (active only in certain cancer
types). If block PRESENCE (which cancers have any covariate from a block) is
itself tissue-clustered, then the per-cancer spike-and-slab PIPs inherit that
tissue structure trivially — the "structural decomposition" is circular through
BIDIFAC+.

Test: For each BIDIFAC+ block (integer prefix of cov_id), construct the
binary "presence vector" over 29 cancers: which cancers have any covariate
from this block. Then run ANOVA of presence-vector vs each hierarchy.
If block-presence eta_squared is COMPARABLE to the screening eta_squared on
the same hierarchy, the structure was baked in by BIDIFAC+ before any
spike-and-slab fit.

We compare for each block, the block-presence eta_squared against the average
PIP-screening eta_squared on covariates from that block.

Output: falsification_block_presence.csv
"""
import csv
from pathlib import Path
import numpy as np
import pandas as pd

PIP_CSV    = Path(r"C:/FkCancer/runs/run_paper_a_gibbs_2026-05-09/per_cancer_beta_pip.csv")
HIER_CSV   = Path(r"C:/FkCancer/runs/run_screening_2026-05-09/hierarchies_long.csv")
SCREEN_CSV = Path(r"C:/FkCancer/runs/run_screening_2026-05-09/screening_results.csv")
OUT_CSV    = Path(r"C:/FkCancer/runs/run_screening_2026-05-09/falsification_block_presence.csv")


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


def block_id(cov_id):
    s = str(cov_id)
    return int(float(s))


def main():
    pip_df = pd.read_csv(PIP_CSV)
    pip_df["cov_id"] = pip_df["cov_id"].astype(str)
    pip_df["block"] = pip_df["cov_id"].apply(block_id)

    hier_df = pd.read_csv(HIER_CSV)
    if hier_df["contested"].dtype != bool:
        hier_df["contested"] = (hier_df["contested"].astype(str).str.upper()
                                == "TRUE")

    screen = pd.read_csv(SCREEN_CSV)
    screen["cov_id"] = screen["cov_id"].astype(str)
    screen["block"] = screen["cov_id"].apply(block_id)
    screen["eta_squared"] = pd.to_numeric(screen["eta_squared"], errors="coerce")

    cancers = sorted(pip_df["cancer"].unique())
    blocks = sorted(pip_df["block"].unique())
    print(f"Cancers: {len(cancers)}, BIDIFAC+ blocks: {len(blocks)}")

    # 1. Build cancer-x-block presence matrix (binary)
    presence = pd.DataFrame(0, index=cancers, columns=blocks)
    for _, r in pip_df.iterrows():
        presence.loc[r["cancer"], r["block"]] = 1

    print("\nBlock presence (n cancers per block):")
    print(presence.sum(axis=0).to_string())

    # 2. For each block × hierarchy × granularity, run eta_sq on presence vector
    cells = list(hier_df.groupby(["hierarchy", "granularity"]).groups.keys())
    rows = []
    for block in blocks:
        present_vec = presence[block].to_numpy()  # 0/1 per cancer
        cancer_idx = {c: i for i, c in enumerate(cancers)}
        for hier, gran in cells:
            sub = hier_df[(hier_df["hierarchy"] == hier)
                          & (hier_df["granularity"] == gran)
                          & (~hier_df["contested"])]
            if len(sub) == 0:
                continue
            cancers_in_test = [c for c in sub["cancer"] if c in cancer_idx]
            if len(cancers_in_test) < 3:
                continue
            values = np.array([present_vec[cancer_idx[c]] for c in cancers_in_test])
            labels = np.array([
                sub[sub["cancer"] == c]["group"].iloc[0]
                for c in cancers_in_test])
            if len(np.unique(labels)) < 2:
                continue
            eta = eta_sq(values.astype(float), labels)
            # Also: average PIP-screening eta_squared on covariates from this block,
            # at the same hierarchy/granularity
            sub_screen = screen[(screen["block"] == block)
                                & (screen["hierarchy"] == hier)
                                & (screen["granularity"] == gran)]
            avg_pip_eta = sub_screen["eta_squared"].dropna().mean()
            n_covs_in_block = sub_screen["eta_squared"].dropna().count()
            rows.append({
                "block": int(block),
                "hierarchy": hier,
                "granularity": gran,
                "n_cancers_with_block": int(present_vec.sum()),
                "n_cancers_in_test": int(len(cancers_in_test)),
                "block_presence_eta": eta,
                "avg_pip_screening_eta": avg_pip_eta,
                "n_covs_in_block": int(n_covs_in_block),
            })

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows)
    df.to_csv(OUT_CSV, index=False)
    print(f"\nWrote {OUT_CSV} ({len(df)} rows)")

    # Diagnostics
    print("\n=== Block-presence eta vs PIP-screening eta (by hierarchy) ===")
    for hier in df["hierarchy"].unique():
        sub = df[df["hierarchy"] == hier].dropna(
            subset=["block_presence_eta", "avg_pip_screening_eta"])
        if len(sub) == 0:
            continue
        # Spearman correlation between the two etas across blocks
        if len(sub) >= 5:
            rho, p = (sub[["block_presence_eta", "avg_pip_screening_eta"]]
                      .corr(method="spearman").iloc[0, 1], None)
        else:
            rho = float("nan")
        n_coupled = ((sub["block_presence_eta"] > 0.30)
                     & (sub["avg_pip_screening_eta"] > 0.30)).sum()
        print(f"  {hier:20s}  rows={len(sub):3d}  rho={rho:+.3f}  "
              f"n_block_eta>0.30 AND avg_pip_eta>0.30: {n_coupled}")

    # The smoking gun: are the high-PIP-eta hits exactly the high-block-eta blocks?
    print("\n=== Top-15 (block, hier, gran) by block-presence eta ===")
    top = df.sort_values("block_presence_eta", ascending=False).head(15)
    print(top[["block", "hierarchy", "granularity", "block_presence_eta",
               "avg_pip_screening_eta", "n_covs_in_block"]].to_string(index=False))


if __name__ == "__main__":
    main()
