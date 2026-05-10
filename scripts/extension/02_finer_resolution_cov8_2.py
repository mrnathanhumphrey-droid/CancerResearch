"""Finer-resolution prediction test for cov 8.2.
DISCOVERY work, not validation. ANOVA-only. No Gibbs.

Test design:
  1. Pull cov 8.2 per-cancer beta_mean from per_cancer_beta_pip.csv
  2. Compute group mean within each epi/2a_3bin group
  3. Compute residuals (beta_mean - group_mean) per cancer
  4. Drop cancers contested at epi/2a (MESO, UCS)
  5. PRIMARY (Option A): within each epi/2a group, ANOVA residuals
     vs tissue_origin/1b_11bin label
  6. SENSITIVITY (Option C): within each epi/2a group, ANOVA residuals
     vs hybrid label (subtype for multi-cancer tissues; tissue otherwise)
  7. Permutation null per group (B=1000)

Subtype assignments are PRE-COMMITTED before any computation.
"""
import csv
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

PIP_CSV = Path(r"C:/FkCancer/runs/run_paper_a_gibbs_2026-05-09/per_cancer_beta_pip.csv")
OUT_DIR = Path(r"C:/FkCancer/runs/run_extension_2026-05-10")
OUT_DIR.mkdir(parents=True, exist_ok=True)

B = 1000
SEED = 20260510

# ============================================================
# PRE-COMMITTED ASSIGNMENTS (locked before computing)
# ============================================================

# epi/2a_3bin
EPI_2A = {
    "ACC":"Epithelial","BLCA":"Epithelial","BRCA":"Epithelial","CESC":"Epithelial",
    "CHOL":"Epithelial","CORE":"Epithelial","DLBC":"Hematological","ESCA":"Epithelial",
    "HNSC":"Epithelial","KICH":"Epithelial","KIRC":"Epithelial","KIRP":"Epithelial",
    "LGG":"Non-epithelial","LIHC":"Epithelial","LUAD":"Epithelial","LUSC":"Epithelial",
    "MESO":None,"OV":"Epithelial","PAAD":"Epithelial",
    "PCPG":"Non-epithelial","PRAD":"Epithelial","SARC":"Non-epithelial",
    "SKCM":"Non-epithelial","STAD":"Epithelial","TGCT":"Non-epithelial",
    "THCA":"Epithelial","THYM":"Epithelial","UCEC":"Epithelial","UCS":None,
}

# tissue_origin/1b_11bin
TISSUE_1B = {
    "ACC":"Endocrine","BLCA":"Male repro+bladder","BRCA":"Breast",
    "CESC":"Gynecological","CHOL":"Hepatobiliary+pancreas","CORE":"GI tract",
    "DLBC":"Soft+heme+H&N","ESCA":"GI tract","HNSC":"Soft+heme+H&N",
    "KICH":"Kidney","KIRC":"Kidney","KIRP":"Kidney","LGG":"CNS",
    "LIHC":"Hepatobiliary+pancreas","LUAD":"Lung+pleura","LUSC":"Lung+pleura",
    "MESO":"Lung+pleura","OV":"Gynecological","PAAD":"Hepatobiliary+pancreas",
    "PCPG":"Endocrine","PRAD":"Male repro+bladder","SARC":"Soft+heme+H&N",
    "SKCM":"Skin","STAD":"GI tract","TGCT":"Male repro+bladder",
    "THCA":"Endocrine","THYM":"Endocrine","UCEC":"Gynecological","UCS":"Gynecological",
}

# Modal molecular subtypes — PRE-COMMITTED
# Sources: TCGA flagship papers per cancer; modal = dominant subtype in TCGA cohort.
# Where modal is uncertain or no canonical TCGA subtypes exist, the cancer name
# is used as the subtype label (effectively makes that cancer its own bin).
MODAL_SUBTYPE = {
    # Lung+pleura
    "LUAD": "TRU",                   # Terminal Respiratory Unit (TCGA LUAD modal)
    "LUSC": "Classical",             # Wilkerson 2010 modal
    "MESO": "MESO",                  # MESO contested for epi/2a anyway; not used
    # Kidney
    "KICH": "Chromophobe",           # Only one chromophobe subtype
    "KIRC": "ccA",                   # ClearCode34 ccA modal
    "KIRP": "Type1",                 # Papillary Type 1 modal
    # Hepatobiliary+pancreas
    "LIHC": "LIHC",                  # iCluster labels heterogeneous; use cancer name
    "CHOL": "CHOL",                  # No canonical subtype distinction at TCGA size
    "PAAD": "Classical",             # Bailey/Moffitt Classical modal
    # GI tract
    "ESCA": "ESCC1",                 # ESCA TCGA mostly squamous; ESCC1 modal
    "STAD": "CIN",                   # Chromosomal instability modal
    "CORE": "CMS2",                  # Consensus Molecular Subtype 2 (canonical) modal
    # Gynecological
    "OV":   "Mesenchymal",           # Verhaak/Tothill modal
    "UCEC": "Endometrioid",          # Endometrioid majority
    "UCS":  "UCS",                   # UCS contested for epi/2a anyway
    "CESC": "Squamous_HPV+",         # CESC TCGA modal
    # Male repro+bladder
    "PRAD": "ERG_fusion+",           # ERG fusion modal
    "TGCT": "Seminoma",              # Seminoma vs nonseminoma; seminoma modal
    "BLCA": "Luminal",               # TCGA BLCA luminal modal
    # Endocrine
    "ACC":  "ACC",                   # CIMP clusters; uncertain modal — use cancer name
    "PCPG": "Pseudohypoxia",         # Cluster 1 (pseudohypoxia) modal
    "THCA": "BRAF_papillary",        # BRAF-mutant papillary modal
    "THYM": "B2",                    # WHO type B2 modal
    # Soft+heme+H&N
    "SARC": "LMS",                   # Leiomyosarcoma modal
    "DLBC": "GCB",                   # Germinal Center B-cell modal
    "HNSC": "Classical",             # HPV-negative classical modal
    # Single-cancer tissues — subtype not used in Option C
    "BRCA": "LumA",                  # PAM50 Luminal A modal (would be used if not single-tissue)
    "LGG":  "IDH1mut",
    "SKCM": "BRAF_mut",
}

# Compose Option C (hybrid) label per cancer:
#   If tissue contains >1 TCGA cancer: use modal subtype as label
#   Else (single-cancer tissue): use tissue label
def build_hybrid_labels(epi_2a, tissue_1b, modal_subtype):
    tissue_counts = {}
    for c, t in tissue_1b.items():
        tissue_counts[t] = tissue_counts.get(t, 0) + 1
    hybrid = {}
    for c in epi_2a:
        t = tissue_1b[c]
        if tissue_counts[t] > 1:
            hybrid[c] = modal_subtype[c]
        else:
            hybrid[c] = t
    return hybrid


HYBRID = build_hybrid_labels(EPI_2A, TISSUE_1B, MODAL_SUBTYPE)

# ============================================================
# ANOVA + permutation null
# ============================================================

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


def perm_null_eta(values, labels, B, rng):
    out = np.empty(B)
    for b in range(B):
        out[b] = eta_sq(rng.permutation(values), labels)
    out = out[~np.isnan(out)]
    return out


def main():
    rng = np.random.default_rng(SEED)

    # 1. Load cov 8.2 per-cancer beta_means
    pip_df = pd.read_csv(PIP_CSV)
    pip_df["cov_id"] = pip_df["cov_id"].astype(str)
    cov82 = pip_df[pip_df["cov_id"] == "8.2"][["cancer", "beta_mean"]].copy()
    cov82 = cov82.sort_values("cancer").reset_index(drop=True)
    print(f"Cov 8.2 per-cancer rows: {len(cov82)}")

    # 2. Annotate epi/2a, tissue, hybrid labels; drop contested
    cov82["epi_2a"]    = cov82["cancer"].map(EPI_2A)
    cov82["tissue_1b"] = cov82["cancer"].map(TISSUE_1B)
    cov82["hybrid_C"]  = cov82["cancer"].map(HYBRID)
    contested = cov82[cov82["epi_2a"].isna()]
    print(f"Dropped contested for epi/2a: {list(contested['cancer'])}")
    cov82 = cov82[cov82["epi_2a"].notna()].reset_index(drop=True)
    print(f"Cancers retained: {len(cov82)}")

    # 3. Group means and residuals
    cov82["group_mean"] = cov82.groupby("epi_2a")["beta_mean"].transform("mean")
    cov82["residual"]   = cov82["beta_mean"] - cov82["group_mean"]
    cov82.to_csv(OUT_DIR / "finer_resolution_residuals.csv", index=False)
    print(f"Wrote finer_resolution_residuals.csv ({len(cov82)} rows)")

    # 4. Per-group: list of cancers + tissue labels + hybrid labels
    print("\n=== Per epi/2a group ===")
    for g in cov82["epi_2a"].unique():
        sub = cov82[cov82["epi_2a"] == g]
        print(f"\n  Group: {g} (n={len(sub)})")
        print(f"    group_mean(beta) = {sub['group_mean'].iloc[0]:+.4f}")
        for _, r in sub.iterrows():
            print(f"      {r['cancer']:5s}  beta={r['beta_mean']:+.4f}  "
                  f"resid={r['residual']:+.4f}  tissue={r['tissue_1b']:25s}  "
                  f"hybrid={r['hybrid_C']}")

    # ============================================================
    # PRIMARY (Option A): residuals ~ tissue, within each epi/2a group
    # ============================================================
    print("\n\n=== PRIMARY (Option A): residuals ~ tissue_origin/1b_11bin within each epi/2a group ===")
    rows_A = []
    for g in cov82["epi_2a"].unique():
        sub = cov82[cov82["epi_2a"] == g]
        n = len(sub)
        labels = sub["tissue_1b"].to_numpy()
        values = sub["residual"].to_numpy()
        # n_singletons = #tissue labels appearing exactly once
        from collections import Counter
        cnt = Counter(labels)
        n_groups = len(cnt)
        n_single = sum(1 for v in cnt.values() if v == 1)
        n_multi  = sum(1 for v in cnt.values() if v >= 2)
        eta_obs = eta_sq(values, labels) if n_groups >= 2 else np.nan
        if n_groups >= 2 and n >= 3:
            null_etas = perm_null_eta(values, labels, B, rng)
            null_med = float(np.median(null_etas))
            null_p95 = float(np.quantile(null_etas, 0.95))
            p_emp   = float(np.mean(null_etas >= eta_obs)) if not np.isnan(eta_obs) else np.nan
            try:
                F, pval = stats.f_oneway(*[values[labels == lab] for lab in cnt])
            except Exception:
                F, pval = (np.nan, np.nan)
        else:
            null_med = null_p95 = p_emp = np.nan
            F = pval = np.nan
        underpowered = (n < 4) or (n_single > n_groups / 2)
        rows_A.append({
            "epi_2a_group": g,
            "n_cancers": n,
            "n_tissues": n_groups,
            "n_singleton_tissues": n_single,
            "n_multi_tissues": n_multi,
            "eta_obs": eta_obs,
            "eta_perm_median": null_med,
            "eta_perm_p95": null_p95,
            "p_empirical": p_emp,
            "F": F, "p_F": pval,
            "underpowered": "TRUE" if underpowered else "FALSE",
        })
        print(f"\n  Group: {g} (n={n})")
        print(f"    n_tissues={n_groups}  n_singleton={n_single}  n_multi={n_multi}")
        print(f"    eta_obs={eta_obs:.4f}  perm_median={null_med:.4f}  "
              f"perm_p95={null_p95:.4f}  p_emp={p_emp:.4f}")
        print(f"    F={F if F == F else 'NA'}  p_F={pval if pval == pval else 'NA'}")
        print(f"    underpowered={'YES' if underpowered else 'no'}")
    pd.DataFrame(rows_A).to_csv(OUT_DIR / "option_a_anova.csv", index=False)
    print("\nWrote option_a_anova.csv")

    # ============================================================
    # SENSITIVITY (Option C): residuals ~ hybrid label
    # ============================================================
    print("\n\n=== SENSITIVITY (Option C): residuals ~ hybrid label within each epi/2a group ===")
    rows_C = []
    for g in cov82["epi_2a"].unique():
        sub = cov82[cov82["epi_2a"] == g]
        n = len(sub)
        labels = sub["hybrid_C"].to_numpy()
        values = sub["residual"].to_numpy()
        from collections import Counter
        cnt = Counter(labels)
        n_groups = len(cnt)
        n_single = sum(1 for v in cnt.values() if v == 1)
        n_multi  = sum(1 for v in cnt.values() if v >= 2)
        eta_obs = eta_sq(values, labels) if n_groups >= 2 else np.nan
        if n_groups >= 2 and n >= 3:
            null_etas = perm_null_eta(values, labels, B, rng)
            null_med = float(np.median(null_etas))
            null_p95 = float(np.quantile(null_etas, 0.95))
            p_emp   = float(np.mean(null_etas >= eta_obs)) if not np.isnan(eta_obs) else np.nan
            try:
                F, pval = stats.f_oneway(*[values[labels == lab] for lab in cnt])
            except Exception:
                F, pval = (np.nan, np.nan)
        else:
            null_med = null_p95 = p_emp = np.nan
            F = pval = np.nan
        underpowered = (n < 4) or (n_single > n_groups / 2)
        rows_C.append({
            "epi_2a_group": g,
            "n_cancers": n,
            "n_labels": n_groups,
            "n_singleton_labels": n_single,
            "n_multi_labels": n_multi,
            "eta_obs": eta_obs,
            "eta_perm_median": null_med,
            "eta_perm_p95": null_p95,
            "p_empirical": p_emp,
            "F": F, "p_F": pval,
            "underpowered": "TRUE" if underpowered else "FALSE",
        })
        print(f"\n  Group: {g} (n={n})")
        print(f"    n_labels={n_groups}  n_singleton={n_single}  n_multi={n_multi}")
        print(f"    eta_obs={eta_obs:.4f}  perm_median={null_med:.4f}  "
              f"perm_p95={null_p95:.4f}  p_emp={p_emp:.4f}")
        print(f"    F={F if F == F else 'NA'}  p_F={pval if pval == pval else 'NA'}")
        print(f"    underpowered={'YES' if underpowered else 'no'}")
    pd.DataFrame(rows_C).to_csv(OUT_DIR / "option_c_anova.csv", index=False)
    print("\nWrote option_c_anova.csv")

    # ============================================================
    # Cross-method comparison
    # ============================================================
    print("\n\n=== Option A vs Option C ===")
    A_df = pd.DataFrame(rows_A)
    C_df = pd.DataFrame(rows_C)
    cmp = pd.merge(A_df[["epi_2a_group", "eta_obs", "p_empirical", "underpowered"]]
                       .rename(columns={"eta_obs": "etaA", "p_empirical": "pA",
                                        "underpowered": "underpowered_A"}),
                   C_df[["epi_2a_group", "eta_obs", "p_empirical", "underpowered"]]
                       .rename(columns={"eta_obs": "etaC", "p_empirical": "pC",
                                        "underpowered": "underpowered_C"}),
                   on="epi_2a_group")
    cmp["delta_eta_C_minus_A"] = cmp["etaC"] - cmp["etaA"]
    print(cmp.to_string(index=False))
    cmp.to_csv(OUT_DIR / "option_a_vs_c_comparison.csv", index=False)
    print("\nWrote option_a_vs_c_comparison.csv")


if __name__ == "__main__":
    main()
