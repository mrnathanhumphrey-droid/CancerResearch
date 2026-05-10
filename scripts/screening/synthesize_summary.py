"""Synthesize SCREENING_SUMMARY.md from screening_results.csv.

Sections (per user brief):
  1. Top-20 (cov, hier, gran) triples by eta_squared, split Hoadley vs non-Hoadley
  2. eta_squared distribution: non-Hoadley vs Hoadley
  3. Broadly-informative hierarchies (>=1 granularity producing eta>0.30 on >=5 covs)
  4. Structurally-explainable covariates (>=1 (hier, gran) producing eta>0.30)
  5. Hoadley vs best-non-Hoadley eta gap per covariate (both > 0.20)

Interpretation framing:
  - Candidate identification only; validation requires reduced-parameter Gibbs run.
  - Flag Hoadley dominance (median eta gap > 0.15) as suspicious -> TCGA-specific.
  - Granularity-robustness check.
"""
import numpy as np
import pandas as pd
from pathlib import Path

IN_CSV  = Path(r"C:/FkCancer/runs/run_screening_2026-05-09/screening_results.csv")
OUT_MD  = Path(r"C:/FkCancer/runs/run_screening_2026-05-09/SCREENING_SUMMARY.md")
ETA_THR = 0.30
ETA_MIN = 0.20

# Use non-degenerate results only for headline claims; full results retained
# in the CSV.

def load():
    df = pd.read_csv(IN_CSV)
    df["F"] = pd.to_numeric(df["F"], errors="coerce")
    df["p_value"] = pd.to_numeric(df["p_value"], errors="coerce")
    df["eta_squared"] = pd.to_numeric(df["eta_squared"], errors="coerce")
    # pandas auto-sniffs "TRUE"/"FALSE" → bool, so coerce robustly
    if df["hoadley_flag"].dtype != bool:
        df["hoadley_flag"] = df["hoadley_flag"].astype(str).str.upper() == "TRUE"
    if df["degenerate_majority_singletons"].dtype != bool:
        df["degenerate_majority_singletons"] = (
            df["degenerate_majority_singletons"].astype(str).str.upper() == "TRUE")
    df["degenerate"] = df["degenerate_majority_singletons"]
    # cov_id may have been auto-converted to float; preserve as string label
    df["cov_id"] = df["cov_id"].astype(str)
    return df


def section_1_top20(df, ndg):
    out = ["## 1. Top-20 (covariate × hierarchy × granularity) triples by η²\n",
           "All valid rows (degenerate ones flagged with † in the deg column). "
           "Per the brief, degenerate rows are not dropped — they still "
           "contribute to between-group variance.\n"]
    valid = df[df["eta_squared"].notna()].copy()
    out.append("\n### Non-Hoadley (external biological hierarchies)\n")
    nh = valid[~valid["hoadley_flag"]].nlargest(20, "eta_squared")
    out.append(_format_table(nh))
    out.append("\n### Hoadley (circularity-flagged ⚠)\n")
    h = valid[valid["hoadley_flag"]].nlargest(20, "eta_squared")
    out.append(_format_table(h))
    return "\n".join(out)


def _format_table(rows):
    lines = ["| cov_id | hierarchy | granularity | n_groups | n_cancers | "
             "n_single | η² | F | p | deg |",
             "|---|---|---|---|---|---|---|---|---|---|"]
    for _, r in rows.iterrows():
        f_str = f"{r['F']:.2f}" if pd.notna(r['F']) else "—"
        p_str = f"{r['p_value']:.3g}" if pd.notna(r['p_value']) else "—"
        deg = "†" if bool(r['degenerate']) else ""
        lines.append(
            f"| {r['cov_id']} | {r['hierarchy']} | {r['granularity']} | "
            f"{int(r['n_groups'])} | {int(r['n_cancers_used'])} | "
            f"{int(r['n_singleton_groups'])} | "
            f"{r['eta_squared']:.3f} | {f_str} | {p_str} | {deg} |"
        )
    return "\n".join(lines) + "\n"


def section_2_distribution(df, ndg):
    out = ["\n## 2. η² distribution: Hoadley vs non-Hoadley\n",
           "Reported on (a) all valid rows and (b) non-degenerate rows.\n"]
    valid = df[df["eta_squared"].notna()]
    for label, sub in [("ALL valid rows", valid),
                       ("Non-degenerate only", ndg)]:
        out.append(f"\n### {label}\n")
        nh = sub[~sub["hoadley_flag"]]["eta_squared"].dropna()
        h  = sub[sub["hoadley_flag"]]["eta_squared"].dropna()
        rows = []
        for q in [0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0]:
            rows.append((f"{int(q*100)}th pct",
                         nh.quantile(q) if len(nh) else float("nan"),
                         h.quantile(q) if len(h) else float("nan")))
        rows.append(("mean", nh.mean(), h.mean()))
        rows.append(("count", len(nh), len(h)))
        out.append("| stat | non-Hoadley | Hoadley | gap |\n|---|---|---|---|")
        for lab, a, b in rows:
            if isinstance(a, (int, np.integer)):
                out.append(f"| {lab} | {a} | {b} | — |")
            elif pd.notna(a) and pd.notna(b):
                out.append(f"| {lab} | {a:.3f} | {b:.3f} | {b - a:+.3f} |")
            else:
                a_s = f"{a:.3f}" if pd.notna(a) else "—"
                b_s = f"{b:.3f}" if pd.notna(b) else "—"
                out.append(f"| {lab} | {a_s} | {b_s} | — |")
        if len(nh) and len(h) and pd.notna(nh.median()) and pd.notna(h.median()):
            median_gap = h.median() - nh.median()
            flag = "⚠ FLAG" if median_gap > 0.15 else "OK"
            out.append(f"\n**Median η² gap (Hoadley − non-Hoadley) = "
                       f"{median_gap:+.3f}** ({flag} per user's >0.15 "
                       f"suspicion threshold)\n")
    return "\n".join(out)


def section_3_broadly_informative(df, ndg):
    """For each hierarchy: list granularities where eta>0.30 hits >=5 covs."""
    out = ["\n## 3. Broadly-informative hierarchies\n",
           "Definition: at least one granularity producing η² > 0.30 on at "
           "least 5 covariates (non-degenerate rows only).\n"]
    out.append("| hierarchy | granularity | n_covs with η²>0.30 | broadly_informative |")
    out.append("|---|---|---|---|")
    summary = []
    for (h, g), sub in ndg.groupby(["hierarchy", "granularity"]):
        n_hits = int((sub["eta_squared"] > ETA_THR).sum())
        bi = "✓" if n_hits >= 5 else ""
        summary.append((h, g, n_hits, bi))
    summary.sort(key=lambda r: (-r[2], r[0], r[1]))
    for h, g, n, bi in summary:
        out.append(f"| {h} | {g} | {n} | {bi} |")
    out.append("")
    informative_hiers = sorted({h for h, g, n, bi in summary if bi == "✓"})
    out.append(f"**Broadly-informative hierarchies**: "
               f"{', '.join(informative_hiers) if informative_hiers else 'none'}\n")
    return "\n".join(out)


def section_4_structurally_explainable(df, ndg):
    out = ["\n## 4. Structurally-explainable covariates\n",
           "Definition: at least one (hierarchy, granularity) producing "
           "η² > 0.30 (non-degenerate, non-Hoadley).\n"]
    nh = ndg[~ndg["hoadley_flag"]]
    by_cov = nh[nh["eta_squared"] > ETA_THR].groupby("cov_id").size().sort_values(ascending=False)
    out.append(f"**{len(by_cov)} covariates** have at least one non-Hoadley "
               f"(hierarchy, granularity) hit at η² > 0.30.\n")
    out.append("| cov_id | n_hits | best (hier, gran) | best η² |\n|---|---|---|---|")
    for cov, n_hits in by_cov.items():
        sub = nh[(nh["cov_id"] == cov) & (nh["eta_squared"] > ETA_THR)]
        best = sub.nlargest(1, "eta_squared").iloc[0]
        out.append(
            f"| {cov} | {n_hits} | "
            f"{best['hierarchy']} / {best['granularity']} | "
            f"{best['eta_squared']:.3f} |")
    out.append("")
    return "\n".join(out)


def section_5_hoadley_vs_best_nonhoadley(df, ndg):
    out = ["\n## 5. Hoadley vs best non-Hoadley η² per covariate\n",
           "For each covariate where BOTH Hoadley and non-Hoadley produce "
           "η² > 0.20 (using ALL valid rows). Median gap > 0 = Hoadley wins; "
           "Hoadley dominance suggests TCGA-circularity rather than external "
           "biology.\n"]
    valid = df[df["eta_squared"].notna()]
    nh = valid[~valid["hoadley_flag"]]
    h  = valid[valid["hoadley_flag"]]
    rows = []
    for cov in sorted(valid["cov_id"].unique(), key=lambda x: float(x)):
        nh_sub = nh[(nh["cov_id"] == cov) & (nh["eta_squared"] > ETA_MIN)]
        h_sub  = h[(h["cov_id"] == cov) & (h["eta_squared"] > ETA_MIN)]
        if len(nh_sub) == 0 or len(h_sub) == 0:
            continue
        nh_best = nh_sub["eta_squared"].max()
        nh_best_row = nh_sub.loc[nh_sub["eta_squared"].idxmax()]
        h_best  = h_sub["eta_squared"].max()
        gap = h_best - nh_best
        rows.append({
            "cov_id": cov,
            "best_nonhoadley_hier_gran": f"{nh_best_row['hierarchy']}/{nh_best_row['granularity']}",
            "nh_eta": nh_best,
            "h_eta": h_best,
            "gap": gap,
        })
    out.append("| cov_id | best non-Hoadley | nh η² | Hoadley η² | gap |")
    out.append("|---|---|---|---|---|")
    rows.sort(key=lambda r: -r["gap"])
    for r in rows:
        out.append(f"| {r['cov_id']} | {r['best_nonhoadley_hier_gran']} | "
                   f"{r['nh_eta']:.3f} | {r['h_eta']:.3f} | "
                   f"{r['gap']:+.3f} |")
    out.append("")
    if rows:
        gaps = [r["gap"] for r in rows]
        out.append(f"**Median gap across {len(rows)} covariates: "
                   f"{np.median(gaps):+.3f}** "
                   f"(positive = Hoadley dominates; negative = external "
                   f"hierarchy outperforms).\n")
    return "\n".join(out)


def section_6_granularity_robustness(df, ndg):
    """For each hierarchy with multiple granularities, check whether the
    ranking of "most-explainable covariates" is consistent across granularities.
    Spearman rho between (cov -> max eta²) at each granularity."""
    out = ["\n## 6. Granularity robustness\n",
           "Spearman ρ between covariate η² rankings at different granularities "
           "of the same hierarchy. High ρ = same covariates land on top "
           "regardless of bin choice (real structure). Low ρ = granularity is "
           "load-bearing (k=4 vs k=6 patterns).\n"]
    nh = ndg[~ndg["hoadley_flag"]]
    out.append("| hierarchy | (g1, g2) | Spearman ρ | n_covs_compared |")
    out.append("|---|---|---|---|")
    for hier, sub in nh.groupby("hierarchy"):
        gs = sorted(sub["granularity"].unique())
        for i, g1 in enumerate(gs):
            for g2 in gs[i+1:]:
                a = sub[sub["granularity"] == g1].set_index("cov_id")["eta_squared"]
                b = sub[sub["granularity"] == g2].set_index("cov_id")["eta_squared"]
                joint = pd.concat([a, b], axis=1, join="inner").dropna()
                if len(joint) < 5:
                    continue
                rho, _ = (joint.corr(method="spearman").iloc[0, 1], None)
                out.append(f"| {hier} | ({g1}, {g2}) | {rho:+.3f} | {len(joint)} |")
    out.append("")
    return "\n".join(out)


def main():
    df = load()
    ndg = df[~df["degenerate"] & df["eta_squared"].notna()].copy()
    print(f"Total rows: {len(df)}; non-degenerate: {len(ndg)}")
    parts = []
    parts.append("# Structural-Hierarchy Screening — SUMMARY\n")
    parts.append("**Run date**: 2026-05-09\n")
    parts.append("**Source PIPs**: "
                 "`runs/run_paper_a_gibbs_2026-05-09/per_cancer_beta_pip.csv` "
                 "(551 rows, 29 cancers × variable covariate sets, 68 unique "
                 "covariate IDs).\n")
    parts.append("**Hierarchies tested**: 5 (tissue_origin, epithelial_class, "
                 "germ_layer, sex_composition, hoadley) × 10 granularities. "
                 "Granularity 1c (29 singleton bins) excluded as degenerate.\n")
    parts.append("**Tests run**: " + str(len(df)) +
                 " (cov × hierarchy × granularity), "
                 + str(len(ndg)) + " non-degenerate after the "
                 "n_singleton_groups > n_groups/2 filter.\n")
    parts.append("**Statistical metric**: η² = SS_between / SS_total per "
                 "one-way ANOVA. F and p reported but η² is the primary "
                 "screening metric (n=29 with mostly small group counts).\n")

    parts.append("\n---\n")
    parts.append("## INTERPRETATION FRAMING (READ FIRST)\n")
    parts.append(
        "This is **screening only**, not validation. Triples flagged here are "
        "**candidates** for validation; validation requires a reduced-parameter "
        "Gibbs replacement using structural priors at the chosen "
        "(hierarchy, granularity), held-out predictive performance vs the "
        "baseline replication. That is a separate downstream phase.\n\n"
        "**Do not** read these results as 'the closed-form decomposition is "
        "valid.' They identify the (hierarchy, granularity, covariate) cells "
        "with sufficient between-group variance to be **worth** validating.\n"
    )
    parts.append("---\n")

    parts.append(section_1_top20(df, ndg))
    parts.append(section_2_distribution(df, ndg))
    parts.append(section_3_broadly_informative(df, ndg))
    parts.append(section_4_structurally_explainable(df, ndg))
    parts.append(section_5_hoadley_vs_best_nonhoadley(df, ndg))
    parts.append(section_6_granularity_robustness(df, ndg))

    parts.append("\n---\n## Files\n")
    parts.append("- `screening_results.csv` — long-format results (every "
                 "(cov, hier, gran) test, including degenerate rows)")
    parts.append("- `hierarchies_long.csv` — cancer × hierarchy × granularity "
                 "× group (with contested flag)")
    parts.append("- `data_handling_log.md` — design choices: contested-cell "
                 "drops, canonical group placements for split histologies, "
                 "Hoadley iCluster primary-assignment choices")
    parts.append("- `run_screening.py` — ANOVA pipeline; "
                 "`hierarchies_long.py` — hierarchy CSV builder; "
                 "`synthesize_summary.py` — this report.\n")

    OUT_MD.write_text("\n".join(parts), encoding="utf-8")
    print(f"Wrote {OUT_MD} ({len(parts)} sections)")


if __name__ == "__main__":
    main()
