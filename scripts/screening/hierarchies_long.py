"""Builds the hierarchies long-format CSV from cancer_type_hierarchies_2026-05-09.md
hand-coded per the master table + canonical-choice notes in the data-handling log.

Output: hierarchies_long.csv with columns:
  cancer, hierarchy, granularity, group, contested

`contested` = True if this (cancer, hierarchy) is flagged in the
"Contested or ambiguous assignments" section of the reference table.
"""
import csv
from pathlib import Path

OUT = Path(r"C:/FkCancer/runs/run_screening_2026-05-09/hierarchies_long.csv")

CANCERS = [
    "ACC", "BLCA", "BRCA", "CESC", "CHOL", "CORE", "DLBC", "ESCA", "HNSC",
    "KICH", "KIRC", "KIRP", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV",
    "PAAD", "PCPG", "PRAD", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM",
    "UCEC", "UCS",
]
assert len(CANCERS) == 29

# -----------------------------------------------------------------------------
# Hierarchy 1 — Tissue of origin
# -----------------------------------------------------------------------------
T1A_5BIN = {  # Granularity 1a (coarsest, 5 bins)
    "Digestive": {"ESCA", "STAD", "CORE", "LIHC", "CHOL", "PAAD"},
    "Urogen+Repro": {"BLCA", "KICH", "KIRC", "KIRP", "PRAD", "TGCT", "BRCA",
                     "CESC", "OV", "UCEC", "UCS"},
    "Endocrine": {"ACC", "PCPG", "THCA", "THYM"},
    "Thoracic": {"LUAD", "LUSC", "MESO"},
    "Other": {"LGG", "DLBC", "SKCM", "SARC", "HNSC"},
}
T1B_11BIN = {  # Granularity 1b (mid, 11 bins)
    "Lung+pleura": {"LUAD", "LUSC", "MESO"},
    "Kidney": {"KICH", "KIRC", "KIRP"},
    "Hepatobiliary+pancreas": {"LIHC", "CHOL", "PAAD"},
    "GI tract": {"ESCA", "STAD", "CORE"},
    "Breast": {"BRCA"},
    "Gynecological": {"OV", "UCEC", "UCS", "CESC"},
    "Male repro+bladder": {"PRAD", "TGCT", "BLCA"},
    "Endocrine": {"ACC", "PCPG", "THCA", "THYM"},
    "CNS": {"LGG"},
    "Skin": {"SKCM"},
    "Soft+heme+H&N": {"SARC", "DLBC", "HNSC"},
}
# Granularity 1c (29 singletons) — DEGENERATE, skip from screening

# -----------------------------------------------------------------------------
# Hierarchy 2 — Epithelial classification
# -----------------------------------------------------------------------------
T2A_3BIN = {  # Granularity 2a (3 bins). MESO and UCS contested per file.
    "Epithelial": {"ACC", "BLCA", "BRCA", "CESC", "CHOL", "CORE", "ESCA",
                   "HNSC", "KICH", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC",
                   "OV", "PAAD", "PRAD", "STAD", "THCA", "THYM", "UCEC"},
    "Non-epithelial": {"LGG", "SARC", "SKCM", "TGCT", "PCPG"},
    "Hematological": {"DLBC"},
    # Contested: MESO, UCS (boundary cases — drop in screening)
}
T2B_6BIN = {  # Granularity 2b (6 bins).
    # ESCA / CESC: file says "split" — primary assignment by Hoadley HPV+squamous
    # for CESC and Pan-GI/Pan-squamous-mixed for ESCA. Place CESC in squamous,
    # ESCA in squamous (HPV-negative ESCC dominant in TCGA per published cohort).
    # NOTE: choice flagged in data_handling_log.md.
    "Squamous": {"HNSC", "LUSC", "CESC", "ESCA"},
    "Adenocarcinoma": {"BRCA", "CORE", "LUAD", "PRAD", "PAAD", "STAD", "UCEC",
                      "OV", "ACC", "THCA", "LIHC", "CHOL"},
    "Transitional": {"BLCA"},
    "Glial": {"LGG"},
    "Mesenchymal/germ-cell": {"SARC", "MESO", "TGCT", "UCS"},
    "Heme/neural-crest/sui-generis": {"DLBC", "SKCM", "PCPG", "KICH", "KIRC",
                                       "KIRP", "THYM"},
}
T2C_10BIN = {  # Granularity 2c (10 bins).
    "Squamous": {"HNSC", "LUSC", "CESC", "ESCA"},
    "Glandular adeno": {"BRCA", "CORE", "LUAD", "PRAD", "PAAD", "STAD", "ACC",
                        "THCA"},
    "HCC": {"LIHC"},
    "Cholangiocarcinoma": {"CHOL"},
    "RCC": {"KIRC", "KIRP", "KICH"},
    "Urothelial": {"BLCA"},
    "Endometrioid+serous": {"UCEC", "OV"},
    "Glial": {"LGG"},
    "Mesench/mesothelial/sarcoma": {"SARC", "MESO", "UCS"},
    "Other (germ/melan/neuroend/heme/thymic)": {"TGCT", "SKCM", "PCPG", "DLBC",
                                                 "THYM"},
}

# -----------------------------------------------------------------------------
# Hierarchy 3 — Embryonic germ layer (most contested)
# -----------------------------------------------------------------------------
T3A_3BIN = {  # Granularity 3a (3 bins). Many contested.
    "Endoderm": {"ESCA", "STAD", "CORE", "LIHC", "CHOL", "PAAD", "LUAD",
                 "LUSC", "THCA"},
    "Mesoderm": {"KICH", "KIRC", "KIRP", "OV", "ACC", "SARC", "DLBC", "MESO",
                 "UCEC", "UCS"},
    "Ectoderm": {"SKCM", "LGG", "PCPG"},
    # Contested per file: BLCA, CESC, HNSC, PRAD, THYM (germ-layer mixed)
    # Plus BRCA (ectoderm with stromal mesoderm), TGCT (PGC pre-germ-layer)
}
T3B_MID = {  # Granularity 3b (sub-divisions).
    "Foregut endoderm": {"ESCA", "STAD", "LIHC", "CHOL", "PAAD", "LUAD",
                         "LUSC", "THCA"},
    "Hindgut endoderm": {"CORE"},
    "Neural ectoderm": {"LGG"},
    "Neural crest": {"PCPG", "SKCM"},
    "Intermediate mesoderm": {"ACC", "KICH", "KIRC", "KIRP", "OV"},
    "Lateral plate mesoderm": {"MESO", "UCEC", "UCS"},
    "Hematopoietic mesoderm": {"DLBC"},
    "Paraxial+lateral mesoderm": {"SARC"},
    # Contested same as 3a; SKCM placed in neural crest per "most authorities"
    # (file note); BRCA dropped (was alone in surface ectoderm with SKCM).
}

# -----------------------------------------------------------------------------
# Hierarchy 4 — Sex composition (1 granularity, ternary or quaternary)
# -----------------------------------------------------------------------------
T4_SEX = {  # 4 bins (Both / Female / Male / Female-dominant).
    "Both": {"ACC", "BLCA", "CHOL", "CORE", "DLBC", "ESCA", "HNSC", "KICH",
             "KIRC", "KIRP", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "PAAD",
             "PCPG", "SARC", "SKCM", "STAD", "THCA", "THYM"},
    "Female-only": {"CESC", "OV", "UCEC", "UCS"},
    "Male-only": {"PRAD", "TGCT"},
    "Female-dominant": {"BRCA"},
}

# -----------------------------------------------------------------------------
# Hierarchy 5 — Hoadley 2018 (CIRCULARITY-FLAGGED)
# -----------------------------------------------------------------------------
T5A_4BIN_PLUS_UNASSIGNED = {  # Granularity 5a (4 patterns + Unassigned).
    "Pan-GI": {"CORE", "STAD", "CHOL", "PAAD"},
    "Pan-gyn": {"BRCA", "OV", "UCEC"},
    "Pan-kidney": {"KIRC", "KIRP"},
    "Pan-squamous": {"HNSC", "LUSC"},
    "Unassigned": {"ACC", "DLBC", "KICH", "LIHC", "LUAD", "MESO", "PCPG",
                   "PRAD", "SARC", "SKCM", "TGCT", "THCA", "UCS"},
    # Dropped contested-for-Hoadley: LGG (split C11/C23), THYM (Table S6 amb.),
    # ESCA (Pan-GI/Pan-sq mixed), CESC (Pan-gyn/Pan-sq split), BLCA (subset).
}
# Granularity 5b — does NOT exist per Hoadley 2018 (no canonical mid-level).

T5C_28BIN = {  # Granularity 5c (primary iCluster assignment).
    # Cancers Hoadley flagged "Diverse" or "<50% in one cluster" → DROP.
    "C3": {"DLBC"},
    "C4": {"CORE"},
    "C5": {"PCPG"},
    "C6": {"OV"},
    "C8": {"UCEC"},
    "C9": {"ACC", "KICH"},
    "C12": {"THCA"},
    "C14": {"LUAD"},
    "C15": {"SKCM"},
    "C16": {"PRAD"},
    "C19": {"BRCA"},  # luminal subset; HER2-amp goes to C2 — flagged
    "C20": {"SARC"},
    "C26": {"LIHC"},
    "C28": {"KIRC", "KIRP"},
    # Dropped: BLCA, CHOL, ESCA, HNSC, MESO, PAAD, STAD, TGCT, UCS, THYM, LGG
    # (Hoadley "<50% in any" or contested per file).
}

# -----------------------------------------------------------------------------
# Contested cancers per hierarchy (per "Contested or ambiguous assignments")
# -----------------------------------------------------------------------------
CONTESTED = {
    "tissue_origin": set(),  # No contested at tissue level per file
    "epithelial_class": {"MESO", "UCS"},  # boundary/biphasic
    "germ_layer": {"BLCA", "BRCA", "CESC", "HNSC", "PRAD", "TGCT", "THYM"},
    "sex_composition": set(),
    "hoadley": {"LGG", "THYM", "ESCA", "CESC", "BLCA"},  # see 5a notes
}
# Granularity-specific contested for epithelial 2c (KIRC/KIRP/KICH RCC at finest)
CONTESTED_BY_GRAN = {
    ("epithelial_class", "2c_10bin"): {"MESO", "UCS"},  # same as 2a/2b
}


def emit(cancer, hierarchy, granularity, group, contested):
    return {
        "cancer": cancer,
        "hierarchy": hierarchy,
        "granularity": granularity,
        "group": group,
        "contested": "TRUE" if contested else "FALSE",
    }


def expand(mapping, hierarchy, granularity, contested_set):
    rows = []
    seen = set()
    for group, members in mapping.items():
        for c in members:
            assert c in CANCERS, f"unknown cancer: {c}"
            assert c not in seen, f"duplicate {c} in {hierarchy} {granularity}"
            seen.add(c)
            rows.append(emit(c, hierarchy, granularity, group,
                             c in contested_set))
    # Mark cancers not in mapping but in CANCERS as 'NA' (not assigned at this
    # granularity) — they're effectively dropped for screening anyway.
    return rows


def main():
    rows = []

    # Hierarchy 1 — tissue
    rows += expand(T1A_5BIN, "tissue_origin", "1a_5bin",
                   CONTESTED["tissue_origin"])
    rows += expand(T1B_11BIN, "tissue_origin", "1b_11bin",
                   CONTESTED["tissue_origin"])
    # 1c (29 singletons) skipped (degenerate)

    # Hierarchy 2 — epithelial
    rows += expand(T2A_3BIN, "epithelial_class", "2a_3bin",
                   CONTESTED["epithelial_class"])
    rows += expand(T2B_6BIN, "epithelial_class", "2b_6bin",
                   CONTESTED["epithelial_class"])
    rows += expand(T2C_10BIN, "epithelial_class", "2c_10bin",
                   CONTESTED["epithelial_class"])

    # Hierarchy 3 — germ layer
    rows += expand(T3A_3BIN, "germ_layer", "3a_3bin",
                   CONTESTED["germ_layer"])
    rows += expand(T3B_MID, "germ_layer", "3b_mid",
                   CONTESTED["germ_layer"])

    # Hierarchy 4 — sex
    rows += expand(T4_SEX, "sex_composition", "4_4bin",
                   CONTESTED["sex_composition"])

    # Hierarchy 5 — Hoadley
    rows += expand(T5A_4BIN_PLUS_UNASSIGNED, "hoadley", "5a_4plus_unassigned",
                   CONTESTED["hoadley"])
    rows += expand(T5C_28BIN, "hoadley", "5c_iclusters",
                   CONTESTED["hoadley"])

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(
            f, fieldnames=["cancer", "hierarchy", "granularity", "group",
                           "contested"])
        w.writeheader()
        w.writerows(rows)
    print(f"Wrote {OUT} ({len(rows)} rows)")
    # Sanity: per (hierarchy, granularity), how many cancers assigned?
    from collections import Counter
    cnt = Counter((r["hierarchy"], r["granularity"]) for r in rows)
    for k, v in cnt.items():
        contested_count = sum(
            1 for r in rows if (r["hierarchy"], r["granularity"]) == k
            and r["contested"] == "TRUE")
        print(f"  {k[0]:20s} {k[1]:25s}  n={v:3d}  contested={contested_count}")


if __name__ == "__main__":
    main()
