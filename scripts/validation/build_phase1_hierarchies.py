"""Build per-cancer-type hierarchy assignments for Phase 1 (Option B').

Inputs (already on disk):
  C:/FkCancer/data/PanTCGA/Thorsson_mmc2.xlsx      sheet 'PanImmune_MS'
    -> per-tumor "Immune Subtype" C1-C6 + TCGA Study cancer label
  C:/FkCancer/data/PanTCGA/Malta_mmc1.xlsx          sheet 'StemnessScores_DNAmeth'
    -> per-tumor mDNAsi + cancer.type
  C:/FkCancer/data/PanTCGA/SanchezVega_TableS4.xlsx sheet 'Pathway level'
    -> per-tumor 10 binary pathway alteration flags + SAMPLE_BARCODE (cancer derived from barcode prefix lookup)

Outputs:
  C:/Cancer Research/reference/phase1_hierarchies_2026-05-15.md
  C:/Cancer Research/reference/phase1_hierarchy_assignments.csv

The 29 BIDIFAC+ cancers (CORE = COAD + READ merged):
"""
import pathlib, openpyxl, pandas as pd
from collections import Counter

BIDIFAC_29 = ["ACC","BLCA","BRCA","CESC","CHOL","CORE","DLBC","ESCA","HNSC","KICH",
              "KIRC","KIRP","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
              "PRAD","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS"]

DATA = pathlib.Path(r"C:/FkCancer/data/PanTCGA")
OUT  = pathlib.Path(r"C:/Cancer Research/reference")
OUT.mkdir(exist_ok=True)

# -----------------------------------------------------------
# 1) Thorsson immune subtypes (C1-C6)
# -----------------------------------------------------------
print("[1] Thorsson immune subtypes ...")
wb = openpyxl.load_workbook(DATA/"Thorsson_mmc2.xlsx", read_only=True, data_only=True)
ws = wb["PanImmune_MS"]
rows = ws.iter_rows(values_only=True)
hdr = next(rows)
idx_barcode = hdr.index("TCGA Participant Barcode")
idx_study   = hdr.index("TCGA Study")
idx_imm     = hdr.index("Immune Subtype")

per_cancer_immune = {c: Counter() for c in BIDIFAC_29}
n_total = 0
for r in rows:
    study = r[idx_study]; imm = r[idx_imm]
    if study is None or imm is None or imm == "NA": continue
    # CORE = COAD + READ
    if study in ("COAD","READ"): study = "CORE"
    if study not in per_cancer_immune: continue
    per_cancer_immune[study][str(imm)] += 1
    n_total += 1
wb.close()
print(f"  read {n_total} tumors across {len([c for c in per_cancer_immune if sum(per_cancer_immune[c].values())])} cancers")

# Dominant subtype per cancer
immune_assignment = {}
for c in BIDIFAC_29:
    counter = per_cancer_immune[c]
    if not counter:
        immune_assignment[c] = None
        continue
    dominant, n = counter.most_common(1)[0]
    total = sum(counter.values())
    dom_str = str(dominant)
    if not dom_str.startswith("C"):
        dom_str = f"C{dom_str}"
    immune_assignment[c] = dom_str
    print(f"  {c:<6} n={total:<4} dominant={dom_str} ({n}/{total} = {100*n/total:.0f}%)   full: {dict(counter)}")

# -----------------------------------------------------------
# 2) Malta stemness — tertile split on per-cancer mean mDNAsi
# -----------------------------------------------------------
print("\n[2] Malta stemness ...")
wb = openpyxl.load_workbook(DATA/"Malta_mmc1.xlsx", read_only=True, data_only=True)
ws = wb["StemnessScores_DNAmeth"]
rows = ws.iter_rows(values_only=True)
hdr = next(rows)
idx_id = hdr.index("TCGAlong.id")
idx_ct = hdr.index("cancer.type")
idx_st = hdr.index("sample.type")  # 1 = primary tumor; we keep all primary
idx_dna = hdr.index("mDNAsi")

per_cancer_mdnasi = {c: [] for c in BIDIFAC_29}
for r in rows:
    ct = r[idx_ct]; st = r[idx_st]; dnasi = r[idx_dna]
    if ct is None or dnasi is None: continue
    if str(st) != "1": continue  # primary tumor only
    if ct in ("COAD","READ"): ct = "CORE"
    if ct not in per_cancer_mdnasi: continue
    try:
        per_cancer_mdnasi[ct].append(float(dnasi))
    except (TypeError, ValueError):
        continue
wb.close()

mean_mdnasi = {c: (sum(v)/len(v) if v else None) for c, v in per_cancer_mdnasi.items()}
print("  per-cancer mean mDNAsi:")
for c in BIDIFAC_29:
    n = len(per_cancer_mdnasi[c])
    mn = mean_mdnasi[c]
    print(f"    {c:<6} n={n:<5} mean_mDNAsi={mn:.4f}" if mn is not None else f"    {c:<6} n=0     NO DATA")

# Tertile-split
present = [(c, mean_mdnasi[c]) for c in BIDIFAC_29 if mean_mdnasi[c] is not None]
present.sort(key=lambda x: x[1])
n = len(present)
t1, t2 = n // 3, 2 * n // 3
stemness_assignment = {}
for i, (c, _) in enumerate(present):
    if i < t1:   stemness_assignment[c] = "low"
    elif i < t2: stemness_assignment[c] = "med"
    else:        stemness_assignment[c] = "high"
for c in BIDIFAC_29:
    if c not in stemness_assignment:
        stemness_assignment[c] = None
print(f"  tertile cutoffs: low<={present[t1-1][1]:.3f}, med<={present[t2-1][1]:.3f}, high>{present[t2-1][1]:.3f}")

# -----------------------------------------------------------
# 3) Sanchez-Vega oncogenic pathways - dominant pathway per cancer
# -----------------------------------------------------------
print("\n[3] Sanchez-Vega oncogenic pathways ...")
# Need cancer type per sample barcode. Use the Sanchez-Vega "Alteration level" sheet
# header to inspect, or use known TCGA barcode -> cancer mapping. Easier: derive cancer
# from the SAMPLE_BARCODE prefix... actually TCGA barcodes don't encode cancer.
# Use the Thorsson file's barcode -> study map we already loaded.

# Re-load Thorsson barcode -> study map
print("  re-loading Thorsson barcode -> cancer map ...")
wb = openpyxl.load_workbook(DATA/"Thorsson_mmc2.xlsx", read_only=True, data_only=True)
ws = wb["PanImmune_MS"]
rows = ws.iter_rows(values_only=True)
hdr = next(rows)
idx_barcode = hdr.index("TCGA Participant Barcode")
idx_study = hdr.index("TCGA Study")
barcode_to_cancer = {}
for r in rows:
    if r[idx_barcode] is None or r[idx_study] is None: continue
    bc = r[idx_barcode]  # TCGA-XX-XXXX (12 chars, no -01 suffix)
    study = r[idx_study]
    if study in ("COAD","READ"): study = "CORE"
    barcode_to_cancer[bc] = study
wb.close()
print(f"  {len(barcode_to_cancer)} barcode->cancer entries")

# Load Sanchez-Vega pathway alterations
wb = openpyxl.load_workbook(DATA/"SanchezVega_TableS4.xlsx", read_only=True, data_only=True)
ws = wb["Pathway level"]
rows = ws.iter_rows(values_only=True)
hdr = list(next(rows))
print(f"  Pathway level header: {hdr}")
pathway_cols = hdr[1:]   # everything except SAMPLE_BARCODE

per_cancer_pathway_alts = {c: {p: [0, 0] for p in pathway_cols} for c in BIDIFAC_29}
# value [0]=count_altered, [1]=total_samples
n_unmapped = 0
for r in rows:
    bc15 = r[0]   # e.g., TCGA-OR-A5J1-01
    if bc15 is None: continue
    bc12 = bc15[:12]
    cancer = barcode_to_cancer.get(bc12)
    if cancer is None or cancer not in per_cancer_pathway_alts:
        n_unmapped += 1; continue
    for j, p in enumerate(pathway_cols):
        v = r[j+1]
        if v is None or v == "NA": continue
        try:
            vi = int(v)
        except (TypeError, ValueError):
            continue
        per_cancer_pathway_alts[cancer][p][1] += 1
        if vi == 1:
            per_cancer_pathway_alts[cancer][p][0] += 1
wb.close()
print(f"  unmapped barcodes: {n_unmapped}")

# Dominant pathway per cancer = highest alteration rate
pathway_assignment = {}
print("  per-cancer top pathway by alteration rate:")
for c in BIDIFAC_29:
    pathways = per_cancer_pathway_alts[c]
    rates = []
    for p, (alt, tot) in pathways.items():
        if tot > 0: rates.append((p, alt/tot, alt, tot))
    if not rates:
        pathway_assignment[c] = None
        print(f"    {c:<6} NO DATA")
        continue
    rates.sort(key=lambda x: -x[1])
    top_p, top_rate, alt, tot = rates[0]
    pathway_assignment[c] = top_p
    print(f"    {c:<6} top={top_p:<12} alt={alt}/{tot} ({100*top_rate:.0f}%)")

# -----------------------------------------------------------
# Write combined assignment table
# -----------------------------------------------------------
rows = []
for c in BIDIFAC_29:
    rows.append({
        "cancer": c,
        "thorsson_immune": immune_assignment[c],
        "malta_stemness_tertile": stemness_assignment[c],
        "sanchez_vega_dominant_pathway": pathway_assignment[c],
    })
df = pd.DataFrame(rows)
out_csv = OUT / "phase1_hierarchy_assignments.csv"
df.to_csv(out_csv, index=False)
print(f"\nWrote {out_csv}")
print(df.to_string(index=False))
