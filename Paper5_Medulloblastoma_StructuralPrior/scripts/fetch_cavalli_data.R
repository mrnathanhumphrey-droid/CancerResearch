# Paper 5 — fetch Cavalli 2017 medulloblastoma data from GEO + supplementary.
#
# Downloads GSE85217 (expression) + GSE85212 (methylation) via GEOquery and
# reads clinical + 12-subtype labels + survival from Cavalli supplementary
# Table S1 (Cell Press mmc2.xlsx, already placed at data/cavalli_supplementary_S1.xlsx).
#
# Supplementary mmc files placed in data/ at pre-reg time (2026-05-15):
#   data/cavalli_supplementary_S1.xlsx (Table S1, clinical + subgroup + 12-subtype + OS)
#   data/cavalli_supplementary_S2_subgroup_specific.xlsx
#   data/cavalli_supplementary_S3_topgenes.xlsx
#   data/cavalli_supplementary_S4_subtype_specific.xlsx
# Source: Elsevier CDN ars.els-cdn.com pattern for PII S1535610817302015 (Cavalli 2017).
#
# Outputs:
#   data/cavalli_expr.rds          — Affymetrix Hu Gene 1.1 ST expression (probes × 763)
#   data/cavalli_meth.rds          — Illumina 450k methylation β-values (CpG × 763)
#   data/cavalli_expr_pdata.rds    — GEO sample metadata for expression (763 rows)
#   data/cavalli_meth_pdata.rds    — GEO sample metadata for methylation (763 rows)
#   data/cavalli_clinical.rds      — clinical + subgroup (L1) + subtype (L2) + survival (763 rows)
#   data/cavalli_supp_sha256.txt   — SHA-256 of mmc2.xlsx at fetch time

user_lib <- file.path(Sys.getenv("LOCALAPPDATA"), "R/win-library/4.6")
.libPaths(c(user_lib, .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"))

repo_root <- "C:/Cancer Research/Paper5_Medulloblastoma_StructuralPrior"
data_dir <- file.path(repo_root, "data")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

# --- Install GEOquery + readxl on first run (idempotent) -----------------
required <- c("GEOquery", "readxl", "digest", "BiocManager")
for (pkg in required) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg == "GEOquery") {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", lib = user_lib)
      BiocManager::install("GEOquery", lib = user_lib, ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg, lib = user_lib)
    }
  }
}
library(GEOquery)
library(readxl)
library(digest)

# --- 1. Expression: GSE85217 (Affy Hu Gene 1.1 ST, n=763) ----------------
cat(sprintf("[%s] Fetching GSE85217 (Cavalli expression)\n", format(Sys.time())))
gse_expr <- getGEO("GSE85217", GSEMatrix = TRUE, destdir = data_dir, AnnotGPL = FALSE)
gse_expr <- gse_expr[[1]]
expr_mat <- exprs(gse_expr)
expr_pdata <- pData(gse_expr)
cat(sprintf("  expression matrix: %d probes × %d samples\n",
            nrow(expr_mat), ncol(expr_mat)))
stopifnot(ncol(expr_mat) == 763)
saveRDS(expr_mat,   file.path(data_dir, "cavalli_expr.rds"))
saveRDS(expr_pdata, file.path(data_dir, "cavalli_expr_pdata.rds"))

# --- 2. Methylation: GSE85212 supplementary β-values (1.9 GB) -----------
# The GSE85212 series_matrix.txt.gz is metadata-only (36 KB) — actual β-values
# are in the supplementary file GSE85212_Methylation_763samples_*.txt.gz (1.9 GB).
# We pull pData from the series_matrix and the matrix itself from supplementary.
cat(sprintf("[%s] Fetching GSE85212 pData (series_matrix)\n", format(Sys.time())))
gse_meth <- getGEO("GSE85212", GSEMatrix = TRUE, destdir = data_dir, AnnotGPL = FALSE)
gse_meth <- gse_meth[[1]]
meth_pdata <- pData(gse_meth)
saveRDS(meth_pdata, file.path(data_dir, "cavalli_meth_pdata.rds"))

supp_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE85nnn/GSE85212/suppl/GSE85212_Methylation_763samples_SubtypeStudy_TaylorLab_beta_values.txt.gz"
supp_meth_gz <- file.path(data_dir, "GSE85212_meth_beta.txt.gz")
if (!file.exists(supp_meth_gz) || file.info(supp_meth_gz)$size < 1e9) {
  cat(sprintf("[%s] Downloading methylation β-values from GEO FTP (1.9 GB)\n",
              format(Sys.time())))
  options(timeout = 3600)
  download.file(supp_url, destfile = supp_meth_gz, mode = "wb", quiet = FALSE)
} else {
  cat(sprintf("  methylation supplementary file present at %s (skipping download)\n",
              supp_meth_gz))
}

if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages("data.table", lib = user_lib)
}
library(data.table)

cat(sprintf("[%s] Parsing methylation β-values (may take several minutes)\n",
            format(Sys.time())))
# Read the raw first line directly to preserve sample-ID names exactly as
# they appear in the file (fread's check.names / as.matrix can mangle
# underscores). The file's row 1 carries 763 sample IDs (no separate CpG-ID
# header in column 1), so the column-1 cell of row 2 is the first CpG.
con <- gzfile(supp_meth_gz, "r")
header_line <- readLines(con, n = 1L)
close(con)
sample_ids_raw <- strsplit(header_line, "\t")[[1]]
cat(sprintf("  header line has %d fields; first 3 = %s\n",
            length(sample_ids_raw), paste(head(sample_ids_raw, 3), collapse = " ")))
# If the header line has 763 sample IDs (no leading ID-column header),
# then data rows 2+ have CpG-ID in column 1 + 763 β-values.
# fread auto-handles the mismatch by assigning V1 to the first column.
meth_dt <- fread(file = supp_meth_gz, sep = "\t", header = TRUE,
                 na.strings = c("NA", "", "?"),
                 check.names = FALSE)
cat(sprintf("  fread result: %d rows × %d cols\n", nrow(meth_dt), ncol(meth_dt)))

# Decide column-name source: prefer the raw header line over fread's parsing.
cpg_ids <- meth_dt[[1]]
meth_mat <- as.matrix(meth_dt[, -1, with = FALSE])
rownames(meth_mat) <- cpg_ids
if (length(sample_ids_raw) == ncol(meth_mat)) {
  colnames(meth_mat) <- sample_ids_raw  # exact preservation from file
} else if (length(sample_ids_raw) == ncol(meth_mat) + 1) {
  # leading column had a header (e.g., "ID_REF"); drop it
  colnames(meth_mat) <- sample_ids_raw[-1]
} else {
  warning(sprintf("Header line has %d fields but matrix has %d columns; using fread column names",
                  length(sample_ids_raw), ncol(meth_mat)))
}
cat(sprintf("  methylation matrix: %d CpGs × %d samples\n",
            nrow(meth_mat), ncol(meth_mat)))
stopifnot(ncol(meth_mat) == 763)

cn <- colnames(meth_mat)
n_mb <- sum(grepl("^MB_SubtypeStudy_", cn))
cat(sprintf("  column names matching MB_SubtypeStudy_*: %d / %d\n", n_mb, length(cn)))
stopifnot(n_mb >= 700)

saveRDS(meth_mat, file.path(data_dir, "cavalli_meth.rds"))
rm(meth_dt, meth_mat); invisible(gc())

# --- 3. Clinical + subtype from mmc2.xlsx (Cavalli Table S1) -------------
# Verified format (Cell Press mmc2.xlsx fetched 2026-05-15):
#   Sheet name: "Sheet1"
#   Row 1: title text "Table S1, related to Figure 1. ..."
#   Row 2: column headers
#   Rows 3-765: 763 sample rows with Study_ID = MB_SubtypeStudy_NNNNN
#   Header columns of interest:
#     col 0: Study_ID            → sample_id (matches GEO Sample_title)
#     col 1: Age                 → age_at_dx
#     col 3: Gender              → sex
#     col 4: histology           → histology
#     col 5: Met status (1 Met, 0 M0)  → metastatic
#     col 6: Dead                → os_event (1=death, 0=censored)
#     col 7: OS (years)          → os_time_years
#     col 8: Subgroup            → subgroup (L1: WNT/SHH/Group3/Group4)
#     col 9: Subtype             → subtype (L2: WNT_alpha/beta, SHH_alpha-delta, Group3_alpha-gamma, Group4_alpha-gamma)

supp_path <- file.path(data_dir, "cavalli_supplementary_S1.xlsx")
stopifnot(file.exists(supp_path))
supp_sha <- digest(file = supp_path, algo = "sha256")
writeLines(sprintf("file: cavalli_supplementary_S1.xlsx\nsha256: %s\nfetched: %s",
                   supp_sha, format(Sys.time())),
           file.path(data_dir, "cavalli_supp_sha256.txt"))

# Read mmc2 starting at row 2 (skip title row) — first row of read_excel becomes header
raw <- read_excel(supp_path, sheet = "Sheet1", skip = 1, na = c("NA", ""))
cat(sprintf("  raw clinical rows: %d × %d cols\n", nrow(raw), ncol(raw)))

# Explicit column rename per the verified Cell Press mmc2 layout
required_cols <- c(
  `Study_ID`                       = "sample_id",
  `Age`                            = "age_at_dx",
  `AgeGroup`                       = "age_group",
  `Gender`                         = "sex",
  `histology`                      = "histology",
  `Met status (1 Met, 0 M0)`       = "metastatic",
  `Dead`                           = "os_event",
  `OS (years)`                     = "os_time_years",
  `Subgroup`                       = "subgroup",
  `Subtype`                        = "subtype"
)
missing_in_supp <- setdiff(names(required_cols), names(raw))
if (length(missing_in_supp) > 0) {
  stop(sprintf("mmc2.xlsx is missing expected columns: %s\nInspect supp table and align with the script's verified layout (rows 1-2 of mmc2).",
               paste(missing_in_supp, collapse = ", ")))
}

clinical <- raw[, names(required_cols), drop = FALSE]
names(clinical) <- unname(required_cols)

# Coerce types
clinical$age_at_dx    <- suppressWarnings(as.numeric(clinical$age_at_dx))
clinical$os_event     <- suppressWarnings(as.integer(clinical$os_event))
clinical$os_time_years <- suppressWarnings(as.numeric(clinical$os_time_years))
clinical$sample_id    <- as.character(clinical$sample_id)
clinical$subgroup     <- as.character(clinical$subgroup)
clinical$subtype      <- as.character(clinical$subtype)

# Drop rows that don't have a valid Study_ID (defensive — should be all 763 OK)
keep <- !is.na(clinical$sample_id) &
        grepl("^MB_SubtypeStudy_", clinical$sample_id)
clinical <- clinical[keep, , drop = FALSE]
cat(sprintf("  valid clinical rows (MB_SubtypeStudy_*): %d (expecting 763)\n", nrow(clinical)))
stopifnot(nrow(clinical) == 763)

# Summary
cat("\n  Subgroup distribution:\n")
print(table(clinical$subgroup, useNA = "ifany"))
cat("\n  12-subtype distribution:\n")
print(table(clinical$subtype, useNA = "ifany"))

with_surv <- !is.na(clinical$os_time_years) & !is.na(clinical$os_event)
cat(sprintf("\n  With-survival n: %d / %d (%.1f%%); events = %d\n",
            sum(with_surv), nrow(clinical), 100 * mean(with_surv),
            sum(clinical$os_event[with_surv], na.rm = TRUE)))

saveRDS(clinical, file.path(data_dir, "cavalli_clinical.rds"))
cat(sprintf("\n[%s] Done. Saved cavalli_clinical.rds (763 rows × %d cols).\n",
            format(Sys.time()), ncol(clinical)))
cat("Next: scripts/build_cavalli_metadata.R\n")
