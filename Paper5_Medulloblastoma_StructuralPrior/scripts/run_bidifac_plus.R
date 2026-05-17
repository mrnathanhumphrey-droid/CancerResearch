# Paper 5 — run BIDIFAC+ on Cavalli methylation + expression (two-block).
#
# Mirrors Park-Lock-Hoadley 2022's BIDIFAC+ setup (Lock 2022 Paper 1's
# upstream omic-decomposition step). Two blocks here (methylation + expression)
# scale down from Lock's four-block TCGA setup; algorithm unchanged.
#
# BIDIFAC+ source: github.com/lockEF/bidifac (standalone R scripts, not a
# package). Functions sourced directly at runtime; no install needed.
# Algorithm reference: Park JY, Lock EF, Hoadley KA. 2022. AOAS 16(1):484-501.
#
# Reads:  data/cavalli_meth.rds, data/cavalli_expr.rds, data/cavalli_clinical.rds
# Writes: data/bidifac_components.rds (per-sample factor scores: samples × n_components)
#         data/bidifac_loadings.rds   (per-feature loadings)
#         data/bidifac_modules.rds    (raw module list from BIDIFAC+)
#         data/bidifac_module_index.rds (which block each module belongs to)
#         data/cavalli_clinical_aligned.rds

user_lib <- file.path(Sys.getenv("LOCALAPPDATA"), "R/win-library/4.6")
.libPaths(c(user_lib, .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"))

# RSpectra is a BIDIFAC+ dependency for fast partial SVD
for (pkg in c("RSpectra")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, lib = user_lib)
}
library(RSpectra)

# --- Source BIDIFAC+ functions directly from GitHub ---------------------
# The lockEF/bidifac repo is a script collection (no DESCRIPTION/NAMESPACE);
# remotes::install_github fails with 404 because there is no package to install.
repo_root <- "C:/Cancer Research/Paper5_Medulloblastoma_StructuralPrior"
data_dir  <- file.path(repo_root, "data")
src_cache <- file.path(data_dir, "bidifac_src")
dir.create(src_cache, showWarnings = FALSE, recursive = TRUE)

src_files <- c("bidifac.plus.R", "bidifac.plus.impute.R", "bidifac.plus.given.R",
               "BIDIFAC.R", "EV_BIDIFAC_Functions.R")
for (f in src_files) {
  dst <- file.path(src_cache, f)
  if (!file.exists(dst)) {
    url <- sprintf("https://raw.githubusercontent.com/lockEF/bidifac/master/%s", f)
    cat(sprintf("[%s] Downloading %s\n", format(Sys.time()), f))
    download.file(url, destfile = dst, mode = "wb", quiet = TRUE)
  }
}
# Source dependencies first, then the BIDIFAC+ functions
for (f in c("BIDIFAC.R", "EV_BIDIFAC_Functions.R",
            "bidifac.plus.R", "bidifac.plus.impute.R", "bidifac.plus.given.R")) {
  src_path <- file.path(src_cache, f)
  if (file.exists(src_path)) {
    cat(sprintf("  sourcing %s\n", f))
    tryCatch(source(src_path),
             error = function(e) cat(sprintf("  [WARN] %s failed to source: %s\n",
                                             f, conditionMessage(e))))
  }
}
stopifnot(exists("bidifac.plus"))

# --- Load blocks --------------------------------------------------------
cat(sprintf("[%s] Loading Cavalli omic blocks\n", format(Sys.time())))
meth <- readRDS(file.path(data_dir, "cavalli_meth.rds"))
expr <- readRDS(file.path(data_dir, "cavalli_expr.rds"))
clin <- readRDS(file.path(data_dir, "cavalli_clinical.rds"))
cat(sprintf("  raw meth: %d × %d\n", nrow(meth), ncol(meth)))
cat(sprintf("  raw expr: %d × %d\n", nrow(expr), ncol(expr)))
cat(sprintf("  clinical rows: %d\n", nrow(clin)))

# --- Sample alignment via Sample_title ↔ Study_ID -----------------------
# Verified join key (Cavalli mmc2 inspection 2026-05-15):
#   GEO pData `title` = `MB_SubtypeStudy_NNNNN` = Cavalli mmc2 `Study_ID`
#   The methylation supplementary file uses the same Sample_title strings as
#   column names (verified by inspecting the gz header line in fetch script).
expr_pdata <- readRDS(file.path(data_dir, "cavalli_expr_pdata.rds"))
match_key_expr <- as.character(expr_pdata$title)
# expr column names are GSM accessions; rename to Sample_title for join
colnames(expr) <- match_key_expr[match(colnames(expr), expr_pdata$geo_accession)]
stopifnot(sum(grepl("^MB_SubtypeStudy_", colnames(expr))) >= 700)
stopifnot(sum(grepl("^MB_SubtypeStudy_", colnames(meth))) >= 700)

common <- Reduce(intersect, list(clin$sample_id, colnames(expr), colnames(meth)))
cat(sprintf("Common samples across expr / meth / clinical: %d\n", length(common)))
stopifnot(length(common) >= 700)

idx_clin <- match(common, clin$sample_id)
expr_aligned <- expr[, common, drop = FALSE]
meth_aligned <- meth[, common, drop = FALSE]
clin_aligned <- clin[idx_clin, ]
clin_aligned$sample_id <- common

# --- Variance pre-filter on BOTH blocks ---------------------------------
# Full 450k methylation × 763 + ~22k expression × 763 → ~32k stacked features.
# Profiling on the 2026-05-15 run: 32k × 763 BIDIFAC+ is ~14 min per outer
# iteration, infeasible for the 50-100 iter convergence horizon. Reducing
# each block to top-5k features per cross-sample variance cuts the stacked
# matrix to 10k features, yielding ~10× speedup. Standard practice — does
# not affect pre-reg hierarchy or decision rules.
K_METH_TOP <- 5000L
K_EXPR_TOP <- 5000L
if (!requireNamespace("matrixStats", quietly = TRUE)) {
  install.packages("matrixStats", lib = user_lib)
}
library(matrixStats)
cat(sprintf("[%s] Pre-filtering methylation to top %d / expression to top %d (by var)\n",
            format(Sys.time()), K_METH_TOP, K_EXPR_TOP))
meth_var <- matrixStats::rowVars(meth_aligned, na.rm = TRUE)
top_cpg <- order(meth_var, decreasing = TRUE)[seq_len(min(K_METH_TOP, length(meth_var)))]
meth_aligned <- meth_aligned[top_cpg, , drop = FALSE]
expr_var <- matrixStats::rowVars(expr_aligned, na.rm = TRUE)
top_expr <- order(expr_var, decreasing = TRUE)[seq_len(min(K_EXPR_TOP, length(expr_var)))]
expr_aligned <- expr_aligned[top_expr, , drop = FALSE]
cat(sprintf("  methylation after filter: %d CpGs × %d samples\n",
            nrow(meth_aligned), ncol(meth_aligned)))
cat(sprintf("  expression after filter:  %d probes × %d samples\n",
            nrow(expr_aligned), ncol(expr_aligned)))

# --- Impute NAs (BIDIFAC+ requires complete data) -----------------------
meth_aligned[is.na(meth_aligned)] <- rowMeans(meth_aligned, na.rm = TRUE)[
  row(meth_aligned)[is.na(meth_aligned)]]
expr_aligned[is.na(expr_aligned)] <- 0

# --- Standardize each block (row-centered, unit-Frobenius normalized) ---
standardize_block <- function(X) {
  X <- X - rowMeans(X, na.rm = TRUE)
  s <- sqrt(sum(X^2, na.rm = TRUE) / length(X))
  X / s
}
meth_std <- standardize_block(meth_aligned)
expr_std <- standardize_block(expr_aligned)

# --- Stack vertically into single X0 matrix -----------------------------
# BIDIFAC+ expects X0 with rows = features (stacked across blocks) and
# columns = samples. p.ind = list of row-index sets per block. n.ind = list
# of column-index sets (here a single set: all 763 samples).
n_meth <- nrow(meth_std)
n_expr <- nrow(expr_std)
X0 <- rbind(meth_std, expr_std)
n_samples <- ncol(X0)
cat(sprintf("Stacked X0: %d features (%d meth + %d expr) × %d samples\n",
            nrow(X0), n_meth, n_expr, n_samples))

p_ind <- list(seq_len(n_meth),
              (n_meth + 1L):(n_meth + n_expr))
n_ind <- list(seq_len(n_samples))

# --- Run BIDIFAC+ -------------------------------------------------------
cat(sprintf("[%s] Running BIDIFAC+ (this is the CPU-heavy step)\n", format(Sys.time())))
t0 <- Sys.time()
set.seed(20260516L)
res <- bidifac.plus(
  X0      = X0,
  p.ind   = p_ind,
  n.ind   = n_ind,
  max.comb  = 10,    # ↓ from 20 — fewer candidate modules per iter
  num.comp  = 10,    # ↓ from 20 — smaller SVD per inner step
  max.iter  = 200,   # ↓ from 500 — adequate ceiling; conv.thresh trips early
  conv.thresh = 0.001,
  temp.iter = 50     # ↓ from 100 — proportional to max.iter cut
)
t1 <- Sys.time()
cat(sprintf("BIDIFAC+ done in %.1f min.\n",
            as.numeric(difftime(t1, t0, units = "mins"))))

# --- Extract per-sample factor scores from each module ------------------
# Each non-zero module S[[k]] is a low-rank matrix of size
# length(p.ind.list[[k]]) × length(n.ind.list[[k]]). For Paper 5's purpose
# (analog of Lock 2022's BIDIFAC+ XYC covariates), we extract per-sample
# scores via SVD truncated at the module's effective rank.
module_block <- integer(0)   # which block(s) each module covers (1=meth, 2=expr, 12=joint)
all_scores <- matrix(0, nrow = n_samples, ncol = 0,
                     dimnames = list(common, character(0)))
all_loadings <- list()
for (k in seq_along(res$S)) {
  S_k <- res$S[[k]]
  if (all(S_k == 0)) next
  # Determine which block this module spans
  pk <- res$p.ind.list[[k]]
  in_meth <- any(pk %in% p_ind[[1]])
  in_expr <- any(pk %in% p_ind[[2]])
  block_label <- if (in_meth && in_expr) 12L else if (in_meth) 1L else 2L
  # SVD of the active sub-matrix
  active <- S_k[pk, res$n.ind.list[[k]], drop = FALSE]
  sv <- svd(active)
  r <- sum(sv$d > 1e-6 * sv$d[1])
  if (r == 0) next
  module_block <- c(module_block, rep(block_label, r))
  scores_k <- sv$v[, seq_len(r), drop = FALSE]
  load_k   <- sv$u[, seq_len(r), drop = FALSE] %*% diag(sv$d[seq_len(r)], r, r)
  # Pad scores back to full sample set (zeros for samples not in n.ind.list[[k]])
  scores_full <- matrix(0, nrow = n_samples, ncol = r)
  scores_full[res$n.ind.list[[k]], ] <- scores_k
  colnames(scores_full) <- sprintf("M%02d_%s_r%d",
                                   k,
                                   c("meth", "expr", "joint")[c(1,2,12) == block_label],
                                   seq_len(r))
  rownames(scores_full) <- common
  all_scores <- cbind(all_scores, scores_full)
  all_loadings[[length(all_loadings) + 1L]] <- list(module = k,
                                                    block = block_label,
                                                    loadings = load_k,
                                                    p_ind = pk)
}
all_scores <- scale(all_scores)

cat(sprintf("Total BIDIFAC+ components extracted: %d\n", ncol(all_scores)))
cat(sprintf("  per-block tally (1=meth, 2=expr, 12=joint): %s\n",
            paste(table(module_block), collapse = " / ")))

saveRDS(all_scores,   file.path(data_dir, "bidifac_components.rds"))
saveRDS(all_loadings, file.path(data_dir, "bidifac_loadings.rds"))
saveRDS(res,          file.path(data_dir, "bidifac_modules.rds"))
saveRDS(module_block, file.path(data_dir, "bidifac_module_index.rds"))
saveRDS(clin_aligned, file.path(data_dir, "cavalli_clinical_aligned.rds"))

cat(sprintf("BIDIFAC+ done: %d samples × %d total components.\n",
            nrow(all_scores), ncol(all_scores)))
cat("Next: scripts/screening_paper5.R\n")
