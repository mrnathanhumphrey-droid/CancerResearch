# Paper 5 — run two-block BIDIFAC+ on Cavalli methylation + expression.
#
# Mirrors Park-Lock-Hoadley 2022's BIDIFAC+ setup (Lock 2022 Paper 1's
# upstream omic-decomposition step). Two blocks here (methylation + expression)
# scale down from Lock's four-block TCGA setup; BIDIFAC+ algorithm is unchanged.
#
# Reads:  data/cavalli_meth.rds, data/cavalli_expr.rds, data/cavalli_clinical.rds
# Writes: data/bidifac_components.rds (factor scores: samples × n_components)
#         data/bidifac_loadings.rds   (per-block loadings)
#         data/bidifac_diagnostics.rds (rank selection trace)

user_lib <- file.path(Sys.getenv("LOCALAPPDATA"), "R/win-library/4.6")
.libPaths(c(user_lib, .libPaths()))

# Install lockEF/bidifac if not present (R package from GitHub)
if (!requireNamespace("bidifac", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", lib = user_lib)
  remotes::install_github("lockEF/bidifac", lib = user_lib, upgrade = "never")
}
library(bidifac)

repo_root <- "C:/Cancer Research/Paper5_Medulloblastoma_StructuralPrior"
data_dir <- file.path(repo_root, "data")

cat(sprintf("[%s] Loading Cavalli omic blocks\n", format(Sys.time())))
meth <- readRDS(file.path(data_dir, "cavalli_meth.rds"))
expr <- readRDS(file.path(data_dir, "cavalli_expr.rds"))
clin <- readRDS(file.path(data_dir, "cavalli_clinical.rds"))

# --- Sample alignment -----------------------------------------------------
# Verified join key (Cavalli mmc2 inspection 2026-05-15):
#   GEO pData column `title` carries `MB_SubtypeStudy_NNNNN` strings (Sample_title
#   from GSE85217 / GSE85212 series_matrix), which match the `Study_ID` column
#   in Cavalli's mmc2.xlsx and are saved as `sample_id` in cavalli_clinical.rds.
expr_pdata <- readRDS(file.path(data_dir, "cavalli_expr_pdata.rds"))
meth_pdata <- readRDS(file.path(data_dir, "cavalli_meth_pdata.rds"))

match_key_expr <- as.character(expr_pdata$title)
match_key_meth <- as.character(meth_pdata$title)
stopifnot(any(grepl("^MB_SubtypeStudy_", match_key_expr)))
stopifnot(any(grepl("^MB_SubtypeStudy_", match_key_meth)))

common <- Reduce(intersect, list(clin$sample_id, match_key_expr, match_key_meth))
cat(sprintf("Common samples across expr / meth / clinical: %d\n", length(common)))
stopifnot(length(common) >= 700)  # sanity floor; Cavalli should be 763

idx_expr <- match(common, match_key_expr)
idx_meth <- match(common, match_key_meth)
idx_clin <- match(common, clin$sample_id)

expr_aligned <- expr[, idx_expr]
meth_aligned <- meth[, idx_meth]
clin_aligned <- clin[idx_clin, ]
colnames(expr_aligned) <- common
colnames(meth_aligned) <- common
clin_aligned$sample_id <- common

# --- Pre-standardize each block (row-centered, unit Frobenius) -----------
standardize_block <- function(X) {
  X <- X - rowMeans(X, na.rm = TRUE)                       # center genes/CpGs
  s <- sqrt(sum(X^2, na.rm = TRUE) / length(X))
  X / s
}
expr_std <- standardize_block(expr_aligned)
meth_std <- standardize_block(meth_aligned)
cat(sprintf("expr_std: %d × %d ; meth_std: %d × %d\n",
            nrow(expr_std), ncol(expr_std),
            nrow(meth_std), ncol(meth_std)))

# --- Two-block BIDIFAC+ ---------------------------------------------------
cat(sprintf("[%s] Running bidifacplus on two blocks\n", format(Sys.time())))
# bidifacplus expects a list of matrices, all with matched columns (samples).
blocks <- list(expr = expr_std, meth = meth_std)
fit <- bidifacplus(
  X_list   = blocks,
  R_max    = 30,                       # max rank per block; algorithm prunes via penalty
  verbose  = TRUE,
  tol      = 1e-4,
  max_iter = 200
)

# Output: fit produces factor matrices analogous to Lock 2022's XYC covariates.
# The shared + per-block components together form the candidate covariates
# for the per-subtype screening pass.
saveRDS(fit$V_shared, file.path(data_dir, "bidifac_components_shared.rds"))
saveRDS(fit$V_block, file.path(data_dir, "bidifac_components_per_block.rds"))
saveRDS(fit$U, file.path(data_dir, "bidifac_loadings.rds"))
saveRDS(fit$diagnostics, file.path(data_dir, "bidifac_diagnostics.rds"))

# Combined per-sample score matrix (samples × all components from all blocks)
# This is the Cavalli analog of Lock's XYC matrix (samples × 68 covariates).
all_components <- if (is.list(fit$V_shared) && is.list(fit$V_block)) {
  cbind(fit$V_shared, do.call(cbind, fit$V_block))
} else {
  cbind(fit$V_shared, fit$V_block)
}
# Standardize each component column
all_components <- scale(all_components)
rownames(all_components) <- common
saveRDS(all_components, file.path(data_dir, "bidifac_components.rds"))
saveRDS(clin_aligned,   file.path(data_dir, "cavalli_clinical_aligned.rds"))

cat(sprintf("BIDIFAC+ done: %d samples × %d total components.\n",
            nrow(all_components), ncol(all_components)))
cat(sprintf("Next: scripts/screening_paper5.R\n"))
