# Paper 5 — leakage-clean BIDIFAC+ run.
#
# Diagnostic finding (2026-05-17, post-hoc): The original BIDIFAC+ fit on
# the full Cavalli cohort (n=763) including the held-out 20% test fold.
# Components therefore encode test-set covariance, contaminating the
# held-out LPPD verdict (+349 nats, 60× per-patient Lock's Paper 1).
#
# This script re-runs BIDIFAC+ on training fold ONLY (n=494), then projects
# test fold (n=118) onto the training-derived module loadings via
# least-squares. Pre-registration deviation: DEVIATIONS.md Entry 001.
#
# Reads:  data/cavalli_meth.rds, data/cavalli_expr.rds, data/cavalli_clinical.rds,
#         data/cavalli_expr_pdata.rds,
#         reference/paper5_split_indices.rda  (locked train/test split)
# Writes: data/bidifac_components_clean.rds      (samples × n_comp, train+test scores
#                                                  combined, aligned to common)
#         data/bidifac_loadings_clean.rds       (training-derived loadings + p_ind)
#         data/bidifac_modules_clean.rds        (raw module list)
#         data/bidifac_module_index_clean.rds
#         data/bidifac_norm_clean.rds           (row-means + Frobenius scales,
#                                                  TRAINING ONLY, for downstream
#                                                  audit)
#         data/cavalli_clinical_aligned.rds     (re-saved with explicit train/test col)
#         logs/bidifac_clean.log                 (progress + timing)

user_lib <- file.path(Sys.getenv("LOCALAPPDATA"), "R/win-library/4.6")
.libPaths(c(user_lib, .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"))

for (pkg in c("RSpectra", "matrixStats")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, lib = user_lib)
}
library(RSpectra)
library(matrixStats)

repo_root <- "C:/Cancer Research/Paper5_Medulloblastoma_StructuralPrior"
data_dir  <- file.path(repo_root, "data")
ref_dir   <- file.path(repo_root, "reference")
log_dir   <- file.path(repo_root, "logs")
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)

log_path <- file.path(log_dir, "bidifac_clean.log")
log_msg <- function(...) {
  msg <- sprintf("[%s] %s", format(Sys.time()), sprintf(...))
  cat(msg, "\n")
  cat(msg, "\n", file = log_path, append = TRUE)
  flush.console()
}

# --- Source BIDIFAC+ scripts (cached from prior run) --------------------
src_cache <- file.path(data_dir, "bidifac_src")
stopifnot(dir.exists(src_cache))
for (f in c("BIDIFAC.R", "EV_BIDIFAC_Functions.R",
            "bidifac.plus.R", "bidifac.plus.impute.R", "bidifac.plus.given.R")) {
  src_path <- file.path(src_cache, f)
  if (file.exists(src_path)) {
    log_msg("sourcing %s", f)
    tryCatch(source(src_path),
             error = function(e) log_msg("[WARN] %s failed: %s", f,
                                         conditionMessage(e)))
  }
}
stopifnot(exists("bidifac.plus"))

# --- Load and align ----------------------------------------------------
log_msg("Loading Cavalli omic blocks")
meth <- readRDS(file.path(data_dir, "cavalli_meth.rds"))
expr <- readRDS(file.path(data_dir, "cavalli_expr.rds"))
clin <- readRDS(file.path(data_dir, "cavalli_clinical.rds"))
expr_pdata <- readRDS(file.path(data_dir, "cavalli_expr_pdata.rds"))

match_key_expr <- as.character(expr_pdata$title)
colnames(expr) <- match_key_expr[match(colnames(expr), expr_pdata$geo_accession)]

common <- Reduce(intersect, list(clin$sample_id, colnames(expr), colnames(meth)))
log_msg("Common samples across expr / meth / clinical: %d", length(common))
stopifnot(length(common) >= 700)

idx_clin <- match(common, clin$sample_id)
expr_aligned <- expr[, common, drop = FALSE]
meth_aligned <- meth[, common, drop = FALSE]
clin_aligned <- clin[idx_clin, ]
clin_aligned$sample_id <- common

# --- Map training/test indices to `common` ordering --------------------
# `paper5_split_indices.rda` indices are positions into clin_s (the
# with-survival subset). Need to map those back to `common` positions.
load(file.path(ref_dir, "paper5_split_indices.rda"))

with_surv <- !is.na(clin$os_time_years) & !is.na(clin$os_event)
clin_s <- clin[with_surv, ]
n_surv <- nrow(clin_s)
log_msg("With-survival n: %d (train + test = %d + %d = %d)",
        n_surv,
        sum(sapply(train_idx, length)),
        sum(sapply(test_idx, length)),
        sum(sapply(train_idx, length)) + sum(sapply(test_idx, length)))

# Map clin_s indices to common positions
train_sample_ids <- unique(unlist(lapply(train_idx, function(idx) clin_s$sample_id[idx])))
test_sample_ids  <- unique(unlist(lapply(test_idx,  function(idx) clin_s$sample_id[idx])))
stopifnot(length(intersect(train_sample_ids, test_sample_ids)) == 0L)

# Which samples in `common` are train vs test (some may be neither if
# they were dropped from the with-survival subset)
train_cols <- match(train_sample_ids, common)
test_cols  <- match(test_sample_ids,  common)
train_cols <- train_cols[!is.na(train_cols)]
test_cols  <- test_cols[!is.na(test_cols)]
log_msg("BIDIFAC+ train cols: %d ; test cols: %d ; remainder (no-surv): %d",
        length(train_cols), length(test_cols),
        length(common) - length(train_cols) - length(test_cols))

# --- Variance pre-filter on BOTH blocks ---------------------------------
# CRITICAL: variance is computed on TRAINING SAMPLES ONLY to keep the
# feature selection itself leakage-free.
K_METH_TOP <- 5000L
K_EXPR_TOP <- 5000L
log_msg("Variance pre-filter (train-fold variance only): top %d meth / top %d expr",
        K_METH_TOP, K_EXPR_TOP)
meth_var <- matrixStats::rowVars(meth_aligned[, train_cols, drop = FALSE], na.rm = TRUE)
top_cpg <- order(meth_var, decreasing = TRUE)[seq_len(min(K_METH_TOP, length(meth_var)))]
meth_aligned <- meth_aligned[top_cpg, , drop = FALSE]
expr_var <- matrixStats::rowVars(expr_aligned[, train_cols, drop = FALSE], na.rm = TRUE)
top_expr <- order(expr_var, decreasing = TRUE)[seq_len(min(K_EXPR_TOP, length(expr_var)))]
expr_aligned <- expr_aligned[top_expr, , drop = FALSE]
log_msg("After filter: meth %d × %d ; expr %d × %d",
        nrow(meth_aligned), ncol(meth_aligned),
        nrow(expr_aligned), ncol(expr_aligned))

# --- NA imputation: row means from TRAINING SAMPLES only ----------------
impute_train_mean <- function(X, train_cols) {
  rm_train <- rowMeans(X[, train_cols, drop = FALSE], na.rm = TRUE)
  na_idx <- which(is.na(X), arr.ind = TRUE)
  if (nrow(na_idx) > 0L) X[na_idx] <- rm_train[na_idx[, 1L]]
  X
}
meth_aligned <- impute_train_mean(meth_aligned, train_cols)
expr_aligned[is.na(expr_aligned)] <- 0  # expr was already imputed in prior pass

# --- Standardize using TRAINING fold stats only -------------------------
standardize_train_apply <- function(X, train_cols) {
  rm_train <- rowMeans(X[, train_cols, drop = FALSE], na.rm = TRUE)
  Xc <- X - rm_train
  s_train <- sqrt(sum(Xc[, train_cols]^2, na.rm = TRUE) /
                  (nrow(X) * length(train_cols)))
  list(X = Xc / s_train, row_mean = rm_train, scale = s_train)
}
meth_norm <- standardize_train_apply(meth_aligned, train_cols)
expr_norm <- standardize_train_apply(expr_aligned, train_cols)
meth_std <- meth_norm$X
expr_std <- expr_norm$X
log_msg("Standardization: meth scale=%.6g ; expr scale=%.6g (train-only)",
        meth_norm$scale, expr_norm$scale)

# Save normalization params for downstream audit
saveRDS(list(meth_row_mean = meth_norm$row_mean, meth_scale = meth_norm$scale,
             expr_row_mean = expr_norm$row_mean, expr_scale = expr_norm$scale,
             top_cpg = top_cpg, top_expr = top_expr,
             train_cols = train_cols, test_cols = test_cols,
             common = common),
        file.path(data_dir, "bidifac_norm_clean.rds"))

# --- Build X0_train (features × n_train) --------------------------------
n_meth <- nrow(meth_std)
n_expr <- nrow(expr_std)
X0_full <- rbind(meth_std, expr_std)
X0_train <- X0_full[, train_cols, drop = FALSE]
X0_test  <- X0_full[, test_cols,  drop = FALSE]
n_train <- ncol(X0_train)
n_test  <- ncol(X0_test)
log_msg("X0_train: %d feats × %d train samples ; X0_test: %d × %d",
        nrow(X0_train), n_train, nrow(X0_test), n_test)

p_ind <- list(seq_len(n_meth), (n_meth + 1L):(n_meth + n_expr))
n_ind_train <- list(seq_len(n_train))

# --- Run BIDIFAC+ on TRAINING FOLD ONLY ---------------------------------
log_msg("Starting BIDIFAC+ on training fold (n=%d). Expect ~7-8h wall.", n_train)
t0 <- Sys.time()
set.seed(20260517L)   # distinct from full-data run (seed 20260516L)
res <- bidifac.plus(
  X0      = X0_train,
  p.ind   = p_ind,
  n.ind   = n_ind_train,
  max.comb  = 10,
  num.comp  = 10,
  max.iter  = 200,
  conv.thresh = 0.001,
  temp.iter = 50
)
t1 <- Sys.time()
wall_min <- as.numeric(difftime(t1, t0, units = "mins"))
log_msg("BIDIFAC+ training-fold fit done in %.1f min", wall_min)

# --- Extract per-module loadings + training scores; project test --------
module_block <- integer(0)
all_train_scores <- matrix(0, nrow = n_train, ncol = 0,
                           dimnames = list(common[train_cols], character(0)))
all_test_scores  <- matrix(0, nrow = n_test,  ncol = 0,
                           dimnames = list(common[test_cols],  character(0)))
all_loadings <- list()

for (k in seq_along(res$S)) {
  S_k <- res$S[[k]]
  if (all(S_k == 0)) next
  pk <- res$p.ind.list[[k]]
  in_meth <- any(pk %in% p_ind[[1]])
  in_expr <- any(pk %in% p_ind[[2]])
  block_label <- if (in_meth && in_expr) 12L else if (in_meth) 1L else 2L

  # SVD of training-active sub-matrix
  active <- S_k[pk, res$n.ind.list[[k]], drop = FALSE]
  sv <- svd(active)
  r <- sum(sv$d > 1e-6 * sv$d[1])
  if (r == 0) next
  L_k <- sv$u[, seq_len(r), drop = FALSE] %*% diag(sv$d[seq_len(r)], r, r)  # length(pk) × r
  V_k_train <- sv$v[, seq_len(r), drop = FALSE]                              # n_train_active × r

  # Pad training scores back to full train sample set (zeros for non-active)
  scores_train_full <- matrix(0, nrow = n_train, ncol = r)
  scores_train_full[res$n.ind.list[[k]], ] <- V_k_train

  # Project test fold onto L_k via least-squares
  # X0_test[pk, ] is length(pk) × n_test
  # Want V_k_test (n_test × r) such that X0_test[pk, ] ≈ L_k %*% t(V_k_test)
  # → V_k_test = t(X0_test[pk, ]) %*% L_k %*% solve(t(L_k) %*% L_k)
  X_test_block <- X0_test[pk, , drop = FALSE]
  gram <- crossprod(L_k)                # r × r
  # Regularize tiny diagonal in case of near-rank-deficient gram
  diag(gram) <- diag(gram) + 1e-10 * max(diag(gram), 1e-10)
  V_k_test <- t(X_test_block) %*% L_k %*% solve(gram)   # n_test × r

  col_names_k <- sprintf("M%02d_%s_r%d",
                         k,
                         c("meth", "expr", "joint")[c(1, 2, 12) == block_label],
                         seq_len(r))
  colnames(scores_train_full) <- col_names_k
  colnames(V_k_test) <- col_names_k
  rownames(scores_train_full) <- common[train_cols]
  rownames(V_k_test) <- common[test_cols]

  module_block <- c(module_block, rep(block_label, r))
  all_train_scores <- cbind(all_train_scores, scores_train_full)
  all_test_scores  <- cbind(all_test_scores,  V_k_test)
  all_loadings[[length(all_loadings) + 1L]] <- list(module = k,
                                                    block = block_label,
                                                    loadings = L_k,
                                                    p_ind = pk)
}

# --- Combine train + test into common-aligned score matrix ---------------
# Standardize using TRAINING fold mean + sd (NOT global scale()), so test
# scores live in the same transformed space without seeing test data.
train_means <- colMeans(all_train_scores)
train_sds   <- apply(all_train_scores, 2, sd)
train_sds[train_sds < 1e-12] <- 1   # avoid 0/0 for zero-variance columns
scale_train <- function(M) sweep(sweep(M, 2, train_means, "-"), 2, train_sds, "/")
all_train_scaled <- scale_train(all_train_scores)
all_test_scaled  <- scale_train(all_test_scores)

# Allocate a master matrix indexed by `common`
all_scores <- matrix(NA_real_, nrow = length(common), ncol = ncol(all_train_scaled),
                     dimnames = list(common, colnames(all_train_scaled)))
all_scores[train_cols, ] <- all_train_scaled
all_scores[test_cols,  ] <- all_test_scaled
# Non-survival samples (neither train nor test) stay NA — they are not used
# downstream by Gibbs/LPPD.

log_msg("BIDIFAC+ components (train-only fit): %d", ncol(all_scores))
log_msg("  per-block tally (1=meth, 2=expr, 12=joint): %s",
        paste(table(module_block), collapse = " / "))
log_msg("  train scores: mean=%.4f sd=%.4f", mean(all_train_scaled, na.rm=TRUE),
        sd(all_train_scaled, na.rm=TRUE))
log_msg("  test  scores (projected): mean=%.4f sd=%.4f",
        mean(all_test_scaled, na.rm=TRUE), sd(all_test_scaled, na.rm=TRUE))

saveRDS(all_scores,    file.path(data_dir, "bidifac_components_clean.rds"))
saveRDS(all_loadings,  file.path(data_dir, "bidifac_loadings_clean.rds"))
saveRDS(res,           file.path(data_dir, "bidifac_modules_clean.rds"))
saveRDS(module_block,  file.path(data_dir, "bidifac_module_index_clean.rds"))

# Re-save cavalli_clinical_aligned with train/test column for audit
clin_aligned$bidifac_fold <- "neither"
clin_aligned$bidifac_fold[train_cols] <- "train"
clin_aligned$bidifac_fold[test_cols]  <- "test"
saveRDS(clin_aligned, file.path(data_dir, "cavalli_clinical_aligned.rds"))

log_msg("DONE. Components written. Train: %d/%d, Test: %d/%d.",
        length(train_cols), length(common),
        length(test_cols),  length(common))
