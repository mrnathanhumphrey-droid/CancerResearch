# Paper 5 — diagnose the projected-test-score scale issue.
#
# Symptom: BIDIFAC+ log reported "train scores: sd=1.0 / test scores
# (projected): sd=98.3". We need to know:
#   1. Per-column scale ratio: which columns have train_sd ≈ test_sd vs
#      blown up?
#   2. Are the locked target covariates {43, 25, 23} in the blown-up set
#      or the stable set?
#   3. Plot the singular-value distribution per module — small d_min within
#      a module is the most likely cause.
#   4. Sanity-check downstream Gibbs fits: did the spike-and-slab zero out
#      the blown-up columns?

user_lib <- file.path(Sys.getenv("LOCALAPPDATA"), "R/win-library/4.6")
.libPaths(c(user_lib, .libPaths()))

repo_root <- "C:/Cancer Research/Paper5_Medulloblastoma_StructuralPrior"
data_dir  <- file.path(repo_root, "data")
res_dir   <- file.path(repo_root, "results")
ref_dir   <- file.path(repo_root, "reference")

X_all <- readRDS(file.path(data_dir, "bidifac_components_clean.rds"))
clin  <- readRDS(file.path(data_dir, "cavalli_clinical_aligned.rds"))
loadings <- readRDS(file.path(data_dir, "bidifac_loadings_clean.rds"))
target_table <- read.csv(file.path(ref_dir, "paper5_target_covariates_clean.csv"),
                         stringsAsFactors = FALSE)

# Train and test masks
train_mask <- clin$bidifac_fold == "train"
test_mask  <- clin$bidifac_fold == "test"
cat(sprintf("Train samples: %d ; Test samples: %d ; total cols: %d\n",
            sum(train_mask), sum(test_mask), ncol(X_all)))

# Per-column scale (already standardized by train-fold stats, so train_sd=1
# by construction; test_sd reveals the ratio)
train_sd <- apply(X_all[train_mask, , drop = FALSE], 2, sd)
test_sd  <- apply(X_all[test_mask,  , drop = FALSE], 2, sd)
train_mean <- colMeans(X_all[train_mask, , drop = FALSE])
test_mean  <- colMeans(X_all[test_mask,  , drop = FALSE])

cat("\n=== Per-column train_sd / test_sd ratio ===\n")
ratio <- test_sd / train_sd
cat(sprintf("median(test_sd / train_sd) = %.3f\n", median(ratio)))
cat(sprintf("max(test_sd / train_sd)    = %.3f (column %d)\n",
            max(ratio), which.max(ratio)))
cat(sprintf("min(test_sd / train_sd)    = %.3f (column %d)\n",
            min(ratio), which.min(ratio)))

# Distribution: how many columns are stable vs blown up?
breaks <- c(0, 0.5, 1.5, 5, 20, 100, 500, Inf)
tab <- cut(ratio, breaks = breaks, include.lowest = TRUE)
cat("\n=== Ratio bins ===\n")
print(table(tab))

cat("\n=== Top 10 worst columns (largest test_sd / train_sd) ===\n")
ord <- order(-ratio)
df_worst <- data.frame(col_idx = ord[1:10],
                       col_name = colnames(X_all)[ord[1:10]],
                       train_sd = round(train_sd[ord[1:10]], 4),
                       test_sd  = round(test_sd[ord[1:10]], 4),
                       ratio    = round(ratio[ord[1:10]], 3),
                       train_mean = round(train_mean[ord[1:10]], 4),
                       test_mean  = round(test_mean[ord[1:10]], 4))
print(df_worst, row.names = FALSE)

cat("\n=== Top 10 stable columns (ratio nearest 1) ===\n")
ord_stable <- order(abs(log(ratio)))
df_stable <- data.frame(col_idx = ord_stable[1:10],
                       col_name = colnames(X_all)[ord_stable[1:10]],
                       train_sd = round(train_sd[ord_stable[1:10]], 4),
                       test_sd  = round(test_sd[ord_stable[1:10]], 4),
                       ratio    = round(ratio[ord_stable[1:10]], 3))
print(df_stable, row.names = FALSE)

# === Are the locked covariates in the blown-up set? ===
# Note: covariates_in_model are 1-indexed component positions matching
# colnames(X_all) sequence 1..55 after intercept/age (those are cov_ids 0, 0.5).
# cov_id N corresponds to col_idx N of X_all (since X_all has 55 components).
cat("\n=== LOCKED TARGET COVARIATES — scale check ===\n")
for (i in seq_len(nrow(target_table))) {
  cid <- target_table$cov_id[i]
  cat(sprintf("  cov_id=%d (col %d, %s): train_sd=%.3f test_sd=%.3f ratio=%.3f\n",
              cid, cid, colnames(X_all)[cid],
              train_sd[cid], test_sd[cid], ratio[cid]))
}

# === Module-wise singular value spectrum ===
cat("\n=== Per-module SVD spectrum (showing smallest d_min within modules with r>1) ===\n")
for (k in seq_along(loadings)) {
  ld <- loadings[[k]]
  L_k <- ld$loadings
  # L_k columns = u %*% diag(d), so per-column norm = d. We can recover d.
  d_k <- sqrt(colSums(L_k^2))
  if (length(d_k) > 1L) {
    cat(sprintf("  Module %d (block %d): r=%d, d=%s, d_min/d_max=%.4g\n",
                ld$module, ld$block, length(d_k),
                paste(sprintf("%.3g", d_k), collapse = " "),
                min(d_k)/max(d_k)))
  }
}

# === Match modules to component columns ===
# Each row of loadings list has $module (k), $block, $loadings (length(pk) x r_k).
# The columns in X_all are named "M%02d_<block>_r<j>" — we can map module k <-> col.
mod_cols <- list()
for (k in seq_along(loadings)) {
  mod <- loadings[[k]]$module
  block <- loadings[[k]]$block
  block_tag <- c("meth", "expr", "joint")[c(1, 2, 12) == block]
  r_k <- ncol(loadings[[k]]$loadings)
  prefix <- sprintf("M%02d_%s_", mod, block_tag)
  cols <- which(startsWith(colnames(X_all), prefix))
  mod_cols[[as.character(mod)]] <- cols
  if (length(cols) != r_k) {
    cat(sprintf("  [WARN] module %d has r=%d but %d columns matched\n",
                mod, r_k, length(cols)))
  }
}

# For each locked covariate, find its module and that module's d-spectrum
cat("\n=== LOCKED COVS → their module's singular-value spectrum ===\n")
for (i in seq_len(nrow(target_table))) {
  cid <- target_table$cov_id[i]
  cn  <- colnames(X_all)[cid]
  for (k in seq_along(loadings)) {
    ld <- loadings[[k]]
    block_tag <- c("meth", "expr", "joint")[c(1, 2, 12) == ld$block]
    prefix <- sprintf("M%02d_%s_", ld$module, block_tag)
    if (startsWith(cn, prefix)) {
      d_k <- sqrt(colSums(ld$loadings^2))
      cat(sprintf("  cov_id=%d (%s) ← module %d (block %d): d=%s\n",
                  cid, cn, ld$module, ld$block,
                  paste(sprintf("%.3g", d_k), collapse = " ")))
      break
    }
  }
}

# === Spike-and-slab: did Gibbs zero out the blown-up columns? ===
# Load primary clean Gibbs chain 1 and compute PIPs (γ posterior mean).
cat("\n=== PIPs from clean primary Gibbs chain 1 (per cancer-type) ===\n")
chain_file <- file.path(res_dir, "gibbs_primary_clean", "chain_1.rda")
if (file.exists(chain_file)) {
  env <- new.env(); load(chain_file, envir = env); out <- env$out
  burn <- 50001L
  iters <- length(out$sigma2) - 1L
  draws <- (burn + 1L):iters
  n_cancer <- length(out$gamma)
  cat(sprintf("  n_cancer (axis units) = %d ; iters = %d\n", n_cancer, iters))

  # Worst-10 blown-up columns
  worst_cids <- ord[1:10]
  for (i in seq_along(worst_cids)) {
    cid <- worst_cids[i]
    cn <- colnames(X_all)[cid]
    # PIP averaged across all cancer-types
    pips <- numeric(n_cancer)
    for (j in seq_len(n_cancer)) {
      gg <- out$gamma[[j]]
      ref_names <- names(gg[[1]])
      cid_str <- as.character(cid)
      if (cid_str %in% ref_names) {
        # Need to extract the indicator at position matching cid_str
        col_pos <- match(cid_str, ref_names)
        pip_d <- vapply(gg[draws], function(g) g[col_pos], numeric(1))
        pips[j] <- mean(pip_d)
      } else {
        pips[j] <- NA_real_
      }
    }
    cat(sprintf("  worst#%d cov_id=%d (%s, ratio=%.1f): mean PIP across cancers = %.3f (max %.3f)\n",
                i, cid, cn, ratio[cid], mean(pips, na.rm = TRUE),
                max(pips, na.rm = TRUE)))
  }

  # Locked covs
  cat("\n  --- locked target covs ---\n")
  for (i in seq_len(nrow(target_table))) {
    cid <- target_table$cov_id[i]
    cn  <- colnames(X_all)[cid]
    pips <- numeric(n_cancer)
    for (j in seq_len(n_cancer)) {
      gg <- out$gamma[[j]]
      ref_names <- names(gg[[1]])
      cid_str <- as.character(cid)
      if (cid_str %in% ref_names) {
        col_pos <- match(cid_str, ref_names)
        pip_d <- vapply(gg[draws], function(g) g[col_pos], numeric(1))
        pips[j] <- mean(pip_d)
      } else {
        pips[j] <- NA_real_
      }
    }
    cat(sprintf("    cov_id=%d (%s, ratio=%.2f): mean PIP = %.3f (range %.3f-%.3f)\n",
                cid, cn, ratio[cid], mean(pips, na.rm = TRUE),
                min(pips, na.rm = TRUE), max(pips, na.rm = TRUE)))
  }
} else {
  cat(sprintf("  Chain file missing: %s\n", chain_file))
}

cat("\n=== Summary verdict ===\n")
cat(sprintf("  Stable columns (ratio in [0.5, 1.5]): %d / %d\n",
            sum(ratio >= 0.5 & ratio <= 1.5), length(ratio)))
cat(sprintf("  Locked covs all stable? %s\n",
            all(ratio[target_table$cov_id] >= 0.5 & ratio[target_table$cov_id] <= 1.5)))
