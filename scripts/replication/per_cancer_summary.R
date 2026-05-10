user_lib <- file.path(Sys.getenv("LOCALAPPDATA"), "R/win-library/4.6")
.libPaths(c(user_lib, .libPaths()))

run_dir <- "C:/FkCancer/runs/run_paper_a_gibbs_2026-05-09"
files <- file.path(run_dir, sprintf(
  "GibbsSamplingResults_100kiters_Chain%d_WithAge_StandardizedPredictors.rda", 1:4))
stopifnot(all(file.exists(files)))

load("C:/FkCancer/repos/HierarchicalSS_PanCanPanOmics/XYC_V2_WithAge_StandardizedPredictors.rda")
cancer_types <- names(Covariates)
n_cancer <- length(cancer_types)
cov_per_cancer <- lapply(Covariates, function(x) as.numeric(colnames(x)))

n_iter <- 100001L
burn   <- 50001L
keep_idx <- (burn + 1L):n_iter

# Pre-allocate per-cancer matrices: list of 29 entries, each a 4-chain stacked
# matrix [4*n_keep, n_cov_for_cancer] for betas and gammas
n_keep <- length(keep_idx)
betas_stack  <- vector("list", n_cancer)
gamma_stack  <- vector("list", n_cancer)
for (c in seq_len(n_cancer)) {
  k <- length(cov_per_cancer[[c]])
  betas_stack[[c]]  <- matrix(NA_real_, 4 * n_keep, k,
                              dimnames = list(NULL, as.character(cov_per_cancer[[c]])))
  gamma_stack[[c]]  <- matrix(NA_integer_, 4 * n_keep, k,
                              dimnames = list(NULL, as.character(cov_per_cancer[[c]])))
}

for (j in 1:4) {
  cat(sprintf("[%s] Loading chain %d ...\n", format(Sys.time()), j))
  e <- new.env()
  load(files[j], envir = e)
  obj <- e[[ls(e)[1]]]
  for (c in seq_len(n_cancer)) {
    k <- length(cov_per_cancer[[c]])
    # betas: list of 100000 iterations (no +1)
    bb <- obj$betas[[c]]
    n_b <- length(bb)               # 100000
    # gammas: list of 100001 iterations, each a list with 1 element
    gg <- obj$gamma[[c]]
    n_g <- length(gg)               # 100001
    # Use post-burn slice consistent with parameters; betas idx (50001:100000) -> 50000 draws
    # gamma idx (50002:100001) -> 50000 draws (skip the very first init element)
    keep_b <- (n_b - n_keep + 1L):n_b
    keep_g <- (n_g - n_keep + 1L):n_g
    rng <- ((j - 1L) * n_keep + 1L):(j * n_keep)
    betas_stack[[c]][rng, ]  <- do.call(rbind, bb[keep_b])
    gamma_stack[[c]][rng, ]  <- do.call(rbind, gg[keep_g])
  }
  rm(obj, e); gc(verbose = FALSE)
}
cat(sprintf("[%s] All 4 chains stacked.\n", format(Sys.time())))

# Per-cancer × covariate summary
rows <- list()
for (c in seq_len(n_cancer)) {
  bb <- betas_stack[[c]]
  gg <- gamma_stack[[c]]
  k  <- ncol(bb)
  for (j in 1:k) {
    rows[[length(rows) + 1L]] <- data.frame(
      cancer    = cancer_types[c],
      cov_id    = colnames(bb)[j],
      beta_mean = mean(bb[, j]),
      beta_sd   = sd(bb[, j]),
      beta_q025 = quantile(bb[, j], 0.025),
      beta_q500 = quantile(bb[, j], 0.5),
      beta_q975 = quantile(bb[, j], 0.975),
      pip       = mean(gg[, j] == 1L),
      stringsAsFactors = FALSE
    )
  }
}
out <- do.call(rbind, rows)
rownames(out) <- NULL
csv_path <- file.path(run_dir, "per_cancer_beta_pip.csv")
write.csv(out, csv_path, row.names = FALSE)
cat(sprintf("Saved per-cancer table: %s   (%d rows)\n", csv_path, nrow(out)))

# Summary by cancer
cat("\n=== Per-cancer summary (k = #covariates, mean PIP, n PIP > 0.5, n PIP > 0.9) ===\n")
agg <- aggregate(cbind(pip, abs_beta = abs(beta_mean)) ~ cancer, data = out,
                 FUN = function(v) c(mean = mean(v), n_hi5 = sum(v > 0.5), n_hi9 = sum(v > 0.9)))
agg <- do.call(data.frame, agg)
agg$k <- as.integer(table(out$cancer)[agg$cancer])
agg <- agg[, c("cancer", "k", "pip.mean", "pip.n_hi5", "pip.n_hi9", "abs_beta.mean")]
agg <- agg[order(-agg$pip.mean), ]
print(agg, row.names = FALSE, digits = 4)

# Top-15 strongest cancer×covariate effects (by PIP, then by |beta|)
cat("\n=== Top-15 cancer x covariate by PIP (PIP > 0.7) ===\n")
strong <- out[out$pip > 0.7, ]
strong <- strong[order(-strong$pip, -abs(strong$beta_mean)), ]
print(head(strong, 15), row.names = FALSE, digits = 4)

# Show the cancer×cov rows where the 95% CI excludes 0 AND PIP > 0.5
cat("\n=== Cancer x covariate rows: PIP>0.5 AND 95% CI excludes 0 ===\n")
sig <- out[out$pip > 0.5 & (out$beta_q025 > 0 | out$beta_q975 < 0), ]
sig <- sig[order(-sig$pip, -abs(sig$beta_mean)), ]
print(sig, row.names = FALSE, digits = 4)

cat("\nDone.\n")
