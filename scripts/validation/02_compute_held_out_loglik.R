# Compute held-out log-likelihood for each fitted model on the 20% test set.
#
# For each test patient i in cancer c with covariates X_i, observed Y_i (or
# censored at C_i), and given posterior samples θ_s = (β_c^s, σ²^s):
#
#   uncensored:  log p(Y_i | θ_s) = dnorm(log Y_i; μ_i^s, σ^s, log = TRUE)
#   censored:    log p(Y_i ≥ C_i | θ_s) = pnorm( (μ_i^s − log C_i) / σ^s, log.p = TRUE )
#
# where μ_i^s = X_i^T β_c^s. Then:
#   pointwise predictive density: LPPD_i = log mean_s exp(log p(Y_i | θ_s))
#   total log-likelihood: sum_i LPPD_i
#
# Posterior thinning: post-burn 50k iters, thin every 400 → 125 per chain × 4
# chains = 500 samples per model.
#
# Outputs:
#   loglik_per_patient_<spec>.csv — per-patient LPPD with cancer label
#   loglik_summary.csv            — per-spec total + per-cancer aggregate

user_lib <- file.path(Sys.getenv("LOCALAPPDATA"), "R/win-library/4.6")
.libPaths(c(user_lib, .libPaths()))

run_dir <- "C:/FkCancer/runs/run_validation_2026-05-09"
load(file.path(run_dir, "data/test_data.rda"))

# Recover covariates_in_model and covariates_by_cancer from the training data
load(file.path(run_dir, "data/train_data.rda"))
covariates_by_cancer <- lapply(Covariates_train, function(c) as.numeric(colnames(c)))
cancer_types <- names(Covariates_train)
n_cancer <- length(cancer_types)

BURN <- 50001L
THIN <- 400L

logsumexp <- function(x) {
  m <- max(x); m + log(sum(exp(x - m)))
}

compute_lppd_for_spec <- function(spec) {
  cat(sprintf("\n=== Computing LPPD for spec: %s ===\n", spec))
  out_dir <- file.path(run_dir, sprintf("gibbs_%s", spec))
  chain_files <- file.path(out_dir, sprintf("chain_%d.rda", 1:4))
  stopifnot(all(file.exists(chain_files)))

  # Build a list of (per-cancer) beta-sample matrices: rows = posterior samples,
  # cols = covariates_by_cancer[[c]]. Plus a sigma^2 vector aligned with rows.
  betas_per_cancer <- vector("list", n_cancer)
  for (c in seq_len(n_cancer)) {
    betas_per_cancer[[c]] <- list()
  }
  sigma2_samples <- numeric(0)

  for (j in 1:4) {
    cat(sprintf("  loading chain %d ...\n", j))
    e <- new.env(); load(chain_files[j], envir = e)
    out_obj <- e$out
    # Post-burn betas: list per cancer of 100000 iters; thin every 400
    keep <- seq(BURN, 100000L, by = THIN)  # ~125 indices per chain
    s2 <- out_obj$sigma2[keep]
    sigma2_samples <- c(sigma2_samples, s2)
    for (c in seq_len(n_cancer)) {
      bc <- out_obj$betas[[c]]
      mat <- do.call(rbind, bc[keep])     # n_keep x n_cov_c
      betas_per_cancer[[c]] <- rbind(
        if (is.null(betas_per_cancer[[c]]) ||
            (is.list(betas_per_cancer[[c]]) && length(betas_per_cancer[[c]]) == 0)) NULL
        else betas_per_cancer[[c]],
        mat)
    }
    rm(out_obj, e); gc(verbose = FALSE)
  }
  S <- length(sigma2_samples)
  cat(sprintf("  total posterior samples (thinned): %d\n", S))

  # Compute per-patient LPPD for the test set
  rows <- list()
  total_lppd <- 0
  for (c in seq_len(n_cancer)) {
    Xte <- as.matrix(Covariates_test[[c]])
    Yte <- Survival_test[[c]]
    Cte <- Censored_test[[c]]
    n_test <- length(Yte)
    if (n_test == 0) next
    Bm <- betas_per_cancer[[c]]    # S x n_cov_c
    stopifnot(ncol(Xte) == ncol(Bm))
    mu_mat <- Xte %*% t(Bm)        # n_test x S; each col is a sample's mu vector
    sd_mat <- matrix(sqrt(sigma2_samples), n_test, S, byrow = TRUE)
    censored_flag <- is.na(Yte)
    log_dens_mat <- matrix(NA_real_, n_test, S)
    for (i in seq_len(n_test)) {
      if (!censored_flag[i]) {
        log_dens_mat[i,] <- dnorm(log(Yte[i]), mean = mu_mat[i,], sd = sd_mat[i,],
                                  log = TRUE)
      } else {
        # log P(Y >= C) for log-normal: log(1 - Phi((log C - mu)/sigma))
        # = pnorm((mu - log C)/sigma, log.p = TRUE)
        z <- (mu_mat[i,] - log(Cte[i])) / sd_mat[i,]
        log_dens_mat[i,] <- pnorm(z, log.p = TRUE)
      }
    }
    # Per-patient LPPD = logsumexp(row) - log(S)
    lppd_per_patient <- apply(log_dens_mat, 1, logsumexp) - log(S)
    cancer_total <- sum(lppd_per_patient)
    total_lppd <- total_lppd + cancer_total
    cat(sprintf("    %-5s n_test=%3d censored=%2d  lppd=%9.2f\n",
                cancer_types[c], n_test, sum(censored_flag), cancer_total))
    for (i in seq_len(n_test)) {
      rows[[length(rows) + 1L]] <- data.frame(
        cancer  = cancer_types[c],
        patient = i,
        censored = censored_flag[i],
        lppd    = lppd_per_patient[i],
        stringsAsFactors = FALSE
      )
    }
  }
  cat(sprintf("  TOTAL test LPPD (%s): %.2f\n", spec, total_lppd))
  per_pat <- do.call(rbind, rows)
  out_csv <- file.path(run_dir, sprintf("loglik_per_patient_%s.csv", spec))
  write.csv(per_pat, out_csv, row.names = FALSE)
  cat(sprintf("  Saved %s\n", out_csv))
  list(spec = spec, total_lppd = total_lppd, per_patient = per_pat)
}

specs <- c("baseline_train", "primary", "secondary")
results <- list()
for (sp in specs) {
  spec_dir <- file.path(run_dir, sprintf("gibbs_%s", sp))
  if (!all(file.exists(file.path(spec_dir, sprintf("chain_%d.rda", 1:4))))) {
    cat(sprintf("[SKIP %s] not all chains present\n", sp))
    next
  }
  results[[sp]] <- compute_lppd_for_spec(sp)
}

if (length(results) >= 1) {
  summary_df <- data.frame(
    spec = names(results),
    total_lppd = sapply(results, function(r) r$total_lppd),
    n_test_patients = sapply(results, function(r) nrow(r$per_patient)),
    stringsAsFactors = FALSE
  )
  if ("baseline_train" %in% names(results)) {
    summary_df$delta_vs_baseline <- summary_df$total_lppd -
      summary_df$total_lppd[summary_df$spec == "baseline_train"]
  }
  write.csv(summary_df, file.path(run_dir, "loglik_summary.csv"), row.names = FALSE)
  cat("\n=== Summary ===\n"); print(summary_df, row.names = FALSE)
}
