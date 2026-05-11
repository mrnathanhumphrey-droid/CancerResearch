# Held-out LPPD for each of the 18 extension specs.
#
# Reuses Paper 1's held-out 80/20 split (train_data + test_data from
# run_validation_2026-05-09/data/). Posterior thinning: post-burn 50k iters,
# thin every 400 -> 125 per chain x 4 chains = 500 samples per spec.
#
# Output:
#   loglik_per_patient_<spec>.csv  (one per spec)
#   loglik_summary_extension.csv   (all 18 specs + Paper 1 anchors)

user_lib <- file.path(Sys.getenv("LOCALAPPDATA"), "R/win-library/4.6")
.libPaths(c(user_lib, .libPaths()))

run_dir <- "C:/FkCancer/runs/run_extension_2026-05-10"
val_dir <- "C:/FkCancer/runs/run_validation_2026-05-09"

load(file.path(val_dir, "data/test_data.rda"))
load(file.path(val_dir, "data/train_data.rda"))
load(file.path(run_dir, "data/extension_specs.rda"))

covariates_by_cancer <- lapply(Covariates_train, function(c) as.numeric(colnames(c)))
cancer_types <- names(Covariates_train)
n_cancer <- length(cancer_types)

BURN <- 50001L
THIN <- 400L

logsumexp <- function(x) { m <- max(x); m + log(sum(exp(x - m))) }

compute_lppd_for_spec <- function(spec) {
  out_dir <- file.path(run_dir, sprintf("gibbs_%s", spec))
  chain_files <- file.path(out_dir, sprintf("chain_%d.rda", 1:4))
  if (!all(file.exists(chain_files))) {
    cat(sprintf("[SKIP %s] chain files missing\n", spec)); return(NULL)
  }
  cat(sprintf("\n=== %s ===\n", spec))

  betas_per_cancer <- vector("list", n_cancer)
  for (c in seq_len(n_cancer)) betas_per_cancer[[c]] <- list()
  sigma2_samples <- numeric(0)

  for (j in 1:4) {
    cat(sprintf("  loading chain %d ...\n", j))
    e <- new.env(); load(chain_files[j], envir = e); o <- e$out
    keep <- seq(BURN, 100000L, by = THIN)
    sigma2_samples <- c(sigma2_samples, o$sigma2[keep])
    for (c in seq_len(n_cancer)) {
      mat <- do.call(rbind, o$betas[[c]][keep])
      betas_per_cancer[[c]] <- if (length(betas_per_cancer[[c]]) == 0) mat
                                else rbind(betas_per_cancer[[c]], mat)
    }
    rm(o, e); gc(verbose = FALSE)
  }
  S <- length(sigma2_samples)
  cat(sprintf("  total posterior samples: %d\n", S))

  rows <- list()
  total_lppd <- 0
  for (c in seq_len(n_cancer)) {
    Xte <- as.matrix(Covariates_test[[c]])
    Yte <- Survival_test[[c]]
    Cte <- Censored_test[[c]]
    n_test <- length(Yte); if (n_test == 0) next
    Bm <- betas_per_cancer[[c]]
    stopifnot(ncol(Xte) == ncol(Bm))
    mu_mat <- Xte %*% t(Bm)
    sd_mat <- matrix(sqrt(sigma2_samples), n_test, S, byrow = TRUE)
    censored_flag <- is.na(Yte)
    log_dens_mat <- matrix(NA_real_, n_test, S)
    for (i in seq_len(n_test)) {
      if (!censored_flag[i]) {
        log_dens_mat[i,] <- dnorm(log(Yte[i]), mean = mu_mat[i,], sd = sd_mat[i,], log = TRUE)
      } else {
        z <- (mu_mat[i,] - log(Cte[i])) / sd_mat[i,]
        log_dens_mat[i,] <- pnorm(z, log.p = TRUE)
      }
    }
    lppd_per_patient <- apply(log_dens_mat, 1, logsumexp) - log(S)
    cancer_total <- sum(lppd_per_patient)
    total_lppd <- total_lppd + cancer_total
    cat(sprintf("    %-5s n_test=%3d censored=%2d  lppd=%9.2f\n",
                cancer_types[c], n_test, sum(censored_flag), cancer_total))
    for (i in seq_len(n_test)) {
      rows[[length(rows) + 1L]] <- data.frame(
        cancer = cancer_types[c], patient = i,
        censored = censored_flag[i], lppd = lppd_per_patient[i],
        stringsAsFactors = FALSE)
    }
  }
  cat(sprintf("  TOTAL test LPPD (%s): %.2f\n", spec, total_lppd))
  per_pat <- do.call(rbind, rows)
  out_csv <- file.path(run_dir, sprintf("loglik_per_patient_%s.csv", spec))
  write.csv(per_pat, out_csv, row.names = FALSE)
  cat(sprintf("  Saved %s\n", out_csv))
  list(spec = spec, total_lppd = total_lppd, per_patient = per_pat)
}

specs <- names(extension_specs)
results <- list()
for (sp in specs) {
  res <- compute_lppd_for_spec(sp)
  if (!is.null(res)) results[[sp]] <- res
}

if (length(results) >= 1) {
  summary_df <- do.call(rbind, lapply(names(results), function(sp) {
    data.frame(
      spec = sp,
      cov_id = extension_specs[[sp]]$cov_id,
      rule = extension_specs[[sp]]$rule,
      pair_key = extension_specs[[sp]]$pair_key,
      total_lppd = results[[sp]]$total_lppd,
      n_test_patients = nrow(results[[sp]]$per_patient),
      stringsAsFactors = FALSE)
  }))
  # Paper 1 anchors
  PAPER1_PRIMARY_LPPD <- -1032.43
  BASELINE_TRAIN_LPPD <- -1038.69
  summary_df$delta_vs_paper1_primary <- summary_df$total_lppd - PAPER1_PRIMARY_LPPD
  summary_df$delta_vs_baseline_train <- summary_df$total_lppd - BASELINE_TRAIN_LPPD
  write.csv(summary_df, file.path(run_dir, "loglik_summary_extension.csv"),
            row.names = FALSE)
  cat("\n=== Summary ===\n"); print(summary_df, row.names = FALSE)
}
