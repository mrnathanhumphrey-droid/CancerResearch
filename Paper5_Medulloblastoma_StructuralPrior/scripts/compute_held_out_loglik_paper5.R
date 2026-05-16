# Paper 5 — compute held-out log predictive density (LPPD) per spec.
#
# For each spec in {baseline_train, primary, sensitivity, primary_drop_t1..t3}:
#   - Load 4 chains
#   - Thin to 500 posterior samples (every 400 post-burn iters)
#   - Compute per-patient log-likelihood under the log-normal AFT spike-and-slab
#     model on the held-out 20% test set
#   - Sum to total LPPD per spec
#
# Writes:
#   results/loglik_summary.csv          — spec, total_lppd, n_test_patients, delta_vs_baseline
#   results/loglik_per_patient_<spec>.csv  — per-patient LPPD breakdowns

user_lib <- file.path(Sys.getenv("LOCALAPPDATA"), "R/win-library/4.6")
.libPaths(c(user_lib, .libPaths()))

repo_root <- "C:/Cancer Research/Paper5_Medulloblastoma_StructuralPrior"
data_dir  <- file.path(repo_root, "data")
res_dir   <- file.path(repo_root, "results")

load(file.path(data_dir, "test_data.rda"))
subtypes <- names(Covariates_test)

BURN  <- 50001L
THIN  <- 400L
ITERS <- 1e5L
N_SAMPLES <- 500L
keep_iters <- seq(BURN, ITERS, by = THIN)
stopifnot(length(keep_iters) >= N_SAMPLES)
keep_iters <- tail(keep_iters, N_SAMPLES)

specs <- c("baseline_train", "primary", "sensitivity",
           "primary_drop_t1", "primary_drop_t2", "primary_drop_t3")

# Log-normal AFT likelihood with right-censoring:
#   uncensored:  -0.5 log(2π σ²) - (y - μ)² / (2σ²)
#   censored:    log(1 - Φ((y - μ) / σ))    [survival beyond y]
loglik_per_patient_one_draw <- function(X_test, y_test, c_test, beta_c, sigma2) {
  mu <- as.numeric(X_test %*% beta_c)
  sigma <- sqrt(sigma2)
  uncens <- c_test == 0L
  ll <- numeric(length(y_test))
  ll[uncens]  <- dnorm(y_test[uncens], mean = mu[uncens], sd = sigma, log = TRUE)
  ll[!uncens] <- pnorm(y_test[!uncens], mean = mu[!uncens], sd = sigma,
                       lower.tail = FALSE, log.p = TRUE)
  ll
}

compute_lppd_for_spec <- function(spec) {
  cat(sprintf("[%s] Loading chains for spec '%s'\n", format(Sys.time()), spec))
  spec_dir <- file.path(res_dir, sprintf("gibbs_%s", spec))
  chain_files <- list.files(spec_dir, pattern = "^chain_[1-4]\\.rda$", full.names = TRUE)
  if (length(chain_files) < 4) {
    cat(sprintf("  [WARN] only %d of 4 chains present; skipping spec\n",
                length(chain_files)))
    return(NULL)
  }

  per_pat_ll_summed <- list()    # accumulate log-mean-exp over all draws per patient
  for (s in subtypes) per_pat_ll_summed[[s]] <- rep(NA_real_, length(Survival_test[[s]]))

  # For numerical stability across draws, accumulate via log-sum-exp.
  # Total draws across 4 chains = 4 × N_SAMPLES.
  # Per-patient LPPD = log(mean_d exp(loglik_d))
  ll_acc <- list()
  for (s in subtypes) ll_acc[[s]] <- matrix(NA_real_, nrow = length(Survival_test[[s]]),
                                            ncol = 0)

  for (cf in chain_files) {
    env <- new.env(); load(cf, envir = env)
    out <- env$out
    for (s in subtypes) {
      X_test <- Covariates_test[[s]]
      y_test <- Survival_test[[s]]
      c_test <- Censored_test[[s]]
      cov_for_s <- as.numeric(colnames(X_test))

      ll_draws <- matrix(NA_real_, nrow = length(y_test), ncol = N_SAMPLES)
      for (d in seq_len(N_SAMPLES)) {
        it <- keep_iters[d]
        # Extract per-subtype β for iter `it`
        if (is.null(out$betas)) stop("Chain RDA missing $betas")
        beta_full <- out$betas[[which(names(out$betas) == s)]][[it]]
        # Restrict β to covariates present in this subtype (in correct order)
        beta_c <- beta_full[as.character(cov_for_s)]
        if (any(is.na(beta_c))) beta_c[is.na(beta_c)] <- 0
        sigma2 <- out$sigma2[[it]]
        ll_draws[, d] <- loglik_per_patient_one_draw(X_test, y_test, c_test, beta_c, sigma2)
      }
      ll_acc[[s]] <- cbind(ll_acc[[s]], ll_draws)
    }
  }

  # Per-patient LPPD = log(mean over draws of exp(loglik))
  per_patient <- list()
  for (s in subtypes) {
    M <- ll_acc[[s]]
    log_mean <- apply(M, 1, function(v) {
      mx <- max(v); mx + log(mean(exp(v - mx)))
    })
    per_patient[[s]] <- log_mean
  }

  total_lppd <- sum(unlist(per_patient))
  per_pat_df <- do.call(rbind, lapply(subtypes, function(s) {
    data.frame(subtype = s,
               patient_idx = seq_along(per_patient[[s]]),
               lppd = per_patient[[s]],
               stringsAsFactors = FALSE)
  }))
  out_path <- file.path(res_dir, sprintf("loglik_per_patient_%s.csv", spec))
  write.csv(per_pat_df, out_path, row.names = FALSE)
  cat(sprintf("  spec %s: total LPPD = %.3f, n_test = %d, saved %s\n",
              spec, total_lppd, nrow(per_pat_df), out_path))

  list(spec = spec, total_lppd = total_lppd, n_test = nrow(per_pat_df))
}

summary_rows <- list()
baseline_lppd <- NA_real_
for (sp in specs) {
  r <- compute_lppd_for_spec(sp)
  if (is.null(r)) next
  summary_rows[[sp]] <- r
  if (sp == "baseline_train") baseline_lppd <- r$total_lppd
}

summary_df <- do.call(rbind, lapply(summary_rows, function(r) {
  data.frame(spec = r$spec,
             total_lppd = r$total_lppd,
             n_test_patients = r$n_test,
             delta_vs_baseline = r$total_lppd - baseline_lppd,
             stringsAsFactors = FALSE)
}))
write.csv(summary_df, file.path(res_dir, "loglik_summary.csv"), row.names = FALSE)
cat("\nFinal loglik summary:\n")
print(summary_df, row.names = FALSE)
