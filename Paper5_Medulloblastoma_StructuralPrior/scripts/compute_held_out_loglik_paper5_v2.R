# Paper 5 — leakage-clean held-out LPPD (v2).
#
# Differences from compute_held_out_loglik_paper5.R:
#   1. Loads leakage-clean test data per spec:
#        - L2 12-subtype: test_data_clean.rda for {baseline_train, primary,
#          primary_drop_t1/t2/t3}
#        - L1 4-subgroup: test_data_L1_clean.rda for {sensitivity}
#   2. Loads chains from gibbs_<spec>_clean/ (not contaminated gibbs_<spec>/).
#   3. Output suffixed _clean to keep contaminated vs clean LPPDs side by side.
#
# Writes:
#   results/loglik_summary_clean.csv
#   results/loglik_per_patient_<spec>_clean.csv

user_lib <- file.path(Sys.getenv("LOCALAPPDATA"), "R/win-library/4.6")
.libPaths(c(user_lib, .libPaths()))

repo_root <- "C:/Cancer Research/Paper5_Medulloblastoma_StructuralPrior"
data_dir  <- file.path(repo_root, "data")
res_dir   <- file.path(repo_root, "results")

BURN  <- 50001L
THIN  <- 400L
ITERS <- 1e5L
keep_iters <- seq(BURN, ITERS, by = THIN)
N_SAMPLES <- length(keep_iters)
cat(sprintf("Thinned post-burn samples per chain: %d (BURN=%d, THIN=%d, ITERS=%d)\n",
            N_SAMPLES, BURN, THIN, ITERS))

specs <- c("baseline_train", "primary", "sensitivity",
           "primary_drop_t1", "primary_drop_t2", "primary_drop_t3")

# Resolver: returns (Covariates_test, Survival_test, Censored_test, axis_names)
load_test_for_spec <- function(spec) {
  if (spec == "sensitivity") {
    load(file.path(data_dir, "test_data_L1_clean.rda"))
    list(Cov = Covariates_test_L1, Surv = Survival_test_L1,
         Cens = Censored_test_L1, axis = names(Covariates_test_L1))
  } else {
    load(file.path(data_dir, "test_data_clean.rda"))
    list(Cov = Covariates_test, Surv = Survival_test,
         Cens = Censored_test, axis = names(Covariates_test))
  }
}

loglik_per_patient_one_draw <- function(X_test, y_test, c_test, beta_c, sigma2) {
  mu <- as.numeric(X_test %*% beta_c)
  sigma <- sqrt(sigma2)
  uncens <- c_test == 0L
  ll <- numeric(length(y_test))
  ll[uncens]  <- dnorm(log(y_test[uncens]), mean = mu[uncens], sd = sigma, log = TRUE)
  ll[!uncens] <- pnorm(log(y_test[!uncens]), mean = mu[!uncens], sd = sigma,
                       lower.tail = FALSE, log.p = TRUE)
  # Lock sampler is log-normal AFT: y is observed in years, model fits on log(y).
  ll
}

compute_lppd_for_spec <- function(spec) {
  cat(sprintf("[%s] %s\n", format(Sys.time()), spec))
  spec_dir <- file.path(res_dir, sprintf("gibbs_%s_clean", spec))
  chain_files <- list.files(spec_dir, pattern = "^chain_[1-4]\\.rda$",
                            full.names = TRUE)
  if (length(chain_files) < 4) {
    cat(sprintf("  [WARN] only %d of 4 chains; skipping\n", length(chain_files)))
    return(NULL)
  }

  td <- load_test_for_spec(spec)
  axis_names <- td$axis
  cat(sprintf("  Test axis: %d cancer-types (%s)\n",
              length(axis_names), paste(axis_names, collapse = ", ")))

  ll_acc <- list()
  for (a in axis_names) ll_acc[[a]] <- matrix(NA_real_,
                                              nrow = length(td$Surv[[a]]),
                                              ncol = 0)

  for (cf in chain_files) {
    env <- new.env(); load(cf, envir = env)
    out <- env$out
    if (is.null(out$betas)) stop("Chain missing $betas: ", cf)

    for (a_idx in seq_along(axis_names)) {
      a <- axis_names[a_idx]
      X_test <- td$Cov[[a]]
      y_test <- td$Surv[[a]]
      c_test <- td$Cens[[a]]
      cov_for_a <- as.numeric(colnames(X_test))

      ll_draws <- matrix(NA_real_, nrow = length(y_test), ncol = N_SAMPLES)
      for (d in seq_len(N_SAMPLES)) {
        it <- keep_iters[d]
        beta_full <- out$betas[[a_idx]][[it]]
        beta_c <- beta_full[as.character(cov_for_a)]
        if (any(is.na(beta_c))) beta_c[is.na(beta_c)] <- 0
        sigma2 <- out$sigma2[[it]]
        ll_draws[, d] <- loglik_per_patient_one_draw(X_test, y_test, c_test,
                                                    beta_c, sigma2)
      }
      ll_acc[[a]] <- cbind(ll_acc[[a]], ll_draws)
    }
  }

  per_patient <- list()
  for (a in axis_names) {
    M <- ll_acc[[a]]
    log_mean <- apply(M, 1, function(v) {
      mx <- max(v); mx + log(mean(exp(v - mx)))
    })
    per_patient[[a]] <- log_mean
  }

  total_lppd <- sum(unlist(per_patient))
  per_pat_df <- do.call(rbind, lapply(axis_names, function(a) {
    data.frame(unit = a,
               patient_idx = seq_along(per_patient[[a]]),
               lppd = per_patient[[a]],
               stringsAsFactors = FALSE)
  }))
  out_path <- file.path(res_dir, sprintf("loglik_per_patient_%s_clean.csv", spec))
  write.csv(per_pat_df, out_path, row.names = FALSE)
  cat(sprintf("  %s: total LPPD = %.3f, n_test = %d → %s\n",
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
write.csv(summary_df, file.path(res_dir, "loglik_summary_clean.csv"),
          row.names = FALSE)
cat("\nFinal LEAKAGE-CLEAN loglik summary:\n")
print(summary_df, row.names = FALSE)
