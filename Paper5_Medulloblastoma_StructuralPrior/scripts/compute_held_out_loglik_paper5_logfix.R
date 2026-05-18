# Paper 5 — recompute LPPD with log() fix on the CONTAMINATED v1 chains.
#
# Purpose: disentangle the v1 +349 nat delta into:
#   (a) units bug:  v1 evaluated dnorm(y, mu, sigma) with y in raw years but
#                    mu/sigma in log-time (Lock sampler internally logs Y).
#   (b) leakage:    BIDIFAC+ was fit on n=763 including test fold.
#
# This script applies the units fix while reusing the v1 contaminated
# components + v1 contaminated chains. The delta-vs-baseline here tells us
# how much of +349 was units bug. Subtracting that from the eventual
# leakage-clean v2 delta isolates the pure leakage contribution.
#
# Reads:  data/test_data.rda  (v1 contaminated, L2 12-subtype)
#         results/gibbs_<spec>/chain_[1-4].rda  (v1 contaminated chains)
# Writes: results/loglik_summary_logfix.csv
#         results/loglik_per_patient_<spec>_logfix.csv

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
keep_iters <- seq(BURN, ITERS, by = THIN)
N_SAMPLES <- length(keep_iters)
cat(sprintf("Thinned post-burn samples per chain: %d\n", N_SAMPLES))

specs <- c("baseline_train", "primary", "sensitivity",
           "primary_drop_t1", "primary_drop_t2", "primary_drop_t3")

# log() fix: pass log(y_test) into dnorm/pnorm since mu/sigma are in log-time.
loglik_per_patient_one_draw <- function(X_test, y_test, c_test, beta_c, sigma2) {
  mu <- as.numeric(X_test %*% beta_c)
  sigma <- sqrt(sigma2)
  ly <- log(y_test)
  uncens <- c_test == 0L
  ll <- numeric(length(y_test))
  ll[uncens]  <- dnorm(ly[uncens], mean = mu[uncens], sd = sigma, log = TRUE)
  ll[!uncens] <- pnorm(ly[!uncens], mean = mu[!uncens], sd = sigma,
                       lower.tail = FALSE, log.p = TRUE)
  ll
}

compute_lppd_for_spec <- function(spec) {
  cat(sprintf("[%s] %s (logfix on v1 chains)\n", format(Sys.time()), spec))
  spec_dir <- file.path(res_dir, sprintf("gibbs_%s", spec))
  chain_files <- list.files(spec_dir, pattern = "^chain_[1-4]\\.rda$",
                            full.names = TRUE)
  if (length(chain_files) < 4) {
    cat(sprintf("  [WARN] only %d/4 chains; skipping\n", length(chain_files)))
    return(NULL)
  }

  ll_acc <- list()
  for (s in subtypes) ll_acc[[s]] <- matrix(NA_real_,
                                            nrow = length(Survival_test[[s]]),
                                            ncol = 0)

  for (cf in chain_files) {
    env <- new.env(); load(cf, envir = env)
    out <- env$out
    if (is.null(out$betas)) stop("Chain missing $betas: ", cf)
    for (s_idx in seq_along(subtypes)) {
      s <- subtypes[s_idx]
      X_test <- Covariates_test[[s]]
      y_test <- Survival_test[[s]]
      c_test <- Censored_test[[s]]
      cov_for_s <- as.numeric(colnames(X_test))

      ll_draws <- matrix(NA_real_, nrow = length(y_test), ncol = N_SAMPLES)
      for (d in seq_len(N_SAMPLES)) {
        it <- keep_iters[d]
        beta_full <- out$betas[[s_idx]][[it]]
        beta_c <- beta_full[as.character(cov_for_s)]
        if (any(is.na(beta_c))) beta_c[is.na(beta_c)] <- 0
        sigma2 <- out$sigma2[[it]]
        ll_draws[, d] <- loglik_per_patient_one_draw(X_test, y_test, c_test,
                                                    beta_c, sigma2)
      }
      ll_acc[[s]] <- cbind(ll_acc[[s]], ll_draws)
    }
  }

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
  out_path <- file.path(res_dir, sprintf("loglik_per_patient_%s_logfix.csv", spec))
  write.csv(per_pat_df, out_path, row.names = FALSE)
  cat(sprintf("  %s logfix: total LPPD = %.3f, n_test = %d → %s\n",
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
write.csv(summary_df, file.path(res_dir, "loglik_summary_logfix.csv"),
          row.names = FALSE)
cat("\n=== LOGFIX summary (v1 contaminated chains, units bug corrected) ===\n")
print(summary_df, row.names = FALSE)

# Also print v1 numbers for side-by-side comparison
v1 <- tryCatch(read.csv(file.path(res_dir, "loglik_summary.csv")), error = function(e) NULL)
if (!is.null(v1)) {
  cat("\n=== v1 (buggy units) for comparison ===\n")
  print(v1, row.names = FALSE)
  cat("\n=== Δ-vs-baseline shifts: v1 → logfix ===\n")
  cmp <- merge(v1[, c("spec", "delta_vs_baseline")],
               summary_df[, c("spec", "delta_vs_baseline")],
               by = "spec", suffixes = c(".v1", ".logfix"))
  cmp$shift <- cmp$delta_vs_baseline.logfix - cmp$delta_vs_baseline.v1
  print(cmp, row.names = FALSE)
}
