# Paper 5 — CHECK 2: patient-level bootstrap of per-covariate marginal.
#
# For each target covariate i, the marginal contribution to held-out LPPD is:
#     marginal_i = total_LPPD(primary) − total_LPPD(primary_drop_t<i>)
# CHECK 2 wraps a B=1000 patient-level bootstrap (paired indices) around this
# difference and reports the 95% CI.
#
# Decision rule (PRE_REGISTRATION.md §4 CHECK 2):
#   Verdict A if 95% CI brackets zero → covariate i is noise-tier
#   Verdict B if 95% CI excludes zero in the direction of contribution → load-bearing
#
# Reads:  results/loglik_per_patient_primary.csv
#         results/loglik_per_patient_primary_drop_t{1,2,3}.csv
# Writes: results/check2_bootstrap.csv          (full B=1000 per target)
#         results/check2_summary.csv             (CI + verdict per target)

user_lib <- file.path(Sys.getenv("LOCALAPPDATA"), "R/win-library/4.6")
.libPaths(c(user_lib, .libPaths()))

repo_root <- "C:/Cancer Research/Paper5_Medulloblastoma_StructuralPrior"
res_dir   <- file.path(repo_root, "results")
ref_dir   <- file.path(repo_root, "reference")

target_table <- read.csv(file.path(ref_dir, "paper5_target_covariates.csv"),
                         stringsAsFactors = FALSE)
target_cov_ids <- target_table$cov_id
B <- 1000L
boot_seed <- 20260516L

primary_lppd <- read.csv(file.path(res_dir, "loglik_per_patient_primary.csv"),
                         stringsAsFactors = FALSE)
n_test <- nrow(primary_lppd)

set.seed(boot_seed)
boot_idx <- matrix(sample.int(n_test, size = n_test * B, replace = TRUE),
                   nrow = n_test, ncol = B)

summary_rows <- list()
boot_rows <- list()
for (i in seq_along(target_cov_ids)) {
  spec_drop <- sprintf("primary_drop_t%d", i)
  drop_path <- file.path(res_dir, sprintf("loglik_per_patient_%s.csv", spec_drop))
  if (!file.exists(drop_path)) {
    cat(sprintf("[WARN] %s not found; skipping target %d\n", drop_path, i))
    next
  }
  drop_lppd <- read.csv(drop_path, stringsAsFactors = FALSE)
  stopifnot(nrow(drop_lppd) == n_test)
  stopifnot(all(primary_lppd$subtype == drop_lppd$subtype))
  stopifnot(all(primary_lppd$patient_idx == drop_lppd$patient_idx))

  delta_per_patient <- primary_lppd$lppd - drop_lppd$lppd
  obs <- sum(delta_per_patient)

  boot_dist <- apply(boot_idx, 2, function(ix) sum(delta_per_patient[ix]))
  ci_lo <- quantile(boot_dist, 0.025)
  ci_md <- quantile(boot_dist, 0.5)
  ci_hi <- quantile(boot_dist, 0.975)
  excludes_zero <- (ci_lo > 0) || (ci_hi < 0)
  direction <- sign(obs)
  verdict <- if (!excludes_zero) {
    "A (noise-tier — drop candidate)"
  } else if (excludes_zero && direction > 0) {
    "B (load-bearing positive — keep)"
  } else {
    "B-neg (drop, hurts performance)"
  }

  summary_rows[[i]] <- data.frame(
    target = sprintf("t%d", i),
    cov_id = target_cov_ids[i],
    observed_marginal_nats = obs,
    ci_lo = ci_lo, ci_med = ci_md, ci_hi = ci_hi,
    excludes_zero = excludes_zero,
    verdict = verdict,
    stringsAsFactors = FALSE
  )
  boot_rows[[i]] <- data.frame(
    target = sprintf("t%d", i),
    cov_id = target_cov_ids[i],
    boot_idx = seq_len(B),
    marginal_nats = boot_dist,
    stringsAsFactors = FALSE
  )

  cat(sprintf("Target t%d (cov_id %s): obs = %.3f nats, 95%% CI [%.3f, %.3f] → %s\n",
              i, target_cov_ids[i], obs, ci_lo, ci_hi, verdict))
}

write.csv(do.call(rbind, summary_rows),
          file.path(res_dir, "check2_summary.csv"), row.names = FALSE)
write.csv(do.call(rbind, boot_rows),
          file.path(res_dir, "check2_bootstrap.csv"), row.names = FALSE)
cat(sprintf("\nSaved check2_summary.csv and check2_bootstrap.csv\n"))
