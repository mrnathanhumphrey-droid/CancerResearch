# Paper 5 — B=1000 paired patient-level bootstrap on joint LPPD differences.
#
# Computes bootstrap 95% CI for each pairwise spec comparison among
# {baseline_train, primary, sensitivity}, mirroring Paper 1's bootstrap_ci.R.
#
# Reads:  results/loglik_per_patient_{baseline_train,primary,sensitivity}.csv
# Writes: results/bootstrap_ci_summary.csv          — pair × {obs_diff, CI, verdict}
#         results/bootstrap_loglik_difference.csv   — full B=1000 distribution per pair

user_lib <- file.path(Sys.getenv("LOCALAPPDATA"), "R/win-library/4.6")
.libPaths(c(user_lib, .libPaths()))

repo_root <- "C:/Cancer Research/Paper5_Medulloblastoma_StructuralPrior"
res_dir   <- file.path(repo_root, "results")

specs <- c("baseline_train", "primary", "sensitivity")
B <- 1000L
boot_seed <- 20260515L
DECISION_NAT <- 2.0

lppd_list <- list()
for (sp in specs) {
  p <- file.path(res_dir, sprintf("loglik_per_patient_%s.csv", sp))
  if (!file.exists(p)) {
    cat(sprintf("[WARN] %s not found; skipping spec %s\n", p, sp))
    next
  }
  lppd_list[[sp]] <- read.csv(p, stringsAsFactors = FALSE)
}

# Reference patient order (use first available spec's ordering)
ref_spec <- names(lppd_list)[1]
ref <- lppd_list[[ref_spec]]
n_test <- nrow(ref)

# Verify all loaded specs share the same patient ordering
for (sp in names(lppd_list)) {
  ll <- lppd_list[[sp]]
  stopifnot(nrow(ll) == n_test)
  stopifnot(all(ll$subtype == ref$subtype))
  stopifnot(all(ll$patient_idx == ref$patient_idx))
}

set.seed(boot_seed)
boot_idx <- matrix(sample.int(n_test, size = n_test * B, replace = TRUE),
                   nrow = n_test, ncol = B)

# All pairwise comparisons
pairs <- combn(names(lppd_list), 2, simplify = FALSE)

summary_rows <- list()
diffs_rows   <- list()
for (pr in pairs) {
  a <- pr[1]; b <- pr[2]
  d_per_patient <- lppd_list[[b]]$lppd - lppd_list[[a]]$lppd
  obs_diff <- sum(d_per_patient)

  boot_dist <- apply(boot_idx, 2, function(ix) sum(d_per_patient[ix]))
  ci_lo <- quantile(boot_dist, 0.025)
  ci_md <- quantile(boot_dist, 0.5)
  ci_hi <- quantile(boot_dist, 0.975)
  excludes_zero <- (ci_lo > 0) || (ci_hi < 0)

  verdict <- if (abs(obs_diff) <= DECISION_NAT || !excludes_zero) {
    "MATCHES"
  } else if (obs_diff > DECISION_NAT && excludes_zero && ci_lo > 0) {
    "OUTPERFORMS"
  } else if (obs_diff < -DECISION_NAT && excludes_zero && ci_hi < 0) {
    "UNDERPERFORMS"
  } else {
    "INCONCLUSIVE"
  }

  summary_rows[[paste(a, b, sep = "_vs_")]] <- data.frame(
    pair = sprintf("%s_vs_%s", b, a),
    obs_diff = obs_diff,
    ci_lo = ci_lo, ci_med = ci_md, ci_hi = ci_hi,
    excludes_zero = excludes_zero,
    verdict = verdict,
    stringsAsFactors = FALSE
  )
  diffs_rows[[paste(a, b, sep = "_vs_")]] <- data.frame(
    pair = sprintf("%s_vs_%s", b, a),
    boot_idx = seq_len(B),
    diff = boot_dist,
    stringsAsFactors = FALSE
  )

  cat(sprintf("%s: obs = %.3f, 95%% CI [%.3f, %.3f] → %s\n",
              sprintf("%s vs %s", b, a), obs_diff, ci_lo, ci_hi, verdict))
}

summary_df <- do.call(rbind, summary_rows)
diffs_df   <- do.call(rbind, diffs_rows)
write.csv(summary_df, file.path(res_dir, "bootstrap_ci_summary.csv"),
          row.names = FALSE)
write.csv(diffs_df, file.path(res_dir, "bootstrap_loglik_difference.csv"),
          row.names = FALSE)
cat(sprintf("\nSaved bootstrap_ci_summary.csv and bootstrap_loglik_difference.csv\n"))
