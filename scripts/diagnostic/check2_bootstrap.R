# CHECK 2: Bootstrap distribution of cov 10.1's marginal contribution
# at the spec level (4-cov LPPD minus 3-cov LPPD), per patient-level resample
# of the held-out test set. B = 1000, seed = 20260510 (matches existing bootstrap).
#
# Decision rule (locked in LOCK_DIAGNOSTIC_PRE_REGISTRATION.md):
#   If bootstrap 95% CI on marginal excludes zero (positive direction): Reading B.
#   If CI brackets zero: Reading A.
#
# Inputs:  results/validation/loglik_per_patient_{baseline_train,primary,primary_3cov}.csv
# Outputs: results/diagnostic/cov_10_1_bootstrap.csv          (B=1000 rows)
#          results/diagnostic/cov_10_1_bootstrap_summary.csv  (CI + verdict)

user_lib <- file.path(Sys.getenv("LOCALAPPDATA"), "R/win-library/4.6")
if (dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))

# Resolve repo root from this script's location (scripts/diagnostic/check2_bootstrap.R)
this_script <- tryCatch(
  normalizePath(sys.frame(1)$ofile),
  error = function(e) NA_character_
)
if (is.na(this_script) || !file.exists(this_script)) {
  # Fallback: assume CWD is repo root
  repo_root <- normalizePath(getwd())
} else {
  repo_root <- normalizePath(file.path(dirname(this_script), "..", ".."))
}
cat(sprintf("Repo root: %s\n", repo_root))

results_val  <- file.path(repo_root, "results", "validation")
results_diag <- file.path(repo_root, "results", "diagnostic")
dir.create(results_diag, showWarnings = FALSE, recursive = TRUE)

B    <- 1000L
SEED <- 20260510L
set.seed(SEED)

read_pp <- function(spec) {
  fp <- file.path(results_val, sprintf("loglik_per_patient_%s.csv", spec))
  if (!file.exists(fp)) stop(sprintf("Missing %s", fp))
  read.csv(fp, stringsAsFactors = FALSE)
}

pp_4cov     <- read_pp("primary")
pp_3cov     <- read_pp("primary_3cov")
pp_baseline <- read_pp("baseline_train")

key   <- function(df) paste(df$cancer, df$patient, sep = "|")
ref   <- key(pp_baseline)
align <- function(df) df[match(ref, key(df)), ]
pp_4cov_a <- align(pp_4cov)
pp_3cov_a <- align(pp_3cov)
stopifnot(!any(is.na(pp_4cov_a$lppd)), !any(is.na(pp_3cov_a$lppd)))

n_pat <- nrow(pp_baseline)
cat(sprintf("Bootstrap CHECK 2: B=%d, n_pat=%d, seed=%d\n", B, n_pat, SEED))

# Observed marginal (cov 10.1 contribution = 4-cov LPPD - 3-cov LPPD, summed over patients)
obs_marginal <- sum(pp_4cov_a$lppd) - sum(pp_3cov_a$lppd)
cat(sprintf("Observed marginal (4-cov - 3-cov): %+.4f nats\n", obs_marginal))

# Bootstrap
diffs <- numeric(B)
for (b in 1:B) {
  idx <- sample.int(n_pat, n_pat, replace = TRUE)
  diffs[b] <- sum(pp_4cov_a$lppd[idx]) - sum(pp_3cov_a$lppd[idx])
}
q <- quantile(diffs, c(0.025, 0.50, 0.975))
ci_lo <- q[1]; ci_med <- q[2]; ci_hi <- q[3]
ci_excludes_zero_positive <- (ci_lo > 0)
ci_excludes_zero_negative <- (ci_hi < 0)

cat(sprintf("\nBootstrap 95%% CI: [%+.4f, %+.4f]  median %+.4f\n",
            ci_lo, ci_hi, ci_med))
cat(sprintf("Excludes zero (positive): %s\n", ci_excludes_zero_positive))
cat(sprintf("Excludes zero (negative): %s\n", ci_excludes_zero_negative))

if (ci_excludes_zero_positive) {
  verdict <- "READING B (cov 10.1 marginal CI excludes zero in positive direction)"
} else if (ci_excludes_zero_negative) {
  verdict <- "AMBIGUOUS (CI excludes zero in NEGATIVE direction - cov 10.1 hurts)"
} else {
  verdict <- "READING A (CI brackets zero - cov 10.1 marginal indistinguishable from zero)"
}
cat(sprintf("CHECK 2 VERDICT: %s\n", verdict))

# Save full distribution
out_df <- data.frame(b = seq_len(B), marginal_4cov_minus_3cov = diffs)
write.csv(out_df, file.path(results_diag, "cov_10_1_bootstrap.csv"), row.names = FALSE)
cat(sprintf("Saved %s (B=%d rows)\n",
            file.path(results_diag, "cov_10_1_bootstrap.csv"), B))

summary_df <- data.frame(
  metric = c("obs_marginal", "ci_lo_2.5", "ci_median", "ci_hi_97.5",
             "ci_excludes_zero_positive", "ci_excludes_zero_negative",
             "B", "n_patients", "seed", "verdict"),
  value  = c(sprintf("%.6f", obs_marginal),
             sprintf("%.6f", ci_lo),
             sprintf("%.6f", ci_med),
             sprintf("%.6f", ci_hi),
             as.character(ci_excludes_zero_positive),
             as.character(ci_excludes_zero_negative),
             as.character(B),
             as.character(n_pat),
             as.character(SEED),
             verdict),
  stringsAsFactors = FALSE)
write.csv(summary_df, file.path(results_diag, "cov_10_1_bootstrap_summary.csv"),
          row.names = FALSE)
cat(sprintf("Saved %s\n",
            file.path(results_diag, "cov_10_1_bootstrap_summary.csv")))
