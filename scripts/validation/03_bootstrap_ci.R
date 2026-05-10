# Bootstrap CI on held-out log-likelihood differences.
#
# For each pair (primary vs baseline_train, secondary vs baseline_train, primary
# vs secondary): resample patients with replacement from the test set, compute
# the difference in summed LPPD per bootstrap, B=1000 reps. Report 2.5/50/97.5
# percentiles and decision verdict per pre-registered ±2-nat rule.

user_lib <- file.path(Sys.getenv("LOCALAPPDATA"), "R/win-library/4.6")
.libPaths(c(user_lib, .libPaths()))

run_dir <- "C:/FkCancer/runs/run_validation_2026-05-09"

B <- 1000L
SEED <- 20260510L
DECISION_NAT <- 2.0

set.seed(SEED)

read_pp <- function(spec) {
  fp <- file.path(run_dir, sprintf("loglik_per_patient_%s.csv", spec))
  if (!file.exists(fp)) return(NULL)
  read.csv(fp, stringsAsFactors = FALSE)
}

baseline <- read_pp("baseline_train")
primary  <- read_pp("primary")
secondary <- read_pp("secondary")

if (is.null(baseline)) stop("baseline_train per-patient LPPD missing")

# Stack rows per patient (keyed by cancer + patient)
make_key <- function(df) paste(df$cancer, df$patient, sep = "|")
b_key <- make_key(baseline)

# Align by key
align_to <- function(other, key_ref) {
  if (is.null(other)) return(NULL)
  k <- make_key(other)
  m <- match(key_ref, k)
  stopifnot(!any(is.na(m)))
  other[m, ]
}
primary_aligned   <- align_to(primary,   b_key)
secondary_aligned <- align_to(secondary, b_key)

n_pat <- nrow(baseline)
cat(sprintf("Bootstrap: B=%d, n_pat=%d, decision threshold = %.1f nats\n",
            B, n_pat, DECISION_NAT))

bootstrap_diff <- function(lppd_a, lppd_b, B = 1000L) {
  diffs <- numeric(B)
  for (b in 1:B) {
    idx <- sample.int(n_pat, n_pat, replace = TRUE)
    diffs[b] <- sum(lppd_a[idx]) - sum(lppd_b[idx])
  }
  diffs
}

results <- list()
all_diff_rows <- list()
total <- function(x) sum(x$lppd)

for (pair in list(
  list(name = "primary_vs_baseline",   a = primary_aligned,   b = baseline),
  list(name = "secondary_vs_baseline", a = secondary_aligned, b = baseline),
  list(name = "primary_vs_secondary",  a = primary_aligned,   b = secondary_aligned)
)) {
  if (is.null(pair$a) || is.null(pair$b)) next
  obs_diff <- total(pair$a) - total(pair$b)
  diffs <- bootstrap_diff(pair$a$lppd, pair$b$lppd, B)
  q <- quantile(diffs, c(0.025, 0.5, 0.975))
  ci_lo <- q[1]; ci_med <- q[2]; ci_hi <- q[3]
  ci_excludes_zero <- (ci_lo > 0) || (ci_hi < 0)
  if (abs(obs_diff) <= DECISION_NAT && !ci_excludes_zero) {
    verdict <- "MATCHES"
  } else if (obs_diff > DECISION_NAT && ci_lo > 0) {
    verdict <- "OUTPERFORMS"
  } else if (obs_diff < -DECISION_NAT && ci_hi < 0) {
    verdict <- "UNDERPERFORMS"
  } else {
    verdict <- "INCONCLUSIVE"  # in the buffer / mixed signal
  }
  cat(sprintf("\n=== %s ===\n", pair$name))
  cat(sprintf("  observed Δ = %+9.2f nats\n", obs_diff))
  cat(sprintf("  bootstrap CI 95%% = [%+9.2f, %+9.2f]\n", ci_lo, ci_hi))
  cat(sprintf("  bootstrap median = %+9.2f\n", ci_med))
  cat(sprintf("  CI excludes zero: %s\n", ci_excludes_zero))
  cat(sprintf("  VERDICT: %s\n", verdict))
  results[[pair$name]] <- list(
    pair = pair$name, obs_diff = obs_diff,
    ci_lo = ci_lo, ci_med = ci_med, ci_hi = ci_hi,
    excludes_zero = ci_excludes_zero, verdict = verdict)
  all_diff_rows[[length(all_diff_rows) + 1L]] <- data.frame(
    pair = pair$name, b = 1:B, diff = diffs, stringsAsFactors = FALSE)
}

if (length(results) > 0) {
  out_summary <- do.call(rbind, lapply(results, as.data.frame))
  write.csv(out_summary, file.path(run_dir, "bootstrap_ci_summary.csv"),
            row.names = FALSE)
  if (length(all_diff_rows) > 0) {
    write.csv(do.call(rbind, all_diff_rows),
              file.path(run_dir, "bootstrap_loglik_difference.csv"),
              row.names = FALSE)
  }
  cat("\nSaved bootstrap_ci_summary.csv and bootstrap_loglik_difference.csv\n")
}
