# Bootstrap CI + Bonferroni p<0.00278 per spec, vs Paper 1 anchors.
#
# Paper 2 pre-reg decision rule (locked):
#   - Improves Paper 1 anchor: Delta > +2 nats AND bootstrap p < 0.00278
#                              (Bonferroni 0.05/18) AND CI excludes 0 below
#   - Matches:      |Delta| <= 2 AND nominal CI includes 0
#   - Underperforms: Delta < -2 OR CI excludes 0 below
#   - Nominal-only: passes nominal 95% CI but fails Bonferroni (NOT counted toward validation)
#
# Anchors:
#   Paper 1 primary 4-cov spec: total LPPD -1032.43 (per-patient at
#     run_validation_2026-05-09/loglik_per_patient_primary.csv)
#   baseline-on-train:          total LPPD -1038.69 (per-patient at
#     run_validation_2026-05-09/loglik_per_patient_baseline_train.csv)
#
# Bootstrap p-value: two-sided proportion of resamples with Delta on the
# opposite side of zero from the observed (i.e., 2 * min(P(diff<=0), P(diff>=0))).
#
# Output:
#   bootstrap_extension_summary.csv   (one row per spec)
#   bootstrap_extension_diffs.csv     (B rows per spec; long format)

user_lib <- file.path(Sys.getenv("LOCALAPPDATA"), "R/win-library/4.6")
.libPaths(c(user_lib, .libPaths()))

run_dir <- "C:/FkCancer/runs/run_extension_2026-05-10"
val_dir <- "C:/FkCancer/runs/run_validation_2026-05-09"

load(file.path(run_dir, "data/extension_specs.rda"))

B <- 1000L
SEED <- 20260510L
DECISION_NAT <- 2.0
BONF_ALPHA  <- 0.05 / 18  # 0.002778
NOMINAL_ALPHA <- 0.05

set.seed(SEED)

read_pp <- function(dir, spec) {
  fp <- file.path(dir, sprintf("loglik_per_patient_%s.csv", spec))
  if (!file.exists(fp)) return(NULL)
  read.csv(fp, stringsAsFactors = FALSE)
}

# Paper 1 anchors
anchor_primary <- read_pp(val_dir, "primary")
anchor_baseline <- read_pp(val_dir, "baseline_train")
stopifnot(!is.null(anchor_primary), !is.null(anchor_baseline))

make_key <- function(df) paste(df$cancer, df$patient, sep = "|")
ref_key <- make_key(anchor_baseline)

align_to <- function(other, key_ref) {
  if (is.null(other)) return(NULL)
  k <- make_key(other)
  m <- match(key_ref, k)
  stopifnot(!any(is.na(m)))
  other[m, ]
}
anchor_primary_a  <- align_to(anchor_primary,  ref_key)
anchor_baseline_a <- anchor_baseline  # already in ref order

n_pat <- nrow(anchor_baseline_a)
cat(sprintf("Bootstrap: B=%d, n_pat=%d, decision_nat=%.1f, Bonferroni alpha=%.5f\n",
            B, n_pat, DECISION_NAT, BONF_ALPHA))

# Pre-generate bootstrap indices (shared across all 18 specs for paired comparisons)
boot_idx <- matrix(sample.int(n_pat, n_pat * B, replace = TRUE), n_pat, B)

bonf_p_two_sided <- function(diffs) {
  # Empirical two-sided p-value: smallest tail mass extended to both sides.
  # NOTE: this is a percentile-bootstrap p analog. With B=1000, smallest p
  # achievable is ~0.002 (one-sided 0/1000 -> use mid-p adjustment to avoid 0).
  p_le <- (sum(diffs <= 0) + 1) / (B + 1)  # mid-p-style adjustment
  p_ge <- (sum(diffs >= 0) + 1) / (B + 1)
  min(2 * min(p_le, p_ge), 1.0)
}

verdict_decide <- function(obs_diff, ci_lo, ci_hi, bonf_p) {
  ci_excludes_zero <- (ci_lo > 0) || (ci_hi < 0)
  bonf_pass <- bonf_p < BONF_ALPHA
  if (obs_diff > DECISION_NAT && ci_lo > 0 && bonf_pass) {
    "IMPROVES"
  } else if (obs_diff > DECISION_NAT && ci_lo > 0 && !bonf_pass) {
    "NOMINAL-ONLY"
  } else if (abs(obs_diff) <= DECISION_NAT && !ci_excludes_zero) {
    "MATCHES"
  } else if (obs_diff < -DECISION_NAT && ci_hi < 0) {
    "UNDERPERFORMS"
  } else {
    "INCONCLUSIVE"
  }
}

specs <- names(extension_specs)
summary_rows <- list()
diffs_rows <- list()

for (sp in specs) {
  pp <- read_pp(run_dir, sp)
  if (is.null(pp)) {
    cat(sprintf("[SKIP %s] per-patient LPPD missing\n", sp))
    next
  }
  pp_a <- align_to(pp, ref_key)
  cat(sprintf("\n=== %s ===\n", sp))

  for (anchor_name in c("paper1_primary", "baseline_train")) {
    anchor_a <- if (anchor_name == "paper1_primary") anchor_primary_a else anchor_baseline_a
    obs_diff <- sum(pp_a$lppd) - sum(anchor_a$lppd)
    diffs <- numeric(B)
    for (b in 1:B) {
      idx <- boot_idx[, b]
      diffs[b] <- sum(pp_a$lppd[idx]) - sum(anchor_a$lppd[idx])
    }
    q <- quantile(diffs, c(0.025, 0.5, 0.975))
    ci_lo <- q[1]; ci_med <- q[2]; ci_hi <- q[3]
    bonf_p <- bonf_p_two_sided(diffs)
    verdict <- verdict_decide(obs_diff, ci_lo, ci_hi, bonf_p)
    cat(sprintf("  vs %-15s  obs=%+8.2f  CI=[%+7.2f,%+7.2f]  bonf_p=%.4f  %s\n",
                anchor_name, obs_diff, ci_lo, ci_hi, bonf_p, verdict))

    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      spec = sp,
      cov_id = extension_specs[[sp]]$cov_id,
      rule = extension_specs[[sp]]$rule,
      pair_key = extension_specs[[sp]]$pair_key,
      anchor = anchor_name,
      obs_diff = obs_diff,
      ci_lo = ci_lo, ci_med = ci_med, ci_hi = ci_hi,
      bonf_p = bonf_p,
      ci_excludes_zero = (ci_lo > 0) || (ci_hi < 0),
      bonf_pass = bonf_p < BONF_ALPHA,
      verdict = verdict,
      stringsAsFactors = FALSE)
    diffs_rows[[length(diffs_rows) + 1L]] <- data.frame(
      spec = sp, anchor = anchor_name, b = 1:B, diff = diffs,
      stringsAsFactors = FALSE)
  }
}

if (length(summary_rows) > 0) {
  out <- do.call(rbind, summary_rows)
  write.csv(out, file.path(run_dir, "bootstrap_extension_summary.csv"),
            row.names = FALSE)
  write.csv(do.call(rbind, diffs_rows),
            file.path(run_dir, "bootstrap_extension_diffs.csv"),
            row.names = FALSE)
  cat("\nSaved bootstrap_extension_summary.csv and bootstrap_extension_diffs.csv\n")
}
