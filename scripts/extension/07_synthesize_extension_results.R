# Synthesize Paper 2 extension results into:
#   - operator_composition_summary.csv (flat one-row-per-spec table)
#   - extension_decision_log.txt        (verdicts + counts per pre-reg)
#   - EXTENSION_VALIDATION_RESULTS.md   (umbrella report with one filled section)
#
# Pre-reg decision rule (locked):
#   - Validates if >=2 of 18 specs IMPROVE (Bonferroni p<0.00278) AND >=1 of those
#     is RULE 2 or RULE 3.
#   - Convergence-failed specs (R-hat>1.05 anywhere monitored) are EXCLUDED from the
#     improvement count per the halt-don't-replace rule.
#   - All anchor comparisons in this synthesis use the Paper 1 primary 4-cov spec
#     (-1032.43). The vs-baseline comparisons are diagnostic, not part of the
#     pre-registered validation/falsification rule.

user_lib <- file.path(Sys.getenv("LOCALAPPDATA"), "R/win-library/4.6")
.libPaths(c(user_lib, .libPaths()))

run_dir <- "C:/FkCancer/runs/run_extension_2026-05-10"
load(file.path(run_dir, "data/extension_specs.rda"))

PAPER1_PRIMARY_LPPD <- -1032.43
BASELINE_TRAIN_LPPD <- -1038.69
BONF_ALPHA  <- 0.05 / 18  # 0.002778
DECISION_NAT <- 2.0

specs <- names(extension_specs)

conv <- read.csv(file.path(run_dir, "convergence_summary_extension.csv"),
                 stringsAsFactors = FALSE)
loglik <- read.csv(file.path(run_dir, "loglik_summary_extension.csv"),
                   stringsAsFactors = FALSE)
boot <- read.csv(file.path(run_dir, "bootstrap_extension_summary.csv"),
                 stringsAsFactors = FALSE)

# Pull the vs paper1_primary rows for the main verdict
boot_p1 <- boot[boot$anchor == "paper1_primary", ]
boot_bl <- boot[boot$anchor == "baseline_train", ]

# Build the joined per-spec summary
summary_rows <- list()
for (sp in specs) {
  cv_row <- conv[conv$spec == sp, , drop = FALSE]
  ll_row <- loglik[loglik$spec == sp, , drop = FALSE]
  bp_row <- boot_p1[boot_p1$spec == sp, , drop = FALSE]
  bb_row <- boot_bl[boot_bl$spec == sp, , drop = FALSE]
  if (!nrow(cv_row) || !nrow(ll_row) || !nrow(bp_row)) {
    cat(sprintf("[SKIP %s] missing one or more upstream rows\n", sp)); next
  }
  converged <- as.logical(cv_row$converged)
  verdict_p1 <- if (!converged) "CONVERGENCE-FAILED" else bp_row$verdict
  summary_rows[[length(summary_rows) + 1L]] <- data.frame(
    spec = sp,
    cov_id = extension_specs[[sp]]$cov_id,
    rule = extension_specs[[sp]]$rule,
    pair_key = extension_specs[[sp]]$pair_key,
    converged = converged,
    rhat_max = cv_row$rhat_max,
    bulk_ess_min = cv_row$bulk_ess_min,
    total_lppd = ll_row$total_lppd,
    delta_paper1_primary = bp_row$obs_diff,
    ci_lo_paper1 = bp_row$ci_lo,
    ci_hi_paper1 = bp_row$ci_hi,
    bonf_p_paper1 = bp_row$bonf_p,
    bonf_pass = bp_row$bonf_pass,
    delta_baseline_train = if (nrow(bb_row)) bb_row$obs_diff else NA_real_,
    ci_lo_baseline = if (nrow(bb_row)) bb_row$ci_lo else NA_real_,
    ci_hi_baseline = if (nrow(bb_row)) bb_row$ci_hi else NA_real_,
    verdict = verdict_p1,
    stringsAsFactors = FALSE)
}

out <- do.call(rbind, summary_rows)
out <- out[order(out$cov_id, out$rule, out$pair_key), ]
write.csv(out, file.path(run_dir, "operator_composition_summary.csv"),
          row.names = FALSE)
cat("Saved operator_composition_summary.csv\n")

# ---- Apply decision rule ----
n_total <- nrow(out)
n_conv_failed <- sum(out$verdict == "CONVERGENCE-FAILED")
n_improves <- sum(out$verdict == "IMPROVES")
n_nominal_only <- sum(out$verdict == "NOMINAL-ONLY")
n_matches <- sum(out$verdict == "MATCHES")
n_underperforms <- sum(out$verdict == "UNDERPERFORMS")
n_inconclusive <- sum(out$verdict == "INCONCLUSIVE")
improvers <- out[out$verdict == "IMPROVES", , drop = FALSE]
n_r2r3_improvers <- sum(improvers$rule %in% c(2L, 3L))

validates <- (n_improves >= 2L) && (n_r2r3_improvers >= 1L)

cat("\n=== Pre-registered Decision Rule (Paper 2) ===\n")
cat(sprintf("Total specs: %d\n", n_total))
cat(sprintf("  CONVERGENCE-FAILED (excluded per pre-reg): %d\n", n_conv_failed))
cat(sprintf("  IMPROVES (Bonferroni p<0.00278): %d\n", n_improves))
cat(sprintf("  NOMINAL-ONLY (passes nominal CI, fails Bonferroni): %d\n", n_nominal_only))
cat(sprintf("  MATCHES:        %d\n", n_matches))
cat(sprintf("  UNDERPERFORMS:  %d\n", n_underperforms))
cat(sprintf("  INCONCLUSIVE:   %d\n", n_inconclusive))
cat(sprintf("  Of IMPROVES, RULE 2 or RULE 3: %d\n", n_r2r3_improvers))
cat(sprintf("\nValidation condition: >=2 improve AND >=1 from R2/R3 -> %s\n",
            if (validates) "VALIDATES" else "FALSIFIES"))

# Save decision log
log_lines <- c(
  sprintf("# Paper 2 extension decision log"),
  sprintf("Generated: %s", format(Sys.time())),
  "",
  "## Pre-registered decision rule",
  sprintf("- Bonferroni alpha: %.6f (= 0.05 / 18)", BONF_ALPHA),
  sprintf("- Decision threshold: |Delta| > %.1f nats", DECISION_NAT),
  sprintf("- Validates iff: >=2 specs IMPROVE AND >=1 of those is RULE 2 or RULE 3"),
  sprintf("- Convergence halt: R-hat > 1.05 -> CONVERGENCE-FAILED, excluded from count"),
  "",
  "## Counts",
  sprintf("  Total specs:                    %d", n_total),
  sprintf("  CONVERGENCE-FAILED (excluded):  %d", n_conv_failed),
  sprintf("  IMPROVES:                       %d", n_improves),
  sprintf("  NOMINAL-ONLY:                   %d", n_nominal_only),
  sprintf("  MATCHES:                        %d", n_matches),
  sprintf("  UNDERPERFORMS:                  %d", n_underperforms),
  sprintf("  INCONCLUSIVE:                   %d", n_inconclusive),
  sprintf("  IMPROVES (RULE 2 or RULE 3):    %d", n_r2r3_improvers),
  "",
  sprintf("## Verdict: %s", if (validates) "VALIDATES" else "FALSIFIES"),
  ""
)
writeLines(log_lines, file.path(run_dir, "extension_decision_log.txt"))
cat(sprintf("\nSaved extension_decision_log.txt\n"))

# ---- Build the EXTENSION_VALIDATION_RESULTS.md report ----
fmt_row <- function(r) {
  sprintf("| %s | %s | R%d | %s | %.2f | %s | %s | %+.2f | [%+.2f, %+.2f] | %.4f | %s | %s |",
          r$spec, r$cov_id, r$rule, r$pair_key, r$rhat_max,
          if (r$converged) "Y" else "N",
          if (r$bonf_pass) "Y" else "N",
          r$delta_paper1_primary, r$ci_lo_paper1, r$ci_hi_paper1,
          r$bonf_p_paper1,
          format(r$total_lppd, nsmall = 2),
          r$verdict)
}
out_lines <- list(out)
md_lines <- c(
  "# Paper 2 — Operator-composition extension results",
  "",
  sprintf("Generated: %s", format(Sys.time())),
  "",
  sprintf("Anchor: Paper 1 primary 4-cov spec, total LPPD = %.2f.", PAPER1_PRIMARY_LPPD),
  sprintf("Baseline-on-train (diagnostic anchor) total LPPD = %.2f.", BASELINE_TRAIN_LPPD),
  sprintf("Bootstrap B = 1000, seed 20260510, Bonferroni alpha = %.6f (= 0.05 / 18).",
          BONF_ALPHA),
  "",
  "## Pre-registered decision rule",
  "",
  "Per `PRE_REGISTRATION_PAPER2.md` (committed before compute fired):",
  "",
  "- **IMPROVES**: Delta > +2 nats AND CI excludes 0 below AND Bonferroni p < 0.00278",
  "- **NOMINAL-ONLY**: passes nominal 95% CI but fails Bonferroni (NOT counted toward validation)",
  "- **MATCHES**: |Delta| <= 2 AND CI includes 0",
  "- **UNDERPERFORMS**: Delta < -2 OR CI excludes 0 below",
  "- **CONVERGENCE-FAILED**: R-hat > 1.05 anywhere monitored (halt-don't-replace, excluded)",
  "",
  "**Validation iff**: >=2 specs IMPROVE AND >=1 of those is RULE 2 or RULE 3.",
  "",
  "## Verdict",
  "",
  sprintf("**%s**", if (validates) "VALIDATES" else "FALSIFIES"),
  "",
  sprintf("- Total specs: %d", n_total),
  sprintf("- CONVERGENCE-FAILED (excluded per pre-reg): %d", n_conv_failed),
  sprintf("- IMPROVES (Bonferroni p<0.00278): %d", n_improves),
  sprintf("- of those, RULE 2 or RULE 3: %d", n_r2r3_improvers),
  sprintf("- NOMINAL-ONLY: %d", n_nominal_only),
  sprintf("- MATCHES: %d", n_matches),
  sprintf("- UNDERPERFORMS: %d", n_underperforms),
  sprintf("- INCONCLUSIVE: %d", n_inconclusive),
  "",
  "## Per-spec results table",
  "",
  "| spec | cov | rule | pair | R-hat max | conv | bonf | Delta vs P1 | 95% CI | bonf p | total LPPD | verdict |",
  "|---|---|---|---|---|---|---|---|---|---|---|---|"
)
for (i in seq_len(nrow(out))) {
  md_lines <- c(md_lines, fmt_row(out[i, , drop = FALSE]))
}
md_lines <- c(md_lines,
  "",
  "## Anchors used",
  "",
  sprintf("- Paper 1 primary 4-cov LPPD: %.2f (per-patient at run_validation_2026-05-09/loglik_per_patient_primary.csv)",
          PAPER1_PRIMARY_LPPD),
  sprintf("- baseline-on-train LPPD: %.2f (per-patient at run_validation_2026-05-09/loglik_per_patient_baseline_train.csv)",
          BASELINE_TRAIN_LPPD),
  "",
  "## Caveats",
  "",
  "- Bootstrap p-values use mid-p adjustment ((B_tail + 1) / (B + 1)). Smallest achievable p with B=1000 is ~0.002.",
  "- Anchor for all 18 specs is the Paper 1 4-covariate joint-fit total. Per-cov anchors would be different but no individual cov was fit alone in Paper 1.",
  "- Per the pre-reg, no spec is replaced or re-fit if convergence fails."
)
writeLines(md_lines, file.path(run_dir, "EXTENSION_VALIDATION_RESULTS.md"))
cat(sprintf("\nSaved EXTENSION_VALIDATION_RESULTS.md (%d lines)\n", length(md_lines)))
