# Paper 5 — CHECK 1: per-subgroup LPPD contribution per target covariate.
#
# For each of the 3 target covariates (t1, t2, t3):
#   contribution[t, subgroup] = (primary K=3 LPPD by subgroup)
#                             − (primary_drop_t<i> K=2 LPPD by subgroup)
#
# Decision rule (PRE_REGISTRATION.md §4 CHECK 1):
#   Verdict A if: no subgroup crosses |1 nat| AND binomial test on 4 directional
#                 contributions is at chance (two-sided p > 0.05).
#   Verdict B if: ≥ 1 subgroup crosses |1 nat| AND binomial p < 0.05.
#
# Reads:  results/loglik_per_patient_primary.csv
#         results/loglik_per_patient_primary_drop_t{1,2,3}.csv
#         data/cavalli_clinical_aligned.rds (for subtype→L1 mapping)
#         reference/paper5_target_covariates.csv
# Writes: results/check1_summary.csv

user_lib <- file.path(Sys.getenv("LOCALAPPDATA"), "R/win-library/4.6")
.libPaths(c(user_lib, .libPaths()))

repo_root <- "C:/Cancer Research/Paper5_Medulloblastoma_StructuralPrior"
res_dir   <- file.path(repo_root, "results")
data_dir  <- file.path(repo_root, "data")
ref_dir   <- file.path(repo_root, "reference")

clin <- readRDS(file.path(data_dir, "cavalli_clinical_aligned.rds"))
target_table <- read.csv(file.path(ref_dir, "paper5_target_covariates.csv"),
                         stringsAsFactors = FALSE)
target_cov_ids <- target_table$cov_id

primary_lppd <- read.csv(file.path(res_dir, "loglik_per_patient_primary.csv"),
                         stringsAsFactors = FALSE)

# Map L2 subtype → L1 subgroup using clinical table
sub_to_grp <- tapply(as.character(clin$subgroup), clin$subtype, function(x) x[1])

check1_rows <- list()
for (i in seq_along(target_cov_ids)) {
  spec_drop <- sprintf("primary_drop_t%d", i)
  drop_path <- file.path(res_dir, sprintf("loglik_per_patient_%s.csv", spec_drop))
  if (!file.exists(drop_path)) {
    cat(sprintf("[WARN] %s not found; skipping target %d\n", drop_path, i))
    next
  }
  drop_lppd <- read.csv(drop_path, stringsAsFactors = FALSE)

  # Align by (subtype, patient_idx) — must match by construction
  stopifnot(nrow(primary_lppd) == nrow(drop_lppd))
  stopifnot(all(primary_lppd$subtype == drop_lppd$subtype))
  stopifnot(all(primary_lppd$patient_idx == drop_lppd$patient_idx))

  # Per-patient marginal contribution of target covariate i
  delta_per_patient <- primary_lppd$lppd - drop_lppd$lppd

  # Aggregate to L1 subgroup
  l1 <- sub_to_grp[primary_lppd$subtype]
  agg <- aggregate(delta_per_patient,
                   by = list(subgroup = l1),
                   FUN = function(v) c(sum = sum(v), n = length(v)))
  contributions <- data.frame(
    target = sprintf("t%d (cov_id %s)", i, target_cov_ids[i]),
    subgroup = agg$subgroup,
    n_patients = agg$x[, "n"],
    sum_delta_nats = agg$x[, "sum"],
    abs_above_1_nat = abs(agg$x[, "sum"]) > 1,
    sign = sign(agg$x[, "sum"]),
    stringsAsFactors = FALSE
  )

  # Binomial test on the 4 directional contributions
  signs <- contributions$sign
  signs <- signs[signs != 0]
  pos <- sum(signs > 0); neg <- sum(signs < 0)
  if (pos + neg > 0) {
    bn <- binom.test(pos, pos + neg, p = 0.5)
    p_val <- bn$p.value
  } else {
    p_val <- 1.0
  }

  n_cross_1 <- sum(contributions$abs_above_1_nat)
  verdict <- if (n_cross_1 == 0 && p_val > 0.05) {
    "A (noise-tier — drop candidate)"
  } else if (n_cross_1 >= 1 && p_val < 0.05) {
    "B (load-bearing — keep)"
  } else {
    "mixed (conservative fallback to ablated spec)"
  }

  contributions$target_summary <- sprintf("n_cross_1=%d  binom_p=%.3f  verdict=%s",
                                          n_cross_1, p_val, verdict)
  check1_rows[[i]] <- contributions

  cat(sprintf("\nTarget t%d (cov_id %s): n_subgroups_crossing_1nat = %d, binom p = %.3f → %s\n",
              i, target_cov_ids[i], n_cross_1, p_val, verdict))
  print(contributions[, c("subgroup", "n_patients", "sum_delta_nats", "sign")],
        row.names = FALSE)
}

check1_df <- do.call(rbind, check1_rows)
write.csv(check1_df, file.path(res_dir, "check1_summary.csv"), row.names = FALSE)
cat(sprintf("\nSaved %s\n", file.path(res_dir, "check1_summary.csv")))
