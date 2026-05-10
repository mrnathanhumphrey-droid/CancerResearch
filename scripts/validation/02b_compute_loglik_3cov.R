# Held-out LPPD for the 3-covariate primary spec (Reading A vs B test).
# Reuses the same posterior-sampling thinning + scoring as 02_compute_held_out_loglik.R
# but processes only the gibbs_primary_3cov chains.

user_lib <- file.path(Sys.getenv("LOCALAPPDATA"), "R/win-library/4.6")
.libPaths(c(user_lib, .libPaths()))

run_dir <- "C:/FkCancer/runs/run_validation_2026-05-09"
load(file.path(run_dir, "data/test_data.rda"))
load(file.path(run_dir, "data/train_data.rda"))
covariates_by_cancer <- lapply(Covariates_train, function(c) as.numeric(colnames(c)))
cancer_types <- names(Covariates_train)
n_cancer <- length(cancer_types)

BURN <- 50001L
THIN <- 400L

logsumexp <- function(x) { m <- max(x); m + log(sum(exp(x - m))) }

compute_lppd_for_spec <- function(spec) {
  out_dir <- file.path(run_dir, sprintf("gibbs_%s", spec))
  chain_files <- file.path(out_dir, sprintf("chain_%d.rda", 1:4))
  if (!all(file.exists(chain_files))) {
    cat(sprintf("[SKIP %s] chain files missing\n", spec)); return(NULL)
  }
  cat(sprintf("\n=== %s ===\n", spec))

  betas_per_cancer <- vector("list", n_cancer)
  for (c in seq_len(n_cancer)) betas_per_cancer[[c]] <- list()
  sigma2_samples <- numeric(0)

  for (j in 1:4) {
    cat(sprintf("  loading chain %d ...\n", j))
    e <- new.env(); load(chain_files[j], envir = e); o <- e$out
    keep <- seq(BURN, 100000L, by = THIN)
    sigma2_samples <- c(sigma2_samples, o$sigma2[keep])
    for (c in seq_len(n_cancer)) {
      mat <- do.call(rbind, o$betas[[c]][keep])
      betas_per_cancer[[c]] <- if (length(betas_per_cancer[[c]]) == 0) mat
                                else rbind(betas_per_cancer[[c]], mat)
    }
    rm(o, e); gc(verbose = FALSE)
  }
  S <- length(sigma2_samples)
  cat(sprintf("  total posterior samples: %d\n", S))

  rows <- list()
  total_lppd <- 0
  for (c in seq_len(n_cancer)) {
    Xte <- as.matrix(Covariates_test[[c]])
    Yte <- Survival_test[[c]]
    Cte <- Censored_test[[c]]
    n_test <- length(Yte); if (n_test == 0) next
    Bm <- betas_per_cancer[[c]]
    stopifnot(ncol(Xte) == ncol(Bm))
    mu_mat <- Xte %*% t(Bm)
    sd_mat <- matrix(sqrt(sigma2_samples), n_test, S, byrow = TRUE)
    censored_flag <- is.na(Yte)
    log_dens_mat <- matrix(NA_real_, n_test, S)
    for (i in seq_len(n_test)) {
      if (!censored_flag[i]) {
        log_dens_mat[i,] <- dnorm(log(Yte[i]), mean = mu_mat[i,], sd = sd_mat[i,], log = TRUE)
      } else {
        z <- (mu_mat[i,] - log(Cte[i])) / sd_mat[i,]
        log_dens_mat[i,] <- pnorm(z, log.p = TRUE)
      }
    }
    lppd_per_patient <- apply(log_dens_mat, 1, logsumexp) - log(S)
    cancer_total <- sum(lppd_per_patient)
    total_lppd <- total_lppd + cancer_total
    cat(sprintf("    %-5s n_test=%3d censored=%2d  lppd=%9.2f\n",
                cancer_types[c], n_test, sum(censored_flag), cancer_total))
    for (i in seq_len(n_test)) {
      rows[[length(rows) + 1L]] <- data.frame(
        cancer = cancer_types[c], patient = i,
        censored = censored_flag[i], lppd = lppd_per_patient[i],
        stringsAsFactors = FALSE)
    }
  }
  cat(sprintf("  TOTAL test LPPD (%s): %.2f\n", spec, total_lppd))
  per_pat <- do.call(rbind, rows)
  out_csv <- file.path(run_dir, sprintf("loglik_per_patient_%s.csv", spec))
  write.csv(per_pat, out_csv, row.names = FALSE)
  cat(sprintf("  Saved %s\n", out_csv))
  list(spec = spec, total_lppd = total_lppd, per_patient = per_pat)
}

# Compute for the 3-cov spec
res_3cov <- compute_lppd_for_spec("primary_3cov")
if (is.null(res_3cov)) stop("primary_3cov chains not found.")

# Read existing primary + baseline_train per-patient LPPDs for comparison
read_pp <- function(spec) {
  fp <- file.path(run_dir, sprintf("loglik_per_patient_%s.csv", spec))
  if (!file.exists(fp)) return(NULL)
  read.csv(fp, stringsAsFactors = FALSE)
}
pp_4cov     <- read_pp("primary")
pp_baseline <- read_pp("baseline_train")
pp_3cov     <- res_3cov$per_patient

stopifnot(nrow(pp_3cov) == nrow(pp_4cov))
stopifnot(nrow(pp_3cov) == nrow(pp_baseline))

key <- function(df) paste(df$cancer, df$patient, sep = "|")
align_to <- function(other, ref_key) other[match(ref_key, key(other)), ]
ref_key <- key(pp_baseline)
pp_3cov_a <- align_to(pp_3cov, ref_key)
pp_4cov_a <- align_to(pp_4cov, ref_key)

total_baseline <- sum(pp_baseline$lppd)
total_3cov     <- sum(pp_3cov_a$lppd)
total_4cov     <- sum(pp_4cov_a$lppd)
cat(sprintf("\n=== Total LPPD ===\n  baseline: %.2f\n  4-cov primary: %.2f (Δ vs baseline = %+.2f)\n  3-cov primary: %.2f (Δ vs baseline = %+.2f)\n  Δ (3-cov - 4-cov) = %+.2f\n",
            total_baseline, total_4cov, total_4cov-total_baseline,
            total_3cov, total_3cov-total_baseline,
            total_3cov - total_4cov))

# Bootstrap CI on (3-cov - 4-cov) and (3-cov - baseline)
B <- 1000L
set.seed(20260510L)
n_pat <- nrow(pp_baseline)
diffs_3cov_vs_4cov <- numeric(B)
diffs_3cov_vs_base <- numeric(B)
diffs_4cov_vs_base <- numeric(B)
for (b in 1:B) {
  idx <- sample.int(n_pat, n_pat, replace = TRUE)
  diffs_3cov_vs_4cov[b] <- sum(pp_3cov_a$lppd[idx]) - sum(pp_4cov_a$lppd[idx])
  diffs_3cov_vs_base[b] <- sum(pp_3cov_a$lppd[idx]) - sum(pp_baseline$lppd[idx])
  diffs_4cov_vs_base[b] <- sum(pp_4cov_a$lppd[idx]) - sum(pp_baseline$lppd[idx])
}
fmt_ci <- function(v) {
  q <- quantile(v, c(0.025, 0.5, 0.975))
  sprintf("[%+.2f, %+.2f] (median %+.2f)", q[1], q[3], q[2])
}
cat(sprintf("\n=== Bootstrap 95%% CI (B=%d) ===\n", B))
cat(sprintf("  3-cov primary vs 4-cov primary: %s, exclude 0? %s\n",
            fmt_ci(diffs_3cov_vs_4cov),
            (quantile(diffs_3cov_vs_4cov, .025) > 0) ||
            (quantile(diffs_3cov_vs_4cov, .975) < 0)))
cat(sprintf("  3-cov primary vs baseline:      %s, exclude 0? %s\n",
            fmt_ci(diffs_3cov_vs_base),
            (quantile(diffs_3cov_vs_base, .025) > 0) ||
            (quantile(diffs_3cov_vs_base, .975) < 0)))
cat(sprintf("  4-cov primary vs baseline:      %s, exclude 0? %s\n",
            fmt_ci(diffs_4cov_vs_base),
            (quantile(diffs_4cov_vs_base, .025) > 0) ||
            (quantile(diffs_4cov_vs_base, .975) < 0)))

# Reading verdict
delta_3v4 <- total_3cov - total_4cov
ci_3v4 <- quantile(diffs_3cov_vs_4cov, c(0.025, 0.975))
verdict <- if (abs(delta_3v4) < 2 && ci_3v4[1] < 0 && ci_3v4[2] > 0) {
  "READING A: cov 10.1 not contributing meaningful signal (3-cov ≈ 4-cov within bootstrap noise)"
} else if (delta_3v4 < -2 && ci_3v4[2] < 0) {
  "READING B: cov 10.1 contributing real signal (3-cov significantly worse than 4-cov)"
} else {
  sprintf("AMBIGUOUS: Δ=%.2f, CI=[%.2f, %.2f]", delta_3v4, ci_3v4[1], ci_3v4[2])
}
cat(sprintf("\n=== VERDICT ===\n%s\n", verdict))

# Save summary
out <- data.frame(
  spec = c("baseline_train", "primary_4cov", "primary_3cov"),
  total_lppd = c(total_baseline, total_4cov, total_3cov),
  delta_vs_baseline = c(0, total_4cov - total_baseline, total_3cov - total_baseline),
  delta_vs_4cov_primary = c(NA, 0, delta_3v4),
  stringsAsFactors = FALSE)
write.csv(out, file.path(run_dir, "loglik_3cov_summary.csv"), row.names = FALSE)
cat(sprintf("Saved %s\n", file.path(run_dir, "loglik_3cov_summary.csv")))

# Per-cancer breakdown comparison
cat("\n=== Per-cancer Δ_3cov vs Δ_4cov (primary) ===\n")
per_cancer <- data.frame(cancer = cancer_types, stringsAsFactors = FALSE)
per_cancer$baseline <- sapply(cancer_types, function(c) sum(pp_baseline$lppd[pp_baseline$cancer == c]))
per_cancer$primary_4cov <- sapply(cancer_types, function(c) sum(pp_4cov_a$lppd[pp_4cov_a$cancer == c]))
per_cancer$primary_3cov <- sapply(cancer_types, function(c) sum(pp_3cov_a$lppd[pp_3cov_a$cancer == c]))
per_cancer$delta_4cov <- per_cancer$primary_4cov - per_cancer$baseline
per_cancer$delta_3cov <- per_cancer$primary_3cov - per_cancer$baseline
per_cancer$delta_3cov_minus_4cov <- per_cancer$primary_3cov - per_cancer$primary_4cov
per_cancer <- per_cancer[order(per_cancer$delta_3cov_minus_4cov), ]
print(per_cancer, row.names = FALSE, digits = 4)
write.csv(per_cancer, file.path(run_dir, "per_cancer_loglik_3cov_comparison.csv"),
          row.names = FALSE)
cat(sprintf("Saved %s\n", file.path(run_dir, "per_cancer_loglik_3cov_comparison.csv")))
