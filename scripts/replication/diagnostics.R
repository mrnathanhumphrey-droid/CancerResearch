user_lib <- file.path(Sys.getenv("LOCALAPPDATA"), "R/win-library/4.6")
.libPaths(c(user_lib, .libPaths()))
suppressPackageStartupMessages({
  library(posterior)
})

run_dir <- "C:/FkCancer/runs/run_paper_a_gibbs_2026-05-09"
files <- file.path(run_dir, sprintf(
  "GibbsSamplingResults_100kiters_Chain%d_WithAge_StandardizedPredictors.rda", 1:4))
stopifnot(all(file.exists(files)))

# --- Load XYC inputs to recover covariate IDs and cancer-type names ---
load("C:/FkCancer/repos/HierarchicalSS_PanCanPanOmics/XYC_V2_WithAge_StandardizedPredictors.rda")
covariates_in_model <- as.numeric(names(table(unlist(
  lapply(Covariates, function(x) as.numeric(colnames(x)))
))))
p <- length(covariates_in_model)
cancer_types <- names(Covariates)
n_cancer <- length(cancer_types)
cat(sprintf("p = %d, n_cancer = %d\n", p, n_cancer))

# --- Load each chain, extract beta_tilde / sigma2 / pi only ---
n_iter <- 100001L
burn   <- 50001L
keep   <- (burn + 1L):n_iter   # 50,000 post-burn draws
n_keep <- length(keep)

beta_tilde_arr <- array(NA_real_, dim = c(n_keep, 4, p),
                        dimnames = list(NULL, NULL, paste0("beta_tilde_", covariates_in_model)))
sigma2_arr     <- array(NA_real_, dim = c(n_keep, 4, 1),
                        dimnames = list(NULL, NULL, "sigma2"))
pi_arr         <- array(NA_real_, dim = c(n_keep, 4, p),
                        dimnames = list(NULL, NULL, paste0("pi_", covariates_in_model)))

for (j in 1:4) {
  cat(sprintf("[%s] Loading chain %d ...\n", format(Sys.time()), j))
  e <- new.env()
  load(files[j], envir = e)
  obj <- e[[ls(e)[1]]]
  cat(sprintf("    obj name: %s\n", ls(e)[1]))
  beta_tilde_arr[, j, ] <- obj$beta_tilde[keep, ]
  sigma2_arr    [, j, ] <- obj$sigma2[keep]
  pi_arr        [, j, ] <- obj$pi[keep, ]
  rm(obj, e); gc(verbose = FALSE)
}
cat(sprintf("[%s] All 4 chains loaded.\n", format(Sys.time())))

# --- summarise_draws ---
cat("\n=== beta_tilde diagnostics (post-burn 50k x 4 chains) ===\n")
draws_bt <- as_draws_array(beta_tilde_arr)
sum_bt <- summarise_draws(draws_bt,
  mean, sd,
  ~quantile(.x, probs = c(0.025, 0.5, 0.975)),
  rhat = rhat,
  ess_bulk = ess_bulk,
  ess_tail = ess_tail
)
print(sum_bt, n = Inf)

cat("\n=== sigma2 diagnostics ===\n")
draws_s2 <- as_draws_array(sigma2_arr)
sum_s2 <- summarise_draws(draws_s2,
  mean, sd,
  ~quantile(.x, probs = c(0.025, 0.5, 0.975)),
  rhat = rhat,
  ess_bulk = ess_bulk,
  ess_tail = ess_tail
)
print(sum_s2)

cat("\n=== pi (inclusion probability) diagnostics ===\n")
draws_pi <- as_draws_array(pi_arr)
sum_pi <- summarise_draws(draws_pi,
  mean, sd,
  ~quantile(.x, probs = c(0.025, 0.5, 0.975)),
  rhat = rhat,
  ess_bulk = ess_bulk,
  ess_tail = ess_tail
)
print(sum_pi, n = Inf)

# --- Convergence acceptance summary ---
cat("\n=== Convergence acceptance summary ===\n")
all_rhats <- c(sum_bt$rhat, sum_s2$rhat, sum_pi$rhat)
all_bulks <- c(sum_bt$ess_bulk, sum_s2$ess_bulk, sum_pi$ess_bulk)
all_tails <- c(sum_bt$ess_tail, sum_s2$ess_tail, sum_pi$ess_tail)
cat(sprintf("Total monitored params: %d\n", length(all_rhats)))
cat(sprintf("R-hat max: %.4f, > 1.01: %d, > 1.05: %d, > 1.10: %d\n",
            max(all_rhats, na.rm = TRUE),
            sum(all_rhats > 1.01, na.rm = TRUE),
            sum(all_rhats > 1.05, na.rm = TRUE),
            sum(all_rhats > 1.10, na.rm = TRUE)))
cat(sprintf("Bulk ESS min: %.0f, < 400: %d, < 1000: %d\n",
            min(all_bulks, na.rm = TRUE),
            sum(all_bulks < 400, na.rm = TRUE),
            sum(all_bulks < 1000, na.rm = TRUE)))
cat(sprintf("Tail ESS min: %.0f, < 400: %d, < 1000: %d\n",
            min(all_tails, na.rm = TRUE),
            sum(all_tails < 400, na.rm = TRUE),
            sum(all_tails < 1000, na.rm = TRUE)))

# --- Per-chain mean cross-check (should be similar across chains) ---
cat("\n=== Per-chain means (sigma2, mean(pi), mean(|beta_tilde|)) ===\n")
per_chain <- data.frame(
  chain = 1:4,
  sigma2_mean        = apply(sigma2_arr, 2, mean),
  beta_tilde_absmean = apply(abs(beta_tilde_arr), 2, mean),
  pi_mean            = apply(pi_arr, 2, mean)
)
print(per_chain, row.names = FALSE, digits = 6)

# --- Save full summaries ---
out_csv_bt <- file.path(run_dir, "diagnostics_beta_tilde.csv")
out_csv_pi <- file.path(run_dir, "diagnostics_pi.csv")
out_csv_s2 <- file.path(run_dir, "diagnostics_sigma2.csv")
write.csv(as.data.frame(sum_bt), out_csv_bt, row.names = FALSE)
write.csv(as.data.frame(sum_pi), out_csv_pi, row.names = FALSE)
write.csv(as.data.frame(sum_s2), out_csv_s2, row.names = FALSE)
cat(sprintf("\nSaved CSVs:\n  %s\n  %s\n  %s\n", out_csv_bt, out_csv_pi, out_csv_s2))

cat("\nDone.\n")
