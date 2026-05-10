# Convergence diagnostics for each spec (baseline_train, primary, secondary).
# Halt-eligible per the brief: if any monitored param has R-hat > 1.05, report
# and halt downstream interpretation.

user_lib <- file.path(Sys.getenv("LOCALAPPDATA"), "R/win-library/4.6")
.libPaths(c(user_lib, .libPaths()))
library(posterior)

run_dir <- "C:/FkCancer/runs/run_validation_2026-05-09"
load(file.path(run_dir, "data/train_data.rda"))
covariates_in_model <- as.numeric(names(table(unlist(
  lapply(Covariates_train, function(c) as.numeric(colnames(c)))))))
p <- length(covariates_in_model)

n_iter <- 100001L
burn   <- 50001L
keep   <- (burn + 1L):n_iter
n_keep <- length(keep)

specs <- c("baseline_train", "primary", "secondary")
all_summary <- list()

for (sp in specs) {
  spec_dir <- file.path(run_dir, sprintf("gibbs_%s", sp))
  files <- file.path(spec_dir, sprintf("chain_%d.rda", 1:4))
  if (!all(file.exists(files))) {
    cat(sprintf("[SKIP %s] not all chains present\n", sp))
    next
  }
  cat(sprintf("\n=== %s ===\n", sp))
  beta_tilde_arr <- array(NA_real_, dim = c(n_keep, 4, p),
                          dimnames = list(NULL, NULL,
                                          paste0("beta_tilde_", covariates_in_model)))
  sigma2_arr <- array(NA_real_, dim = c(n_keep, 4, 1),
                      dimnames = list(NULL, NULL, "sigma2"))
  pi_arr <- array(NA_real_, dim = c(n_keep, 4, p),
                  dimnames = list(NULL, NULL, paste0("pi_", covariates_in_model)))
  for (j in 1:4) {
    e <- new.env(); load(files[j], envir = e); o <- e$out
    beta_tilde_arr[, j, ] <- o$beta_tilde[keep, ]
    sigma2_arr[, j, ]     <- o$sigma2[keep]
    pi_arr[, j, ]         <- o$pi[keep, ]
    rm(o, e); gc(verbose = FALSE)
  }
  s_bt <- summarise_draws(as_draws_array(beta_tilde_arr),
                           rhat = rhat, ess_bulk = ess_bulk, ess_tail = ess_tail)
  s_s2 <- summarise_draws(as_draws_array(sigma2_arr),
                           rhat = rhat, ess_bulk = ess_bulk, ess_tail = ess_tail)
  s_pi <- summarise_draws(as_draws_array(pi_arr),
                           rhat = rhat, ess_bulk = ess_bulk, ess_tail = ess_tail)
  rh <- c(s_bt$rhat, s_s2$rhat, s_pi$rhat)
  bk <- c(s_bt$ess_bulk, s_s2$ess_bulk, s_pi$ess_bulk)
  tk <- c(s_bt$ess_tail, s_s2$ess_tail, s_pi$ess_tail)
  rh_max <- max(rh, na.rm = TRUE)
  rh_above_105 <- sum(rh > 1.05, na.rm = TRUE)
  rh_above_101 <- sum(rh > 1.01, na.rm = TRUE)
  bk_min <- min(bk, na.rm = TRUE)
  bk_below_400 <- sum(bk < 400, na.rm = TRUE)
  cat(sprintf("  R-hat max:   %.4f  (>1.01: %d, >1.05: %d)\n",
              rh_max, rh_above_101, rh_above_105))
  cat(sprintf("  Bulk ESS min: %.0f  (<400: %d)\n", bk_min, bk_below_400))
  cat(sprintf("  Tail ESS min: %.0f\n", min(tk, na.rm = TRUE)))
  if (rh_above_105 > 0) {
    cat(sprintf("  *** WARNING: %d params have R-hat > 1.05 — bad sampling ***\n",
                rh_above_105))
  }
  all_summary[[sp]] <- data.frame(
    spec = sp,
    rhat_max = rh_max,
    rhat_above_101 = rh_above_101,
    rhat_above_105 = rh_above_105,
    bulk_ess_min = bk_min,
    bulk_ess_below_400 = bk_below_400,
    tail_ess_min = min(tk, na.rm = TRUE),
    stringsAsFactors = FALSE)
  # Save full per-param diagnostics
  write.csv(as.data.frame(s_bt),
            file.path(run_dir, sprintf("convergence_beta_tilde_%s.csv", sp)),
            row.names = FALSE)
  write.csv(as.data.frame(s_s2),
            file.path(run_dir, sprintf("convergence_sigma2_%s.csv", sp)),
            row.names = FALSE)
  write.csv(as.data.frame(s_pi),
            file.path(run_dir, sprintf("convergence_pi_%s.csv", sp)),
            row.names = FALSE)
}

if (length(all_summary) > 0) {
  out <- do.call(rbind, all_summary)
  write.csv(out, file.path(run_dir, "convergence_summary.csv"), row.names = FALSE)
  cat("\n=== Summary ===\n"); print(out, row.names = FALSE)
}
