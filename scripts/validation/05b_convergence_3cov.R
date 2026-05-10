user_lib <- file.path(Sys.getenv("LOCALAPPDATA"), "R/win-library/4.6")
.libPaths(c(user_lib, .libPaths()))
suppressPackageStartupMessages(library(posterior))

run_dir <- "C:/FkCancer/runs/run_validation_2026-05-09"
load(file.path(run_dir, "data/train_data.rda"))
covariates_in_model <- as.numeric(names(table(unlist(
  lapply(Covariates_train, function(c) as.numeric(colnames(c)))))))
p <- length(covariates_in_model)

n_iter <- 100001L; burn <- 50001L
keep <- (burn + 1L):n_iter; n_keep <- length(keep)

spec_dir <- file.path(run_dir, "gibbs_primary_3cov")
files <- file.path(spec_dir, sprintf("chain_%d.rda", 1:4))
stopifnot(all(file.exists(files)))

beta_tilde_arr <- array(NA_real_, dim = c(n_keep, 4, p))
sigma2_arr <- array(NA_real_, dim = c(n_keep, 4, 1))
pi_arr <- array(NA_real_, dim = c(n_keep, 4, p))
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
cat(sprintf("=== primary_3cov convergence ===\n"))
cat(sprintf("  R-hat max:   %.4f  (>1.01: %d, >1.05: %d)\n",
            max(rh, na.rm=TRUE),
            sum(rh > 1.01, na.rm=TRUE), sum(rh > 1.05, na.rm=TRUE)))
cat(sprintf("  Bulk ESS min: %.0f  (<400: %d)\n",
            min(bk, na.rm=TRUE), sum(bk < 400, na.rm=TRUE)))
cat(sprintf("  Tail ESS min: %.0f\n", min(tk, na.rm=TRUE)))
if (sum(rh > 1.05, na.rm=TRUE) > 0)
  cat("*** WARNING: R-hat > 1.05 — halt rule triggered ***\n")
