# Paper 5 — Gibbs convergence diagnostics per spec.
#
# Loads 4 chains per spec, extracts beta_tilde (p × iters), sigma2, pi.
# Uses posterior::summarise_draws for R-hat / bulk-ESS / tail-ESS.
# Applies the halt rule: any monitored param with R-hat > 1.05 → halt that spec.
#
# Reads:  results/gibbs_<spec>/chain_[1-4].rda for spec ∈ {all 6 specs}
# Writes: results/convergence_<spec>.csv per spec
#         results/convergence_summary.csv  (max R-hat, min bulk-ESS, verdict per spec)

user_lib <- file.path(Sys.getenv("LOCALAPPDATA"), "R/win-library/4.6")
.libPaths(c(user_lib, .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Ensure posterior package is available
for (pkg in c("posterior", "backports", "abind", "checkmate", "tensorA",
              "distributional", "matrixStats")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, lib = user_lib)
  }
}
library(posterior)

repo_root <- "C:/Cancer Research/Paper5_Medulloblastoma_StructuralPrior"
res_dir   <- file.path(repo_root, "results")

specs <- c("baseline_train", "primary", "sensitivity",
           "primary_drop_t1", "primary_drop_t2", "primary_drop_t3")
BURN <- 50001L
ITERS <- 1e5L
RHAT_HALT <- 1.05

summarize_one_spec <- function(spec) {
  spec_dir <- file.path(res_dir, sprintf("gibbs_%s", spec))
  chain_files <- sort(list.files(spec_dir, pattern = "^chain_[1-4]\\.rda$",
                                 full.names = TRUE))
  if (length(chain_files) < 4) {
    cat(sprintf("[%s] only %d/4 chains landed — skip\n", spec, length(chain_files)))
    return(NULL)
  }

  # Load and stack post-burn iterations across 4 chains
  bt <- list(); s2 <- list(); pi <- list()
  for (cf in chain_files) {
    env <- new.env(); load(cf, envir = env); out <- env$out
    keep <- (BURN + 1L):ITERS
    bt_mat <- do.call(rbind, out$beta_tilde[keep])
    bt[[length(bt) + 1]] <- bt_mat
    s2[[length(s2) + 1]] <- as.numeric(out$sigma2[keep])
    pi[[length(pi) + 1]] <- do.call(rbind, out$pi[keep])
  }

  # Stack into 3D arrays [iter × chain × param]
  n_iter <- nrow(bt[[1]])
  p_bt <- ncol(bt[[1]])
  p_pi <- ncol(pi[[1]])

  bt_arr <- array(NA_real_, dim = c(n_iter, 4, p_bt))
  pi_arr <- array(NA_real_, dim = c(n_iter, 4, p_pi))
  s2_arr <- array(NA_real_, dim = c(n_iter, 4, 1))
  for (k in 1:4) {
    bt_arr[, k, ] <- bt[[k]]
    pi_arr[, k, ] <- pi[[k]]
    s2_arr[, k, 1] <- s2[[k]]
  }
  dimnames(bt_arr) <- list(NULL, NULL, paste0("beta_tilde_", seq_len(p_bt)))
  dimnames(pi_arr) <- list(NULL, NULL, paste0("pi_", seq_len(p_pi)))
  dimnames(s2_arr) <- list(NULL, NULL, "sigma2")

  bt_draws <- as_draws_array(bt_arr)
  pi_draws <- as_draws_array(pi_arr)
  s2_draws <- as_draws_array(s2_arr)

  bt_sum <- summarise_draws(bt_draws, "rhat", "ess_bulk", "ess_tail")
  pi_sum <- summarise_draws(pi_draws, "rhat", "ess_bulk", "ess_tail")
  s2_sum <- summarise_draws(s2_draws, "rhat", "ess_bulk", "ess_tail")
  full <- rbind(as.data.frame(bt_sum),
                as.data.frame(pi_sum),
                as.data.frame(s2_sum))
  out_path <- file.path(res_dir, sprintf("convergence_%s.csv", spec))
  write.csv(full, out_path, row.names = FALSE)

  max_rhat <- max(full$rhat, na.rm = TRUE)
  min_ess_bulk <- min(full$ess_bulk, na.rm = TRUE)
  min_ess_tail <- min(full$ess_tail, na.rm = TRUE)
  n_above_halt <- sum(full$rhat > RHAT_HALT, na.rm = TRUE)
  verdict <- if (n_above_halt > 0) "HALT" else "PASS"

  cat(sprintf("[%s] max R-hat = %.4f (n>%.2f = %d), min bulk-ESS = %.0f → %s\n",
              spec, max_rhat, RHAT_HALT, n_above_halt, min_ess_bulk, verdict))

  data.frame(spec = spec,
             max_rhat = max_rhat,
             n_params_above_halt = n_above_halt,
             min_ess_bulk = min_ess_bulk,
             min_ess_tail = min_ess_tail,
             verdict = verdict,
             stringsAsFactors = FALSE)
}

summary_rows <- list()
for (sp in specs) {
  r <- summarize_one_spec(sp)
  if (!is.null(r)) summary_rows[[sp]] <- r
}
summary_df <- do.call(rbind, summary_rows)
write.csv(summary_df, file.path(res_dir, "convergence_summary.csv"),
          row.names = FALSE)
cat("\nConvergence summary:\n")
print(summary_df, row.names = FALSE)
