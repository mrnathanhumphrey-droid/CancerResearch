# Paper 5 — leakage-clean single-chain Gibbs runner (v2).
#
# Differences from run_gibbs_paper5.R:
#   1. Loads leakage-clean training data (train_data_clean.rda for L2 specs;
#      train_data_L1_clean.rda for sensitivity spec).
#   2. Sensitivity now uses 4 L1 subgroups as the cancer-type axis WITH
#      identity hierarchy (vs primary's 12 L2 with L1 pooling). Fixes Flag 1.
#   3. Output dir suffixed _clean to avoid clobbering contaminated chains.
#
# Usage:
#   Rscript run_gibbs_paper5_v2.R <spec> <chain_id>
#     spec ∈ {baseline_train, primary, sensitivity,
#             primary_drop_t1, primary_drop_t2, primary_drop_t3}
#     chain_id ∈ {1..4}
# Seeds: 20260518 + (chain_id - 1)   [+3 vs v1 to mark distinct run]

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript run_gibbs_paper5_v2.R <spec> <chain_id>")
spec <- args[[1]]
chain_id <- as.integer(args[[2]])
allowed_specs <- c("baseline_train", "primary", "sensitivity",
                   "primary_drop_t1", "primary_drop_t2", "primary_drop_t3")
stopifnot(spec %in% allowed_specs)

user_lib <- file.path(Sys.getenv("LOCALAPPDATA"), "R/win-library/4.6")
.libPaths(c(user_lib, .libPaths()))

repo_root <- "C:/Cancer Research/Paper5_Medulloblastoma_StructuralPrior"
data_dir  <- file.path(repo_root, "data")
sampler_path <- "C:/Cancer Research/scripts/validation/sampler_structural_prior.R"
source(sampler_path)

# Resolve training data + target_covs per spec
if (spec == "sensitivity") {
  load(file.path(data_dir, "train_data_L1_clean.rda"))
  Covariates <- Covariates_train_L1
  Survival   <- Survival_train_L1
  Censored   <- Censored_train_L1
  load(file.path(data_dir, "target_covs_sensitivity_clean.rda"))
  target_covs <- target_covs_sensitivity
  cat(sprintf("Spec '%s': L1 4-subgroup, identity hierarchy, %d cancer-types\n",
              spec, length(Covariates)))
} else {
  # L2 12-subtype train data for all other specs
  load(file.path(data_dir, "train_data_clean.rda"))
  Covariates <- Covariates_train
  Survival   <- Survival_train
  Censored   <- Censored_train

  if (spec == "primary") {
    load(file.path(data_dir, "target_covs_primary_clean.rda"))
    target_covs <- target_covs_primary
  } else if (grepl("^primary_drop_t[1-3]$", spec)) {
    load(file.path(data_dir, "target_covs_primary_clean.rda"))
    drop_idx <- as.integer(sub("primary_drop_t", "", spec))
    target_covs <- target_covs_primary[-drop_idx]
    cat(sprintf("Ablation %s: dropping cov %s\n",
                spec, target_covs_primary[[drop_idx]]$cov_id))
  } else {
    target_covs <- list()   # baseline_train
  }
  cat(sprintf("Spec '%s': L2 12-subtype, %d cancer-types\n",
              spec, length(Covariates)))
}

covariates_by_cancer <- lapply(Covariates, function(c) as.numeric(colnames(c)))
covariates_in_model  <- as.numeric(names(table(unlist(covariates_by_cancer))))
p <- length(covariates_in_model)
n_cancer <- length(Covariates)
iters <- 1e5

priors <- list(betatilde_priorvar_intercept = 10^2,
               betatilde_priorvar_coefficient = 1,
               lambda2_priorshape_intercept = 1,
               lambda2_priorrate_intercept = 1,
               lambda2_priorshape_coefficient = 5,
               lambda2_priorrate_coefficient = 1,
               sigma2_priorshape = .01,
               sigma2_priorrate = .01,
               spike_priorvar = 1/10000)

seed_val <- 20260518L + (chain_id - 1L)
set.seed(seed_val)
cat(sprintf("=== CLEAN %s chain %d, seed %d, iters %d, p=%d, n_cancer=%d ===\n",
            spec, chain_id, seed_val, iters, p, n_cancer))

beta_tilde_start <- rnorm(p, 0, 1)
sigma2_start <- 1/rgamma(1, 1, 1)
lambda2_start <- 1/rgamma(p, 1, 1)
pi_start <- rep(0.5, p)
gamma_start <- lapply(1:n_cancer, function(type) {
  ind <- rbinom(length(covariates_by_cancer[[type]]), 1, 1)
  ind[1] <- 1
  names(ind) <- covariates_by_cancer[[type]]
  list(ind)
})
starting_values <- list(sigma2_start = sigma2_start,
                        beta_tilde_start = beta_tilde_start,
                        lambda2_start = lambda2_start,
                        gamma_start = gamma_start,
                        pi_start = pi_start)

cat(sprintf("[%s] Starting Gibbs\n", format(Sys.time())))
t0 <- Sys.time()
out <- HierarchicalLogNormalSpikeSlab_StructuralPrior(
  Covariates, Survival, Censored,
  starting_values, iters, covariates_in_model,
  covariates_by_cancer, priors,
  pi_generation = "shared_across_cancers",
  target_covs = target_covs,
  progress = "model_comparison"
)
t1 <- Sys.time()
cat(sprintf("[%s] Finished in %.2f min\n",
            format(t1), as.numeric(difftime(t1, t0, units = "mins"))))

out_dir <- file.path(repo_root, "results", sprintf("gibbs_%s_clean", spec))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_file <- file.path(out_dir, sprintf("chain_%d.rda", chain_id))
save(out, file = out_file)
cat(sprintf("Saved %s\n", out_file))
