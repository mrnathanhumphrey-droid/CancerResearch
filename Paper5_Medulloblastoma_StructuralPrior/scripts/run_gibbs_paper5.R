# Paper 5 — single-chain Gibbs runner.
#
# Usage:
#   Rscript run_gibbs_paper5.R <spec> <chain_id>
#     spec ∈ {baseline_train, primary, sensitivity,
#             primary_drop_t1, primary_drop_t2, primary_drop_t3}
#     chain_id ∈ {1, 2, 3, 4}
# Seeds: 20260515 + (chain_id - 1)
#
# Outputs: gibbs_<spec>/chain_<id>.rda containing 'out' (full posterior list)
#
# Reuses the same sampler as Paper 1 (sampler_structural_prior.R, sourced
# from C:/Cancer Research/scripts/validation/).
#
# Ablation specs primary_drop_t1/t2/t3 drop one target covariate at a time
# from the primary spec, enabling the CHECK 1 + CHECK 2 marginal-contribution
# diagnostic committed in PRE_REGISTRATION.md §4.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript run_gibbs_paper5.R <spec> <chain_id>")
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

# Load training fold
load(file.path(data_dir, "train_data.rda"))
Covariates <- Covariates_train
Survival   <- Survival_train
Censored   <- Censored_train

# Resolve target_covs per spec
target_covs <- list()
if (spec == "primary") {
  load(file.path(data_dir, "target_covs_primary.rda"))
  target_covs <- target_covs_primary
} else if (spec == "sensitivity") {
  load(file.path(data_dir, "target_covs_sensitivity.rda"))
  target_covs <- target_covs_sensitivity
} else if (grepl("^primary_drop_t[1-3]$", spec)) {
  load(file.path(data_dir, "target_covs_primary.rda"))
  drop_idx <- as.integer(sub("primary_drop_t", "", spec))
  target_covs <- target_covs_primary[-drop_idx]
  cat(sprintf("Ablation spec %s: dropping target covariate %d (cov_id %s)\n",
              spec, drop_idx, target_covs_primary[[drop_idx]]$cov_id))
}
# baseline_train: target_covs = list() (empty)

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

seed_val <- 20260515L + (chain_id - 1L)
set.seed(seed_val)
cat(sprintf("=== %s Chain %d, seed %d, iters %d ===\n",
            spec, chain_id, seed_val, iters))

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
cat(sprintf("[%s] Finished\n", format(t1)))
cat(sprintf("Elapsed: %.2f minutes\n",
            as.numeric(difftime(t1, t0, units = "mins"))))

out_dir <- file.path(repo_root, "results", sprintf("gibbs_%s", spec))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_file <- file.path(out_dir, sprintf("chain_%d.rda", chain_id))
save(out, file = out_file)
cat(sprintf("Saved %s\n", out_file))
