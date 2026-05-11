# Single-chain runner for one of the 18 extension specifications.
#
# Usage:
#   Rscript run_extension_chain.R <spec_name> <chain_id>
# spec_name must match one of the 18 names in extension_specs.rda
# chain_id in {1, 2, 3, 4}; seeds at 20260509 + chain_id - 1
#
# Outputs:  gibbs_<spec_name>/chain_<chain_id>.rda

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript run_extension_chain.R <spec_name> <chain_id>")
spec_name <- args[[1]]
chain_id  <- as.integer(args[[2]])

user_lib <- file.path(Sys.getenv("LOCALAPPDATA"), "R/win-library/4.6")
.libPaths(c(user_lib, .libPaths()))

run_dir <- "C:/FkCancer/runs/run_extension_2026-05-10"
val_dir <- "C:/FkCancer/runs/run_validation_2026-05-09"
source(file.path(run_dir, "sampler_composition.R"))

load(file.path(val_dir, "data/train_data.rda"))
load(file.path(run_dir, "data/extension_specs.rda"))

if (!spec_name %in% names(extension_specs))
  stop(sprintf("Unknown spec name '%s'. Available: %s",
               spec_name, paste(names(extension_specs), collapse=", ")))
spec <- extension_specs[[spec_name]]

Covariates <- Covariates_train
Survival   <- Survival_train
Censored   <- Censored_train
covariates_by_cancer_full <- lapply(Covariates, function(c) as.numeric(colnames(c)))
covariates_table <- table(unlist(covariates_by_cancer_full))
covariates_in_model_full <- as.numeric(names(covariates_table))
covariates_in_model <- covariates_in_model_full
p <- length(covariates_in_model)
covariates_by_cancer <- covariates_by_cancer_full
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

seed_val <- 20260509L + (chain_id - 1L)
set.seed(seed_val)
cat(sprintf("=== %s chain %d, seed %d, iters %d ===\n",
            spec_name, chain_id, seed_val, iters))
cat(sprintf("    cov %s, RULE %d, pair %s\n",
            spec$cov_id, spec$rule, spec$pair_key))

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
out <- HierarchicalLogNormalSpikeSlab_Composition(
  Covariates, Survival, Censored,
  starting_values, iters, covariates_in_model,
  covariates_by_cancer, priors,
  pi_generation = "shared_across_cancers",
  target_covs = spec$target_covs,
  progress = "model_comparison"
)
t1 <- Sys.time()
cat(sprintf("[%s] Finished\n", format(t1)))
cat(sprintf("Elapsed: %.2f minutes\n",
            as.numeric(difftime(t1, t0, units = "mins"))))

out_dir <- file.path(run_dir, sprintf("gibbs_%s", spec_name))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_file <- file.path(out_dir, sprintf("chain_%d.rda", chain_id))
save(out, file = out_file)
cat(sprintf("Saved %s\n", out_file))
