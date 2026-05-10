# Single-chain runner for the validation Gibbs.
#
# Usage:
#   Rscript run_chain.R <spec> <chain_id>
#     spec ∈ {baseline_train, primary, secondary}
#     chain_id ∈ {1, 2, 3, 4}
# Seeds: 20260509 + (chain_id - 1) [matches baseline pattern]
#
# Outputs: spec_chain_<id>.rda containing 'out' (full posterior list)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript run_chain.R <spec> <chain_id>")
spec <- args[[1]]
chain_id <- as.integer(args[[2]])
stopifnot(spec %in% c("baseline_train", "primary", "secondary"))

user_lib <- file.path(Sys.getenv("LOCALAPPDATA"), "R/win-library/4.6")
.libPaths(c(user_lib, .libPaths()))

run_dir <- "C:/FkCancer/runs/run_validation_2026-05-09"
source(file.path(run_dir, "sampler_structural_prior.R"))

load(file.path(run_dir, "data/train_data.rda"))
Covariates <- Covariates_train
Survival   <- Survival_train
Censored   <- Censored_train

target_covs <- list()
if (spec == "primary") {
  load(file.path(run_dir, "data/target_covs_primary.rda"))
  target_covs <- target_covs_primary
} else if (spec == "secondary") {
  load(file.path(run_dir, "data/target_covs_secondary.rda"))
  target_covs <- target_covs_secondary
}

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

out_dir <- file.path(run_dir, sprintf("gibbs_%s", spec))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_file <- file.path(out_dir, sprintf("chain_%d.rda", chain_id))
save(out, file = out_file)
cat(sprintf("Saved %s\n", out_file))
