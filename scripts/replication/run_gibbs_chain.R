# Parameterized Gibbs runner — accepts a chain ID via Rscript args.
# Usage: Rscript run_gibbs_chain.R <chain_id>
# Seeds at 20260509 + chain_id - 1 (chain 1 → 20260509, matching the original).

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript run_gibbs_chain.R <chain_id>")
chain_id <- as.integer(args[[1]])
seed_val <- 20260509L + (chain_id - 1L)

cat(sprintf("=== Chain %d, seed %d ===\n", chain_id, seed_val))

user_lib <- file.path(Sys.getenv("LOCALAPPDATA"), "R/win-library/4.6")
.libPaths(c(user_lib, .libPaths()))
repo <- "C:/FkCancer/repos/HierarchicalSS_PanCanPanOmics"
setwd(repo)

load("XYC_V2_WithAge_StandardizedPredictors.rda")
source("ExtendedHierarchicalModelLogNormalSpikeAndSlabInterceptOutsideSS.R")
source("HelperFunctions.R")

covariates_by_cancer_full <- lapply(Covariates, function(cancer) as.numeric(colnames(cancer)))
covariates_table <- table(unlist(covariates_by_cancer_full))
covariates_in_model_full <- as.numeric(names(covariates_table))
covariates_in_model <- covariates_in_model_full
p <- length(covariates_in_model)
covariates_by_cancer <- covariates_by_cancer_full
iters <- 1e5
n_cancer <- length(Covariates)
cancer_types <- names(Covariates)

priors <- list(betatilde_priorvar_intercept = 10^2, betatilde_priorvar_coefficient = 1,
               lambda2_priorshape_intercept = 1, lambda2_priorrate_intercept = 1,
               lambda2_priorshape_coefficient = 5, lambda2_priorrate_coefficient = 1,
               sigma2_priorshape = .01, sigma2_priorrate = .01, spike_priorvar = 1/10000)

set.seed(seed_val)
beta_tilde_start <- rnorm(p, 0, 1)
sigma2_start <- 1/rgamma(1, 1, 1)
lambda2_start <- 1/rgamma(p, 1, 1)
pi_start <- rep(0.5, p)
gamma_start <- lapply(1:n_cancer, function(type) {
  indicators <- rbinom(length(covariates_by_cancer[[type]]), size = 1, prob = 1)
  indicators[1] <- 1
  names(indicators) <- covariates_by_cancer[[type]]
  list(indicators)
})
starting_values <- list(sigma2_start = sigma2_start, beta_tilde_start = beta_tilde_start,
                        lambda2_start = lambda2_start, gamma_start = gamma_start, pi_start = pi_start)

cat(sprintf("Chain %d starting Gibbs at %s\n", chain_id, format(Sys.time())))
t0 <- Sys.time()
out <- HierarchicalLogNormalSpikeSlab(
  Covariates, Survival, Censored,
  starting_values, iters, covariates_in_model,
  covariates_by_cancer, priors,
  pi_generation = "shared_across_cancers"
)
t1 <- Sys.time()
cat(sprintf("Chain %d finished at %s\n", chain_id, format(t1)))
cat("Elapsed:\n"); print(t1 - t0)

# Save with chain-specific filenames
out_run  <- sprintf("C:/FkCancer/runs/run_paper_a_gibbs_2026-05-09/GibbsSamplingResults_100kiters_Chain%d_WithAge_StandardizedPredictors.rda", chain_id)
out_repo <- sprintf("%s/GibbsSamplingResults_100kiters_Chain%d_WithAge_StandardizedPredictors.rda", repo, chain_id)
# Save under the canonical name the upstream code expects
assign(sprintf("Posteriors_100k_WithAge_StandardizedPredictors_Chain%d", chain_id), out)
save(list = sprintf("Posteriors_100k_WithAge_StandardizedPredictors_Chain%d", chain_id),
     file = out_run)
save(list = sprintf("Posteriors_100k_WithAge_StandardizedPredictors_Chain%d", chain_id),
     file = out_repo)
cat(sprintf("Chain %d saved to %s\n", chain_id, out_run))
