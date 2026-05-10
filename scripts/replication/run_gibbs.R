# Runner for Lock-lab pan-cancer hierarchical spike-and-slab Gibbs sampler.
# Fires the 100k-iteration Gibbs fit on the pre-computed XYC inputs.
#
# Output: GibbsSamplingResults_100kiters_Chain1_WithAge_StandardizedPredictors.rda
# Expected wall time: hours (likely 2-8 h on this hardware, depends on per-iter cost).
#
# Run via:
#   Rscript C:/FkCancer/runs/run_paper_a_gibbs_2026-05-09/run_gibbs.R 2>&1 | tee run.log

# --- Library setup ---
user_lib <- file.path(Sys.getenv("LOCALAPPDATA"), "R/win-library/4.6")
.libPaths(c(user_lib, .libPaths()))

# --- Working dir = repo (so relative paths in source/load work) ---
repo <- "C:/FkCancer/repos/HierarchicalSS_PanCanPanOmics"
setwd(repo)
cat("Working dir:", getwd(), "\n")

# --- Load the pre-computed model inputs FIRST ---
# (HelperFunctions.R references Covariates at parse time; must exist before sourcing.)
load("XYC_V2_WithAge_StandardizedPredictors.rda")

# --- Source the sampler code + helpers ---
# (HelperFunctions.R is the patched copy without the hardcoded load.)
source("ExtendedHierarchicalModelLogNormalSpikeAndSlabInterceptOutsideSS.R")
source("HelperFunctions.R")
cat("Loaded model inputs. Cancer types:", length(Covariates), "\n")
cat("Covariates list names:", paste(head(names(Covariates), 5), collapse=", "), "...\n")

# --- Setup (verbatim from upstream Gibbs script) ---
covariates_by_cancer_full <- lapply(Covariates, function(cancer) as.numeric(colnames(cancer)))
covariates_table <- table(unlist(covariates_by_cancer_full))
covariates_in_model_full <- as.numeric(names(covariates_table))
covariates_in_model <- covariates_in_model_full
p <- length(covariates_in_model)
covariates_by_cancer <- covariates_by_cancer_full
iters <- 1e5
n_cancer <- length(Covariates)
cancer_types <- names(Covariates)

cat("p (total PCs across cancer types) =", p, "\n")
cat("iters =", iters, "\n")
cat("n_cancer =", n_cancer, "\n")

# --- Priors (verbatim) ---
priors <- list(betatilde_priorvar_intercept = 10^2,
               betatilde_priorvar_coefficient = 1,
               lambda2_priorshape_intercept = 1,
               lambda2_priorrate_intercept = 1,
               lambda2_priorshape_coefficient = 5,
               lambda2_priorrate_coefficient = 1,
               sigma2_priorshape = .01,
               sigma2_priorrate = .01,
               spike_priorvar = 1/10000)

# --- Starting values (verbatim) ---
set.seed(20260509)  # reproducibility
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

starting_values <- list(sigma2_start = sigma2_start,
                        beta_tilde_start = beta_tilde_start,
                        lambda2_start = lambda2_start,
                        gamma_start = gamma_start,
                        pi_start = pi_start)

# --- Fire the sampler ---
cat("\n=== Starting Gibbs sampler at", format(Sys.time()), "===\n")
start <- Sys.time()
Posteriors_100k_WithAge_StandardizedPredictors <- HierarchicalLogNormalSpikeSlab(
  Covariates, Survival, Censored,
  starting_values, iters, covariates_in_model,
  covariates_by_cancer, priors,
  pi_generation = "shared_across_cancers"
)
end <- Sys.time()
cat("\n=== Sampler finished at", format(end), "===\n")
cat("Elapsed:\n"); print(end - start)

# --- Save output to the run dir AND to the repo (for upstream-script compat) ---
out_run <- "C:/FkCancer/runs/run_paper_a_gibbs_2026-05-09/GibbsSamplingResults_100kiters_Chain1_WithAge_StandardizedPredictors.rda"
out_repo <- file.path(repo, "GibbsSamplingResults_100kiters_Chain1_WithAge_StandardizedPredictors.rda")
save(Posteriors_100k_WithAge_StandardizedPredictors, file = out_run)
save(Posteriors_100k_WithAge_StandardizedPredictors, file = out_repo)
cat("\nSaved:\n  ", out_run, "\n  ", out_repo, "\n")

cat("\nSession info:\n"); print(sessionInfo())
