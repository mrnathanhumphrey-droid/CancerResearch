# Smoke test: 1000 iters on training data for one Rule 1, one Rule 2, one
# Rule 3 spec. Verifies the composition sampler produces sensible values.

user_lib <- file.path(Sys.getenv("LOCALAPPDATA"), "R/win-library/4.6")
.libPaths(c(user_lib, .libPaths()))

run_dir <- "C:/FkCancer/runs/run_extension_2026-05-10"
val_dir <- "C:/FkCancer/runs/run_validation_2026-05-09"
source(file.path(run_dir, "sampler_composition.R"))
load(file.path(val_dir, "data/train_data.rda"))
load(file.path(run_dir, "data/extension_specs.rda"))

# Pick three representative specs (one per rule) on cov 8.2 pair A
test_specs <- c(
  "spec_01_cov8_2_rule1_pairA",  # cov 8.2, RULE 1, pair A
  "spec_03_cov8_2_rule2_pairA",  # cov 8.2, RULE 2, pair A
  "spec_05_cov8_2_rule3_pairA"   # cov 8.2, RULE 3, pair A
)

Covariates <- Covariates_train
Survival   <- Survival_train
Censored   <- Censored_train
covariates_by_cancer_full <- lapply(Covariates, function(c) as.numeric(colnames(c)))
covariates_in_model <- as.numeric(names(table(unlist(covariates_by_cancer_full))))
p <- length(covariates_in_model)
covariates_by_cancer <- covariates_by_cancer_full
n_cancer <- length(Covariates)

iters <- 500L
priors <- list(betatilde_priorvar_intercept = 10^2,
               betatilde_priorvar_coefficient = 1,
               lambda2_priorshape_intercept = 1,
               lambda2_priorrate_intercept = 1,
               lambda2_priorshape_coefficient = 5,
               lambda2_priorrate_coefficient = 1,
               sigma2_priorshape = .01,
               sigma2_priorrate = .01,
               spike_priorvar = 1/10000)

set.seed(20260509L)
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

for (sp_name in test_specs) {
  spec <- extension_specs[[sp_name]]
  cat(sprintf("\n========== %s (RULE %d) ==========\n", sp_name, spec$rule))
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
  cat(sprintf("Smoke test elapsed: %.2f sec (%.4f s/iter)\n",
              as.numeric(difftime(t1, t0, units = "secs")),
              as.numeric(difftime(t1, t0, units = "secs"))/iters))
  cat(sprintf("sigma2 last 100 iters: mean=%.3f range=[%.3f, %.3f]\n",
              mean(out$sigma2[(iters-100):iters], na.rm=TRUE),
              min(out$sigma2[(iters-100):iters], na.rm=TRUE),
              max(out$sigma2[(iters-100):iters], na.rm=TRUE)))
  cat(sprintf("any NA in sigma2: %s\n", any(is.na(out$sigma2))))

  # Check the structural override: cancers in the same target group should share beta
  tm_full <- spec$target_covs[[1]]
  cov_id <- tm_full$cov_id
  rule <- tm_full$rule
  gpc_h1 <- tm_full$pair$gpc_h1
  gpc_h2 <- tm_full$pair$gpc_h2

  last_iter <- iters
  betas_per_cancer <- sapply(1:n_cancer, function(c) {
    k_idx <- which(covariates_by_cancer[[c]] == cov_id)
    if (length(k_idx) == 0) return(NA_real_)
    out$betas[[c]][[last_iter]][k_idx]
  })
  names(betas_per_cancer) <- names(Covariates_train)
  cat(sprintf("Per-cancer beta for cov %s at iter %d (NA = cov absent):\n",
              cov_id, last_iter))
  print(round(betas_per_cancer, 4))

  # Check group-shared structure
  if (rule == 1L) {
    cat("\nRule 1 check: each cancer's beta should equal (h1_g + h2_g)/2.\n")
    # For each cancer with both labels non-NA, compare beta to expected
  } else if (rule == 2L) {
    cat("\nRule 2 check: cancers in same joint (h1, h2) group should share beta.\n")
    joint <- ifelse(is.na(gpc_h1) | is.na(gpc_h2), NA,
                    paste(gpc_h1, gpc_h2, sep = "__"))
    for (g in unique(joint[!is.na(joint)])) {
      cancers <- which(joint == g)
      bs <- betas_per_cancer[cancers]
      bs <- bs[!is.na(bs)]
      if (length(bs) >= 2) {
        ok <- abs(diff(range(bs))) < 1e-9
        cat(sprintf("  joint '%s': n=%d, share? %s, betas: %s\n",
                    g, length(bs), if (ok) "OK" else "FAIL",
                    paste(round(bs, 4), collapse = ",")))
      }
    }
  } else if (rule == 3L) {
    cat("\nRule 3 check: joint groups |g|>=2 share; singletons fall back to h1 mean.\n")
  }
}

cat("\nSmoke test complete.\n")
