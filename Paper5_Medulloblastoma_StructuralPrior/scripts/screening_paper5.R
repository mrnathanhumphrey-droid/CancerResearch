# Paper 5 — screening pass: ANOVA of per-L2-subtype PIPs against L2 hierarchy.
#
# For each BIDIFAC+ component, fits a baseline Gibbs (Lock 2022 sampler with
# empty target_covs) and extracts per-L2-subtype posterior inclusion probabilities
# (PIPs). ANOVA of those PIPs against L2 subtype label gives an η² ranking;
# the top 3 components by effect size are pre-registered as the target covariates
# for Paper 5 (PRE_REGISTRATION.md §4).
#
# Reads:  data/bidifac_components.rds, data/cavalli_clinical_aligned.rds
# Writes: results/paper5_screening_results.csv (all components: η², F, p)
#         reference/paper5_target_covariates.csv (top 3 by η²)
#
# This is a separate Gibbs fire from the validation; runs on the full data
# (not just training fold) since screening must precede the train/test split.

user_lib <- file.path(Sys.getenv("LOCALAPPDATA"), "R/win-library/4.6")
.libPaths(c(user_lib, .libPaths()))

repo_root <- "C:/Cancer Research/Paper5_Medulloblastoma_StructuralPrior"
data_dir <- file.path(repo_root, "data")
ref_dir  <- file.path(repo_root, "reference")
res_dir  <- file.path(repo_root, "results")
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

# Reuse Lock sampler from Paper 1's validation folder
sampler_path <- "C:/Cancer Research/scripts/validation/sampler_structural_prior.R"
stopifnot(file.exists(sampler_path))
source(sampler_path)

# --- Load aligned BIDIFAC+ components + clinical ------------------------
X_all <- readRDS(file.path(data_dir, "bidifac_components.rds"))  # samples × n_comp
clin  <- readRDS(file.path(data_dir, "cavalli_clinical_aligned.rds"))
stopifnot(nrow(X_all) == nrow(clin))
n_comp <- ncol(X_all)
cat(sprintf("Screening on %d samples × %d BIDIFAC+ components\n",
            nrow(X_all), n_comp))

# --- Build Lock-style Covariates / Survival / Censored per L2 subtype ----
# "Cancer types" in Paper 1's terminology -> L2 subtypes here.
subtype <- as.character(clin$subtype)
with_surv <- !is.na(clin$os_time_years) & !is.na(clin$os_event)
clin_s <- clin[with_surv, ]
X_s <- X_all[with_surv, , drop = FALSE]
sub_s <- subtype[with_surv]

subtypes <- sort(unique(sub_s))
n_subtype <- length(subtypes)
cat(sprintf("L2 subtypes: %d\n", n_subtype))

# Add intercept (cov_id 0) and age-standardized (cov_id 0.5) per Lock convention
age_z <- if ("age_at_dx" %in% names(clin_s)) {
  as.numeric(scale(clin_s$age_at_dx))
} else {
  rep(0, nrow(clin_s))  # if absent, age effect is zeroed
}

build_lists <- function(target_cov_idx = NULL) {
  Covariates <- list(); Survival <- list(); Censored <- list()
  for (s in subtypes) {
    keep <- which(sub_s == s)
    X_sub <- cbind(intercept = 1, age = age_z[keep], X_s[keep, , drop = FALSE])
    colnames(X_sub) <- c("0", "0.5", as.character(seq_len(ncol(X_s))))
    Covariates[[s]] <- X_sub
    Survival[[s]]   <- log(pmax(clin_s$os_time_years[keep], 1/365))  # log-AFT
    Censored[[s]]   <- 1L - as.integer(clin_s$os_event[keep])  # Lock: 1=censored
  }
  names(Covariates) <- subtypes
  list(Covariates = Covariates, Survival = Survival, Censored = Censored)
}

dat <- build_lists()
Covariates <- dat$Covariates
Survival   <- dat$Survival
Censored   <- dat$Censored

# --- Run baseline Gibbs (target_covs empty) ------------------------------
covariates_by_cancer <- lapply(Covariates, function(c) as.numeric(colnames(c)))
covariates_in_model <- as.numeric(names(table(unlist(covariates_by_cancer))))
p <- length(covariates_in_model)
n_cancer <- length(Covariates)

priors <- list(betatilde_priorvar_intercept = 10^2,
               betatilde_priorvar_coefficient = 1,
               lambda2_priorshape_intercept = 1,
               lambda2_priorrate_intercept = 1,
               lambda2_priorshape_coefficient = 5,
               lambda2_priorrate_coefficient = 1,
               sigma2_priorshape = .01,
               sigma2_priorrate = .01,
               spike_priorvar = 1/10000)

iters <- 1e5
set.seed(20260515L)
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

cat(sprintf("[%s] Starting baseline Gibbs (screening fit)\n", format(Sys.time())))
out <- HierarchicalLogNormalSpikeSlab_StructuralPrior(
  Covariates, Survival, Censored,
  starting_values, iters, covariates_in_model,
  covariates_by_cancer, priors,
  pi_generation = "shared_across_cancers",
  target_covs = list(),
  progress = "normal"
)
saveRDS(out, file.path(data_dir, "screening_baseline_gibbs.rds"))

# --- Extract per-subtype PIPs and ANOVA against L2 -----------------------
burn <- 50001L
draws <- (burn + 1L):iters
pip <- matrix(NA_real_, nrow = n_subtype, ncol = p,
              dimnames = list(subtypes, as.character(covariates_in_model)))
for (i in seq_len(n_subtype)) {
  gg <- out$gamma[[i]]
  M <- do.call(rbind, gg[draws])
  pip[i, colnames(M)] <- colMeans(M)
}

screening <- data.frame(
  cov_id = covariates_in_model,
  eta_sq = NA_real_,
  F_stat = NA_real_,
  p_val  = NA_real_
)
# ANOVA: per-subtype PIP grouped by subtype label (1 row per subtype, so single-cell
# ANOVA doesn't make sense at the subtype level — instead, group at the L1 subgroup
# level using subtype→subgroup mapping). Pull the L1 mapping from clinical data.
sub_to_grp <- tapply(as.character(clin_s$subgroup), sub_s, function(x) x[1])
grp_per_subtype <- sub_to_grp[subtypes]
stopifnot(!any(is.na(grp_per_subtype)))

for (j in seq_len(p)) {
  y <- pip[, j]
  if (sd(y, na.rm = TRUE) < 1e-12) next
  ss_within  <- sum(tapply(y, grp_per_subtype, function(z) sum((z - mean(z))^2)),
                    na.rm = TRUE)
  ss_between <- sum(tapply(y, grp_per_subtype, function(z) length(z) * (mean(z) - mean(y))^2),
                    na.rm = TRUE)
  ss_total   <- ss_within + ss_between
  screening$eta_sq[j] <- ss_between / ss_total
  k <- length(unique(grp_per_subtype))
  n <- length(y)
  if (k > 1 && n > k && ss_within > 0) {
    screening$F_stat[j] <- (ss_between / (k - 1)) / (ss_within / (n - k))
    screening$p_val[j]  <- pf(screening$F_stat[j], k - 1, n - k, lower.tail = FALSE)
  }
}

# Exclude intercept (cov_id 0) and age (cov_id 0.5) from target-covariate ranking
candidates <- screening[!(screening$cov_id %in% c(0, 0.5)) &
                          !is.na(screening$eta_sq), ]
candidates <- candidates[order(-candidates$eta_sq), ]
write.csv(candidates, file.path(res_dir, "paper5_screening_results.csv"),
          row.names = FALSE)

# --- Pick top-K by η² (K = 3 per pre-reg §4) -----------------------------
K <- 3L
target <- head(candidates, K)
target$rank <- seq_len(nrow(target))
write.csv(target, file.path(ref_dir, "paper5_target_covariates.csv"),
          row.names = FALSE)

cat(sprintf("\nTop %d target covariates (locked from this output):\n", K))
print(target, row.names = FALSE)
cat(sprintf("\nReference frozen at %s\n",
            file.path(ref_dir, "paper5_target_covariates.csv")))
cat(sprintf("Commit this file to git before run_gibbs_paper5.R fires.\n"))
