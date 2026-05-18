# Paper 5 — leakage-clean screening: ANOVA of per-L2-subtype PIPs on
# TRAINING FOLD ONLY (n=494).
#
# Differs from screening_paper5.R in two ways:
#   1. Uses leakage-clean BIDIFAC+ components (data/bidifac_components_clean.rds),
#      i.e. components derived from training-fold BIDIFAC+ + projected test.
#   2. Restricts the screening Gibbs to training rows (those with
#      bidifac_fold == "train" in cavalli_clinical_aligned.rds).
#
# This makes the covariate selection itself leakage-free, as required by
# DEVIATIONS.md Entry 001 (the original screening peeked at test fold).
#
# Reads:  data/bidifac_components_clean.rds, data/cavalli_clinical_aligned.rds
# Writes: results/paper5_screening_results_clean.csv (all components: η², F, p)
#         reference/paper5_target_covariates_clean.csv (top 3 by η²)
#         data/screening_baseline_gibbs_clean.rds (cached Gibbs fit)

user_lib <- file.path(Sys.getenv("LOCALAPPDATA"), "R/win-library/4.6")
.libPaths(c(user_lib, .libPaths()))

repo_root <- "C:/Cancer Research/Paper5_Medulloblastoma_StructuralPrior"
data_dir <- file.path(repo_root, "data")
ref_dir  <- file.path(repo_root, "reference")
res_dir  <- file.path(repo_root, "results")
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

sampler_path <- "C:/Cancer Research/scripts/validation/sampler_structural_prior.R"
stopifnot(file.exists(sampler_path))
source(sampler_path)

# --- Load leakage-clean components + clinical ----------------------------
X_all <- readRDS(file.path(data_dir, "bidifac_components_clean.rds"))
clin  <- readRDS(file.path(data_dir, "cavalli_clinical_aligned.rds"))
stopifnot(nrow(X_all) == nrow(clin))
n_comp <- ncol(X_all)
cat(sprintf("Clean components matrix: %d samples × %d comps\n",
            nrow(X_all), n_comp))

# Restrict to TRAINING fold (rows with bidifac_fold == "train")
train_mask <- clin$bidifac_fold == "train"
stopifnot(sum(train_mask) > 400)
clin_train <- clin[train_mask, ]
X_train_full <- X_all[train_mask, , drop = FALSE]
cat(sprintf("Training rows: %d\n", nrow(clin_train)))
stopifnot(!anyNA(X_train_full))

# Also require non-NA OS for the screening Gibbs (should be all of them)
with_surv <- !is.na(clin_train$os_time_years) & !is.na(clin_train$os_event)
clin_s <- clin_train[with_surv, ]
X_s <- X_train_full[with_surv, , drop = FALSE]
cat(sprintf("Train rows with OS: %d\n", nrow(clin_s)))

subtype <- as.character(clin_s$subtype)
subtypes <- sort(unique(subtype))
n_subtype <- length(subtypes)
cat(sprintf("L2 subtypes (train): %d\n", n_subtype))

if ("age_at_dx" %in% names(clin_s)) {
  age_raw <- as.numeric(clin_s$age_at_dx)
  age_raw[is.na(age_raw)] <- median(age_raw, na.rm = TRUE)
  age_z <- as.numeric(scale(age_raw))
} else {
  age_z <- rep(0, nrow(clin_s))
}
stopifnot(!any(is.na(age_z)))

build_lists <- function() {
  Covariates <- list(); Survival <- list(); Censored <- list()
  for (s in subtypes) {
    keep <- which(subtype == s)
    X_sub <- cbind(intercept = 1, age = age_z[keep], X_s[keep, , drop = FALSE])
    colnames(X_sub) <- c("0", "0.5", as.character(seq_len(ncol(X_s))))
    Covariates[[s]] <- X_sub
    Survival[[s]]   <- pmax(clin_s$os_time_years[keep], 1/365)
    Censored[[s]]   <- 1L - as.integer(clin_s$os_event[keep])
  }
  names(Covariates) <- subtypes
  list(Covariates = Covariates, Survival = Survival, Censored = Censored)
}

dat <- build_lists()
Covariates <- dat$Covariates
Survival   <- dat$Survival
Censored   <- dat$Censored

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
set.seed(20260517L)
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

gibbs_cache <- file.path(data_dir, "screening_baseline_gibbs_clean.rds")
if (file.exists(gibbs_cache) && file.info(gibbs_cache)$size > 1e6) {
  cat(sprintf("[%s] Loading cached clean baseline Gibbs from %s\n",
              format(Sys.time()), gibbs_cache))
  out <- readRDS(gibbs_cache)
} else {
  cat(sprintf("[%s] Starting clean baseline Gibbs (training-fold only)\n",
              format(Sys.time())))
  out <- HierarchicalLogNormalSpikeSlab_StructuralPrior(
    Covariates, Survival, Censored,
    starting_values, iters, covariates_in_model,
    covariates_by_cancer, priors,
    pi_generation = "shared_across_cancers",
    target_covs = list(),
    progress = "normal"
  )
  saveRDS(out, gibbs_cache)
  cat(sprintf("[%s] Saved screening Gibbs to %s\n",
              format(Sys.time()), gibbs_cache))
}

# --- Extract per-subtype PIPs and ANOVA against L1 -----------------------
burn <- 50001L
draws <- (burn + 1L):iters
pip <- matrix(NA_real_, nrow = n_subtype, ncol = p,
              dimnames = list(subtypes, as.character(covariates_in_model)))
for (i in seq_len(n_subtype)) {
  gg <- out$gamma[[i]]
  ref_names <- names(gg[[1]])
  if (is.null(ref_names)) {
    cat(sprintf("[WARN] subtype %s: gamma iter 1 unnamed; skipping\n", subtypes[i]))
    next
  }
  M <- do.call(rbind, gg[draws])
  colnames(M) <- ref_names
  pip_means <- colMeans(M)
  in_pip <- intersect(ref_names, colnames(pip))
  pip[i, in_pip] <- pip_means[in_pip]
}

screening <- data.frame(
  cov_id = covariates_in_model,
  eta_sq = NA_real_,
  F_stat = NA_real_,
  p_val  = NA_real_
)

sub_to_grp <- tapply(as.character(clin_s$subgroup), subtype, function(x) x[1])
grp_per_subtype <- sub_to_grp[subtypes]
stopifnot(!any(is.na(grp_per_subtype)))

cat(sprintf("\nPIP matrix: %d subtypes × %d covariates\n", nrow(pip), ncol(pip)))
cat(sprintf("  fully-NA columns: %d / %d\n",
            sum(apply(pip, 2, function(v) all(is.na(v)))), ncol(pip)))

for (j in seq_len(p)) {
  y <- pip[, j]
  if (all(is.na(y))) next
  s_y <- sd(y, na.rm = TRUE)
  if (is.na(s_y) || s_y < 1e-12) next
  keep_idx <- !is.na(y)
  y_k <- y[keep_idx]
  g_k <- grp_per_subtype[keep_idx]
  if (length(unique(g_k)) < 2L) next
  ss_within  <- sum(tapply(y_k, g_k, function(z) sum((z - mean(z))^2)),
                    na.rm = TRUE)
  ss_between <- sum(tapply(y_k, g_k, function(z) length(z) * (mean(z) - mean(y_k))^2),
                    na.rm = TRUE)
  ss_total   <- ss_within + ss_between
  if (ss_total <= 0) next
  screening$eta_sq[j] <- ss_between / ss_total
  k <- length(unique(g_k))
  n <- length(y_k)
  if (k > 1 && n > k && ss_within > 0) {
    screening$F_stat[j] <- (ss_between / (k - 1)) / (ss_within / (n - k))
    screening$p_val[j]  <- pf(screening$F_stat[j], k - 1, n - k, lower.tail = FALSE)
  }
}

candidates <- screening[!(screening$cov_id %in% c(0, 0.5)) &
                          !is.na(screening$eta_sq), ]
candidates <- candidates[order(-candidates$eta_sq), ]
write.csv(candidates, file.path(res_dir, "paper5_screening_results_clean.csv"),
          row.names = FALSE)

K <- 3L
target <- head(candidates, K)
target$rank <- seq_len(nrow(target))
write.csv(target, file.path(ref_dir, "paper5_target_covariates_clean.csv"),
          row.names = FALSE)

cat(sprintf("\nTop %d target covariates (LEAKAGE-CLEAN screening):\n", K))
print(target, row.names = FALSE)
cat(sprintf("\nReference frozen at %s\n",
            file.path(ref_dir, "paper5_target_covariates_clean.csv")))
cat(sprintf("Commit this file before firing leakage-clean Gibbs.\n"))
