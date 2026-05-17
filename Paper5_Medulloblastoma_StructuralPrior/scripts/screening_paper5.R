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

# Add intercept (cov_id 0) and age-standardized (cov_id 0.5) per Lock convention.
# Some patients have OS data but missing age_at_dx; impute to cohort median
# before standardizing so X has no NA (else eigen(Sigma) inside the sampler
# blows up at iter 0).
if ("age_at_dx" %in% names(clin_s)) {
  age_raw <- as.numeric(clin_s$age_at_dx)
  age_raw[is.na(age_raw)] <- median(age_raw, na.rm = TRUE)
  age_z <- as.numeric(scale(age_raw))
  cat(sprintf("age_at_dx: imputed %d NAs with cohort median %.2f\n",
              sum(is.na(clin_s$age_at_dx)), median(age_raw, na.rm = TRUE)))
} else {
  age_z <- rep(0, nrow(clin_s))
}
stopifnot(!any(is.na(age_z)))

build_lists <- function(target_cov_idx = NULL) {
  Covariates <- list(); Survival <- list(); Censored <- list()
  for (s in subtypes) {
    keep <- which(sub_s == s)
    X_sub <- cbind(intercept = 1, age = age_z[keep], X_s[keep, , drop = FALSE])
    colnames(X_sub) <- c("0", "0.5", as.character(seq_len(ncol(X_s))))
    Covariates[[s]] <- X_sub
    # Lock sampler is log-normal AFT — takes raw survival in years and logs
    # internally (`log(Y[[k]])`). Do NOT pre-log here. Clip to ≥ 1 day to
    # avoid log(0) inside the sampler.
    Survival[[s]]   <- pmax(clin_s$os_time_years[keep], 1/365)
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

gibbs_cache <- file.path(data_dir, "screening_baseline_gibbs.rds")
if (file.exists(gibbs_cache) && file.info(gibbs_cache)$size > 1e6) {
  cat(sprintf("[%s] Loading cached baseline Gibbs from %s\n",
              format(Sys.time()), gibbs_cache))
  out <- readRDS(gibbs_cache)
} else {
  cat(sprintf("[%s] Starting baseline Gibbs (screening fit)\n", format(Sys.time())))
  out <- HierarchicalLogNormalSpikeSlab_StructuralPrior(
    Covariates, Survival, Censored,
    starting_values, iters, covariates_in_model,
    covariates_by_cancer, priors,
    pi_generation = "shared_across_cancers",
    target_covs = list(),
    progress = "normal"
  )
  saveRDS(out, gibbs_cache)
}

# --- Extract per-subtype PIPs and ANOVA against L2 -----------------------
# Lock sampler quirk: out$gamma[[i]][[1]] is a named numeric vector with
# cov_id labels, but post-burn iterations drop the names. Take names from
# iter 1 and apply by column position after rbind.
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
  # Align by name (some cov_ids may not be present per subtype if covariate
  # restriction were applied; here all 73 should be present).
  in_pip <- intersect(ref_names, colnames(pip))
  pip[i, in_pip] <- pip_means[in_pip]
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

cat(sprintf("\nPIP matrix: %d subtypes × %d covariates\n", nrow(pip), ncol(pip)))
cat(sprintf("  fully-NA columns: %d / %d\n",
            sum(apply(pip, 2, function(v) all(is.na(v)))), ncol(pip)))
cat(sprintf("  partial-NA columns: %d / %d\n",
            sum(apply(pip, 2, function(v) any(is.na(v)) && !all(is.na(v)))), ncol(pip)))

for (j in seq_len(p)) {
  y <- pip[, j]
  # Skip if entirely NA or has zero variance (degenerate covariate)
  if (all(is.na(y))) next
  s_y <- sd(y, na.rm = TRUE)
  if (is.na(s_y) || s_y < 1e-12) next
  # Pair against grp labels; drop subtypes with NA PIP
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
