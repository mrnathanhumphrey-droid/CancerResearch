# Paper 5 — build held-out 80/20 split, stratified at L2 (12-subtype).
# Build Lock-style Covariates/Survival/Censored lists for train and test.
# Build target_covs_primary (L2 hierarchy) + target_covs_sensitivity (L1).
#
# Reads:  data/bidifac_components.rds, data/cavalli_clinical_aligned.rds,
#         reference/paper5_target_covariates.csv
# Writes: reference/paper5_split_indices.rda
#         data/train_data.rda, data/test_data.rda
#         data/target_covs_primary.rda, data/target_covs_sensitivity.rda

user_lib <- file.path(Sys.getenv("LOCALAPPDATA"), "R/win-library/4.6")
.libPaths(c(user_lib, .libPaths()))

repo_root <- "C:/Cancer Research/Paper5_Medulloblastoma_StructuralPrior"
data_dir <- file.path(repo_root, "data")
ref_dir  <- file.path(repo_root, "reference")

X_all <- readRDS(file.path(data_dir, "bidifac_components.rds"))
clin  <- readRDS(file.path(data_dir, "cavalli_clinical_aligned.rds"))
target_table <- read.csv(file.path(ref_dir, "paper5_target_covariates.csv"),
                         stringsAsFactors = FALSE)
target_cov_ids <- target_table$cov_id  # length K=3
stopifnot(length(target_cov_ids) == 3)

with_surv <- !is.na(clin$os_time_years) & !is.na(clin$os_event)
clin_s <- clin[with_surv, ]
X_s <- X_all[with_surv, , drop = FALSE]

subtype <- as.character(clin_s$subtype)
subgroup <- as.character(clin_s$subgroup)
subtypes <- sort(unique(subtype))
subgroups <- sort(unique(subgroup))
cat(sprintf("With-survival n: %d. L2 subtypes: %d. L1 subgroups: %d.\n",
            nrow(clin_s), length(subtypes), length(subgroups)))

# --- Stratified 80/20 at L2 ----------------------------------------------
set.seed(20260515L)
train_idx <- list(); test_idx <- list()
for (s in subtypes) {
  idx <- which(subtype == s)
  n_s <- length(idx)
  n_test <- max(1L, floor(0.20 * n_s))
  perm <- sample(idx, n_s)
  test_idx[[s]] <- sort(perm[seq_len(n_test)])
  train_idx[[s]] <- sort(perm[(n_test + 1L):n_s])
}
cat("\nL2-stratified split summary (train_n / test_n):\n")
for (s in subtypes) {
  cat(sprintf("  %-12s  train=%4d  test=%3d\n",
              s, length(train_idx[[s]]), length(test_idx[[s]])))
}

save(train_idx, test_idx,
     file = file.path(ref_dir, "paper5_split_indices.rda"))

# --- Build Lock-style Covariates / Survival / Censored ------------------
age_z <- if ("age_at_dx" %in% names(clin_s)) {
  as.numeric(scale(clin_s$age_at_dx))
} else { rep(0, nrow(clin_s)) }

build_cov_for_subtype <- function(idx_set) {
  out <- list(); surv <- list(); cens <- list()
  for (s in subtypes) {
    keep <- intersect(idx_set[[s]], which(subtype == s))
    # Subset by relative index since idx_set returned absolute indices into clin_s
    X_sub <- cbind(intercept = 1, age = age_z[keep], X_s[keep, , drop = FALSE])
    colnames(X_sub) <- c("0", "0.5", as.character(seq_len(ncol(X_s))))
    out[[s]]  <- X_sub
    surv[[s]] <- log(pmax(clin_s$os_time_years[keep], 1/365))
    cens[[s]] <- 1L - as.integer(clin_s$os_event[keep])
  }
  names(out) <- subtypes; names(surv) <- subtypes; names(cens) <- subtypes
  list(Covariates = out, Survival = surv, Censored = cens)
}

train <- build_cov_for_subtype(train_idx)
test  <- build_cov_for_subtype(test_idx)
Covariates_train <- train$Covariates
Survival_train   <- train$Survival
Censored_train   <- train$Censored
Covariates_test  <- test$Covariates
Survival_test    <- test$Survival
Censored_test    <- test$Censored

save(Covariates_train, Survival_train, Censored_train,
     file = file.path(data_dir, "train_data.rda"))
save(Covariates_test, Survival_test, Censored_test,
     file = file.path(data_dir, "test_data.rda"))

# --- L1 (4-subgroup) mapping per subtype --------------------------------
sub_to_grp <- tapply(subgroup, subtype, function(x) x[1])
l1_per_subtype <- sub_to_grp[subtypes]
l2_per_subtype <- subtypes  # each subtype is its own L2 bin
names(l1_per_subtype) <- subtypes

# --- Build target_covs lists --------------------------------------------
build_target <- function(cov_id, hier_label, group_vec) {
  list(
    cov_id = cov_id,
    hierarchy_label = hier_label,
    group_per_cancer = unname(group_vec[subtypes])
  )
}

# Primary: L2 hierarchy means each subtype is its own structural-prior bin
# WITH the partition being identity at L2 — that's degenerate (singleton bins).
# CORRECTED PER PRE-REG §3: L2 means each subtype shares its β WITH OTHER SUBTYPES IN ITS L1 SUBGROUP.
# The "12-subtype hierarchy at granularity 12" means the per-subtype β is pooled
# at the parent L1 subgroup level (so WNT-α and WNT-β share, etc.).
# In the Lock sampler, the structural-prior group is given by the hierarchy mapping;
# we're using L1 subgroup as the group label for L2-level pooling.

target_covs_primary <- lapply(target_cov_ids, function(cid) {
  build_target(cid, "L2_12subtype_pooled_by_L1_subgroup", l1_per_subtype)
})

target_covs_sensitivity <- lapply(target_cov_ids, function(cid) {
  build_target(cid, "L1_4subgroup", l1_per_subtype)
})
# Note: when "cancer type" = L2 subtype, the L1 mapping above also functions as
# the L1 sensitivity hierarchy. The difference between primary and sensitivity
# in Paper 5 is the level at which "cancer type" lives (L2 vs L1), not the
# group structure. See PRE_REGISTRATION.md §3. The sensitivity spec aggregates
# data to the 4 L1 subgroups (separate train/test build below, NOT this file).

save(target_covs_primary, subtypes,
     file = file.path(data_dir, "target_covs_primary.rda"))
save(target_covs_sensitivity, subtypes,
     file = file.path(data_dir, "target_covs_sensitivity.rda"))

cat("\nPRIMARY target_covs (L2 cancer-types, L1 hierarchy bins for pooling):\n")
for (t in target_covs_primary) {
  cat(sprintf("  cov %s @ %s: bins = %s\n",
              t$cov_id, t$hierarchy_label,
              paste(unique(t$group_per_cancer), collapse = " / ")))
}
cat(sprintf("\nSaved split + train/test + target_covs.\n"))
