# Paper 5 — leakage-clean train/test build (v2).
#
# Differences from build_paper5_split.R:
#   1. Uses leakage-clean BIDIFAC+ components (bidifac_components_clean.rds)
#   2. Uses LEAKAGE-CLEAN target covariates (paper5_target_covariates_clean.csv)
#   3. FIXES FLAG 1 SPEC COLLAPSE: builds a real L1 sensitivity spec at
#      cancer-level = 4 L1 subgroups with identity hierarchy (vs primary's
#      12 L2 subtypes pooled by L1). The original v1 had both specs sharing
#      the same Covariates/Survival/Censored + group_per_cancer, producing
#      bit-identical LPPD.
#
# Reads:  data/bidifac_components_clean.rds, data/cavalli_clinical_aligned.rds,
#         reference/paper5_target_covariates_clean.csv,
#         reference/paper5_split_indices.rda    (LOCKED split, unchanged)
# Writes:
#   data/train_data_clean.rda           — L2-level (12 subtypes) train lists
#   data/test_data_clean.rda            — L2-level (12 subtypes) test lists
#   data/train_data_L1_clean.rda        — L1-level (4 subgroups) train lists
#   data/test_data_L1_clean.rda         — L1-level (4 subgroups) test lists
#   data/target_covs_primary_clean.rda  — L1-pooling structural prior over 12 L2
#   data/target_covs_sensitivity_clean.rda — identity over 4 L1 (no pooling)

user_lib <- file.path(Sys.getenv("LOCALAPPDATA"), "R/win-library/4.6")
.libPaths(c(user_lib, .libPaths()))

repo_root <- "C:/Cancer Research/Paper5_Medulloblastoma_StructuralPrior"
data_dir <- file.path(repo_root, "data")
ref_dir  <- file.path(repo_root, "reference")

X_all <- readRDS(file.path(data_dir, "bidifac_components_clean.rds"))
clin  <- readRDS(file.path(data_dir, "cavalli_clinical_aligned.rds"))
stopifnot(nrow(X_all) == nrow(clin))

target_table <- read.csv(file.path(ref_dir, "paper5_target_covariates_clean.csv"),
                         stringsAsFactors = FALSE)
target_cov_ids <- target_table$cov_id
stopifnot(length(target_cov_ids) == 3L)
cat(sprintf("Leakage-clean target covariates: %s\n",
            paste(target_cov_ids, collapse = ", ")))

# --- Restrict to with-OS rows, preserving locked split ------------------
with_surv <- !is.na(clin$os_time_years) & !is.na(clin$os_event)
clin_s <- clin[with_surv, ]
X_s <- X_all[with_surv, , drop = FALSE]

# Locked split — re-load and verify same n
load(file.path(ref_dir, "paper5_split_indices.rda"))
n_train <- sum(sapply(train_idx, length))
n_test  <- sum(sapply(test_idx,  length))
stopifnot(n_train + n_test == nrow(clin_s))

subtype  <- as.character(clin_s$subtype)
subgroup <- as.character(clin_s$subgroup)
subtypes  <- sort(unique(subtype))   # 12 L2
subgroups <- sort(unique(subgroup))  # 4 L1
cat(sprintf("With-survival n: %d. L2 subtypes: %d. L1 subgroups: %d.\n",
            nrow(clin_s), length(subtypes), length(subgroups)))

# Verify NA pattern — train/test rows in clean components should be non-NA
stopifnot(!anyNA(X_s[unlist(train_idx), ]))
stopifnot(!anyNA(X_s[unlist(test_idx),  ]))

# Age standardization: cohort-median impute, then z-score (TRAINING ROWS ONLY
# for scale params, applied to all)
age_raw <- as.numeric(clin_s$age_at_dx)
age_raw[is.na(age_raw)] <- median(age_raw, na.rm = TRUE)
age_train <- age_raw[unlist(train_idx)]
age_mean <- mean(age_train); age_sd <- sd(age_train)
age_z <- (age_raw - age_mean) / age_sd
stopifnot(!any(is.na(age_z)))
cat(sprintf("Age standardization (train-only): mean=%.2f sd=%.2f\n",
            age_mean, age_sd))

# --- L2 (12-subtype) builder --------------------------------------------
build_L2 <- function(idx_set) {
  Covariates <- list(); Survival <- list(); Censored <- list()
  for (s in subtypes) {
    keep <- intersect(idx_set[[s]], which(subtype == s))
    X_sub <- cbind(intercept = 1, age = age_z[keep], X_s[keep, , drop = FALSE])
    colnames(X_sub) <- c("0", "0.5", as.character(seq_len(ncol(X_s))))
    Covariates[[s]]  <- X_sub
    Survival[[s]]    <- pmax(clin_s$os_time_years[keep], 1/365)
    Censored[[s]]    <- 1L - as.integer(clin_s$os_event[keep])
  }
  names(Covariates) <- subtypes; names(Survival) <- subtypes
  names(Censored) <- subtypes
  list(Covariates = Covariates, Survival = Survival, Censored = Censored)
}

# --- L1 (4-subgroup) builder --------------------------------------------
# For Flag 1 fix: aggregate all L2 subtypes within each L1 subgroup into a
# single cancer-type unit. Sensitivity spec runs at this granularity.
build_L1 <- function(idx_set) {
  Covariates <- list(); Survival <- list(); Censored <- list()
  for (g in subgroups) {
    # All L2 subtypes whose subgroup == g, indices from idx_set
    sub_in_g <- subtypes[
      vapply(subtypes, function(s) {
        unique(subgroup[subtype == s])[1] == g
      }, logical(1))
    ]
    keep <- unique(unlist(idx_set[sub_in_g]))
    keep <- intersect(keep, which(subgroup == g))
    X_sub <- cbind(intercept = 1, age = age_z[keep], X_s[keep, , drop = FALSE])
    colnames(X_sub) <- c("0", "0.5", as.character(seq_len(ncol(X_s))))
    Covariates[[g]] <- X_sub
    Survival[[g]]   <- pmax(clin_s$os_time_years[keep], 1/365)
    Censored[[g]]   <- 1L - as.integer(clin_s$os_event[keep])
  }
  names(Covariates) <- subgroups; names(Survival) <- subgroups
  names(Censored) <- subgroups
  list(Covariates = Covariates, Survival = Survival, Censored = Censored)
}

# --- Build all four sets ------------------------------------------------
train_L2 <- build_L2(train_idx); test_L2 <- build_L2(test_idx)
train_L1 <- build_L1(train_idx); test_L1 <- build_L1(test_idx)

# L2 train sanity
cat("\nL2 (12-subtype) — train/test row counts per subtype:\n")
for (s in subtypes) {
  cat(sprintf("  %-14s  train=%4d  test=%3d\n", s,
              length(train_L2$Survival[[s]]),
              length(test_L2$Survival[[s]])))
}
cat("\nL1 (4-subgroup) — train/test row counts per subgroup:\n")
for (g in subgroups) {
  cat(sprintf("  %-10s  train=%4d  test=%3d\n", g,
              length(train_L1$Survival[[g]]),
              length(test_L1$Survival[[g]])))
}

# Save L2 (primary)
Covariates_train <- train_L2$Covariates
Survival_train   <- train_L2$Survival
Censored_train   <- train_L2$Censored
Covariates_test  <- test_L2$Covariates
Survival_test    <- test_L2$Survival
Censored_test    <- test_L2$Censored
save(Covariates_train, Survival_train, Censored_train,
     file = file.path(data_dir, "train_data_clean.rda"))
save(Covariates_test, Survival_test, Censored_test,
     file = file.path(data_dir, "test_data_clean.rda"))

# Save L1 (sensitivity)
Covariates_train_L1 <- train_L1$Covariates
Survival_train_L1   <- train_L1$Survival
Censored_train_L1   <- train_L1$Censored
Covariates_test_L1  <- test_L1$Covariates
Survival_test_L1    <- test_L1$Survival
Censored_test_L1    <- test_L1$Censored
save(Covariates_train_L1, Survival_train_L1, Censored_train_L1,
     file = file.path(data_dir, "train_data_L1_clean.rda"))
save(Covariates_test_L1, Survival_test_L1, Censored_test_L1,
     file = file.path(data_dir, "test_data_L1_clean.rda"))

# --- L1 mapping per L2 subtype (for primary's structural prior) ---------
sub_to_grp <- tapply(subgroup, subtype, function(x) x[1])
l1_per_subtype <- sub_to_grp[subtypes]
stopifnot(!any(is.na(l1_per_subtype)))

# --- target_covs builders -----------------------------------------------
# Primary: 12 L2 cancer-types, structural prior pools β by L1 parent.
target_covs_primary <- lapply(target_cov_ids, function(cid) {
  list(cov_id          = cid,
       hierarchy_label = "L2_12subtype_pooled_by_L1_subgroup",
       group_per_cancer = unname(l1_per_subtype[subtypes]))
})

# Sensitivity: 4 L1 cancer-types, IDENTITY hierarchy (each L1 is its own bin).
# group_per_cancer length = n_cancer = 4 = subgroups (1:1 with cancer-type axis).
# This is what makes sensitivity meaningfully different from primary —
# coarser cancer-type unit + no further pooling.
target_covs_sensitivity <- lapply(target_cov_ids, function(cid) {
  list(cov_id          = cid,
       hierarchy_label = "L1_4subgroup_identity",
       group_per_cancer = subgroups)
})

save(target_covs_primary, subtypes, l1_per_subtype,
     file = file.path(data_dir, "target_covs_primary_clean.rda"))
save(target_covs_sensitivity, subgroups,
     file = file.path(data_dir, "target_covs_sensitivity_clean.rda"))

cat("\n=== PRIMARY (12 L2 → pool by L1) ===\n")
for (t in target_covs_primary) {
  cat(sprintf("  cov %s : bins = %s\n", t$cov_id,
              paste(unique(t$group_per_cancer), collapse = " / ")))
}
cat("\n=== SENSITIVITY (4 L1 → identity, no pooling) ===\n")
for (t in target_covs_sensitivity) {
  cat(sprintf("  cov %s : bins = %s\n", t$cov_id,
              paste(t$group_per_cancer, collapse = " / ")))
}

cat(sprintf("\nSaved clean train/test + target_covs (primary + L1 sensitivity).\n"))
