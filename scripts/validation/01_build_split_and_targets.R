# Build held-out split (20% test per cancer, stratified, seed=42) and the
# structural-prior target tables for the 4 pre-specified covariates.
#
# Outputs:
#   data/split_indices.rda — list per cancer of train/test integer indices
#   data/train_data.rda    — Covariates_train, Survival_train, Censored_train
#   data/test_data.rda     — Covariates_test, Survival_test, Censored_test
#   data/target_covs_primary.rda    — primary spec target_covs object
#   data/target_covs_secondary.rda  — secondary spec target_covs object

user_lib <- file.path(Sys.getenv("LOCALAPPDATA"), "R/win-library/4.6")
.libPaths(c(user_lib, .libPaths()))

run_dir <- "C:/FkCancer/runs/run_validation_2026-05-09"
dir.create(file.path(run_dir, "data"), showWarnings = FALSE, recursive = TRUE)

# Load full XYC data
load("C:/FkCancer/repos/HierarchicalSS_PanCanPanOmics/XYC_V2_WithAge_StandardizedPredictors.rda")
cancer_types <- names(Covariates)
n_cancer <- length(cancer_types)
cat(sprintf("Loaded XYC: %d cancers\n", n_cancer))

# 80/20 split with seed 42, stratified per cancer
set.seed(42L)
train_idx <- list()
test_idx  <- list()
for (c in seq_len(n_cancer)) {
  n_c <- nrow(Covariates[[c]])
  n_test <- max(1L, floor(0.20 * n_c))   # at least 1 held out
  perm <- sample.int(n_c)
  test_idx[[c]]  <- sort(perm[seq_len(n_test)])
  train_idx[[c]] <- sort(perm[(n_test + 1L):n_c])
}
names(train_idx) <- cancer_types
names(test_idx)  <- cancer_types

cat("Split summary (train_n / test_n):\n")
for (c in seq_len(n_cancer)) {
  cat(sprintf("  %-5s  train=%4d  test=%3d\n",
              cancer_types[c], length(train_idx[[c]]), length(test_idx[[c]])))
}

# Build train and test versions of the data
Covariates_train <- list(); Survival_train <- list(); Censored_train <- list()
Covariates_test  <- list(); Survival_test  <- list(); Censored_test  <- list()
for (c in seq_len(n_cancer)) {
  X_c <- Covariates[[c]]
  Y_c <- Survival[[c]]
  Cn_c <- Censored[[c]]
  Covariates_train[[c]] <- X_c[train_idx[[c]], , drop = FALSE]
  Survival_train[[c]]   <- Y_c[train_idx[[c]]]
  Censored_train[[c]]   <- Cn_c[train_idx[[c]]]
  Covariates_test[[c]]  <- X_c[test_idx[[c]], , drop = FALSE]
  Survival_test[[c]]    <- Y_c[test_idx[[c]]]
  Censored_test[[c]]    <- Cn_c[test_idx[[c]]]
}
names(Covariates_train) <- cancer_types
names(Covariates_test)  <- cancer_types
names(Survival_train)   <- cancer_types
names(Survival_test)    <- cancer_types
names(Censored_train)   <- cancer_types
names(Censored_test)    <- cancer_types

save(train_idx, test_idx,
     file = file.path(run_dir, "data/split_indices.rda"))
save(Covariates_train, Survival_train, Censored_train,
     file = file.path(run_dir, "data/train_data.rda"))
save(Covariates_test, Survival_test, Censored_test,
     file = file.path(run_dir, "data/test_data.rda"))
cat(sprintf("Saved split_indices.rda, train_data.rda, test_data.rda\n"))

# ---- Build target_covs for primary and secondary specifications ----
# Hierarchy assignments (mirroring hierarchies_long.csv from the screening run)
# Primary spec target covariates and their (hierarchy, granularity):
#   cov 10.1 -> epithelial_class / 2c_10bin
#   cov 1.1  -> epithelial_class / 2c_10bin
#   cov 8.2  -> epithelial_class / 2a_3bin   (PRIMARY)
#   cov 20.1 -> germ_layer / 3b_mid
# Secondary spec: same except cov 8.2 -> germ_layer / 3b_mid

# Hierarchy assignments per cancer (from cancer_type_hierarchies_2026-05-09.md)
# Returns a named character vector: cancer -> group, NA for contested/excluded
hier_2a_3bin <- c(
  ACC="Epithelial", BLCA="Epithelial", BRCA="Epithelial", CESC="Epithelial",
  CHOL="Epithelial", CORE="Epithelial", DLBC="Hematological", ESCA="Epithelial",
  HNSC="Epithelial", KICH="Epithelial", KIRC="Epithelial", KIRP="Epithelial",
  LGG="Non-epithelial", LIHC="Epithelial", LUAD="Epithelial", LUSC="Epithelial",
  MESO=NA_character_, OV="Epithelial", PAAD="Epithelial",
  PCPG="Non-epithelial", PRAD="Epithelial", SARC="Non-epithelial",
  SKCM="Non-epithelial", STAD="Epithelial", TGCT="Non-epithelial",
  THCA="Epithelial", THYM="Epithelial", UCEC="Epithelial",
  UCS=NA_character_  # contested (biphasic)
)

hier_2c_10bin <- c(
  ACC="Glandular adeno", BLCA="Urothelial", BRCA="Glandular adeno",
  CESC="Squamous", CHOL="Cholangiocarcinoma", CORE="Glandular adeno",
  DLBC="Other (germ/melan/neuroend/heme/thymic)", ESCA="Squamous",
  HNSC="Squamous", KICH="RCC", KIRC="RCC", KIRP="RCC",
  LGG="Glial", LIHC="HCC", LUAD="Glandular adeno", LUSC="Squamous",
  MESO=NA_character_, OV="Endometrioid+serous", PAAD="Glandular adeno",
  PCPG="Other (germ/melan/neuroend/heme/thymic)", PRAD="Glandular adeno",
  SARC="Mesench/mesothelial/sarcoma",
  SKCM="Other (germ/melan/neuroend/heme/thymic)", STAD="Glandular adeno",
  TGCT="Other (germ/melan/neuroend/heme/thymic)", THCA="Glandular adeno",
  THYM="Other (germ/melan/neuroend/heme/thymic)", UCEC="Endometrioid+serous",
  UCS=NA_character_   # contested
)

hier_3b_mid <- c(
  ACC="Intermediate mesoderm", BLCA=NA_character_, BRCA=NA_character_,
  CESC=NA_character_, CHOL="Foregut endoderm", CORE="Hindgut endoderm",
  DLBC="Hematopoietic mesoderm", ESCA="Foregut endoderm",
  HNSC=NA_character_, KICH="Intermediate mesoderm", KIRC="Intermediate mesoderm",
  KIRP="Intermediate mesoderm", LGG="Neural ectoderm",
  LIHC="Foregut endoderm", LUAD="Foregut endoderm", LUSC="Foregut endoderm",
  MESO="Lateral plate mesoderm", OV="Intermediate mesoderm",
  PAAD="Foregut endoderm", PCPG="Neural crest", PRAD=NA_character_,
  SARC="Paraxial+lateral mesoderm", SKCM="Neural crest",
  STAD="Foregut endoderm", TGCT=NA_character_, THCA="Foregut endoderm",
  THYM=NA_character_, UCEC="Lateral plate mesoderm",
  UCS="Lateral plate mesoderm"
)

stopifnot(all(cancer_types %in% names(hier_2a_3bin)))
stopifnot(all(cancer_types %in% names(hier_2c_10bin)))
stopifnot(all(cancer_types %in% names(hier_3b_mid)))

# target_covs format: list of {cov_id, hierarchy_label, group_per_cancer}
# group_per_cancer is a vector of length n_cancer in cancer_types order, NA for excluded
build_target <- function(cov_id, hier_label, group_vec) {
  list(
    cov_id = cov_id,
    hierarchy_label = hier_label,
    group_per_cancer = unname(group_vec[cancer_types])
  )
}

target_covs_primary <- list(
  build_target(10.1, "epithelial_class/2c_10bin", hier_2c_10bin),
  build_target(1.1,  "epithelial_class/2c_10bin", hier_2c_10bin),
  build_target(8.2,  "epithelial_class/2a_3bin",  hier_2a_3bin),
  build_target(20.1, "germ_layer/3b_mid",         hier_3b_mid)
)

target_covs_secondary <- list(
  build_target(10.1, "epithelial_class/2c_10bin", hier_2c_10bin),
  build_target(1.1,  "epithelial_class/2c_10bin", hier_2c_10bin),
  build_target(8.2,  "germ_layer/3b_mid",         hier_3b_mid),
  build_target(20.1, "germ_layer/3b_mid",         hier_3b_mid)
)

# Sanity check the structural-prior tables
cat("\nPRIMARY target_covs:\n")
for (t in target_covs_primary) {
  cat(sprintf("  cov %s @ %s: groups = %s\n",
              t$cov_id, t$hierarchy_label,
              paste(unique(t$group_per_cancer[!is.na(t$group_per_cancer)]),
                    collapse = " / ")))
  n_in <- sum(!is.na(t$group_per_cancer))
  n_out <- sum(is.na(t$group_per_cancer))
  cat(sprintf("    %d cancers in non-contested group, %d contested/excluded\n",
              n_in, n_out))
}

cat("\nSECONDARY target_covs:\n")
for (t in target_covs_secondary) {
  cat(sprintf("  cov %s @ %s: groups = %s\n",
              t$cov_id, t$hierarchy_label,
              paste(unique(t$group_per_cancer[!is.na(t$group_per_cancer)]),
                    collapse = " / ")))
  n_in <- sum(!is.na(t$group_per_cancer))
  n_out <- sum(is.na(t$group_per_cancer))
  cat(sprintf("    %d cancers in non-contested group, %d contested/excluded\n",
              n_in, n_out))
}

save(target_covs_primary, cancer_types,
     file = file.path(run_dir, "data/target_covs_primary.rda"))
save(target_covs_secondary, cancer_types,
     file = file.path(run_dir, "data/target_covs_secondary.rda"))
cat("\nSaved target_covs_primary.rda and target_covs_secondary.rda\n")
