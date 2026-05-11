# Build 18 target_covs objects for the operator-composition extension.
# 3 covariates × 3 rules × 2 hierarchy pairs.
#
# Output:  data/extension_specs.rda  containing list `extension_specs` of 18
#          named entries, each a target_covs object suitable for the
#          composition sampler.

user_lib <- file.path(Sys.getenv("LOCALAPPDATA"), "R/win-library/4.6")
.libPaths(c(user_lib, .libPaths()))

run_dir <- "C:/FkCancer/runs/run_extension_2026-05-10"
dir.create(file.path(run_dir, "data"), showWarnings = FALSE, recursive = TRUE)

# Recover cancer types
load("C:/FkCancer/runs/run_validation_2026-05-09/data/train_data.rda")
cancer_types <- names(Covariates_train)
n_cancer <- length(cancer_types)
stopifnot(n_cancer == 29)

# --- Hierarchy assignments (one per cancer; NA for contested) ---
hier_2a_3bin <- c(
  ACC="Epithelial", BLCA="Epithelial", BRCA="Epithelial", CESC="Epithelial",
  CHOL="Epithelial", CORE="Epithelial", DLBC="Hematological", ESCA="Epithelial",
  HNSC="Epithelial", KICH="Epithelial", KIRC="Epithelial", KIRP="Epithelial",
  LGG="Non-epithelial", LIHC="Epithelial", LUAD="Epithelial", LUSC="Epithelial",
  MESO=NA_character_, OV="Epithelial", PAAD="Epithelial",
  PCPG="Non-epithelial", PRAD="Epithelial", SARC="Non-epithelial",
  SKCM="Non-epithelial", STAD="Epithelial", TGCT="Non-epithelial",
  THCA="Epithelial", THYM="Epithelial", UCEC="Epithelial",
  UCS=NA_character_)
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
  UCS=NA_character_)
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
  UCS="Lateral plate mesoderm")
hier_1b_11bin <- c(
  ACC="Endocrine", BLCA="Male repro+bladder", BRCA="Breast",
  CESC="Gynecological", CHOL="Hepatobiliary+pancreas", CORE="GI tract",
  DLBC="Soft+heme+H&N", ESCA="GI tract", HNSC="Soft+heme+H&N",
  KICH="Kidney", KIRC="Kidney", KIRP="Kidney", LGG="CNS",
  LIHC="Hepatobiliary+pancreas", LUAD="Lung+pleura", LUSC="Lung+pleura",
  MESO="Lung+pleura", OV="Gynecological", PAAD="Hepatobiliary+pancreas",
  PCPG="Endocrine", PRAD="Male repro+bladder", SARC="Soft+heme+H&N",
  SKCM="Skin", STAD="GI tract", TGCT="Male repro+bladder",
  THCA="Endocrine", THYM="Endocrine", UCEC="Gynecological",
  UCS="Gynecological")

stopifnot(all(cancer_types %in% names(hier_2a_3bin)))
stopifnot(all(cancer_types %in% names(hier_2c_10bin)))
stopifnot(all(cancer_types %in% names(hier_3b_mid)))
stopifnot(all(cancer_types %in% names(hier_1b_11bin)))

pairs <- list(
  A = list(h1_label = "epi/2a_3bin",  gpc_h1 = unname(hier_2a_3bin[cancer_types]),
           h2_label = "germ/3b_mid",  gpc_h2 = unname(hier_3b_mid[cancer_types])),
  B = list(h1_label = "epi/2c_10bin", gpc_h1 = unname(hier_2c_10bin[cancer_types]),
           h2_label = "tissue/1b_11bin", gpc_h2 = unname(hier_1b_11bin[cancer_types]))
)

target_cov_ids <- c(8.2, 10.1, 1.1)
rules <- 1:3
pair_keys <- c("A", "B")

extension_specs <- list()
spec_id <- 0L
for (cov in target_cov_ids) {
  for (r in rules) {
    for (pkey in pair_keys) {
      spec_id <- spec_id + 1L
      cov_str <- gsub("\\.", "_", as.character(cov))
      spec_name <- sprintf("spec_%02d_cov%s_rule%d_pair%s",
                           spec_id, cov_str, r, pkey)
      extension_specs[[spec_name]] <- list(
        spec_id = spec_id,
        spec_name = spec_name,
        cov_id = cov,
        rule = as.integer(r),
        pair_key = pkey,
        target_covs = list(list(
          cov_id = cov,
          rule = as.integer(r),
          pair = pairs[[pkey]]
        ))
      )
      tm_meta <- pairs[[pkey]]
      n_h1 <- sum(!is.na(tm_meta$gpc_h1))
      n_h2 <- sum(!is.na(tm_meta$gpc_h2))
      joint <- ifelse(is.na(tm_meta$gpc_h1) | is.na(tm_meta$gpc_h2),
                      NA_character_,
                      paste(tm_meta$gpc_h1, tm_meta$gpc_h2, sep="__"))
      n_jt <- sum(!is.na(joint))
      cat(sprintf("  %s: cov %s, RULE %d, %s x %s | h1 n=%d, h2 n=%d, joint n=%d\n",
                  spec_name, cov, r, tm_meta$h1_label, tm_meta$h2_label,
                  n_h1, n_h2, n_jt))
    }
  }
}

stopifnot(length(extension_specs) == 18L)
save(extension_specs, cancer_types,
     file = file.path(run_dir, "data/extension_specs.rda"))
cat(sprintf("\nSaved %s\n", file.path(run_dir, "data/extension_specs.rda")))

# Print the joint-group distribution per pair (so we can see how many singleton joint groups exist)
cat("\n=== Joint-group sizes per pair ===\n")
for (pkey in pair_keys) {
  joint <- ifelse(is.na(pairs[[pkey]]$gpc_h1) | is.na(pairs[[pkey]]$gpc_h2),
                  NA_character_,
                  paste(pairs[[pkey]]$gpc_h1, pairs[[pkey]]$gpc_h2, sep = "__"))
  tab <- table(joint[!is.na(joint)])
  cat(sprintf("\nPair %s (%s x %s): %d joint groups, sizes:\n",
              pkey, pairs[[pkey]]$h1_label, pairs[[pkey]]$h2_label, length(tab)))
  print(tab)
  cat(sprintf("  n_singleton_joint_groups: %d\n", sum(tab == 1)))
  cat(sprintf("  n_multi_joint_groups (>=2): %d\n", sum(tab >= 2)))
  cat(sprintf("  n_cancers in non-singleton joint groups: %d\n", sum(tab[tab >= 2])))
}
