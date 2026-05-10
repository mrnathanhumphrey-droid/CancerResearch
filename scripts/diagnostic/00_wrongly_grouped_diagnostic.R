# Wrongly-grouped diagnostic — Paper 2 prerequisite (b).
#
# For each (target covariate, hierarchy/granularity) cell that is touched by
# either Paper 1 or Paper 2, examine the per-cancer posterior beta from the
# Paper 1 baseline-on-train chains. Within each hierarchy group, compute:
#   - n_cancers (in group with the covariate present)
#   - mean of per-cancer baseline β posterior means
#   - sd of per-cancer baseline β posterior means (between-cancer-within-group)
#   - sign agreement (fraction of cancers with same sign as group mean)
#   - max-min range across cancers' posterior means
#
# A "wrongly grouped" cell has high within-group between-cancer variance and
# low sign agreement — i.e., the structural assignment forces together cancers
# whose baseline-fitted betas point in different directions.
#
# Output:
#   wrongly_grouped_diagnostic.csv  — per (cov, hier, gran, group) row
#   diag_within_group_summary.csv   — per (cov, hier, gran) aggregate
#   per_cancer_baseline_betas.csv   — per (cov, cancer) baseline posterior

user_lib <- file.path(Sys.getenv("LOCALAPPDATA"), "R/win-library/4.6")
.libPaths(c(user_lib, .libPaths()))

run_dir <- "C:/FkCancer/runs/run_extension_2026-05-10"
val_dir <- "C:/FkCancer/runs/run_validation_2026-05-09"

# Paper 1 baseline-on-train chains (per-cancer fitted betas)
chain_files <- file.path(val_dir, "gibbs_baseline_train",
                         sprintf("chain_%d.rda", 1:4))
stopifnot(all(file.exists(chain_files)))

# Recover covariate structure from training data
load(file.path(val_dir, "data/train_data.rda"))
covariates_by_cancer <- lapply(Covariates_train, function(c) as.numeric(colnames(c)))
cancer_types <- names(Covariates_train)
n_cancer <- length(cancer_types)

# Target covariates of interest
target_covs <- c(8.2, 10.1, 1.1, 20.1)

# Hierarchy assignments (mirror the Paper 1 build_split_and_targets.R)
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

hierarchies <- list(
  list(label = "epi/2a_3bin",   gpc = hier_2a_3bin),
  list(label = "epi/2c_10bin",  gpc = hier_2c_10bin),
  list(label = "germ/3b_mid",   gpc = hier_3b_mid),
  list(label = "tissue/1b_11bin", gpc = hier_1b_11bin)
)

# --- Extract per-cancer posterior means for each target cov ---
BURN <- 50001L
N_ITER <- 100001L
keep_iters <- BURN:N_ITER

cat(sprintf("Loading %d Paper-1 baseline-on-train chains and extracting per-cancer betas for target covs %s ...\n",
            length(chain_files), paste(target_covs, collapse=", ")))

# beta_samples[[cov_id_str]][[cancer]] = vector of post-burn samples
beta_samples <- list()
for (cov in target_covs) {
  beta_samples[[as.character(cov)]] <- vector("list", n_cancer)
  names(beta_samples[[as.character(cov)]]) <- cancer_types
}
for (j in seq_along(chain_files)) {
  cat(sprintf("  chain %d ...\n", j))
  e <- new.env(); load(chain_files[j], envir = e); o <- e$out
  # baseline_train uses post-burn iters of obj$betas (length 100000)
  # Note: betas list has 100000 entries (no +1) in this sampler version
  n_b <- length(o$betas[[1]])
  keep_b <- (n_b - (N_ITER - BURN) + 1L):n_b   # last 50,000 iters
  for (c in seq_len(n_cancer)) {
    cov_idx_in_c <- match(target_covs, covariates_by_cancer[[c]])
    bb <- o$betas[[c]]
    for (kk in seq_along(target_covs)) {
      if (is.na(cov_idx_in_c[kk])) next
      cov_str <- as.character(target_covs[kk])
      vec <- sapply(bb[keep_b], function(v) v[cov_idx_in_c[kk]])
      beta_samples[[cov_str]][[c]] <- c(beta_samples[[cov_str]][[c]], vec)
    }
  }
  rm(o, e); gc(verbose = FALSE)
}

# --- Per-cancer posterior summary ---
per_cancer_rows <- list()
for (cov in target_covs) {
  cov_str <- as.character(cov)
  for (c in seq_len(n_cancer)) {
    samp <- beta_samples[[cov_str]][[c]]
    if (length(samp) == 0) next
    per_cancer_rows[[length(per_cancer_rows) + 1L]] <- data.frame(
      cov_id = cov_str, cancer = cancer_types[c],
      n_samples = length(samp),
      beta_mean = mean(samp),
      beta_sd = sd(samp),
      beta_q025 = quantile(samp, 0.025),
      beta_q500 = quantile(samp, 0.5),
      beta_q975 = quantile(samp, 0.975),
      stringsAsFactors = FALSE
    )
  }
}
per_cancer_df <- do.call(rbind, per_cancer_rows)
write.csv(per_cancer_df,
          file.path(run_dir, "per_cancer_baseline_betas.csv"), row.names = FALSE)
cat(sprintf("Wrote per_cancer_baseline_betas.csv (%d rows)\n", nrow(per_cancer_df)))

# --- Within-group diagnostic per (cov, hier, group) ---
group_rows <- list()
agg_rows <- list()
for (cov in target_covs) {
  cov_str <- as.character(cov)
  cov_pc <- per_cancer_df[per_cancer_df$cov_id == cov_str, ]
  for (h in hierarchies) {
    label <- h$label
    gpc <- h$gpc[cancer_types]
    # Annotate
    cov_pc_h <- cov_pc
    cov_pc_h$group <- gpc[match(cov_pc_h$cancer, cancer_types)]
    cov_pc_h <- cov_pc_h[!is.na(cov_pc_h$group), ]
    if (nrow(cov_pc_h) == 0) next

    # per-(cov, hier) aggregate
    # within-group SD pooled vs between-group SD
    grand_mean <- mean(cov_pc_h$beta_mean)
    ss_total <- sum((cov_pc_h$beta_mean - grand_mean)^2)
    ss_between <- 0
    ss_within  <- 0
    n_groups <- 0
    n_singleton <- 0
    sign_agreement_total <- numeric(0)
    for (g in unique(cov_pc_h$group)) {
      sub <- cov_pc_h[cov_pc_h$group == g, ]
      n_g <- nrow(sub)
      n_groups <- n_groups + 1
      if (n_g == 1) n_singleton <- n_singleton + 1
      g_mean <- mean(sub$beta_mean)
      g_sd <- if (n_g > 1) sd(sub$beta_mean) else NA_real_
      sign_g <- sign(g_mean)
      same_sign <- sum(sign(sub$beta_mean) == sign_g & sign_g != 0)
      sign_agreement <- if (n_g > 0) same_sign / n_g else NA_real_
      sign_agreement_total <- c(sign_agreement_total, rep(sign_agreement, n_g))
      ss_between <- ss_between + n_g * (g_mean - grand_mean)^2
      if (n_g > 1) ss_within <- ss_within + sum((sub$beta_mean - g_mean)^2)
      group_rows[[length(group_rows) + 1L]] <- data.frame(
        cov_id = cov_str, hierarchy = label, group = g,
        n_cancers = n_g,
        group_mean = g_mean,
        within_group_sd = g_sd,
        within_group_range = if (n_g > 1) max(sub$beta_mean) - min(sub$beta_mean) else 0,
        sign_agreement = sign_agreement,
        cancer_list = paste(sub$cancer, collapse = ","),
        cancer_betas = paste(sprintf("%s=%.3f", sub$cancer, sub$beta_mean),
                             collapse = "; "),
        stringsAsFactors = FALSE
      )
    }
    eta_squared_baseline <- if (ss_total > 0) ss_between / ss_total else NA_real_
    pooled_within_sd <- if (n_groups < nrow(cov_pc_h))
      sqrt(ss_within / (nrow(cov_pc_h) - n_groups)) else NA_real_
    agg_rows[[length(agg_rows) + 1L]] <- data.frame(
      cov_id = cov_str, hierarchy = label,
      n_cancers_used = nrow(cov_pc_h),
      n_groups = n_groups,
      n_singleton_groups = n_singleton,
      eta_squared_baseline_betas = eta_squared_baseline,
      pooled_within_sd = pooled_within_sd,
      mean_sign_agreement = mean(sign_agreement_total),
      stringsAsFactors = FALSE)
  }
}
group_df <- do.call(rbind, group_rows)
agg_df <- do.call(rbind, agg_rows)
write.csv(group_df,
          file.path(run_dir, "wrongly_grouped_diagnostic.csv"), row.names = FALSE)
write.csv(agg_df,
          file.path(run_dir, "diag_within_group_summary.csv"), row.names = FALSE)

cat(sprintf("\nWrote wrongly_grouped_diagnostic.csv (%d rows)\n", nrow(group_df)))
cat(sprintf("Wrote diag_within_group_summary.csv (%d rows)\n", nrow(agg_df)))

cat("\n=== Per-(cov, hier) aggregate (sorted by mean_sign_agreement asc) ===\n")
print(agg_df[order(agg_df$mean_sign_agreement), ], row.names = FALSE, digits = 3)

cat("\n=== Cells with low sign agreement (<0.7) and >=2 non-singleton groups ===\n")
flagged <- agg_df[agg_df$mean_sign_agreement < 0.7 &
                  (agg_df$n_groups - agg_df$n_singleton_groups) >= 2, ]
print(flagged, row.names = FALSE, digits = 3)

cat("\n=== For each flagged cell, the worst within-group sign disagreement ===\n")
if (nrow(flagged) > 0) {
  for (i in seq_len(nrow(flagged))) {
    cov_id <- flagged$cov_id[i]; hier <- flagged$hierarchy[i]
    sub <- group_df[group_df$cov_id == cov_id & group_df$hierarchy == hier &
                    group_df$n_cancers > 1, ]
    sub <- sub[order(sub$sign_agreement), ]
    cat(sprintf("\n  cov %s @ %s — worst groups:\n", cov_id, hier))
    print(head(sub[, c("group", "n_cancers", "group_mean",
                       "within_group_sd", "sign_agreement", "cancer_betas")], 3),
          row.names = FALSE, digits = 3)
  }
}

cat("\nDone.\n")
