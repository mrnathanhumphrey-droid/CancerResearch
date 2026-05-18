# Paper 5 — diagnose LPPD sensitivity to test-fold rescaling.
#
# Hypothesis: the test scores have inflated sd (median 5.13× train_sd, up
# to 714× for some columns) due to projection through small singular
# values. If this scale mismatch is actively corrupting the LPPD
# comparison, re-standardizing X_test by its own fold stats (instead of
# train fold stats) should shift the LPPD substantially. If LPPD is
# stable, the scale issue is cosmetic.
#
# Tests three rescaling schemes:
#   (orig)  X_test as-is (standardized by train stats; current pipeline)
#   (own)   X_test z-scored by test fold stats (test_mean → 0, test_sd → 1)
#   (clip)  X_test as-is but clipped to |z|<5 per column
# Recomputes LPPD on the existing _clean Gibbs chains; β unchanged.

user_lib <- file.path(Sys.getenv("LOCALAPPDATA"), "R/win-library/4.6")
.libPaths(c(user_lib, .libPaths()))

repo_root <- "C:/Cancer Research/Paper5_Medulloblastoma_StructuralPrior"
data_dir  <- file.path(repo_root, "data")
res_dir   <- file.path(repo_root, "results")
ref_dir   <- file.path(repo_root, "reference")

X_all <- readRDS(file.path(data_dir, "bidifac_components_clean.rds"))
clin  <- readRDS(file.path(data_dir, "cavalli_clinical_aligned.rds"))
load(file.path(ref_dir, "paper5_split_indices.rda"))

# Reuse the same L2/L1 build pipeline as build_paper5_split_v2.R but with
# the alternative rescaling applied to test rows.
with_surv <- !is.na(clin$os_time_years) & !is.na(clin$os_event)
clin_s <- clin[with_surv, ]
X_s <- X_all[with_surv, , drop = FALSE]

subtype  <- as.character(clin_s$subtype)
subgroup <- as.character(clin_s$subgroup)
subtypes  <- sort(unique(subtype))
subgroups <- sort(unique(subgroup))

train_rows <- unlist(train_idx); test_rows <- unlist(test_idx)
train_sd  <- apply(X_s[train_rows, , drop = FALSE], 2, sd)
train_mn  <- colMeans(X_s[train_rows, , drop = FALSE])
test_sd   <- apply(X_s[test_rows,  , drop = FALSE], 2, sd)
test_mn   <- colMeans(X_s[test_rows,  , drop = FALSE])
ratio <- test_sd / train_sd
cat(sprintf("Sample sizes: train=%d test=%d\n", length(train_rows), length(test_rows)))
cat(sprintf("Scale ratio: median=%.3f max=%.3f min=%.3f\n",
            median(ratio), max(ratio), min(ratio)))

# === Rescaling variants applied to TEST rows only =======================
make_X_variant <- function(variant) {
  X_v <- X_s
  if (variant == "orig") return(X_v)
  if (variant == "own") {
    # Z-score test rows by test fold stats
    test_block <- X_v[test_rows, , drop = FALSE]
    test_block <- sweep(test_block, 2, test_mn, "-")
    test_block <- sweep(test_block, 2, test_sd, "/")
    X_v[test_rows, ] <- test_block
    return(X_v)
  }
  if (variant == "clip") {
    test_block <- X_v[test_rows, , drop = FALSE]
    test_block[test_block > 5]  <- 5
    test_block[test_block < -5] <- -5
    X_v[test_rows, ] <- test_block
    return(X_v)
  }
  if (variant == "drop_worst") {
    # Zero out columns with ratio > 5 (worst-scale ones).
    bad <- which(ratio > 5)
    cat(sprintf("  drop_worst: zeroing %d cols (ratio>5) in TEST rows only\n", length(bad)))
    X_v[test_rows, bad] <- 0
    return(X_v)
  }
  stop("unknown variant: ", variant)
}

# Build L2 and L1 test sets for a given X variant
build_L2_test <- function(X_var) {
  X_sv <- X_var   # X_var is already with-survival rows (was built from X_s)
  # age standardization, same as build_paper5_split_v2.R
  age_raw <- as.numeric(clin_s$age_at_dx)
  age_raw[is.na(age_raw)] <- median(age_raw, na.rm = TRUE)
  age_z <- (age_raw - mean(age_raw[unlist(train_idx)])) / sd(age_raw[unlist(train_idx)])

  Cov <- list(); Surv <- list(); Cens <- list()
  for (s in subtypes) {
    keep <- intersect(test_idx[[s]], which(subtype == s))
    X_sub <- cbind(intercept = 1, age = age_z[keep], X_sv[keep, , drop = FALSE])
    colnames(X_sub) <- c("0", "0.5", as.character(seq_len(ncol(X_sv))))
    Cov[[s]]  <- X_sub
    Surv[[s]] <- pmax(clin_s$os_time_years[keep], 1/365)
    Cens[[s]] <- 1L - as.integer(clin_s$os_event[keep])
  }
  list(Cov = Cov, Surv = Surv, Cens = Cens, axis = subtypes)
}
build_L1_test <- function(X_var) {
  X_sv <- X_var
  age_raw <- as.numeric(clin_s$age_at_dx)
  age_raw[is.na(age_raw)] <- median(age_raw, na.rm = TRUE)
  age_z <- (age_raw - mean(age_raw[unlist(train_idx)])) / sd(age_raw[unlist(train_idx)])

  Cov <- list(); Surv <- list(); Cens <- list()
  for (g in subgroups) {
    sub_in_g <- subtypes[
      vapply(subtypes, function(s) unique(subgroup[subtype == s])[1] == g,
             logical(1))
    ]
    keep <- unique(unlist(test_idx[sub_in_g]))
    keep <- intersect(keep, which(subgroup == g))
    X_sub <- cbind(intercept = 1, age = age_z[keep], X_sv[keep, , drop = FALSE])
    colnames(X_sub) <- c("0", "0.5", as.character(seq_len(ncol(X_sv))))
    Cov[[g]] <- X_sub
    Surv[[g]] <- pmax(clin_s$os_time_years[keep], 1/365)
    Cens[[g]] <- 1L - as.integer(clin_s$os_event[keep])
  }
  list(Cov = Cov, Surv = Surv, Cens = Cens, axis = subgroups)
}

# === LPPD compute using existing Gibbs chains ==========================
BURN  <- 50001L; THIN <- 400L; ITERS <- 1e5L
keep_iters <- seq(BURN, ITERS, by = THIN)
N_SAMPLES <- length(keep_iters)

loglik_per_patient_one_draw <- function(X_test, y_test, c_test, beta_c, sigma2) {
  mu <- as.numeric(X_test %*% beta_c)
  sigma <- sqrt(sigma2)
  ly <- log(y_test)
  uncens <- c_test == 0L
  ll <- numeric(length(y_test))
  ll[uncens]  <- dnorm(ly[uncens], mean = mu[uncens], sd = sigma, log = TRUE)
  ll[!uncens] <- pnorm(ly[!uncens], mean = mu[!uncens], sd = sigma,
                       lower.tail = FALSE, log.p = TRUE)
  ll
}

compute_lppd <- function(spec, test_data) {
  spec_dir <- file.path(res_dir, sprintf("gibbs_%s_clean", spec))
  chain_files <- list.files(spec_dir, pattern = "^chain_[1-4]\\.rda$",
                            full.names = TRUE)
  axis_names <- test_data$axis
  ll_acc <- list()
  for (a in axis_names) ll_acc[[a]] <- matrix(NA_real_,
                                              nrow = length(test_data$Surv[[a]]),
                                              ncol = 0)
  for (cf in chain_files) {
    env <- new.env(); load(cf, envir = env); out <- env$out
    for (a_idx in seq_along(axis_names)) {
      a <- axis_names[a_idx]
      X_test <- test_data$Cov[[a]]
      y_test <- test_data$Surv[[a]]; c_test <- test_data$Cens[[a]]
      cov_for_a <- as.numeric(colnames(X_test))
      ll_draws <- matrix(NA_real_, nrow = length(y_test), ncol = N_SAMPLES)
      for (d in seq_len(N_SAMPLES)) {
        it <- keep_iters[d]
        beta_full <- out$betas[[a_idx]][[it]]
        beta_c <- beta_full[as.character(cov_for_a)]
        if (any(is.na(beta_c))) beta_c[is.na(beta_c)] <- 0
        sigma2 <- out$sigma2[[it]]
        ll_draws[, d] <- loglik_per_patient_one_draw(X_test, y_test, c_test,
                                                    beta_c, sigma2)
      }
      ll_acc[[a]] <- cbind(ll_acc[[a]], ll_draws)
    }
  }
  total <- 0
  for (a in axis_names) {
    M <- ll_acc[[a]]
    log_mean <- apply(M, 1, function(v) {
      mx <- max(v); mx + log(mean(exp(v - mx)))
    })
    total <- total + sum(log_mean)
  }
  total
}

# === Run across variants ================================================
specs <- c("baseline_train", "primary", "sensitivity",
           "primary_drop_t1", "primary_drop_t2", "primary_drop_t3")
variants <- c("orig", "own", "clip", "drop_worst")

results <- list()
for (v in variants) {
  cat(sprintf("\n=== Variant: %s ===\n", v))
  X_var <- make_X_variant(v)
  td_L2 <- build_L2_test(X_var)
  td_L1 <- build_L1_test(X_var)
  for (sp in specs) {
    td <- if (sp == "sensitivity") td_L1 else td_L2
    lppd <- compute_lppd(sp, td)
    cat(sprintf("  %s: LPPD = %.3f\n", sp, lppd))
    results[[paste0(v, "_", sp)]] <- list(variant = v, spec = sp, lppd = lppd)
  }
}

# === Summary table ======================================================
df <- do.call(rbind, lapply(results, function(r) {
  data.frame(variant = r$variant, spec = r$spec, lppd = r$lppd,
             stringsAsFactors = FALSE)
}))
# Delta vs baseline within each variant
df$delta <- NA_real_
for (v in variants) {
  baseline <- df$lppd[df$variant == v & df$spec == "baseline_train"]
  df$delta[df$variant == v] <- df$lppd[df$variant == v] - baseline
}
out_path <- file.path(res_dir, "lppd_rescaling_diagnostic.csv")
write.csv(df, out_path, row.names = FALSE)
cat(sprintf("\nWrote %s\n\nSummary (LPPD by variant × spec):\n", out_path))

# Pivot for readability
pivot <- reshape(df[, c("variant", "spec", "delta")],
                 idvar = "spec", timevar = "variant", direction = "wide")
print(pivot, row.names = FALSE)

cat("\n=== Interpretation ===\n")
prim_shift <- df$delta[df$variant == "own" & df$spec == "primary"] -
              df$delta[df$variant == "orig" & df$spec == "primary"]
sens_shift <- df$delta[df$variant == "own" & df$spec == "sensitivity"] -
              df$delta[df$variant == "orig" & df$spec == "sensitivity"]
cat(sprintf("  primary    delta shift (own vs orig): %+.2f nats\n", prim_shift))
cat(sprintf("  sensitivity delta shift (own vs orig): %+.2f nats\n", sens_shift))
if (abs(prim_shift) < 2) {
  cat("  → primary delta robust to test-fold rescaling. +6.68 nats likely real.\n")
} else {
  cat("  → primary delta sensitive to test-fold rescaling. Suggests projection artifact.\n")
}
