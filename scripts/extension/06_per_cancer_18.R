# Per-cancer LPPD breakdown for each of the 18 extension specs.
#
# Output:
#   per_cancer_extension_long.csv    long format: cancer x spec x lppd x delta_paper1
#   per_cancer_extension_wide.csv    wide: rows=cancers, cols=specs (delta vs Paper 1 primary)

user_lib <- file.path(Sys.getenv("LOCALAPPDATA"), "R/win-library/4.6")
.libPaths(c(user_lib, .libPaths()))

run_dir <- "C:/FkCancer/runs/run_extension_2026-05-10"
val_dir <- "C:/FkCancer/runs/run_validation_2026-05-09"

load(file.path(run_dir, "data/extension_specs.rda"))
load(file.path(val_dir, "data/train_data.rda"))
cancer_types <- names(Covariates_train)

read_pp <- function(dir, spec) {
  fp <- file.path(dir, sprintf("loglik_per_patient_%s.csv", spec))
  if (!file.exists(fp)) return(NULL)
  read.csv(fp, stringsAsFactors = FALSE)
}

per_cancer_total <- function(pp) {
  if (is.null(pp)) return(NULL)
  agg <- aggregate(pp$lppd, by = list(cancer = pp$cancer), FUN = sum)
  setNames(agg$x, agg$cancer)
}

# Anchors
anchor_primary <- read_pp(val_dir, "primary")
anchor_baseline <- read_pp(val_dir, "baseline_train")
stopifnot(!is.null(anchor_primary), !is.null(anchor_baseline))
pc_primary <- per_cancer_total(anchor_primary)
pc_baseline <- per_cancer_total(anchor_baseline)

specs <- names(extension_specs)
long_rows <- list()
wide_mat <- matrix(NA_real_, length(cancer_types), length(specs),
                    dimnames = list(cancer_types, specs))

for (sp in specs) {
  pp <- read_pp(run_dir, sp)
  if (is.null(pp)) next
  pc <- per_cancer_total(pp)
  for (cancer in cancer_types) {
    lppd_spec <- if (cancer %in% names(pc)) pc[cancer] else NA_real_
    lppd_p1   <- if (cancer %in% names(pc_primary)) pc_primary[cancer] else NA_real_
    lppd_bl   <- if (cancer %in% names(pc_baseline)) pc_baseline[cancer] else NA_real_
    long_rows[[length(long_rows) + 1L]] <- data.frame(
      spec = sp,
      cov_id = extension_specs[[sp]]$cov_id,
      rule = extension_specs[[sp]]$rule,
      pair_key = extension_specs[[sp]]$pair_key,
      cancer = cancer,
      lppd_spec = unname(lppd_spec),
      lppd_paper1_primary = unname(lppd_p1),
      lppd_baseline_train = unname(lppd_bl),
      delta_vs_paper1 = unname(lppd_spec - lppd_p1),
      delta_vs_baseline = unname(lppd_spec - lppd_bl),
      stringsAsFactors = FALSE)
    if (cancer %in% names(pc) && cancer %in% names(pc_primary)) {
      wide_mat[cancer, sp] <- pc[cancer] - pc_primary[cancer]
    }
  }
}

if (length(long_rows) > 0) {
  long_df <- do.call(rbind, long_rows)
  write.csv(long_df, file.path(run_dir, "per_cancer_extension_long.csv"),
            row.names = FALSE)
  wide_df <- data.frame(cancer = rownames(wide_mat), wide_mat,
                         check.names = FALSE, stringsAsFactors = FALSE)
  write.csv(wide_df, file.path(run_dir, "per_cancer_extension_wide.csv"),
            row.names = FALSE)
  cat(sprintf("Saved per_cancer_extension_long.csv (%d rows) and per_cancer_extension_wide.csv\n",
              nrow(long_df)))
}
