# Per-cancer log-likelihood breakdown.
#
# Aggregate per-patient LPPD by cancer type, for each spec. Compute per-cancer
# differences (primary - baseline, secondary - baseline) to identify whether
# structural priors help/hurt uniformly or vary by cancer.

run_dir <- "C:/FkCancer/runs/run_validation_2026-05-09"

read_pp <- function(spec) {
  fp <- file.path(run_dir, sprintf("loglik_per_patient_%s.csv", spec))
  if (!file.exists(fp)) return(NULL)
  read.csv(fp, stringsAsFactors = FALSE)
}

baseline  <- read_pp("baseline_train")
primary   <- read_pp("primary")
secondary <- read_pp("secondary")

if (is.null(baseline)) stop("baseline missing")

agg <- function(df, label) {
  if (is.null(df)) return(NULL)
  data.frame(
    cancer = unique(df$cancer),
    spec   = label,
    n_test = sapply(unique(df$cancer), function(c) sum(df$cancer == c)),
    n_censored = sapply(unique(df$cancer), function(c) sum(df$cancer == c & df$censored)),
    sum_lppd = sapply(unique(df$cancer), function(c) sum(df$lppd[df$cancer == c])),
    mean_lppd = sapply(unique(df$cancer), function(c) mean(df$lppd[df$cancer == c])),
    stringsAsFactors = FALSE
  )
}

per_cancer <- rbind(
  agg(baseline,  "baseline_train"),
  agg(primary,   "primary"),
  agg(secondary, "secondary"))
per_cancer <- per_cancer[!is.null(per_cancer), ]
write.csv(per_cancer, file.path(run_dir, "per_cancer_loglik_breakdown.csv"),
          row.names = FALSE)
cat(sprintf("Saved per_cancer_loglik_breakdown.csv (%d rows)\n", nrow(per_cancer)))

# Wide format for easy inspection
specs <- intersect(c("baseline_train", "primary", "secondary"), unique(per_cancer$spec))
cancers <- unique(per_cancer$cancer)
wide <- data.frame(cancer = cancers, stringsAsFactors = FALSE)
for (sp in specs) {
  sub <- per_cancer[per_cancer$spec == sp, c("cancer", "sum_lppd")]
  m <- match(wide$cancer, sub$cancer)
  wide[[sp]] <- sub$sum_lppd[m]
}
if (all(c("baseline_train", "primary") %in% specs)) {
  wide$delta_primary  <- wide$primary - wide$baseline_train
}
if (all(c("baseline_train", "secondary") %in% specs)) {
  wide$delta_secondary <- wide$secondary - wide$baseline_train
}
cat("\n=== Per-cancer LPPD by spec ===\n")
print(wide, row.names = FALSE, digits = 4)

write.csv(wide, file.path(run_dir, "per_cancer_loglik_wide.csv"), row.names = FALSE)
cat("\nSaved per_cancer_loglik_wide.csv\n")
