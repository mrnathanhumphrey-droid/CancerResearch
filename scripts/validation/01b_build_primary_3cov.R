# Build target_covs_primary_3cov.rda — Paper 1 primary spec MINUS cov 10.1.
# Tests whether cov 10.1 is contributing real signal or fitting noise.
# (Reading A vs Reading B test from the 2026-05-10 follow-up brief.)

user_lib <- file.path(Sys.getenv("LOCALAPPDATA"), "R/win-library/4.6")
.libPaths(c(user_lib, .libPaths()))

run_dir <- "C:/FkCancer/runs/run_validation_2026-05-09"
load(file.path(run_dir, "data/target_covs_primary.rda"))

# target_covs_primary has 4 entries. Drop the cov_id == 10.1 entry.
target_covs_primary_3cov <- target_covs_primary[
  sapply(target_covs_primary, function(t) t$cov_id != 10.1)
]
stopifnot(length(target_covs_primary_3cov) == 3)
cat("Retained 3 covariates:\n")
for (t in target_covs_primary_3cov) {
  cat(sprintf("  cov %s @ %s\n", t$cov_id, t$hierarchy_label))
}

save(target_covs_primary_3cov, cancer_types,
     file = file.path(run_dir, "data/target_covs_primary_3cov.rda"))
cat(sprintf("Saved %s\n", file.path(run_dir, "data/target_covs_primary_3cov.rda")))
