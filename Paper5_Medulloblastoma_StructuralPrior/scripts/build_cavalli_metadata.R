# Paper 5 — build per-subgroup + per-subtype distribution and event-count table.
#
# Reads:  data/cavalli_clinical.rds
# Writes: reference/cavalli_subtype_distribution.csv
#
# This file resolves the TBD numbers cited in PRE_REGISTRATION.md §13.
# Run AFTER fetch_cavalli_data.R has populated cavalli_clinical.rds.

user_lib <- file.path(Sys.getenv("LOCALAPPDATA"), "R/win-library/4.6")
.libPaths(c(user_lib, .libPaths()))

repo_root <- "C:/Cancer Research/Paper5_Medulloblastoma_StructuralPrior"
data_dir <- file.path(repo_root, "data")
ref_dir  <- file.path(repo_root, "reference")
dir.create(ref_dir, showWarnings = FALSE, recursive = TRUE)

clinical <- readRDS(file.path(data_dir, "cavalli_clinical.rds"))
cat(sprintf("Loaded clinical: %d rows × %d cols\n",
            nrow(clinical), ncol(clinical)))

# Required columns
stopifnot(all(c("sample_id", "subgroup", "subtype",
                "os_time_years", "os_event") %in% names(clinical)))

with_surv <- !is.na(clinical$os_time_years) & !is.na(clinical$os_event)
cat(sprintf("With-survival n: %d / %d (%.1f%%)\n",
            sum(with_surv), nrow(clinical), 100 * mean(with_surv)))

# --- Per-subgroup (L1) summary -------------------------------------------
agg_l1 <- aggregate(
  cbind(
    n = rep(1, nrow(clinical)),
    with_surv = as.integer(with_surv),
    events = ifelse(with_surv, clinical$os_event, NA_integer_)
  ),
  by = list(level = "L1_subgroup", bin = clinical$subgroup),
  FUN = function(x) sum(x, na.rm = TRUE)
)

# Per-subgroup median follow-up
median_fu_l1 <- aggregate(
  clinical$os_time_years[with_surv],
  by = list(level = "L1_subgroup", bin = clinical$subgroup[with_surv]),
  FUN = function(x) median(x, na.rm = TRUE)
)
names(median_fu_l1)[3] <- "median_followup_years"
agg_l1 <- merge(agg_l1, median_fu_l1, by = c("level", "bin"))

# --- Per-subtype (L2) summary --------------------------------------------
agg_l2 <- aggregate(
  cbind(
    n = rep(1, nrow(clinical)),
    with_surv = as.integer(with_surv),
    events = ifelse(with_surv, clinical$os_event, NA_integer_)
  ),
  by = list(level = "L2_subtype", bin = clinical$subtype),
  FUN = function(x) sum(x, na.rm = TRUE)
)

median_fu_l2 <- aggregate(
  clinical$os_time_years[with_surv],
  by = list(level = "L2_subtype", bin = clinical$subtype[with_surv]),
  FUN = function(x) median(x, na.rm = TRUE)
)
names(median_fu_l2)[3] <- "median_followup_years"
agg_l2 <- merge(agg_l2, median_fu_l2, by = c("level", "bin"))

dist <- rbind(agg_l1, agg_l2)
dist$event_rate <- dist$events / pmax(dist$with_surv, 1)

# Print summary
cat("\nL1 subgroup distribution:\n")
print(agg_l1, row.names = FALSE)
cat("\nL2 subtype distribution:\n")
print(agg_l2, row.names = FALSE)

cat(sprintf("\nTotal with-survival n: %d, total events: %d, overall event rate: %.3f\n",
            sum(with_surv),
            sum(clinical$os_event[with_surv], na.rm = TRUE),
            mean(clinical$os_event[with_surv], na.rm = TRUE)))

# Power-aware framing check (PRE_REGISTRATION.md §13)
total_events <- sum(clinical$os_event[with_surv], na.rm = TRUE)
if (total_events < 85) {
  cat(sprintf("\n[POWER WARNING] Total events %d < 85 — pre-reg power-aware framing clause triggers; any 'matches' verdict will be labelled 'underpowered'.\n",
              total_events))
} else {
  cat(sprintf("\n[POWER OK] Total events %d ≥ 85 — pre-reg power-aware framing clause does NOT trigger.\n",
              total_events))
}

out_path <- file.path(ref_dir, "cavalli_subtype_distribution.csv")
write.csv(dist, out_path, row.names = FALSE)
cat(sprintf("\nSaved %s\n", out_path))
