# Paper 5 — build per-subgroup + per-subtype distribution and event-count table.
#
# Reads:  data/cavalli_clinical.rds
# Writes: reference/cavalli_subtype_distribution.csv
#
# Resolves the TBD numbers cited in PRE_REGISTRATION.md §13.
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
stopifnot(all(c("sample_id", "subgroup", "subtype",
                "os_time_years", "os_event") %in% names(clinical)))

with_surv <- !is.na(clinical$os_time_years) & !is.na(clinical$os_event)
cat(sprintf("With-survival n: %d / %d (%.1f%%); events = %d\n",
            sum(with_surv), nrow(clinical), 100 * mean(with_surv),
            sum(clinical$os_event[with_surv], na.rm = TRUE)))

# Helper: build one summary table for a grouping vector
summarize_by <- function(group_vec, level_label, has_surv, events, fu_years) {
  bins <- sort(unique(group_vec[!is.na(group_vec)]))
  do.call(rbind, lapply(bins, function(b) {
    rows <- which(group_vec == b)
    n_total <- length(rows)
    rows_surv <- rows[has_surv[rows]]
    n_surv <- length(rows_surv)
    n_events <- sum(events[rows_surv], na.rm = TRUE)
    median_fu <- if (n_surv > 0) median(fu_years[rows_surv], na.rm = TRUE) else NA_real_
    event_rate <- if (n_surv > 0) n_events / n_surv else NA_real_
    data.frame(level = level_label, bin = b,
               n = n_total, with_surv = n_surv,
               events = n_events,
               event_rate = event_rate,
               median_followup_years = median_fu,
               stringsAsFactors = FALSE)
  }))
}

agg_l1 <- summarize_by(as.character(clinical$subgroup), "L1_subgroup",
                       with_surv, clinical$os_event, clinical$os_time_years)
agg_l2 <- summarize_by(as.character(clinical$subtype), "L2_subtype",
                       with_surv, clinical$os_event, clinical$os_time_years)

cat("\nL1 subgroup distribution:\n")
print(agg_l1, row.names = FALSE)
cat("\nL2 subtype distribution:\n")
print(agg_l2, row.names = FALSE)

total_events <- sum(clinical$os_event[with_surv], na.rm = TRUE)
cat(sprintf("\nTotal with-survival n: %d, total events: %d, overall event rate: %.3f\n",
            sum(with_surv), total_events,
            mean(clinical$os_event[with_surv], na.rm = TRUE)))

# Power-aware framing check (PRE_REGISTRATION.md §13)
if (total_events < 85) {
  cat(sprintf("\n[POWER WARNING] Total events %d < 85 — pre-reg power-aware framing clause triggers; any 'matches' verdict will be labelled 'underpowered'.\n",
              total_events))
} else {
  cat(sprintf("\n[POWER OK] Total events %d ≥ 85 — pre-reg power-aware framing clause does NOT trigger.\n",
              total_events))
}

dist <- rbind(agg_l1, agg_l2)
out_path <- file.path(ref_dir, "cavalli_subtype_distribution.csv")
write.csv(dist, out_path, row.names = FALSE)
cat(sprintf("\nSaved %s\n", out_path))
