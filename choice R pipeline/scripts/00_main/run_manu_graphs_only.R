# ============================================================
# run_manu_graphs_only.R
# Regenerates manuscript figures using the existing STEP5 folder.
# Pre-sets STEP5_OUT so the main script reuses the folder.
# Stats models will re-fit (same data → same results).
# ============================================================

setwd(file.path(PROJECT_ROOT, "choice R pipeline"))

# Auto-detect the most recent run unless the caller pre-sets STEP5_OUT.
# Hard-coding a timestamp here is how this file silently went stale before
# (it still pointed at STEP5_stats_20260808_123821 on 2026-08-18, ten days and
# several re-runs later).
if (!exists("STEP5_OUT") || !nzchar(STEP5_OUT) || !dir.exists(STEP5_OUT)) {
  .p <- file.path(PROJECT_ROOT, "choice R pipeline/output/STEP5_stats")
  .r <- list.dirs(.p, full.names = TRUE, recursive = FALSE)
  .r <- .r[grepl("^STEP5_stats_\\d{8}_\\d{6}$", basename(.r))]
  if (!length(.r)) stop("No STEP5_stats_* run folder found under ", .p)
  STEP5_OUT <- .r[which.max(file.mtime(.r))]
}
STEP5_OUT <- normalizePath(STEP5_OUT, mustWork = FALSE)
cat("Reusing run:", STEP5_OUT, "\n")

# Pre-load trial data from the reorganised output folder so the STATS script
# doesn't need to search for STEP2_output at the pipeline root.
.step2_latest <- sort(
  list.dirs(file.path(PROJECT_ROOT, "choice R pipeline/output/STEP2_output"),
            recursive = FALSE), decreasing = TRUE)[1]
if (!exists("trial_activity_summary", envir = .GlobalEnv, inherits = FALSE)) {
  .f <- file.path(.step2_latest, "trial_activity_summary.csv")
  if (file.exists(.f)) {
    trial_activity_summary <- readr::read_csv(.f, show_col_types = FALSE)
    assign("trial_activity_summary", trial_activity_summary, envir = .GlobalEnv)
    message("Pre-loaded trial_activity_summary (", nrow(trial_activity_summary), " rows)")
  }
}
if (!exists("trial_occupancy_long", envir = .GlobalEnv, inherits = FALSE)) {
  .f <- file.path(.step2_latest, "trial_occupancy_long.csv")
  if (file.exists(.f)) {
    trial_occupancy_long <- readr::read_csv(.f, show_col_types = FALSE)
    assign("trial_occupancy_long", trial_occupancy_long, envir = .GlobalEnv)
    message("Pre-loaded trial_occupancy_long (", nrow(trial_occupancy_long), " rows)")
  }
}
if (!exists("group_dynamics_summary", envir = .GlobalEnv, inherits = FALSE)) {
  .step2b_latest <- sort(
    list.dirs(file.path(PROJECT_ROOT, "choice R pipeline/output/STEP2b_output"),
              recursive = FALSE), decreasing = TRUE)[1]
  .f <- file.path(.step2b_latest, "group_dynamics_summary.csv")
  if (file.exists(.f)) {
    group_dynamics_summary <- readr::read_csv(.f, show_col_types = FALSE)
    assign("group_dynamics_summary", group_dynamics_summary, envir = .GlobalEnv)
    message("Pre-loaded group_dynamics_summary (", nrow(group_dynamics_summary), " rows)")
  }
}

source(file.path(PROJECT_ROOT, "choice R pipeline/scripts/01_pipeline_analysis/activity_analysis_STATS_choice_exp.R"))
