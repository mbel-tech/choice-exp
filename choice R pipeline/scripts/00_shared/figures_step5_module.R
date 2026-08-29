# =============================================================================
# figures_step5_module.R
#
# 2026-08-28 (RN review): axis labels carry their unit and spell out their acronym.
#   "Flow<->calm crossings" -> "... (N)"        RN #532: it is a count
#   "ALR(x vs calm)"        -> "ALR(x vs. calm)" RN #129: one form of vs. throughout
#   "Mean NND/IID (cm)"     -> spelled out       RN #334/#418: define before use
#   "School speed (cm/s)"   -> "(cm s-1)"        RN #321: one form for compound units
# =============================================================================
# Standalone, MODEL-FIT-FREE regeneration of the manuscript's STEP5 figures
# (Figure 7 zone/ALR -> "Figure_8", Figure 8 collective -> "Figure_9",
# Figure 9 interval -> "Figure_10", Figure 7bis). Companion to
# figures_seq_bout_module.R, which already does this for the sequence/bout
# figures; this is the same idea for STEP5's own figures.
#
# WHY: activity_analysis_STATS_choice_exp.R is monolithic -- it fits ~40
# models (several with N=1000 parametric-bootstrap LRTs) AND draws every
# figure in one ~50-minute pass. A purely cosmetic figure change (caption
# font size, star size, legend placement, panel order) previously forced a
# full re-fit just to see the result. Every data frame the figure code
# actually touches is pure data wrangling from THREE already-saved CSVs --
# no model object is required except for the two cell-mean panel captions
# (Figure 8 C/D), which are now served from a small cache
# (anova_statement_table_cells.csv, built by
# 00_shared/rebuild_figure_caption_table_cells.R) instead of a live res$anova
# table, closing the last gap.
#
# INPUTS (all already written by a completed STEP5 run -- point STEP5_OUT at
# it; nothing here re-fits anything):
#   STEP5_OUT/analysis_ready.csv                              (= df)
#   STEP2_output/*/trial_occupancy_long.csv                    (= df_long, via .find_latest_csv)
#   STEP2b_output/*/group_dynamics_summary.csv                 (= df_gd, via .find_latest_csv)
#   easy_scripts/standalone/anova_statement_table_trial.csv    (trial-level captions)
#   easy_scripts/standalone/anova_statement_table_interval.csv (interval captions)
#   easy_scripts/standalone/anova_statement_table_cells.csv    (cell-mean captions)
#   STEP5_OUT/<label_dir>/cld_treatmentxzone.csv               (CLD, per model)
#   STEP5_OUT/<label_dir>/cld_treatmentxtimepoint(_f).csv       (CLD, per model)
#   all_manu_graphs/sequence_and_patterns_graphs/*.png          (raster-embedded panels;
#                                                                 run figures_seq_bout_module.R first)
#
# OUT OF SCOPE: Figure S7 (Jacobs' D) is NOT regenerated here -- its caption
# is the one remaining panel family built directly from a live res$anova
# object (.mg_make_jacobs_panel(res_list=...)) with no CSV cache. Extending
# rebuild_figure_caption_table_cells.R's pattern to the four Jacobs models
# would close this the same way, if ever needed.
#
# USAGE:
#   STEP5_OUT <- file.path(PROJECT_ROOT, "choice R pipeline/output/STEP5_stats/STEP5_stats_<ts>")
#   source(file.path(PROJECT_ROOT, "choice R pipeline/scripts/00_shared/figures_step5_module.R"))
#   (or just source() it with STEP5_OUT unset -- it auto-finds the latest run)
#
# Runtime: seconds, not ~50 minutes -- confirmed 2026-08-09.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggtext)
  library(stringr)
  library(patchwork)
  library(png)
  library(readr)
})

`%||%` <- function(a, b) if (is.null(a)) b else a
ts_msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n")

# Canonical CLD letter rule (control/interval 1 == "a"), shared with
# export_tukey_interval_module.R so figure letters and exported letters agree.
source(file.path(file.path(PROJECT_ROOT, "choice R pipeline/scripts/00_shared"),
                 "cld_letter_policy.R"))

# ---- Locate the target STEP5 run ------------------------------------------
if (!exists("STEP5_OUT") || !nzchar(STEP5_OUT) || !dir.exists(STEP5_OUT)) {
  .step5_parent <- file.path(PROJECT_ROOT, "choice R pipeline/output/STEP5_stats")
  .runs <- list.dirs(.step5_parent, full.names = TRUE, recursive = FALSE)
  .runs <- .runs[grepl("^STEP5_stats_\\d{8}_\\d{6}$", basename(.runs))]
  # Only runs that reached a LATE model count as usable. Without this the module
  # will happily latch onto a run that is still being written and render panels
  # from half-present output -- which happened on 2026-08-17.
  .ok <- vapply(.runs, function(x)
    file.exists(file.path(x, "centroid_speed_timepoint", "anova.csv")), logical(1))
  if (any(!.ok))
    ts_msg("Ignoring ", sum(!.ok), " incomplete STEP5 run(s), newest: ",
           basename(.runs[!.ok][which.max(file.mtime(.runs[!.ok]))]))
  .runs <- .runs[.ok]
  if (!length(.runs)) stop("No COMPLETE STEP5_stats_* run found under ", .step5_parent)
  STEP5_OUT <- .runs[which.max(file.mtime(.runs))]
}
ts_msg("Target STEP5 run: ", STEP5_OUT)

.find_latest_csv_root <- file.path(PROJECT_ROOT, "choice R pipeline/output")
.find_latest_csv <- function(step_name, csv_filename) {
  parent <- file.path(.find_latest_csv_root, step_name)
  if (!dir.exists(parent)) return(NULL)
  subdirs <- list.dirs(parent, full.names = TRUE, recursive = FALSE)
  subdirs <- subdirs[grepl(paste0("^", step_name, "_\\d{8}_\\d{6}$"), basename(subdirs))]
  if (!length(subdirs)) return(NULL)
  latest <- subdirs[which.max(file.mtime(subdirs))]
  candidate <- file.path(latest, csv_filename)
  if (file.exists(candidate)) candidate else NULL
}

# ---- Reconstruct the plotting data frames (pure data wrangling; verbatim
#      logic from activity_analysis_STATS_choice_exp.R lines 245-476, minus
#      the parts already baked into analysis_ready.csv) -------------------
TREATMENT_LEVELS_g  <- c("control", "exercise choice")
DENSITY_LEVELS_g    <- c(4L, 8L, 12L, 16L)
DENSITY_AS_FACTOR_g <- TRUE
TIMEPOINT_LEVELS_g  <- 1:3

df <- readr::read_csv(file.path(STEP5_OUT, "analysis_ready.csv"), show_col_types = FALSE)
df <- df %>%
  dplyr::mutate(
    treatment      = factor(trimws(tolower(as.character(treatment))), levels = TREATMENT_LEVELS_g),
    fish_density_f = factor(fish_density, levels = DENSITY_LEVELS_g, ordered = isTRUE(DENSITY_AS_FACTOR_g)),
    timepoint      = as.integer(timepoint),
    timepoint_f    = factor(timepoint, levels = TIMEPOINT_LEVELS_g),
    trial_id       = as.integer(trial_id),
    tank           = as.character(tank)
  )
ts_msg("df (analysis_ready.csv): ", nrow(df), " rows")

.occ_path <- .find_latest_csv("STEP2_output", "trial_occupancy_long.csv")
if (is.null(.occ_path)) stop("trial_occupancy_long.csv not found under STEP2_output")
df_long <- readr::read_csv(.occ_path, show_col_types = FALSE) %>%
  dplyr::mutate(
    treatment   = factor(trimws(tolower(as.character(treatment))), levels = TREATMENT_LEVELS_g),
    timepoint   = as.integer(timepoint),
    timepoint_f = factor(timepoint, levels = TIMEPOINT_LEVELS_g),
    trial_id    = as.integer(trial_id),
    tank        = as.character(tank),
    zone_level  = factor(zone_level, levels = c("main", "sec")),
    zone        = factor(zone, levels = c("flow", "calm", "high", "medium", "low"))
  )
ts_msg("df_long (trial_occupancy_long.csv): ", nrow(df_long), " rows")

# ---- Safety net: define res_* fallback symbols BEFORE any figure code runs.
# Several panels guard construction on exists("res_xxx_tp") (e.g. Figure 10's
# NND/IID/area interval panels) -- since this module never fits models, those
# names must exist (as NULL) or the guard silently skips the whole panel.
# Defining them NULL is safe: every caption/CLD path here reads from a CSV
# cache first and only falls back to res_obj$anova when the cache is absent
# (R's lazy argument evaluation means that fallback is never forced when the
# cache is present).
for (.rn in c("res_cspd_agg","res_hull_agg","res_hull_tp","res_iid_agg","res_iid_tp",
              "res_nnd_agg","res_nnd_tp",
              "res_zone_flow_logit_agg","res_zone_main","res_zone_main_cells_agg",
              "res_zone_sec_cells_agg","res_zone_sec_high_agg")) {
  if (!exists(.rn, inherits = FALSE)) assign(.rn, NULL)
}
# =============================================================================
# ==== 2) BUILD LONG-FORMAT DATASETS ==========================================
# =============================================================================
# trial_occupancy_long has columns: trial_id, tank, treatment, timepoint,
#   zone_level (main/sec), zone (flow/calm/high/medium/low), prop_time.
# Split into main- and sec-zone subsets for separate fits.

# Ensure metadata columns exist on df_long; backfill from df where missing.
.meta_cols <- intersect(c("fish_density", "fish_density_f", "trial_date"), names(df))
if (length(.meta_cols) > 0 &&
    !all(.meta_cols %in% names(df_long))) {
  .meta <- dplyr::distinct(
    df[, c("trial_id", "timepoint", .meta_cols)],
    trial_id, timepoint, .keep_all = TRUE
  )
  df_long <- dplyr::left_join(df_long, .meta, by = c("trial_id", "timepoint"))
}
if ("fish_density" %in% names(df_long)) {
  df_long$fish_density_f <- factor(df_long$fish_density, levels = DENSITY_LEVELS_g)
} else {
  df_long$fish_density_f <- NA
}

df_main_long <- dplyr::filter(df_long, zone_level == "main") %>%
  dplyr::mutate(zone = factor(as.character(zone), levels = c("flow", "calm")))

df_sec_long <- dplyr::filter(df_long, zone_level == "sec") %>%
  dplyr::mutate(zone = factor(as.character(zone),
                              levels = c("high", "medium", "low", "calm")))

# Propagate phys_trial_id to long format
if ("phys_trial_id" %in% names(df)) {
  .phys_map <- dplyr::distinct(df[, c("trial_id", "phys_trial_id")])
  df_long <- dplyr::left_join(df_long, .phys_map, by = "trial_id")
  df_main_long <- dplyr::left_join(df_main_long, .phys_map, by = "trial_id")
  df_sec_long  <- dplyr::left_join(df_sec_long,  .phys_map, by = "trial_id")
}

ts_msg("Long-format: main=", nrow(df_main_long),
       " | sec=", nrow(df_sec_long))

# Aggregated long-format: average prop_time across timepoints per physical trial × zone.
# Physical trial = tank × trial_date × treatment (N = 16: 4 tanks × 2 dates × 2 treatments).
# Used for A1_agg / A2_agg (treatment × zone, no timepoint term).
.agg_long_cols <- c("phys_trial_id", "zone", "treatment", "tank",
                    "fish_density_f", "trial_date")
df_main_long_agg <- df_main_long %>%
  dplyr::mutate(phys_trial_id = as.integer(factor(paste(tank, trial_date, treatment)))) %>%
  dplyr::group_by(dplyr::across(dplyr::any_of(.agg_long_cols))) %>%
  dplyr::summarise(prop_time = mean(prop_time, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(zone = factor(as.character(zone), levels = levels(df_main_long$zone)))

df_sec_long_agg <- df_sec_long %>%
  dplyr::mutate(phys_trial_id = as.integer(factor(paste(tank, trial_date, treatment)))) %>%
  dplyr::group_by(dplyr::across(dplyr::any_of(.agg_long_cols))) %>%
  dplyr::summarise(prop_time = mean(prop_time, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(zone = factor(as.character(zone), levels = levels(df_sec_long$zone)))

# =============================================================================
# COMPOSITIONALLY CORRECT ZONE DATASETS
# =============================================================================
# Why: prop_flow + prop_calm = 1 (main zones) and sum(sub-zone props) = 1.
# Fitting beta GLMMs or Gaussian LMMs on the long-format stacked data treats
# zone rows within a session as independent — they are not (algebraic dependency).
# Solution:
#   Main zones  → scalar logit(prop_flow) per session / physical trial.
#   Sub-zones   → three log-ratios vs calm (log(prop_z / prop_calm)) per session /
#                 physical trial. Calm is the reference simplex component and is
#                 implicitly modelled. ILR basis would be formally correct but
#                 harder to interpret; log-ratios (pairwise vs calm) are sufficient
#                 when calm is a natural reference and the principal research
#                 question is directional (treatment pushes fish away from calm).
#   Zeros: a small Haldane-type offset (.eps_lr = 1e-4, ~0.01% of session time)
#          is added to all four sub-zone proportions before log-ratio computation
#          to stabilise the ratio when any zone is never used. No offset is needed
#          for the main-zone logit because rows with prop_flow ∈ {0,1} are removed
#          (they arise from sessions where fish never used one zone at all; n very
#          small, removing them prevents infinite logit and is documented below).

.eps_lr <- 1e-4   # Haldane offset for sub-zone log-ratios

# — Session-level main zone (N = 48) ——————————————————————————————————————————
df_main_wide <- tryCatch({
  df_main_long %>%
    dplyr::filter(zone %in% c("flow", "calm")) %>%
    dplyr::select(dplyr::any_of(c("trial_id", "phys_trial_id", "tank", "trial_date",
                                   "treatment", "timepoint", "timepoint_f",
                                   "fish_density", "fish_density_f", "zone", "prop_time"))) %>%
    tidyr::pivot_wider(names_from = zone, values_from = prop_time, names_prefix = "prop_") %>%
    dplyr::filter(is.finite(prop_flow)) %>%
    dplyr::mutate(
      prop_flow_sq = pmin(pmax(prop_flow, 1e-4), 1 - 1e-4),
      logit_flow   = log(prop_flow_sq / (1 - prop_flow_sq))
    )
}, error = function(e) { warning("df_main_wide failed: ", e$message); NULL })

# — Physical-trial-level main zone (N = 16) ———————————————————————————————————
df_main_wide_agg <- tryCatch({
  df_main_long_agg %>%
    dplyr::filter(zone %in% c("flow", "calm")) %>%
    dplyr::select(dplyr::any_of(c("phys_trial_id", "tank", "trial_date", "treatment",
                                   "fish_density_f", "zone", "prop_time"))) %>%
    tidyr::pivot_wider(names_from = zone, values_from = prop_time, names_prefix = "prop_") %>%
    dplyr::filter(is.finite(prop_flow) & prop_flow > 0 & prop_flow < 1) %>%
    dplyr::mutate(logit_flow = log(prop_flow / (1 - prop_flow)))
}, error = function(e) { warning("df_main_wide_agg failed: ", e$message); NULL })

# — Session-level sub-zones (N = 48) ——————————————————————————————————————————
df_sec_wide <- tryCatch({
  df_sec_long %>%
    dplyr::filter(zone %in% c("high", "medium", "low", "calm")) %>%
    dplyr::select(dplyr::any_of(c("trial_id", "phys_trial_id", "tank", "trial_date",
                                   "treatment", "timepoint", "timepoint_f",
                                   "fish_density", "fish_density_f", "zone", "prop_time"))) %>%
    tidyr::pivot_wider(names_from = zone, values_from = prop_time, names_prefix = "prop_") %>%
    dplyr::mutate(
      # DO NOT area-normalize here. The sub-zone proportions arrive from STEP2
      # section 6b ALREADY divided by their area constants (high/10, medium/18,
      # low/13, calm/42) and rescaled to sum to 1 (prop_time_ac_*), which is what
      # trial_occupancy_long carries as prop_time for zone_level == "sec".
      #
      # THE 2026-08-18 DOUBLE-AREA-NORMALISATION FIX
      #   Until 2026-08-18 this block divided by the area constants a SECOND
      #   time, so every sub-zone log-ratio carried a zone-specific offset of
      #   log(A_calm / A_zone): high +1.435, medium +0.847, low +1.173.
      #   Consequences of the double correction:
      #     - treatment, interval and treatment x interval tests WITHIN a zone
      #       were unaffected (a constant per zone cancels in every contrast);
      #     - the zero point did NOT mean "used exactly in proportion to area",
      #       so the sign of an emmean was not interpretable as preference;
      #     - the stacked zone-vs-zone contrasts (A2_pairs) WERE biased, by the
      #       differences between those offsets (e.g. high vs medium by 0.588);
      #     - the Haldane floor did not mean what its comment says: applied to
      #       prop/A it clamped at A x 1e-4 of occupancy (0.42% of session time
      #       for calm), not at 1e-4.
      #   The area correction belongs in STEP2 and is applied there once.
      prop_high_an   = pmax(prop_high,   .eps_lr),
      prop_medium_an = pmax(prop_medium, .eps_lr),
      prop_low_an    = pmax(prop_low,    .eps_lr),
      prop_calm_an   = pmax(prop_calm,   .eps_lr),
      lr_high   = log(prop_high_an   / prop_calm_an),
      lr_medium = log(prop_medium_an / prop_calm_an),
      lr_low    = log(prop_low_an    / prop_calm_an)
    )
}, error = function(e) { warning("df_sec_wide failed: ", e$message); NULL })

# — Physical-trial-level sub-zones (N = 16) ———————————————————————————————————
df_sec_wide_agg <- tryCatch({
  df_sec_long_agg %>%
    dplyr::filter(zone %in% c("high", "medium", "low", "calm")) %>%
    dplyr::select(dplyr::any_of(c("phys_trial_id", "tank", "trial_date", "treatment",
                                   "fish_density_f", "zone", "prop_time"))) %>%
    tidyr::pivot_wider(names_from = zone, values_from = prop_time, names_prefix = "prop_") %>%
    dplyr::mutate(
      # DO NOT area-normalize here. The sub-zone proportions arrive from STEP2
      # section 6b ALREADY divided by their area constants (high/10, medium/18,
      # low/13, calm/42) and rescaled to sum to 1 (prop_time_ac_*), which is what
      # trial_occupancy_long carries as prop_time for zone_level == "sec".
      #
      # THE 2026-08-18 DOUBLE-AREA-NORMALISATION FIX
      #   Until 2026-08-18 this block divided by the area constants a SECOND
      #   time, so every sub-zone log-ratio carried a zone-specific offset of
      #   log(A_calm / A_zone): high +1.435, medium +0.847, low +1.173.
      #   Consequences of the double correction:
      #     - treatment, interval and treatment x interval tests WITHIN a zone
      #       were unaffected (a constant per zone cancels in every contrast);
      #     - the zero point did NOT mean "used exactly in proportion to area",
      #       so the sign of an emmean was not interpretable as preference;
      #     - the stacked zone-vs-zone contrasts (A2_pairs) WERE biased, by the
      #       differences between those offsets (e.g. high vs medium by 0.588);
      #     - the Haldane floor did not mean what its comment says: applied to
      #       prop/A it clamped at A x 1e-4 of occupancy (0.42% of session time
      #       for calm), not at 1e-4.
      #   The area correction belongs in STEP2 and is applied there once.
      prop_high_an   = pmax(prop_high,   .eps_lr),
      prop_medium_an = pmax(prop_medium, .eps_lr),
      prop_low_an    = pmax(prop_low,    .eps_lr),
      prop_calm_an   = pmax(prop_calm,   .eps_lr),
      lr_high   = log(prop_high_an   / prop_calm_an),
      lr_medium = log(prop_medium_an / prop_calm_an),
      lr_low    = log(prop_low_an    / prop_calm_an)
    )
}, error = function(e) { warning("df_sec_wide_agg failed: ", e$message); NULL })

ts_msg("Wide-format datasets: main_wide=", if (!is.null(df_main_wide)) nrow(df_main_wide) else "NULL",
       " | main_wide_agg=", if (!is.null(df_main_wide_agg)) nrow(df_main_wide_agg) else "NULL",
       " | sec_wide=", if (!is.null(df_sec_wide)) nrow(df_sec_wide) else "NULL",
       " | sec_wide_agg=", if (!is.null(df_sec_wide_agg)) nrow(df_sec_wide_agg) else "NULL")

# Wide-format aggregated: one row per physical trial (tank × trial_date × treatment).
# Averages rate/proportion columns across timepoints; sums count/time columns.
# N = 16 (4 tanks × 2 recording dates × 2 treatments).
# Used for all A*_agg analyses (no timepoint term).
{
  .agg_sum_cols  <- intersect(c("n_main_switches", "total_time_s"), names(df))
  .agg_mean_cols <- intersect(c("switches_per_session", "prop_active",
                                 "zone_flux_per_session"), names(df))
  df_agg <- df %>%
    dplyr::mutate(phys_trial_id = as.integer(factor(paste(tank, trial_date, treatment)))) %>%
    dplyr::group_by(phys_trial_id, tank, trial_date, treatment, fish_density, fish_density_f) %>%
    dplyr::summarise(
      dplyr::across(dplyr::all_of(.agg_sum_cols),  ~ sum(.x,  na.rm = TRUE)),
      dplyr::across(dplyr::all_of(.agg_mean_cols), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      treatment      = factor(treatment,      levels = levels(df$treatment)),
      fish_density_f = factor(fish_density_f, levels = levels(df$fish_density_f))
    )
  ts_msg("df_agg: ", nrow(df_agg), " physical trials | ",
         dplyr::n_distinct(df_agg$treatment), " treatments")
}

# ---- Group dynamics (NND / IID / hull area / centroid speed): the four
#      collective indicators computable from an unlabelled per-frame position
#      set, and therefore unaffected by identity error. No alignment/heading
#      measure is loaded or plotted -- adapted from lines 3182-3257:
#      loads group_dynamics_summary.csv directly instead of via .get_global(),
#      since this module never runs the STEP2b step that populates that global. ----
.gd_path <- .find_latest_csv("STEP2b_output", "group_dynamics_summary.csv")
df_gd     <- if (!is.null(.gd_path)) readr::read_csv(.gd_path, show_col_types = FALSE) else NULL
df_gd_agg <- NULL

if (!is.null(df_gd) && nrow(df_gd) > 0) {
  .gd_want    <- c("trial_id", "trial", "tank", "treatment",
                   "fish_density", "fish_density_f", "timepoint",
                   "timepoint_f", "trial_date", "phys_trial_id")
  .gd_missing <- setdiff(.gd_want, names(df_gd))
  if (length(.gd_missing) > 0) {
    .gd_meta <- dplyr::distinct(
      df[, intersect(c("trial_id", "timepoint", .gd_missing), names(df))],
      trial_id, timepoint, .keep_all = TRUE
    )
    df_gd <- dplyr::left_join(df_gd, .gd_meta, by = c("trial_id", "timepoint"))
  }
  df_gd <- df_gd %>%
    dplyr::mutate(
      treatment      = factor(trimws(tolower(as.character(treatment))), levels = TREATMENT_LEVELS_g),
      timepoint_f    = factor(timepoint, levels = TIMEPOINT_LEVELS_g),
      fish_density_f = factor(fish_density, levels = DENSITY_LEVELS_g)
    )
  ts_msg("Group dynamics data: ", nrow(df_gd), " rows | ",
         dplyr::n_distinct(df_gd$trial_id), " sessions")

  .gd_resp_cols  <- intersect(
    c("mean_nnd_cm", "mean_iid_cm",
      "mean_hull_area_cm2", "mean_centroid_spd_cm"),
    names(df_gd))
  .gd_group_cols <- intersect(
    c("tank", "trial_date", "treatment", "fish_density", "fish_density_f"),
    names(df_gd))
  df_gd_agg <- df_gd %>%
    dplyr::mutate(phys_trial_id = as.integer(factor(paste(tank, trial_date, treatment)))) %>%
    dplyr::group_by(phys_trial_id, dplyr::across(dplyr::all_of(.gd_group_cols))) %>%
    dplyr::summarise(
      dplyr::across(dplyr::all_of(.gd_resp_cols), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      treatment      = factor(treatment,      levels = levels(df_gd$treatment)),
      fish_density_f = factor(fish_density_f, levels = levels(df_gd$fish_density_f))
    )
  ts_msg("df_gd_agg: ", nrow(df_gd_agg), " physical trials")
} else {
  ts_msg("group_dynamics_summary.csv not found -- NND/IID/area/speed panels will be skipped.")
}

# ============ get_global (169) ============
.get_global <- function(name, default = NULL) {
  if (exists(name, envir = .GlobalEnv, inherits = FALSE))
    get(name, envir = .GlobalEnv, inherits = FALSE)
  else default
}

# find_latest_csv (196): NOT redefined here. An earlier extraction pass had
# duplicated it with a getwd()-relative search (for the production script's
# "run from subdirectory" context), which silently shadowed the correct
# .find_latest_csv_root-based definition above for every call site after this
# point -- broke the new BOUT_output lookup below with no error, just a
# silently-NULL path (found & fixed 2026-08-09). Removed rather than fixed in
# place since the root-based version above is already correct and this
# module is never run from the production pipeline's own working directory.

# ============ LEVEL constants (184-194) ============
TREATMENT_LEVELS_g  <- .get_global("TREATMENT_LEVELS",  c("control", "exercise choice"))
DENSITY_LEVELS_g    <- .get_global("DENSITY_LEVELS",    c(4L, 8L, 12L, 16L))
DENSITY_AS_FACTOR_g <- .get_global("DENSITY_AS_FACTOR", TRUE)
TIMEPOINT_LEVELS_g  <- .get_global("TIMEPOINT_LEVELS",  1:3)

# Parametric-bootstrap LRT replicates for beta-GLMM treatment tests (df_method_memo
# section 4: KR/Satterthwaite are undefined for GLMMs; pbkrtest::PBmodcomp has no
# glmmTMB method, confirmed 2026-08-08, so this project implements the same idea
# manually via simulate.glmmTMB). Lower this for a faster exploratory run.
N_PB_DEFAULT <- .get_global("N_PB", 1000L)
PB_SEED      <- .get_global("PB_SEED", 20260808L)

# ============ fmt_F (655) ============
fmt_F <- function(x) {
  if (!is.finite(x)) return("NA")
  s <- format(signif(x, 4), scientific = FALSE, trim = TRUE)
  # Only trim trailing zeros when a decimal point is present -- otherwise an
  # integer denominator df like 30 is mangled to "3" (see project number-
  # formatting defect log).
  if (grepl(".", s, fixed = TRUE)) {
    s <- sub("0+$", "", s)
    s <- sub("\\.$", "", s)
  }
  s
}

# ============ fmt3 (667) ============
fmt3 <- function(x) {
  if (!is.finite(x)) return("NA")
  s <- sprintf("%.3f", x)
  if (grepl("0$", s)) s <- substr(s, 1, nchar(s) - 1)
  s
}

# ============ fmt_Fstat (679) ============
fmt_Fstat <- function(x) {
  if (!is.finite(x)) return("NA")
  if (abs(x) < 0.001) return("< 0.001")
  fmt_F(x)
}

# ============ fmt_es3 (684) ============
fmt_es3 <- function(x) {
  if (!is.finite(x)) return("NA")
  if (abs(x) < 0.001) return("< 0.001")
  fmt3(x)
}

# ============ fmt_p (690) ============
fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (round(p, 3) == 0) return("p < 0.001")
  paste0("p = ", fmt3(p))
}

# ============ CELL-MEAN pivots (2415-2480) ============
ts_msg("=== CELL-MEAN: zone x treatment cross-comparisons (main + sub) ===")

# --- Main-zone cell-mean long data (4 cells: flow|calm x 2 treatments) -------
df_main_cells_agg <- if (!is.null(df_main_wide_agg)) {
  df_main_wide_agg %>%
    dplyr::select(dplyr::any_of(c("phys_trial_id","tank","trial_date",
                                   "treatment","fish_density_f")),
                  prop_flow, prop_calm) %>%
    tidyr::pivot_longer(c(prop_flow, prop_calm),
                        names_to = "zone", values_to = "p",
                        names_prefix = "prop_") %>%
    dplyr::mutate(zone = factor(zone, levels = c("flow","calm")))
} else NULL

df_main_cells_tp <- if (!is.null(df_main_wide)) {
  df_main_wide %>%
    dplyr::select(dplyr::any_of(c("phys_trial_id","tank","trial_date",
                                   "treatment","timepoint_f","fish_density_f")),
                  prop_flow, prop_calm) %>%
    tidyr::pivot_longer(c(prop_flow, prop_calm),
                        names_to = "zone", values_to = "p",
                        names_prefix = "prop_") %>%
    dplyr::mutate(zone = factor(zone, levels = c("flow","calm")))
} else NULL

# --- Sub-zone CLR long data (8 cells: 4 zones x 2 treatments) ----------------
# CLR is the geometric-mean-centred log of the 4-part area-normalised
# composition. Computed within each composition (per trial for agg; per
# trial x timepoint for tp) so calm enters as a real cell with its own value.
.clr_long <- function(d, group_keys) {
  d %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_keys))) %>%
    dplyr::mutate(clr = log(p_an) - mean(log(p_an), na.rm = TRUE)) %>%
    dplyr::ungroup()
}

df_sec_cells_agg <- if (!is.null(df_sec_wide_agg)) {
  .d_agg <- df_sec_wide_agg %>%
    dplyr::select(dplyr::any_of(c("phys_trial_id","tank","trial_date",
                                   "treatment","fish_density_f")),
                  prop_high_an, prop_medium_an, prop_low_an, prop_calm_an) %>%
    tidyr::pivot_longer(c(prop_high_an, prop_medium_an, prop_low_an,
                          prop_calm_an),
                        names_to = "zone", values_to = "p_an") %>%
    dplyr::mutate(zone = factor(sub("^prop_(.+)_an$", "\\1", zone),
                                levels = c("high","medium","low","calm")))
  .clr_long(.d_agg, "phys_trial_id")
} else NULL

df_sec_cells_tp <- if (!is.null(df_sec_wide)) {
  .d_tp <- df_sec_wide %>%
    dplyr::select(dplyr::any_of(c("phys_trial_id","tank","trial_date",
                                   "treatment","timepoint_f","fish_density_f")),
                  prop_high_an, prop_medium_an, prop_low_an, prop_calm_an) %>%
    tidyr::pivot_longer(c(prop_high_an, prop_medium_an, prop_low_an,
                          prop_calm_an),
                        names_to = "zone", values_to = "p_an") %>%
    dplyr::mutate(zone = factor(sub("^prop_(.+)_an$", "\\1", zone),
                                levels = c("high","medium","low","calm")))
  .clr_long(.d_tp, c("phys_trial_id","timepoint_f"))
} else NULL

# Balance check
if (!is.null(df_sec_cells_agg)) {
  .bal_cells <- df_sec_cells_agg %>% dplyr::count(treatment, zone)
  ts_msg("  Cell-mean agg balance: ",
         paste(sprintf("%s/%s=%d", .bal_cells$treatment, .bal_cells$zone,
                       .bal_cells$n), collapse = "; "))
}

# ============ STYLE CONSTANTS + BASE_THEME + color maps (4200-4270) ============
.sg_dir <- file.path(STEP5_OUT, "skinny_graphs")
.sg_png <- file.path(.sg_dir, "png_preview")
dir.create(.sg_png, recursive = TRUE, showWarnings = FALSE)

# ---- Plot settings -----------------------------------------------------------
TREATMENT_COLORS <- c("control" = "#2166AC", "exercise choice" = "#D6604D")
TREATMENT_SHAPES <- c("control" = 17L,       "exercise choice" = 16L)
ZONE_COLORS_MAIN <- c("flow" = "#4DAC26",    "calm" = "#8073AC")

JITTER_W             <- 0.12
DODGE_W              <- 0.70
PT_ALPHA             <- 0.80
PT_SIZE              <- 2.4
MEAN_W               <- 0.40
LW_MEAN              <- 0.85
LW_ERR               <- 0.65
ERR_W                <- 0.13
CLD_SIZE             <- 4.0
CLD_NUDGE_SEM_MULT   <- 1.6
CLD_NUDGE_RANGE_FRAC <- 0.05

BASE_THEME <- ggplot2::theme_minimal(base_size = 13) +
  ggplot2::theme(
    panel.grid        = ggplot2::element_blank(),
    panel.background  = ggplot2::element_blank(),
    axis.line.x       = ggplot2::element_line(color = "black", linewidth = 0.85),
    axis.line.y       = ggplot2::element_line(color = "black", linewidth = 0.85),
    axis.ticks        = ggplot2::element_line(color = "black", linewidth = 0.7),
    axis.text.x       = ggplot2::element_text(size = 13, face = "bold", color = "black"),
    axis.text.y       = ggplot2::element_text(size = 11, color = "black"),
    axis.title        = ggplot2::element_text(size = 13, face = "bold", color = "black"),
    axis.title.x      = ggplot2::element_blank(),
    legend.title      = ggplot2::element_text(size = 11, face = "bold"),
    legend.text       = ggplot2::element_text(size = 10),
    strip.text        = ggplot2::element_text(size = 12, face = "bold"),
    strip.background  = ggplot2::element_blank(),
    plot.margin       = ggplot2::margin(8, 10, 16, 10),
    plot.caption      = ggtext::element_markdown(size = 14, hjust = 0,
                                                 face = "italic", lineheight = 1.3,
                                                 margin = ggplot2::margin(t = 10)),
    plot.subtitle     = ggplot2::element_text(size = 9,  hjust = 0, face = "italic",
                                              colour = "grey40")
  )

TIMEPOINT_LABELS <- c("1" = "5\u201325 min",
                       "2" = "45\u201365 min",
                       "3" = "85\u2013105 min")

# Zone color maps (user-defined)
color_map       <- c(flow = "red4", calm = "lightskyblue3", high = "tomato3",
                     medium = "mediumturquoise", low = "goldenrod3")
color_map_broad <- c(flow = "red4", calm = "lightskyblue3")

# Okabe palette for treatment-colored indicator line plots
.pal_tmp <- tryCatch(as.character(cols4all::c4a("okabe")),
                     error = function(e) c("#E69F00","#56B4E9","#009E73","#F0E442",
                                           "#0072B2","#D55E00","#CC79A7","#000000"))
TREATMENT_COLORS_PAL <- c("control" = .pal_tmp[6], "exercise choice" = .pal_tmp[5])
rm(.pal_tmp)

LINE_TP_LABELS <- c("1" = "5\u201325 min", "2" = "45\u201365 min", "3" = "85\u2013105 min")

# ---- Summary helpers ---------------------------------------------------------
.sem <- function(x) { x <- x[is.finite(x)]; if (length(x) < 2) NA_real_ else sd(x)/sqrt(length(x)) }

.smry_zone <- function(d, y_col, groups = c("zone","treatment")) {
  d %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(groups))) %>%
    dplyr::summarise(mean_y = mean(.data[[y_col]], na.rm=TRUE),
                     sem_y  = .sem(.data[[y_col]]),
                     .groups = "drop")
}

.smry_std <- function(d, y_col, groups = "treatment") {
  d %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(groups))) %>%
    dplyr::summarise(mean_y = mean(.data[[y_col]], na.rm=TRUE),
                     sem_y  = .sem(.data[[y_col]]),
                     .groups = "drop")
}

# ============ .make_std_plot (4374) ============
.make_std_plot <- function(df, y_col, y_label, title, caption_txt = NULL,
                            cld_df = NULL, is_sig = FALSE, facet_tp = FALSE,
                            n_txt = NULL, sig_textsize = 4) {
  groups_smry <- if (facet_tp) c("treatment","timepoint_f") else "treatment"
  smry <- .smry_std(df, y_col, groups_smry)

  # Map treatment to numeric x so annotations can use x = 1.5 (rvg-safe).
  .trt_lvls <- levels(factor(df$treatment))
  .n_trt    <- length(.trt_lvls)
  .trt_map  <- stats::setNames(seq_len(.n_trt), .trt_lvls)
  df_plot   <- df[is.finite(df[[y_col]]), ]
  df_plot$.x_pos <- unname(.trt_map[as.character(df_plot$treatment)])
  smry$.x_pos    <- unname(.trt_map[as.character(smry$treatment)])

  POS_JITTER <- ggplot2::position_jitter(width = JITTER_W, height = 0)

  p <- ggplot2::ggplot(df_plot, ggplot2::aes(
      x = .x_pos, y = .data[[y_col]],
      colour = treatment, shape = treatment)) +
    ggplot2::geom_point(position = POS_JITTER, alpha = PT_ALPHA, size = PT_SIZE) +
    ggplot2::geom_crossbar(
      data = smry,
      ggplot2::aes(x = .x_pos, y = mean_y, ymin = mean_y, ymax = mean_y),
      colour = "black",
      width = MEAN_W, linewidth = LW_MEAN, show.legend = FALSE) +
    ggplot2::geom_errorbar(
      data = smry,
      ggplot2::aes(x = .x_pos, y = mean_y, ymin = mean_y - sem_y, ymax = mean_y + sem_y),
      colour = "black",
      width = ERR_W, linewidth = LW_ERR, show.legend = FALSE) +
    ggplot2::scale_x_continuous(
      breaks = seq_len(.n_trt),
      labels = stringr::str_to_title(.trt_lvls),
      limits = c(0.5, .n_trt + 0.5),
      expand = c(0, 0)) +
    ggplot2::scale_colour_manual(values = TREATMENT_COLORS_PAL, name = "Treatment") +
    ggplot2::scale_shape_manual(values  = TREATMENT_SHAPES, name = "Treatment") +
    BASE_THEME +
    ggplot2::labs(title = title, x = "Treatment", y = y_label,
                  caption = .compose_caption(n_txt, caption_txt))

  # Significance bracket for two-group comparisons via ggsignif.
  # Bracket sits 10% of the data range above the highest plotted point.
  if (isTRUE(is_sig) && .n_trt == 2L) {
    .raw_y  <- df[[y_col]][is.finite(df[[y_col]])]
    .y_max  <- max(.raw_y, na.rm = TRUE)
    .y_rng  <- max(.y_max - min(.raw_y, na.rm = TRUE), 1e-10)
    .bar_y  <- .y_max + 0.10 * .y_rng

    .sig_df <- data.frame(xmin = 1, xmax = 2,
                          y_position = .bar_y, annotations = "*")
    p <- p +
      ggplot2::expand_limits(y = .y_max + 0.22 * .y_rng) +
      ggsignif::geom_signif(
        data        = .sig_df,
        mapping     = ggplot2::aes(xmin = xmin, xmax = xmax,
                                   y_position = y_position,
                                   annotations = annotations),
        manual      = TRUE,
        tip_length  = 0.02,
        size        = 0.5,
        textsize    = sig_textsize,
        vjust       = 0.5,
        color       = "black",
        inherit.aes = FALSE
      )
  }

  # Suppress y-axis ticks below zero (values are always non-negative).
  p <- p + ggplot2::scale_y_continuous(
    breaks = function(lims) { b <- pretty(lims); b[b >= 0] }
  )

  if (facet_tp)
    p <- p + ggplot2::facet_wrap(
      ~timepoint_f, ncol = 1,
      labeller = ggplot2::labeller(timepoint_f = TIMEPOINT_LABELS))
  p
}

# ============ .cap_all_sig (4568) ============
.cap_all_sig <- function(res) {
  if (is.null(res) || is.null(res$anova_caps)) return(NULL)
  caps <- vapply(res$anova_caps, function(v) {
    if (is.null(v) || (length(v) == 1L && is.na(v))) NA_character_ else as.character(v)
  }, character(1))
  caps <- caps[!is.na(caps) & nzchar(caps)]
  if (length(caps) == 0) return(NULL)
  paste(unname(caps), collapse = "\n")
}

# ============ .compose_caption (4579) ============
.compose_caption <- function(n_txt, anova_txt) {
  parts <- c(n_txt, anova_txt)
  parts <- parts[!is.null(parts)]
  parts <- parts[!is.na(parts) & nzchar(parts)]
  if (length(parts) == 0) NULL else paste(parts, collapse = "\n")
}

# ============ .is_sig (4591) ============
.is_sig <- function(res, key) {
  if (is.null(res) || is.null(res$anova_caps)) return(FALSE)
  val <- res$anova_caps[[key]]
  !is.null(val) && !(length(val) == 1L && is.na(val))
}

# ============ .cld_letters_differ (4598) ============
.cld_letters_differ <- function(cld_df) {
  if (is.null(cld_df) || !".group" %in% names(cld_df)) return(FALSE)
  u <- unique(trimws(as.character(cld_df$.group[!is.na(cld_df$.group)])))
  length(u) > 1L
}

# ============ .cap_zone (4662) ============
.cap_zone <- function(res, term_key, cld_df) {
  base <- .cap_all_sig(res)
  if (is.null(res) || is.null(res$anova) ||
      is.null(term_key) || !.cld_letters_differ(cld_df)) return(base)
  is_term_sig <- .is_sig(res, term_key)
  if (isTRUE(is_term_sig)) return(base)
  # Map the focal_terms key (e.g. "Treatment:Zone:Timepoint") to a regex
  # against the anova table's term column. Underscores / colons separate parts.
  pattern <- paste0(strsplit(term_key, ":", fixed = TRUE)[[1]],
                    collapse = ".*")
  stat_str <- .term_str_anyp(res$anova, pattern, res$stat_type)
  if (is.na(stat_str)) return(base)
  note <- paste0(.pretty_term(term_key),
                 ": ANOVA n.s. (", stat_str,
                 "); post-hoc CLD letters differ — see panel.")
  .compose_caption(base, note)
}

# ============ .n_cap (121) ============
.n_cap <- function(df, type = c("agg", "tp")) {
  type <- match.arg(type)
  if (is.null(df) || !"treatment" %in% names(df)) return(NULL)
  tryCatch({
    if (type == "agg") {
      ns <- df %>%
        dplyr::group_by(treatment) %>%
        dplyr::summarise(
          n = if ("phys_trial_id" %in% names(df))
                dplyr::n_distinct(phys_trial_id)
              else
                dplyr::n(),
          .groups = "drop")
      unit <- "Trial"
    } else {
      ns <- df %>%
        dplyr::group_by(treatment) %>%
        dplyr::summarise(
          n = dplyr::n_distinct(interaction(trial_id, timepoint_f, drop = TRUE)),
          .groups = "drop")
      unit <- "Trial × Timepoint combination"
    }
    n_uniq <- dplyr::n_distinct(ns$n)
    n_str  <- if (n_uniq == 1) {
      paste0("N = ", ns$n[1], " per group")
    } else {
      paste(sprintf("%s: N = %d", ns$treatment, ns$n), collapse = "; ")
    }
    paste0(n_str, "; one datapoint = one ", unit, ".")
  }, error = function(e) NULL)
}

# ============ mg_dir setup + MG_TAG_THEME + mg_frame + mg_treatment_xaxis + mg_save (5418-5510) ============
.mg_dir <- file.path(STEP5_OUT, "manu_graphs")
dir.create(.mg_dir, recursive = TRUE, showWarnings = FALSE)
.mg_dir_new     <- file.path(PROJECT_ROOT, "manu_graphs_26.05.2026")
.mg_dir_new_pdf <- file.path(.mg_dir_new, "PDF")
dir.create(.mg_dir_new,     recursive = TRUE, showWarnings = FALSE)
dir.create(.mg_dir_new_pdf, recursive = TRUE, showWarnings = FALSE)
ts_msg("  Dated output folder (PNG): ", .mg_dir_new)
ts_msg("  Dated output folder (PDF): ", .mg_dir_new_pdf)

# Common patchwork tag theme. Tag positioned in the top-right margin of each
# panel so that it sits inside the 10 mm margin without colliding with data.
.MG_TAG_THEME <- ggplot2::theme(
  plot.tag          = ggplot2::element_text(face = "bold", size = 18),
  plot.tag.position = "topleft"
)

# Helper: drop title, white background, generous top margin so the bold
# panel tag has room without colliding with axis labels.
.mg_frame <- function(p) {
  if (is.null(p)) return(NULL)
  p + ggplot2::theme(
    plot.title       = ggplot2::element_blank(),
    plot.margin      = ggplot2::margin(10, 10, 10, 10, unit = "mm"),
    plot.background  = ggplot2::element_rect(fill = "white", colour = NA)
  )
}

# Helper: tighten x-axis labels for treatment scatter panels. The underlying
# .make_std_plot uses scale_x_continuous with str_to_title labels ("Control",
# "Exercise Choice"); when panels are narrow the two-word "Exercise Choice"
# label crowds against "Control". This helper replaces the scale labels with
# a stacked two-line version so they fit at any sensible panel width.
.mg_treatment_xaxis <- function(p) {
  if (is.null(p)) return(NULL)
  p + ggplot2::scale_x_continuous(
        breaks = c(1, 2),
        labels = c("Control", "Exercise\nchoice"),
        limits = c(0.5, 2.5), expand = c(0, 0)
      ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 12, face = "bold",
                                          lineheight = 0.9, hjust = 0.5)
    )
}

# Save helper:
#   (1) PNG + PDF  → STEP5_OUT/manu_graphs/   (standard pipeline location)
#   (2) PNG        → D:/CHOICE R SCRIPTS/manu_graphs_26.05.2026/
#   (3) PDF (font-embedded, cairo_pdf device)
#              → D:/CHOICE R SCRIPTS/manu_graphs_26.05.2026/PDF/
.mg_save <- function(p, filename, width_mm, height_mm, dpi = 300) {
  if (is.null(p)) { ts_msg("  manu_graphs: skipped ", filename, " (NULL)"); return(invisible()) }

  # --- (1) Standard pipeline location: PNG + PDF ---
  for (.ext in c("png", "pdf")) {
    .out <- file.path(.mg_dir, paste0(filename, ".", .ext))
    tryCatch(
      ggplot2::ggsave(.out, p, width = width_mm, height = height_mm,
                      units = "mm", dpi = dpi, bg = "white"),
      error = function(e) warning(filename, " ", .ext, " save failed: ", e$message)
    )
  }

  # --- (2) Dated folder: PNG ---
  if (exists(".mg_dir_new") && dir.exists(.mg_dir_new)) {
    .out_png <- file.path(.mg_dir_new, paste0(filename, ".png"))
    tryCatch(
      ggplot2::ggsave(.out_png, p, width = width_mm, height = height_mm,
                      units = "mm", dpi = dpi, bg = "white"),
      error = function(e) warning(filename, " dated PNG save failed: ", e$message)
    )
  }

  # --- (3) Dated PDF subfolder: font-embedded PDF via cairo_pdf device ---
  if (exists(".mg_dir_new_pdf") && dir.exists(.mg_dir_new_pdf)) {
    .out_pdf <- file.path(.mg_dir_new_pdf, paste0(filename, ".pdf"))
    tryCatch(
      ggplot2::ggsave(.out_pdf, p,
                      width  = width_mm / 25.4,   # ggsave with cairo_pdf expects inches
                      height = height_mm / 25.4,
                      units  = "in",
                      device = cairo_pdf,
                      bg     = "white"),
      error = function(e) warning(filename, " dated PDF save failed: ", e$message)
    )
  }

  ts_msg("  manu_graphs: saved ", filename, " (.png + .pdf | dated .png + cairo .pdf)")
}

# Helper: add per-group (N = X) labels just above each x-axis tick.
# Works for 2-group scatter panels that use numeric x positions 1 and 2.
# Expands the lower y-limit to create visible padding below the data region.

# ============ .mg_add_n_labels (5511) ============
.mg_add_n_labels <- function(p, dat, y_col, tmap) {
  if (is.null(p) || is.null(dat) || !y_col %in% names(dat)) return(p)
  tryCatch({
    .n_df <- dat %>%
      dplyr::group_by(treatment) %>%
      dplyr::summarise(
        n = if ("phys_trial_id" %in% names(dat))
              dplyr::n_distinct(phys_trial_id) else dplyr::n(),
        .groups = "drop") %>%
      dplyr::mutate(
        .x    = unname(tmap[as.character(treatment)]),
        label = paste0("(N = ", n, ")"))
    .raw   <- dat[[y_col]][is.finite(dat[[y_col]])]
    .y_min <- min(.raw, na.rm = TRUE)
    .y_rng <- max(diff(range(.raw, na.rm = TRUE)), 1e-10)
    # Place label between x-axis line and y = 0 (or below y_min for panels
    # that already span negative values).
    if (.y_min >= 0) {
      .y_n      <- -0.08 * .y_rng
      .y_expand <- -0.18 * .y_rng
    } else {
      .y_n      <- .y_min - 0.10 * .y_rng
      .y_expand <- .y_min - 0.20 * .y_rng
    }
    .n_df$lbl_y <- .y_n
    p +
      ggplot2::expand_limits(y = .y_expand) +
      ggplot2::geom_text(data = .n_df,
        ggplot2::aes(x = .x, y = lbl_y, label = label),
        inherit.aes = FALSE, size = 3.0, colour = "grey35", vjust = 0.5,
        show.legend = FALSE)
  }, error = function(e) p)
}

# ============ .mg_make_alr_scatter (5582) ============
.mg_make_alr_scatter <- function(df_in, y_col, y_label, res_obj, n_cap_df = NULL,
                                  show_n = TRUE, cap_csv = NULL, sig_from_csv = NULL) {
  if (is.null(df_in) || !y_col %in% names(df_in)) return(NULL)
  tryCatch({
    .d    <- df_in[is.finite(df_in[[y_col]]), ]
    if (nrow(.d) == 0) return(NULL)
    .lvls <- levels(factor(.d$treatment))
    .tmap <- stats::setNames(seq_along(.lvls), .lvls)
    .d$.x <- unname(.tmap[as.character(.d$treatment)])
    .smry <- .d %>%
      dplyr::group_by(treatment) %>%
      dplyr::summarise(mean_y = mean(.data[[y_col]], na.rm = TRUE),
                       se_y   = .sem(.data[[y_col]]),
                       .groups = "drop") %>%
      dplyr::mutate(.x = unname(.tmap[as.character(treatment)]))
    .n_df <- if (!show_n) NULL else if (is.null(n_cap_df)) df_in else n_cap_df
    .cap  <- if (!is.null(cap_csv)) .compose_caption(.n_cap(.n_df, "agg"), cap_csv) else
             .compose_caption(.n_cap(.n_df, "agg"), .cap_all_sig(res_obj))

    p <- ggplot2::ggplot(.d, ggplot2::aes(x = .x, y = .data[[y_col]],
                                           colour = treatment, shape = treatment)) +
      ggplot2::geom_point(
        position = ggplot2::position_jitter(width = JITTER_W, height = 0),
        alpha = PT_ALPHA, size = PT_SIZE) +
      ggplot2::geom_crossbar(data = .smry,
        ggplot2::aes(x = .x, y = mean_y, ymin = mean_y, ymax = mean_y),
        colour = "black", width = MEAN_W, linewidth = LW_MEAN, show.legend = FALSE) +
      ggplot2::geom_errorbar(data = .smry,
        ggplot2::aes(x = .x, y = mean_y, ymin = mean_y - se_y, ymax = mean_y + se_y),
        colour = "black", width = ERR_W, linewidth = LW_ERR, show.legend = FALSE) +
      ggplot2::scale_x_continuous(
        breaks = seq_along(.lvls), labels = stringr::str_to_title(.lvls),
        limits = c(0.5, length(.lvls) + 0.5), expand = c(0, 0)) +
      ggplot2::scale_colour_manual(values = TREATMENT_COLORS_PAL, name = "Treatment") +
      ggplot2::scale_shape_manual(values  = TREATMENT_SHAPES,     name = "Treatment") +
      BASE_THEME +
      ggplot2::labs(x = "Treatment", y = y_label, caption = .cap) +
      ggplot2::theme(legend.position = "bottom")

    # Asterisk significance bracket for 2-group comparisons via ggsignif.
    # Use CSV-parsed significance when provided, otherwise fall back to LMM result.
    .is_sig_val <- if (!is.null(sig_from_csv)) sig_from_csv else .is_sig(res_obj, "Treatment")
    if (.is_sig_val && length(.lvls) == 2L) {
      .raw  <- .d[[y_col]]
      .ym   <- max(.raw, na.rm = TRUE)
      .rng  <- max(diff(range(.raw, na.rm = TRUE)), 1e-10)
      .bary <- .ym + 0.10 * .rng
      .sig_df <- data.frame(xmin = 1, xmax = 2,
                            y_position = .bary, annotations = "*")
      p <- p +
        ggplot2::expand_limits(y = .ym + 0.22 * .rng) +
        ggsignif::geom_signif(
          data        = .sig_df,
          mapping     = ggplot2::aes(xmin = xmin, xmax = xmax,
                                     y_position = y_position,
                                     annotations = annotations),
          manual      = TRUE,
          tip_length  = 0.02,
          size        = 0.5,
          textsize    = 8,
          vjust       = 0.5,
          color       = "black",
          inherit.aes = FALSE
        )
    }
    # Per-group (N = X) labels above x-axis ticks — only when show_n = TRUE.
    if (show_n) p <- .mg_add_n_labels(p, .d, y_col, .tmap)
    p
  }, error = function(e) { warning("alr scatter '", y_col, "' failed: ", e$message); NULL })
}

# ============ .mg_wrap_caption forward-decl (5817-5828) ============
# ---- Cell-mean panel helpers (used by new Figure 7) -------------------------
# Helper: wrap caption text to a per-panel character width (forward-declared
# here because .mg_make_cells_panel below uses it; the canonical definition
# at line ~5234 is identical).
if (!exists(".mg_wrap_caption", inherits = FALSE)) {
  .mg_wrap_caption <- function(text, panel_width_chars = 26) {
    # No longer wraps (2026-08-09): the ANOVA caption must render on one
    # line. panel_width_chars is kept as a no-op parameter so call sites
    # don't need editing; text is returned unchanged.
    text
  }
}

# ============ .mg_read_cld_treatxzone (5832) ============
.mg_read_cld_treatxzone <- function(label_dir) {
  fp <- file.path(STEP5_OUT, label_dir, "cld_treatmentxzone.csv")
  if (!file.exists(fp)) return(NULL)
  d <- tryCatch(read.csv(fp, stringsAsFactors = FALSE),
                error = function(e) NULL)
  if (is.null(d) || !all(c("treatment","zone",".group") %in% names(d)))
    return(NULL)
  d$.group    <- substr(trimws(as.character(d$.group)), 1, 2)
  d$treatment <- trimws(as.character(d$treatment))
  d$zone      <- trimws(as.character(d$zone))
  d <- d[nzchar(d$.group), , drop = FALSE]
  if (nrow(d) == 0) return(NULL)
  d
}

# ============ .fig7_anova_cap (5849) ============
.fig7_anova_cap <- function(res) {
  if (is.null(res) || is.null(res$anova)) return(NULL)
  av  <- res$anova
  row <- av[grepl("treatment.*zone$|^zone.*treatment$|treatment:zone|zone:treatment",
                   av$term, ignore.case = TRUE, perl = TRUE), ]
  if (nrow(row) == 0L) return(NULL)
  row <- row[1L, ]
  p   <- suppressWarnings(as.numeric(row$p_value))
  if (!is.finite(p)) return(NULL)
  p_str <- fmt_p(p)
  chi   <- fmt_Fstat(as.numeric(row$chisq))
  st    <- attr(av, "stat_type") %||% "Wald-chisq"
  is_f  <- isTRUE(grepl("^F", st))
  if (is_f) {
    df1 <- as.integer(row$df)
    df2 <- suppressWarnings(as.numeric(
             if (!is.null(row$df_denom)) row$df_denom else NA_real_))
    if (is.finite(df2)) {
      eta <- fmt_es3((as.numeric(row$chisq) * df1) / (as.numeric(row$chisq) * df1 + df2))
      sprintf("Treatment × Zone: F<sub>%d, %s</sub> = %s, %s, &eta;<sup>2</sup><sub>p</sub> = %s",
              df1, fmt_F(df2), chi, p_str, eta)
    } else
      sprintf("Treatment × Zone: F<sub>%d</sub> = %s, %s", df1, chi, p_str)
  } else {
    sprintf("Treatment × Zone: χ²<sub>%d</sub> = %s, %s",
            as.integer(row$df), chi, p_str)
  }
}

# ============ .mg_remap_cld_zone (5881) ============
.mg_remap_cld_zone <- function(cld_df, zones_order,
                                treatments_order = c("control","exercise choice")) {
  if (is.null(cld_df) || !".group" %in% names(cld_df)) return(cld_df)
  prec_idx <- order(
    match(cld_df$zone,      zones_order),
    match(cld_df$treatment, treatments_order)
  )
  seen <- character(0)
  for (i in prec_idx) {
    lets <- strsplit(cld_df$.group[i], "")[[1]]
    for (l in lets) if (!l %in% seen) seen <- c(seen, l)
  }
  if (!length(seen)) return(cld_df)
  letter_map <- stats::setNames(letters[seq_along(seen)], seen)
  cld_df$.group <- vapply(cld_df$.group, function(g) {
    lets     <- strsplit(g, "")[[1]]
    remapped <- letter_map[lets]
    remapped <- remapped[!is.na(remapped)]
    if (!length(remapped)) return(g)
    paste(sort(remapped), collapse = "")
  }, character(1L), USE.NAMES = FALSE)
  cld_df
}

# ============ .mg_make_cells_panel (5909) ============
.mg_make_cells_panel <- function(df_long, y_col, y_label, zones_order,
                                  label_dir, res_obj,
                                  hline_at = NULL,
                                  panel_width_chars = 30,
                                  treatments_order = c("control","exercise choice"),
                                  remap_cld = FALSE,
                                  always_show_cld = FALSE,
                                  show_zone_legend = TRUE,
                                  show_treatment_legend = TRUE,
                                  cap_csv = NULL) {
  if (is.null(df_long) || !y_col %in% names(df_long) ||
      !all(c("treatment","zone") %in% names(df_long))) return(NULL)
  tryCatch({
    .d <- df_long
    .d$treatment <- factor(as.character(.d$treatment), levels = treatments_order)
    .d$zone      <- factor(as.character(.d$zone),      levels = zones_order)
    .d <- .d[is.finite(.d[[y_col]]) & !is.na(.d$treatment) & !is.na(.d$zone), ]
    if (nrow(.d) == 0) return(NULL)

    # Build cell coordinates: zones on integer x; treatments dodged ±0.18.
    .n_trts <- length(treatments_order)
    .dodge  <- 0.36 / .n_trts
    .tx_off <- stats::setNames(
      ((seq_along(treatments_order) - mean(seq_along(treatments_order))) * 2) * .dodge,
      treatments_order)
    .d$cell_x <- as.integer(.d$zone) + unname(.tx_off[as.character(.d$treatment)])

    .smry <- .d %>%
      dplyr::group_by(treatment, zone) %>%
      dplyr::summarise(mean_y = mean(.data[[y_col]], na.rm = TRUE),
                       sem_y  = .sem(.data[[y_col]]),
                       max_y  = max(.data[[y_col]], na.rm = TRUE),
                       .groups = "drop") %>%
      dplyr::mutate(cell_x = as.integer(factor(zone, levels = zones_order)) +
                              unname(.tx_off[as.character(treatment)]))

    ymin <- min(c(.smry$mean_y - .smry$sem_y, .d[[y_col]]), na.rm = TRUE)
    ymax <- max(c(.smry$mean_y + .smry$sem_y, .d[[y_col]]), na.rm = TRUE)
    yrange <- ymax - ymin
    if (!is.finite(yrange) || yrange == 0) yrange <- 1

    p <- ggplot2::ggplot(.d, ggplot2::aes(x = cell_x, y = .data[[y_col]],
                                          colour = zone, shape = treatment))
    if (!is.null(hline_at))
      p <- p + ggplot2::geom_hline(yintercept = hline_at, linetype = "dashed",
                                   colour = "grey50", linewidth = 0.4)
    p <- p +
      ggplot2::geom_point(
        position = ggplot2::position_jitter(width = 0.05, height = 0),
        alpha = PT_ALPHA, size = PT_SIZE) +
      ggplot2::geom_crossbar(
        data = .smry,
        ggplot2::aes(x = cell_x, y = mean_y, ymin = mean_y, ymax = mean_y),
        colour = "black", width = 0.20, linewidth = LW_MEAN, show.legend = FALSE,
        inherit.aes = FALSE) +
      ggplot2::geom_errorbar(
        data = .smry,
        ggplot2::aes(x = cell_x, y = mean_y, ymin = mean_y - sem_y,
                     ymax = mean_y + sem_y),
        colour = "black", width = 0.10, linewidth = LW_ERR, show.legend = FALSE,
        inherit.aes = FALSE) +
      ggplot2::scale_x_continuous(
        breaks = seq_along(zones_order),
        labels = stringr::str_to_title(zones_order),
        limits = c(0.5, length(zones_order) + 0.5),
        expand = c(0, 0)) +
      ggplot2::scale_colour_manual(values = color_map[zones_order], name = "Zone",
        guide = if (show_zone_legend) "legend" else "none") +
      ggplot2::scale_shape_manual(values  = TREATMENT_SHAPES,    name = "Treatment",
        guide = if (show_treatment_legend) "legend" else "none") +
      BASE_THEME +
      ggplot2::labs(x = "Zone", y = y_label,
                    caption = .mg_wrap_caption(
                      if (!is.null(cap_csv)) cap_csv else .fig7_anova_cap(res_obj),
                      panel_width_chars)) +
      ggplot2::theme(
        # 14 -> 11 on 2026-08-09: these cell-mean panels (Figure 8C/8D) carry
        # the treatment x zone caption, and once eta^2_p joined it on a single
        # line the sub-zone version ("Treatment x Zone: F3, 42 = 8.131,
        # p < 0.001, eta^2_p = 0.367", ~55 chars) ran past the panel edge and
        # lost the effect size. Confirmed by inspecting the rendered PNG.
        plot.caption = ggtext::element_markdown(size = 11, hjust = 0,
                                                face = "italic", lineheight = 1.3,
                                                margin = ggplot2::margin(t = 10))
      )

    # CLD letters per cell (one per treatment x zone), black, bold, size 6,
    # centred on the cell and placed 15% of the panel's data range above that
    # cell's highest plotted point (author decision, 2026-08-09) — the previous
    # mean+SE anchor put letters inside the point cloud wherever a cell's
    # scatter extended well above its own error bar.
    cld_df <- .mg_read_cld_treatxzone(label_dir)
    if (remap_cld && !is.null(cld_df))
      cld_df <- .mg_remap_cld_zone(cld_df, zones_order, treatments_order)
    if (!is.null(cld_df)) {
      cld_join <- merge(cld_df,
                        .smry[, c("treatment","zone","mean_y","sem_y","max_y","cell_x")],
                        by = c("treatment","zone"), all.x = TRUE)
      cld_join <- cld_join[!is.na(cld_join$mean_y), , drop = FALSE]
      if (nrow(cld_join) > 0 && (always_show_cld || length(unique(cld_join$.group)) > 1L)) {
        cld_join$y <- cld_join$max_y + 0.15 * yrange
        p <- p + ggplot2::geom_text(
          data = cld_join,
          ggplot2::aes(x = cell_x, y = y, label = .group),
          inherit.aes = FALSE, fontface = "bold", size = 6,
          colour = "black", show.legend = FALSE)
        p <- p + ggplot2::coord_cartesian(
          ylim = c(ymin - 0.05 * yrange, ymax + 0.26 * yrange))
      }
    }
    p
  }, error = function(e) {
    warning("cells panel '", y_col, "' failed: ", e$message); NULL
  })
}

# ============ MG_BODY_H/LEG_H/SPC_H + MG_GHOST_THEME (6028-6054) ============
# .MG_SPC_H = 0.5 × .MG_LEG_H gives 0.5× caption-to-bottom padding for
# panels that carry no legend.
.MG_BODY_H <- 10
.MG_LEG_H  <- 1.5
.MG_SPC_H  <- 0.75

# Ghost theme: applied to a copy of a legend-carrying panel to produce a
# legend-only strip (all plot elements blanked, caption suppressed).
# Stacking body / ghost gives the order: plot → ANOVA caption → legend.
.MG_GHOST_THEME <- ggplot2::theme(
  panel.background = ggplot2::element_blank(),
  panel.border     = ggplot2::element_blank(),
  panel.grid.major = ggplot2::element_blank(),
  panel.grid.minor = ggplot2::element_blank(),
  axis.line        = ggplot2::element_blank(),
  axis.title       = ggplot2::element_blank(),
  axis.text        = ggplot2::element_blank(),
  axis.ticks       = ggplot2::element_blank(),
  plot.background  = ggplot2::element_blank(),
  plot.caption     = ggplot2::element_blank(),
  plot.margin      = ggplot2::margin(0, 0, 0, 0),
  legend.margin    = ggplot2::margin(2, 0, 0, 0)
)

# Caption theme override for Figure_10 (interval line plots): black text.
# Was 13 pt to match legend text; dropped to 11 on 2026-08-09 when partial
# eta-squared was appended to these captions (longest statement grew 43 -> 56

# ============ .FIG10_CAP_OVR (6055) ============
# characters, which at 13 pt would have run to the panel edge in this 2-column
# layout). Measured against the regenerated caption table, not guessed.
.FIG10_CAP_OVR <- ggplot2::theme(
  plot.caption = ggtext::element_markdown(
    size   = 11, hjust  = 0, face   = "italic",
    colour = "black", margin = ggplot2::margin(t = 10)
  )
)

# ---------------------------------------------------------------------------
# Parse p-value from an ANOVA statement string (e.g. "Treatment F₁,₁₁ = 10.348, p = 0.008"
# or "Treatment χ²₁ = 19.004, p < 0.001") and return TRUE if p < 0.05.
# Returns FALSE if the string is NULL, NA, or the p-value cannot be parsed.
# Convert Unicode subscript digits (U+2080–U+2089) to <sub>...</sub> HTML so
# that ggtext::element_markdown renders them correctly in both PNG and PDF.
# Standard PDF fonts (Helvetica, Arial) lack these codepoints, causing blank
# squares; <sub> tags are handled by ggtext's own renderer and are device-agnostic.
.sub_to_html <- function(s) {
  if (is.null(s) || is.na(s) || !nzchar(s)) return(s)
  sub_chars <- "₀₁₂₃₄₅₆₇₈₉"
  # Match runs of subscript digits (with optional embedded commas for df1,df2).
  pat <- paste0("[", sub_chars, "][", sub_chars, ",]*[", sub_chars, "]|[", sub_chars, "]")
  stringr::str_replace_all(s, pat, function(m) {
    paste0("<sub>", chartr(sub_chars, "0123456789", m), "</sub>")
  })
}

# ============ .sub_to_html (6070) ============
# Standard PDF fonts (Helvetica, Arial) lack these codepoints, causing blank
# squares; <sub> tags are handled by ggtext's own renderer and are device-agnostic.
.sub_to_html <- function(s) {
  if (is.null(s) || is.na(s) || !nzchar(s)) return(s)
  sub_chars <- "₀₁₂₃₄₅₆₇₈₉"
  # Match runs of subscript digits (with optional embedded commas for df1,df2).
  pat <- paste0("[", sub_chars, "][", sub_chars, ",]*[", sub_chars, "]|[", sub_chars, "]")
  stringr::str_replace_all(s, pat, function(m) {
    paste0("<sub>", chartr(sub_chars, "0123456789", m), "</sub>")
  })
}

# ============ Line-with-Tukey panel helpers (moved up from FIGURE 9 section,
# 2026-08-09) ============
# .mg_make_line_with_tukey() and its dependents were originally defined only
# where Figure 9 (interval) first needed them, further down. Figure 7bis
# (below) now also needs it, moved earlier so it can double as Figure_8 panel
# F -- so the whole self-contained helper bundle (pure function definitions,
# no side effects) moved up here, ahead of every place that calls it. Note:
# .mg_wrap_caption is redefined identically here to the forward-declared
# no-op version above (both just return text unchanged) -- harmless.
# Helper: build raw cell means/SE table from a long df by treatment x timepoint.
.mg_cell_summary <- function(df_in, y_col) {
  if (is.null(df_in) || !y_col %in% names(df_in) ||
      !"timepoint_f" %in% names(df_in)) return(NULL)
  df_in %>%
    dplyr::filter(is.finite(.data[[y_col]]), !is.na(timepoint_f)) %>%
    dplyr::group_by(treatment, timepoint_f) %>%
    dplyr::summarise(mean_y = mean(.data[[y_col]], na.rm = TRUE),
                     sem_y  = .sem(.data[[y_col]]),
                     .groups = "drop") %>%
    dplyr::mutate(
      tp_num = as.integer(as.character(timepoint_f)),
      cell   = paste(as.character(treatment), as.character(timepoint_f))
    )
}

# Helper: parse Tukey contrast label "control timepoint1 - exercise choice timepoint3"
# into a row of (cell1_treatment, cell1_tp, cell2_treatment, cell2_tp, p).
.mg_parse_contrast <- function(contr_row) {
  txt <- as.character(contr_row$contrast)
  parts <- strsplit(txt, " - ", fixed = TRUE)[[1]]
  if (length(parts) != 2) return(NULL)
  pat_tp <- "^(.+?)\\s+timepoint_?f?(\\d)$"
  m1 <- regmatches(parts[1], regexec(pat_tp, parts[1]))[[1]]
  m2 <- regmatches(parts[2], regexec(pat_tp, parts[2]))[[1]]
  if (length(m1) != 3 || length(m2) != 3) return(NULL)
  data.frame(
    tx1 = m1[2], tp1 = as.integer(m1[3]),
    tx2 = m2[2], tp2 = as.integer(m2[3]),
    p   = contr_row$p.value,
    stringsAsFactors = FALSE
  )
}

.mg_sig_stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  ""
}

# Helper: wrap subtitle text to a per-panel character width.
.mg_wrap_caption <- function(text, panel_width_chars = 26) {
  # No longer wraps (2026-08-09): the ANOVA caption must render on one
  # line. panel_width_chars is kept as a no-op parameter so call sites
  # don't need editing; text is returned unchanged.
  text
}

# Helper: read the (treatment x timepoint) CLD table written by
# run_lmm_analysis. Trims whitespace, truncates letters to max 2 chars,
# coerces timepoint to integer.
.mg_read_cld_tp <- function(label_dir, cld_path = NULL) {
  # cld_path: explicit CSV path, for panels whose CLD is not produced by a
  # STEP5 model directory (commitment_index comes from the bout engine).
  if (!is.null(cld_path) && nzchar(cld_path) && file.exists(cld_path)) {
    use_fp <- cld_path
  } else {
  if (is.null(label_dir) || !nzchar(label_dir)) return(NULL)
  fp   <- file.path(STEP5_OUT, label_dir, "cld_treatmentxtimepoint.csv")
  fp_f <- file.path(STEP5_OUT, label_dir, "cld_treatmentxtimepoint_f.csv")
  # Two possible filenames because the CLD tag follows the interaction term name:
  # cld_treatmentxtimepoint.csv for a continuous timepoint,
  # cld_treatmentxtimepoint_f.csv for the 3-level factor.
  # Under the uniform 2-df rule (.ALLOW_CONTINUOUS_TP = FALSE in the STATS
  # script) every model now writes the _f variant, so that is the branch taken
  # in practice; the non-_f branch is retained only so archived output from
  # before the revision still renders.
  use_fp <- if (file.exists(fp)) fp else if (file.exists(fp_f)) fp_f else return(NULL)
  }
  d <- tryCatch(read.csv(use_fp, stringsAsFactors = FALSE),
                error = function(e) NULL)
  if (is.null(d) || !all(c("treatment",".group") %in% names(d))) return(NULL)
  # Support "timepoint" (continuous model) or "timepoint_f" (factor model) as tp source
  tp_col <- if ("timepoint" %in% names(d)) "timepoint" else
            if ("timepoint_f" %in% names(d)) "timepoint_f" else return(NULL)
  # NOT truncated (fixed 2026-08-16). This was substr(..., 1, 2), which silently
  # destroyed the third letter of any 3-group cell and thereby changed which
  # cells appeared to differ -- IID's true groups are a|b|c|abc|abc|abc, and the
  # cut turned the three "abc" cells into "ab", removing the letter they shared
  # with the "c" cell. See cld_letter_policy.R for the full worked example.
  d$.group    <- trimws(as.character(d$.group))
  d$treatment <- trimws(as.character(d$treatment))
  d$tp_num    <- suppressWarnings(as.integer(as.character(d[[tp_col]])))
  d <- d[!is.na(d$tp_num) & nzchar(d$.group), , drop = FALSE]
  if (nrow(d) == 0) return(NULL)
  d
}

# Build the per-cell label table for plotting. Places a label just above the
# +SE error bar tip (7% of the data range gap) per (tp, treatment); when both
# treatments at the same tp share the same letter string, emits a single
# merged label above the higher of the two SE tips. Each row carries x_adj = 0
# (callers may override specific cells via cld_nudges).
# collapse_tp: integer vector of tp values where a single label is shown above
# max(both tops)+off, using the first letter of the Control treatment's group.
.mg_cld_label_rows <- function(cld_df, smry, yrange, force_split_tp = integer(0),
                                collapse_tp = integer(0)) {
  # House rule (2026-08-10): letters sit centred on the cell's x position at the
  # upper limit of its error bar plus 10% of the plotted data range.
  off <- 0.10 * yrange
  out <- list()
  for (tp in sort(unique(cld_df$tp_num))) {
    rows <- cld_df[cld_df$tp_num == tp, , drop = FALSE]
    sm   <- smry[smry$tp_num == tp, , drop = FALSE]
    if (!nrow(rows) || !nrow(sm)) next
    sm$top <- sm$mean_y + sm$sem_y
    rows <- merge(rows, sm[, c("treatment","top")],
                  by = "treatment", all.x = TRUE)
    rows <- rows[!is.na(rows$top), , drop = FALSE]
    if (!nrow(rows)) next
    if (tp %in% collapse_tp) {
      # One label for the interval, carrying the CONTROL cell's letters, drawn
      # above the higher of the two error-bar tips. Used where the two
      # treatments' letter strings differ textually but overlap (e.g. IID
      # interval 1: control "abc" vs exercise choice "c" -- they share "c", so
      # they do not differ, and two stacked labels misrepresent that).
      prim <- rows[rows$treatment == "control", , drop = FALSE]
      if (!nrow(prim)) prim <- rows[1L, , drop = FALSE]
      out[[length(out) + 1L]] <- data.frame(
        tp_num    = tp,
        treatment = NA_character_,
        y         = max(rows$top, na.rm = TRUE) + off,
        label     = prim$.group[1],
        merged    = TRUE,
        x_adj     = 0,
        stringsAsFactors = FALSE
      )
    } else if (length(unique(rows$.group)) == 1L && !tp %in% force_split_tp) {
      out[[length(out) + 1L]] <- data.frame(
        tp_num    = tp,
        treatment = NA_character_,
        y         = max(rows$top, na.rm = TRUE) + off,
        label     = rows$.group[1],
        merged    = TRUE,
        x_adj     = 0,
        stringsAsFactors = FALSE
      )
    } else {
      for (i in seq_len(nrow(rows))) {
        out[[length(out) + 1L]] <- data.frame(
          tp_num    = tp,
          treatment = rows$treatment[i],
          y         = rows$top[i] + off,
          label     = rows$.group[i],
          merged    = FALSE,
          x_adj     = 0,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (!length(out)) return(NULL)
  do.call(rbind, out)
}

# Remap CLD letter assignments so the first group encountered under
# precedence order (control → EC; interval 1→2→3) is labelled "a",
# the second new letter "b", etc.  Letters within each .group string
# are re-sorted alphabetically after remapping.
# Delegates to the canonical rule in cld_letter_policy.R so the letters drawn
# here and the letters tabulated by export_tukey_interval_module.R cannot drift
# apart. Kept as a thin wrapper because ~10 call sites already use this name.
.mg_remap_cld <- function(cld_df) .cld_remap_control_first(cld_df)

.mg_make_line_with_tukey <- function(df_src, y_col, y_label, label_dir,
                                      anova_caption,
                                      panel_width_chars = 20,
                                      cld_nudges = NULL,
                                      # Default flipped to TRUE (2026-08-16):
                                      # the control/interval-1 == "a" rule is a
                                      # house convention that must hold in EVERY
                                      # panel, so opting in per call site was
                                      # the wrong default -- a panel that forgot
                                      # it silently carried inconsistent
                                      # letters. Remapping is idempotent, so
                                      # applying it to already-correct input is
                                      # a no-op.
                                      remap_cld = TRUE,
                                      # "control_first" (house rule 1) or
                                      # "lone_first" -- see cld_letter_policy.R.
                                      # Only use "lone_first" on a panel whose
                                      # control/interval-1 cell carries more
                                      # than one letter, so the display would
                                      # otherwise have no plain "a" anywhere.
                                      remap_mode = "control_first",
                                      collapse_tp = integer(0),
                                      cld_path = NULL) {
  smry <- .mg_cell_summary(df_src, y_col)
  if (is.null(smry) || nrow(smry) == 0) return(NULL)

  cld_df <- .mg_read_cld_tp(label_dir, cld_path = cld_path)
  if (remap_cld && !is.null(cld_df))
    cld_df <- if (identical(remap_mode, "lone_first"))
                .cld_remap_lone_first(cld_df) else .mg_remap_cld(cld_df)
  # Post-condition: warn loudly if the invariant does not hold.
  .cld_check_control_first(cld_df, label = y_label, mode = remap_mode)
  # HOUSE RULE 2: when every cell shares one group, nothing differs and the
  # letters are noise -- six identical "a"s make a reader hunt for a contrast
  # that is not there. Drop them.
  if (.cld_is_uniform(cld_df)) {
    ts_msg("    CLD suppressed (all cells share one group): ", y_label)
    cld_df <- NULL
  }

  ymin <- min(smry$mean_y - smry$sem_y, na.rm = TRUE)
  ymax <- max(smry$mean_y + smry$sem_y, na.rm = TRUE)
  yrange <- ymax - ymin
  if (!is.finite(yrange) || yrange == 0) yrange <- 1

  p <- ggplot2::ggplot(
      smry,
      ggplot2::aes(x = tp_num, y = mean_y,
                   colour   = treatment,
                   shape    = treatment,
                   linetype = treatment,
                   group    = treatment)
    ) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = PT_SIZE * 1.6) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = mean_y - sem_y, ymax = mean_y + sem_y),
      colour = "black", width = 0.10, linewidth = LW_ERR,
      show.legend = FALSE
    ) +
    ggplot2::scale_colour_manual(values = TREATMENT_COLORS_PAL, name = "Treatment") +
    ggplot2::scale_shape_manual(values  = TREATMENT_SHAPES,    name = "Treatment") +
    ggplot2::scale_linetype_manual(
      values = c("exercise choice" = "solid", "control" = "dashed"),
      name   = "Treatment"
    ) +
    ggplot2::scale_x_continuous(
      breaks = 1:3,
      labels = c("5–25", "45–65", "85–105"),
      limits = c(0.6, 3.4)
    ) +
    BASE_THEME +
    ggplot2::labs(x = "Interval (min)", y = y_label,
                  caption = .mg_wrap_caption(anova_caption,
                                             panel_width_chars)) +
    ggplot2::theme(
      plot.caption = ggtext::element_markdown(size = 11, hjust = 0,
                                              face = "italic", colour = "grey25",
                                              margin = ggplot2::margin(t = 10)),
      axis.title.x = ggplot2::element_text(size = 13, face = "bold", color = "black")
    )

  # CLD letters per (treatment x timepoint) cell, all black, size 6, positioned
  # just above the +SE error bar tip (7% of yrange gap). Shared-letter cells at
  # the same interval render as a single label above the higher SE tip.
  # Per-cell x nudges are applied via cld_nudges (data frame: treatment, tp_num, x_adj,
  # optionally y_extra_frac, y_at_mean, y_at_lower_se, suppress, force_separate).
  if (!is.null(cld_df)) {
    # force_separate: tp_num values where merging is prevented even when letters match
    force_split_tp <- if (!is.null(cld_nudges) && "force_separate" %in% names(cld_nudges))
      cld_nudges$tp_num[!is.na(cld_nudges$force_separate) & cld_nudges$force_separate]
    else integer(0)
    lab <- .mg_cld_label_rows(cld_df, smry, yrange, force_split_tp = force_split_tp,
                               collapse_tp = collapse_tp)
    if (!is.null(lab) && nrow(lab) > 0) {
      if (!is.null(cld_nudges) && nrow(cld_nudges) > 0 &&
          all(c("treatment","tp_num") %in% names(cld_nudges))) {
        lab <- merge(lab, cld_nudges, by = c("treatment","tp_num"), all.x = TRUE,
                     suffixes = c("", ".nudge"))
        # x_adj: replace default 0 with nudge value where provided
        if ("x_adj.nudge" %in% names(lab)) {
          lab$x_adj <- ifelse(is.na(lab$x_adj.nudge), lab$x_adj, lab$x_adj.nudge)
          lab$x_adj.nudge <- NULL
        }
        # y_extra_frac: additional vertical offset (fraction of yrange) for nudged letters
        if ("y_extra_frac" %in% names(lab)) {
          lab$y <- lab$y + ifelse(is.na(lab$y_extra_frac), 0, lab$y_extra_frac * yrange)
          lab$y_extra_frac <- NULL
        }
        # y_at_mean: override y to the treatment x timepoint cell mean (not above SE bar)
        if ("y_at_mean" %in% names(lab)) {
          need_y <- !is.na(lab$y_at_mean) & lab$y_at_mean & !is.na(lab$treatment)
          if (any(need_y)) {
            mean_lkp <- stats::setNames(smry$mean_y,
                                        paste(smry$treatment, smry$tp_num))
            for (i in which(need_y)) {
              mv <- mean_lkp[paste(lab$treatment[i], lab$tp_num[i])]
              if (!is.na(mv) && is.finite(mv)) lab$y[i] <- mv
            }
          }
          lab$y_at_mean <- NULL
        }
        # y_at_lower_se: put the letter level with the BOTTOM of its own error
        # bar (mean - SE) instead of above the top of it. For a cell whose
        # marker sits low and close to the other series, the default position
        # lands the letter in the middle of the other group's error bar.
        if ("y_at_lower_se" %in% names(lab)) {
          need_lo <- !is.na(lab$y_at_lower_se) & lab$y_at_lower_se &
                     !is.na(lab$treatment)
          if (any(need_lo)) {
            lo_lkp <- stats::setNames(smry$mean_y - smry$sem_y,
                                      paste(smry$treatment, smry$tp_num))
            for (i in which(need_lo)) {
              lv <- lo_lkp[paste(lab$treatment[i], lab$tp_num[i])]
              if (!is.na(lv) && is.finite(lv)) lab$y[i] <- lv
            }
          }
          lab$y_at_lower_se <- NULL
        }
        # y_at_treat_avg: set y to grand mean of ALL treatment cell means at that
        # tp_num (from smry) plus the standard 0.11*yrange offset.  Used for a
        # single shared label positioned vertically at the midpoint between groups.
        if ("y_at_treat_avg" %in% names(lab)) {
          need_avg <- !is.na(lab$y_at_treat_avg) & lab$y_at_treat_avg
          if (any(need_avg)) {
            for (i in which(need_avg)) {
              tp_means <- smry$mean_y[smry$tp_num == lab$tp_num[i]]
              tp_means <- tp_means[is.finite(tp_means)]
              if (length(tp_means) >= 1)
                lab$y[i] <- mean(tp_means) + 0.11 * yrange
            }
          }
          lab$y_at_treat_avg <- NULL
        }
        # y_centered_gap: centre the letter vertically in the empty band between
        # the two treatments' error bars at that interval (own SE tip facing the
        # other group, to the other group's facing SE tip). Used where the
        # default "above my own error bar" position would collide with the other
        # treatment's marker.
        if ("y_centered_gap" %in% names(lab)) {
          need_gap <- !is.na(lab$y_centered_gap) & lab$y_centered_gap &
                      !is.na(lab$treatment)
          if (any(need_gap)) {
            for (i in which(need_gap)) {
              sm_tp <- smry[smry$tp_num == lab$tp_num[i], , drop = FALSE]
              own   <- sm_tp[sm_tp$treatment == lab$treatment[i], , drop = FALSE]
              oth   <- sm_tp[sm_tp$treatment != lab$treatment[i], , drop = FALSE]
              if (nrow(own) == 1L && nrow(oth) == 1L) {
                if (oth$mean_y > own$mean_y) {
                  lab$y[i] <- mean(c(own$mean_y + own$sem_y,
                                     oth$mean_y - oth$sem_y))
                } else {
                  lab$y[i] <- mean(c(own$mean_y - own$sem_y,
                                     oth$mean_y + oth$sem_y))
                }
              }
            }
          }
          lab$y_centered_gap <- NULL
        }
        # label_override: replace computed label text for specific cells
        if ("label_override" %in% names(lab)) {
          has_ovr <- !is.na(lab$label_override) & nzchar(lab$label_override)
          if (any(has_ovr)) lab$label[has_ovr] <- lab$label_override[has_ovr]
          lab$label_override <- NULL
        }
        # suppress: drop specific cells from the label set
        if ("suppress" %in% names(lab)) {
          lab <- lab[is.na(lab$suppress) | !lab$suppress, , drop = FALSE]
          if ("suppress" %in% names(lab)) lab$suppress <- NULL
        }
        # force_separate: already consumed above; remove from lab if present
        if ("force_separate" %in% names(lab)) lab$force_separate <- NULL
      }
      lab_split <- !lab$merged
      if (any(lab_split)) {
        p <- p + ggplot2::geom_text(
          data = lab[lab_split, , drop = FALSE],
          ggplot2::aes(x = tp_num + x_adj, y = y, label = label),
          inherit.aes = FALSE, fontface = "bold", size = 6,
          colour = "black", show.legend = FALSE
        )
      }
      if (any(lab$merged)) {
        p <- p + ggplot2::geom_text(
          data = lab[lab$merged, , drop = FALSE],
          ggplot2::aes(x = tp_num + x_adj, y = y, label = label),
          inherit.aes = FALSE, fontface = "bold", size = 6,
          colour = "black"
        )
      }
      # Expand top margin: base offset 11% + max nudge extra 4% + buffer → 24%.
      p <- p + ggplot2::coord_cartesian(
        ylim = c(ymin - 0.05 * yrange, ymax + 0.24 * yrange))
    }
  }
  p
}

# Build a short ANOVA subtitle for each panel: Treatment × Interval interaction
# term only, with full-word labels to avoid abbreviation in the figure.
# Wrapping is applied at render time via .mg_wrap_caption(panel_width_chars).
.mg_anova_cap <- function(res_tp_obj, label_main_term) {
  if (is.null(res_tp_obj) || is.null(res_tp_obj$anova)) return(NULL)
  av <- res_tp_obj$anova
  row <- av[grepl("treatment.*timepoint|timepoint.*treatment", av$term,
                   ignore.case = TRUE, perl = TRUE), ]
  if (nrow(row) == 0) return(NULL)
  p <- suppressWarnings(as.numeric(row$p_value[1]))
  if (!is.finite(p)) return(NULL)
  p_str <- fmt_p(p)
  df1 <- if ("df" %in% names(row)) as.integer(round(row$df[1])) else NA_integer_
  df2_raw <- if ("df_denom" %in% names(row)) as.numeric(row$df_denom[1]) else NA_real_
  chi <- fmt_F(as.numeric(row$chisq[1]))
  if (is.finite(df2_raw)) {
    eta <- fmt3((as.numeric(row$chisq[1]) * df1) / (as.numeric(row$chisq[1]) * df1 + df2_raw))
    sprintf("Treatment × Interval: F<sub>%d, %s</sub> = %s, %s, &eta;<sup>2</sup><sub>p</sub> = %s",
            df1, fmt_F(df2_raw), chi, p_str, eta)
  } else {
    sprintf("Treatment × Interval: F<sub>%d</sub> = %s, %s", df1, chi, p_str)
  }
}

# ---------------------------------------------------------------------------
# Figure_10 ANOVA captions — sourced from anova_statement_table_interval.csv.
# The helper returns the "statement" string for a given indicator label,
# falling back to NULL (which suppresses the caption) if the file or row
# is missing.
.fig10_cap_csv <- local({
  .csv_path <- file.path(PROJECT_ROOT, "choice R pipeline/easy_scripts/standalone/anova_statement_table_interval.csv")
  .tbl <- tryCatch(read.csv(.csv_path, stringsAsFactors = FALSE, encoding = "UTF-8"),
                   error = function(e) NULL)
  function(indicator) {
    if (is.null(.tbl)) return(NULL)
    row <- .tbl[trimws(.tbl$indicator) == indicator, , drop = FALSE]
    if (nrow(row) == 0) return(NULL)
    .sub_to_html(trimws(row$statement[1]))
  }
})
# ---------------------------------------------------------------------------

# ============ FIGURE 7 section (csv_is_sig through Figure_8 save) (6083-6308) ============
.csv_is_sig <- function(statement) {
  if (is.null(statement) || is.na(statement) || !nzchar(statement)) return(FALSE)
  # "p < 0.001" → always significant
  if (grepl("p\\s*<\\s*0\\.0", statement)) return(TRUE)
  # "p = 0.NNN" → extract numeric value
  m <- regmatches(statement, regexpr("p\\s*=\\s*([0-9]+\\.[0-9]+)", statement))
  if (length(m) == 0 || !nzchar(m)) return(FALSE)
  pval <- suppressWarnings(as.numeric(sub("p\\s*=\\s*", "", m)))
  is.finite(pval) && pval < 0.05
}
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Trial-level ANOVA captions — sourced from anova_statement_table_trial.csv.
# Lookup by variable name (column "variable"); returns the "statement" string,
# or NULL if the file or row is missing. Used for Figure_8 (zone/ALR) and
# Figure_9 (collective trial-level) panels.
.fig9_cap_csv <- local({
  .csv_path <- file.path(PROJECT_ROOT, "choice R pipeline/easy_scripts/standalone/anova_statement_table_trial.csv")
  .tbl <- tryCatch(read.csv(.csv_path, stringsAsFactors = FALSE, encoding = "UTF-8"),
                   error = function(e) NULL)
  function(variable) {
    if (is.null(.tbl)) return(NULL)
    row <- .tbl[trimws(.tbl$variable) == variable, , drop = FALSE]
    if (nrow(row) == 0) return(NULL)
    .sub_to_html(trimws(row$statement[1]))
  }
})
# ---------------------------------------------------------------------------

# ---- Figure 7: 2 x 2 — alr panels (top) + cell-mean CLD panels (bottom) ----
.fig7_A <- .mg_frame(.mg_make_alr_scatter(
  df_in        = if (!is.null(df_main_wide_agg) && "logit_flow" %in% names(df_main_wide_agg))
                   df_main_wide_agg else NULL,
  y_col        = "logit_flow",
  y_label      = "ALR(flow vs calm)",
  res_obj      = if (exists("res_zone_flow_logit_agg")) res_zone_flow_logit_agg else NULL,
  n_cap_df     = df_main_wide_agg,
  show_n       = FALSE,
  cap_csv      = .fig9_cap_csv("logit_flow"),
  sig_from_csv = .csv_is_sig(.fig9_cap_csv("logit_flow"))))

.fig7_B <- .mg_frame(.mg_make_alr_scatter(
  df_in        = if (!is.null(df_sec_wide_agg) && "lr_high" %in% names(df_sec_wide_agg))
                   df_sec_wide_agg else NULL,
  y_col        = "lr_high",
  y_label      = "ALR(high vs calm)",
  res_obj      = if (exists("res_zone_sec_high_agg")) res_zone_sec_high_agg else NULL,
  n_cap_df     = df_sec_wide_agg,
  show_n       = FALSE,
  cap_csv      = .fig9_cap_csv("lr_high"),
  sig_from_csv = .csv_is_sig(.fig9_cap_csv("lr_high"))))

# Trial-level cell-mean ANOVA captions -- sourced from
# anova_statement_table_cells.csv (rebuild_figure_caption_table_cells.R).
# Added 2026-08-09 alongside figures_step5_module.R: this is the one panel
# family whose caption previously came ONLY from a live res_obj$anova table
# (.fig7_anova_cap()), with no CSV cache like every sibling panel -- the last
# obstacle to a model-fit-free figure-only regeneration path. res_obj is kept
# as a fallback so the live pipeline is unaffected if the cache is stale/missing.
.fig7_cells_cap_csv <- local({
  .csv_path <- file.path(PROJECT_ROOT, "choice R pipeline/easy_scripts/standalone/anova_statement_table_cells.csv")
  .tbl <- tryCatch(read.csv(.csv_path, stringsAsFactors = FALSE, encoding = "UTF-8"),
                   error = function(e) NULL)
  function(variable) {
    if (is.null(.tbl)) return(NULL)
    row <- .tbl[trimws(.tbl$variable) == variable, , drop = FALSE]
    if (nrow(row) == 0 || is.na(row$statement[1])) return(NULL)
    trimws(row$statement[1])
  }
})

.fig7_C <- .mg_frame(.mg_make_cells_panel(
  df_long               = df_main_cells_agg,
  y_col                 = "p",
  y_label               = "Zone occupancy",
  zones_order           = c("flow","calm"),
  label_dir             = "zone_main_cells_aggregated",
  res_obj               = res_zone_main_cells_agg,
  cap_csv               = .fig7_cells_cap_csv("zone_main"),
  hline_at              = NULL,
  panel_width_chars     = 32,
  remap_cld             = TRUE,
  show_zone_legend      = FALSE,
  show_treatment_legend = TRUE))

.fig7_D <- .mg_frame(.mg_make_cells_panel(
  df_long               = df_sec_cells_agg,
  y_col                 = "clr",
  y_label               = "area normalized CLR",
  zones_order           = c("high","medium","low","calm"),
  label_dir             = "zone_sec_cells_aggregated",
  res_obj               = res_zone_sec_cells_agg,
  cap_csv               = .fig7_cells_cap_csv("zone_sec"),
  hline_at              = NULL,
  panel_width_chars     = 32,
  remap_cld             = TRUE,
  show_zone_legend      = FALSE,
  show_treatment_legend = FALSE))

# ---- Figure 7bis panel: alr(flow) by treatment x interval (single panel) ----
# Built here (moved up from its original standalone-only position, 2026-08-09)
# so the same native ggplot object can double as Figure_8 panel F below, in
# addition to still being saved on its own further down. Already
# .mg_frame()-wrapped inside .mg_make_line_with_tukey()'s caller just like
# .fig7_A-.fig7_D, so it needs no raster embedding and no scale compensation
# to match its neighbours -- same theme, same wrapper, same native panel path.
.fig7bis_cap <- .fig7_cells_cap_csv("zone_main_timepoint") %||%
                (if (exists("res_zone_main")) .mg_anova_cap(res_zone_main) else NULL)
.fig7bis <- if (!is.null(df_main_wide) && "logit_flow" %in% names(df_main_wide)) {
  .mg_frame(.mg_make_line_with_tukey(
    df_src        = df_main_wide,
    y_col         = "logit_flow",
    y_label       = "ALR(flow vs calm)",
    label_dir     = "zone_main_timepoint",
    anova_caption = .fig7bis_cap,
    remap_cld     = TRUE,
    # Control at interval 2 sits well below exercise choice; its default label
    # position (own SE tip + 11% of range) landed on the exercise-choice marker.
    # Centre it in the gap between the two groups' error bars instead.
    cld_nudges    = data.frame(
      treatment      = "control",
      tp_num         = 2L,
      y_centered_gap = TRUE,
      stringsAsFactors = FALSE)
  ))
} else NULL

# ---- Legend images for Figures 8–10 (external PNG files) ----------------
# Loaded once here; re-used by Figure 9 and Figure 10 assembly below.
# Hard-stop if either file is missing — do not silently substitute.
{
  .leg_dir  <- file.path(PROJECT_ROOT, "manu_graphs_26.05.2026/Legends")
  .ab_path  <- file.path(.leg_dir, "Figure_10_legend_panelAB.PNG")
  .trt_path <- file.path(.leg_dir, "treatment_legend.PNG")
  for (.p in c(.ab_path, .trt_path)) {
    if (!file.exists(.p))
      stop("[Figure assembly] Legend file not found: ", .p,
           "\n  Place the file there and re-run.")
  }
  # Wrap a PNG as a ggplot panel that preserves the original pixel aspect ratio.
  # coord_equal() enforces asp = px_height/px_width; any extra cell space is
  # filled by the white plot background — no stretching of legend text/symbols.
  .gg_legend <- function(img) {
    asp <- dim(img)[1] / dim(img)[2]        # px_height / px_width
    ggplot2::ggplot() +
      ggplot2::annotation_custom(
        grid::rasterGrob(img, interpolate = TRUE),
        xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
      ) +
      ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
      ggplot2::scale_y_continuous(limits = c(0, asp), expand = c(0, 0)) +
      ggplot2::coord_equal() +
      ggplot2::theme_void() +
      ggplot2::theme(
        plot.background  = ggplot2::element_rect(fill = "white", colour = NA),
        plot.margin      = ggplot2::margin(0, 0, 0, 0)
      )
  }
  .legend_ab_elem  <- .gg_legend(png::readPNG(.ab_path))
  .legend_trt_elem <- .gg_legend(png::readPNG(.trt_path))
  ts_msg("[legends] loaded: Figure_10_legend_panelAB.PNG + treatment_legend.PNG")
}

# ---- Panels E, F: sequence-module metrics (raster-embedded, 2026-08-09) ----
# Panel E: p_stay_Flow (raster-embedded, sequence-module engine). Re-added
# 2026-08-09 after an earlier removal; the size mismatch against the native
# A-D panels traced to the destination-cell aspect ratio (see
# fig_seq_manuscript_embeds() in figures_seq_bout_module.R) is compensated via
# its scale = 1.25 render parameter, the same fix validated on Figure_9's A/F.
.seq_fig_dir     <- file.path(PROJECT_ROOT, "all_manu_graphs/sequence_and_patterns_graphs")
.pstay_flow_path <- file.path(.seq_fig_dir, "seq_metric_binary_p_stay_Flow.png")
if (!file.exists(.pstay_flow_path))
  stop("[Figure assembly] Sequence-module figure not found: ", .pstay_flow_path,
       "\n  Run figures_seq_bout_module.R first.")
.fig7_E <- .mg_frame(.gg_legend(png::readPNG(.pstay_flow_path)))

# Panel F: Figure 7bis (alr(flow) by treatment x interval), reusing the native
# ggplot object built above -- same theme/wrapper as A-D, so no raster
# embedding or scale compensation is needed for this panel.
.fig7_F <- .fig7bis

if (!any(sapply(list(.fig7_A, .fig7_B, .fig7_C, .fig7_D, .fig7_E, .fig7_F), is.null))) {
  # Legend-free panels + panel tags
  .f7_body_A <- .fig7_A + ggplot2::theme(legend.position = "none") + .MG_TAG_THEME + ggplot2::labs(tag = "A")
  .f7_body_B <- .fig7_B + ggplot2::theme(legend.position = "none") + .MG_TAG_THEME + ggplot2::labs(tag = "B")
  .f7_body_C <- .fig7_C + ggplot2::theme(legend.position = "none") + .MG_TAG_THEME + ggplot2::labs(tag = "C")
  .f7_body_D <- .fig7_D + ggplot2::theme(legend.position = "none") + .MG_TAG_THEME + ggplot2::labs(tag = "D")
  .f7_body_E <- .fig7_E + .MG_TAG_THEME + ggplot2::labs(tag = "E")
  .f7_body_F <- .fig7_F + ggplot2::theme(legend.position = "none") + .MG_TAG_THEME + ggplot2::labs(tag = "F")
  # Single shared legend centred below the lowest panel row (E/F), matching
  # Figure 9's disposition. Shape encodes treatment in A/B/E/F; colour maps to
  # treatment in A/B/E/F but to zone in C/D, so only the shape (treatment)
  # legend is shared across the whole grid.
  .f7_top <- (patchwork::wrap_elements(full = .f7_body_A) |
               patchwork::wrap_elements(full = .f7_body_B))
  .f7_mid <- (patchwork::wrap_elements(full = .f7_body_C) |
               patchwork::wrap_elements(full = .f7_body_D))
  .f7_bot <- (patchwork::wrap_elements(full = .f7_body_E) |
               patchwork::wrap_elements(full = .f7_body_F))
  .fig7 <- (patchwork::wrap_elements(full = .f7_top) /
             patchwork::wrap_elements(full = .f7_mid) /
             patchwork::wrap_elements(full = .f7_bot) /
             .legend_trt_elem) +
             # Legend row enlarged 0.5 -> 1.8 (2026-08-10). The legend is an
             # aspect-locked raster (1369x167 px), so at 0.5 it was height-
             # limited to ~9 mm in a 320 mm figure and barely legible; 1.8
             # gives ~29 mm, close to the width-limited maximum.
             patchwork::plot_layout(heights = c(6, 6, 6, 1.8)) +
             patchwork::plot_annotation(
               theme = ggplot2::theme(
                 plot.background = ggplot2::element_rect(fill = "white", colour = NA),
                 plot.margin     = ggplot2::margin(5, 5, 5, 5, unit = "mm")
               )
             )
  .mg_save(.fig7, "Figure_8", width_mm = 280, height_mm = 320)
} else {
  ts_msg("  Figure 8 skipped — missing component plot(s)")
}
writeLines(c(
  "Figure 7. Zone preference by treatment (N = 8 trials per treatment).",
  "",
  "(A) Additive log-ratio of flow vs calm occupancy (alr(flow)) per trial.",
  "    Points jittered around treatment group. Crossbar = mean; error bar = +/- 1 SE.",
  "    Asterisk (*) indicates significant treatment effect (LMM, Kenward-Roger).",
  "",
  "(B) Additive log-ratio of high sub-zone vs calm occupancy (alr(high)) per trial.",
  "    As panel A.",
  "",
  "(C) Cell-mean proportion of observation time in each main zone (flow, calm)",
  "    by treatment. Cells from a beta generalised linear mixed model with logit",
  "    link (p ~ treatment * zone + (1|trial)). Horizontal dashed line at 0.5",
  "    indicates equal flow:calm occupancy. CLD letters denote Tukey-adjusted",
  "    pairwise cell groupings (all cells group a: no pairwise differences).",
  "",
  "(D) Cell-mean centred log-ratio (CLR) of the area-normalised 4-part composition",
  "    (high, medium, low, calm) by treatment. Cells from a Gaussian linear mixed",
  "    model (clr ~ treatment * zone + (1|trial) + (1|zone:trial)). Horizontal",
  "    dashed line at 0 indicates use proportional to the geometric mean of the",
  "    composition. CLD letters denote Tukey-adjusted pairwise cell groupings.",
  "",
  "(E) Probability the school remains in the flow zone from one 1-s bin to the",
  "    next (binary alphabet), per trial. Gaussian LMM, treatment fixed effect.",
  "    Sourced from the sequence (bin-state) analysis, not refitted here.",
  "    DESCRIPTIVE: shares 76% of its school-level variance with ALR(flow)",
  "    (r = 0.87; Table S19) and is reported as a re-expression of that result,",
  "    not independent confirmation.",
  "",
  "(F) Additive log-ratio of flow vs calm occupancy (alr(flow)) by treatment",
  "    across observation intervals -- the interval-resolved view of panel A",
  "    (N = 8 trials per treatment per interval; intervals: 5-25 min, 45-65 min,",
  "    85-105 min). Lines connect treatment group means; error bars = +/- 1 SE;",
  "    points are trial-level values. CLD letters denote Tukey-adjusted pairwise",
  "    groupings across treatment x interval cells.",
  "",
  "Bars show cell mean +/- 1 SE; points are trial-level values (jittered).",
  "ANOVA captions emitted per panel only when the treatment, treatment x zone,",
  "or treatment x interval term is significant (p < 0.05; 3 d.p.). Cell-mean",
  "analyses (C, D) are exploratory — Tukey-adjusted within-model only; no BH",
  "correction applied."
), file.path(.mg_dir, "figure7_caption.txt"))
ts_msg("  Figure 7 caption written.")


# ============ FIGURE 8 section -> saves Figure_9 (6353-6578) ============
# ---- Figure 8 — Trial-level collective + transitions (6 panels A-F, 3x2) ----
# Panel A: flow-calm transitions (moved here from standalone Figure 9)
# Panels B-F: NND, IID, Convex hull, Centroid speed
# CLD letters instead of asterisks; ANOVA statement + N label per panel.

# ANOVA caption formatter for Figure 8: emits F / χ² stat only when p < 0.05;
# p formatted to exactly 3 decimal places (e.g., p = 0.019, p < 0.001).
.fig8_anova_cap <- function(res) {
  if (is.null(res) || is.null(res$anova)) return(NULL)
  av  <- res$anova
  row <- av[grepl("^treatment$", av$term, ignore.case = TRUE, perl = TRUE), ]
  if (nrow(row) == 0L) return(NULL)
  row <- row[1L, ]
  p   <- suppressWarnings(as.numeric(row$p_value))
  if (!is.finite(p)) return(NULL)
  p_str  <- fmt_p(p)
  chi    <- fmt_F(as.numeric(row$chisq))
  st     <- attr(av, "stat_type") %||% "Wald-chisq"
  is_f   <- isTRUE(grepl("^F", st))
  if (is_f) {
    df1 <- as.integer(row$df)
    df2 <- suppressWarnings(as.numeric(
             if (!is.null(row$df_denom)) row$df_denom else NA_real_))
    if (is.finite(df2)) {
      eta <- fmt_es3((as.numeric(row$chisq) * df1) / (as.numeric(row$chisq) * df1 + df2))
      return(sprintf("F<sub>%d, %s</sub> = %s, %s, &eta;<sup>2</sup><sub>p</sub> = %s",
                     df1, fmt_F(df2), chi, p_str, eta))
    } else
      return(sprintf("F<sub>%d</sub> = %s, %s", df1, chi, p_str))
  } else {
    return(sprintf("χ²<sub>%d</sub> = %s, %s", as.integer(row$df), chi, p_str))
  }
}

# Figure 9 only: split the one-line ANOVA statement so the F/df/p part fills
# the first line and the effect size drops to a second. The six panels of this
# figure are the narrowest in the manuscript; on one line the caption had to
# run at 7.5 pt to fit, which is below the journal's minimum readable size.
# Every other figure keeps the single-line form (author decision, 2026-08-09).
# Two spellings of the effect size reach here: the caption CSV carries a literal
# U+03B7 ("η²ₚ"), the in-script fallback .fig8_anova_cap() emits the HTML entity
# "&eta;". Match either so the break lands regardless of the caption's source.
.mg_break_before_eta <- function(txt) {
  if (is.null(txt) || is.na(txt) || !nzchar(txt)) return(txt)
  sub(",\\s*(&eta;|η)", ",<br>\\1", txt)
}

.mg_mk_collective_panel <- function(res, y_col, y_label, dat = NULL,
                                     y_limits = NULL,
                                     show_colour_legend = TRUE,
                                     show_n = TRUE,
                                     cap_csv = NULL,
                                     sig_from_csv = NULL) {
  .dat <- if (!is.null(dat)) dat else {
    if (!is.null(df_gd_agg) && y_col %in% names(df_gd_agg)) df_gd_agg else NULL
  }
  if (is.null(.dat)) return(NULL)
  tryCatch({
    .lvls <- levels(factor(.dat$treatment))
    .tmap <- stats::setNames(seq_along(.lvls), .lvls)
    # Determine significance for bracket: use CSV-parsed p-value when provided,
    # falling back to the LMM result object. This ensures the bracket matches
    # the caption statement from anova_statement_table_trial.csv.
    .is_sig_val <- if (!is.null(sig_from_csv)) sig_from_csv else .is_sig(res, "Treatment")
    p <- .make_std_plot(
      df          = .dat,
      y_col       = y_col,
      y_label     = y_label,
      title       = "",
      caption_txt = .mg_break_before_eta(
                      if (!is.null(cap_csv)) cap_csv else .fig8_anova_cap(res)),
      cld_df      = NULL,
      is_sig      = .is_sig_val,
      facet_tp    = FALSE,
      n_txt       = NULL,
      sig_textsize = 8
    )
    if (show_n) p <- .mg_add_n_labels(p, .dat, y_col, .tmap)
    # When show_colour_legend = FALSE, suppress colour from the legend so only
    # shape remains (colours are still rendered on the points themselves).
    if (!show_colour_legend)
      p <- p + ggplot2::scale_colour_manual(
                 values = TREATMENT_COLORS_PAL, name = "Treatment", guide = "none")
    if (!is.null(y_limits))
      p <- p + ggplot2::coord_cartesian(ylim = y_limits)
    .mg_frame(p)
  }, error = function(e) {
    warning("Figure 8 panel '", y_col, "' failed: ", e$message); NULL
  })
}

# PANEL A REPLACED 2026-08-09: switches_per_session (the STEP5 collective-
# metrics engine's own raw flow<->calm transition count, F(1,8) = 0.105,
# p = 0.754, ns) is swapped for switch_rate (the sequence engine's per-minute
# transition RATE from the 1-s bin state series, F(1,10) = 59.67, p < 0.001,
# eta2p = 0.856 -- see seq_anova_trial_level.csv). Same underlying behaviour
# (flow<->calm switching), a materially stronger and better-normalised
# measurement of it. Raster-embedded from figures_seq_bout_module.R's
# fig_seq_manuscript_embeds() for the same reason as Figure 7 panels E/F: the
# sequence engine's own model is the source of truth, never re-fitted here.
.seq_fig_dir2      <- file.path(PROJECT_ROOT, "all_manu_graphs/sequence_and_patterns_graphs")
.switch_rate_path  <- file.path(.seq_fig_dir2, "seq_switch_rate_manuscript.png")
.entropy_gr_path   <- file.path(.seq_fig_dir2, "seq_metric_graded_entropy_rate.png")
for (.p in c(.switch_rate_path, .entropy_gr_path)) {
  if (!file.exists(.p))
    stop("[Figure assembly] Sequence-module figure not found: ", .p,
         "\n  Run figures_seq_bout_module.R first.")
}
.fig8_A <- .mg_frame(.gg_legend(png::readPNG(.switch_rate_path)))
.fig8_F <- .mg_frame(.gg_legend(png::readPNG(.entropy_gr_path)))
.fig8_B <- .mg_treatment_xaxis(.mg_mk_collective_panel(res_nnd_agg,
  "mean_nnd_cm",          "Mean nearest-neighbour distance (cm)",                                  show_n = FALSE,
  cap_csv = .fig9_cap_csv("mean_nnd_cm"),
  sig_from_csv = .csv_is_sig(.fig9_cap_csv("mean_nnd_cm"))))
# Every indicator in this figure (transitions, NND, IID, school area, school
# speed) is a permutation-invariant function of the per-frame position set and
# is therefore immune to identity error. Alignment/heading-based measures are
# not, and are not computed or plotted anywhere in this module -- see
# METHODS_CHANGES.md section 7.
# Panels renumbered 2026-08-09: old D/E/F -> new C/D/E, giving a 5-panel figure.
.fig8_C <- .mg_treatment_xaxis(.mg_mk_collective_panel(res_iid_agg,
  "mean_iid_cm",          "Mean inter-individual distance (cm)",          show_colour_legend = FALSE, show_n = FALSE,
  cap_csv = .fig9_cap_csv("mean_iid_cm"),
  sig_from_csv = .csv_is_sig(.fig9_cap_csv("mean_iid_cm"))))
.fig8_D <- .mg_treatment_xaxis(.mg_mk_collective_panel(res_hull_agg,
  "mean_hull_area_cm2",   "School area (cm²)",       show_colour_legend = FALSE, show_n = FALSE,
  cap_csv = .fig9_cap_csv("mean_school_area_cm2"),
  sig_from_csv = .csv_is_sig(.fig9_cap_csv("mean_school_area_cm2"))))
.fig8_E <- .mg_treatment_xaxis(.mg_mk_collective_panel(res_cspd_agg,
  "mean_centroid_spd_cm", "School speed (cm/s)",
  show_colour_legend = FALSE, show_n = FALSE,
  cap_csv = .fig9_cap_csv("mean_school_speed_cm_s"),
  sig_from_csv = .csv_is_sig(.fig9_cap_csv("mean_school_speed_cm_s"))))

if (!any(sapply(list(.fig8_A, .fig8_B, .fig8_C, .fig8_D, .fig8_E, .fig8_F), is.null))) {
  # Panel tag theme for Figure 9: tags placed in the left margin (further left
  # than "topleft" which anchors to the panel boundary).
  .f8_tag <- .MG_TAG_THEME + ggplot2::theme(plot.tag.position = c(0, 1.0))
  # Caption size reduced from BASE_THEME (14 pt) so all 5 captions fit within
  # each narrow panel without truncation. Raised 7.5 -> 9 on 2026-08-09 once
  # .mg_break_before_eta() moved the effect size onto a second line: line 1 is
  # then at most 35 characters ("Treatment F1,8 = 0.1048, p = 0.754,") vs the
  # 47-character single-line string that was CONFIRMED to fit at 7.5 pt in
  # this exact layout, so 7.5 * (47/35) * 0.9 safety margin = 9.06 -> 9 pt.
  # (11 pt was tried first and clipped panel A/B on the rendered PNG -- do
  # not re-attempt a value above ~9.5 without re-verifying on real output;
  # a full STEP5 re-run costs ~50 min, so this is deliberately conservative.)
  .f8_cap_theme <- ggplot2::theme(
    plot.caption = ggtext::element_markdown(size = 9, hjust = 0, face = "italic",
                                            lineheight = 1.3,
                                            margin = ggplot2::margin(t = 10))
  )
  # Legend-free panels + tags
  .f8_body_A <- .fig8_A + ggplot2::theme(legend.position = "none") + .f8_tag + .f8_cap_theme + ggplot2::labs(tag = "A")
  .f8_body_B <- .fig8_B + ggplot2::theme(legend.position = "none") + .f8_tag + .f8_cap_theme + ggplot2::labs(tag = "B")
  .f8_body_C <- .fig8_C + ggplot2::theme(legend.position = "none") + .f8_tag + .f8_cap_theme + ggplot2::labs(tag = "C")
  # Bottom-row panels: halve bottom margin (10→5 mm) to reduce gap to legend
  .f8_body_D <- .fig8_D + ggplot2::theme(legend.position = "none",
    plot.margin = ggplot2::margin(10, 10, 5, 10, unit = "mm")) + .f8_tag + .f8_cap_theme + ggplot2::labs(tag = "D")
  .f8_body_E <- .fig8_E + ggplot2::theme(legend.position = "none",
    plot.margin = ggplot2::margin(10, 10, 5, 10, unit = "mm")) + .f8_tag + .f8_cap_theme + ggplot2::labs(tag = "E")
  # Panel F (2026-08-09): entropy_rate, graded alphabet, raster-embedded from
  # the sequence module -- fills the slot left empty when the fifth collective
  # old panel C) was removed, restoring the 3x2 grid.
  .f8_body_F <- .fig8_F + ggplot2::theme(legend.position = "none",
    plot.margin = ggplot2::margin(10, 10, 5, 10, unit = "mm")) + .f8_tag + .f8_cap_theme + ggplot2::labs(tag = "F")
  # 3x2 panel grid (restored 2026-08-09 -- see panel F above).
  .f8_grid <- (patchwork::wrap_elements(full = .f8_body_A) |
               patchwork::wrap_elements(full = .f8_body_B) |
               patchwork::wrap_elements(full = .f8_body_C)) /
              (patchwork::wrap_elements(full = .f8_body_D) |
               patchwork::wrap_elements(full = .f8_body_E) |
               patchwork::wrap_elements(full = .f8_body_F))
  # Figure_10_legend_panelAB centred below grid (0.5x reduced height); 0.5 cm white frame
  .fig8 <- (patchwork::wrap_elements(full = .f8_grid) /
             .legend_ab_elem) +
             patchwork::plot_layout(heights = c(6, 0.5)) +
             patchwork::plot_annotation(
               theme = ggplot2::theme(
                 plot.background = ggplot2::element_rect(fill = "white", colour = NA),
                 plot.margin     = ggplot2::margin(5, 5, 5, 5, unit = "mm")
               )
             )
  .mg_save(.fig8, "Figure_9", width_mm = 320, height_mm = 220)
} else {
  ts_msg("  Figure 9 skipped — missing component plot(s)")
}
writeLines(c(
  "Figure 8. Trial-level collective movement and transition indicators",
  "(N = 8 trials per treatment).",
  "",
  "(A) School transition rate between main zones (switch_rate; transitions",
  "    per minute of the flow/calm bin-state sequence), per trial. Gaussian",
  "    LMM, treatment fixed effect. Sourced from the sequence (bin-state)",
  "    analysis (behavioural_sequence_analysis_choice_exp.R), not refitted",
  "    here; replaces the earlier switches_per_session count (F(1,8) = 0.105,",
  "    p = 0.754, ns), a cruder, non-rate-normalised measure of the same",
  "    behaviour.",
  "(B) Mean nearest-neighbour distance (NND, cm). Gaussian LMM.",
  "(C) Mean inter-individual distance (IID, cm). Gaussian LMM.",
  "(D) School area (cm2). Gaussian LMM.",
  "(E) School speed (cm/s). Gaussian LMM.",
  "(F) Normalised entropy rate of the graded (Low/Med/High) engagement-state",
  "    sequence (0-1), per trial. Gaussian LMM, treatment fixed effect.",
  "    Sourced from the sequence analysis, as panel A.",
  "",
  "Bars show group mean +/- 1 SE. Points are raw trial-level values (jittered).",
  "CLD letters denote Tukey-adjusted pairwise groupings (max 2 characters).",
  "p-values from LMMs with treatment as fixed effect and fish density as",
  "random effect (RE-selection by AICc); Kenward-Roger df correction applied.",
  "* p < 0.05; ** p < 0.01; *** p < 0.001."
), file.path(.mg_dir, "figure8_caption.txt"))
ts_msg("  Figure 8 caption written.")


# ============ FIGURE 9 section -> saves Figure_10 (6579-7108) ============
# ---- Figure 9 — Interval-resolved indicators with Tukey brackets -----------
# (was Figure 10; standalone transitions Figure 9 deleted)
# Panels: A NND, B IID, C Hull. Centroid speed dropped
# (no significant ANOVA term and no significant Tukey contrast).
# Each panel is a line plot (mean ± SE per timepoint × treatment); significant
# Tukey contrasts (p < 0.05) are annotated with brackets and asterisks.

# (helper function bundle moved up above -- see 'Line-with-Tukey panel helpers' near FIGURE 7 section)
.fig9_A <- if (exists("res_nnd_tp")) .mg_make_line_with_tukey(
  df_src        = if (!is.null(df_gd) && "mean_nnd_cm" %in% names(df_gd)) df_gd else NULL,
  y_col         = "mean_nnd_cm",
  y_label       = "Mean nearest-neighbour distance (cm)",
  label_dir     = "nnd_timepoint",
  anova_caption = .fig10_cap_csv("Mean NND"),
  remap_cld     = TRUE,
  # POSITION nudges only. The label_override column was removed 2026-08-16: it
  # hard-coded the letter text, so the panel stopped tracking its own CLD. It
  # happened to still agree with the refitted model here, but the identical
  # construct on the school-area panel did not (see .fig9_C), and a letter that
  # is typed by hand rather than read from cld_treatmentxtimepoint_f.csv cannot
  # be trusted to survive a re-fit.
  cld_nudges    = data.frame(
    treatment      = c("control",  "control",  "control",
                       "exercise choice", "exercise choice", "exercise choice"),
    tp_num         = c(1L, 2L, 3L,
                       1L, 2L, 3L),
    x_adj          = c(NA_real_, NA_real_, NA_real_,
                       -0.22,    NA_real_, NA_real_),
    y_at_mean      = c(NA, NA, NA,
                       TRUE, NA, NA),
    stringsAsFactors = FALSE)
) else NULL

# Panel removed 2026-08-09 on the same methodological ground as the Figure 9
# removal above: an alignment/heading indicator is the only interval-resolved
# whose computation requires PERSISTENT IDENTITY TRACKING (heading = per-
# fish_id frame-to-frame displacement), so it is dropped from the manuscript.
# res_pol_tp is still fitted upstream and remains in the pipeline's own CSV/
# Word outputs. Panels renumbered: old C/D (IID/school area) -> new B/C.

# label_dir was NULL until 2026-08-09, which suppressed this panel's CLD letters
# even though iid_timepoint/contrasts_treatmentxtimepoint.csv holds three
# significant within-exercise-choice Tukey contrasts (intervals 1-2, 1-3 and
# 2-3, all p = 0.046) — the very contrasts the figure legend already describes.
# The letters are now drawn so the panel matches its own caption and the
# contrast table.
.fig9_B <- if (exists("res_iid_tp")) .mg_make_line_with_tukey(
  df_src        = if (!is.null(df_gd) && "mean_iid_cm" %in% names(df_gd)) df_gd else NULL,
  y_col         = "mean_iid_cm",
  y_label       = "Mean inter-individual distance (cm)",
  label_dir     = "iid_timepoint",
  anova_caption = .fig10_cap_csv("Mean IID"),
  remap_cld     = TRUE
  # collapse_tp = 1L REMOVED 2026-08-16. It existed because interval 1 rendered
  # as "c" above "ab" and read as a difference -- but that "ab" was the
  # TRUNCATED form of "abc", and the truncation is the bug now fixed in
  # .mg_read_cld_tp(). With full letters the cell shows "abc", which shares "c"
  # with the other treatment and is read correctly by the ordinary CLD
  # convention. Collapsing to a single label was also hiding the
  # exercise-choice group outright, which loses more than it fixed.
) else NULL

.fig9_C <- if (exists("res_hull_tp")) .mg_make_line_with_tukey(
  df_src        = if (!is.null(df_gd) && "mean_hull_area_cm2" %in% names(df_gd)) df_gd else NULL,
  y_col         = "mean_hull_area_cm2",
  y_label       = "School area (cm²)",
  label_dir     = "hull_area_timepoint",
  anova_caption = .fig10_cap_csv("Mean school area"),
  remap_cld     = TRUE,
  # POSITION nudges only.
  #   Interval 1: single label at mean(ctrl, ec) height + 0.11*yrange, offset
  #               left; the control label is suppressed and the exercise-choice
  #               one stands for both (they carry the same group).
  #   Intervals 2-3: per-treatment labels at default error-bar-top positions.
  #
  # label_override REMOVED 2026-08-16 -- it was asserting a result the model does
  # not support. The hand-typed letters were "a" (control) vs "b" (exercise) at
  # 85-105 min, i.e. a difference; the actual CLD gives "ab" and "a", which SHARE
  # a letter. The family-adjusted contrast for that cell pair is p = 0.053, so
  # the two do not differ once the 15 comparisons are accounted for -- the
  # override was displaying the unadjusted simple-effect conclusion (p = 0.005)
  # under a caption that says letters are the compact letter display. Letters now
  # come from cld_treatmentxtimepoint_f.csv, as the caption claims.
  cld_nudges    = data.frame(
    treatment      = c("control",  "exercise choice",
                       "control",  "exercise choice",
                       "control",  "exercise choice"),
    tp_num         = c(1L, 1L,
                       2L, 2L,
                       3L, 3L),
    x_adj          = c(NA_real_, -0.22,
                       NA_real_, NA_real_,
                       NA_real_, NA_real_),
    force_separate = c(TRUE, TRUE,
                       TRUE, TRUE,
                       TRUE, TRUE),
    suppress       = c(TRUE, NA,
                       NA,   NA,
                       NA,   NA),
    y_at_treat_avg = c(NA, TRUE,
                       NA, NA,
                       NA, NA),
    stringsAsFactors = FALSE)
) else NULL


# ---- Panel D: commitment_index by treatment x interval (bout-structure
# engine, 2026-08-09). Sourced from that engine's own session-level export
# (STEP5's df_gd only holds NND/IID/area/speed) but rendered with the SAME
# .mg_make_line_with_tukey() function as A-C, so the aesthetic matches
# exactly -- same theme, point/line sizes, colour/shape scales, panel frame.
# commitment_index is EXCLUDED from formal hypothesis testing by the bout
# engine's own collinearity gate (school-level |r| ~ 0.98 with both
# switch_rate and entropy_rate, already reported elsewhere -- see
# metric_collinearity.csv), so no CLD letters are drawn (label_dir = NULL)
# and the caption states that plainly instead of reporting a claimed F/p.
.bout_session_path <- .find_latest_csv("BOUT_output", "bout_metrics_per_session.csv")
.df_bout_session <- if (!is.null(.bout_session_path)) {
  d <- readr::read_csv(.bout_session_path, show_col_types = FALSE)
  d$treatment   <- factor(trimws(tolower(as.character(d$treatment))), levels = TREATMENT_LEVELS_g)
  d$timepoint_f <- factor(d$timepoint, levels = TIMEPOINT_LEVELS_g)
  d
} else NULL

# 2026-08-10: panel D now carries CLD letters and its own Treatment x Interval
# ANOVA statement, matching panels A-C. The letters come from a dedicated
# emmeans refit of the SAME Level-B model the bout engine reports
# (val ~ treatment * timepoint_f + fish_density + (1|trial)), written to
# easy_scripts/standalone/cld_commitment_index_treatmentxtimepoint.csv with
# letters remapped so control/interval 1 carries "a". NOTE: emmeans substitutes
# Sidak for Tukey on this 6-cell family -- but VERIFIED 2026-08-17 that the
# substitution applies to the CONFIDENCE-INTERVAL adjustment only. cld() reports
# both: "Conf-level adjustment: sidak" and "P value adjustment: tukey". The
# p-values that decide the letter groupings are Tukey; the letters are therefore
# Tukey-adjusted, not Sidak-adjusted.
.cld_commit_path <- file.path(file.path(PROJECT_ROOT, "choice R pipeline/easy_scripts/standalone"),
                              "cld_commitment_index_treatmentxtimepoint.csv")
.bout_tp_path <- .find_latest_csv("BOUT_output", "bout_anova_trial_x_timepoint.csv")
.commit_cap <- local({
  if (is.null(.bout_tp_path)) return("")
  d <- readr::read_csv(.bout_tp_path, show_col_types = FALSE)
  r <- d[d$metric == "commitment_index" & d$term == "treatment:timepoint_f", ]
  if (!nrow(r)) return("")
  p_txt <- if (r$p[1] < 0.001) "p < 0.001" else sprintf("p = %.3f", r$p[1])
  sprintf("Treatment\u00d7Interval *F*<sub>%s,%s</sub> = %s, %s, \u03b7\u00b2<sub>p</sub> = %s",
          format(r$df1[1]), format(round(r$df2[1], 2)),
          format(round(r$F[1], 3)), p_txt, format(round(r$eta2_p[1], 3)))
})

.fig9_D <- if (!is.null(.df_bout_session)) .mg_make_line_with_tukey(
  df_src        = .df_bout_session,
  y_col         = "commitment_index",
  y_label       = "Commitment index",
  label_dir     = NULL,
  cld_path      = .cld_commit_path,
  anova_caption = .commit_cap,
  # Was FALSE, relying on the static CSV having been hand-remapped. Now applied
  # in code so the invariant is enforced rather than trusted (idempotent if the
  # file was already correct).
  remap_cld     = TRUE
) else NULL

.fig9_panels <- list(.fig9_A, .fig9_B, .fig9_C, .fig9_D)
.fig9_panels <- lapply(.fig9_panels, function(p) if (is.null(p)) NULL else .mg_frame(p))
if (!any(sapply(.fig9_panels, is.null))) {
  # Legend-free panels + tags; panel B (IID) retains scale overrides (suppress duplicate legend)
  .f9_body_A <- .fig9_panels[[1]] + ggplot2::theme(legend.position = "none") + .FIG10_CAP_OVR + .MG_TAG_THEME + ggplot2::labs(tag = "A")
  .f9_body_B <- .fig9_panels[[2]] + ggplot2::theme(legend.position = "none") + .FIG10_CAP_OVR +
    ggplot2::scale_colour_manual(values = TREATMENT_COLORS_PAL, name = "Treatment", guide = "none") +
    ggplot2::scale_linetype_manual(values = c("exercise choice" = "solid", "control" = "dashed"),
      name = "Treatment", guide = "none") +
    .MG_TAG_THEME + ggplot2::labs(tag = "B")
  .f9_body_C <- .fig9_panels[[3]] + ggplot2::theme(legend.position = "none") + .FIG10_CAP_OVR + .MG_TAG_THEME + ggplot2::labs(tag = "C")
  .f9_body_D <- .fig9_panels[[4]] + ggplot2::theme(legend.position = "none") + .FIG10_CAP_OVR + .MG_TAG_THEME + ggplot2::labs(tag = "D")
  # True 2x2 grid (2026-08-09): the trailing plot_spacer() that used to hold
  # this cell (kept from when the fifth collective panel was removed, to preserve
  # each panel's ~130mm width and the already-verified 11pt caption size) is
  # now filled by panel D -- same cell, same width, no re-verification needed.
  .f9_grid <- (patchwork::wrap_elements(full = .f9_body_A) |
               patchwork::wrap_elements(full = .f9_body_B)) /
              (patchwork::wrap_elements(full = .f9_body_C) |
               patchwork::wrap_elements(full = .f9_body_D))
  # Figure_10_legend_panelAB centred below grid (height decreased by 70% from 0.7); 0.5 cm white frame
  .fig9_iv <- (patchwork::wrap_elements(full = .f9_grid) /
                .legend_ab_elem) +
                patchwork::plot_layout(heights = c(6, 0.5)) +
                patchwork::plot_annotation(
                  theme = ggplot2::theme(
                    plot.background = ggplot2::element_rect(fill = "white", colour = NA),
                    plot.margin     = ggplot2::margin(5, 5, 5, 5, unit = "mm")
                  )
                )
  .mg_save(.fig9_iv, "Figure_10", width_mm = 260, height_mm = 220)
  writeLines(c(
    "Figure 9. Interval-resolved collective movement indicators with significant",
    "treatment-related terms or Tukey contrasts (N = 8 trials per treatment",
    "per interval; intervals: 5-25 min, 45-65 min, 85-105 min).",
    "",
    "(A) Mean nearest-neighbour distance (NND, cm).",
    "    Significant Treatment x Interval interaction (LMM with F-KR).",
    "(B) Mean inter-individual distance (IID, cm).",
    "    Significant within-exercise-choice Tukey contrasts across intervals.",
    "(C) School area (cm2).",
    "    Significant Tukey contrast (control Interval 2 vs exercise choice",
    "    Interval 3).",
    "(D) Commitment index (longest flow bout / interval duration), from the",
    "    bout-structure analysis. Collinear with switch_rate and entropy_rate",
    "    at school level (|r| ~ 0.98, already reported in Figure 9 panels A",
    "    and F), so it is read as a re-expression of those outcomes rather",
    "    than as independent evidence; its Treatment x Interval term and CLD",
    "    letters are shown for completeness. Letters are Tukey-adjusted across",
    "    the six treatment x interval cells, as in panels A-C.",
    "",
    "Lines show group mean +/- 1 SE per interval. Compact letter display",
    "(CLD) labels above each cell denote Tukey-adjusted pairwise groupings",
    "(treatment x interval); cells sharing a letter do not differ at p < 0.05.",
    "When both treatments at the same interval share the same letter(s),",
    "a single label is shown. Letters are shown in full and are relabelled so",
    "that control at interval 1 always carries 'a' (a relabelling only: which",
    "cells share a letter is unchanged). Subtitle shows the Treatment x Interval",
    "F-statistic from the LMM for every panel, including panel D."
  ), file.path(.mg_dir, "figure9_caption.txt"))
  ts_msg("  Figure 9 (interval) caption written.")
} else {
  ts_msg("  Figure 9 (interval) skipped — missing component plot(s)")
}


# ============ FIGURE 7bis section (7109-7152) ============
# ---- Figure 7bis: alr(flow) by treatment × interval (single panel) ----------
# .fig7bis itself is built earlier now (with the Figure 7/8 panels, so it can
# double as Figure_8 panel F) -- this just saves it standalone and writes its
# own caption, unchanged from before.
if (!is.null(.fig7bis)) {
  .mg_save(.fig7bis, "figure7bis_alr_flow_timepoint", width_mm = 140, height_mm = 130)
  writeLines(c(
    "Figure 7bis. Additive log-ratio of flow vs calm occupancy (alr(flow)) by",
    "treatment across observation intervals (N = 8 trials per treatment per",
    "interval; intervals: 5-25 min, 45-65 min, 85-105 min).",
    "",
    "Lines connect treatment group means; error bars = +/- 1 SE. Points are",
    "trial-level values. CLD letters denote Tukey-adjusted pairwise groupings",
    "across all treatment x interval cells (precedence: control > exercise choice;",
    "interval 1 > 2 > 3). ANOVA caption emitted only when the treatment x",
    "interval interaction term is significant (p < 0.05; 3 d.p.)."
  ), file.path(.mg_dir, "figure7bis_caption.txt"))
  ts_msg("  Figure 7bis saved and caption written.")
} else {
  ts_msg("  Figure 7bis skipped — df_main_wide not available")
}


# =============================================================================
# ==== MANUSCRIPT INTERVAL FIGURES -- final indicator set =====================
# =============================================================================
# The paper reports behavioural outcomes at the trial x interval level ONLY
# (N = 48 sessions; 16 trials x 3 intervals). The aggregated / trial-level module
# is not reported. Interval is a three-level factor throughout, so every
# interval and treatment x interval term is a 2-df test.
#
# Three figures, matching the three questions the paper asks:
#
#   FIG_A  Zone preference          2x2  ALR flow / high / medium / low vs calm
#          -> "do fish spend time in the flow rather than the calm?" and is the
#             preference graded across the velocity gradient?
#
#   FIG_B  Pattern of engagement    2x3  flow<->calm crossings, entropy rate,
#          dwell time (Flow and Calm together), flow bouts per minute, longest
#          flow bout, commitment index
#          -> "long uninterrupted bouts, or short bouts spaced by calm?"
#
#   FIG_C  Collective movement      2x2  NND, IID, school area, school speed
#
# Figure NUMBERS here are placeholders (A/B/C). The manuscript numbering is
# assigned when the figures are placed; do not hard-code 1/2/3 anywhere.
#
# Denominator df are printed as integers in every caption (the fractional
# Kenward-Roger / Satterthwaite value stays in the engine's anova.csv).
#
# Every panel is drawn through the STEP5 panel builders even when the data come
# from the SEQ or BOUT engines, so point sizes, line widths and axis themes
# cannot differ within a figure.
# =============================================================================

# ---- Caption builder for the SEQ / BOUT anova tables ------------------------
# `terms` may name more than one ANOVA row; each becomes its own caption line.
.mg_seqbout_cap <- function(csv_path, metric_name, alphabet = NULL,
                            terms = "treatment:timepoint_f") {
  if (is.null(csv_path) || !file.exists(csv_path)) return("")
  d <- tryCatch(readr::read_csv(csv_path, show_col_types = FALSE),
                error = function(e) NULL)
  if (is.null(d) || !all(c("metric", "term") %in% names(d))) return("")
  lbl <- c(treatment = "Treatment", timepoint_f = "Interval",
           `treatment:timepoint_f` = "Treatment×Interval")
  lines <- vapply(terms, function(tm) {
    r <- d[d$metric == metric_name & d$term == tm, , drop = FALSE]
    if (!is.null(alphabet) && "alphabet" %in% names(r))
      r <- r[r$alphabet == alphabet, , drop = FALSE]
    if (!nrow(r)) return("")
    r <- r[1, , drop = FALSE]
    if (!is.finite(r$p[1])) return("")
    p_txt <- if (r$p[1] < 0.001) "p < 0.001" else sprintf("p = %.3f", r$p[1])
    eta <- if ("eta2_p" %in% names(r) && is.finite(r$eta2_p[1]))
      sprintf(", η²<sub>p</sub> = %s", format(round(r$eta2_p[1], 3))) else ""
    # Denominator df as an INTEGER; fractional value retained in the anova CSV.
    sprintf("%s *F*<sub>%s,%s</sub> = %s, %s%s",
            lbl[[tm]] %||% tm,
            format(round(r$df1[1])), format(round(r$df2[1])),
            format(round(r$F[1], 3)), p_txt, eta)
  }, character(1))
  paste(lines[nzchar(lines)], collapse = "<br>")
}

# ---- Caption builder reading a STEP5 model's anova.csv DIRECTLY --------------
# The static anova_statement_table_interval.csv carries only the interaction
# term, and only for the eight indicators listed in its MAP. This reads whatever
# terms are asked for straight from the model directory, so a panel can show the
# treatment main effect beside the interaction without the caption table needing
# an entry. Denominator df printed as integers, matching .mg_seqbout_cap().
.mg_step5_cap <- function(label_dir, terms = "treatment:timepoint_f") {
  f <- file.path(STEP5_OUT, label_dir, "anova.csv")
  if (!file.exists(f)) return("")
  a <- tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(a)) return("")
  lbl <- c(treatment = "Treatment", timepoint_f = "Interval",
           `treatment:timepoint_f` = "Treatment×Interval")
  is_F <- length(a$stat_type) && grepl("^F", a$stat_type[1])
  lines <- vapply(terms, function(tm) {
    r <- a[a$term == tm, , drop = FALSE]
    if (!nrow(r) || !is.finite(r$p_value[1])) return("")
    p_txt <- if (r$p_value[1] < 0.001) "p < 0.001" else sprintf("p = %.3f", r$p_value[1])
    d1 <- r$df[1]; d2 <- r$df_denom[1]; st <- r$chisq[1]
    eta <- if (is.finite(d2)) {
      e <- (st * d1) / (st * d1 + d2)
      sprintf(", η²<sub>p</sub> = %s", format(round(e, 3)))
    } else ""
    if (is_F && is.finite(d2))
      sprintf("%s *F*<sub>%s,%s</sub> = %s, %s%s", lbl[[tm]] %||% tm,
              format(round(d1)), format(round(d2)), format(round(st, 3)), p_txt, eta)
    else
      sprintf("%s χ²<sub>%s</sub> = %s, %s", lbl[[tm]] %||% tm,
              format(round(d1)), format(round(st, 3)), p_txt)
  }, character(1))
  paste(lines[nzchar(lines)], collapse = "<br>")
}

# NOTE: a combined two-state dwell panel (Flow and Calm on shared axes, colour
# = state) was built here on 2026-08-17 and removed the same day. Dwell time is
# now shown as two separate panels using the standard treatment palette, so
# colour means the same thing in every panel of every figure.

# ---- data sources -----------------------------------------------------------
.seq_bin_path <- .find_latest_csv("SEQ_output", "seq_metrics_per_trial_binary.csv")
.df_seq_bin <- if (!is.null(.seq_bin_path)) {
  d <- readr::read_csv(.seq_bin_path, show_col_types = FALSE)
  d$treatment   <- factor(trimws(tolower(as.character(d$treatment))),
                          levels = TREATMENT_LEVELS_g)
  d$timepoint_f <- factor(d$timepoint, levels = TIMEPOINT_LEVELS_g)
  d
} else NULL
.seq_tp_path   <- .find_latest_csv("SEQ_output",  "seq_anova_trial_x_timepoint.csv")
.bout_tp_path2 <- if (exists(".bout_tp_path") && !is.null(.bout_tp_path)) .bout_tp_path else
  .find_latest_csv("BOUT_output", "bout_anova_trial_x_timepoint.csv")

# Every panel caption reports the treatment main effect and the treatment x
# interval interaction, on two lines.
.TERMS2 <- c("treatment", "treatment:timepoint_f")

# The compact-letter-display convention, written once and spliced into all three
# figure captions so they cannot describe the letters differently.
.CLD_CAPTION_NOTE <- c(
  "Compact letter display: cells sharing a letter do not differ at p < 0.05,",
  "across the six treatment x interval cells. Letters are relabelled so that",
  "control at the first interval always begins with 'a'; this is a relabelling",
  "only, and which cells share a letter is unchanged. A cell carrying two",
  "letters (e.g. 'ab') belongs to two non-difference groups -- it differs from",
  "neither -- and cannot be reduced to one letter without asserting a",
  "difference the model does not support.",
  "",
  "PANELS WITHOUT LETTERS have no significant pairwise differences at all: every",
  "cell falls in the same group, so letters would carry no information and are",
  "omitted rather than printed as six identical marks."
)

.mg_body <- function(p, tg, keep_legend = FALSE)
  p + (if (keep_legend) ggplot2::theme() else ggplot2::theme(legend.position = "none")) +
    .FIG10_CAP_OVR + .MG_TAG_THEME + ggplot2::labs(tag = tg)

.mg_finish <- function(grid, legend_elem = .legend_ab_elem)
  (patchwork::wrap_elements(full = grid) / legend_elem) +
    patchwork::plot_layout(heights = c(6, 0.5)) +
    patchwork::plot_annotation(theme = ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.margin     = ggplot2::margin(5, 5, 5, 5, unit = "mm")))

# =============================================================================
# FIGURE A -- zone preference across the velocity gradient (2 x 2)
# =============================================================================
# Captions come from .mg_step5_cap() for ALL FOUR panels, not from the static
# anova_statement_table_interval.csv. That table holds only the interaction term
# and only for indicators in its MAP, which would have given panels A/B a
# one-line caption and panels C/D a two-line one inside a single figure.
.figA_specs <- list(
  # Panel A: at interval 2 control is "a" and exercise choice "ab", so two
  # separate labels are drawn. Control's default position (its own +SE tip +
  # 0.10*yrange) lands on the exercise-choice marker and error bar, which
  # partially covers the letter. Shift it left only -- the vertical position is
  # left alone so it stays consistent with every other label in the panel.
  list(df = "df_main_wide", y = "logit_flow", lab = "ALR(flow vs. calm)",
       dir = "zone_main_timepoint",
       nudge = data.frame(treatment = "control", tp_num = 2L, x_adj = -0.22,
                          stringsAsFactors = FALSE)),
  # Panel B: both treatments carry "a" at interval 1, so .mg_cld_label_rows()
  # MERGES them into one label with treatment = NA. A merged row cannot be
  # targeted (y_at_mean is guarded by !is.na(treatment)), so force_separate
  # splits interval 1 first -- it is keyed per interval, hence set on BOTH rows.
  # Control then moves left of its mean point at the mean's own height;
  # exercise choice keeps the default. Interval 2 is not merged (its two cells
  # carry different strings); it only needed a nudge -- see the note below.
  list(df = "df_sec_wide",  y = "lr_high",   lab = "ALR(high vs. calm)",
       dir = "zone_sec_high_timepoint",
       # 2026-08-18: control/interval 2 letter dropped at the author's request
       # (suppress). 2026-08-21: RESTORED -- MV read the gap as a missing letter
       # ("Figure B is missing a letter in the middle interval"), and the CLD file
       # does letter that cell ("ab"). It was only ever suppressed because it
       # collided with the exercise-choice marker, so it now gets the same remedy
       # as interval 1: shifted left and dropped to its own mean height. No
       # force_separate is needed at interval 2 -- control "ab" and exercise
       # choice "a" are different strings, so .mg_cld_label_rows() never merged
       # them; only the collision hid it.
       nudge = data.frame(
         treatment      = c("control", "exercise choice", "control"),
         tp_num         = c(1L, 1L, 2L),
         force_separate = c(TRUE, TRUE, NA),
         x_adj          = c(-0.22, NA_real_, -0.22),
         y_at_mean      = c(TRUE,  NA, TRUE),
         stringsAsFactors = FALSE)),
  list(df = "df_sec_wide",  y = "lr_medium", lab = "ALR(medium vs. calm)",
       dir = "zone_sec_medium_timepoint"),
  list(df = "df_sec_wide",  y = "lr_low",    lab = "ALR(low vs. calm)",
       # 2026-08-18: at interval 2 the two cell means are close (control 0.46,
       # exercise choice 0.22), so both letters ("ab" and "a") were drawn at
       # essentially the same point and overprinted into an illegible glyph --
       # it read as a single label, i.e. one letter apparently missing. Same
       # remedy as panel B interval 1: shift control left and drop it to its own
       # mean height so the two letters separate.
       dir = "zone_sec_low_timepoint",
       # 2026-08-21: interval 1 added for the same reason, but a DIFFERENT cause --
       # here both cells carry the identical string "ab", which trips the
       # length(unique(rows$.group)) == 1L branch and collapses them into a single
       # centred label. MV read that as "figure D is missing a letter in the first
       # interval". force_separate splits the pair; control then moves left at its
       # own mean height, exercise choice keeps the default position.
       # 2026-08-21 (round 2): the exercise-choice letter at interval 2 moves down
       # to the BOTTOM of its own error bar. Its marker (0.22) sits just under
       # control's (0.46), so the default "above my own SE tip" position put the
       # letter inside control's error bar.
       nudge = data.frame(
         treatment      = c("control", "exercise choice", "control", "exercise choice"),
         tp_num         = c(2L, 2L, 1L, 1L),
         force_separate = c(TRUE, TRUE, TRUE, TRUE),
         x_adj          = c(-0.22, NA_real_, -0.22, NA_real_),
         y_at_mean      = c(TRUE,  NA, TRUE,  NA),
         y_at_lower_se  = c(NA,    TRUE, NA,   NA),
         stringsAsFactors = FALSE))
)
.figA_panels <- lapply(.figA_specs, function(s) {
  d <- get0(s$df, ifnotfound = NULL)
  if (is.null(d) || !s$y %in% names(d)) return(NULL)
  .mg_make_line_with_tukey(df_src = d, y_col = s$y, y_label = s$lab,
                           label_dir = s$dir,
                           anova_caption = .mg_step5_cap(s$dir, terms = .TERMS2),
                           remap_cld = TRUE,
                           cld_nudges = s$nudge)
})
if (!any(vapply(.figA_panels, is.null, logical(1)))) {
  b <- Map(.mg_body, lapply(.figA_panels, .mg_frame), c("A", "B", "C", "D"))
  .figA_grid <- (patchwork::wrap_elements(full = b[[1]]) |
                 patchwork::wrap_elements(full = b[[2]])) /
                (patchwork::wrap_elements(full = b[[3]]) |
                 patchwork::wrap_elements(full = b[[4]]))
  .mg_save(.mg_finish(.figA_grid), "Figure_A_zone_preference",
           width_mm = 280, height_mm = 240)
  writeLines(c(
    "Figure A. Zone preference across the velocity gradient, by treatment and",
    "observation interval (N = 8 trials per treatment per interval; intervals",
    "5-25, 45-65 and 85-105 min after trial onset).",
    "",
    "(A) ALR(flow vs calm)   (B) ALR(high vs calm)",
    "(C) ALR(medium vs calm) (D) ALR(low vs calm)",
    "",
    "Each panel is an additive log-ratio of area-normalised occupancy: the log",
    "of time in the named zone over time in the calm zone, after dividing each",
    "by that zone's share of arena area. Zero means the zone is used exactly in",
    "proportion to its area; positive means over-used relative to calm. The",
    "log-ratio is used rather than a raw proportion because time in one zone",
    "necessarily denies time to the other.",
    "",
    "Panels B-D partition the flowing side of panel A into its three velocity",
    "bands, so together they show WHERE on the gradient any preference sits.",
    "",
    "Interval is fitted as a three-level factor, so every Treatment x Interval",
    "term is a 2-df test. Denominator degrees of freedom are printed as",
    "integers; the fractional Kenward-Roger values are retained in the",
    "analysis output.",
    "",
    "Lines join treatment group means; error bars = +/- 1 SE.",
    "",
    .CLD_CAPTION_NOTE
  ), file.path(.mg_dir, "figure_A_caption.txt"))
  ts_msg("  Figure A (zone preference, 2x2) written.")
} else {
  ts_msg("  Figure A skipped -- missing panel(s): ",
         paste(which(vapply(.figA_panels, is.null, logical(1))), collapse = ", "))
}

# =============================================================================
# FIGURE B -- pattern of engagement with the flow (2 x 2)
# =============================================================================
# Layout, in reading order:
#   A  flow<->calm crossings      B  flow bouts per minute  <- how OFTEN
#   C  longest flow bout          D  longest calm bout      <- how LONG, at the extreme
#
# Rebuilt 2026-08-17 (was 2 wide x 3 tall). Entropy rate and the commitment
# index left the paper entirely; both dwell-time panels were dropped from the
# figure but stay in Table 1 and the Results text, so the figure now carries one
# rate measure and one extreme measure per state.
#
# WHY BOTH A MEAN AND A MAXIMUM ARE REPORTED (the mean lives in the text, the
# maximum in panels C/D). Dwell time IS mean bout duration -- bout_structure_
# analysis_choice_exp.R asserts mean_bout_flow == dwell_Flow to 1e-8 -- so the
# two are statistics of ONE duration distribution, and they differ in what they
# can see. The mean and the per-minute count are computed on UNCENSORED INTERIOR
# bouts only (the first and last run of each segment are dropped as boundary-
# censored); the maximum is taken over the FULL inventory. A visit still running
# when the interval ends is therefore invisible by construction to the mean and
# the count, and visible only to the maximum. That is why the flow result --
# mean flat (p = 0.458), maximum p < 0.001 -- is not a contradiction: it is the
# signature of engagement consolidating into a few sustained visits.
#
# Panel D is null on both terms and is kept anyway. It is the calm-side
# counterpart to C, and the PAIR is what licenses the claim that the lengthening
# is specific to the flow rather than a general lengthening of bouts. It is also
# defined for all 48 sessions (8/8/8 in both arms), unlike the bout-shape
# measures.
#
# CLD on every panel. The SEQ and BOUT engines write no compact letter display,
# so theirs is built by build_seqbout_cld.R -- a refit of the engine's own
# Level-B model, lettered under the shared letter policy -- and read here through
# the ordinary `cld_path` argument.
.SEQBOUT_CLD <- file.path(.find_latest_csv_root, "SEQBOUT_cld")
.seqbout_cld <- function(metric) {
  p <- file.path(.SEQBOUT_CLD, metric, "cld_treatmentxtimepoint_f.csv")
  if (file.exists(p)) p else NULL
}

.figB <- list()

# B-A flow<->calm crossings (zone flux)
.figB$crossings <- if (!is.null(df) && "zone_flux_per_session" %in% names(df))
  .mg_make_line_with_tukey(
    df_src = df, y_col = "zone_flux_per_session",
    y_label = "Flow↔calm crossings (n)", label_dir = "flux_timepoint",
    anova_caption = .mg_step5_cap("flux_timepoint", terms = .TERMS2),
    # 2026-08-21 (round 2): control and exercise choice both carry "a" at
    # interval 2, so .mg_cld_label_rows() merged them into one centred label and
    # the panel looked like the exercise-choice cell had lost its letter. Every
    # cell gets its own label. The two markers are far apart vertically here
    # (91.7 vs 50.8), so no x nudge is needed to keep them legible.
    remap_cld = TRUE,
    cld_nudges = data.frame(
      treatment      = c("control", "exercise choice"),
      tp_num         = c(2L, 2L),
      force_separate = c(TRUE, TRUE),
      stringsAsFactors = FALSE)) else NULL

# B-B flow bouts per minute -- how OFTEN the flow is entered
#
# At interval 1 the exercise-choice letter ("ad") would default to its own +SE
# tip + 0.10*yrange = ~0.62, which is exactly where the CONTROL group's lower
# error-bar whisker ends (~0.66): the letter is drawn on top of the other
# series' error bar. Control's "ab" sits clear above its own marker and is left
# alone. Same remedy as Figure A panel B and Figure C panel A -- left of its own
# mean point, at that mean's height.
.figB$rate <- if (!is.null(.df_bout_session) && "bout_rate_flow" %in% names(.df_bout_session))
  .mg_make_line_with_tukey(
    df_src = .df_bout_session, y_col = "bout_rate_flow",
    y_label = "Flow bouts per minute", label_dir = NULL,
    cld_path = .seqbout_cld("bout_rate_flow"),
    anova_caption = .mg_seqbout_cap(.bout_tp_path2, "bout_rate_flow",
                                    terms = .TERMS2),
    remap_cld = TRUE,
    # 2026-08-21 (round 2): control/interval 1 belongs to two non-difference
    # groups here, so no labelling can give it a plain "a" and the panel had no
    # single "a" anywhere. "lone_first" puts "a" on the nearest cell that does
    # carry one letter -- control/interval 2 -- and renames the rest in reading
    # order. Cosmetic only; which cells share a letter is unchanged.
    remap_mode = "lone_first",
    cld_nudges = data.frame(
      treatment = "exercise choice", tp_num = 1L,
      x_adj = -0.22, y_at_mean = TRUE,
      stringsAsFactors = FALSE)) else NULL

# B-C longest flow bout -- the single most sustained visit to the flow
.figB$max_flow <- if (!is.null(.df_bout_session) && "max_flow_bout_s" %in% names(.df_bout_session))
  .mg_make_line_with_tukey(
    df_src = .df_bout_session, y_col = "max_flow_bout_s",
    y_label = "Longest flow bout (s)", label_dir = NULL,
    cld_path = .seqbout_cld("max_flow_bout_s"),
    anova_caption = .mg_seqbout_cap(.bout_tp_path2, "max_flow_bout_s",
                                    terms = .TERMS2),
    remap_cld = TRUE) else NULL

# B-D longest calm bout -- the calm-state counterpart to C. Null on both terms,
# and kept for exactly that reason: it is what makes C's effect flow-SPECIFIC
# rather than a general lengthening of every bout.
.figB$max_calm <- if (!is.null(.df_bout_session) && "max_calm_bout_s" %in% names(.df_bout_session))
  .mg_make_line_with_tukey(
    df_src = .df_bout_session, y_col = "max_calm_bout_s",
    y_label = "Longest calm bout (s)", label_dir = NULL,
    cld_path = .seqbout_cld("max_calm_bout_s"),
    anova_caption = .mg_seqbout_cap(.bout_tp_path2, "max_calm_bout_s",
                                    terms = .TERMS2),
    remap_cld = TRUE) else NULL

if (!any(vapply(.figB, is.null, logical(1))) && length(.figB) == 4L) {
  b <- Map(.mg_body, lapply(.figB, .mg_frame), c("A", "B", "C", "D"))
  .figB_grid <- (patchwork::wrap_elements(full = b[[1]]) |
                 patchwork::wrap_elements(full = b[[2]])) /
                (patchwork::wrap_elements(full = b[[3]]) |
                 patchwork::wrap_elements(full = b[[4]]))
  .mg_save(.mg_finish(.figB_grid), "Figure_B_engagement_pattern",
           width_mm = 280, height_mm = 240)
  writeLines(c(
    "Figure B. Pattern of engagement with the flow, by treatment and",
    "observation interval (N = 8 trials per treatment per interval; intervals",
    "5-25, 45-65 and 85-105 min after trial onset).",
    "",
    "(A) Flow<->calm crossings per session -- traffic across the zone boundary.",
    "(B) Flow bouts per minute -- how OFTEN the school enters the flow.",
    "(C) Longest flow bout (s) -- the single most sustained visit to the flow",
    "    within an interval.",
    "(D) Longest calm bout (s) -- the same measure for rest, shown as the",
    "    calm-state counterpart to C.",
    "",
    "Panels B-D answer the shape of engagement directly: B counts how often the",
    "flow is entered, C how long the most sustained visit lasts. Few entries",
    "with a long maximum indicates engagement consolidated into a few sustained",
    "visits; many entries with a short maximum indicates fragmented engagement.",
    "D is the same question for rest, and its null result is what makes C's",
    "effect specific to the flow rather than a general lengthening of bouts.",
    "",
    "MEAN AND MAXIMUM ARE BOTH REPORTED, and only the maximum is plotted here;",
    "mean bout duration (dwell time) is given in Table 1 and the Results text.",
    "The two are statistics of the same duration distribution but cannot see the",
    "same things: the mean and the per-minute count are computed on uncensored",
    "interior bouts only, whereas the maximum uses the full inventory, so a",
    "visit still running when the interval ends is invisible to the mean by",
    "construction. In the flow the two dissociate -- the typical completed visit",
    "is unchanged while the longest is far longer -- which is the signature of",
    "engagement consolidating rather than of every visit lengthening.",
    "",
    "Bout-SHAPE measures (burstiness, CV of bout duration) are deliberately",
    "absent. They require several flow bouts to estimate a dispersion, and",
    "exercise-choice schools average 1.75 flow bouts by the final interval, so",
    "those metrics are undefined for that entire arm. That is a property of the",
    "behaviour rather than of the measurement: by the last interval most",
    "exercise-choice schools never left the flow zone, so engagement had",
    "stopped consisting of separable bouts.",
    "",
    "Interval is fitted as a three-level factor throughout, so every",
    "Treatment x Interval term is a 2-df test; denominator degrees of freedom",
    "are printed as integers. Lines join treatment group means; error bars =",
    "+/- 1 SE.",
    "",
    .CLD_CAPTION_NOTE,
    "",
    "For panels sourced from the sequence and bout engines, which write no",
    "compact letter display of their own, the letters come from a refit of the",
    "engine's own Level-B model."
  ), file.path(.mg_dir, "figure_B_caption.txt"))
  ts_msg("  Figure B (engagement pattern, 2x2) written.")
} else {
  ts_msg("  Figure B skipped -- missing: ",
         paste(names(.figB)[vapply(.figB, is.null, logical(1))], collapse = ", "))
}

# =============================================================================
# FIGURE C -- collective movement (2 x 2)
# =============================================================================
.figC_specs <- list(
  # Panel A: at interval 1 control is "a" and exercise choice "ab". The letter
  # that collides is the EXERCISE-CHOICE one, not control's: control's "a" sits
  # above its own error-bar tip and is clear, while "ab" -- drawn above the
  # exercise-choice tip -- lands on the control marker, which is the higher of
  # the two groups here. Moving control's "a" down to its mean height would have
  # placed it on "ab" and created a second collision, so the exercise-choice
  # label is the one that moves: left of its own mean point, at that mean's
  # height. Same construct already in production on .fig9_A.
  list(y = "mean_nnd_cm",         lab = "Mean nearest-neighbour distance (cm)",
       dir = "nnd_timepoint",             cap = "Mean NND",
       nudge = data.frame(treatment = "exercise choice", tp_num = 1L,
                          x_adj = -0.22, y_at_mean = TRUE,
                          stringsAsFactors = FALSE)),
  list(y = "mean_iid_cm",         lab = "Mean inter-individual distance (cm)",
       dir = "iid_timepoint",             cap = "Mean IID"),
  list(y = "mean_hull_area_cm2",  lab = "School area (cm²)",
       dir = "hull_area_timepoint",       cap = "Mean school area"),
  list(y = "mean_centroid_spd_cm", lab = "School speed (cm/s)",
       dir = "centroid_speed_timepoint",  cap = "Mean school speed")
)
.figC_panels <- lapply(.figC_specs, function(s) {
  if (is.null(df_gd) || !s$y %in% names(df_gd)) return(NULL)
  .mg_make_line_with_tukey(
    df_src = df_gd, y_col = s$y, y_label = s$lab, label_dir = s$dir,
    anova_caption = .mg_step5_cap(s$dir, terms = .TERMS2), remap_cld = TRUE,
    cld_nudges = s$nudge)
})
if (!any(vapply(.figC_panels, is.null, logical(1)))) {
  b <- Map(.mg_body, lapply(.figC_panels, .mg_frame), c("A", "B", "C", "D"))
  .figC_grid <- (patchwork::wrap_elements(full = b[[1]]) |
                 patchwork::wrap_elements(full = b[[2]])) /
                (patchwork::wrap_elements(full = b[[3]]) |
                 patchwork::wrap_elements(full = b[[4]]))
  .mg_save(.mg_finish(.figC_grid), "Figure_C_collective_movement",
           width_mm = 280, height_mm = 240)
  writeLines(c(
    "Figure C. Collective movement, by treatment and observation interval",
    "(N = 8 trials per treatment per interval).",
    "",
    "(A) Mean nearest-neighbour distance (cm) -- local packing.",
    "(B) Mean inter-individual distance (cm) -- global spread.",
    "(C) School area (cm2) -- convex-hull footprint.",
    "(D) School speed (cm/s) -- centroid displacement.",
    "",
    "All four are computed from the unlabelled set of per-frame positions and",
    "are therefore unaffected by identity error. Panels A-C describe how",
    "tightly the group is arranged; panel D separates that from how fast it",
    "moves, which forecloses the reading that any zone preference is simply",
    "faster-swimming fish being carried into the current.",
    "",
    "Read this figure as a family rather than as four independent tests: all",
    "four collective indicators are reported together, whatever each shows,",
    "because they were selected on construct coverage rather than on their",
    "results.",
    "",
    "Interval is a three-level factor; denominator degrees of freedom are",
    "printed as integers. Lines join treatment group means; error bars =",
    "+/- 1 SE.",
    "",
    .CLD_CAPTION_NOTE
  ), file.path(.mg_dir, "figure_C_caption.txt"))
  ts_msg("  Figure C (collective movement, 2x2) written.")
} else {
  ts_msg("  Figure C skipped -- missing panel(s): ",
         paste(which(vapply(.figC_panels, is.null, logical(1))), collapse = ", "))
}

# ============ MIRROR to all_manu_graphs (7156-7174) ============
# Mirror every figure from manu_graphs to D:/CHOICE R SCRIPTS/all_manu_graphs
.all_manu_dir <- file.path(PROJECT_ROOT, "all_manu_graphs")
dir.create(.all_manu_dir, recursive = TRUE, showWarnings = FALSE)
.mg_files <- list.files(.mg_dir, pattern = "\\.(png|pdf|txt)$",
                        full.names = TRUE, recursive = FALSE)
invisible(file.copy(.mg_files, file.path(.all_manu_dir, basename(.mg_files)),
                    overwrite = TRUE))
ts_msg("Mirrored ", length(.mg_files), " manu_graphs file(s) → ", .all_manu_dir)
# Copy TXT caption files to the dated PNG folder and PDF subfolder.
.mg_txt_files <- list.files(.mg_dir, pattern = "\\.txt$", full.names = TRUE, recursive = FALSE)
if (length(.mg_txt_files) && dir.exists(.mg_dir_new)) {
  invisible(file.copy(.mg_txt_files, file.path(.mg_dir_new, basename(.mg_txt_files)),
                      overwrite = TRUE))
}
if (length(.mg_txt_files) && exists(".mg_dir_new_pdf") && dir.exists(.mg_dir_new_pdf)) {
  invisible(file.copy(.mg_txt_files, file.path(.mg_dir_new_pdf, basename(.mg_txt_files)),
                      overwrite = TRUE))
}
ts_msg("Captions copied to dated folders (", length(.mg_txt_files), " txt files)")

# =============================================================================
# ==== SAFETY NET: define any res_* fallback symbols the extracted figure
#      code still references (only reached if the corresponding cap_csv/CLD
#      cache is missing -- R's lazy argument evaluation means these never
#      actually get forced when the cache is present, but defining them keeps
#      the module self-contained and independent of that subtlety). ==========
# =============================================================================
for (.rn in c("res_cspd_agg","res_hull_agg","res_hull_tp","res_iid_agg","res_iid_tp",
              "res_nnd_agg","res_nnd_tp",
              "res_zone_flow_logit_agg","res_zone_main","res_zone_main_cells_agg",
              "res_zone_sec_cells_agg","res_zone_sec_high_agg")) {
  if (!exists(.rn, inherits = FALSE)) assign(.rn, NULL)
}

ts_msg("=== figures_step5_module: figure regeneration complete ===")
ts_msg("Figure_8 (zone/ALR, 6 panels), Figure_9 (collective, 6 panels),")
ts_msg("Figure_10 (interval, 4 panels), figure7bis all written to: ", .mg_dir)
ts_msg("Revised interval set: Figure_INT1 (3 panels), Figure_INT2 (2 panels),")
ts_msg("Figure_INT3 (4 panels) — the three construct families.")
ts_msg("Mirrored to: ", .all_manu_dir)
ts_msg("NOTE: Figure S7 (Jacobs D) was NOT regenerated -- out of scope, see header.")
