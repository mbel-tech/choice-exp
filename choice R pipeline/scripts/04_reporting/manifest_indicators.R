# =============================================================================
# manifest_indicators.R
# -----------------------------------------------------------------------------
# Single source of truth: one row per presented indicator across Pipeline A
# (activity_analysis_STATS_choice_exp.R) and Pipeline B (analysis_b3.R).
#
# IN SCOPE
#   - Primary inferential models reported in the manuscript
#   - Descriptive demoted models (still inferentially fit, reported alongside)
#   - Pipeline B: cortisol + every populated mono_cell_models entry
#
# OUT OF SCOPE (intentionally excluded; do NOT add these rows)
#   - Jacobs atanh sensitivity models  (res_jac_<zone>_<path>)
#   - Jacobs joint models              (res_jac_beta_main_joint_agg,
#                                       res_jac_beta_sec_joint_agg)
#   - Monoamine omnibus §5 exploratory (mono_models)
#   - Any reference to BH correction (do not name BH families, do not validate
#     against bh_*; the table reports the FIT, not the post-hoc correction)
# =============================================================================

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
})

.PATH_A <- file.path(PROJECT_ROOT, "choice R pipeline/scripts/01_pipeline_analysis/activity_analysis_STATS_choice_exp.R")
.PATH_B <- file.path(PROJECT_ROOT, "choice exp for claude mono and cortisol/claude output/Scripts/analysis_b3_REVISED.R")
.PATH_C <- file.path(PROJECT_ROOT, "choice R pipeline/scripts/01_pipeline_analysis/bout_structure_analysis_choice_exp.R")
.PATH_D <- file.path(PROJECT_ROOT, "choice R pipeline/scripts/01_pipeline_analysis/behavioural_sequence_analysis_choice_exp.R")

# ---- Helper to build manifest rows tersely ----------------------------------
.row <- function(block, indicator, tier, object_expr, meta_expr,
                 source_path, source_line, response_pretty, unit) {
  tibble(
    block           = block,
    indicator       = indicator,
    tier            = tier,
    object_expr     = object_expr,
    meta_expr       = meta_expr,
    source_path     = source_path,
    source_line     = as.character(source_line),
    response_pretty = response_pretty,
    unit            = unit
  )
}

# ---- Pipeline A: PRIMARY (collective movement & switches) -------------------
.A_primary_movement <- bind_rows(
  .row("Collective movement", "School transitions (aggregated)",
       "primary", "res_switch_agg$model", "res_switch_agg",
       .PATH_A, 2717, "Switches per session", "Trial (aggregated)"),
  .row("Collective movement", "School transitions (interval)",
       "primary", "res_switch_tp$model", "res_switch_tp",
       .PATH_A, 2702, "Switches per session", "Session"),
  .row("Collective movement", "Nearest-neighbour distance (aggregated)",
       "primary", "res_nnd_agg$model", "res_nnd_agg",
       .PATH_A, 2874, "Mean NND (cm)", "Trial (aggregated)"),
  .row("Collective movement", "Nearest-neighbour distance (interval)",
       "primary", "res_nnd_tp$model", "res_nnd_tp",
       .PATH_A, 2859, "Mean NND (cm)", "Session"),
  # Polarisation removed 2026-08-19 ("drop dead/unreported outcomes"). It was
  # deleted from the manuscript on 2026-08-09 -- it is the only indicator that
  # requires persistent identity tracking, so an ID swap corrupts it -- and
  # res_pol_agg / res_pol_tp no longer exist in the STATS script. The matching
  # rows were removed from indicator_column_map.R the same day but NOT from
  # here, which left validate_indicator_map() aborting every workbook rebuild
  # with "Manifest indicators NOT in map (2)". Removed here 2026-08-29 to close
  # that. STEP2b still computes mean_polarisation if it is ever wanted again.
  .row("Collective movement", "School centroid speed (aggregated)",
       "primary", "res_cspd_agg$model", "res_cspd_agg",
       .PATH_A, 2974, "School speed (cm/s)", "Trial (aggregated)"),
  .row("Collective movement", "School centroid speed (interval)",
       "primary", "res_cspd_tp$model", "res_cspd_tp",
       .PATH_A, 2959, "School speed (cm/s)", "Session")
)

# ---- Pipeline A: PRIMARY (Jacobs β-GLMM per zone) ---------------------------
# 4 zones × 2 strata (agg + tp) = 8 rows. Jacobs joints + atanh sensitivity
# are EXCLUDED by design (see file header).
.A_primary_jacobs <- bind_rows(
  .row("Jacobs D — per zone", "Jacobs D — flow zone (aggregated)",
       "primary", "res_jac_beta_flow_agg$model", "res_jac_beta_flow_agg",
       .PATH_A, 2536, "Jacobs D (flow)", "Trial (aggregated)"),
  .row("Jacobs D — per zone", "Jacobs D — flow zone (interval)",
       "primary", "res_jac_beta_flow_tp$model", "res_jac_beta_flow_tp",
       .PATH_A, 2541, "Jacobs D (flow)", "Session"),
  .row("Jacobs D — per zone", "Jacobs D — high zone (aggregated)",
       "primary", "res_jac_beta_high_agg$model", "res_jac_beta_high_agg",
       .PATH_A, 2537, "Jacobs D (high)", "Trial (aggregated)"),
  .row("Jacobs D — per zone", "Jacobs D — high zone (interval)",
       "primary", "res_jac_beta_high_tp$model", "res_jac_beta_high_tp",
       .PATH_A, 2543, "Jacobs D (high)", "Session"),
  .row("Jacobs D — per zone", "Jacobs D — medium zone (aggregated)",
       "primary", "res_jac_beta_med_agg$model", "res_jac_beta_med_agg",
       .PATH_A, 2538, "Jacobs D (medium)", "Trial (aggregated)"),
  .row("Jacobs D — per zone", "Jacobs D — medium zone (interval)",
       "primary", "res_jac_beta_med_tp$model", "res_jac_beta_med_tp",
       .PATH_A, 2545, "Jacobs D (medium)", "Session"),
  .row("Jacobs D — per zone", "Jacobs D — low zone (aggregated)",
       "primary", "res_jac_beta_low_agg$model", "res_jac_beta_low_agg",
       .PATH_A, 2539, "Jacobs D (low)", "Trial (aggregated)"),
  .row("Jacobs D — per zone", "Jacobs D — low zone (interval)",
       "primary", "res_jac_beta_low_tp$model", "res_jac_beta_low_tp",
       .PATH_A, 2547, "Jacobs D (low)", "Session")
)

# ---- Pipeline A: DESCRIPTIVE (zone preference + cells + activity/flux/IID/hull) ----
.A_descriptive <- bind_rows(
  # Zone preference — main
  .row("Zone preference — main", "Main-zone preference (interval, alr)",
       "descriptive", "res_zone_main$model", "res_zone_main",
       .PATH_A, 1792, "ALR(flow vs calm)", "Session"),
  .row("Zone preference — main", "Main-zone preference (aggregated, β-GLMM)",
       "descriptive", "res_zone_main_agg$model", "res_zone_main_agg",
       .PATH_A, 1864, "Proportion in flow", "Trial (aggregated)"),
  .row("Zone preference — main", "Main-zone preference (aggregated, logit-LMM)",
       "descriptive", "res_zone_flow_logit_agg$model", "res_zone_flow_logit_agg",
       .PATH_A, 1880, "ALR(flow vs calm)", "Trial (aggregated)"),
  # Zone preference — sub (alr per stratum)
  .row("Zone preference — sub", "Sub-zone high (interval, alr)",
       "descriptive", "res_zone_sec_high_tp$model", "res_zone_sec_high_tp",
       .PATH_A, 1847, "ALR(high vs calm)", "Session"),
  .row("Zone preference — sub", "Sub-zone medium (interval, alr)",
       "descriptive", "res_zone_sec_med_tp$model", "res_zone_sec_med_tp",
       .PATH_A, 1848, "ALR(medium vs calm)", "Session"),
  .row("Zone preference — sub", "Sub-zone low (interval, alr)",
       "descriptive", "res_zone_sec_low_tp$model", "res_zone_sec_low_tp",
       .PATH_A, 1849, "ALR(low vs calm)", "Session"),
  .row("Zone preference — sub", "Sub-zone high (aggregated, alr)",
       "descriptive", "res_zone_sec_high_agg$model", "res_zone_sec_high_agg",
       .PATH_A, 1914, "ALR(high vs calm)", "Trial (aggregated)"),
  .row("Zone preference — sub", "Sub-zone medium (aggregated, alr)",
       "descriptive", "res_zone_sec_med_agg$model", "res_zone_sec_med_agg",
       .PATH_A, 1915, "ALR(medium vs calm)", "Trial (aggregated)"),
  .row("Zone preference — sub", "Sub-zone low (aggregated, alr)",
       "descriptive", "res_zone_sec_low_agg$model", "res_zone_sec_low_agg",
       .PATH_A, 1916, "ALR(low vs calm)", "Trial (aggregated)"),
  # Sub-zone alr pairwise comparisons (long format, all zones in one model)
  .row("Zone preference — sub", "Sub-zone alr pairs (aggregated, long)",
       "descriptive", "res_zone_sec_lr_pairs_agg$model", "res_zone_sec_lr_pairs_agg",
       .PATH_A, 1964, "ALR (long format)", "Trial (aggregated)"),
  .row("Zone preference — sub", "Sub-zone alr pairs (interval, long)",
       "descriptive", "res_zone_sec_lr_pairs_tp$model", "res_zone_sec_lr_pairs_tp",
       .PATH_A, 1985, "ALR (long format)", "Session"),
  # Cell means (zone × treatment)
  .row("Zone occupancy cells — main", "Main cells (aggregated, β-GLMM)",
       "descriptive", "res_zone_main_cells_agg$model", "res_zone_main_cells_agg",
       .PATH_A, 2097, "Cell proportion (main)", "Trial × zone"),
  .row("Zone occupancy cells — main", "Main cells (interval, β-GLMM)",
       "descriptive", "res_zone_main_cells_tp$model", "res_zone_main_cells_tp",
       .PATH_A, 2119, "Cell proportion (main)", "Session × zone"),
  .row("Zone occupancy cells — sub", "Sub cells CLR (aggregated)",
       "descriptive", "res_zone_sec_cells_agg$model", "res_zone_sec_cells_agg",
       .PATH_A, 2284, "CLR(sub-zone)", "Trial × zone"),
  .row("Zone occupancy cells — sub", "Sub cells CLR (interval)",
       "descriptive", "res_zone_sec_cells_tp$model", "res_zone_sec_cells_tp",
       .PATH_A, 2322, "CLR(sub-zone)", "Session × zone"),
  # Activity / flux / IID / hull (β-GLMM or LMM)
  .row("Activity & flux", "Proportion active (aggregated, β-GLMM)",
       "descriptive", "res_active_agg$model", "res_active_agg",
       .PATH_A, 2743, "Proportion active", "Trial (aggregated)"),
  .row("Activity & flux", "Proportion active (interval, β-GLMM)",
       "descriptive", "res_active_tp$model", "res_active_tp",
       .PATH_A, 2728, "Proportion active", "Session"),
  .row("Activity & flux", "Zone flux (aggregated)",
       "descriptive", "res_flux_agg$model", "res_flux_agg",
       .PATH_A, 2768, "Flux per session", "Trial (aggregated)"),
  .row("Activity & flux", "Zone flux (interval)",
       "descriptive", "res_flux_tp$model", "res_flux_tp",
       .PATH_A, 2753, "Flux per session", "Session"),
  .row("Collective movement (descriptive)", "Inter-individual distance (aggregated)",
       "descriptive", "res_iid_agg$model", "res_iid_agg",
       .PATH_A, 2924, "Mean IID (cm)", "Trial (aggregated)"),
  .row("Collective movement (descriptive)", "Inter-individual distance (interval)",
       "descriptive", "res_iid_tp$model", "res_iid_tp",
       .PATH_A, 2909, "Mean IID (cm)", "Session"),
  .row("Collective movement (descriptive)", "Convex hull area (aggregated)",
       "descriptive", "res_hull_agg$model", "res_hull_agg",
       .PATH_A, 2949, "School area (cm²)", "Trial (aggregated)"),
  .row("Collective movement (descriptive)", "Convex hull area (interval)",
       "descriptive", "res_hull_tp$model", "res_hull_tp",
       .PATH_A, 2934, "School area (cm²)", "Session")
)

# ---- Pipeline A: optional fish-level social indicators ----------------------
# These are conditional on data availability; if NULL, the orchestrator
# emits a row marked "Object NULL or does not exist" rather than failing.
.A_optional_social <- bind_rows(
  .row("Social (fish-level)", "Fish-fish NND (aggregated)",
       "descriptive", "res_fnnd_agg$model", "res_fnnd_agg",
       .PATH_A, 3160, "Mean fish-fish NND (px)", "Trial (aggregated)"),
  .row("Social (fish-level)", "Fish-fish NND (interval)",
       "descriptive", "res_fnnd_tp$model", "res_fnnd_tp",
       .PATH_A, 3145, "Mean fish-fish NND (px)", "Session"),
  .row("Social (fish-level)", "Centroid distance (aggregated)",
       "descriptive", "res_cdist_agg$model", "res_cdist_agg",
       .PATH_A, 3185, "Mean centroid distance (px)", "Trial (aggregated)"),
  .row("Social (fish-level)", "Centroid distance (interval)",
       "descriptive", "res_cdist_tp$model", "res_cdist_tp",
       .PATH_A, 3170, "Mean centroid distance (px)", "Session"),
  .row("Social (fish-level)", "Turning rate (aggregated)",
       "descriptive", "res_turn_agg$model", "res_turn_agg",
       .PATH_A, 3210, "Mean turning rate", "Trial (aggregated)"),
  .row("Social (fish-level)", "Turning rate (interval)",
       "descriptive", "res_turn_tp$model", "res_turn_tp",
       .PATH_A, 3195, "Mean turning rate", "Session")
)

# ---- Pipeline B: cortisol ---------------------------------------------------
.B_cortisol <- .row(
  "Endocrine — cortisol", "Plasma cortisol",
  "primary", "fit_cort", NA_character_,
  .PATH_B, 747, "Plasma cortisol (ng/mL)", "Sample"
)

# ---- Pipeline B: monoamine per-cell (7 analytes × 4 areas, up to 28) --------
# cell_id format: paste(analyte_key, area, sep="__")
# NE removed 2026-08-29: noradrenaline is quantified but is not entered into any
# reported model, so it is not a presented indicator. The analysed grid is the
# 24 analyte x region cells below, not 28; monoamine_long.csv (the analysis
# source) already excludes NE.
.B_analyte_keys <- c("ht_5","hiaa_5","hiaa_5_ratio","da","dopac","dopac_da_ratio")
.B_analyte_label <- c(ht_5 = "5-HT", hiaa_5 = "5-HIAA",
                     hiaa_5_ratio = "5-HIAA/5-HT",
                     da = "DA", dopac = "DOPAC",
                     dopac_da_ratio = "DOPAC/DA")
.B_areas <- c("DM","POA","VV","VD")

.B_per_cell <- do.call(bind_rows, lapply(.B_analyte_keys, function(ak) {
  do.call(bind_rows, lapply(.B_areas, function(ar) {
    cid <- paste(ak, ar, sep = "__")
    pretty <- sprintf("%s — %s", .B_analyte_label[[ak]], ar)
    .row(
      block        = sprintf("Endocrine — %s", .B_analyte_label[[ak]]),
      indicator    = pretty,
      tier         = "primary",
      object_expr  = sprintf('mono_cell_models[["%s"]]$fit', cid),
      meta_expr    = sprintf('mono_cell_models[["%s"]]', cid),
      source_path  = .PATH_B,
      source_line  = "1145–1182",
      response_pretty = sprintf("%s (per-cell)", .B_analyte_label[[ak]]),
      unit         = "Sample"
    )
  }))
}))

# ---- Pipeline C: bout structure (BOUT engine) -------------------------------
# Reported outcomes 6-8 of the fourteen in
# scripts/00_shared/export_tukey_interval_module.R. Fit at trial x interval
# (n = 48) as val ~ treatment * timepoint_f [+ fish_density] + (1|trial);
# results in output/BOUT_output/<run>/bout_anova_trial_x_timepoint.csv.
# Added 2026-08-29: these were reported in the manuscript but had no manifest
# row, so they never reached the distributed workbook.
.C_bouts <- bind_rows(
  .row("Engagement pattern", "Flow bouts per minute (interval)",
       "primary", 'bout_anova_trial_x_timepoint[metric == "bout_rate_flow"]',
       NA_character_, .PATH_C, 942, "Flow bouts per minute", "Session"),
  .row("Engagement pattern", "Longest flow bout (interval)",
       "primary", 'bout_anova_trial_x_timepoint[metric == "max_flow_bout_s"]',
       NA_character_, .PATH_C, 942, "Longest flow bout (s)", "Session"),
  .row("Engagement pattern", "Longest calm bout (interval)",
       "primary", 'bout_anova_trial_x_timepoint[metric == "max_calm_bout_s"]',
       NA_character_, .PATH_C, 942, "Longest calm bout (s)", "Session")
)

# ---- Pipeline D: state sequence (SEQ engine, binary alphabet) ---------------
# Reported outcomes 9-10. Dropped from Figure B on 2026-08-17 but retained in
# Table 1 and the Results text: they are the MEAN of the same bout-duration
# distribution whose MAXIMUM the BOUT engine reports.
.D_seq <- bind_rows(
  .row("Engagement pattern", "Mean flow bout duration (interval)",
       "primary", 'seq_anova_trial_x_timepoint[metric == "dwell_Flow"]',
       NA_character_, .PATH_D, 545, "Mean flow bout duration (s)", "Session"),
  .row("Engagement pattern", "Mean calm bout duration (interval)",
       "primary", 'seq_anova_trial_x_timepoint[metric == "dwell_Calm"]',
       NA_character_, .PATH_D, 545, "Mean calm bout duration (s)", "Session")
)

# ---- Bind everything --------------------------------------------------------
MANIFEST <- bind_rows(
  .A_primary_movement,
  .A_primary_jacobs,
  .A_descriptive,
  .A_optional_social,
  .C_bouts,
  .D_seq,
  .B_cortisol,
  .B_per_cell
)

message(sprintf("[manifest] %d rows: %d primary, %d descriptive",
                nrow(MANIFEST),
                sum(MANIFEST$tier == "primary"),
                sum(MANIFEST$tier == "descriptive")))
