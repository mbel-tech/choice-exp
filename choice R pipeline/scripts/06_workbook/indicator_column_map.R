# =============================================================================
# indicator_column_map.R — manifest → CSV column → workbook column
# =============================================================================
# Contract between manifest_indicators.R (74 rows) and the consolidated workbook.
#
# Every primary+descriptive manifest row maps to exactly one workbook column.
#
# Source columns reference:
#   - behaviour: D:/CHOICE R SCRIPTS/choice R pipeline/easy_scripts/easy_scripts_dataset.csv
#   - endo:      D:/CHOICE R SCRIPTS/choice R pipeline/easy_scripts/easy_scripts_endo_dataset.csv
#
# Derivation modes:
#   "direct"          — CSV column is taken as-is (already at correct grain)
#   "interval"        — CSV column is per-interval; copied directly
#   "broadcast_mean"  — aggregate the per-interval CSV column by phys_trial mean
#   "representative"  — CSV column is a representative cell of a multi-cell model;
#                       flagged in legend with note explaining the choice
#   "missing"         — no CSV column available; emit NA column + legend note
#   "endo_cell"       — endocrine cell, lives in endo_wide (not behaviour_wide)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

# ---- Behaviour map ----------------------------------------------------------
# 48 behavioural rows from manifest_indicators.R
INDICATOR_MAP_BEH <- tibble::tribble(
  ~manifest_label,                                       ~es_col,                  ~wb_col,                          ~derive,           ~unit,
  # ---- Collective movement (primary, 8) ----
  "School transitions (interval)",                       "switches_per_session",   "switches_per_session",           "interval",        "count / session",
  "School transitions (aggregated)",                     "switches_per_session",   "switches_per_session_agg",       "broadcast_mean",  "count / session",
  "Nearest-neighbour distance (interval)",               "mean_nnd_cm",            "mean_nnd_cm",                    "interval",        "cm",
  "Nearest-neighbour distance (aggregated)",             "mean_nnd_cm",            "mean_nnd_cm_agg",                "broadcast_mean",  "cm",
  "School centroid speed (interval)",                    "mean_centroid_spd_cm",   "mean_centroid_spd_cm",           "interval",        "cm / s",
  "School centroid speed (aggregated)",                  "mean_centroid_spd_cm",   "mean_centroid_spd_cm_agg",       "broadcast_mean",  "cm / s",

  # ---- Jacobs D per zone (primary, 8) ----
  "Jacobs D — flow zone (interval)",                "D_main_flow",            "jacobs_d_flow",                  "interval",        "Jacobs D (-1..1)",
  "Jacobs D — flow zone (aggregated)",              "D_main_flow_agg",        "jacobs_d_flow_agg",              "direct",          "Jacobs D (-1..1)",
  "Jacobs D — high zone (interval)",                "D_sec_high",             "jacobs_d_high",                  "interval",        "Jacobs D (-1..1)",
  "Jacobs D — high zone (aggregated)",              "D_sec_high_agg",         "jacobs_d_high_agg",              "direct",          "Jacobs D (-1..1)",
  "Jacobs D — medium zone (interval)",              "D_sec_medium",           "jacobs_d_medium",                "interval",        "Jacobs D (-1..1)",
  "Jacobs D — medium zone (aggregated)",            "D_sec_medium_agg",       "jacobs_d_medium_agg",            "direct",          "Jacobs D (-1..1)",
  "Jacobs D — low zone (interval)",                 "D_sec_low",              "jacobs_d_low",                   "interval",        "Jacobs D (-1..1)",
  "Jacobs D — low zone (aggregated)",               "D_sec_low_agg",          "jacobs_d_low_agg",               "direct",          "Jacobs D (-1..1)",

  # ---- Zone preference — main (descriptive, 3) ----
  "Main-zone preference (interval, alr)",                "logit_flow",             "logit_flow",                     "interval",        "logit",
  "Main-zone preference (aggregated, β-GLMM)",      "prop_flow",              "prop_flow_agg",                  "broadcast_mean",  "proportion (0-1)",
  "Main-zone preference (aggregated, logit-LMM)",        "logit_flow",             "logit_flow_agg",                 "broadcast_mean",  "logit",

  # ---- Zone preference — sub (descriptive, 8) ----
  "Sub-zone high (interval, alr)",                       "lr_high",                "lr_high",                        "interval",        "log-ratio",
  "Sub-zone medium (interval, alr)",                     "lr_medium",              "lr_medium",                      "interval",        "log-ratio",
  "Sub-zone low (interval, alr)",                        "lr_low",                 "lr_low",                         "interval",        "log-ratio",
  "Sub-zone high (aggregated, alr)",                     "lr_high",                "lr_high_agg",                    "broadcast_mean",  "log-ratio",
  "Sub-zone medium (aggregated, alr)",                   "lr_medium",              "lr_medium_agg",                  "broadcast_mean",  "log-ratio",
  "Sub-zone low (aggregated, alr)",                      "lr_low",                 "lr_low_agg",                     "broadcast_mean",  "log-ratio",
  "Sub-zone alr pairs (aggregated, long)",               NA,                       "subzone_alr_pairs_agg",          "missing",         NA,
  "Sub-zone alr pairs (interval, long)",                 NA,                       "subzone_alr_pairs_tp",           "missing",         NA,

  # ---- Zone occupancy cells (descriptive, 4) ----
  # These are multi-cell models — picking a representative "headline" cell:
  "Main cells (aggregated, β-GLMM)",                "p_main_flow_agg",        "main_cells_flow_agg",            "representative",  "proportion (0-1)",
  "Main cells (interval, β-GLMM)",                  "p_main_flow_tp",         "main_cells_flow_tp",             "representative",  "proportion (0-1)",
  "Sub cells CLR (aggregated)",                          "clr_sec_high_agg",       "sub_cells_high_clr_agg",         "representative",  "CLR",
  "Sub cells CLR (interval)",                            "clr_sec_high_tp",        "sub_cells_high_clr_tp",          "representative",  "CLR",

  # ---- Activity & flux (descriptive, 4) ----
  "Proportion active (aggregated, β-GLMM)",         "prop_active",            "prop_active_agg",                "broadcast_mean",  "proportion (0-1)",
  "Proportion active (interval, β-GLMM)",           "prop_active",            "prop_active",                    "interval",        "proportion (0-1)",
  "Zone flux (aggregated)",                              "zone_flux_per_session",  "zone_flux_per_session_agg",      "broadcast_mean",  "count / session",
  "Zone flux (interval)",                                "zone_flux_per_session",  "zone_flux_per_session",          "interval",        "count / session",

  # ---- Collective movement — descriptive (4) ----
  "Inter-individual distance (aggregated)",              "mean_iid_cm",            "mean_iid_cm_agg",                "broadcast_mean",  "cm",
  "Inter-individual distance (interval)",                "mean_iid_cm",            "mean_iid_cm",                    "interval",        "cm",
  "Convex hull area (aggregated)",                       "mean_hull_area_cm2",     "mean_hull_area_cm2_agg",         "broadcast_mean",  "cm²",
  "Convex hull area (interval)",                         "mean_hull_area_cm2",     "mean_hull_area_cm2",             "interval",        "cm²",

  # ---- Social (fish-level) (descriptive, 6) — source data not in current run
  "Fish-fish NND (aggregated)",                          NA,                       "fish_fish_nnd_agg",              "missing",         "px",
  "Fish-fish NND (interval)",                            NA,                       "fish_fish_nnd",                  "missing",         "px",
  "Centroid distance (aggregated)",                      NA,                       "centroid_distance_agg",          "missing",         "px",
  "Centroid distance (interval)",                        NA,                       "centroid_distance",              "missing",         "px",
  "Turning rate (aggregated)",                           NA,                       "turning_rate_agg",               "missing",         "deg / s",
  "Turning rate (interval)",                             NA,                       "turning_rate",                   "missing",         "deg / s",

  # ---- Engagement pattern: bout structure & state sequence (primary, 5) ----
  # Added 2026-08-29. Reported outcomes 6-10 of the fourteen in
  # export_tukey_interval_module.R; they had no map row, so the distributed
  # workbook shipped without five of the outcomes the manuscript reports.
  # es_col names come from BOUT bout_metrics_per_session.csv (bout_*, max_*)
  # and SEQ seq_metrics_per_trial_binary.csv (dwell_*), which the builder
  # joins onto the behaviour frame on (tank, fish_density, interval) -- both
  # engines spell the interval `timepoint` in their own headers.
  "Flow bouts per minute (interval)",                    "bout_rate_flow",         "bouts_per_min",                  "interval",        "count / min",
  "Longest flow bout (interval)",                        "max_flow_bout_s",        "longest_flow_bout_s",            "interval",        "s",
  "Longest calm bout (interval)",                        "max_calm_bout_s",        "longest_calm_bout_s",            "interval",        "s",
  "Mean flow bout duration (interval)",                  "dwell_Flow",             "mean_flow_bout_s",               "interval",        "s",
  "Mean calm bout duration (interval)",                  "dwell_Calm",             "mean_calm_bout_s",               "interval",        "s"
)

# ---- Endocrine map ----------------------------------------------------------
# 25 endocrine rows from manifest_indicators.R (1 cortisol + 6 analytes x 4 areas)
# NE dropped 2026-08-29. Noradrenaline was quantified but is not entered into
# any reported model -- the analysed grid is 24 analyte x region cells, not 28
# -- and monoamine_long.csv (the analysis source) already excludes it.
.endo_analytes <- c("ht_5", "hiaa_5", "hiaa_5_ratio", "da", "dopac", "dopac_da_ratio")
.endo_areas    <- c("DM", "POA", "VV", "VD")
.endo_pretty <- c(
  ht_5            = "5-HT",
  hiaa_5          = "5-HIAA",
  hiaa_5_ratio    = "5-HIAA/5-HT",
  da              = "DA",
  dopac           = "DOPAC",
  dopac_da_ratio  = "DOPAC/DA"

)
.endo_units <- c(
  ht_5            = "ng/mg tissue",
  hiaa_5          = "ng/mg tissue",
  hiaa_5_ratio    = "ratio (unitless)",
  da              = "ng/mg tissue",
  dopac           = "ng/mg tissue",
  dopac_da_ratio  = "ratio (unitless)"

)

.endo_rows <- do.call(rbind, lapply(.endo_analytes, function(a) {
  do.call(rbind, lapply(.endo_areas, function(ar) {
    data.frame(
      manifest_label = sprintf("%s — %s", .endo_pretty[a], ar),
      endo_cell      = sprintf("%s__%s", a, ar),
      wb_col         = sprintf("%s__%s", a, ar),
      derive         = "endo_cell",
      unit           = unname(.endo_units[a]),
      stringsAsFactors = FALSE
    )
  }))
}))

INDICATOR_MAP_ENDO <- tibble::tibble(
  manifest_label = c("Plasma cortisol", .endo_rows$manifest_label),
  endo_cell      = c("cort",                .endo_rows$endo_cell),
  wb_col         = c("cort_plasma",         .endo_rows$wb_col),
  derive         = c("endo_cell",           .endo_rows$derive),
  unit           = c("ng/mL plasma",       .endo_rows$unit)
)

# ---- Validator --------------------------------------------------------------
validate_indicator_map <- function(manifest, beh_map = INDICATOR_MAP_BEH,
                                    endo_map = INDICATOR_MAP_ENDO,
                                    verbose = TRUE) {
  stopifnot(is.data.frame(manifest), "indicator" %in% names(manifest),
            "tier" %in% names(manifest))

  mf_in <- manifest[manifest$tier %in% c("primary", "descriptive"), ]
  all_labels  <- mf_in$indicator
  mapped_beh  <- beh_map$manifest_label
  mapped_endo <- endo_map$manifest_label
  mapped_all  <- c(mapped_beh, mapped_endo)

  miss_in_map <- setdiff(all_labels,   mapped_all)
  orphan_map  <- setdiff(mapped_all,   all_labels)

  if (length(miss_in_map) > 0) {
    stop("[indicator_map] Manifest indicators NOT in map (", length(miss_in_map), "):\n  - ",
         paste(miss_in_map, collapse = "\n  - "), call. = FALSE)
  }
  if (length(orphan_map) > 0) {
    stop("[indicator_map] Map entries NOT in manifest (", length(orphan_map), "):\n  - ",
         paste(orphan_map, collapse = "\n  - "), call. = FALSE)
  }

  # Duplicate workbook columns?
  dup_beh  <- beh_map$wb_col[duplicated(beh_map$wb_col)]
  dup_endo <- endo_map$wb_col[duplicated(endo_map$wb_col)]
  if (length(dup_beh)  > 0) stop("[indicator_map] Duplicate beh wb_col: ",
                                 paste(dup_beh, collapse = ", "), call. = FALSE)
  if (length(dup_endo) > 0) stop("[indicator_map] Duplicate endo wb_col: ",
                                 paste(dup_endo, collapse = ", "), call. = FALSE)

  if (verbose) {
    cat(sprintf("[indicator_map] OK: %d manifest indicators -> %d beh + %d endo workbook columns\n",
                length(all_labels), nrow(beh_map), nrow(endo_map)))
    .miss <- sum(beh_map$derive == "missing")
    .rep  <- sum(beh_map$derive == "representative")
    .bro  <- sum(beh_map$derive == "broadcast_mean")
    cat(sprintf("  - broadcast_mean: %d  |  representative: %d  |  missing (NA cols): %d\n",
                .bro, .rep, .miss))
  }
  invisible(TRUE)
}
