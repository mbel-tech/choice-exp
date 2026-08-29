###### NORMALIZATION ##########################################################################

# ---------------------------
# FUNCTION 1: Fine-Grained Normalization
# ---------------------------
# NOTE (2026-08-18): `df` must carry RAW occupancy in high/medium/low/calm.
# Passing STEP2's prop_time_ac_* values (already area-corrected and rescaled)
# applies the area correction twice and distorts the composition -- that was
# the bug fixed in activity_analysis_GRAPHS_choice_exp.R on 2026-08-18.
normalize_fine_zones <- function(df) {
  df %>%
    dplyr::mutate(
      high_norm   = high / 10,
      medium_norm = medium / 18,
      low_norm    = low / 13,
      calm_norm   = calm / 42
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      total_fine = sum(c_across(c(high_norm, medium_norm, low_norm, calm_norm)), na.rm = TRUE),
      high_pct   = high_norm / total_fine * 100,
      medium_pct = medium_norm / total_fine * 100,
      low_pct    = low_norm / total_fine * 100,
      calm_pct   = calm_norm / total_fine * 100
    ) %>%
    dplyr::ungroup()
}

#####################################

color_map <- c(flow="red4", calm="lightskyblue3", high="tomato3",
               medium="mediumturquoise", low="goldenrod3")
color_map_broad <- c(flow="red4", calm="lightskyblue3")

##################################


# Get misc / okabe palette (older versions return a palette object)
pal_full <- cols4all::c4a("okabe")

# Convert palette to actual hex colors
pal_vec <- as.character(pal_full)

# We want: color 5 = Exercise choice, color 6 = Control
pal <- c(
  "Control"         = pal_vec[6],
  "Exercise choice" = pal_vec[5]
)

