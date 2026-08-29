# =============================================================================
# STEP 4 — CHOICE EXPERIMENT
# Build scatter plots and zone occupancy plots; save as editable PowerPoint.
#
# TWO PARALLEL OUTPUT SETS:
#
#   A) TIMEPOINT SET  (one row per fish × timepoint; faceted by timepoint)
#      STEP4_OUT/timepoint/
#        standard_indicator_plots_tp.pptx  — one faceted slide per indicator
#        zone_plots_tp.pptx                — slide 1: main zones; slide 2: sec zones
#        png_preview/
#      Caption: "n (fish \u00d7 timepoint) = X"
#
#   B) AGGREGATED SET  (means per fish across timepoints; plain scatter)
#      STEP4_OUT/aggregated/
#        standard_indicator_plots_agg.pptx
#        zone_plots_agg.pptx
#        png_preview/
#      Caption: "n (fish) = X  (mean across timepoints)"
#
# PLOT TYPES:
#   Standard indicator:  x = treatment; colour = treatment; shape = treatment
#   Main zone:           x = zone (flow, calm); colour = zone; shape = treatment;
#                        treatments dodged within zone
#   Secondary zone:      x = zone (high, medium, low, calm); colour = zone;
#                        shape = treatment; treatments dodged within zone
#
# INPUT:
#   fish_activity_summary   (one row per fish × timepoint, from GlobalEnv or disk)
#   group_dynamics_summary  (one row per trial × timepoint, from STEP2b; optional)
#
# OUTPUTS:
#   STEP4_graphs/STEP4_graphs_<timestamp>/  (two subdirectories as above)
#     timepoint/
#       standard_indicator_plots_tp.pptx
#       zone_plots_tp.pptx
#       group_dynamics_plots_tp.pptx     (if group_dynamics_summary available)
#       social_context_plots_tp.pptx     (if social columns in fish_activity_summary)
#     aggregated/
#       standard_indicator_plots_agg.pptx
#       zone_plots_agg.pptx
#       group_dynamics_plots_agg.pptx
#       social_context_plots_agg.pptx
# =============================================================================


# =============================================================================
# ==== 0) SETUP ===============================================================
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(officer)
  library(rvg)
  library(readr)
  library(purrr)
  library(stringr)
  library(cols4all)
  library(grid)        # rasterGrob()
  if (requireNamespace("patchwork", quietly = TRUE)) library(patchwork)
  if (requireNamespace("magick",    quietly = TRUE)) library(magick)
})

ts_msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), " | ", ..., "\n", sep = "")

.get_global <- function(name, default = NULL) {
  if (exists(name, envir = .GlobalEnv, inherits = FALSE))
    get(name, envir = .GlobalEnv, inherits = FALSE)
  else default
}

.step4_dir <- tryCatch({
  frames     <- sys.frames()
  ofile_envs <- Filter(function(f) exists("ofile", envir = f, inherits = FALSE), frames)
  if (length(ofile_envs) > 0)
    normalizePath(dirname(get("ofile", envir = ofile_envs[[length(ofile_envs)]])),
                  winslash = "/", mustWork = FALSE)
  else getwd()
}, error = function(e) getwd())

.find_latest_csv <- function(step_name, csv_filename) {
  parent <- file.path(getwd(), step_name)
  if (!dir.exists(parent)) return(NULL)
  subdirs <- list.dirs(parent, full.names = TRUE, recursive = FALSE)
  subdirs <- subdirs[grepl(paste0("^", step_name, "_\\d{8}_\\d{6}$"), basename(subdirs))]
  if (length(subdirs) == 0) return(NULL)
  latest    <- subdirs[which.max(file.mtime(subdirs))]
  candidate <- file.path(latest, csv_filename)
  if (file.exists(candidate)) candidate else NULL
}

TREATMENT_LEVELS_g <- .get_global("TREATMENT_LEVELS", c("control", "exercise choice"))
DENSITY_LEVELS_g   <- .get_global("DENSITY_LEVELS",   c(4L, 8L, 12L, 16L))
TIMEPOINT_LEVELS_g <- .get_global("TIMEPOINT_LEVELS", 1:3)

.step4_parent <- file.path(getwd(), "STEP4_graphs")
if (!dir.exists(.step4_parent))
  dir.create(.step4_parent, recursive = TRUE, showWarnings = FALSE)
STEP4_OUT <- file.path(.step4_parent,
                        paste0("STEP4_graphs_", format(Sys.time(), "%Y%m%d_%H%M%S")))
dir.create(file.path(STEP4_OUT, "timepoint",   "png_preview"), recursive = TRUE)
dir.create(file.path(STEP4_OUT, "aggregated",  "png_preview"), recursive = TRUE)
ts_msg("Output directory: ", STEP4_OUT)


# =============================================================================
# ==== 1) LOAD DATA ===========================================================
# =============================================================================

if (exists("fish_activity_summary", envir = .GlobalEnv, inherits = FALSE)) {
  ts_msg("Using fish_activity_summary from GlobalEnv")
  df_sum <- get("fish_activity_summary", envir = .GlobalEnv, inherits = FALSE)
} else {
  ts_msg("fish_activity_summary not in GlobalEnv; searching disk...")
  # STEP2 writes trial_activity_summary.csv; also check legacy name
  .csv_path <- .find_latest_csv("STEP2_output", "trial_activity_summary.csv")
  if (is.null(.csv_path))
    .csv_path <- .find_latest_csv("STEP2_output", "fish_activity_summary.csv")
  if (!is.null(.csv_path)) {
    df_sum <- readr::read_csv(.csv_path, show_col_types = FALSE)
    assign("fish_activity_summary", df_sum, envir = .GlobalEnv)
    ts_msg("Loaded from: ", .csv_path)
  } else {
    stop("Cannot find trial_activity_summary / fish_activity_summary. Run Steps 1-2 first.", call. = FALSE)
  }
}

# Check timepoint
if (!"timepoint" %in% names(df_sum)) {
  warning("'timepoint' column not found — timepoint plot set will be skipped.", call. = FALSE)
  df_sum$timepoint <- NA_integer_
}

# Add identifiers and factors
df_sum <- df_sum %>%
  dplyr::mutate(
    treatment    = factor(trimws(tolower(as.character(treatment))),
                          levels = TREATMENT_LEVELS_g),
    fish_density_f = factor(fish_density, levels = DENSITY_LEVELS_g, ordered = TRUE),
    timepoint    = as.integer(timepoint),
    timepoint_f  = factor(timepoint, levels = TIMEPOINT_LEVELS_g),
    trial_id     = as.integer(trial_id),
    tank         = as.character(tank),
    # One observation per trial × timepoint (school-level data).
    # fish_uid is trial_id as character — used as the scatter-plot point ID.
    fish_uid     = as.character(trial_id)
  )

.n_tanks   <- dplyr::n_distinct(df_sum$tank)
.n_trials  <- dplyr::n_distinct(df_sum$trial_id)
ts_msg("Graph data: ", nrow(df_sum), " rows | ",
       .n_tanks, " tanks | ", .n_trials, " trials")
readr::write_csv(df_sum, file.path(STEP4_OUT, "fish_activity_summary_for_graphs.csv"))

# The aggregated (non-faceted) plot set uses df_sum directly — all observations
# from all timepoints shown as individual points.  No per-fish averaging is
# possible because fish identities are not preserved across segments.

# ---- Load group_dynamics_summary (one row per trial x timepoint) ----
if (exists("group_dynamics_summary", envir = .GlobalEnv, inherits = FALSE)) {
  ts_msg("Using group_dynamics_summary from GlobalEnv")
  df_gd <- get("group_dynamics_summary", envir = .GlobalEnv, inherits = FALSE)
} else {
  ts_msg("group_dynamics_summary not in GlobalEnv; searching disk...")
  .gd_csv <- .find_latest_csv("STEP2b_output", "group_dynamics_summary.csv")
  if (!is.null(.gd_csv)) {
    df_gd <- readr::read_csv(.gd_csv, show_col_types = FALSE)
    assign("group_dynamics_summary", df_gd, envir = .GlobalEnv)
    ts_msg("Loaded from: ", .gd_csv)
  } else {
    ts_msg("group_dynamics_summary not found — group dynamics plots will be skipped.")
    df_gd <- NULL
  }
}

if (!is.null(df_gd) && nrow(df_gd) > 0) {
  # Attach metadata columns that df_gd may be missing.
  # Only join columns not already present to avoid .x/.y suffix conflicts.
  .gd_want <- c("trial_id","tank","treatment","fish_density","fish_density_f",
                "timepoint","timepoint_f")
  .gd_missing <- setdiff(.gd_want, names(df_gd))
  if (length(.gd_missing) > 0) {
    .gd_meta <- dplyr::distinct(
      df_sum[, intersect(c("trial_id", "timepoint", .gd_missing), names(df_sum))],
      trial_id, timepoint, .keep_all = TRUE
    )
    df_gd <- dplyr::left_join(df_gd, .gd_meta, by = c("trial_id", "timepoint"))
  }
  df_gd <- df_gd %>%
    dplyr::mutate(
      treatment   = factor(trimws(tolower(as.character(treatment))),
                           levels = TREATMENT_LEVELS_g),
      timepoint_f = factor(timepoint, levels = TIMEPOINT_LEVELS_g)
    )
  ts_msg("Group dynamics graph data: ", nrow(df_gd), " rows")
}


# =============================================================================
# ==== 3) PLOT SETTINGS =======================================================
# =============================================================================

PT_ALPHA       <- 0.85
PT_SIZE        <- 2.8    # standard indicator plots
ZONE_PT_SIZE   <- 4      # zone plots — 1/3 reduction from reference size 6
JITTER_WIDTH   <- 0.12
DODGE_WIDTH    <- 0.70
MEAN_WIDTH     <- 0.45
MEAN_LINE_SIZE <- 0.85
ERRORBAR_WIDTH <- 0.14
ERRORBAR_SIZE  <- 0.65

# Treatment: Okabe-Ito palette (colorblind-friendly; color 5 = exercise choice, 6 = control)
pal_vec          <- as.character(cols4all::c4a("okabe"))
TREATMENT_COLORS <- c("control" = pal_vec[6], "exercise choice" = pal_vec[5])
TREATMENT_SHAPES <- c("control" = 17L,        "exercise choice" = 16L)

# Zone colour maps — broad (main zones) and full (secondary zones)
# Matches reference: flow="red4", calm="lightskyblue3", high="tomato3",
#                    medium="mediumturquoise", low="goldenrod3"
color_map_broad <- c(flow = "red4", calm = "lightskyblue3")
color_map       <- c(flow = "red4", calm = "lightskyblue3",
                     high = "tomato3", medium = "mediumturquoise", low = "goldenrod3")

MAIN_ZONE_COLORS <- color_map_broad
SEC_ZONE_COLORS  <- color_map[c("high", "medium", "low", "calm")]

# BASE_THEME: used for standard indicator / group-dynamics scatter plots.
BASE_THEME <- ggplot2::theme_minimal(base_size = 13) +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    axis.line.x      = ggplot2::element_line(color = "black", linewidth = 0.85),
    axis.line.y      = ggplot2::element_line(color = "black", linewidth = 0.85),
    axis.ticks       = ggplot2::element_line(color = "black", linewidth = 0.7),
    axis.text.x      = ggplot2::element_text(size = 13, face = "bold", color = "black"),
    axis.text.y      = ggplot2::element_text(size = 11, color = "black"),
    axis.title       = ggplot2::element_text(size = 14, face = "bold", color = "black"),
    legend.title     = ggplot2::element_text(size = 12, face = "bold"),
    legend.text      = ggplot2::element_text(size = 11),
    plot.margin      = ggplot2::margin(t = 8, r = 10, b = 16, l = 10),
    strip.text       = ggplot2::element_text(size = 12, face = "bold", color = "black"),
    strip.background = ggplot2::element_blank()
  )

# ZONE_THEME: used for all zone occupancy plots (main and secondary, agg and TP).
# Aesthetic parameters matched to reference script "choice behavioral def june 2025":
#   axis.text.x  size=20, bold  (was 13 in BASE_THEME)
#   axis.text.y  size=20        (was 11)
#   axis.title.x element_blank() (reference hides x-axis title for zone plots)
#   axis.title.y size=26, bold  (was 14 combined)
#   plot.caption size=22, lineheight=1.15, margin(t=24)
#   plot.caption.position = "plot"
#   plot.margin  margin(25,25,25,25)  (was margin(8,10,16,10))
ZONE_THEME <- ggplot2::theme_minimal(base_size = 13) +
  ggplot2::theme(
    panel.grid       = ggplot2::element_blank(),
    axis.line.x      = ggplot2::element_line(),
    axis.line.y      = ggplot2::element_line(),
    axis.title.x     = ggplot2::element_blank(),
    axis.text.x      = ggplot2::element_text(size = 20, face = "bold"),
    axis.text.y      = ggplot2::element_text(size = 20),
    axis.title.y     = ggplot2::element_text(size = 26, face = "bold"),
    legend.title     = ggplot2::element_text(size = 12, face = "bold"),
    legend.text      = ggplot2::element_text(size = 11),
    plot.caption.position = "plot",
    plot.caption     = ggplot2::element_text(size = 22, lineheight = 1.15,
                                             margin = ggplot2::margin(t = 24)),
    plot.margin      = ggplot2::margin(25, 25, 25, 25),
    strip.text       = ggplot2::element_text(size = 14, face = "bold"),
    strip.background = ggplot2::element_blank()
  )

# =============================================================================
# Picture4 logo — upper-left inset for secondary zone occupancy plots.
# Placement matches reference: patchwork::inset_element(), align_to = "full",
# rel_width = 0.35, top-left corner (left=0, right=rel_width, top=1).
# File: D:\CHOICE R SCRIPTS\Picture4.jpg
# =============================================================================
.LOGO_PATH <- file.path(PROJECT_ROOT, "Picture4.jpg")
.logo_grob <- NULL   # populated below if file exists

.HAS_MAGICK    <- requireNamespace("magick",    quietly = TRUE)
.HAS_PATCHWORK <- requireNamespace("patchwork", quietly = TRUE)
.logo_aspect   <- 1   # fallback if image not loaded

if (.HAS_MAGICK && file.exists(.LOGO_PATH)) {
  tryCatch({
    .logo_img    <- magick::image_read(.LOGO_PATH)
    .logo_info   <- magick::image_info(.logo_img)
    .logo_aspect <- .logo_info$height / .logo_info$width   # h/w
    .logo_grob   <- grid::rasterGrob(as.raster(.logo_img), interpolate = TRUE)
    ts_msg("Logo loaded: ", .LOGO_PATH, " (", .logo_info$width, "x", .logo_info$height, ")")
  }, error = function(e) {
    warning("Logo load failed: ", conditionMessage(e)); .logo_grob <- NULL
  })
} else {
  if (!.HAS_MAGICK)
    ts_msg("NOTE: 'magick' package not installed — Picture4 logo will be omitted. ",
           "Install with: install.packages('magick')")
  else
    warning("Picture4.jpg not found at: ", .LOGO_PATH, " — logo will be omitted.")
}

# Helper: inset logo in the upper-RIGHT of a ggplot using patchwork::inset_element.
# rel_width: fraction of full plot width the image occupies (default 0.35).
# Requires both 'patchwork' and 'magick'; returns p unchanged if either is absent.
add_logo_top_right <- function(p, rel_width = 0.35) {
  if (is.null(.logo_grob) || !.HAS_PATCHWORK) return(p)
  rel_height <- min(rel_width * .logo_aspect, 0.30)   # clamp to avoid covering plot body
  p + patchwork::inset_element(
    .logo_grob,
    left   = 1 - rel_width,
    bottom = 1 - rel_height,
    right  = 1,
    top    = 1,
    align_to = "full"
  )
}


# =============================================================================
# ==== 4) INDICATOR LABELS AND TIMEPOINT LABELS ================================
# =============================================================================

# Pretty y-axis labels using Unicode superscripts
INDICATOR_PRETTY_LABELS <- c(
  # School-level identity-free indicators (current pipeline)
  prop_active            = "Proportion active (> 1 BL s\u207b\u00b9)",
  switches_per_session   = "Flow\u2194calm switches per 20-min session",
  zone_flux_per_session  = "Zone flux per session (|\u0394 prop\u209a\u1d40\u02b3| \u00d7 session; diagnostic)",
  prop_time_in_flow      = "Proportion of time in Flow zone (s)",
  prop_time_in_calm      = "Proportion of time in Calm zone (s)",
  prop_time_ac_high      = "Proportion of time in High zone (s)",
  prop_time_ac_medium    = "Proportion of time in Medium zone (s)",
  prop_time_ac_low       = "Proportion of time in Low zone (s)",
  prop_time_ac_calm      = "Proportion of time in Calm zone (s)",
  bl_cm_trial            = "Body length (cm, trial mean)"
)

# Facet strip labels for timepoint
TIMEPOINT_LABELS <- c(
  "1" = "5\u201325 min",
  "2" = "45\u201365 min",
  "3" = "85\u2013105 min"
)
TIMEPOINT_LABELLER <- ggplot2::labeller(timepoint_f = TIMEPOINT_LABELS)


# =============================================================================
# ==== 5) SUMMARY HELPERS =====================================================
# =============================================================================

sem <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  sd(x, na.rm = TRUE) / sqrt(length(x))
}

make_summary_std <- function(df, y_col) {
  df %>%
    dplyr::group_by(treatment) %>%
    dplyr::summarise(
      mean_y = mean(.data[[y_col]], na.rm = TRUE),
      sem_y  = sem(.data[[y_col]]),
      n      = dplyr::n(),
      .groups = "drop"
    )
}

make_summary_std_tp <- function(df, y_col) {
  df %>%
    dplyr::group_by(treatment, timepoint_f) %>%
    dplyr::summarise(
      mean_y = mean(.data[[y_col]], na.rm = TRUE),
      sem_y  = sem(.data[[y_col]]),
      n      = dplyr::n(),
      .groups = "drop"
    )
}

make_summary_zone <- function(df_long) {
  df_long %>%
    dplyr::group_by(zone, treatment) %>%
    dplyr::summarise(
      mean_y = mean(occupancy_s, na.rm = TRUE),
      sem_y  = sem(occupancy_s),
      n      = dplyr::n(),
      .groups = "drop"
    )
}

make_summary_zone_tp <- function(df_long) {
  df_long %>%
    dplyr::group_by(zone, treatment, timepoint_f) %>%
    dplyr::summarise(
      mean_y = mean(occupancy_s, na.rm = TRUE),
      sem_y  = sem(occupancy_s),
      n      = dplyr::n(),
      .groups = "drop"
    )
}

# Area-normalised secondary zone occupancy.
# Divides RAW occupancy (time or proportion) by the zone area constants, then
# re-proportionalises to 100% across the four sub-zones. This is exactly the
# correction STEP2 section 6b applies, so the output equals 100 * prop_time_ac_*.
#
# Input columns must be named: high, medium, low, calm, and must be RAW.
# Passing the already-corrected prop_time_ac_* values applies the area
# correction TWICE. Because this function re-proportionalises, that is not a
# constant offset but a genuine distortion of the composition (it read as
# high 58.2% / medium 21.8% / low 16.9% / calm 3.2% against the correct
# 48.8 / 26.3 / 15.7 / 9.2). See the 2026-08-18 fix at the call site.
normalize_fine_zones <- function(df) {
  df %>%
    dplyr::mutate(
      high_norm   = high   / 10,
      medium_norm = medium / 18,
      low_norm    = low    / 13,
      calm_norm   = calm   / 42
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      total_fine = sum(c_across(c(high_norm, medium_norm, low_norm, calm_norm)),
                       na.rm = TRUE),
      high_pct   = if (total_fine > 0) high_norm   / total_fine * 100 else NA_real_,
      medium_pct = if (total_fine > 0) medium_norm / total_fine * 100 else NA_real_,
      low_pct    = if (total_fine > 0) low_norm    / total_fine * 100 else NA_real_,
      calm_pct   = if (total_fine > 0) calm_norm   / total_fine * 100 else NA_real_
    ) %>%
    dplyr::ungroup()
}

# Apply normalization to df_sum immediately — adds high_pct/medium_pct/low_pct/calm_pct.
#
# THE 2026-08-18 DOUBLE-AREA-NORMALISATION FIX
#   This block used to pass the prop_time_ac_* columns, which STEP2 section 6b
#   has ALREADY divided by the zone areas and rescaled to sum to 1.
#   normalize_fine_zones() divided by the areas again, so the plotted
#   percentages were a distorted composition rather than the area-corrected
#   one. Feeding the RAW proportions instead is what the function documents,
#   and its single division then reproduces STEP2 exactly:
#   high_pct == 100 * prop_time_ac_high.
if (all(c("prop_time_in_high", "prop_time_in_medium",
          "prop_time_in_low",  "prop_time_in_calm_sec") %in% names(df_sum))) {
  df_sum <- df_sum %>%
    dplyr::mutate(high   = prop_time_in_high,
                  medium = prop_time_in_medium,
                  low    = prop_time_in_low,
                  calm   = prop_time_in_calm_sec) %>%
    normalize_fine_zones() %>%
    dplyr::select(-high, -medium, -low, -calm,
                  -high_norm, -medium_norm, -low_norm, -calm_norm, -total_fine)
  ts_msg("Sub-zone occupancy normalized (high_pct / medium_pct / low_pct / calm_pct added to df_sum)")
}


# =============================================================================
# ==== 6) PLOT BUILDER: STANDARD INDICATOR (AGGREGATED) =======================
# =============================================================================

POS_JITTER <- ggplot2::position_jitter(width = JITTER_WIDTH, height = 0)

make_scatter_plot <- function(df, y_col, y_label, title,
                               caption_note = NULL) {
  smry        <- make_summary_std(df, y_col)
  n_obs       <- sum(is.finite(df[[y_col]]))
  caption_txt <- paste0("n (trials) = ", n_obs)
  if (!is.null(caption_note)) caption_txt <- paste0(caption_txt, "  ", caption_note)

  # N-label y position: just below data minimum, ~10% of range below zero-floor.
  # Matches reference style: paste0("(N = ", N, ")"), size=9, vjust=1.5, colour="black".
  .y_vals  <- df[[y_col]][is.finite(df[[y_col]])]
  .y_rng   <- if (length(.y_vals) > 1) diff(range(.y_vals)) else max(abs(.y_vals), 1e-6)
  .n_y_pos <- min(.y_vals, na.rm = TRUE) - 0.10 * max(.y_rng, 1e-6)
  smry$.n_y <- .n_y_pos

  ggplot2::ggplot(df, ggplot2::aes(
    x      = treatment,
    y      = .data[[y_col]],
    colour = treatment,
    shape  = treatment
  )) +
    ggplot2::geom_point(position = POS_JITTER, alpha = PT_ALPHA, size = PT_SIZE) +
    ggplot2::geom_crossbar(
      data      = smry,
      ggplot2::aes(y = mean_y, ymin = mean_y, ymax = mean_y),
      colour    = "black",
      width     = MEAN_WIDTH,
      linewidth = MEAN_LINE_SIZE,
      show.legend = FALSE
    ) +
    ggplot2::geom_errorbar(
      data      = smry,
      ggplot2::aes(y = mean_y, ymin = mean_y - sem_y, ymax = mean_y + sem_y),
      colour    = "black",
      width     = ERRORBAR_WIDTH,
      linewidth = ERRORBAR_SIZE,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data        = smry,
      ggplot2::aes(x = treatment, y = .n_y, label = paste0("(N = ", n, ")")),
      size        = 9,
      vjust       = 1.5,
      colour      = "black",
      inherit.aes = FALSE
    ) +
    ggplot2::scale_colour_manual(values = TREATMENT_COLORS, name = "Treatment") +
    ggplot2::scale_shape_manual(values  = TREATMENT_SHAPES, name = "Treatment") +
    ggplot2::scale_x_discrete(labels = function(x) stringr::str_to_title(x)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.18, 0.05))) +
    BASE_THEME +
    ggplot2::labs(title = title, x = "Treatment", y = y_label, caption = caption_txt)
}


# =============================================================================
# ==== 7) PLOT BUILDER: STANDARD INDICATOR (TIMEPOINT, FACETED) ===============
# =============================================================================

make_scatter_plot_tp <- function(df, y_col, y_label, title) {
  smry        <- make_summary_std_tp(df, y_col)
  n_obs       <- sum(is.finite(df[[y_col]]))
  caption_txt <- paste0("n (trials \u00d7 timepoint) = ", n_obs)

  ggplot2::ggplot(df, ggplot2::aes(
    x      = treatment,
    y      = .data[[y_col]],
    colour = treatment,
    shape  = treatment
  )) +
    ggplot2::geom_point(position = POS_JITTER, alpha = PT_ALPHA, size = PT_SIZE) +
    ggplot2::geom_crossbar(
      data      = smry,
      ggplot2::aes(y = mean_y, ymin = mean_y, ymax = mean_y),
      colour    = "black",
      width     = MEAN_WIDTH,
      linewidth = MEAN_LINE_SIZE,
      show.legend = FALSE
    ) +
    ggplot2::geom_errorbar(
      data      = smry,
      ggplot2::aes(y = mean_y, ymin = mean_y - sem_y, ymax = mean_y + sem_y),
      colour    = "black",
      width     = ERRORBAR_WIDTH,
      linewidth = ERRORBAR_SIZE,
      show.legend = FALSE
    ) +
    ggplot2::facet_wrap(~ timepoint_f, nrow = 1, labeller = TIMEPOINT_LABELLER) +
    ggplot2::scale_colour_manual(values = TREATMENT_COLORS, name = "Treatment") +
    ggplot2::scale_shape_manual(values  = TREATMENT_SHAPES, name = "Treatment") +
    ggplot2::scale_x_discrete(labels = function(x) stringr::str_to_title(x)) +
    BASE_THEME +
    ggplot2::labs(title = title, x = "Treatment", y = y_label, caption = caption_txt)
}


# =============================================================================
# ==== 8) PLOT BUILDER: ZONE PLOT (AGGREGATED) ================================
# =============================================================================

POS_JD    <- ggplot2::position_jitterdodge(jitter.width  = JITTER_WIDTH,
                                            jitter.height = 0,
                                            dodge.width   = DODGE_WIDTH)
POS_DODGE <- ggplot2::position_dodge(width = DODGE_WIDTH)

make_zone_plot <- function(df_long, zone_colors, y_label, title,
                            caption_note = NULL) {
  smry   <- make_summary_zone(df_long)

  ggplot2::ggplot(df_long, ggplot2::aes(
    x      = zone,
    y      = occupancy_s,
    colour = zone,
    shape  = treatment,
    group  = treatment
  )) +
    ggplot2::geom_point(position = POS_JD, alpha = PT_ALPHA, size = ZONE_PT_SIZE) +
    ggplot2::geom_crossbar(
      data      = smry,
      ggplot2::aes(y = mean_y, ymin = mean_y, ymax = mean_y,
                   colour = zone, group = treatment),
      position  = POS_DODGE,
      width     = MEAN_WIDTH,
      linewidth = MEAN_LINE_SIZE,
      show.legend = FALSE
    ) +
    ggplot2::geom_errorbar(
      data      = smry,
      ggplot2::aes(y = mean_y, ymin = mean_y - sem_y, ymax = mean_y + sem_y,
                   colour = zone, group = treatment),
      position  = POS_DODGE,
      width     = ERRORBAR_WIDTH,
      linewidth = ERRORBAR_SIZE,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data        = smry,
      ggplot2::aes(x = zone, y = -8, label = paste0("(N = ", n, ")"),
                   group = treatment),
      position    = POS_DODGE,
      size        = 9,
      vjust       = 1.5,
      colour      = "black",
      inherit.aes = FALSE
    ) +
    # Zone colour legend removed (guide = "none"); Treatment shape legend kept
    ggplot2::scale_colour_manual(values = zone_colors, guide = "none") +
    ggplot2::scale_shape_manual(values  = TREATMENT_SHAPES, name = "Treatment") +
    # Lower limit extended to -20 so (N = X) labels sit clearly below y = 0
    ggplot2::scale_y_continuous(limits = c(-20, 100),
                                breaks = seq(0, 100, by = 25)) +
    ZONE_THEME +
    # Restore x-axis title (ZONE_THEME blanks it; override here for zone plots)
    ggplot2::theme(
      axis.title.x  = ggplot2::element_text(size = 20, face = "bold"),
      legend.position = "bottom"
    ) +
    # No caption — (N = X) per group shown inline via geom_text instead
    ggplot2::labs(title = title, y = y_label, x = "% of time in zone")
}


# =============================================================================
# ==== 9) PLOT BUILDER: ZONE PLOT (TIMEPOINT, FACETED) ========================
# =============================================================================

make_zone_plot_tp <- function(df_long, zone_colors, y_label, title) {
  smry        <- make_summary_zone_tp(df_long)
  n_obs       <- dplyr::n_distinct(df_long$fish_uid[is.finite(df_long$occupancy_s)])
  caption_txt <- paste0("n (trials \u00d7 timepoint) = ", n_obs)

  ggplot2::ggplot(df_long, ggplot2::aes(
    x      = zone,
    y      = occupancy_s,
    colour = zone,
    shape  = treatment,
    group  = treatment
  )) +
    ggplot2::geom_point(position = POS_JD, alpha = PT_ALPHA, size = ZONE_PT_SIZE) +
    ggplot2::geom_crossbar(
      data      = smry,
      ggplot2::aes(y = mean_y, ymin = mean_y, ymax = mean_y,
                   colour = zone, group = treatment),
      position  = POS_DODGE,
      width     = MEAN_WIDTH,
      linewidth = MEAN_LINE_SIZE,
      show.legend = FALSE
    ) +
    ggplot2::geom_errorbar(
      data      = smry,
      ggplot2::aes(y = mean_y, ymin = mean_y - sem_y, ymax = mean_y + sem_y,
                   colour = zone, group = treatment),
      position  = POS_DODGE,
      width     = ERRORBAR_WIDTH,
      linewidth = ERRORBAR_SIZE,
      show.legend = FALSE
    ) +
    ggplot2::facet_wrap(~ timepoint_f, nrow = 1, labeller = TIMEPOINT_LABELLER) +
    # Zone colour legend removed; Treatment shape legend kept
    ggplot2::scale_colour_manual(values = zone_colors, guide = "none") +
    ggplot2::scale_shape_manual(values  = TREATMENT_SHAPES, name = "Treatment") +
    ggplot2::scale_y_continuous(limits = c(-20, 100),
                                breaks = seq(0, 100, by = 25)) +
    ZONE_THEME +
    ggplot2::theme(
      axis.title.x    = ggplot2::element_text(size = 20, face = "bold"),
      legend.position = "bottom"
    ) +
    ggplot2::labs(title = title, y = y_label, x = "% of time in zone")
}


# =============================================================================
# ==== 10) PPT EXPORT HELPER ==================================================
# =============================================================================

save_pptx <- function(plot_list, pptx_path, slide_w = 10, slide_h = 6) {
  prs <- officer::read_pptx()
  for (nm in names(plot_list)) {
    prs <- officer::add_slide(prs, layout = "Blank", master = "Office Theme")
    prs <- officer::ph_with(
      prs,
      value    = rvg::dml(ggobj = plot_list[[nm]]),
      location = officer::ph_location(left = 0, top = 0, width = slide_w, height = slide_h)
    )
  }
  print(prs, target = pptx_path)
  ts_msg("Saved: ", basename(pptx_path))
}


# =============================================================================
# ==== 11) DEFINE INDICATORS ==================================================
# =============================================================================

INDICATORS <- list(
  prop_active          = list(label = INDICATOR_PRETTY_LABELS[["prop_active"]]),
  switches_per_session = list(label = INDICATOR_PRETTY_LABELS[["switches_per_session"]]),
  zone_flux_per_session = list(label = INDICATOR_PRETTY_LABELS[["zone_flux_per_session"]])
)

# Group-level dynamics (one row per trial x timepoint; from group_dynamics_summary)
GROUP_INDICATORS <- list(
  mean_nnd_cm          = list(label = "Mean nearest-neighbour distance (cm)"),
  mean_iid_cm          = list(label = "Mean inter-individual distance (cm)"),
  mean_polarisation    = list(label = "Group polarisation (0\u20131)"),
  mean_hull_area_cm2   = list(label = "School area (cm\u00b2)"),
  mean_centroid_spd_cm = list(label = "School speed (cm s\u207b\u00b9)")
)

# Individual social context (columns in fish_activity_summary merged from STEP 2b)
SOCIAL_INDICATORS <- list(
  mean_nnd_fish_px  = list(label = "Fish nearest-neighbour distance (px)"),
  mean_centdist_px  = list(label = "Distance to group centroid (px)"),
  mean_turning_rate = list(label = "Mean turning rate (rad frame\u207b\u00b9)")
)


# =============================================================================
# ==== 12) LONG-FORMAT ZONE DATA HELPER =======================================
# =============================================================================

.make_main_long <- function(d) {
  # Prefer prop_time_in_* (already proportions); fall back to time_in_*_s.
  # Multiply proportions by 100 so y-axis is % of session time.
  if (all(c("prop_time_in_flow", "prop_time_in_calm") %in% names(d))) {
    d %>%
      dplyr::select(fish_uid, treatment,
                    dplyr::any_of(c("timepoint", "timepoint_f")),
                    prop_time_in_flow, prop_time_in_calm) %>%
      tidyr::pivot_longer(
        cols      = c(prop_time_in_flow, prop_time_in_calm),
        names_to  = "zone_raw",
        values_to = "occupancy_s"
      ) %>%
      dplyr::mutate(
        occupancy_s = occupancy_s * 100,
        zone = factor(
          dplyr::case_when(
            zone_raw == "prop_time_in_flow" ~ "flow",
            zone_raw == "prop_time_in_calm" ~ "calm"
          ),
          levels = c("flow", "calm")
        )
      )
  } else if (all(c("time_in_flow_s", "time_in_calm_s") %in% names(d))) {
    d %>%
      dplyr::select(fish_uid, treatment,
                    dplyr::any_of(c("timepoint", "timepoint_f")),
                    time_in_flow_s, time_in_calm_s) %>%
      tidyr::pivot_longer(
        cols      = c(time_in_flow_s, time_in_calm_s),
        names_to  = "zone_raw",
        values_to = "occupancy_s"
      ) %>%
      dplyr::mutate(
        zone = factor(
          dplyr::case_when(
            zone_raw == "time_in_flow_s" ~ "flow",
            zone_raw == "time_in_calm_s" ~ "calm"
          ),
          levels = c("flow", "calm")
        )
      )
  } else NULL
}

.make_sec_long <- function(d) {
  # normalize_fine_zones() was already applied to df_sum; _pct columns present.
  .needed <- c("high_pct", "medium_pct", "low_pct", "calm_pct")
  if (!all(.needed %in% names(d))) return(NULL)

  d %>%
    dplyr::select(dplyr::any_of(c("trial_id", "tank", "fish_uid", "treatment",
                                   "timepoint", "timepoint_f")),
                  high_pct, medium_pct, low_pct, calm_pct) %>%
    tidyr::pivot_longer(
      cols      = c(high_pct, medium_pct, low_pct, calm_pct),
      names_to  = "zone_raw",
      values_to = "occupancy_s"
    ) %>%
    dplyr::mutate(
      zone = factor(
        dplyr::case_when(
          zone_raw == "high_pct"   ~ "high",
          zone_raw == "medium_pct" ~ "medium",
          zone_raw == "low_pct"    ~ "low",
          zone_raw == "calm_pct"   ~ "calm"
        ),
        levels = c("high", "medium", "low", "calm")
      )
    )
}


# =============================================================================
# ==== 13) BLOCK A — TIMEPOINT PLOTS ==========================================
# =============================================================================

.has_timepoints <- any(!is.na(df_sum$timepoint))
scatter_plots_tp <- list()
zone_plots_tp    <- list()

if (.has_timepoints) {
  ts_msg("=== BLOCK A: TIMEPOINT PLOTS ===")

  .tp_out <- file.path(STEP4_OUT, "timepoint")

  # Standard indicator plots (faceted by timepoint)
  for (ind_name in names(INDICATORS)) {
    if (!ind_name %in% names(df_sum)) {
      ts_msg("  Skipping (column not found): ", ind_name); next
    }
    ts_msg("  Plotting (tp): ", ind_name)
    p <- tryCatch(
      make_scatter_plot_tp(df_sum, ind_name,
                           INDICATORS[[ind_name]]$label,
                           INDICATORS[[ind_name]]$label),
      error = function(e) {
        warning("Scatter tp plot failed for ", ind_name, ": ", conditionMessage(e)); NULL
      }
    )
    if (!is.null(p)) scatter_plots_tp[[ind_name]] <- p
  }

  # Zone plots (faceted by timepoint)
  df_main_long_tp <- .make_main_long(df_sum)
  if (!is.null(df_main_long_tp)) {
    ts_msg("  Main zones (tp): ", nrow(df_main_long_tp), " rows")
    p_main_tp <- tryCatch(
      make_zone_plot_tp(df_main_long_tp, MAIN_ZONE_COLORS,
                        "Time in zone (% of session)", "Main zone occupancy"),
      error = function(e) { warning("Main zone tp plot failed: ", conditionMessage(e)); NULL }
    )
    if (!is.null(p_main_tp)) zone_plots_tp[["main_zones"]] <- p_main_tp
  }

  df_sec_long_tp <- .make_sec_long(df_sum)
  if (!is.null(df_sec_long_tp)) {
    ts_msg("  Secondary zones (tp): ", nrow(df_sec_long_tp), " rows")
    p_sec_tp <- tryCatch({
      p <- make_zone_plot_tp(df_sec_long_tp, SEC_ZONE_COLORS,
                             "Area-normalised occupancy (% of sub-zone time)", "Sub-zone occupancy")
      # Picture4 logo: upper-right inset
      add_logo_top_right(p)
    }, error = function(e) { warning("Sec zone tp plot failed: ", conditionMessage(e)); NULL })
    if (!is.null(p_sec_tp)) zone_plots_tp[["sec_zones"]] <- p_sec_tp
  }

  # Save PPTX
  if (length(scatter_plots_tp) > 0)
    tryCatch(save_pptx(scatter_plots_tp,
                        file.path(.tp_out, "standard_indicator_plots_tp.pptx"),
                        slide_w = 13, slide_h = 6),
             error = function(e)
               warning("Could not save standard tp PPTX: ", conditionMessage(e), call. = FALSE))

  if (length(zone_plots_tp) > 0)
    tryCatch(save_pptx(zone_plots_tp,
                        file.path(.tp_out, "zone_plots_tp.pptx"),
                        slide_w = 13, slide_h = 6),
             error = function(e)
               warning("Could not save zone tp PPTX: ", conditionMessage(e), call. = FALSE))

  # PNG previews
  .png_tp <- file.path(.tp_out, "png_preview")
  invisible(lapply(names(c(scatter_plots_tp, zone_plots_tp)), function(nm) {
    tryCatch(
      ggplot2::ggsave(file.path(.png_tp, paste0(nm, ".png")),
                      c(scatter_plots_tp, zone_plots_tp)[[nm]],
                      width = 13, height = 6, dpi = 150),
      error = function(e) NULL
    )
  }))
  ts_msg("Timepoint PNG previews: ", .png_tp)

  # ---- Group dynamics plots (timepoint, faceted) ----
  group_plots_tp   <- list()
  social_plots_tp  <- list()

  if (!is.null(df_gd) && nrow(df_gd) > 0 && any(!is.na(df_gd$timepoint))) {
    ts_msg("  Group dynamics plots (tp):")
    for (ind_name in names(GROUP_INDICATORS)) {
      if (!ind_name %in% names(df_gd)) {
        ts_msg("    Skipping (not found): ", ind_name); next
      }
      ts_msg("    ", ind_name)
      p <- tryCatch(
        make_scatter_plot_tp(df_gd, ind_name,
                             GROUP_INDICATORS[[ind_name]]$label,
                             GROUP_INDICATORS[[ind_name]]$label),
        error = function(e) { warning("Group tp plot failed for ", ind_name, ": ",
                                      conditionMessage(e)); NULL }
      )
      if (!is.null(p)) group_plots_tp[[ind_name]] <- p
    }
    if (length(group_plots_tp) > 0)
      tryCatch(save_pptx(group_plots_tp,
                          file.path(.tp_out, "group_dynamics_plots_tp.pptx"),
                          slide_w = 13, slide_h = 6),
               error = function(e)
                 warning("Could not save group tp PPTX: ", conditionMessage(e), call. = FALSE))
  }

  # Individual social context (from df_sum, which has these columns post-STEP2b)
  .soc_avail_sum <- intersect(names(SOCIAL_INDICATORS), names(df_sum))
  if (length(.soc_avail_sum) > 0) {
    ts_msg("  Individual social context plots (tp):")
    for (ind_name in .soc_avail_sum) {
      ts_msg("    ", ind_name)
      p <- tryCatch(
        make_scatter_plot_tp(df_sum, ind_name,
                             SOCIAL_INDICATORS[[ind_name]]$label,
                             SOCIAL_INDICATORS[[ind_name]]$label),
        error = function(e) { warning("Social tp plot failed for ", ind_name, ": ",
                                      conditionMessage(e)); NULL }
      )
      if (!is.null(p)) social_plots_tp[[ind_name]] <- p
    }
    if (length(social_plots_tp) > 0)
      tryCatch(save_pptx(social_plots_tp,
                          file.path(.tp_out, "social_context_plots_tp.pptx"),
                          slide_w = 13, slide_h = 6),
               error = function(e)
                 warning("Could not save social tp PPTX: ", conditionMessage(e), call. = FALSE))
  }

  # PNG previews for group + social
  .all_gd_tp <- c(group_plots_tp, social_plots_tp)
  invisible(lapply(names(.all_gd_tp), function(nm) {
    tryCatch(
      ggplot2::ggsave(file.path(.png_tp, paste0("gd_", nm, ".png")),
                      .all_gd_tp[[nm]], width = 13, height = 6, dpi = 150),
      error = function(e) NULL
    )
  }))

} else {
  ts_msg("All timepoints NA — Block A (timepoint plots) skipped.")
  group_plots_tp  <- list()
  social_plots_tp <- list()
}


# =============================================================================
# ==== 14) BLOCK B — AGGREGATED PLOTS =========================================
# =============================================================================

ts_msg("=== BLOCK B: AGGREGATED PLOTS ===")
.agg_out <- file.path(STEP4_OUT, "aggregated")
.agg_note <- if (.has_timepoints) "(all timepoints combined)" else NULL
scatter_plots_agg <- list()
zone_plots_agg    <- list()

# Standard indicator plots
for (ind_name in names(INDICATORS)) {
  if (!ind_name %in% names(df_sum)) {
    ts_msg("  Skipping (column not found in df_sum): ", ind_name); next
  }
  ts_msg("  Plotting (agg): ", ind_name)
  p <- tryCatch(
    make_scatter_plot(df_sum, ind_name,
                      INDICATORS[[ind_name]]$label,
                      INDICATORS[[ind_name]]$label,
                      caption_note = .agg_note),
    error = function(e) {
      warning("Scatter agg plot failed for ", ind_name, ": ", conditionMessage(e)); NULL
    }
  )
  if (!is.null(p)) scatter_plots_agg[[ind_name]] <- p
}

# Zone plots
df_main_long_agg <- .make_main_long(df_sum)
if (!is.null(df_main_long_agg)) {
  ts_msg("  Main zones (agg): ", nrow(df_main_long_agg), " rows")
  p_main_agg <- tryCatch(
    make_zone_plot(df_main_long_agg, MAIN_ZONE_COLORS,
                   "Time in zone (% of session)", "Main zone occupancy",
                   caption_note = .agg_note),
    error = function(e) { warning("Main zone agg plot failed: ", conditionMessage(e)); NULL }
  )
  if (!is.null(p_main_agg)) zone_plots_agg[["main_zones"]] <- p_main_agg
}

df_sec_long_agg <- .make_sec_long(df_sum)
if (!is.null(df_sec_long_agg)) {
  ts_msg("  Secondary zones (agg): ", nrow(df_sec_long_agg), " rows")
  p_sec_agg <- tryCatch({
    p <- make_zone_plot(df_sec_long_agg, SEC_ZONE_COLORS,
                        "Area-normalised occupancy (% of sub-zone time)", "Sub-zone occupancy",
                        caption_note = .agg_note)
    # Picture4 logo: upper-right inset
    add_logo_top_right(p)
  }, error = function(e) { warning("Sec zone agg plot failed: ", conditionMessage(e)); NULL })
  if (!is.null(p_sec_agg)) zone_plots_agg[["sec_zones"]] <- p_sec_agg
}

# Save PPTX
if (length(scatter_plots_agg) > 0)
  tryCatch(save_pptx(scatter_plots_agg,
                      file.path(.agg_out, "standard_indicator_plots_agg.pptx")),
           error = function(e)
             warning("Could not save standard agg PPTX: ", conditionMessage(e), call. = FALSE))

if (length(zone_plots_agg) > 0)
  tryCatch(save_pptx(zone_plots_agg,
                      file.path(.agg_out, "zone_plots_agg.pptx")),
           error = function(e)
             warning("Could not save zone agg PPTX: ", conditionMessage(e), call. = FALSE))

# PNG previews
.png_agg <- file.path(.agg_out, "png_preview")
invisible(lapply(names(c(scatter_plots_agg, zone_plots_agg)), function(nm) {
  tryCatch(
    ggplot2::ggsave(file.path(.png_agg, paste0(nm, ".png")),
                    c(scatter_plots_agg, zone_plots_agg)[[nm]],
                    width = 10, height = 6, dpi = 150),
    error = function(e) NULL
  )
}))
ts_msg("Aggregated PNG previews: ", .png_agg)

# ---- Group dynamics plots (aggregated) ----
group_plots_agg  <- list()
social_plots_agg <- list()

if (!is.null(df_gd) && nrow(df_gd) > 0) {
  ts_msg("  Group dynamics plots (agg):")
  for (ind_name in names(GROUP_INDICATORS)) {
    if (!ind_name %in% names(df_gd)) {
      ts_msg("    Skipping (not found): ", ind_name); next
    }
    ts_msg("    ", ind_name)
    p <- tryCatch(
      make_scatter_plot(df_gd, ind_name,
                        GROUP_INDICATORS[[ind_name]]$label,
                        GROUP_INDICATORS[[ind_name]]$label,
                        caption_note = if (.has_timepoints) "(all timepoints combined)" else NULL),
      error = function(e) { warning("Group agg plot failed for ", ind_name, ": ",
                                    conditionMessage(e)); NULL }
    )
    if (!is.null(p)) group_plots_agg[[ind_name]] <- p
  }
  if (length(group_plots_agg) > 0)
    tryCatch(save_pptx(group_plots_agg,
                        file.path(.agg_out, "group_dynamics_plots_agg.pptx")),
             error = function(e)
               warning("Could not save group agg PPTX: ", conditionMessage(e), call. = FALSE))
}

.soc_avail_sum2 <- intersect(names(SOCIAL_INDICATORS), names(df_sum))
if (length(.soc_avail_sum2) > 0) {
  ts_msg("  Individual social context plots (agg):")
  for (ind_name in .soc_avail_sum2) {
    ts_msg("    ", ind_name)
    p <- tryCatch(
      make_scatter_plot(df_sum, ind_name,
                        SOCIAL_INDICATORS[[ind_name]]$label,
                        SOCIAL_INDICATORS[[ind_name]]$label,
                        caption_note = if (.has_timepoints) "(all timepoints combined)" else NULL),
      error = function(e) { warning("Social agg plot failed for ", ind_name, ": ",
                                    conditionMessage(e)); NULL }
    )
    if (!is.null(p)) social_plots_agg[[ind_name]] <- p
  }
  if (length(social_plots_agg) > 0)
    tryCatch(save_pptx(social_plots_agg,
                        file.path(.agg_out, "social_context_plots_agg.pptx")),
             error = function(e)
               warning("Could not save social agg PPTX: ", conditionMessage(e), call. = FALSE))
}

# PNG previews for group + social (aggregated)
.all_gd_agg <- c(group_plots_agg, social_plots_agg)
invisible(lapply(names(.all_gd_agg), function(nm) {
  tryCatch(
    ggplot2::ggsave(file.path(.png_agg, paste0("gd_", nm, ".png")),
                    .all_gd_agg[[nm]], width = 10, height = 6, dpi = 150),
    error = function(e) NULL
  )
}))


# =============================================================================
# ==== 15) PUBLISH TO GLOBALENV ===============================================
# =============================================================================

assign("scatter_plots_tp_choice",  scatter_plots_tp,  envir = .GlobalEnv)
assign("zone_plots_tp_choice",     zone_plots_tp,     envir = .GlobalEnv)
assign("scatter_plots_agg_choice", scatter_plots_agg, envir = .GlobalEnv)
assign("zone_plots_agg_choice",    zone_plots_agg,    envir = .GlobalEnv)
assign("group_plots_tp_choice",    group_plots_tp,    envir = .GlobalEnv)
assign("social_plots_tp_choice",   social_plots_tp,   envir = .GlobalEnv)
assign("group_plots_agg_choice",   group_plots_agg,   envir = .GlobalEnv)
assign("social_plots_agg_choice",  social_plots_agg,  envir = .GlobalEnv)
assign("STEP4_OUTPUT_DIR",         STEP4_OUT,         envir = .GlobalEnv)

ts_msg("Step 4 complete.")
ts_msg("  Timepoint set:   ",
       length(scatter_plots_tp), " standard + ",
       length(zone_plots_tp), " zone + ",
       length(group_plots_tp), " group + ",
       length(social_plots_tp), " social plots -> timepoint/")
ts_msg("  Aggregated set:  ",
       length(scatter_plots_agg), " standard + ",
       length(zone_plots_agg), " zone + ",
       length(group_plots_agg), " group + ",
       length(social_plots_agg), " social plots -> aggregated/")
ts_msg("  Output root: ", STEP4_OUT)
