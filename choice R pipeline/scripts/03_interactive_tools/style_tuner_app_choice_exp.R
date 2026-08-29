#!/usr/bin/env Rscript
# =============================================================================
# STYLE TUNER — CHOICE EXPERIMENT PLOTS
# Live-preview app for tuning all visual parameters of the scatter plots
# produced by activity_analysis_GRAPHS_choice_exp.R.
#
# PLOT TYPES:
#   Standard indicator  — x = treatment, colour + shape = treatment
#   Main zones          — x = zone (flow/calm), colour = zone, shape = treatment
#   Secondary zones     — x = zone (high/medium/low/calm), colour = zone, shape = treatment
#
# HOW TO USE:
#   1. Run Steps 1-2 first so fish_activity_summary is in GlobalEnv.
#      OR open this script with the pipeline working directory set — it will
#      search STEP2_output/ automatically.
#   2. Source this file: click Source in RStudio.
#   3. Adjust sliders / inputs in the sidebar — the preview updates live.
#   4. Use the "Plot code" tab to copy the final settings back into the
#      pipeline script.
#   5. Use the Download button to save the current plot as a PNG.
# =============================================================================

options(stringsAsFactors = FALSE)

required_pkgs <- c("shiny", "ggplot2", "dplyr", "tidyr", "readr", "stringr", "colourpicker")
missing_pkgs  <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]
if (length(missing_pkgs) > 0) {
  install.packages(missing_pkgs, repos = "https://cloud.r-project.org")
}
invisible(lapply(required_pkgs, library, character.only = TRUE))


# =============================================================================
# ==== HELPERS ================================================================
# =============================================================================

sem <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  sd(x) / sqrt(length(x))
}

.find_latest_csv <- function(step_name, csv_filename) {
  parent <- file.path(getwd(), step_name)
  if (!dir.exists(parent)) return(NULL)
  subdirs <- list.dirs(parent, full.names = TRUE, recursive = FALSE)
  subdirs <- subdirs[grepl(paste0("^", step_name, "_\\d{8}_\\d{6}$"),
                           basename(subdirs))]
  if (length(subdirs) == 0) return(NULL)
  latest    <- subdirs[which.max(file.mtime(subdirs))]
  candidate <- file.path(latest, csv_filename)
  if (file.exists(candidate)) candidate else NULL
}


# =============================================================================
# ==== LOAD DATA ==============================================================
# =============================================================================

if (exists("fish_activity_summary", envir = .GlobalEnv, inherits = FALSE)) {
  .fas <- get("fish_activity_summary", envir = .GlobalEnv, inherits = FALSE)
  message("Using fish_activity_summary from GlobalEnv (", nrow(.fas), " rows)")
} else {
  message("fish_activity_summary not in GlobalEnv; searching STEP2_output/...")
  .csv <- .find_latest_csv("STEP2_output", "fish_activity_summary.csv")
  if (is.null(.csv))
    stop("Cannot find fish_activity_summary. Run Steps 1-2 first.", call. = FALSE)
  .fas <- readr::read_csv(.csv, show_col_types = FALSE)
  assign("fish_activity_summary", .fas, envir = .GlobalEnv)
  message("Loaded from: ", .csv)
}

TREATMENT_LEVELS_g <- if (exists("TREATMENT_LEVELS", envir = .GlobalEnv, inherits = FALSE))
  get("TREATMENT_LEVELS", envir = .GlobalEnv) else c("control", "exercise choice")

.fas <- .fas %>%
  dplyr::mutate(
    treatment = factor(trimws(tolower(as.character(treatment))),
                       levels = TREATMENT_LEVELS_g)
  )

# Identify numeric indicator columns (exclude ID / grouping cols)
.exclude_cols <- c("trial_id", "trial_date", "fish_id", "fish_in_session",
                   "treatment", "fish_density", "fish_density_f", "tank",
                   "motor_side", "total_distance_cm")
INDICATORS_AVAILABLE <- setdiff(
  names(.fas)[vapply(.fas, is.numeric, logical(1))],
  .exclude_cols
)

INDICATOR_LABELS <- c(
  mean_speed_cm_s         = "Mean speed (cm/s)",
  max_speed_cm_s          = "Max speed (cm/s)",
  total_distance_m        = "Total distance (m)",
  moving_time_TS1_s       = "Time moving — TS1 (s)",
  prop_moving_TS1         = "Proportion time moving — TS1",
  mean_speed_BL_s         = "Mean speed (BL/s)",
  time_in_flow_s          = "Time in flow zone (s)",
  time_in_calm_s          = "Time in calm zone (s)",
  prop_time_in_flow       = "Proportion time in flow zone",
  time_in_flow_high_s     = "Time — high flow (s)",
  time_in_flow_medium_s   = "Time — medium flow (s)",
  time_in_flow_low_s      = "Time — low flow (s)",
  area_corr_s_flow_high   = "Area-corrected time — high flow (s/unit)",
  area_corr_s_flow_medium = "Area-corrected time — medium flow (s/unit)",
  area_corr_s_flow_low    = "Area-corrected time — low flow (s/unit)",
  area_corr_s_calm        = "Area-corrected time — calm (s/unit)",
  norm_pct_flow_high      = "Area-corrected occupancy — high flow (%)",
  norm_pct_flow_medium    = "Area-corrected occupancy — medium flow (%)",
  norm_pct_flow_low       = "Area-corrected occupancy — low flow (%)",
  norm_pct_calm           = "Area-corrected occupancy — calm (%)"
)
INDICATOR_LABELS <- INDICATOR_LABELS[names(INDICATOR_LABELS) %in% INDICATORS_AVAILABLE]
.unlabelled      <- setdiff(INDICATORS_AVAILABLE, names(INDICATOR_LABELS))
INDICATOR_LABELS <- c(INDICATOR_LABELS, setNames(.unlabelled, .unlabelled))
INDICATOR_CHOICES <- setNames(names(INDICATOR_LABELS), INDICATOR_LABELS)

# Pre-build zone long formats (rebuilt inside reactive if needed for filtering)
MAIN_ZONE_LEVELS <- c("flow", "calm")
SEC_ZONE_LEVELS  <- c("high", "medium", "low", "calm")

.has_main_zones <- all(c("time_in_flow_s", "time_in_calm_s") %in% names(.fas))
.has_sec_zones  <- all(c("area_corr_s_flow_high", "area_corr_s_flow_medium",
                          "area_corr_s_flow_low",  "area_corr_s_calm") %in% names(.fas))


# =============================================================================
# ==== STYLE DEFAULTS =========================================================
# =============================================================================

STYLE_DEFAULTS <- list(
  plot_type         = "standard",
  indicator         = if ("mean_speed_cm_s" %in% INDICATORS_AVAILABLE) "mean_speed_cm_s" else INDICATORS_AVAILABLE[1],
  lower_p           = 0.02,
  upper_p           = 0.98,
  # points
  point_size        = 2.8,
  point_alpha       = 0.85,
  jitter_width      = 0.12,
  dodge_width       = 0.70,
  # error bars & crossbar
  errorbar_width    = 0.14,
  errorbar_lwd      = 0.65,
  crossbar_width    = 0.45,
  crossbar_lwd      = 0.85,
  # treatment colours (always shown)
  col_control       = "#2166AC",
  col_exercise      = "#D6604D",
  # main zone colours
  col_zone_flow     = "#4DAC26",
  col_zone_calm     = "#8073AC",
  # secondary zone colours
  col_zone_high     = "#D73027",
  col_zone_medium   = "#FC8D59",
  col_zone_low      = "#91BFDB",
  col_zone_calm_sec = "#8073AC",
  # theme
  base_size         = 13,
  axis_text_size    = 11,
  axis_title_size   = 14,
  plot_title_size   = 14,
  legend_text_size  = 12,
  legend_title_size = 12,
  show_legend       = TRUE,
  legend_position   = "right",
  show_major_grid   = FALSE,
  show_minor_grid   = FALSE,
  axis_lwd          = 0.85,
  margin_t          = 8,
  margin_r          = 10,
  margin_b          = 16,
  margin_l          = 10,
  # export
  plot_width        = 8,
  plot_height       = 6
)


# =============================================================================
# ==== PLOT BUILDER ===========================================================
# =============================================================================

make_tuned_plot <- function(df, input) {

  # ---- shared theme -----------------------------------------------------------
  base_theme <- ggplot2::theme_minimal(base_size = input$base_size) +
    ggplot2::theme(
      panel.grid.major  = if (isTRUE(input$show_major_grid)) ggplot2::element_line() else ggplot2::element_blank(),
      panel.grid.minor  = if (isTRUE(input$show_minor_grid)) ggplot2::element_line() else ggplot2::element_blank(),
      panel.background  = ggplot2::element_blank(),
      axis.line.x       = ggplot2::element_line(color = "black", linewidth = input$axis_lwd),
      axis.line.y       = ggplot2::element_line(color = "black", linewidth = input$axis_lwd),
      axis.ticks        = ggplot2::element_line(color = "black", linewidth = input$axis_lwd * 0.8),
      axis.text         = ggplot2::element_text(size  = input$axis_text_size, color = "black"),
      axis.text.x       = ggplot2::element_text(face  = "bold"),
      axis.title        = ggplot2::element_text(size  = input$axis_title_size, face = "bold", color = "black"),
      plot.title        = ggplot2::element_text(size  = input$plot_title_size, face = "bold"),
      plot.title.position = "plot",
      legend.title      = ggplot2::element_text(size  = input$legend_title_size, face = "bold"),
      legend.text       = ggplot2::element_text(size  = input$legend_text_size),
      legend.position   = if (isTRUE(input$show_legend)) input$legend_position else "none",
      plot.margin       = ggplot2::margin(t = input$margin_t, r = input$margin_r,
                                          b = input$margin_b, l = input$margin_l)
    )

  treatment_colors <- c(
    setNames(input$col_control,  "control"),
    setNames(input$col_exercise, "exercise choice")
  )
  treatment_shapes <- c("control" = 17L, "exercise choice" = 16L)

  pos_jitter <- ggplot2::position_jitter(width = input$jitter_width, height = 0)
  pos_jd     <- ggplot2::position_jitterdodge(jitter.width  = input$jitter_width,
                                               jitter.height = 0,
                                               dodge.width   = input$dodge_width)
  pos_dodge  <- ggplot2::position_dodge(width = input$dodge_width)

  # ---- standard indicator -----------------------------------------------------
  if (input$plot_type == "standard") {
    ind   <- input$indicator
    req(ind %in% names(df))

    y_raw  <- df[[ind]]
    lo     <- quantile(y_raw, input$lower_p, na.rm = TRUE)
    hi     <- quantile(y_raw, input$upper_p, na.rm = TRUE)
    df_plt <- df[!is.na(y_raw) & y_raw >= lo & y_raw <= hi, ]

    if (nrow(df_plt) == 0)
      return(ggplot2::ggplot() +
               ggplot2::annotate("text", x=0.5, y=0.5, label="No data after outlier filter", size=6) +
               ggplot2::theme_void())

    smry <- df_plt %>%
      dplyr::group_by(treatment) %>%
      dplyr::summarise(mean_y = mean(.data[[ind]], na.rm = TRUE),
                       sem_y  = sem(.data[[ind]]),
                       .groups = "drop")

    y_label     <- INDICATOR_LABELS[ind]; if (is.na(y_label)) y_label <- ind
    caption_txt <- paste0("n = ", nrow(df_plt), " fish",
                          if (nrow(df_plt) < sum(!is.na(y_raw)))
                            paste0(" (", sum(!is.na(y_raw)) - nrow(df_plt), " outliers removed)")
                          else "")

    ggplot2::ggplot(df_plt, ggplot2::aes(x = treatment, y = .data[[ind]],
                                          colour = treatment, shape = treatment)) +
      ggplot2::geom_point(position = pos_jitter,
                          size = input$point_size, alpha = input$point_alpha) +
      ggplot2::geom_crossbar(data = smry,
                              ggplot2::aes(y = mean_y, ymin = mean_y, ymax = mean_y),
                              width = input$crossbar_width, linewidth = input$crossbar_lwd,
                              show.legend = FALSE) +
      ggplot2::geom_errorbar(data = smry,
                              ggplot2::aes(y = mean_y, ymin = mean_y - sem_y, ymax = mean_y + sem_y),
                              width = input$errorbar_width, linewidth = input$errorbar_lwd,
                              show.legend = FALSE) +
      ggplot2::scale_colour_manual(values = treatment_colors, name = "Treatment") +
      ggplot2::scale_shape_manual(values  = treatment_shapes, name = "Treatment") +
      ggplot2::scale_x_discrete(labels = function(x) stringr::str_to_title(x)) +
      base_theme +
      ggplot2::labs(title = y_label, x = "Treatment", y = y_label, caption = caption_txt)

  # ---- main zones -------------------------------------------------------------
  } else if (input$plot_type == "main_zones") {
    if (!.has_main_zones)
      return(ggplot2::ggplot() +
               ggplot2::annotate("text", x=0.5, y=0.5,
                                 label="time_in_flow_s / time_in_calm_s not found in data", size=5) +
               ggplot2::theme_void())

    zone_colors <- c(flow = input$col_zone_flow, calm = input$col_zone_calm)

    df_long <- df %>%
      dplyr::mutate(fish_uid = paste(trial_id, fish_id, sep = "_")) %>%
      dplyr::select(fish_uid, treatment, time_in_flow_s, time_in_calm_s) %>%
      tidyr::pivot_longer(c(time_in_flow_s, time_in_calm_s),
                          names_to = "zone_raw", values_to = "occupancy_s") %>%
      dplyr::mutate(zone = factor(
        dplyr::case_when(zone_raw == "time_in_flow_s" ~ "flow",
                         zone_raw == "time_in_calm_s" ~ "calm"),
        levels = MAIN_ZONE_LEVELS
      ))

    lo <- quantile(df_long$occupancy_s, input$lower_p, na.rm = TRUE)
    hi <- quantile(df_long$occupancy_s, input$upper_p, na.rm = TRUE)
    df_long <- df_long[!is.na(df_long$occupancy_s) &
                         df_long$occupancy_s >= lo & df_long$occupancy_s <= hi, ]

    smry <- df_long %>%
      dplyr::group_by(zone, treatment) %>%
      dplyr::summarise(mean_y = mean(occupancy_s, na.rm = TRUE),
                       sem_y  = sem(occupancy_s), .groups = "drop")

    ggplot2::ggplot(df_long, ggplot2::aes(x = zone, y = occupancy_s,
                                           colour = zone, shape = treatment, group = treatment)) +
      ggplot2::geom_point(position = pos_jd, size = input$point_size, alpha = input$point_alpha) +
      ggplot2::geom_crossbar(data = smry,
                              ggplot2::aes(y = mean_y, ymin = mean_y, ymax = mean_y,
                                           colour = zone, group = treatment),
                              position = pos_dodge, width = input$crossbar_width,
                              linewidth = input$crossbar_lwd, show.legend = FALSE) +
      ggplot2::geom_errorbar(data = smry,
                              ggplot2::aes(y = mean_y, ymin = mean_y - sem_y, ymax = mean_y + sem_y,
                                           colour = zone, group = treatment),
                              position = pos_dodge, width = input$errorbar_width,
                              linewidth = input$errorbar_lwd, show.legend = FALSE) +
      ggplot2::scale_colour_manual(values = zone_colors, name = "Zone") +
      ggplot2::scale_shape_manual(values  = treatment_shapes, name = "Treatment") +
      base_theme +
      ggplot2::labs(title = "Main zone occupancy", x = "Zone", y = "Time in zone (s)",
                    caption = paste0("n = ", dplyr::n_distinct(df_long$fish_uid), " fish"))

  # ---- secondary zones --------------------------------------------------------
  } else {
    if (!.has_sec_zones)
      return(ggplot2::ggplot() +
               ggplot2::annotate("text", x=0.5, y=0.5,
                                 label="area_corr_s_* columns not found in data", size=5) +
               ggplot2::theme_void())

    zone_colors <- c(
      high   = input$col_zone_high,
      medium = input$col_zone_medium,
      low    = input$col_zone_low,
      calm   = input$col_zone_calm_sec
    )

    df_long <- df %>%
      dplyr::mutate(fish_uid = paste(trial_id, fish_id, sep = "_")) %>%
      dplyr::select(fish_uid, treatment,
                    area_corr_s_flow_high, area_corr_s_flow_medium,
                    area_corr_s_flow_low,  area_corr_s_calm) %>%
      tidyr::pivot_longer(c(area_corr_s_flow_high, area_corr_s_flow_medium,
                             area_corr_s_flow_low,  area_corr_s_calm),
                          names_to = "zone_raw", values_to = "occupancy_s") %>%
      dplyr::mutate(zone = factor(
        dplyr::case_when(
          zone_raw == "area_corr_s_flow_high"   ~ "high",
          zone_raw == "area_corr_s_flow_medium" ~ "medium",
          zone_raw == "area_corr_s_flow_low"    ~ "low",
          zone_raw == "area_corr_s_calm"        ~ "calm"
        ),
        levels = SEC_ZONE_LEVELS
      ))

    lo <- quantile(df_long$occupancy_s, input$lower_p, na.rm = TRUE)
    hi <- quantile(df_long$occupancy_s, input$upper_p, na.rm = TRUE)
    df_long <- df_long[!is.na(df_long$occupancy_s) &
                         df_long$occupancy_s >= lo & df_long$occupancy_s <= hi, ]

    smry <- df_long %>%
      dplyr::group_by(zone, treatment) %>%
      dplyr::summarise(mean_y = mean(occupancy_s, na.rm = TRUE),
                       sem_y  = sem(occupancy_s), .groups = "drop")

    ggplot2::ggplot(df_long, ggplot2::aes(x = zone, y = occupancy_s,
                                           colour = zone, shape = treatment, group = treatment)) +
      ggplot2::geom_point(position = pos_jd, size = input$point_size, alpha = input$point_alpha) +
      ggplot2::geom_crossbar(data = smry,
                              ggplot2::aes(y = mean_y, ymin = mean_y, ymax = mean_y,
                                           colour = zone, group = treatment),
                              position = pos_dodge, width = input$crossbar_width,
                              linewidth = input$crossbar_lwd, show.legend = FALSE) +
      ggplot2::geom_errorbar(data = smry,
                              ggplot2::aes(y = mean_y, ymin = mean_y - sem_y, ymax = mean_y + sem_y,
                                           colour = zone, group = treatment),
                              position = pos_dodge, width = input$errorbar_width,
                              linewidth = input$errorbar_lwd, show.legend = FALSE) +
      ggplot2::scale_colour_manual(values = zone_colors, name = "Zone") +
      ggplot2::scale_shape_manual(values  = treatment_shapes, name = "Treatment") +
      base_theme +
      ggplot2::labs(title = "Secondary zone occupancy", x = "Zone",
                    y = "Area-corrected time (s / unit area)",
                    caption = paste0("n = ", dplyr::n_distinct(df_long$fish_uid), " fish"))
  }
}


# =============================================================================
# ==== CODE EXPORT ============================================================
# =============================================================================

build_plot_code <- function(input) {
  common <- paste0(
'# ---- Paste these settings into activity_analysis_GRAPHS_choice_exp.R ----

# Treatment colours + shapes
TREATMENT_COLORS <- c(
  "control"         = "', input$col_control,  '",
  "exercise choice" = "', input$col_exercise, '"
)
TREATMENT_SHAPES <- c(
  "control"         = 17L,
  "exercise choice" = 16L
)
'
  )

  zone_sec <- if (input$plot_type == "main_zones") paste0(
'# Main zone colours
MAIN_ZONE_COLORS <- c(
  "flow" = "', input$col_zone_flow, '",
  "calm" = "', input$col_zone_calm, '"
)
'
  ) else if (input$plot_type == "sec_zones") paste0(
'# Secondary zone colours
SEC_ZONE_COLORS <- c(
  "high"   = "', input$col_zone_high,     '",
  "medium" = "', input$col_zone_medium,   '",
  "low"    = "', input$col_zone_low,      '",
  "calm"   = "', input$col_zone_calm_sec, '"
)
'
  ) else ""

  geom_settings <- paste0(
'# Geom settings
PT_ALPHA       <- ', input$point_alpha, '
PT_SIZE        <- ', input$point_size, '
JITTER_WIDTH   <- ', input$jitter_width, '
DODGE_WIDTH    <- ', input$dodge_width, '
MEAN_WIDTH     <- ', input$crossbar_width, '
MEAN_LINE_SIZE <- ', input$crossbar_lwd, '
ERRORBAR_WIDTH <- ', input$errorbar_width, '
ERRORBAR_SIZE  <- ', input$errorbar_lwd, '

# BASE_THEME
BASE_THEME <- ggplot2::theme_minimal(base_size = ', input$base_size, ') +
  ggplot2::theme(
    panel.grid.major  = ', if (isTRUE(input$show_major_grid)) "ggplot2::element_line()" else "ggplot2::element_blank()", ',
    panel.grid.minor  = ', if (isTRUE(input$show_minor_grid)) "ggplot2::element_line()" else "ggplot2::element_blank()", ',
    panel.background  = ggplot2::element_blank(),
    axis.line.x       = ggplot2::element_line(color = "black", linewidth = ', input$axis_lwd, '),
    axis.line.y       = ggplot2::element_line(color = "black", linewidth = ', input$axis_lwd, '),
    axis.ticks        = ggplot2::element_line(color = "black", linewidth = ', round(input$axis_lwd * 0.8, 3), '),
    axis.text         = ggplot2::element_text(size = ', input$axis_text_size, ', color = "black"),
    axis.text.x       = ggplot2::element_text(face = "bold"),
    axis.title        = ggplot2::element_text(size = ', input$axis_title_size, ', face = "bold", color = "black"),
    plot.title        = ggplot2::element_text(size = ', input$plot_title_size, ', face = "bold"),
    plot.title.position = "plot",
    legend.title      = ggplot2::element_text(size = ', input$legend_title_size, ', face = "bold"),
    legend.text       = ggplot2::element_text(size = ', input$legend_text_size, '),
    legend.position   = ', if (isTRUE(input$show_legend)) paste0('"', input$legend_position, '"') else '"none"', ',
    plot.margin       = ggplot2::margin(t = ', input$margin_t, ', r = ', input$margin_r,
                                        ', b = ', input$margin_b, ', l = ', input$margin_l, ')
  )

# Export size (inches, 300 dpi)
# plot_width  <- ', input$plot_width, '
# plot_height <- ', input$plot_height
  )

  paste0(common, zone_sec, geom_settings)
}

build_data_summary <- function(df) {
  df %>%
    dplyr::group_by(treatment) %>%
    dplyr::summarise(n_fish = dplyr::n(), .groups = "drop") %>%
    dplyr::arrange(treatment)
}


# =============================================================================
# ==== UI =====================================================================
# =============================================================================

ui <- shiny::fluidPage(
  shiny::titlePanel("Choice experiment — plot style tuner"),
  shiny::tags$head(shiny::tags$style(shiny::HTML("
    .control-scroll {
      max-height: calc(100vh - 120px);
      overflow-y: auto;
      padding-right: 8px;
    }
    .preview-sticky {
      position: sticky;
      top: 10px;
      align-self: flex-start;
      background: white;
      z-index: 10;
    }
    hr { margin: 10px 0; }
  "))),
  shiny::sidebarLayout(
    shiny::sidebarPanel(
      width = 3,
      div(class = "control-scroll",

        shiny::h5("Plot type", style = "font-weight:bold; margin-top:0;"),
        shiny::selectInput("plot_type", NULL,
                           choices  = c("Standard indicator" = "standard",
                                        "Main zones"         = "main_zones",
                                        "Secondary zones"    = "sec_zones"),
                           selected = STYLE_DEFAULTS$plot_type),

        shiny::conditionalPanel(
          condition = "input.plot_type == 'standard'",
          shiny::h5("Indicator", style = "font-weight:bold;"),
          shiny::selectInput("indicator", NULL,
                             choices  = INDICATOR_CHOICES,
                             selected = STYLE_DEFAULTS$indicator)
        ),

        shiny::hr(),
        shiny::h5("Outlier filter", style = "font-weight:bold;"),
        shiny::sliderInput("lower_p", "Lower percentile cut",
                           min = 0, max = 0.10, value = STYLE_DEFAULTS$lower_p, step = 0.005),
        shiny::sliderInput("upper_p", "Upper percentile cut",
                           min = 0.90, max = 1, value = STYLE_DEFAULTS$upper_p, step = 0.005),

        shiny::hr(),
        shiny::h5("Points", style = "font-weight:bold;"),
        shiny::sliderInput("point_size",   "Point size",   min = 0.5, max = 7,   value = STYLE_DEFAULTS$point_size,   step = 0.1),
        shiny::sliderInput("point_alpha",  "Point alpha",  min = 0.1, max = 1,   value = STYLE_DEFAULTS$point_alpha,  step = 0.05),
        shiny::sliderInput("jitter_width", "Jitter width", min = 0,   max = 0.5, value = STYLE_DEFAULTS$jitter_width, step = 0.01),
        shiny::conditionalPanel(
          condition = "input.plot_type != 'standard'",
          shiny::sliderInput("dodge_width", "Dodge width", min = 0.2, max = 1.2, value = STYLE_DEFAULTS$dodge_width, step = 0.05)
        ),

        shiny::hr(),
        shiny::h5("Mean marker", style = "font-weight:bold;"),
        shiny::sliderInput("crossbar_width", "Crossbar width",     min = 0.05, max = 1,   value = STYLE_DEFAULTS$crossbar_width, step = 0.05),
        shiny::sliderInput("crossbar_lwd",   "Crossbar linewidth", min = 0.1,  max = 3,   value = STYLE_DEFAULTS$crossbar_lwd,   step = 0.05),
        shiny::sliderInput("errorbar_width", "Errorbar width",     min = 0.01, max = 0.5, value = STYLE_DEFAULTS$errorbar_width, step = 0.01),
        shiny::sliderInput("errorbar_lwd",   "Errorbar linewidth", min = 0.1,  max = 3,   value = STYLE_DEFAULTS$errorbar_lwd,   step = 0.05),

        shiny::hr(),
        shiny::h5("Treatment colours", style = "font-weight:bold;"),
        colourpicker::colourInput("col_control",  "Control",         value = STYLE_DEFAULTS$col_control),
        colourpicker::colourInput("col_exercise", "Exercise choice", value = STYLE_DEFAULTS$col_exercise),

        shiny::conditionalPanel(
          condition = "input.plot_type == 'main_zones'",
          shiny::hr(),
          shiny::h5("Zone colours — main zones", style = "font-weight:bold;"),
          colourpicker::colourInput("col_zone_flow", "Flow", value = STYLE_DEFAULTS$col_zone_flow),
          colourpicker::colourInput("col_zone_calm", "Calm", value = STYLE_DEFAULTS$col_zone_calm)
        ),

        shiny::conditionalPanel(
          condition = "input.plot_type == 'sec_zones'",
          shiny::hr(),
          shiny::h5("Zone colours — secondary zones", style = "font-weight:bold;"),
          colourpicker::colourInput("col_zone_high",     "High flow",   value = STYLE_DEFAULTS$col_zone_high),
          colourpicker::colourInput("col_zone_medium",   "Medium flow", value = STYLE_DEFAULTS$col_zone_medium),
          colourpicker::colourInput("col_zone_low",      "Low flow",    value = STYLE_DEFAULTS$col_zone_low),
          colourpicker::colourInput("col_zone_calm_sec", "Calm",        value = STYLE_DEFAULTS$col_zone_calm_sec)
        ),

        shiny::hr(),
        shiny::h5("Typography", style = "font-weight:bold;"),
        shiny::sliderInput("base_size",        "Base size",         min = 8,  max = 24, value = STYLE_DEFAULTS$base_size,        step = 1),
        shiny::sliderInput("axis_text_size",   "Axis text size",    min = 6,  max = 24, value = STYLE_DEFAULTS$axis_text_size,   step = 1),
        shiny::sliderInput("axis_title_size",  "Axis title size",   min = 6,  max = 28, value = STYLE_DEFAULTS$axis_title_size,  step = 1),
        shiny::sliderInput("plot_title_size",  "Plot title size",   min = 8,  max = 30, value = STYLE_DEFAULTS$plot_title_size,  step = 1),
        shiny::sliderInput("legend_text_size", "Legend text size",  min = 6,  max = 24, value = STYLE_DEFAULTS$legend_text_size, step = 1),
        shiny::sliderInput("legend_title_size","Legend title size", min = 6,  max = 24, value = STYLE_DEFAULTS$legend_title_size,step = 1),

        shiny::hr(),
        shiny::h5("Legend & grid", style = "font-weight:bold;"),
        shiny::checkboxInput("show_legend",     "Show legend",     value = STYLE_DEFAULTS$show_legend),
        shiny::selectInput("legend_position", "Legend position",
                           choices  = c("right", "left", "top", "bottom"),
                           selected = STYLE_DEFAULTS$legend_position),
        shiny::checkboxInput("show_major_grid", "Major grid lines", value = STYLE_DEFAULTS$show_major_grid),
        shiny::checkboxInput("show_minor_grid", "Minor grid lines", value = STYLE_DEFAULTS$show_minor_grid),

        shiny::hr(),
        shiny::h5("Axes & margins", style = "font-weight:bold;"),
        shiny::sliderInput("axis_lwd", "Axis line width", min = 0.1, max = 2.5, value = STYLE_DEFAULTS$axis_lwd,  step = 0.05),
        shiny::sliderInput("margin_t", "Top margin",    min = 0, max = 50, value = STYLE_DEFAULTS$margin_t, step = 1),
        shiny::sliderInput("margin_r", "Right margin",  min = 0, max = 50, value = STYLE_DEFAULTS$margin_r, step = 1),
        shiny::sliderInput("margin_b", "Bottom margin", min = 0, max = 50, value = STYLE_DEFAULTS$margin_b, step = 1),
        shiny::sliderInput("margin_l", "Left margin",   min = 0, max = 50, value = STYLE_DEFAULTS$margin_l, step = 1),

        shiny::hr(),
        shiny::h5("Export size (inches)", style = "font-weight:bold;"),
        shiny::numericInput("plot_width",  "Width",  value = STYLE_DEFAULTS$plot_width,  min = 2, max = 20, step = 0.5),
        shiny::numericInput("plot_height", "Height", value = STYLE_DEFAULTS$plot_height, min = 2, max = 20, step = 0.5),
        shiny::downloadButton("download_plot", "Download PNG", style = "width:100%;")
      )
    ),

    shiny::mainPanel(
      width = 9,
      div(class = "preview-sticky",
        shiny::tabsetPanel(
          shiny::tabPanel("Preview",
            shiny::plotOutput("plot_preview", height = "520px")
          ),
          shiny::tabPanel("Plot code",
            shiny::tags$p(style = "margin-top:8px; color:#555;",
              "Copy-paste into activity_analysis_GRAPHS_choice_exp.R"),
            shiny::tags$pre(
              style = "background:#f5f5f5; padding:12px; font-size:12px; overflow-x:auto;",
              shiny::textOutput("plot_code")
            )
          ),
          shiny::tabPanel("Data summary",
            shiny::tags$p(style = "margin-top:8px; color:#555;", "n fish per treatment"),
            shiny::tableOutput("data_summary")
          )
        )
      )
    )
  )
)


# =============================================================================
# ==== SERVER =================================================================
# =============================================================================

server <- function(input, output, session) {

  tuned_plot <- shiny::reactive({
    req(input$plot_type)
    tryCatch(
      make_tuned_plot(.fas, input),
      error = function(e) {
        ggplot2::ggplot() +
          ggplot2::annotate("text", x = 0.5, y = 0.5,
                            label = conditionMessage(e), size = 5) +
          ggplot2::theme_void()
      }
    )
  })

  output$plot_preview <- shiny::renderPlot({ tuned_plot() })

  output$plot_code <- shiny::renderText({ build_plot_code(input) })

  output$data_summary <- shiny::renderTable({
    build_data_summary(.fas) %>%
      dplyr::rename(Treatment = treatment, `N fish` = n_fish)
  })

  output$download_plot <- shiny::downloadHandler(
    filename = function() {
      suffix <- switch(input$plot_type,
                       standard   = input$indicator,
                       main_zones = "main_zones",
                       sec_zones  = "sec_zones")
      paste0("choice_exp_", suffix, "_styled.png")
    },
    content = function(file) {
      ggplot2::ggsave(
        filename = file,
        plot     = tuned_plot(),
        width    = input$plot_width,
        height   = input$plot_height,
        dpi      = 300
      )
    }
  )
}

shiny::shinyApp(ui, server)
