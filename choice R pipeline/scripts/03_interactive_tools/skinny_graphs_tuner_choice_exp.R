#!/usr/bin/env Rscript
# =============================================================================
# SKINNY GRAPHS STYLE TUNER — CHOICE EXPERIMENT
# Live-preview app for tuning all visual parameters of the six skinny graphs
# produced by activity_analysis_STATS_choice_exp.R (Step 5).
#
# PLOT TYPES:
#   1. Zone occupancy — aggregated      (zone on x, treatment dodge)
#   2. Zone occupancy × timepoint       (zone on x, treatment dodge, faceted)
#   3. Zone entries — aggregated        (zone on x, treatment dodge)
#   4. Zone entries × timepoint         (zone on x, treatment dodge, faceted)
#   5. Time moving × timepoint          (treatment on x, faceted)
#   6. Time moving — aggregated         (treatment on x)
#
# FEATURES:
#   - All visual style parameters (colours, shapes, sizes, fonts, margins)
#   - CLD letter annotations: toggle, size, vertical nudge controls
#   - ANOVA caption: toggle, per-plot editable text (auto-loaded from Step 5)
#   - Timepoint strip labels: editable in app
#   - Code export → paste constants into activity_analysis_STATS_choice_exp.R
#   - Download PNG
#
# HOW TO USE:
#   1. Run Steps 1-2 (and optionally Step 5) so fish_activity_summary and
#      res_* result objects are in GlobalEnv.  OR source from the pipeline
#      working directory — the app will search STEP2_output/ and
#      STEP5_stats/ automatically.
#   2. Click Source in RStudio.
#   3. Adjust controls in the sidebar — preview updates live.
#   4. Edit ANOVA captions in the "ANOVA captions" tab (pre-filled from
#      Step 5 results; leave blank to suppress the caption).
#   5. Copy the generated code from "Plot code" → paste into STATS script.
#   6. Use the Download button to save the current plot as a PNG.
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

.sem <- function(x) { x <- x[is.finite(x)]; if (length(x) < 2) NA_real_ else sd(x)/sqrt(length(x)) }

.find_latest_dir <- function(parent, prefix) {
  if (!dir.exists(parent)) return(NULL)
  subs <- list.dirs(parent, full.names = TRUE, recursive = FALSE)
  subs <- subs[grepl(paste0("^", prefix, "_\\d{8}_\\d{6}$"), basename(subs))]
  if (length(subs) == 0) return(NULL)
  subs[which.max(file.mtime(subs))]
}

.find_latest_csv <- function(step_name, csv_filename) {
  d <- .find_latest_dir(file.path(getwd(), step_name), step_name)
  if (is.null(d)) return(NULL)
  f <- file.path(d, csv_filename)
  if (file.exists(f)) f else NULL
}

.get_global <- function(name, default = NULL) {
  if (exists(name, envir = .GlobalEnv, inherits = FALSE))
    get(name, envir = .GlobalEnv, inherits = FALSE)
  else default
}


# =============================================================================
# ==== LOAD FISH ACTIVITY SUMMARY =============================================
# =============================================================================

if (exists("fish_activity_summary", envir = .GlobalEnv, inherits = FALSE)) {
  .fas <- get("fish_activity_summary", envir = .GlobalEnv, inherits = FALSE)
  message("Using fish_activity_summary from GlobalEnv (", nrow(.fas), " rows)")
} else {
  message("fish_activity_summary not in GlobalEnv; searching STEP2_output/...")
  .csv_path <- .find_latest_csv("STEP2_output", "fish_activity_summary.csv")
  if (is.null(.csv_path))
    stop("Cannot find fish_activity_summary.csv. Run Steps 1-2 first.", call. = FALSE)
  .fas <- readr::read_csv(.csv_path, show_col_types = FALSE)
  assign("fish_activity_summary", .fas, envir = .GlobalEnv)
  message("Loaded: ", .csv_path)
}

TREATMENT_LEVELS_g <- .get_global("TREATMENT_LEVELS", c("control", "exercise choice"))
TIMEPOINT_LEVELS_g <- .get_global("TIMEPOINT_LEVELS", 1:3)

.fas <- .fas %>%
  dplyr::mutate(
    treatment   = factor(trimws(tolower(as.character(treatment))),
                         levels = TREATMENT_LEVELS_g),
    timepoint   = suppressWarnings(as.integer(timepoint)),
    timepoint_f = factor(timepoint, levels = TIMEPOINT_LEVELS_g),
    fish_uid    = paste(trial_id, fish_id, sep = "_")
  )

.has_tp       <- any(!is.na(.fas$timepoint))
.has_switches <- all(c("n_flow_entries", "n_calm_entries") %in% names(.fas))
.has_moving   <- "moving_time_TS1_s" %in% names(.fas)
.has_zones    <- all(c("time_in_flow_s", "time_in_calm_s") %in% names(.fas))

# Pre-build long-format datasets (rebuilt inside reactive only if needed)
if (.has_zones) {
  .df_main_long <- .fas %>%
    dplyr::select(fish_uid, trial_id, tank, treatment, timepoint, timepoint_f,
                  time_in_flow_s, time_in_calm_s) %>%
    tidyr::pivot_longer(c(time_in_flow_s, time_in_calm_s),
                        names_to = "zone", values_to = "occupancy_s") %>%
    dplyr::mutate(zone = factor(
      dplyr::recode(zone, time_in_flow_s = "flow", time_in_calm_s = "calm"),
      levels = c("flow", "calm")
    ))
} else {
  .df_main_long <- NULL
}

if (.has_switches) {
  .df_switch_long <- .fas %>%
    dplyr::select(fish_uid, trial_id, tank, treatment, timepoint, timepoint_f,
                  n_flow_entries, n_calm_entries) %>%
    tidyr::pivot_longer(c(n_flow_entries, n_calm_entries),
                        names_to = "zone", values_to = "n_entries") %>%
    dplyr::mutate(zone = factor(
      dplyr::recode(zone, n_flow_entries = "flow", n_calm_entries = "calm"),
      levels = c("flow", "calm")
    ))
} else {
  .df_switch_long <- NULL
}


# =============================================================================
# ==== LOAD STEP 5 RESULTS (CLD + ANOVA CAPTIONS) =============================
# =============================================================================

.step5_latest <- .find_latest_dir(file.path(getwd(), "STEP5_stats"), "STEP5_stats")
if (!is.null(.step5_latest))
  message("Found STEP5 output: ", basename(.step5_latest))

# Get result object from GlobalEnv or NULL
.get_res <- function(name) .get_global(name, NULL)

.res <- list(
  zone_agg   = .get_res("res_zone_agg"),
  zone_tp    = .get_res("res_zone_tp"),
  switch_agg = .get_res("res_switch_agg"),
  switch_tp  = .get_res("res_switch_tp"),
  moving_tp  = .get_res("res_moving_tp"),
  moving_agg = .get_res("res_moving_agg")
)

# ---- ANOVA caption utilities ------------------------------------------------

.anova_cap_from_csv <- function(anova_csv, term_pattern) {
  if (is.null(anova_csv) || !file.exists(anova_csv)) return("")
  av <- tryCatch(readr::read_csv(anova_csv, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(av)) return("")
  # normalise column names (in-memory vs disk may differ)
  names(av) <- dplyr::recode(names(av),
    Chisq = "chisq", Df = "df", `Pr(>Chisq)` = "p_value")
  row <- av[grepl(term_pattern, av$term, ignore.case = TRUE, perl = TRUE), ]
  if (nrow(row) == 0) return("")
  row <- row[1, ]
  p <- row$p_value
  if (is.na(p)) return("")
  p_str <- if (p < 0.001) "p < 0.001" else sprintf("p = %.3f", p)
  sprintf("\u03c7\u00b2(%d) = %.2f, %s", as.integer(row$df), row$chisq, p_str)
}

.get_anova_cap <- function(res, res_key, step5_label, term_pattern) {
  # GlobalEnv result object first
  if (!is.null(res) && !is.null(res$anova_caps[[res_key]])) return(res$anova_caps[[res_key]])
  # Fall back to STEP5 CSV on disk
  if (!is.null(.step5_latest)) {
    f <- file.path(.step5_latest, step5_label, "anova.csv")
    cap <- .anova_cap_from_csv(f, term_pattern)
    if (nchar(cap) > 0) return(cap)
  }
  ""
}

# ---- CLD utilities ----------------------------------------------------------

.get_cld <- function(res, term_re, step5_label) {
  # GlobalEnv result object first
  if (!is.null(res) && !is.null(res$posthoc)) {
    ph_match <- res$posthoc[grepl(term_re, names(res$posthoc), ignore.case = TRUE, perl = TRUE)]
    if (length(ph_match) > 0 && !is.null(ph_match[[1]]$cld)) return(ph_match[[1]]$cld)
  }
  # Fall back to disk: load all cld_*.csv and return the most complex one
  if (!is.null(.step5_latest)) {
    cld_dir <- file.path(.step5_latest, step5_label)
    cld_files <- list.files(cld_dir, pattern = "^cld_.*\\.csv$", full.names = TRUE)
    if (length(cld_files) > 0) {
      # Prefer files whose name matches term_re
      preferred <- cld_files[grepl(term_re, sub("\\.csv$","",basename(cld_files)),
                                   ignore.case = TRUE, perl = TRUE)]
      f <- if (length(preferred) > 0) preferred[1] else cld_files[length(cld_files)]
      return(tryCatch(readr::read_csv(f, show_col_types = FALSE), error = function(e) NULL))
    }
  }
  NULL
}

# ---- Pre-load captions and CLD per plot type --------------------------------

DEFAULT_CAPTIONS <- list(
  zone_agg   = .get_anova_cap(.res$zone_agg,   "Treatment:Zone",
                              "zone_aggregated",   "treatment.*zone|zone.*treatment"),
  zone_tp    = .get_anova_cap(.res$zone_tp,    "Treatment:Zone:Timepoint",
                              "zone_timepoint",    "treatment.*zone.*timepoint"),
  switch_agg = .get_anova_cap(.res$switch_agg, "Treatment:Zone",
                              "switches_aggregated","treatment.*zone|zone.*treatment"),
  switch_tp  = .get_anova_cap(.res$switch_tp,  "Treatment:Zone:Timepoint",
                              "switches_timepoint","treatment.*zone.*timepoint"),
  moving_tp  = .get_anova_cap(.res$moving_tp,  "Treatment:Timepoint",
                              "moving_timepoint",  "treatment.*timepoint|timepoint.*treatment"),
  moving_agg = .get_anova_cap(.res$moving_agg, "Treatment",
                              "moving_aggregated", "^treatment$")
)

DEFAULT_CLD <- list(
  zone_agg   = .get_cld(.res$zone_agg,   "treatment.*zone|zone.*treatment",          "zone_aggregated"),
  zone_tp    = .get_cld(.res$zone_tp,    "treatment.*zone.*timepoint",               "zone_timepoint"),
  switch_agg = .get_cld(.res$switch_agg, "treatment.*zone|zone.*treatment",          "switches_aggregated"),
  switch_tp  = .get_cld(.res$switch_tp,  "treatment.*zone.*timepoint",               "switches_timepoint"),
  moving_tp  = .get_cld(.res$moving_tp,  "treatment.*timepoint|timepoint.*treatment","moving_timepoint"),
  moving_agg = .get_cld(.res$moving_agg, "^treatment$",                              "moving_aggregated")
)

.cld_available  <- vapply(DEFAULT_CLD,      Negate(is.null), logical(1))
.caps_available <- vapply(DEFAULT_CAPTIONS, function(x) nchar(x) > 0, logical(1))

if (any(.cld_available))
  message("CLD loaded for: ",  paste(names(.cld_available)[.cld_available],  collapse = ", "))
if (any(.caps_available))
  message("Captions loaded for: ", paste(names(.caps_available)[.caps_available], collapse = ", "))


# =============================================================================
# ==== STYLE DEFAULTS =========================================================
# =============================================================================

STYLE_DEFAULTS <- list(
  plot_type            = "zone_agg",
  # points
  point_size           = 2.4,
  point_alpha          = 0.80,
  jitter_width         = 0.12,
  dodge_width          = 0.70,
  # mean marker
  errorbar_width       = 0.13,
  errorbar_lwd         = 0.65,
  crossbar_width       = 0.40,
  crossbar_lwd         = 0.85,
  # treatment colours
  col_control          = "#2166AC",
  col_exercise         = "#D6604D",
  # zone colours
  col_zone_flow        = "#4DAC26",
  col_zone_calm        = "#8073AC",
  # typography
  base_size            = 13,
  axis_text_x_size     = 13,
  axis_text_y_size     = 11,
  axis_title_size      = 13,
  strip_text_size      = 12,
  caption_size         = 10,
  legend_text_size     = 10,
  legend_title_size    = 11,
  # legend & grid
  show_legend          = TRUE,
  legend_position      = "right",
  show_major_grid      = FALSE,
  show_minor_grid      = FALSE,
  # axes & margins
  axis_lwd             = 0.85,
  margin_t             = 8,
  margin_r             = 10,
  margin_b             = 16,
  margin_l             = 10,
  # CLD letters
  show_cld             = TRUE,
  cld_size             = 4.0,
  cld_nudge_sem_mult   = 1.6,
  cld_nudge_range_frac = 0.05,
  # ANOVA caption
  show_caption         = TRUE,
  # timepoint strip labels
  tp_label_1           = "Segment 1  (0\u201320 min)",
  tp_label_2           = "Segment 2  (20\u201340 min)",
  tp_label_3           = "Segment 3  (40\u201360 min)",
  # export size
  plot_width           = 10,
  plot_height          = 6
)


# =============================================================================
# ==== PLOT BUILDER ===========================================================
# =============================================================================

.empty_plot <- function(msg) {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0.5, y = 0.5, label = msg, size = 5,
                      hjust = 0.5, vjust = 0.5) +
    ggplot2::theme_void()
}

.build_base_theme <- function(input) {
  ggplot2::theme_minimal(base_size = input$base_size) +
    ggplot2::theme(
      panel.grid.major = if (isTRUE(input$show_major_grid))
        ggplot2::element_line() else ggplot2::element_blank(),
      panel.grid.minor = if (isTRUE(input$show_minor_grid))
        ggplot2::element_line() else ggplot2::element_blank(),
      panel.background  = ggplot2::element_blank(),
      axis.line.x       = ggplot2::element_line(color = "black", linewidth = input$axis_lwd),
      axis.line.y       = ggplot2::element_line(color = "black", linewidth = input$axis_lwd),
      axis.ticks        = ggplot2::element_line(color = "black",
                                                linewidth = input$axis_lwd * 0.8),
      axis.text.x       = ggplot2::element_text(size  = input$axis_text_x_size,
                                                face  = "bold", color = "black"),
      axis.text.y       = ggplot2::element_text(size  = input$axis_text_y_size,
                                                color = "black"),
      axis.title        = ggplot2::element_text(size  = input$axis_title_size,
                                                face  = "bold", color = "black"),
      legend.title      = ggplot2::element_text(size  = input$legend_title_size, face = "bold"),
      legend.text       = ggplot2::element_text(size  = input$legend_text_size),
      strip.text        = ggplot2::element_text(size  = input$strip_text_size, face = "bold"),
      strip.background  = ggplot2::element_blank(),
      legend.position   = if (isTRUE(input$show_legend)) input$legend_position else "none",
      plot.margin       = ggplot2::margin(t = input$margin_t, r = input$margin_r,
                                          b = input$margin_b, l = input$margin_l),
      plot.caption      = ggplot2::element_text(size = input$caption_size,
                                                hjust = 0, face = "italic")
    )
}

.add_cld_zone <- function(p, cld_df, smry, POS_DODGE, input) {
  if (is.null(cld_df) || !".group" %in% names(cld_df)) return(p)
  merge_cols <- intersect(c("zone", "treatment", "timepoint_f"), names(cld_df))
  if (length(merge_cols) == 0) return(p)
  cld_df$.group <- trimws(as.character(cld_df$.group))
  cld_join <- tryCatch(
    dplyr::left_join(smry,
                     dplyr::select(cld_df, dplyr::all_of(c(merge_cols, ".group"))),
                     by = merge_cols),
    error = function(e) NULL
  )
  if (is.null(cld_join) || !".group" %in% names(cld_join)) return(p)
  y_range  <- diff(range(smry$mean_y, na.rm = TRUE))
  cld_join$label_y <- cld_join$mean_y +
    cld_join$sem_y * input$cld_nudge_sem_mult +
    y_range * input$cld_nudge_range_frac
  p + ggplot2::geom_text(
    data = cld_join,
    ggplot2::aes(x = zone, y = label_y, label = .group,
                 group = treatment, colour = zone),
    position  = POS_DODGE,
    size      = input$cld_size,
    fontface  = "bold",
    show.legend = FALSE
  )
}

.add_cld_std <- function(p, cld_df, smry, input) {
  if (is.null(cld_df) || !".group" %in% names(cld_df)) return(p)
  merge_cols <- intersect(c("treatment", "timepoint_f"), names(cld_df))
  cld_df$.group <- trimws(as.character(cld_df$.group))
  cld_join <- tryCatch(
    dplyr::left_join(smry,
                     dplyr::select(cld_df, dplyr::all_of(c(merge_cols, ".group"))),
                     by = merge_cols),
    error = function(e) NULL
  )
  if (is.null(cld_join) || !".group" %in% names(cld_join)) return(p)
  y_range  <- diff(range(smry$mean_y, na.rm = TRUE))
  cld_join$label_y <- cld_join$mean_y +
    cld_join$sem_y * input$cld_nudge_sem_mult +
    y_range * input$cld_nudge_range_frac
  p + ggplot2::geom_text(
    data = cld_join,
    ggplot2::aes(x = treatment, y = label_y, label = .group, colour = treatment),
    size = input$cld_size, fontface = "bold", show.legend = FALSE
  )
}

make_skinny_plot <- function(input) {

  pt     <- input$plot_type
  is_tp  <- grepl("_tp$", pt)
  is_std <- grepl("^moving", pt)

  trt_colors <- c("control" = input$col_control, "exercise choice" = input$col_exercise)
  trt_shapes <- c("control" = 17L, "exercise choice" = 16L)
  zn_colors  <- c("flow" = input$col_zone_flow, "calm" = input$col_zone_calm)

  BASE_THEME <- .build_base_theme(input)

  tp_labels <- stats::setNames(
    c(input$tp_label_1, input$tp_label_2, input$tp_label_3),
    as.character(TIMEPOINT_LEVELS_g)
  )

  # Caption for this plot type
  cap_id  <- paste0("cap_", pt)
  cap_val <- input[[cap_id]]
  caption_txt <- if (isTRUE(input$show_caption) && !is.null(cap_val) &&
                     nchar(trimws(cap_val)) > 0)
    trimws(cap_val) else NULL

  # CLD data
  cld_df <- if (isTRUE(input$show_cld)) DEFAULT_CLD[[pt]] else NULL

  # Positions
  POS_JD    <- ggplot2::position_jitterdodge(jitter.width  = input$jitter_width,
                                              jitter.height = 0,
                                              dodge.width   = input$dodge_width)
  POS_DODGE <- ggplot2::position_dodge(width = input$dodge_width)
  POS_JIT   <- ggplot2::position_jitter(width = input$jitter_width, height = 0)

  # ---------- Zone plots (zone_agg / zone_tp / switch_agg / switch_tp) ---------
  if (!is_std) {
    if (grepl("^zone", pt)) {
      if (is.null(.df_main_long))
        return(.empty_plot("Zone occupancy columns not found.\nRun Steps 1-2 first."))
      d     <- .df_main_long
      y_col <- "occupancy_s"
      y_lbl <- "Time in zone (s)"
      title <- if (is_tp) "Zone occupancy by treatment \u00d7 timepoint"
               else "Zone occupancy by treatment"
    } else {
      if (is.null(.df_switch_long))
        return(.empty_plot("Zone entry counts not found.\nRe-run Step 2 (needs n_flow_entries)."))
      d     <- .df_switch_long
      y_col <- "n_entries"
      y_lbl <- "Zone entries (n)"
      title <- if (is_tp) "Zone entries by treatment \u00d7 timepoint"
               else "Zone entries by treatment"
    }

    grps <- if (is_tp) c("zone","treatment","timepoint_f") else c("zone","treatment")
    smry <- d %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(grps))) %>%
      dplyr::summarise(mean_y = mean(.data[[y_col]], na.rm = TRUE),
                       sem_y  = .sem(.data[[y_col]]), .groups = "drop")

    p <- ggplot2::ggplot(d, ggplot2::aes(x = zone, y = .data[[y_col]],
                                          colour = zone, shape = treatment, group = treatment)) +
      ggplot2::geom_point(position = POS_JD,
                          alpha = input$point_alpha, size = input$point_size) +
      ggplot2::geom_crossbar(
        data = smry,
        ggplot2::aes(y = mean_y, ymin = mean_y, ymax = mean_y,
                     colour = zone, group = treatment),
        position = POS_DODGE, width = input$crossbar_width,
        linewidth = input$crossbar_lwd, show.legend = FALSE) +
      ggplot2::geom_errorbar(
        data = smry,
        ggplot2::aes(y = mean_y, ymin = mean_y - sem_y, ymax = mean_y + sem_y,
                     colour = zone, group = treatment),
        position = POS_DODGE, width = input$errorbar_width,
        linewidth = input$errorbar_lwd, show.legend = FALSE) +
      ggplot2::scale_colour_manual(values = zn_colors,  name = "Zone") +
      ggplot2::scale_shape_manual(values  = trt_shapes, name = "Treatment") +
      BASE_THEME +
      ggplot2::labs(title = title, x = "Zone", y = y_lbl, caption = caption_txt)

    p <- .add_cld_zone(p, cld_df, smry, POS_DODGE, input)
    if (is_tp)
      p <- p + ggplot2::facet_wrap(~timepoint_f, ncol = 1,
                                   labeller = ggplot2::labeller(timepoint_f = tp_labels))
    return(p)
  }

  # ---------- Standard plots (moving_tp / moving_agg) -------------------------
  if (!.has_moving)
    return(.empty_plot("moving_time_TS1_s not found in fish_activity_summary."))

  d <- .fas[is.finite(.fas$moving_time_TS1_s), ]

  grps <- if (is_tp) c("treatment","timepoint_f") else "treatment"
  smry <- d %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(grps))) %>%
    dplyr::summarise(mean_y = mean(moving_time_TS1_s, na.rm = TRUE),
                     sem_y  = .sem(moving_time_TS1_s), .groups = "drop")

  p <- ggplot2::ggplot(d, ggplot2::aes(x = treatment, y = moving_time_TS1_s,
                                        colour = treatment, shape = treatment)) +
    ggplot2::geom_point(position = POS_JIT,
                        alpha = input$point_alpha, size = input$point_size) +
    ggplot2::geom_crossbar(
      data = smry,
      ggplot2::aes(y = mean_y, ymin = mean_y, ymax = mean_y),
      width = input$crossbar_width, linewidth = input$crossbar_lwd,
      show.legend = FALSE) +
    ggplot2::geom_errorbar(
      data = smry,
      ggplot2::aes(y = mean_y, ymin = mean_y - sem_y, ymax = mean_y + sem_y),
      width = input$errorbar_width, linewidth = input$errorbar_lwd,
      show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = trt_colors, name = "Treatment") +
    ggplot2::scale_shape_manual(values  = trt_shapes, name = "Treatment") +
    ggplot2::scale_x_discrete(labels = stringr::str_to_title) +
    BASE_THEME +
    ggplot2::labs(
      title   = if (is_tp) "Time moving by treatment \u00d7 timepoint"
                else "Time moving by treatment",
      x = "Treatment", y = "Time moving (s)", caption = caption_txt
    )

  p <- .add_cld_std(p, cld_df, smry, input)
  if (is_tp)
    p <- p + ggplot2::facet_wrap(~timepoint_f, ncol = 1,
                                 labeller = ggplot2::labeller(timepoint_f = tp_labels))
  p
}


# =============================================================================
# ==== CODE EXPORT ============================================================
# =============================================================================

build_plot_code <- function(input) {
  grid_str <- if (!isTRUE(input$show_major_grid) && !isTRUE(input$show_minor_grid))
    "ggplot2::element_blank()"
  else if (isTRUE(input$show_major_grid) && !isTRUE(input$show_minor_grid))
    "ggplot2::element_line()  # major only; set minor separately"
  else
    "ggplot2::element_line()"

  paste0(
'# ---- Paste into activity_analysis_STATS_choice_exp.R ----
# (Section 6 — SKINNY GRAPHS, "Plot settings" block)

TREATMENT_COLORS <- c("control"         = "', input$col_control,  '",
                       "exercise choice" = "', input$col_exercise, '")
TREATMENT_SHAPES <- c("control" = 17L, "exercise choice" = 16L)
ZONE_COLORS_MAIN <- c("flow" = "', input$col_zone_flow, '",
                       "calm" = "', input$col_zone_calm, '")

JITTER_W             <- ', input$jitter_width, '
DODGE_W              <- ', input$dodge_width, '
PT_ALPHA             <- ', input$point_alpha, '
PT_SIZE              <- ', input$point_size, '
MEAN_W               <- ', input$crossbar_width, '
LW_MEAN              <- ', input$crossbar_lwd, '
LW_ERR               <- ', input$errorbar_lwd, '
ERR_W                <- ', input$errorbar_width, '
CLD_SIZE             <- ', input$cld_size, '
CLD_NUDGE_SEM_MULT   <- ', input$cld_nudge_sem_mult, '
CLD_NUDGE_RANGE_FRAC <- ', input$cld_nudge_range_frac, '

TIMEPOINT_LABELS <- c(
  "1" = "', input$tp_label_1, '",
  "2" = "', input$tp_label_2, '",
  "3" = "', input$tp_label_3, '"
)

BASE_THEME <- ggplot2::theme_minimal(base_size = ', input$base_size, ') +
  ggplot2::theme(
    panel.grid        = ', grid_str, ',
    panel.background  = ggplot2::element_blank(),
    axis.line.x       = ggplot2::element_line(color = "black", linewidth = ', input$axis_lwd, '),
    axis.line.y       = ggplot2::element_line(color = "black", linewidth = ', input$axis_lwd, '),
    axis.ticks        = ggplot2::element_line(color = "black", linewidth = ', round(input$axis_lwd * 0.8, 3), '),
    axis.text.x       = ggplot2::element_text(size = ', input$axis_text_x_size, ', face = "bold", color = "black"),
    axis.text.y       = ggplot2::element_text(size = ', input$axis_text_y_size, ', color = "black"),
    axis.title        = ggplot2::element_text(size = ', input$axis_title_size, ', face = "bold", color = "black"),
    legend.title      = ggplot2::element_text(size = ', input$legend_title_size, ', face = "bold"),
    legend.text       = ggplot2::element_text(size = ', input$legend_text_size, '),
    strip.text        = ggplot2::element_text(size = ', input$strip_text_size, ', face = "bold"),
    strip.background  = ggplot2::element_blank(),
    plot.margin       = ggplot2::margin(', input$margin_t, ', ', input$margin_r, ', ', input$margin_b, ', ', input$margin_l, '),
    plot.caption      = ggplot2::element_text(size = ', input$caption_size, ', hjust = 0, face = "italic")
  )
'
  )
}

build_data_summary <- function() {
  tbl <- .fas %>%
    dplyr::group_by(treatment) %>%
    dplyr::summarise(
      n_fish      = dplyr::n(),
      n_timepoints = dplyr::n_distinct(timepoint[!is.na(timepoint)]),
      n_tanks     = dplyr::n_distinct(tank),
      .groups     = "drop"
    )
  if (.has_moving)
    tbl <- tbl %>% dplyr::left_join(
      .fas %>% dplyr::group_by(treatment) %>%
        dplyr::summarise(mean_moving_s = round(mean(moving_time_TS1_s, na.rm=TRUE), 1),
                         .groups="drop"),
      by = "treatment")
  tbl
}


# =============================================================================
# ==== UI =====================================================================
# =============================================================================

PLOT_TYPE_CHOICES <- c(
  "Zone occupancy — aggregated"  = "zone_agg",
  "Zone occupancy \u00d7 timepoint" = "zone_tp",
  "Zone entries — aggregated"    = "switch_agg",
  "Zone entries \u00d7 timepoint"   = "switch_tp",
  "Time moving \u00d7 timepoint"    = "moving_tp",
  "Time moving — aggregated"     = "moving_agg"
)

ui <- shiny::fluidPage(
  shiny::titlePanel("Choice experiment — skinny graphs style tuner"),
  shiny::tags$head(shiny::tags$style(shiny::HTML("
    .ctrl-scroll {
      max-height: calc(100vh - 120px);
      overflow-y: auto;
      padding-right: 6px;
    }
    .preview-sticky {
      position: sticky;
      top: 10px;
    }
    hr { margin: 8px 0; }
    .note-text { font-size: 11px; color: #777; margin-top: 2px; }
  "))),

  shiny::sidebarLayout(
    shiny::sidebarPanel(
      width = 3,
      shiny::div(class = "ctrl-scroll",

        shiny::h5("Plot type", style = "font-weight:bold; margin-top:0;"),
        shiny::selectInput("plot_type", NULL,
                           choices  = PLOT_TYPE_CHOICES,
                           selected = STYLE_DEFAULTS$plot_type),

        shiny::hr(),
        shiny::h5("Points", style = "font-weight:bold;"),
        shiny::sliderInput("point_size",   "Point size",
                           min=0.5, max=7, value=STYLE_DEFAULTS$point_size, step=0.1),
        shiny::sliderInput("point_alpha",  "Opacity",
                           min=0.1, max=1, value=STYLE_DEFAULTS$point_alpha, step=0.05),
        shiny::sliderInput("jitter_width", "Jitter width",
                           min=0, max=0.5, value=STYLE_DEFAULTS$jitter_width, step=0.01),
        shiny::conditionalPanel(
          condition = "['zone_agg','zone_tp','switch_agg','switch_tp'].indexOf(input.plot_type) >= 0",
          shiny::sliderInput("dodge_width", "Dodge width",
                             min=0.2, max=1.2, value=STYLE_DEFAULTS$dodge_width, step=0.05)
        ),

        shiny::hr(),
        shiny::h5("Mean marker", style = "font-weight:bold;"),
        shiny::sliderInput("crossbar_width", "Crossbar width",
                           min=0.05, max=1, value=STYLE_DEFAULTS$crossbar_width, step=0.05),
        shiny::sliderInput("crossbar_lwd", "Crossbar line width",
                           min=0.1, max=3, value=STYLE_DEFAULTS$crossbar_lwd, step=0.05),
        shiny::sliderInput("errorbar_width", "Errorbar width",
                           min=0.01, max=0.5, value=STYLE_DEFAULTS$errorbar_width, step=0.01),
        shiny::sliderInput("errorbar_lwd", "Errorbar line width",
                           min=0.1, max=3, value=STYLE_DEFAULTS$errorbar_lwd, step=0.05),

        shiny::hr(),
        shiny::h5("Treatment colours", style = "font-weight:bold;"),
        colourpicker::colourInput("col_control",  "Control",
                                  value = STYLE_DEFAULTS$col_control),
        colourpicker::colourInput("col_exercise", "Exercise choice",
                                  value = STYLE_DEFAULTS$col_exercise),

        shiny::hr(),
        shiny::h5("Zone colours", style = "font-weight:bold;"),
        colourpicker::colourInput("col_zone_flow", "Flow",
                                  value = STYLE_DEFAULTS$col_zone_flow),
        colourpicker::colourInput("col_zone_calm", "Calm",
                                  value = STYLE_DEFAULTS$col_zone_calm),

        shiny::hr(),
        shiny::h5("Typography", style = "font-weight:bold;"),
        shiny::sliderInput("base_size",         "Base size",
                           min=8, max=24, value=STYLE_DEFAULTS$base_size, step=1),
        shiny::sliderInput("axis_text_x_size",  "X axis text",
                           min=6, max=24, value=STYLE_DEFAULTS$axis_text_x_size, step=1),
        shiny::sliderInput("axis_text_y_size",  "Y axis text",
                           min=6, max=24, value=STYLE_DEFAULTS$axis_text_y_size, step=1),
        shiny::sliderInput("axis_title_size",   "Axis title",
                           min=6, max=28, value=STYLE_DEFAULTS$axis_title_size, step=1),
        shiny::sliderInput("strip_text_size",   "Facet strip text",
                           min=6, max=24, value=STYLE_DEFAULTS$strip_text_size, step=1),
        shiny::sliderInput("caption_size",      "Caption size",
                           min=6, max=18, value=STYLE_DEFAULTS$caption_size, step=1),
        shiny::sliderInput("legend_text_size",  "Legend text",
                           min=6, max=24, value=STYLE_DEFAULTS$legend_text_size, step=1),
        shiny::sliderInput("legend_title_size", "Legend title",
                           min=6, max=24, value=STYLE_DEFAULTS$legend_title_size, step=1),

        shiny::hr(),
        shiny::h5("Legend & grid", style = "font-weight:bold;"),
        shiny::checkboxInput("show_legend",     "Show legend",      STYLE_DEFAULTS$show_legend),
        shiny::selectInput("legend_position", "Legend position",
                           choices  = c("right","left","top","bottom"),
                           selected = STYLE_DEFAULTS$legend_position),
        shiny::checkboxInput("show_major_grid", "Major grid lines", STYLE_DEFAULTS$show_major_grid),
        shiny::checkboxInput("show_minor_grid", "Minor grid lines", STYLE_DEFAULTS$show_minor_grid),

        shiny::hr(),
        shiny::h5("Axes & margins", style = "font-weight:bold;"),
        shiny::sliderInput("axis_lwd",  "Axis line width",
                           min=0.1, max=2.5, value=STYLE_DEFAULTS$axis_lwd, step=0.05),
        shiny::sliderInput("margin_t", "Top margin",
                           min=0, max=50, value=STYLE_DEFAULTS$margin_t, step=1),
        shiny::sliderInput("margin_r", "Right margin",
                           min=0, max=50, value=STYLE_DEFAULTS$margin_r, step=1),
        shiny::sliderInput("margin_b", "Bottom margin",
                           min=0, max=50, value=STYLE_DEFAULTS$margin_b, step=1),
        shiny::sliderInput("margin_l", "Left margin",
                           min=0, max=50, value=STYLE_DEFAULTS$margin_l, step=1),

        shiny::hr(),
        shiny::h5("CLD letters", style = "font-weight:bold;"),
        shiny::p(class = "note-text",
          if (any(.cld_available))
            paste("Available for:", paste(names(.cld_available)[.cld_available], collapse=", "))
          else
            "No CLD data found. Run Step 5 first."
        ),
        shiny::checkboxInput("show_cld", "Show CLD letters", STYLE_DEFAULTS$show_cld),
        shiny::conditionalPanel(
          condition = "input.show_cld",
          shiny::sliderInput("cld_size", "Letter size",
                             min=2, max=8, value=STYLE_DEFAULTS$cld_size, step=0.25),
          shiny::sliderInput("cld_nudge_sem_mult", "Vertical nudge (× SEM)",
                             min=0, max=5, value=STYLE_DEFAULTS$cld_nudge_sem_mult, step=0.1),
          shiny::sliderInput("cld_nudge_range_frac", "Nudge (fraction of y range)",
                             min=0, max=0.3, value=STYLE_DEFAULTS$cld_nudge_range_frac, step=0.01)
        ),

        shiny::hr(),
        shiny::h5("ANOVA caption", style = "font-weight:bold;"),
        shiny::checkboxInput("show_caption", "Show caption", STYLE_DEFAULTS$show_caption),
        shiny::p(class = "note-text", "Edit captions in the 'ANOVA captions' tab."),

        shiny::hr(),
        shiny::h5("Export size (inches)", style = "font-weight:bold;"),
        shiny::numericInput("plot_width",  "Width",
                            value=STYLE_DEFAULTS$plot_width,  min=2, max=24, step=0.5),
        shiny::numericInput("plot_height", "Height",
                            value=STYLE_DEFAULTS$plot_height, min=2, max=24, step=0.5),
        shiny::downloadButton("download_plot", "Download PNG",
                              style = "width:100%;")
      )
    ),

    shiny::mainPanel(
      width = 9,
      shiny::div(class = "preview-sticky",
        shiny::tabsetPanel(

          shiny::tabPanel("Preview",
            shiny::plotOutput("plot_preview", height = "560px")
          ),

          shiny::tabPanel("ANOVA captions",
            shiny::tags$p(
              style = "margin-top:10px; color:#555;",
              "One caption per plot type. Pre-populated from Step 5 results if available.",
              "Leave blank to suppress. Unicode characters (χ², ×) are supported."
            ),
            shiny::fluidRow(
              shiny::column(6,
                shiny::textAreaInput("cap_zone_agg",   "Zone occupancy — aggregated",
                                     value = DEFAULT_CAPTIONS$zone_agg,   rows = 2),
                shiny::textAreaInput("cap_switch_agg", "Zone entries — aggregated",
                                     value = DEFAULT_CAPTIONS$switch_agg, rows = 2),
                shiny::textAreaInput("cap_moving_tp",  "Time moving \u00d7 timepoint",
                                     value = DEFAULT_CAPTIONS$moving_tp,  rows = 2)
              ),
              shiny::column(6,
                shiny::textAreaInput("cap_zone_tp",    "Zone occupancy \u00d7 timepoint",
                                     value = DEFAULT_CAPTIONS$zone_tp,    rows = 2),
                shiny::textAreaInput("cap_switch_tp",  "Zone entries \u00d7 timepoint",
                                     value = DEFAULT_CAPTIONS$switch_tp,  rows = 2),
                shiny::textAreaInput("cap_moving_agg", "Time moving — aggregated",
                                     value = DEFAULT_CAPTIONS$moving_agg, rows = 2)
              )
            ),
            shiny::hr(),
            shiny::h5("Timepoint strip labels", style = "font-weight:bold;"),
            shiny::fluidRow(
              shiny::column(4,
                shiny::textInput("tp_label_1", "Segment 1",
                                  value = STYLE_DEFAULTS$tp_label_1)),
              shiny::column(4,
                shiny::textInput("tp_label_2", "Segment 2",
                                  value = STYLE_DEFAULTS$tp_label_2)),
              shiny::column(4,
                shiny::textInput("tp_label_3", "Segment 3",
                                  value = STYLE_DEFAULTS$tp_label_3))
            )
          ),

          shiny::tabPanel("Plot code",
            shiny::tags$p(
              style = "margin-top:10px; color:#555;",
              "Copy-paste into activity_analysis_STATS_choice_exp.R",
              "(Section 6 — Plot settings block)."
            ),
            shiny::tags$pre(
              style = "background:#f5f5f5; padding:12px; font-size:12px; overflow-x:auto; white-space:pre-wrap;",
              shiny::textOutput("plot_code")
            )
          ),

          shiny::tabPanel("Data summary",
            shiny::tags$p(style = "margin-top:10px; color:#555;",
                          "Fish counts and data availability."),
            shiny::tableOutput("data_summary"),
            shiny::hr(),
            shiny::tags$p(style="color:#555;",
              paste0(
                if (.has_zones) "\u2713 Zone occupancy data available.  " else "\u2717 Zone data missing.  ",
                if (.has_switches) "\u2713 Switch counts available.  " else "\u2717 Switch counts missing (re-run Step 2).  ",
                if (.has_moving) "\u2713 Moving time available.  " else "\u2717 Moving time missing.  ",
                if (.has_tp) paste0("\u2713 Timepoints: ", paste(TIMEPOINT_LEVELS_g, collapse=", "), ".") else "\u2717 No timepoint data."
              )
            ),
            shiny::hr(),
            shiny::tags$p(style="color:#555;",
              paste0(
                "CLD loaded: ",
                if (any(.cld_available)) paste(names(.cld_available)[.cld_available], collapse=", ") else "none",
                "  |  ANOVA captions loaded: ",
                if (any(.caps_available)) paste(names(.caps_available)[.caps_available], collapse=", ") else "none"
              )
            )
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
    shiny::req(input$plot_type)
    # Ensure all caption inputs are initialized before building
    shiny::req(input$cap_zone_agg, input$cap_zone_tp,
               input$cap_switch_agg, input$cap_switch_tp,
               input$cap_moving_tp, input$cap_moving_agg,
               cancelOutput = TRUE)
    tryCatch(
      make_skinny_plot(input),
      error = function(e) {
        .empty_plot(paste0("Plot error:\n", conditionMessage(e)))
      }
    )
  })

  output$plot_preview <- shiny::renderPlot({ tuned_plot() })

  output$plot_code    <- shiny::renderText({ build_plot_code(input) })

  output$data_summary <- shiny::renderTable({
    build_data_summary() %>%
      dplyr::rename(Treatment = treatment,
                    `N fish`      = n_fish,
                    `Timepoints`  = n_timepoints,
                    `N tanks`     = n_tanks)
  })

  output$download_plot <- shiny::downloadHandler(
    filename = function() {
      paste0("skinny_", input$plot_type, "_styled.png")
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
