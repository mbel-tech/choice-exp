# =============================================================================
# figures_seq_bout_module.R
# =============================================================================
# Standalone FIGURE module for the behavioural sequence analysis and the
# bout/pattern-structure analysis, drawn in the manuscript house style.
#
# WHY THIS EXISTS
# ---------------
# behavioural_sequence_analysis_choice_exp.R and bout_structure_analysis_
# choice_exp.R each fit their models AND draw their figures in one pass, so a
# purely cosmetic figure change forces a full re-fit. This module reads the CSVs
# those engines already export and redraws every figure from them. NO model is
# fitted here, NO statistic is recomputed and NO dataset is regenerated -- every
# number drawn comes from the engines' own exported output.
#
# STYLE
# -----
# Follows all_manu_graphs/design_memo.md. In particular: no grid lines, black
# axis rules only, Okabe treatment colours double-encoded with shape, jittered
# raw points over a BLACK mean crossbar + SE bar (no box plots), black
# inference ink, one-line markdown statistics captions with the "< 0.001" rule,
# and PNG + PDF at 300 dpi. Per memo section 9, these supplementary figures keep
# their titles, carry their own legend, and outline significant grid cells in
# black rather than colour.
#
# USAGE
# -----
#   source(".../figures_seq_bout_module.R")   # auto-discovers the latest runs
#   fig_module_run()                          # everything, both analyses
#   fig_module_run(which = "sequence")        # one analysis only
#   fig_module_run(seq_dir = "...", bout_dir = "...")   # pin specific runs
#   fig_module_run(out_dir = "...")           # write elsewhere
#
# Output goes to a single NON-timestamped folder:
#   D:/CHOICE R SCRIPTS/all_manu_graphs/sequence_and_patterns_graphs/
# The engines' own figures/ folders are never touched.
#
# NOT RECONSTRUCTABLE FROM EXPORTS (needs the parent engine, which holds the
# bin-level state series): the sequence state-ribbon timelines. fig_module_run()
# reports these explicitly rather than silently skipping them.
# =============================================================================

suppressMessages({
  library(data.table)
  library(ggplot2)
  library(ggtext)
  library(ggsignif)
  library(patchwork)
})

# ---- Paths ------------------------------------------------------------------
.FM_PIPE <- file.path(PROJECT_ROOT, "choice R pipeline")
.FM_OUT  <- file.path(PROJECT_ROOT, "all_manu_graphs/sequence_and_patterns_graphs")

.fm_latest_run <- function(kind = c("SEQ", "BOUT")) {
  kind <- match.arg(kind)
  root <- file.path(.FM_PIPE, "output", paste0(kind, "_output"))
  if (!dir.exists(root)) return(NA_character_)
  ds <- list.dirs(root, recursive = FALSE)
  ds <- ds[grepl(paste0(kind, "_output_\\d{8}_\\d{6}$"), ds)]
  if (!length(ds)) return(NA_character_)
  sort(ds, decreasing = TRUE)[1]
}

.fm_read <- function(dir, fname, required = TRUE) {
  fp <- file.path(dir, fname)
  if (!file.exists(fp)) {
    if (required) warning("missing export: ", fp, call. = FALSE)
    return(NULL)
  }
  data.table::fread(fp)
}

# PNG + PDF at 300 dpi, mm units, white ground (design memo section 8). The PDF
# goes through cairo_pdf so the caption's subscript/superscript glyphs embed --
# the standard PDF fonts lack those codepoints and render blank boxes.
.fm_save <- function(p, dir, stem, width, height, dpi = 300) {
  if (is.null(p)) return(invisible(FALSE))
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  ggplot2::ggsave(file.path(dir, paste0(stem, ".png")), p,
                  width = width, height = height, units = "mm",
                  dpi = dpi, bg = "white", limitsize = FALSE)
  tryCatch(
    ggplot2::ggsave(file.path(dir, paste0(stem, ".pdf")), p,
                    width = width, height = height, units = "mm",
                    device = grDevices::cairo_pdf, bg = "white",
                    limitsize = FALSE),
    error = function(e) warning(stem, " pdf save failed: ", e$message, call. = FALSE)
  )
  message("  saved: ", stem, " (png + pdf)")
  invisible(TRUE)
}

# ---- Number formatting ------------------------------------------------------
# Mirrors the project convention exactly, including the rule that a test
# statistic or effect size below 0.001 renders as "< 0.001" rather than a long
# 4-sig-fig decimal ("0.0009258") or a rounded-away zero ("0.00"). Denominator
# df keep the plain 4-sig-fig rule (a df is never < 1).
.fm_F <- function(x) {
  s <- format(signif(x, 4), scientific = FALSE, trim = TRUE)
  has_dot <- grepl(".", s, fixed = TRUE)
  s[has_dot] <- sub("\\.$", "", sub("0+$", "", s[has_dot]))
  ifelse(is.finite(x), s, "NA")
}
.fm_df <- .fm_F
.fm_Fstat <- function(x) {
  ifelse(!is.finite(x), "NA", ifelse(abs(x) < 0.001, "< 0.001", .fm_F(x)))
}
.fm_es <- function(x) {
  s <- sprintf("%.3f", x)
  s <- ifelse(grepl("0$", s), substr(s, 1, nchar(s) - 1), s)
  ifelse(is.finite(x), s, "NA")
}
.fm_es2 <- function(x) {
  ifelse(!is.finite(x), "NA", ifelse(abs(x) < 0.001, "< 0.001", .fm_es(x)))
}
.fm_p <- function(p) {
  ifelse(!is.finite(p), "NA", ifelse(p < 0.001, "< 0.001", .fm_es(p)))
}
.fm_p_txt <- function(p) {
  ifelse(!is.finite(p), "p = NA",
         ifelse(p < 0.001, "p < 0.001", paste0("p = ", .fm_es(p))))
}

# One-line ANOVA caption (design memo section 6):
#   F<sub>df1,df2</sub> = stat, p = p, eta2p = es
.fm_cap <- function(df1, df2, F, p, eta2_p) {
  sprintf("F<sub>%s,%s</sub> = %s, %s, &eta;<sup>2</sup><sub>p</sub> = %s",
          .fm_df(df1), .fm_df(df2), .fm_Fstat(F), .fm_p_txt(p), .fm_es2(eta2_p))
}
# Multi-line variant for heat-map cells, where horizontal room inside a tile is
# the binding constraint rather than line count.
.fm_cell <- function(df1, df2, F, p, eta2_p) {
  sprintf("F<sub>%s,%s</sub> = %s<br>%s<br>&eta;<sup>2</sup><sub>p</sub> = %s",
          .fm_df(df1), .fm_df(df2), .fm_Fstat(F), .fm_p_txt(p), .fm_es2(eta2_p))
}

# ---- House style (design memo sections 1-4) ---------------------------------
FM_TREAT_COLORS <- c("control" = "#D55E00", "exercise choice" = "#0072B2")
FM_TREAT_SHAPES <- c("control" = 17L,       "exercise choice" = 16L)
FM_TREAT_LINES  <- c("control" = "dashed",  "exercise choice" = "solid")
FM_STATE_COLORS <- c("Flow" = "red4", "Calm" = "lightskyblue3",
                     "High" = "tomato3", "Med" = "mediumturquoise",
                     "Low" = "goldenrod3")
FM_RAMP_HI      <- "#0072B2"

FM_JITTER_W <- 0.12
FM_PT_SIZE  <- 2.4
FM_PT_ALPHA <- 0.80
FM_MEAN_W   <- 0.40
FM_LW_MEAN  <- 0.85
FM_ERR_W    <- 0.13
FM_LW_ERR   <- 0.65

.fm_treat_lvl <- function(x) {
  factor(tolower(as.character(x)), levels = names(FM_TREAT_COLORS))
}
.fm_pretty_treat <- function(x) {
  x <- tolower(as.character(x))
  ifelse(x == "control", "Control", "Exercise\nchoice")
}

theme_fm <- function(base_size = 13) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid       = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      axis.line.x      = ggplot2::element_line(colour = "black", linewidth = 0.85),
      axis.line.y      = ggplot2::element_line(colour = "black", linewidth = 0.85),
      axis.ticks       = ggplot2::element_line(colour = "black", linewidth = 0.7),
      axis.text.x      = ggplot2::element_text(size = base_size, face = "bold",
                                               colour = "black", lineheight = 0.9),
      axis.text.y      = ggplot2::element_text(size = base_size - 2, colour = "black"),
      axis.title       = ggplot2::element_text(size = base_size, face = "bold",
                                               colour = "black"),
      strip.text       = ggplot2::element_text(size = base_size - 1, face = "bold"),
      strip.background = ggplot2::element_blank(),
      legend.title     = ggplot2::element_text(size = base_size - 2, face = "bold"),
      legend.text      = ggplot2::element_text(size = base_size - 3),
      legend.position  = "bottom",
      plot.title       = ggplot2::element_text(size = base_size + 1, face = "bold",
                                               colour = "black"),
      plot.subtitle    = ggplot2::element_text(size = base_size - 4, hjust = 0,
                                               face = "italic", colour = "grey40"),
      plot.caption     = ggtext::element_markdown(size = base_size - 2, hjust = 0,
                                                  face = "italic", colour = "black",
                                                  lineheight = 1.3,
                                                  margin = ggplot2::margin(t = 10)),
      plot.margin      = ggplot2::margin(8, 10, 8, 10)
    )
}

.fm_sem <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) NA_real_ else stats::sd(x) / sqrt(length(x))
}

# theme_fm sets axis.line.x / axis.line.y individually, so blanking the parent
# `axis.line` does NOT clear them. Heat maps have no continuous scale to anchor,
# so they drop both children plus the ticks.
.fm_no_axis_rules <- ggplot2::theme(
  axis.line.x = ggplot2::element_blank(),
  axis.line.y = ggplot2::element_blank(),
  axis.ticks  = ggplot2::element_blank()
)

# Word-wrap only (no truncation) -- used for the grid figures' row/axis labels,
# where the full pretty label must stay unique (stripping the parenthetical, as
# an earlier version did, collapsed p_stay_Flow and p_stay_Calm to the same "P"
# and silently merged two different metrics onto one grid row).
.fm_short <- function(lab, width = 30L) {
  vapply(as.character(lab), function(z) paste(strwrap(z, width = width), collapse = "\n"),
         character(1L), USE.NAMES = FALSE)
}

# The house two-group panel (design memo section 4): jittered raw points,
# black mean crossbar, black mean +/- 1 SE bar, optional black significance
# bracket, one-line markdown ANOVA caption. Deliberately NOT a box plot.
.fm_treatment_panel <- function(d, y_label, caption = "", is_sig = FALSE,
                                 show_legend = TRUE, caption_size = NULL,
                                 scale = 1) {
  d <- data.table::as.data.table(d)
  d <- d[is.finite(val)]
  if (!nrow(d)) return(NULL)
  d[, treatment := .fm_treat_lvl(treatment)]
  d <- d[!is.na(treatment)]
  if (!nrow(d)) return(NULL)
  d[, x := as.integer(treatment)]

  smry <- d[, .(mean_y = mean(val), sem_y = .fm_sem(val)), by = .(treatment, x)]

  p <- ggplot2::ggplot(d, ggplot2::aes(x = x, y = val,
                                       colour = treatment, shape = treatment)) +
    ggplot2::geom_point(position = ggplot2::position_jitter(width = FM_JITTER_W,
                                                            height = 0),
                        size = FM_PT_SIZE * scale, alpha = FM_PT_ALPHA) +
    ggplot2::geom_crossbar(data = smry, inherit.aes = FALSE,
                           ggplot2::aes(x = x, y = mean_y,
                                        ymin = mean_y, ymax = mean_y),
                           colour = "black", width = FM_MEAN_W,
                           linewidth = FM_LW_MEAN * scale) +
    ggplot2::geom_errorbar(data = smry, inherit.aes = FALSE,
                           ggplot2::aes(x = x, ymin = mean_y - sem_y,
                                        ymax = mean_y + sem_y),
                           colour = "black", width = FM_ERR_W,
                           linewidth = FM_LW_ERR * scale) +
    ggplot2::scale_x_continuous(breaks = seq_along(levels(d$treatment)),
                                labels = .fm_pretty_treat(levels(d$treatment)),
                                limits = c(0.5, length(levels(d$treatment)) + 0.5),
                                expand = c(0, 0)) +
    ggplot2::scale_colour_manual(values = FM_TREAT_COLORS, name = "Treatment",
                                 labels = c("control", "exercise choice"),
                                 drop = FALSE) +
    ggplot2::scale_shape_manual(values = FM_TREAT_SHAPES, name = "Treatment",
                                labels = c("control", "exercise choice"),
                                drop = FALSE) +
    ggplot2::labs(x = NULL, y = y_label, caption = caption) +
    theme_fm(13 * scale) +
    ggplot2::theme(axis.title.x = ggplot2::element_blank())
  if (!show_legend) p <- p + ggplot2::theme(legend.position = "none")
  # caption_size override: theme_fm()'s own caption size (base_size-2 = 11) is
  # this module's own house style for STANDALONE sequence-analysis figures.
  # Manuscript-embedded panels instead need to match whichever STEP5 caption
  # size is used by the NATIVE (non-embedded) panels in their destination
  # figure, so the embedded panel doesn't read as a different, smaller font
  # next to its neighbours -- passed explicitly per embed (see
  # fig_seq_manuscript_embeds()).
  # caption_size is NOT multiplied by scale: it was already measured/matched
  # directly against the destination figure's native caption size (2026-08-09),
  # and the panel's render width is fixed -- scaling it up the same as the
  # axis/point elements clipped the effect-size term off the end of the line
  # (confirmed on the rendered PNG). scale only compensates the axis text,
  # points and lines, which had visible headroom.
  if (!is.null(caption_size))
    p <- p + ggplot2::theme(
      plot.caption = ggtext::element_markdown(size = caption_size, hjust = 0,
                                              face = "italic", colour = "black",
                                              lineheight = 1.3,
                                              margin = ggplot2::margin(t = 10)))

  if (isTRUE(is_sig) && length(unique(d$treatment)) == 2L) {
    y_max <- max(d$val); y_rng <- max(y_max - min(d$val), 1e-10)
    p <- p +
      ggplot2::expand_limits(y = y_max + 0.22 * y_rng) +
      ggsignif::geom_signif(
        data = data.frame(xmin = 1, xmax = 2,
                          y_position = y_max + 0.10 * y_rng, annotations = "*"),
        mapping = ggplot2::aes(xmin = xmin, xmax = xmax,
                               y_position = y_position, annotations = annotations),
        manual = TRUE, tip_length = 0.02, size = 0.5 * scale, textsize = 8 * scale,
        vjust = 0.5, colour = "black", inherit.aes = FALSE)
  }
  p
}

# =============================================================================
# GRID FIGURE (the headline summary, one per analysis)
# =============================================================================
# `grid_dt` must carry: metric, term, F, df1, df2, p, eta2_p, and the facet
# columns named in `facet_row` / `facet_col`. Significant cells are outlined in
# black (memo section 9.3 -- colour is reserved for group identity, and these
# grids have no group channel).
.fm_results_grid <- function(grid_dt,
                             facet_row = NULL, facet_col = "level_lab",
                             significant_only = FALSE,
                             metric_order = NULL, term_order = NULL,
                             metric_lab = NULL, cell_size = 2.4) {
  d <- data.table::copy(data.table::as.data.table(grid_dt))
  if (!nrow(d)) return(NULL)

  d[, sig := is.finite(p) & p < 0.05]
  if ("tested" %in% names(d)) d[tested %in% c(FALSE, "FALSE"), sig := FALSE]
  if (significant_only) d <- d[sig == TRUE]
  if (!nrow(d)) return(NULL)

  d[, cell_lab := .fm_cell(df1, df2, F, p, eta2_p)]

  if (!is.null(metric_order)) {
    keep <- intersect(metric_order, unique(d$metric))
    d[, metric := factor(metric, levels = rev(keep))]
  } else {
    d[, metric := factor(metric, levels = rev(sort(unique(as.character(metric)))))]
  }
  if (!is.null(term_order)) {
    keep_t <- intersect(term_order, unique(as.character(d$term_lab)))
    d[, term_lab := factor(as.character(term_lab), levels = keep_t)]
  }

  if (!is.null(metric_lab)) {
    lv <- levels(d$metric)
    levels(d$metric) <- .fm_short(ifelse(lv %in% names(metric_lab),
                                         metric_lab[lv], lv), width = 26L)
  }

  p_out <- ggplot2::ggplot(d, ggplot2::aes(x = term_lab, y = metric, fill = eta2_p)) +
    ggplot2::geom_tile(ggplot2::aes(colour = sig, linewidth = sig),
                       width = 0.97, height = 0.97) +
    ggtext::geom_richtext(ggplot2::aes(label = cell_lab), size = cell_size,
                          colour = "black", fill = NA, label.color = NA,
                          label.padding = grid::unit(rep(1, 4), "pt")) +
    ggplot2::scale_fill_gradient(low = "white", high = FM_RAMP_HI,
                                 limits = c(0, 1),
                                 name = expression(eta[p]^2)) +
    ggplot2::scale_colour_manual(values = c("TRUE" = "black", "FALSE" = "grey85"),
                                 guide = "none") +
    ggplot2::scale_linewidth_manual(values = c("TRUE" = 1.1, "FALSE" = 0.4),
                                    guide = "none") +
    ggplot2::labs(x = NULL, y = NULL) +
    theme_fm(12) + .fm_no_axis_rules +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 20, hjust = 1,
                                                       size = 11, face = "bold"),
                   axis.text.y = ggplot2::element_text(size = 10, face = "bold",
                                                       colour = "black"))

  if (!is.null(facet_row)) {
    p_out + ggplot2::facet_grid(stats::reformulate(facet_col, facet_row),
                                scales = "free", space = "free")
  } else {
    p_out + ggplot2::facet_grid(stats::reformulate(facet_col, "."),
                                scales = "free_x", space = "free_x")
  }
}

# =============================================================================
# SEQUENCE ANALYSIS FIGURES
# =============================================================================
SEQ_TERM_LAB <- c(treatment = "treatment", timepoint_f = "interval",
                  `treatment:timepoint_f` = "treatment \u00d7 interval")
SEQ_LEVEL_LAB <- c(trial = "Level A \u2014 trial (n = 16)",
                   trial_x_timepoint = "Level B \u2014 trial \u00d7 interval (n = 48)")
SEQ_METRIC_ORDER <- c("switch_rate", "entropy_rate",
                      "occ_Flow", "occ_High", "occ_Low",
                      "dwell_Flow", "dwell_Calm", "dwell_High", "dwell_Low",
                      "p_stay_Flow", "p_stay_Calm", "p_stay_High", "p_stay_Low",
                      "obs_min")
# Axis labels: the metric column names are analysis shorthand, not reader-facing.
SEQ_METRIC_LAB <- c(
  switch_rate  = "Switch rate (transitions / min)",
  entropy_rate = "Normalised entropy rate (0\u20131)",
  occ_Flow     = "Flow occupancy (time share)",
  occ_High     = "High-flow occupancy (time share)",
  occ_Low      = "Low-flow occupancy (time share)",
  dwell_Flow   = "Mean flow dwell time (s)",
  dwell_Calm   = "Mean calm dwell time (s)",
  dwell_High   = "Mean high-flow dwell time (s)",
  dwell_Low    = "Mean low-flow dwell time (s)",
  p_stay_Flow  = "P(stay in flow)",
  p_stay_Calm  = "P(stay in calm)",
  p_stay_High  = "P(stay in high flow)",
  p_stay_Low   = "P(stay in low flow)",
  obs_min      = "Observation time (min)")
.seq_lab <- function(m) ifelse(m %in% names(SEQ_METRIC_LAB), SEQ_METRIC_LAB[m], m)

fig_seq_all <- function(seq_dir = .fm_latest_run("SEQ"), out_dir = .FM_OUT) {
  stopifnot(!is.na(seq_dir), dir.exists(seq_dir))
  message("SEQUENCE figures <- ", basename(seq_dir))
  made <- character(0); skipped <- character(0)

  # --- (1) Transition-probability heat maps, per alphabet, BY TREATMENT ------
  for (alph in c("binary", "graded")) {
    tm <- .fm_read(seq_dir, sprintf("transition_matrix_pooled_%s.csv", alph),
                   required = FALSE)
    if (is.null(tm) || !nrow(tm) ||
        !all(c("from", "to", "prob") %in% names(tm))) {
      skipped <- c(skipped, paste0("transition_heat_", alph)); next
    }
    setDT(tm)
    lv <- if (alph == "binary") c("Calm", "Flow") else c("Low", "Med", "High")
    tm[, `:=`(from = factor(from, levels = rev(lv)), to = factor(to, levels = lv))]
    if ("treatment" %in% names(tm))
      tm[, treatment := factor(.fm_pretty_treat(treatment),
                               levels = c("Control", "Exercise\nchoice"))]
    p <- ggplot2::ggplot(tm, ggplot2::aes(x = to, y = from, fill = prob)) +
      ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
      ggplot2::geom_text(ggplot2::aes(label = .fm_es(prob)), size = 4,
                         colour = "black") +
      ggplot2::scale_fill_gradient(low = "white", high = FM_RAMP_HI,
                                   limits = c(0, 1), name = "P(to | from)") +
      ggplot2::labs(x = "State at next 5-s bin", y = "State at current 5-s bin") +
      theme_fm(13) + .fm_no_axis_rules
    if ("treatment" %in% names(tm)) p <- p + ggplot2::facet_wrap(~ treatment)
    if (isTRUE(.fm_save(p, out_dir, sprintf("seq_transition_heat_%s", alph),
                        if (alph == "binary") 190 else 220, 150)))
      made <- c(made, paste0("transition_heat_", alph))
  }

  # --- (2) Per-trial metric panels, house idiom + one-line ANOVA caption -----
  av_trial <- .fm_read(seq_dir, "seq_anova_trial_level.csv", required = FALSE)
  for (alph in c("binary", "graded")) {
    mt <- .fm_read(seq_dir, sprintf("seq_metrics_per_trial_%s.csv", alph),
                   required = FALSE)
    if (is.null(mt) || !nrow(mt) || !"treatment" %in% names(mt)) {
      skipped <- c(skipped, paste0("metrics_", alph)); next
    }
    setDT(mt)
    show <- intersect(c("switch_rate", "entropy_rate", "occ_Flow", "occ_High",
                        "dwell_Flow", "dwell_Calm", "dwell_High", "p_stay_Flow",
                        "p_stay_High"),
                      names(mt))
    for (met in show) {
      d <- mt[is.finite(get(met)), .(val = mean(get(met), na.rm = TRUE)),
              by = .(trial, treatment)]
      if (!nrow(d)) next
      cap <- ""; sig <- FALSE
      if (!is.null(av_trial) && nrow(av_trial)) {
        r <- as.data.table(av_trial)[alphabet == alph & metric == met &
                                       term == "treatment"]
        if (nrow(r)) {
          cap <- paste0("Treatment ", .fm_cap(r$df1[1], r$df2[1], r$F[1],
                                              r$p[1], r$eta2_p[1]))
          sig <- is.finite(r$p[1]) && r$p[1] < 0.05
        }
      }
      p <- .fm_treatment_panel(d, y_label = .seq_lab(met),
                               caption = cap, is_sig = sig)
      if (isTRUE(.fm_save(p, out_dir, sprintf("seq_metric_%s_%s", alph, met),
                          140, 130)))
        made <- c(made, sprintf("metric_%s_%s", alph, met))
    }
  }

  # --- (3) The two grid summaries ------------------------------------------
  grid <- .fm_read(seq_dir, "seq_results_grid.csv", required = FALSE)
  if (!is.null(grid) && nrow(grid)) {
    g <- as.data.table(grid)
    g[, term_lab  := factor(SEQ_TERM_LAB[term], levels = unname(SEQ_TERM_LAB))]
    g[, level_lab := factor(SEQ_LEVEL_LAB[level], levels = unname(SEQ_LEVEL_LAB))]
    p_all <- .fm_results_grid(
      g, facet_row = "alphabet", facet_col = "level_lab",
      significant_only = FALSE, metric_order = SEQ_METRIC_ORDER,
      term_order = unname(SEQ_TERM_LAB), metric_lab = SEQ_METRIC_LAB)
    if (isTRUE(.fm_save(p_all, out_dir, "seq_grid_complete", 330, 250)))
      made <- c(made, "grid_complete")

    p_sig <- .fm_results_grid(
      g, facet_row = "alphabet", facet_col = "level_lab",
      significant_only = TRUE, metric_order = SEQ_METRIC_ORDER,
      term_order = unname(SEQ_TERM_LAB), metric_lab = SEQ_METRIC_LAB)
    if (isTRUE(.fm_save(p_sig, out_dir, "seq_grid_significant", 330, 220)))
      made <- c(made, "grid_significant")
  } else skipped <- c(skipped, "grid")

  skipped <- c(skipped, "state_ribbons (needs bin-level data held only inside the engine)")
  list(made = made, skipped = skipped, out_dir = out_dir)
}

# =============================================================================
# MANUSCRIPT-SPECIFIC EMBEDS (2026-08-09)
# =============================================================================
# One-off panels embedded (as raster PNGs, same technique already used for the
# legend strips) into the STEP5 composite figures:
#   - Figure "7"/zone-preference: p_stay_Flow, p_stay_High as panels E/F.
#   - Figure "8"/collective-trial: switch_rate replaces panel A
#     (switches_per_session), with a manuscript-specific y-axis label the
#     generic per-trial loop above does not use; entropy_rate (graded) fills
#     the new panel F.
# Same data, same CSV-sourced caption, same .fm_treatment_panel() house style
# as every other panel in this module. Every one of these four re-renders
# WITHOUT its own legend (show_legend = FALSE, 2026-08-09) and overwrites the
# generic loop's same-named file: the parent STEP5 figure already carries one
# shared legend for the whole composite, so a second legend baked into the
# embedded PNG was redundant. The generic-loop version (with legend) is what
# a reader gets when browsing the standalone sequence-analysis figure set.
fig_seq_manuscript_embeds <- function(seq_dir = .fm_latest_run("SEQ"), out_dir = .FM_OUT) {
  stopifnot(!is.na(seq_dir), dir.exists(seq_dir))
  av_trial <- .fm_read(seq_dir, "seq_anova_trial_level.csv", required = FALSE)
  made <- character(0)

  # dest_w/dest_h: render dimensions matched to the EXACT destination panel
  # cell in the parent STEP5 figure (not a generic 140x130mm square). A
  # raster embedded via annotation_custom()+coord_equal() preserves its own
  # source aspect ratio, so if that aspect doesn't match the cell it gets
  # letterboxed (shrunk to fit, with blank padding) -- which is exactly why
  # these panels looked smaller than their native ggplot neighbours before
  # this fix (measured directly on the rendered PNGs, 2026-08-09):
  #   Figure_8 (.fig7, 280x320mm, 2 cols x rows heights c(6,6,6,0.5)):
  #     E/F cell = (280-10)/2 wide x 6/18.5*(320-10) tall = 135 x 100.5mm.
  #   Figure_9 (.fig8, 320x220mm, 3 cols x rows heights c(6,0.5)):
  #     A/F cell = (320-10)/3 wide x 6/6.5/2*(220-10) tall = 103.3 x 96.9mm.
  # cap_size: matched to the NATIVE panels' own caption theme in that same
  # destination figure -- Figure_8 A/B use BASE_THEME's default (14pt, no
  # override in .mg_make_alr_scatter); Figure_9 B-E use .f8_cap_theme (9pt).
  # scale: compensates for a residual size deficit measured directly on the
  # rendered composite (2026-08-09) -- even after dest_w/dest_h were matched
  # to the destination cell's mm dimensions, Figure_8's E/F axis-box width
  # measured ~933px against panel A/B's ~1167px (80%), UNCHANGED by the
  # dest_w/dest_h fix (ruling out simple letterboxing as the sole cause; some
  # additional shrinkage survives the wrap_elements()/patchwork embedding
  # that was not fully traced to its root cause). Compensated empirically by
  # rendering the source panel's own text/points/lines 1/0.80 = 1.25x larger
  # than the nominal target size, so the FINAL embedded appearance matches.
  .one <- function(alph, met, y_label, stem, dest_w = 140, dest_h = 130, cap_size = 11,
                    scale = 1) {
    mt <- .fm_read(seq_dir, sprintf("seq_metrics_per_trial_%s.csv", alph), required = FALSE)
    if (is.null(mt) || !met %in% names(mt)) return(invisible(FALSE))
    setDT(mt)
    d <- mt[is.finite(get(met)), .(val = mean(get(met), na.rm = TRUE)),
            by = .(trial, treatment)]
    if (!nrow(d)) return(invisible(FALSE))
    cap <- ""; sig <- FALSE
    if (!is.null(av_trial) && nrow(av_trial)) {
      r <- as.data.table(av_trial)[alphabet == alph & metric == met & term == "treatment"]
      if (nrow(r)) {
        cap <- paste0("Treatment ", .fm_cap(r$df1[1], r$df2[1], r$F[1],
                                            r$p[1], r$eta2_p[1]))
        sig <- is.finite(r$p[1]) && r$p[1] < 0.05
      }
    }
    p <- .fm_treatment_panel(d, y_label = y_label, caption = cap, is_sig = sig,
                             show_legend = FALSE, caption_size = cap_size, scale = scale)
    if (isTRUE(.fm_save(p, out_dir, stem, dest_w, dest_h))) made <<- c(made, stem)
    invisible(TRUE)
  }

  # Figure_9 (collective-trial) panels A and F: cell = 103.3 x 96.9mm, caption 9pt.
  # scale = 1.25 to compensate the same residual shrinkage documented above for
  # Figure_8's (now-removed) E/F embeds -- same raster/coord_equal embedding
  # chain, same measured ~80% size deficit against native neighbour panels.
  .one("binary", "switch_rate", "School transition\nbetween main zone",
       "seq_switch_rate_manuscript", dest_w = 103.3, dest_h = 96.9, cap_size = 9,
       scale = 1.25)
  .one("graded", "entropy_rate", "Normalised entropy rate",
       "seq_metric_graded_entropy_rate", dest_w = 103.3, dest_h = 96.9, cap_size = 9,
       scale = 1.25)
  # Figure_8 (zone/ALR) panel E: cell = 135 x 100.5mm, caption 14pt (matching
  # panels A/B, the same simple two-group scatter panel type), scale = 1.25 to
  # compensate the measured residual shrinkage (see above). Re-added 2026-08-09
  # after a brief removal -- panel E is back, panel F is now the native
  # figure7bis panel (built directly in activity_analysis_STATS_choice_exp.R,
  # no raster embedding needed there). p_stay_High is NOT re-added here: it is
  # not used by any current manuscript panel; fig_seq_all()'s generic per-trial
  # loop still produces its standalone (legend-bearing, un-scaled) render.
  # Axis title spelled out (and wrapped after "staying") for the manuscript
  # embed only -- Figure_8 panel E. The compact SEQ_METRIC_LAB form
  # ("P(stay in flow)") is retained for the module's own standalone renders.
  .one("binary", "p_stay_Flow", "Probability of staying\nin flow zone",
       "seq_metric_binary_p_stay_Flow", dest_w = 135, dest_h = 100.5, cap_size = 14,
       scale = 1.25)

  list(made = made, out_dir = out_dir)
}

# =============================================================================
# FIGURE_PATTERN -- SUPERSEDED 2026-08-09, see fig11_temporal_structure() below
# =============================================================================
# An earlier "figure_pattern" composite (p_stay_Flow, dwell_Flow, dwell_High,
# memory_flow) stood here. It predates the collinearity screen against the
# ALR primary endpoint (collinearity_vs_ALR.R): p_stay_Flow and dwell_High
# both clear |r| >= 0.70 against ALR(flow) and are reported as DESCRIPTIVE
# alongside switch_rate/entropy_rate (Figure 9), not as independent temporal
# structure. The correct independent set -- |r| < 0.70 vs every ALR outcome,
# confirmed mutually non-collinear -- is burst_flow, dwell_Flow, dwell_Calm,
# memory_flow, now built as Figure 11 by fig11_temporal_structure(). Its
# stale output (figure_pattern.png/.pdf/_caption.txt) has been deleted from
# all_manu_graphs/; do not regenerate under that name.

# =============================================================================
# PATTERN / BOUT-STRUCTURE FIGURES
# =============================================================================
BOUT_TERM_LAB <- c(treatment = "treatment", timepoint_f = "interval",
                   `treatment:timepoint_f` = "treatment \u00d7 interval")
BOUT_LEVEL_LAB <- c(trial = "Level A \u2014 trial (n = 16)",
                    trial_x_timepoint = "Level B \u2014 trial \u00d7 interval (n = 48)")
BOUT_METRIC_ORDER <- c("commitment_index", "max_flow_bout_s", "max_calm_bout_s",
                       "n_flow_ge30", "long_flow_frac", "flow_centroid",
                       "t_first_sustained_s", "si_markov_max", "engage_depth",
                       "mean_intensity_sel", "med_bout_flow", "med_bout_calm",
                       "burst_flow", "burst_calm", "memory_flow",
                       "bout_rate_flow", "onset_hazard_flow",
                       "mean_bout_flow", "mean_bout_calm",
                       "cv_bout_flow", "cv_bout_calm", "log2_med_ratio")
BOUT_METRIC_LAB <- c(
  commitment_index    = "Commitment index (longest bout / interval)",
  max_flow_bout_s     = "Longest flow bout (s)",
  max_calm_bout_s     = "Longest calm bout (s)",
  n_flow_ge30         = "Flow bouts \u2265 30 s (count)",
  long_flow_frac      = "Time in flow bouts \u2265 60 s (share)",
  flow_centroid       = "Flow centroid (0 = front-loaded, 1 = back-loaded)",
  t_first_sustained_s = "Latency to first sustained flow bout (s)",
  si_markov_max       = "Structure index vs Markov null (log\u2082)",
  engage_depth        = "Engagement depth (mean graded-state rank)",
  mean_intensity_sel  = "Selected intensity (mean intensity rank)",
  med_bout_flow       = "Median flow bout (s, uncensored)",
  med_bout_calm       = "Median calm bout (s, uncensored)",
  burst_flow          = "Burstiness of flow bouts (\u22121\u20131)",
  burst_calm          = "Burstiness of calm bouts (\u22121\u20131)",
  memory_flow         = "Memory of flow bouts (lag-1 r)",
  bout_rate_flow      = "Flow bout rate (bouts / min)",
  onset_hazard_flow   = "Flow onset hazard (per min)",
  mean_bout_flow      = "Mean flow bout (s)",
  mean_bout_calm      = "Mean calm bout (s)",
  cv_bout_flow        = "CV of flow bout durations",
  cv_bout_calm        = "CV of calm bout durations",
  log2_med_ratio      = "Bout asymmetry, log\u2082(flow / calm median)")
# A few exports (repeatability ICC, rank stability) carry sequence-module
# metrics alongside bout metrics for a joint comparison, so fall back to the
# sequence label dictionary before giving up and showing the raw column name.
.bout_lab <- function(m) {
  ifelse(m %in% names(BOUT_METRIC_LAB), BOUT_METRIC_LAB[m],
         ifelse(m %in% names(SEQ_METRIC_LAB), SEQ_METRIC_LAB[m], m))
}

fig_bout_all <- function(bout_dir = .fm_latest_run("BOUT"), out_dir = .FM_OUT) {
  stopifnot(!is.na(bout_dir), dir.exists(bout_dir))
  message("PATTERN (bout) figures <- ", basename(bout_dir))
  made <- character(0); skipped <- character(0)

  # --- (1) Bout raster ------------------------------------------------------
  inv <- .fm_read(bout_dir, "bout_inventory.csv", required = FALSE)
  if (!is.null(inv) && nrow(inv) &&
      all(c("trial", "timepoint", "treatment", "state", "start_bin", "end_bin")
          %in% names(inv))) {
    d <- as.data.table(inv)[!is.na(state)]
    if ("alphabet" %in% names(d)) d <- d[alphabet == "binary"]
    BIN_S <- 5
    d[, lbl := sprintf("%s \u2014 school %s, interval %s",
                       ifelse(tolower(treatment) == "control",
                              "Control", "Exercise choice"), trial, timepoint)]
    d[, lbl := factor(lbl, levels = unique(lbl[order(treatment, trial, timepoint)]))]
    p <- ggplot2::ggplot(d, ggplot2::aes(xmin = start_bin * BIN_S,
                                         xmax = (end_bin + 1) * BIN_S,
                                         ymin = 0, ymax = 1, fill = state)) +
      ggplot2::geom_rect() +
      ggplot2::facet_wrap(~ lbl, ncol = 1, strip.position = "left") +
      ggplot2::scale_fill_manual(values = FM_STATE_COLORS, name = "State") +
      ggplot2::scale_y_continuous(breaks = NULL) +
      ggplot2::labs(x = "Time within interval (s)", y = NULL) +
      theme_fm(11) +
      ggplot2::theme(axis.line.y = ggplot2::element_blank(),
                     strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 1,
                                                               size = 5, face = "plain"))
    if (isTRUE(.fm_save(p, out_dir, "bout_raster_full", 250, 420)))
      made <- c(made, "raster")
  } else skipped <- c(skipped, "raster")

  # --- (2) Kaplan-Meier median bout duration (censoring-aware) -------------
  # km_bout_durations.csv exports the KM SUMMARY (median + CI per state x
  # treatment), not the survival curve, so this is a median-with-CI comparison.
  # The curve itself would need the per-bout survival fit held in the engine.
  km <- .fm_read(bout_dir, "km_bout_durations.csv", required = FALSE)
  if (!is.null(km) && nrow(km) &&
      all(c("state", "treatment", "median_km") %in% names(km))) {
    d <- as.data.table(km)
    d[, treatment := .fm_treat_lvl(treatment)]
    d[, state := factor(state, levels = c("Flow", "Calm"))]
    has_ci <- all(c("lo", "hi") %in% names(d))
    dodge <- ggplot2::position_dodge(width = 0.5)
    p <- ggplot2::ggplot(d, ggplot2::aes(x = state, y = median_km,
                                         colour = treatment, shape = treatment))
    if (has_ci)
      p <- p + ggplot2::geom_errorbar(ggplot2::aes(ymin = lo, ymax = hi),
                                      width = FM_ERR_W, linewidth = FM_LW_ERR,
                                      colour = "black", position = dodge)
    p <- p +
      ggplot2::geom_point(size = 3.8, position = dodge) +
      ggplot2::scale_colour_manual(values = FM_TREAT_COLORS, name = "Treatment") +
      ggplot2::scale_shape_manual(values = FM_TREAT_SHAPES, name = "Treatment") +
      ggplot2::labs(x = NULL,
                    y = "Median bout duration (s, Kaplan\u2013Meier \u00b1 95% CI)") +
      theme_fm(13)
    if (isTRUE(.fm_save(p, out_dir, "bout_km_median_duration", 160, 130)))
      made <- c(made, "km_median_duration")
  } else skipped <- c(skipped, "km_median_duration")

  # --- (3) Per-school metric panels, house idiom + ANOVA caption ------------
  sch <- .fm_read(bout_dir, "bout_metrics_per_school.csv", required = FALSE)
  av  <- .fm_read(bout_dir, "bout_anova_trial_level.csv", required = FALSE)
  if (!is.null(sch) && nrow(sch) && "treatment" %in% names(sch)) {
    d0 <- as.data.table(sch)
    show <- intersect(c("commitment_index", "max_flow_bout_s", "n_flow_ge30",
                        "long_flow_frac", "flow_centroid", "t_first_sustained_s",
                        "engage_depth", "mean_intensity_sel", "si_markov_max",
                        "burst_flow", "memory_flow"),
                      names(d0))
    for (met in show) {
      d <- d0[is.finite(get(met)), .(val = get(met), treatment)]
      if (!nrow(d)) next
      cap <- ""; sig <- FALSE
      if (!is.null(av) && nrow(av)) {
        r <- as.data.table(av)[metric == met & term == "treatment"]
        if (nrow(r)) {
          cap <- paste0("Treatment ", .fm_cap(r$df1[1], r$df2[1], r$F[1],
                                              r$p[1], r$eta2_p[1]))
          sig <- is.finite(r$p[1]) && r$p[1] < 0.05
          # A metric excluded by the collinearity gate is reported but never
          # claimed: the caption says so and the significance bracket is
          # withheld, so the panel cannot be read as an inferential result.
          if ("tested" %in% names(r) && !isTRUE(as.logical(r$tested[1]))) {
            rsn <- if ("exclude_reason" %in% names(r)) r$exclude_reason[1] else ""
            cap <- paste0(cap, "<br>*not tested \u2014 ", rsn, "*")
            sig <- FALSE
          }
        }
      }
      p <- .fm_treatment_panel(d, y_label = .bout_lab(met),
                               caption = cap, is_sig = sig)
      if (isTRUE(.fm_save(p, out_dir, sprintf("bout_metric_%s", met), 140, 135)))
        made <- c(made, paste0("metric_", met))
    }
  } else skipped <- c(skipped, "metric panels")

  # --- (4) Observed vs Markov/shuffle null ---------------------------------
  nullsum <- .fm_read(bout_dir, "bout_null_summary.csv", required = FALSE)
  if (!is.null(nullsum) && nrow(nullsum) &&
      all(c("null_type", "statistic", "treatment", "median_log2_ratio")
          %in% names(nullsum))) {
    d <- as.data.table(nullsum)
    d[, treatment := .fm_treat_lvl(treatment)]
    d[, statistic := factor(statistic, levels = BOUT_METRIC_ORDER[
        BOUT_METRIC_ORDER %in% unique(statistic)])]
    p <- ggplot2::ggplot(d, ggplot2::aes(x = statistic, y = median_log2_ratio,
                                         fill = treatment)) +
      ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8),
                        width = 0.72, colour = "black", linewidth = 0.3) +
      ggplot2::geom_hline(yintercept = 0, linewidth = 0.6, colour = "black") +
      ggplot2::facet_wrap(~ null_type) +
      ggplot2::scale_fill_manual(values = FM_TREAT_COLORS, name = "Treatment") +
      ggplot2::scale_x_discrete(labels = .bout_lab) +
      ggplot2::labs(x = NULL,
                    y = "Median log\u2082(observed / null)\n0 = indistinguishable from null") +
      theme_fm(12) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1,
                                                         size = 9, face = "bold"),
                     plot.margin = ggplot2::margin(8, 10, 8, 55)) +
      # The first bar's rotated label is the widest thing on the canvas and
      # ggsave's declared width otherwise clips it; reserve blank space to its
      # left the same way plot.margin does for every other edge.
      ggplot2::coord_cartesian(clip = "off")
    if (isTRUE(.fm_save(p, out_dir, "bout_observed_vs_null", 290, 155)))
      made <- c(made, "observed_vs_null")
  } else skipped <- c(skipped, "observed_vs_null")

  # --- (5) Repeatability (ICC) forest --------------------------------------
  icc <- .fm_read(bout_dir, "bout_repeatability_icc.csv", required = FALSE)
  if (!is.null(icc) && nrow(icc) && all(c("metric", "icc") %in% names(icc))) {
    d <- as.data.table(icc)
    if ("subset" %in% names(d)) d <- d[subset == "all"]
    if ("adjusted" %in% names(d)) d <- d[adjusted %in% c(TRUE, "TRUE")]
    if (nrow(d)) {
      d[, lab := .bout_lab(metric)]
      has_ci <- all(c("icc_lo", "icc_hi") %in% names(d))
      p <- ggplot2::ggplot(d, ggplot2::aes(x = icc,
                                           y = stats::reorder(lab, icc)))
      if (has_ci)
        p <- p + ggplot2::geom_errorbar(ggplot2::aes(xmin = icc_lo, xmax = icc_hi),
                                        orientation = "y", width = FM_ERR_W,
                                        linewidth = FM_LW_ERR, colour = "black")
      p <- p + ggplot2::geom_point(size = 3.2, colour = FM_RAMP_HI, shape = 16) +
        ggplot2::scale_x_continuous(limits = c(0, 1)) +
        ggplot2::labs(x = "Adjusted ICC (school repeatability ± bootstrap 95% CI)",
                      y = NULL) +
        theme_fm(12) +
        ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10, face = "bold",
                                                           colour = "black"))
      if (isTRUE(.fm_save(p, out_dir, "bout_icc_forest", 200, 140)))
        made <- c(made, "icc_forest")
    } else skipped <- c(skipped, "icc_forest")
  } else skipped <- c(skipped, "icc_forest")

  # --- (6) Rank stability across intervals ---------------------------------
  rk <- .fm_read(bout_dir, "bout_rank_stability.csv", required = FALSE)
  if (!is.null(rk) && nrow(rk) && all(c("metric", "kendall_W", "p") %in% names(rk))) {
    d <- as.data.table(rk)
    if ("subset" %in% names(d)) d <- d[subset == "all"]
    d[, sig := is.finite(p) & p < 0.05]
    d[, lab := .bout_lab(metric)]
    p <- ggplot2::ggplot(d, ggplot2::aes(x = kendall_W,
                                         y = stats::reorder(lab, kendall_W),
                                         fill = sig)) +
      ggplot2::geom_col(width = 0.7, colour = "black", linewidth = 0.3) +
      ggplot2::scale_fill_manual(values = c("TRUE" = FM_RAMP_HI, "FALSE" = "grey85"),
                                 name = "p < 0.05") +
      ggplot2::scale_x_continuous(limits = c(0, 1)) +
      ggplot2::labs(
        x = "Kendall's W across the three intervals\n1 = schools hold their relative order exactly",
        y = NULL) +
      theme_fm(12) +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10, face = "bold",
                                                         colour = "black"))
    if (isTRUE(.fm_save(p, out_dir, "bout_rank_stability", 200, 140)))
      made <- c(made, "rank_stability")
  } else skipped <- c(skipped, "rank_stability")

  # --- (7) Collinearity with the established sequence metrics --------------
  col <- .fm_read(bout_dir, "metric_collinearity.csv", required = FALSE)
  if (!is.null(col) && nrow(col) &&
      all(c("metric_new", "metric_seq", "r_pearson", "flagged") %in% names(col))) {
    d <- as.data.table(col)
    d[, flag := flagged %in% c(TRUE, "TRUE")]
    p <- ggplot2::ggplot(d, ggplot2::aes(x = metric_seq, y = metric_new,
                                         fill = r_pearson)) +
      ggplot2::geom_tile(ggplot2::aes(colour = flag, linewidth = flag),
                         width = 0.95, height = 0.95) +
      ggplot2::geom_text(ggplot2::aes(label = .fm_es(r_pearson)), size = 2.7,
                         colour = "black") +
      ggplot2::scale_fill_gradient2(low = "#B2432F", mid = "white",
                                    high = "#2E5FA3", midpoint = 0,
                                    limits = c(-1, 1), name = "r") +
      ggplot2::scale_colour_manual(values = c("TRUE" = "black", "FALSE" = "grey88"),
                                   guide = "none") +
      ggplot2::scale_linewidth_manual(values = c("TRUE" = 1.1, "FALSE" = 0.4),
                                      guide = "none") +
      ggplot2::labs(x = "Established sequence metric",
                    y = "Candidate bout metric") +
      theme_fm(11) + .fm_no_axis_rules +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 20, hjust = 1,
                                                         size = 10, face = "bold"),
                     axis.text.y = ggplot2::element_text(size = 9, face = "bold",
                                                         colour = "black"))
    if (isTRUE(.fm_save(p, out_dir, "bout_collinearity_heat", 200, 250)))
      made <- c(made, "collinearity_heat")
  } else skipped <- c(skipped, "collinearity_heat")

  # --- (8) The two grid summaries ------------------------------------------
  grid <- .fm_read(bout_dir, "bout_results_grid.csv", required = FALSE)
  if (!is.null(grid) && nrow(grid)) {
    g <- as.data.table(grid)
    g[, term_lab  := factor(BOUT_TERM_LAB[term], levels = unname(BOUT_TERM_LAB))]
    g[, level_lab := factor(BOUT_LEVEL_LAB[level], levels = unname(BOUT_LEVEL_LAB))]
    facet_r <- if ("tier" %in% names(g) && length(unique(g$tier)) > 1) "tier" else NULL
    p_all <- .fm_results_grid(
      g, facet_row = facet_r, facet_col = "level_lab",
      significant_only = FALSE, metric_order = BOUT_METRIC_ORDER,
      term_order = unname(BOUT_TERM_LAB), metric_lab = BOUT_METRIC_LAB)
    if (isTRUE(.fm_save(p_all, out_dir, "bout_grid_complete", 330, 300)))
      made <- c(made, "grid_complete")

    p_sig <- .fm_results_grid(
      g, facet_row = facet_r, facet_col = "level_lab",
      significant_only = TRUE, metric_order = BOUT_METRIC_ORDER,
      term_order = unname(BOUT_TERM_LAB), metric_lab = BOUT_METRIC_LAB)
    if (isTRUE(.fm_save(p_sig, out_dir, "bout_grid_significant", 330, 220)))
      made <- c(made, "grid_significant")
  } else skipped <- c(skipped, "grid")

  list(made = made, skipped = skipped, out_dir = out_dir)
}

# =============================================================================
# FIGURE 11 -- INDEPENDENT TEMPORAL STRUCTURE (2026-08-09)
# =============================================================================
# The 4 metrics that survive the school-level collinearity screen against the
# ALR primary endpoint (collinearity_vs_ALR.R, |r| < 0.70 vs every ALR
# outcome): burst_flow, dwell_Flow, dwell_Calm, memory_flow. All four are
# null at Level A (trial, n = 16). Two of them are NOT null at Level B (trial
# x interval, n = 48): burst_flow's treatment main effect (p = 0.026) and
# dwell_Calm's treatment x interval interaction (p = 0.031) -- an independent
# temporal-structure signal that is only resolvable at session level. Panels
# E-F exist to show exactly that Level A/Level B contrast, not merely to
# repeat A-D at finer grain.
#   A-D (Level A, trial):        .fm_treatment_panel(), same idiom as every
#                                 other per-metric panel in this module.
#   E-F (Level B, trial x interval): a new interval line-plot helper, drawn
#                                 in the same house style as the treatment
#                                 panels (Okabe colours, black inference ink)
#                                 since this module has no existing interval
#                                 line-plot (that idiom lives in the STEP5
#                                 script, for Figure 10, and is not reused
#                                 here to avoid a cross-script dependency).
# =============================================================================

# Two-line caption: Treatment term on line 1, Treatment x Interval on line 2
# (the interval main effect is a nuisance term here, omitted to keep the
# panel legible -- full three-term statistics are in the manuscript text).
.fm_cap_levelB <- function(av_tp, met, alph_filter = NULL) {
  d <- as.data.table(av_tp)
  if (!is.null(alph_filter) && "alphabet" %in% names(d)) d <- d[alphabet == alph_filter]
  r_t  <- d[metric == met & term == "treatment"]
  r_ix <- d[metric == met & term == "treatment:timepoint_f"]
  l1 <- if (nrow(r_t))
    paste0("Treatment ", .fm_cap(r_t$df1[1], r_t$df2[1], r_t$F[1], r_t$p[1], r_t$eta2_p[1]))
  else "Treatment: NA"
  l2 <- if (nrow(r_ix))
    paste0("Treatment×Interval ",
           .fm_cap(r_ix$df1[1], r_ix$df2[1], r_ix$F[1], r_ix$p[1], r_ix$eta2_p[1]))
  else "Treatment×Interval: NA"
  paste0(l1, "<br>", l2)
}

# Session-level (n = 48) interval line plot: mean +/- SE per treatment x
# timepoint, house colours/shapes/linetypes, two-line Level-B caption.
#
# cld: optional data frame (treatment, timepoint, .group) of compact-letter
#      groupings. Letters are drawn CENTRED on each cell's x position, at the
#      upper limit of that cell's error bar plus 10% of the plotted data range
#      (2026-08-10 house rule, applied identically in the STEP5 module).
# dodge: horizontal separation between the two treatments. Needed where the two
#      treatment means nearly coincide at an interval (dwell_Calm interval 3),
#      which otherwise hides one marker completely behind the other.
.fm_interval_panel <- function(d0, y_col, y_label, cap, cld = NULL, dodge = 0) {
  d <- as.data.table(d0)[is.finite(get(y_col))]
  d[, treatment := .fm_treat_lvl(treatment)]
  d <- d[!is.na(treatment)]
  d[, tp := as.integer(timepoint)]
  smry <- d[, .(mean_y = mean(get(y_col)), sem_y = .fm_sem(get(y_col))),
            by = .(treatment, tp)]
  # signed offset: control shifted left, exercise choice right
  smry[, xpos := tp + ifelse(as.character(treatment) == "control", -dodge/2, dodge/2)]
  ymin <- min(smry$mean_y - smry$sem_y, na.rm = TRUE)
  ymax <- max(smry$mean_y + smry$sem_y, na.rm = TRUE)
  yrange <- ymax - ymin
  p <- ggplot2::ggplot(smry, ggplot2::aes(x = xpos, y = mean_y, colour = treatment,
                                          shape = treatment, linetype = treatment,
                                          group = treatment)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = FM_PT_SIZE * 1.6) +
    ggplot2::geom_errorbar(ggplot2::aes(x = xpos, ymin = mean_y - sem_y, ymax = mean_y + sem_y),
                           colour = "black", width = 0.10, linewidth = FM_LW_ERR,
                           inherit.aes = FALSE,
                           data = smry[, .(xpos, mean_y, sem_y, treatment)])
  if (!is.null(cld) && nrow(cld) > 0) {
    lab <- merge(smry, as.data.table(cld)[, .(treatment = .fm_treat_lvl(treatment),
                                              tp = as.integer(timepoint),
                                              label = trimws(.group))],
                 by = c("treatment", "tp"), all.x = TRUE)
    lab <- lab[!is.na(label) & nzchar(label)]
    if (nrow(lab)) {
      lab[, ylab := mean_y + sem_y + 0.10 * yrange]
      p <- p + ggplot2::geom_text(
        data = lab,
        ggplot2::aes(x = xpos, y = ylab, label = label),
        inherit.aes = FALSE, colour = "black", fontface = "bold",
        size = 5.2, hjust = 0.5, vjust = 0)
      ymax <- max(ymax, max(lab$ylab, na.rm = TRUE))
    }
  }
  p +
    ggplot2::scale_colour_manual(values = FM_TREAT_COLORS, name = "Treatment") +
    ggplot2::scale_shape_manual(values = FM_TREAT_SHAPES, name = "Treatment") +
    ggplot2::scale_linetype_manual(values = FM_TREAT_LINES, name = "Treatment") +
    ggplot2::scale_x_continuous(breaks = 1:3,
                                labels = c("5–25", "45–65", "85–105"),
                                limits = c(0.6, 3.4)) +
    ggplot2::expand_limits(y = ymax + 0.06 * yrange) +
    ggplot2::labs(x = "Interval (min)", y = y_label, caption = cap) +
    theme_fm(13)
}

fig11_temporal_structure <- function(seq_dir = .fm_latest_run("SEQ"),
                                      bout_dir = .fm_latest_run("BOUT"),
                                      out_dir = .FM_OUT) {
  stopifnot(!is.na(seq_dir), dir.exists(seq_dir), !is.na(bout_dir), dir.exists(bout_dir))
  made <- character(0)

  av_seq_trial <- .fm_read(seq_dir, "seq_anova_trial_level.csv", required = FALSE)
  av_seq_tp    <- .fm_read(seq_dir, "seq_anova_trial_x_timepoint.csv", required = FALSE)
  av_bout_trial<- .fm_read(bout_dir, "bout_anova_trial_level.csv", required = FALSE)
  av_bout_tp   <- .fm_read(bout_dir, "bout_anova_trial_x_timepoint.csv", required = FALSE)
  seq_bin      <- .fm_read(seq_dir, "seq_metrics_per_trial_binary.csv", required = FALSE)
  bout_sch     <- .fm_read(bout_dir, "bout_metrics_per_school.csv", required = FALSE)
  bout_sess    <- .fm_read(bout_dir, "bout_metrics_per_session.csv", required = FALSE)

  # ---- A, D: bout metrics, Level A (trial) ----------------------------------
  .bout_panel_A <- function(met, y_label) {
    d <- as.data.table(bout_sch)[is.finite(get(met)), .(val = get(met), treatment)]
    r <- as.data.table(av_bout_trial)[metric == met & term == "treatment"]
    cap <- if (nrow(r)) paste0("Treatment ", .fm_cap(r$df1[1], r$df2[1], r$F[1],
                                                      r$p[1], r$eta2_p[1])) else ""
    sig <- nrow(r) > 0 && is.finite(r$p[1]) && r$p[1] < 0.05
    .fm_treatment_panel(d, y_label = y_label, caption = cap, is_sig = sig)
  }
  .fig11_D <- .bout_panel_A("memory_flow", "Flow bout memory")

  # ---- E: burst_flow, Level B (trial x interval) ----------------------------
  .fig11_E <- .fm_interval_panel(
    bout_sess, "burst_flow", "Burstiness of flow bouts",
    .fm_cap_levelB(av_bout_tp, "burst_flow"))

  # ---- F: dwell_Calm, Level B (trial x interval) -----------------------------
  # seq_metrics_per_trial_binary.csv IS the session-level (n=48) table --
  # "trial" here means physical school, one row per (trial, timepoint).
  # CLD letters come from cld_dwell_Calm_treatmentxtimepoint.csv, computed by
  # 00_shared/make_clds-style refit (emmeans on the same Level-B model, letters
  # remapped so control/interval 1 carries "a"). All six cells share "a": the
  # treatment x interval interaction is significant but no individual cell pair
  # survives the multiplicity adjustment.
  .cld_dwell_calm <- local({
    f <- file.path(file.path(PROJECT_ROOT, "choice R pipeline/easy_scripts/standalone"),
                   "cld_dwell_Calm_treatmentxtimepoint.csv")
    if (file.exists(f)) readr::read_csv(f, show_col_types = FALSE) else NULL
  })
  .fig11_F <- .fm_interval_panel(
    seq_bin, "dwell_Calm", "Mean calm dwell time (s)",
    .fm_cap_levelB(av_seq_tp, "dwell_Calm", alph_filter = "binary"),
    cld = .cld_dwell_calm, dodge = 0.14)

  # 2026-08-10: trial-level burst_flow / dwell_Flow / dwell_Calm panels dropped
  # (former A, B, C); the surviving three are promoted D->A, E->B, F->C.
  panels <- list(A = .fig11_D, B = .fig11_E, C = .fig11_F)
  if (any(vapply(panels, is.null, logical(1)))) {
    message("Figure 11 skipped -- missing component panel(s): ",
            paste(names(panels)[vapply(panels, is.null, logical(1))], collapse = ", "))
    return(list(made = made, out_dir = out_dir))
  }

  # Single shared legend: patchwork's automatic guide-collection can't see
  # inside wrap_elements(full=...) (it grob-ifies the whole rendered panel,
  # legend included, before composition -- the same reason the STEP5 script
  # embeds a separate legend PNG rather than relying on guides = "collect").
  # cowplot::get_legend() is the usual workaround but segfaults against this
  # ggplot2 build. Simplest robust fix, no extra dependency: keep the legend
  # on panel F only (bottom-right) and strip it everywhere else.
  # Caption size reduced from theme_fm's default (base_size - 2 = 11 pt): these
  # panels are composed 3-across at 320 mm (~106 mm/panel), narrower than this
  # module's standalone 140 mm single-panel exports the default was tuned for.
  # A-D (one-line, .fm_cap) fit at 9.5 pt; E-F (two-line, .fm_cap_levelB, up to
  # ~58 chars on line 2) needed 7.5 pt -- both confirmed by rendering and
  # visually inspecting the actual PNG (11 pt clipped eta2p off every panel).
  tag_theme <- function(p, tag, keep_legend, cap_size) {
    p <- p + ggplot2::labs(tag = tag) +
      ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 16),
                     plot.tag.position = "topleft",
                     plot.margin = ggplot2::margin(8, 10, 8, 10),
                     plot.caption = ggtext::element_markdown(
                       size = cap_size, hjust = 0, face = "italic",
                       colour = "black", lineheight = 1.3,
                       margin = ggplot2::margin(t = 8)))
    if (!keep_legend) p <- p + ggplot2::theme(legend.position = "none")
    p
  }
  # All three panels keep their guides so patchwork can collect them into a
  # SINGLE legend centred at the bottom of the figure (plot_layout(guides =
  # "collect")). This replaces the old "legend on panel F only" workaround --
  # that left panel F's plotting area visibly smaller than the others, because
  # the in-panel legend ate its vertical space. Guide collection works here
  # only because the panels are composed directly: wrap_elements(full = ...)
  # grob-ifies each panel (legend included) before composition, which is what
  # made collection impossible in the previous 3x2 layout.
  # Panel A's colour/shape merely restate its x axis, so its guide is dropped;
  # B and C carry identical colour/shape/linetype scales and merge cleanly.
  tagged <- Map(tag_theme, panels, names(panels),
               keep_legend = TRUE,
               cap_size = c(A = 9.5, B = 7.5, C = 7.5)[names(panels)])
  tagged$A <- tagged$A +
    ggplot2::guides(colour = "none", shape = "none", linetype = "none", fill = "none")

  # guide_area() gives the collected legend its own full-width row at the foot
  # of the figure, so it is centred on the FIGURE rather than on whichever
  # panel happened to carry it. Using `& theme(legend.position = "bottom")`
  # instead would reserve legend space inside every panel as well, which pushes
  # the per-panel captions away from their plots and leaves a large gap.
  fig <- ((tagged$A | tagged$B | tagged$C) / patchwork::guide_area()) +
    patchwork::plot_layout(guides = "collect", heights = c(1, 0.08)) +
    patchwork::plot_annotation(
      theme = ggplot2::theme(
        plot.background = ggplot2::element_rect(fill = "white", colour = NA),
        plot.margin     = ggplot2::margin(5, 5, 5, 5, unit = "mm"))) &
    ggplot2::theme(legend.direction = "horizontal",
                   legend.justification = "center",
                   legend.box = "horizontal",
                   legend.box.just = "center")

  if (isTRUE(.fm_save(fig, out_dir, "figure11_temporal_structure", 320, 125)))
    made <- c(made, "figure11_temporal_structure")

  list(made = made, out_dir = out_dir)
}

# =============================================================================
# ENTRY POINT
# =============================================================================
fig_module_run <- function(which = c("both", "sequence", "pattern"),
                           seq_dir = .fm_latest_run("SEQ"),
                           bout_dir = .fm_latest_run("BOUT"),
                           out_dir = .FM_OUT) {
  which <- match.arg(which)
  res <- list()
  if (which %in% c("both", "sequence")) res$sequence <- fig_seq_all(seq_dir, out_dir)
  if (which %in% c("both", "sequence")) res$manuscript_embeds <- fig_seq_manuscript_embeds(seq_dir, out_dir)
  if (which %in% c("both", "pattern"))  res$pattern  <- fig_bout_all(bout_dir, out_dir)
  # Figure 11 (temporal structure) needs both engines' exports regardless of
  # `which`, so only runs for "both". Manuscript-cited as "Figure 11"; saved
  # under a distinct filename ("figure11_temporal_structure", NOT
  # "Figure_11") because the STEP5 script's own copy-rename step already
  # claims "Figure_11.png" internally for the (manuscript-numbered) Figure 13
  # serotonin panel -- see activity_analysis_STATS_choice_exp.R's
  # DM_serotonin_group_SEM -> Figure_11 rename.
  if (which == "both") res$temporal_structure <- fig11_temporal_structure(seq_dir, bout_dir, out_dir)

  message("\n=== figure module summary ===")
  for (nm in names(res)) {
    r <- res[[nm]]
    message(sprintf("%-9s %2d figure(s) -> %s", nm, length(r$made), r$out_dir))
    if (length(r$skipped))
      message("          not rebuilt: ", paste(r$skipped, collapse = "; "))
  }
  invisible(res)
}

if (identical(environment(), globalenv()) &&
    !interactive() && sys.nframe() == 0L) {
  fig_module_run()
}
