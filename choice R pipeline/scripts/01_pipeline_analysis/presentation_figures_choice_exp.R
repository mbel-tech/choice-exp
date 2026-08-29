# =============================================================================
# PRESENTATION FIGURES — clean single-panel images of each SIGNIFICANT result
# =============================================================================
# For talks: one panel per significant result, every element enlarged so it
# reads from across a room. Significance is judged on RAW p (never BH/q).
#
# Self-contained: reads already-saved engine outputs; does not re-run any of
# the three source engines. Retains every visual element of the panels already
# coded for (behaviour scatters + Dm jitter plots), only scaled up. The three
# behavioural-sequence visuals are new (drafts for feedback).
#
# SOURCES
#   Behaviour : easy_scripts/easy_scripts_dataset.csv  (aggregated to phys trial)
#               captions from easy_scripts/standalone/anova_statement_table_trial.csv
#   Endocrine : choice exp .../claude output/Data/monoamine_long.csv  (Dm DA/DOPAC)
#               captions from .../Data/cell_cross_summary_bh.csv  (p_raw column)
#   Sequence  : latest output/SEQ_output/SEQ_output_*/  (metrics + pooled matrix)
#               ribbon states recomputed from latest STEP1 master_fish_by_frame.csv
#
# OUTPUT  ->  D:/CHOICE R SCRIPTS/presentation figure/   (PNG only, large fonts)
#
# Run: "C:/Program Files/R/R-4.6.0/bin/Rscript.exe" --vanilla <this file>
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr)
  library(ggplot2); library(ggsignif); library(ggtext); library(stringr)
})

ts_msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), " | ", ..., "\n", sep = "")

# ---- PATHS ------------------------------------------------------------------
PIPE_ROOT <- file.path(PROJECT_ROOT, "choice R pipeline")
ENDO_DATA <- file.path(PROJECT_ROOT, "choice exp for claude mono and cortisol/claude output/Data")
OUT_DIR   <- file.path(PROJECT_ROOT, "presentation figure")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- SHARED AESTHETICS (verbatim from both engines) -------------------------
# Okabe-Ito: control = index 6 (#D55E00), exercise choice = index 5 (#0072B2).
PAL <- c("control" = "#D55E00", "exercise choice" = "#0072B2")
SHP <- c("control" = 17L,       "exercise choice" = 16L)
# Title-case aliases (Dm panels label treatments "Control"/"Exercise choice").
PAL_TC <- c("Control" = "#D55E00", "Exercise choice" = "#0072B2")

# Presentation scale — everything big.
P_JITTER <- 0.12; P_ALPHA <- 0.9; P_SIZE <- 5.0
P_MEANW  <- 0.42; P_LWMEAN <- 1.6; P_ERRW <- 0.16; P_LWERR <- 1.3
STAR_TS  <- 13;   SIG_LW <- 0.9;   N_SIZE <- 7.0
DOT_W_IN <- 8;    DOT_H_IN <- 6;   DPI <- 300     # 4:3 panels @ 2400x1800

theme_pres <- function(base = 24) {
  ggplot2::theme_minimal(base_size = base) +
    ggplot2::theme(
      panel.grid       = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      axis.line.x      = ggplot2::element_line(colour = "black", linewidth = 1.3),
      axis.line.y      = ggplot2::element_line(colour = "black", linewidth = 1.3),
      axis.ticks       = ggplot2::element_line(colour = "black", linewidth = 1.0),
      axis.ticks.length = ggplot2::unit(6, "pt"),
      axis.text.x      = ggplot2::element_text(size = base + 2, face = "bold", colour = "black"),
      axis.text.y      = ggplot2::element_text(size = base,     colour = "black"),
      axis.title.y     = ggplot2::element_text(size = base + 4, face = "bold", colour = "black"),
      axis.title.x     = ggplot2::element_blank(),
      legend.position  = "none",
      plot.margin      = ggplot2::margin(14, 18, 14, 14),
      plot.caption     = ggtext::element_markdown(size = base, hjust = 0.5,
                                                  face = "italic", lineheight = 1.25,
                                                  margin = ggplot2::margin(t = 12))
    )
}

.sem  <- function(x) { x <- x[is.finite(x)]; if (length(x) < 2) NA_real_ else sd(x) / sqrt(length(x)) }
.save <- function(p, name, w = DOT_W_IN, h = DOT_H_IN) {
  if (is.null(p)) { ts_msg("  SKIP (NULL): ", name); return(invisible()) }
  f <- file.path(OUT_DIR, paste0(name, ".png"))
  tryCatch({ ggplot2::ggsave(f, p, width = w, height = h, units = "in", dpi = DPI, bg = "white")
             ts_msg("  wrote ", basename(f)) },
           error = function(e) ts_msg("  FAILED ", name, ": ", conditionMessage(e)))
}

# Two-line x labels matching the manuscript panels.
X_LABELS <- c("control" = "Control", "exercise choice" = "Exercise\nchoice")

# ---- Add per-group (N = X) labels below the data (retained element) ---------
add_n_labels <- function(p, dat, ycol, unit_col = "phys_trial_id", paren = TRUE) {
  nd <- dat %>% dplyr::group_by(treatment) %>%
    dplyr::summarise(n = if (unit_col %in% names(dat))
                           dplyr::n_distinct(.data[[unit_col]]) else dplyr::n(),
                     .groups = "drop") %>%
    dplyr::mutate(.x = as.integer(factor(treatment, levels = names(X_LABELS))),
                  label = if (paren) paste0("(N = ", n, ")") else paste0("N = ", n))
  raw <- dat[[ycol]][is.finite(dat[[ycol]])]
  ymin <- min(raw, na.rm = TRUE); yr <- max(diff(range(raw, na.rm = TRUE)), 1e-10)
  if (ymin >= 0) { y_n <- -0.08 * yr; y_exp <- -0.18 * yr }
  else           { y_n <- ymin - 0.10 * yr; y_exp <- ymin - 0.20 * yr }
  nd$lbl_y <- y_n
  p + ggplot2::expand_limits(y = y_exp) +
    ggplot2::geom_text(data = nd, ggplot2::aes(x = .x, y = lbl_y, label = label),
                       inherit.aes = FALSE, size = N_SIZE, colour = "grey35", vjust = 0.5)
}

# ---- Significance bracket + "*" (manual ggsignif, scaled) --------------------
add_sig_bracket <- function(p, yvals) {
  ym <- max(yvals, na.rm = TRUE); rng <- max(diff(range(yvals, na.rm = TRUE)), 1e-10)
  bary <- ym + 0.10 * rng
  sig_df <- data.frame(xmin = 1, xmax = 2, y_position = bary, annotations = "*")
  p + ggplot2::expand_limits(y = ym + 0.24 * rng) +
    ggsignif::geom_signif(data = sig_df,
      mapping = ggplot2::aes(xmin = xmin, xmax = xmax, y_position = y_position,
                             annotations = annotations),
      manual = TRUE, tip_length = 0.02, size = SIG_LW, textsize = STAR_TS,
      vjust = 0.5, colour = "black", inherit.aes = FALSE)
}

# =============================================================================
# GENERIC DOT PANEL — behaviour style (jitter + black mean CROSSBAR + SEM +
# N labels + optional dashed zero ref + sig bracket + raw-p caption).
# =============================================================================
build_beh_panel <- function(dat, ycol, ylab, cap_txt, is_sig,
                            dashed_zero = FALSE, nonneg_breaks = TRUE) {
  d <- dat[is.finite(dat[[ycol]]), ]
  if (nrow(d) == 0) return(NULL)
  d$.x <- as.integer(factor(d$treatment, levels = names(X_LABELS)))
  smry <- d %>% dplyr::group_by(treatment) %>%
    dplyr::summarise(mean_y = mean(.data[[ycol]], na.rm = TRUE),
                     se_y = .sem(.data[[ycol]]), .groups = "drop") %>%
    dplyr::mutate(.x = as.integer(factor(treatment, levels = names(X_LABELS))))

  p <- ggplot2::ggplot(d, ggplot2::aes(x = .x, y = .data[[ycol]],
                                       colour = treatment, shape = treatment))
  if (dashed_zero)
    p <- p + ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                                 colour = "grey55", linewidth = 0.9)
  p <- p +
    ggplot2::geom_point(position = ggplot2::position_jitter(width = P_JITTER, height = 0),
                        alpha = P_ALPHA, size = P_SIZE) +
    ggplot2::geom_crossbar(data = smry,
      ggplot2::aes(x = .x, y = mean_y, ymin = mean_y, ymax = mean_y),
      colour = "black", width = P_MEANW, linewidth = P_LWMEAN, show.legend = FALSE) +
    ggplot2::geom_errorbar(data = smry,
      ggplot2::aes(x = .x, y = mean_y, ymin = mean_y - se_y, ymax = mean_y + se_y),
      colour = "black", width = P_ERRW, linewidth = P_LWERR, show.legend = FALSE) +
    ggplot2::scale_x_continuous(breaks = c(1, 2), labels = unname(X_LABELS),
                                limits = c(0.5, 2.5), expand = c(0, 0)) +
    ggplot2::scale_colour_manual(values = PAL) +
    ggplot2::scale_shape_manual(values = SHP) +
    ggplot2::labs(y = ylab, caption = cap_txt) +
    theme_pres()
  if (isTRUE(nonneg_breaks))
    p <- p + ggplot2::scale_y_continuous(breaks = function(l) { b <- pretty(l); b[b >= 0] })
  if (isTRUE(is_sig)) p <- add_sig_bracket(p, d[[ycol]])
  p <- add_n_labels(p, d, ycol, unit_col = "phys_trial_id")
  p
}

# =============================================================================
# 1) BEHAVIOUR — 5 trial-level dot panels (raw-p significant)
# =============================================================================
ts_msg("=== Behaviour panels ===")
beh <- data.table::fread(file.path(PIPE_ROOT, "easy_scripts", "easy_scripts_dataset.csv"))
beh[, treatment := tolower(trimws(treatment))]
beh <- beh[treatment %in% names(PAL)]

# Aggregate session rows -> physical trial (N = 16; 8 per treatment), matching
# the aggregated manuscript panels.
beh_resps <- c("logit_flow", "lr_high", "mean_nnd_cm", "mean_polarisation",
               "mean_centroid_spd_cm")
beh_agg <- beh[, c(list(treatment = treatment[1]),
                   lapply(.SD, mean, na.rm = TRUE)),
               by = phys_trial_id, .SDcols = intersect(beh_resps, names(beh))]

# Captions: raw-p ANOVA statements straight from the manuscript CSV.
cap_tbl <- data.table::fread(
  file.path(PIPE_ROOT, "easy_scripts", "standalone", "anova_statement_table_trial.csv"),
  encoding = "UTF-8")
get_cap <- function(var) {
  r <- cap_tbl[variable == var]
  if (nrow(r) == 0) return("")
  paste0("*", r$statement[1], "*")
}

beh_specs <- list(
  list(col = "logit_flow",           lab = "Logit(flow occupancy)",  cap = "logit_flow",       zero = TRUE,  nn = FALSE),
  list(col = "lr_high",              lab = "High-flow log-ratio",    cap = "lr_high",          zero = TRUE,  nn = FALSE),
  list(col = "mean_nnd_cm",          lab = "Mean NND (cm)",          cap = "mean_nnd_cm",           zero = FALSE, nn = TRUE),
  list(col = "mean_polarisation",    lab = "Polarisation",           cap = "mean_polarisation",     zero = FALSE, nn = TRUE),
  list(col = "mean_centroid_spd_cm", lab = "School speed (cm/s)",    cap = "mean_school_speed_cm_s",zero = FALSE, nn = TRUE)
)
beh_names <- c(logit_flow = "beh_logit_flow", lr_high = "beh_lr_high",
               mean_nnd_cm = "beh_nnd", mean_polarisation = "beh_polarisation",
               mean_centroid_spd_cm = "beh_school_speed")
for (s in beh_specs) {
  p <- tryCatch(build_beh_panel(as.data.frame(beh_agg), s$col, s$lab, get_cap(s$cap),
                                is_sig = TRUE, dashed_zero = s$zero, nonneg_breaks = s$nn),
                error = function(e) { ts_msg("  beh ", s$col, " failed: ", conditionMessage(e)); NULL })
  .save(p, beh_names[[s$col]])
}

# =============================================================================
# 2) ENDOCRINE — Dm DA + DOPAC (retain make_mono_area_plot elements: jitter +
#    black mean POINT + SEM + "N = X" + sig bracket + raw-p F caption)
# =============================================================================
ts_msg("=== Dm endocrine panels ===")
mono <- data.table::fread(file.path(ENDO_DATA, "monoamine_long.csv"))
cell_bh <- data.table::fread(file.path(ENDO_DATA, "cell_cross_summary_bh.csv"))

build_dm_panel <- function(akey, analyte_lab, ylab) {
  d <- mono[analyte_key == akey & area == "DM" & is.finite(value) & value > 0 &
              !is.na(treatment)]
  if (nrow(d) == 0) return(NULL)
  d[, treat_lbl := factor(dplyr::recode(as.character(treatment),
                          "control" = "Control", "treat" = "Exercise choice"),
                          levels = c("Control", "Exercise choice"))]
  d[, .x := as.integer(treat_lbl)]
  smry <- d[, .(mean_y = mean(value, na.rm = TRUE), se_y = .sem(value),
                n = .N), by = treat_lbl]
  smry[, .x := as.integer(treat_lbl)]

  # Raw-p F caption (NOT p_BH) from cell_cross_summary_bh.
  row <- cell_bh[area == "DM" & analyte == analyte_lab]
  cap <- if (nrow(row))
    sprintf("*F<sub>%d,%d</sub> = %.2f, p = %s*",
            round(row$df1[1]), round(row$df2[1]), row$F_stat[1],
            ifelse(row$p_raw[1] < 0.001, "< 0.001", sprintf("%.3f", row$p_raw[1])))
  else ""
  is_sig <- nrow(row) && is.finite(row$p_raw[1]) && row$p_raw[1] < 0.05

  yr <- max(diff(range(d$value, na.rm = TRUE)), 1e-10); ymin <- min(d$value, na.rm = TRUE)
  y_n <- if (ymin >= 0) -0.02 * yr else ymin - 0.10 * yr

  p <- ggplot2::ggplot(d, ggplot2::aes(x = .x, y = value, colour = treat_lbl)) +
    ggplot2::geom_point(position = ggplot2::position_jitter(width = P_JITTER, height = 0),
                        size = P_SIZE, alpha = P_ALPHA) +
    ggplot2::geom_errorbar(data = smry,
      ggplot2::aes(x = .x, y = mean_y, ymin = mean_y - se_y, ymax = mean_y + se_y),
      width = P_ERRW, colour = "black", linewidth = P_LWERR, inherit.aes = FALSE) +
    ggplot2::geom_point(data = smry, ggplot2::aes(x = .x, y = mean_y),
                        colour = "black", size = P_SIZE + 1.5, inherit.aes = FALSE) +
    ggplot2::geom_text(data = smry, ggplot2::aes(x = .x, y = y_n, label = paste0("N = ", n)),
                       inherit.aes = FALSE, vjust = 1, size = N_SIZE) +
    ggplot2::scale_x_continuous(breaks = c(1, 2),
                                labels = c("Control", "Exercise\nchoice"),
                                limits = c(0.5, 2.5), expand = c(0, 0)) +
    ggplot2::scale_colour_manual(values = PAL_TC) +
    ggplot2::labs(y = ylab, caption = cap) +
    ggplot2::expand_limits(y = y_n - 0.06 * yr) +
    theme_pres()
  if (is_sig) p <- add_sig_bracket(p, d$value)
  p
}
.save(build_dm_panel("da",    "DA",    "[DA] (ng/mg)"),    "dm_da")
.save(build_dm_panel("dopac", "DOPAC", "[DOPAC] (ng/mg)"), "dm_dopac")

# =============================================================================
# 3) BEHAVIOURAL SEQUENCE — NEW visuals (drafts for feedback)
# =============================================================================
ts_msg("=== Sequence panels (new) ===")
seq_dir <- {
  parent <- file.path(PIPE_ROOT, "output", "SEQ_output")
  subs <- list.dirs(parent, full.names = TRUE, recursive = FALSE)
  subs <- subs[grepl("SEQ_output_\\d{8}_\\d{6}$", basename(subs))]
  if (length(subs)) subs[which.max(file.mtime(subs))] else NA_character_
}
# Ribbon state palette (user-specified): calm + the three flow bands.
STATE_FILL <- c("Calm" = "honeydew4", "Low" = "goldenrod3",
                "Medium" = "cyan3", "High" = "tomato3")

if (!is.na(seq_dir)) {
  seq_bin <- data.table::fread(file.path(seq_dir, "seq_metrics_per_trial_binary.csv"))
  seq_grd <- data.table::fread(file.path(seq_dir, "seq_metrics_per_trial_graded.csv"))
  seq_con <- data.table::fread(file.path(seq_dir, "treatment_contrasts.csv"))

  # LMM STATEMENT: Type III ANOVA F of the linear mixed model (Satterthwaite)
  # primary + Wilcoxon robust, on two separate lines, matching the italic in-panel
  # statistic (F<sub>df1,df2</sub> = …, p = …) used on the other figures.
  .fmt_p <- function(p) ifelse(!is.finite(p), "NA",
                       ifelse(p < 0.001, "< 0.001", sprintf("%.3f", p)))
  .fmt_df <- function(x) ifelse(!is.finite(x), "NA",
                        ifelse(abs(x - round(x)) < 0.05, sprintf("%d", as.integer(round(x))),
                               sprintf("%.1f", x)))
  seq_row <- function(alph, metr) {
    r <- seq_con[alphabet == alph & metric == metr]
    if (nrow(r)) r[1] else NULL
  }
  seq_cap <- function(r) {
    if (is.null(r)) return("")
    sprintf("*treatment: F<sub>%s,%s</sub> = %.2f, p = %s*<br>*Wilcoxon: W = %g, p = %s*",
            .fmt_df(r$df1_lmm[1]), .fmt_df(r$df2_lmm[1]), r$F_lmm[1], .fmt_p(r$p_lmm[1]),
            r$W[1], .fmt_p(r$p_wilcox[1]))
  }

  # --- 3a. Significant metric dot-panels (school as unit; 8 per treatment) ---
  build_seq_metric <- function(tbl, ycol, ylab, r) {
    d <- as.data.frame(tbl)[is.finite(as.data.frame(tbl)[[ycol]]), ]
    d$treatment <- tolower(trimws(d$treatment))
    d <- d[d$treatment %in% names(PAL), ]
    # collapse timepoints -> one value per physical school
    d <- d %>% dplyr::group_by(trial, treatment) %>%
      dplyr::summarise(val = mean(.data[[ycol]], na.rm = TRUE), .groups = "drop")
    names(d)[names(d) == "val"] <- ycol
    p_lmm <- if (is.null(r)) NA_real_ else r$p_lmm[1]
    build_beh_panel(as.data.frame(d), ycol, ylab, seq_cap(r),
                    is_sig = is.finite(p_lmm) && p_lmm < 0.05,
                    dashed_zero = FALSE,
                    nonneg_breaks = FALSE)
  }
  # NOTE: build_beh_panel counts N via phys_trial_id; sequence uses `trial`.
  #       Provide a phys_trial_id alias so N labels read "(N = 8)".
  add_alias <- function(df) { df$phys_trial_id <- df$trial; df }

  .save(build_seq_metric(add_alias(as.data.frame(seq_bin)), "switch_rate",
                         "Switches / min", seq_row("binary", "switch_rate")),
        "seq_switch_rate")
  .save(build_seq_metric(add_alias(as.data.frame(seq_bin)), "entropy_rate",
                         "Sequence entropy", seq_row("binary", "entropy_rate")),
        "seq_entropy")
  .save(build_seq_metric(add_alias(as.data.frame(seq_grd)), "occ_High",
                         "Time in High-flow", seq_row("graded", "occ_High")),
        "seq_occ_high")

  lv <- c("Calm", "Low", "Medium", "High")   # ribbon zone-state order (calm-inclusive)

  # --- 3b. Example state-ribbon timeline (recompute bin states from master) --
  # Hand-picked representative sessions (per user feedback): the control school
  # that switches MOST vs the exercise-choice school that switches LEAST (most
  # committed), chosen by switch_rate across any timepoint. Labelled honestly
  # as representative examples.
  ribbon_ok <- FALSE
  master <- {
    parent <- file.path(PIPE_ROOT, "output", "STEP1_output")
    subs <- list.dirs(parent, full.names = TRUE, recursive = FALSE)
    subs <- subs[grepl("STEP1_output_\\d{8}_\\d{6}$", basename(subs))]
    if (length(subs)) file.path(subs[which.max(file.mtime(subs))], "master_fish_by_frame.csv")
    else NA_character_
  }
  if (!is.na(master) && file.exists(master)) {
    # Choose example sessions by switch rate (any timepoint), from the SEQ
    # per-session metrics: control = MOST switching, exercise = LEAST switching.
    sb <- as.data.frame(seq_bin)
    sb$treatment <- tolower(trimws(sb$treatment))
    sb <- sb[is.finite(sb$switch_rate) & sb$treatment %in% names(PAL), ]
    sb_c <- sb[sb$treatment == "control", ]
    sb_e <- sb[sb$treatment == "exercise choice", ]
    ctrl_pick <- sb_c[which.max(sb_c$switch_rate), c("trial", "timepoint")]
    exer_pick <- sb_e[which.min(sb_e$switch_rate), c("trial", "timepoint")]

    # Stream the master in chunks (memory-safe: avoids fread's whole-file memory
    # map, which fails on very large masters). Each 1-s bin is classified by the
    # DOMINANT sub-zone the school occupied, now INCLUDING the calm zone, using
    # the `sec_zone` field (calm / low / medium / high). Only the two example
    # sessions are retained per chunk, so accumulated memory stays tiny.
    ts_msg("  ribbon: streaming master (", basename(dirname(master)), ") ...")
    zones4 <- c("calm", "low", "medium", "high")
    keep <- list(c(ctrl_pick$trial, ctrl_pick$timepoint),
                 c(exer_pick$trial, exer_pick$timepoint))
    acc <- list()
    .cb <- function(df, pos) {
      dt <- data.table::as.data.table(df)
      dt <- dt[sec_zone %in% zones4 &
               ((trial == keep[[1]][1] & timepoint == keep[[1]][2]) |
                (trial == keep[[2]][1] & timepoint == keep[[2]][2]))]
      if (!nrow(dt)) return(invisible())
      dt[, treatment := tolower(trimws(treatment))]
      dt[, bin := floor(as.numeric(time) / 1.0)]
      acc[[length(acc) + 1L]] <<- dt[, .N,
        by = .(trial, timepoint, treatment, bin, sec_zone)]
    }
    readr::read_csv_chunked(
      master, callback = readr::SideEffectChunkCallback$new(.cb),
      chunk_size = 500000L, progress = FALSE,
      col_types = readr::cols(.default = readr::col_character()))
    cnt <- data.table::rbindlist(acc)[, .(N = sum(N)),
             by = .(trial, timepoint, treatment, bin, sec_zone)]
    # Dominant sub-zone per bin (ties break toward calm via factor order).
    cnt[, sec_zone := factor(sec_zone, levels = zones4)]
    bo <- cnt[order(-N, sec_zone)][, .(sec_zone = sec_zone[1]),
             by = .(trial, timepoint, treatment, bin)]
    bo[, state := factor(c(calm = "Calm", low = "Low", medium = "Medium",
                           high = "High")[as.character(sec_zone)], levels = lv)]
    # Keep only the two hand-picked example sessions.
    ex <- bo[(trial == ctrl_pick$trial & timepoint == ctrl_pick$timepoint) |
             (trial == exer_pick$trial & timepoint == exer_pick$timepoint)]
    ex[, bin0 := bin - min(bin), by = .(trial, timepoint)]   # align each to t=0
    ex[, lbl := fifelse(treatment == "control",
                        "Control\n(high-switching example)",
                        "Exercise choice\n(committed example)")]
    ex[, lbl := factor(lbl, levels = c("Control\n(high-switching example)",
                                       "Exercise choice\n(committed example)"))]
    p_rib <- ggplot2::ggplot(ex, ggplot2::aes(x = bin0, y = 1, fill = state)) +
      ggplot2::geom_tile() +
      ggplot2::facet_wrap(~ lbl, ncol = 1, strip.position = "left") +
      ggplot2::scale_fill_manual(values = STATE_FILL, name = "Collective zone state", drop = FALSE) +
      ggplot2::scale_y_continuous(breaks = NULL) +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.005, 0.02))) +
      ggplot2::labs(x = "Time (s)", y = NULL) +
      ggplot2::theme_minimal(base_size = 24) +
      ggplot2::theme(strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 1, size = 20, face = "bold"),
                     axis.text.x = ggplot2::element_text(size = 22, colour = "black"),
                     axis.title.x = ggplot2::element_text(size = 24, face = "bold"),
                     panel.grid = ggplot2::element_blank(),
                     plot.margin = ggplot2::margin(10, 26, 10, 10),
                     legend.position = "bottom",
                     legend.text = ggplot2::element_text(size = 20),
                     legend.title = ggplot2::element_text(size = 22, face = "bold"))
    .save(p_rib, "seq_state_ribbon", w = 14, h = 5)
    ts_msg("  ribbon examples: control school ", ctrl_pick$trial, " tp", ctrl_pick$timepoint,
           " (switch rate ", round(max(sb_c$switch_rate), 2), "/min) vs exercise school ",
           exer_pick$trial, " tp", exer_pick$timepoint,
           " (switch rate ", round(min(sb_e$switch_rate), 2), "/min)")
    ribbon_ok <- TRUE
  }
  if (!ribbon_ok) ts_msg("  ribbon SKIPPED: master_fish_by_frame.csv not found")
} else {
  ts_msg("  Sequence panels SKIPPED: no SEQ_output folder found")
}

ts_msg("==================== PRESENTATION FIGURES DONE ====================")
ts_msg("Output: ", OUT_DIR)
cat("\nPanels: 5 behaviour (beh_*), 2 Dm endocrine (dm_*), ",
    "3 sequence (seq_switch_rate/entropy/occ_high + seq_state_ribbon). ",
    "All significance by RAW p; captions never show q.\n", sep = "")
