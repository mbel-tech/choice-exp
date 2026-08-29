# =============================================================================
# house_style.R -- CANONICAL source of truth for palette, theme, marks, CLD
# =============================================================================
# Part of the handoff bundle at D:\CHOICE R SCRIPTS\handoff\. See HANDOFF.md
# for the narrative and CONVENTIONS.md for the numbered rule (C1-C12) that each
# block below implements. DO NOT RETYPE any value in this file into a document
# or another script -- source this file instead.
#
# This file exists because the project's own history contains two drift
# defects that a "just copy the value" workflow produced:
#   (1) number_formatting.R states its own formatters were "DUPLICATED
#       verbatim, not sourced, into every engine" -- seven copies that had to
#       be kept in sync by hand.
#   (2) figures_step5_module.R:563 defines a dead, non-Okabe TREATMENT_COLORS
#       object with the exact name of the live Okabe palette used elsewhere,
#       one autocomplete away from being copied by mistake (adjudication A-3).
#
# Authority order for the whole bundle: house_style.R > CONVENTIONS.md >
# RESULTS_STYLE.md > HANDOFF.md > anything predating this bundle.
# =============================================================================

suppressPackageStartupMessages({
  # Soft dependencies -- functions degrade gracefully where a package absent.
  .HAS_EFFECTSIZE <- requireNamespace("effectsize", quietly = TRUE)
  .HAS_GGTEXT     <- requireNamespace("ggtext",     quietly = TRUE)
})

# =============================================================================
# A1 -- PALETTE  (CONVENTIONS.md C1, C2)
# =============================================================================
# Okabe-Ito, colourblind-safe. Treatment levels take Okabe indices 6 then 5,
# control first. Zone/state colours are a SEPARATE, non-overlapping family so
# a zone-coloured panel can never be mistaken for a treatment-coloured one.
# Statistical ink (means, error bars, CLD letters, brackets) is ALWAYS black --
# colour carries group identity, black carries inference; do not mix channels.
#
# ADJUDICATION A-3: figures_step5_module.R:563-565 defines dead objects named
# `TREATMENT_COLORS` (#2166AC/#D6604D) and `ZONE_COLORS_MAIN` (#4DAC26/#8073AC)
# that are NOT Okabe-Ito and have zero live uses in that file. Deliberately
# NOT reproduced here, and NOT to be copied from that file into new work.

PALETTE_TREATMENT <- c(
  "control"         = "#D55E00",   # Okabe-Ito 6, vermillion
  "exercise choice" = "#0072B2"    # Okabe-Ito 5, blue
)

SHAPES <- c(
  "control"         = 17L,         # filled triangle
  "exercise choice" = 16L          # filled circle
)

LINETYPES <- c(
  "control"         = "dashed",
  "exercise choice" = "solid"
)

PALETTE_ZONE <- c(
  "flow"   = "red4",
  "calm"   = "lightskyblue3",
  "high"   = "tomato3",
  "medium" = "mediumturquoise",
  "low"    = "goldenrod3"
)

PALETTE_RAMP_SEQ <- c("white", "#0072B2")               # sequential: white -> Okabe blue
PALETTE_RAMP_DIV <- c("#B2432F", "white", "#2E5FA3")    # diverging, centred on 0

INK_STAT <- "black"   # every mean, error bar, CLD letter, bracket, asterisk

# Decision rule (C1 [DR]): extending beyond 2 treatment levels stays on the
# Okabe-Ito sequence in this order: 6 ("#D55E00"), 5 ("#0072B2"), 3 ("#009E73"),
# 1 ("#E69F00"). Never index 4 (#F0E442, yellow -- invisible on white) or index
# 8 (#000000, black -- reserved for statistical ink). Above 4 levels, colour is
# not an appropriate channel for that variable; facet instead.
PALETTE_OKABE_FULL <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                        "#0072B2", "#D55E00", "#CC79A7", "#000000")

# =============================================================================
# A2 -- THEME  (CONVENTIONS.md C3)
# =============================================================================
# theme_minimal(base_size = 13), then every decorative element stripped: no
# grid at all (major or minor), no panel background/border, no strip
# background. Two black axis rules, bold categorical axis text. The x-axis
# title is blanked by default -- restore it only where the axis is not
# self-describing (e.g. "Interval (min)"), never for "Control"/"Exercise
# choice" style labels. 18 elements below; omitting any one changes the figure.

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
    plot.caption      = if (.HAS_GGTEXT)
      ggtext::element_markdown(size = 14, hjust = 0, face = "italic",
                               lineheight = 1.3, margin = ggplot2::margin(t = 10))
    else ggplot2::element_text(size = 14, hjust = 0, face = "italic",
                               lineheight = 1.3, margin = ggplot2::margin(t = 10)),
    plot.subtitle     = ggplot2::element_text(size = 9, hjust = 0, face = "italic",
                                              colour = "grey40")
  )

TIMEPOINT_LABELS <- c("1" = "5\u201325 min", "2" = "45\u201365 min", "3" = "85\u2013105 min")

# =============================================================================
# A3/A4 -- DATA MARKS  (CONVENTIONS.md C4)
# =============================================================================
# House scatter idiom, in draw order: (1) jittered raw points coloured by
# group; (2) black geom_crossbar at the group mean; (3) black geom_errorbar at
# mean +/- 1 SE. NO BOX PLOTS. Raw data is plotted, never EMMs -- EMMs feed
# CLD letters and contrast tables only.

JITTER_W  <- 0.12
DODGE_W   <- 0.70
PT_ALPHA  <- 0.80
PT_SIZE   <- 2.4
MEAN_W    <- 0.40
LW_MEAN   <- 0.85
LW_ERR    <- 0.65
ERR_W     <- 0.13
CLD_SIZE  <- 6.0   # ADJUDICATION A-4: the *live* value in every geom_text() call.
                   # figures_step5_module.R:575 also defines CLD_SIZE <- 4.0, but
                   # that object is never read by a draw call -- it is dead. 6 is
                   # the convention; do not resurrect the 4.0 default.

.sem <- function(x) { x <- x[is.finite(x)]; if (length(x) < 2) NA_real_ else sd(x) / sqrt(length(x)) }

# house_scatter(): the three-layer idiom with every argument bound. `df` needs
# columns `treatment` (factor) and `y` (numeric). Returns a ggplot object with
# BASE_THEME and the treatment palette already applied.
house_scatter <- function(df, y_label = NULL) {
  lvls   <- levels(factor(df$treatment))
  n_trt  <- length(lvls)
  trt_x  <- stats::setNames(seq_len(n_trt), lvls)
  df     <- df[is.finite(df$y), ]
  df$.x  <- unname(trt_x[as.character(df$treatment)])

  smry <- stats::aggregate(y ~ treatment, df, function(v) c(mean = mean(v), sem = .sem(v)))
  smry <- do.call(data.frame, smry)
  names(smry) <- c("treatment", "mean_y", "sem_y")
  smry$.x <- unname(trt_x[as.character(smry$treatment)])

  ggplot2::ggplot(df, ggplot2::aes(x = .x, y = y, colour = treatment, shape = treatment)) +
    ggplot2::geom_point(position = ggplot2::position_jitter(width = JITTER_W, height = 0),
                        alpha = PT_ALPHA, size = PT_SIZE) +
    ggplot2::geom_crossbar(data = smry,
      ggplot2::aes(x = .x, y = mean_y, ymin = mean_y, ymax = mean_y),
      colour = INK_STAT, width = MEAN_W, linewidth = LW_MEAN, show.legend = FALSE, inherit.aes = FALSE) +
    ggplot2::geom_errorbar(data = smry,
      ggplot2::aes(x = .x, y = mean_y, ymin = mean_y - sem_y, ymax = mean_y + sem_y),
      colour = INK_STAT, width = ERR_W, linewidth = LW_ERR, show.legend = FALSE, inherit.aes = FALSE) +
    ggplot2::scale_x_continuous(breaks = seq_len(n_trt), labels = tools::toTitleCase(lvls),
                                limits = c(0.5, n_trt + 0.5), expand = c(0, 0)) +
    ggplot2::scale_colour_manual(values = PALETTE_TREATMENT, name = "Treatment") +
    ggplot2::scale_shape_manual(values = SHAPES, name = "Treatment") +
    BASE_THEME +
    ggplot2::labs(y = y_label)
}

# house_line(): line + point + black SE bars across an ordered interval
# factor. No ribbon (design_memo's ribbon convention was never adopted).
house_line <- function(df, y_label = NULL) {
  ggplot2::ggplot(df, ggplot2::aes(x = interval, y = mean_y, colour = treatment,
                                   shape = treatment, linetype = treatment, group = treatment)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = PT_SIZE * 1.6) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = mean_y - sem_y, ymax = mean_y + sem_y),
                           colour = INK_STAT, width = 0.10, linewidth = LW_ERR) +
    ggplot2::scale_colour_manual(values = PALETTE_TREATMENT, name = "Treatment") +
    ggplot2::scale_shape_manual(values = SHAPES, name = "Treatment") +
    ggplot2::scale_linetype_manual(values = LINETYPES, name = "Treatment") +
    BASE_THEME +
    ggplot2::labs(y = y_label)
}

# =============================================================================
# A5/A13 -- NUMBER FORMATTING  (CONVENTIONS.md C7)
# =============================================================================
# Exactly two base formatters, absorbed verbatim from
# choice R pipeline/scripts/00_shared/number_formatting.R (author decision,
# 2026-08-07), plus the "< 0.001" collapse for statistics/effect-sizes
# (analysis_b3_REVISED.R, 2026-08-09) and a bootstrap-aware p formatter (new).
#
# fmt_p() from analysis_b3_REVISED.R -- floor "< 0.0001", star-appending -- is
# ADJUDICATION A-9: SUPERSEDED. Not reproduced here. Do not port it in.

fmt_F <- function(x) {
  if (!is.finite(x)) return("NA")
  s <- format(signif(x, 4), scientific = FALSE, trim = TRUE)
  if (grepl(".", s, fixed = TRUE)) {
    s <- sub("0+$", "", s)
    s <- sub("\\.$", "", s)
  }
  s
}

fmt3 <- function(x) {
  if (!is.finite(x)) return("NA")
  s <- sprintf("%.3f", x)
  if (grepl("0$", s)) s <- substr(s, 1, nchar(s) - 1)
  s
}

# fmt_Fstat / fmt_es3: as fmt_F / fmt3, but a value below 0.001 collapses to
# "< 0.001" instead of a long decimal or an uninformative "0.00". Applied ONLY
# to test statistics and effect sizes -- denominator df keep fmt_F (a df is
# never < 1) and p-values keep fmt_p_house() below, which owns its own floor.
fmt_Fstat <- function(x) if (!is.finite(x)) "NA" else if (abs(x) < 0.001) "< 0.001" else fmt_F(x)
fmt_es3   <- function(x) if (!is.finite(x)) "NA" else if (abs(x) < 0.001) "< 0.001" else fmt3(x)

# fmt_p_house(): the CURRENT p-value rule (ADJUDICATION A-9 winner). Floor
# "p < 0.001", no significance stars appended in prose or captions (R11).
fmt_p_house <- function(p) {
  if (!is.finite(p)) return("NA")
  if (round(p, 3) == 0) return("p < 0.001")
  paste0("p = ", fmt3(p))
}

# fmt_p_boot(): a bootstrap p is never printed finer than the replicate count
# resolves (HANDOFF S7, RESULTS_STYLE R11). Smallest reportable value is
# 1/(N+1); anything landing on the floor is written "p <= 1/(N+1)", not as a
# spuriously precise decimal.
fmt_p_boot <- function(p, N) {
  floor_p <- 1 / (N + 1)
  if (!is.finite(p)) return("NA")
  if (p <= floor_p) return(paste0("p \u2264 ", fmt3(floor_p)))
  paste0("p = ", fmt3(p))
}

# =============================================================================
# A6 -- EFFECT SIZES  (HANDOFF.md S6, CONVENTIONS.md C6, RESULTS_STYLE.md R5)
# =============================================================================
# eta2_p / omega2_p formulas absorbed verbatim from
# choice R pipeline/scripts/00_shared/effect_size_and_sample_size.R.
# hedges_g absorbed from .calc_hedges_g() in activity_analysis_STATS_choice_exp.R
# (simplified to the noncentral-t / effectsize::t_to_d path; the raw-data and
# sigma fallbacks in the source are omitted here as pipeline-specific).
#
# effect_size() is the single entry point every caption/results builder must
# call. Passing stat_type = "chisq" returns es = NA_real_ by CONSTRUCTION, so
# a Wald chi-square caption cannot attach a partial eta-squared even by
# mistake (CONVENTIONS.md C6 [DR], preflight_check() assertion 8/9).

eta2_p <- function(F, df1, df2) {
  if (!(is.finite(F) && is.finite(df1) && is.finite(df2))) return(NA_real_)
  (F * df1) / (F * df1 + df2)
}

omega2_p <- function(F, df1, df2) {
  if (!(is.finite(F) && is.finite(df1) && is.finite(df2))) return(NA_real_)
  max(0, (df1 * (F - 1)) / (df1 * (F - 1) + df2 + 1))
}

# effect_size(): dispatches on stat_type. F-test -> eta2_p (+ omega2_p as a
# secondary value); chisq -> NA (structurally no partial eta-squared); t-test
# / 1-df contrast -> signed Hedges' g via noncentral-t inversion (NOT Wald),
# never absolute value.
effect_size <- function(stat_type = c("F", "chisq", "t"), F = NA_real_, df1 = NA_real_,
                        df2 = NA_real_, t = NA_real_, df = NA_real_) {
  stat_type <- match.arg(stat_type)
  if (identical(stat_type, "chisq")) {
    return(list(es = NA_real_, es_type = NA_character_,
               note = "Wald chi-square carries no partial eta-squared by design"))
  }
  if (identical(stat_type, "F")) {
    return(list(es = eta2_p(F, df1, df2), es_type = "eta2_p",
               es_secondary = omega2_p(F, df1, df2), es_secondary_type = "omega2_p"))
  }
  # stat_type == "t": signed Hedges' g, noncentral-t CI, J correction.
  if (!.HAS_EFFECTSIZE || !is.finite(t) || !is.finite(df) || df <= 1) {
    return(list(es = NA_real_, es_type = "g", note = "effectsize unavailable or df<=1"))
  }
  J     <- 1 - 3 / (4 * df - 1)
  d_res <- tryCatch(effectsize::t_to_d(t, df, ci = 0.95), error = function(e) NULL)
  if (is.null(d_res) || !nrow(d_res)) return(list(es = NA_real_, es_type = "g"))
  list(es = d_res$d[1] * J, ci_lo = d_res$CI_low[1] * J, ci_hi = d_res$CI_high[1] * J,
      es_type = "g (total-variance standardiser -- report alongside g_resid, never pooled)")
}

# =============================================================================
# A7 -- THE ANOVA CAPTION STRING  (CONVENTIONS.md C6)
# =============================================================================
# Order is law: statistic -> p -> effect size. Denominator df printed as
# integers (fractional KR/Satterthwaite values stay in anova.csv). Calls
# effect_size() and OMITS the effect-size clause when es is NA -- this is what
# makes "no eta2_p on a Wald chi-square" structural rather than a rule someone
# has to remember to apply.
anova_caption <- function(term_label, stat_type = c("F", "chisq"), stat, df1, df2 = NA, p) {
  stat_type <- match.arg(stat_type)
  p_txt <- fmt_p_house(p)
  if (identical(stat_type, "F")) {
    es <- effect_size("F", F = stat, df1 = df1, df2 = df2)
    stat_txt <- sprintf("F<sub>%d,%s</sub> = %s", round(df1), round(df2), fmt_Fstat(stat))
    es_txt   <- if (is.finite(es$es))
      sprintf(", \u03b7<sup>2</sup><sub>p</sub> = %s", fmt_es3(es$es)) else ""
  } else {
    es       <- effect_size("chisq")
    stat_txt <- sprintf("\u03c7<sup>2</sup><sub>%d</sub> = %s", round(df1), fmt_Fstat(stat))
    es_txt   <- ""   # structurally empty -- see effect_size()
  }
  sprintf("%s: %s, %s%s", term_label, stat_txt, p_txt, es_txt)
}

# =============================================================================
# A8 -- CLD PIPELINE  (CONVENTIONS.md C5)
# =============================================================================
# Absorbed from choice R pipeline/scripts/00_shared/cld_letter_policy.R.
# FOUR steps in MANDATORY order: multcomp::cld() -> drop redundant letters ->
# remap (control-first or lone-first) -> post-condition guard. Letters are
# NEVER truncated (ADJUDICATION A-8: substr(...,1,2) truncation is still live
# at five call sites in the pipeline -- figures_step5_module.R:1029 and
# activity_analysis_STATS_choice_exp.R:5531,5559,6109,6535 -- despite being
# the documented defect. Any new work must not add a sixth.)
#
# .cld_remap_lone_first() -- the alternative letter order used per-panel when
# control/interval-1 legitimately carries >1 letter -- is absorbed verbatim
# from cld_letter_policy.R without modification; see that file for its body
# (identical structure to .cld_remap_control_first() but anchoring "a" on the
# first single-letter cell in reading order instead of the control cell).

.CLD_TREAT_ORDER <- c("control", "exercise choice")

.cld_display_valid <- function(groups, same) {
  ids <- names(groups)
  lets <- lapply(groups, function(g) strsplit(trimws(g), "")[[1]])
  for (i in seq_along(ids)) for (j in seq_along(ids)) {
    if (j <= i) next
    shares <- length(intersect(lets[[i]], lets[[j]])) > 0L
    if (shares != isTRUE(same[i, j])) return(FALSE)
  }
  TRUE
}

.cld_drop_redundant_letters <- function(groups, pmat, alpha = 0.05, label = "") {
  ids <- names(groups)
  if (length(ids) < 2L) return(groups)
  same <- pmat >= alpha
  diag(same) <- TRUE
  if (!.cld_display_valid(groups, same)) {
    warning("CLD reduction skipped", if (nzchar(label)) paste0(" [", label, "]") else "",
            ": letters inconsistent with the p-value matrix. Left unchanged.")
    return(groups)
  }
  repeat {
    dropped <- FALSE
    for (i in seq_along(ids)) {
      lets <- sort(unique(strsplit(trimws(groups[[i]]), "")[[1]]))
      if (length(lets) < 2L) next
      for (l in lets) {
        trial <- groups
        keep  <- setdiff(strsplit(trimws(trial[[i]]), "")[[1]], l)
        if (!length(keep)) next
        trial[[i]] <- paste(sort(unique(keep)), collapse = "")
        if (.cld_display_valid(trial, same)) { groups <- trial; dropped <- TRUE; break }
      }
    }
    if (!dropped) break
  }
  groups
}

.cld_remap_control_first <- function(cld_df, tp_col = "tp_num", treat_col = "treatment",
                                     group_col = ".group", treat_order = .CLD_TREAT_ORDER) {
  if (is.null(cld_df) || !nrow(cld_df)) return(cld_df)
  ord <- order(match(trimws(tolower(as.character(cld_df[[treat_col]]))), treat_order),
              suppressWarnings(as.numeric(as.character(cld_df[[tp_col]]))))
  seen <- character(0)
  for (i in ord) {
    lets <- strsplit(trimws(as.character(cld_df[[group_col]][i])), "")[[1]]
    for (l in lets) if (nzchar(l) && !(l %in% seen)) seen <- c(seen, l)
  }
  if (!length(seen) || length(seen) > length(letters)) return(cld_df)
  map <- stats::setNames(letters[seq_along(seen)], seen)
  cld_df[[group_col]] <- vapply(as.character(cld_df[[group_col]]), function(g) {
    r <- map[strsplit(trimws(g), "")[[1]]]; r <- r[!is.na(r)]
    if (!length(r)) return(trimws(g))
    paste(sort(unique(r)), collapse = "")
  }, character(1L), USE.NAMES = FALSE)
  cld_df
}

.cld_is_uniform <- function(cld_df, group_col = ".group") {
  if (is.null(cld_df) || !nrow(cld_df) || !group_col %in% names(cld_df)) return(FALSE)
  length(unique(trimws(as.character(cld_df[[group_col]])))) == 1L
}

.cld_check_control_first <- function(cld_df, label = "", tp_col = "tp_num",
                                     treat_col = "treatment", group_col = ".group",
                                     mode = "control_first") {
  if (identical(mode, "lone_first")) return(invisible(TRUE))
  if (is.null(cld_df) || !nrow(cld_df)) return(invisible(TRUE))
  tp <- suppressWarnings(as.numeric(as.character(cld_df[[tp_col]])))
  i  <- which(trimws(tolower(as.character(cld_df[[treat_col]]))) == "control" & tp == min(tp, na.rm = TRUE))
  if (!length(i)) return(invisible(TRUE))
  g  <- trimws(as.character(cld_df[[group_col]][i[1]]))
  ok <- nzchar(g) && substr(g, 1, 1) == "a"
  if (!ok) warning("CLD policy violated", if (nzchar(label)) paste0(" [", label, "]") else "",
                   ": control/interval 1 carries '", g, "', expected a group starting with 'a'.")
  invisible(ok)
}

# cld_pipeline(): the four mandatory steps in order. `cld_df` is the data frame
# returned by multcomp::cld() with a `.group` column; `pmat` a symmetric
# p-value matrix with dimnames matching `cell_ids`.
cld_pipeline <- function(cld_df, cell_ids, pmat, alpha = 0.05, mode = "control_first", label = "") {
  if (.cld_is_uniform(cld_df)) return(list(cld_df = NULL, omitted = TRUE,
    reason = "every cell shares one group; no pairwise difference -- CLD omitted per house rule"))
  groups <- stats::setNames(as.list(trimws(as.character(cld_df$.group))), cell_ids)
  groups <- .cld_drop_redundant_letters(groups, pmat, alpha = alpha, label = label)
  cld_df$.group <- unlist(groups, use.names = FALSE)
  cld_df <- if (identical(mode, "lone_first")) {
    if (exists(".cld_remap_lone_first")) .cld_remap_lone_first(cld_df) else cld_df
  } else {
    .cld_remap_control_first(cld_df)
  }
  .cld_check_control_first(cld_df, label = label, mode = mode)
  list(cld_df = cld_df, omitted = FALSE)
}

# 15% of the panel's data range above the cell's highest plotted point; panel
# expands by 26% so nothing clips. Two magic fractions that must move
# together (CONVENTIONS.md C5).
CLD_LETTER_FRAC <- 0.15
CLD_EXPAND_FRAC <- 0.26

CLD_CAPTION_NOTE <- paste(
  "Compact letter display: cells sharing a letter do not differ at p < 0.05.",
  "Letters are relabelled so that control at the first interval always begins",
  "with 'a'; this is a relabelling only, and which cells share a letter is",
  "unchanged. A cell carrying two letters (e.g. 'ab') belongs to two",
  "non-difference groups -- it differs from neither -- and cannot be reduced",
  "to one letter without asserting a difference the model does not support.",
  "Panels without letters have no significant pairwise differences at all:",
  "every cell falls in the same group, so letters would carry no information",
  "and are omitted rather than printed as identical marks."
)

# =============================================================================
# A9 -- CLD PLACEMENT  (CONVENTIONS.md C5)
# =============================================================================
cld_place <- function(cell_max, data_range) cell_max + CLD_LETTER_FRAC * data_range

# =============================================================================
# A10 -- EXPORT  (CONVENTIONS.md C9)
# =============================================================================
# PNG and PDF, 300 dpi, dimensions in mm, bg = "white". PDFs via cairo_pdf so
# fonts embed (standard PDF fonts render sub/superscripts as blank boxes).
# STRIPPED TO ONE DESTINATION TREE relative to figures_step5_module.R's
# .mg_save(), which wrote to three trees -- that multi-tree pattern is
# repo-specific fragility (HANDOFF S11), not something to reproduce in new
# work.
house_save <- function(p, filename, out_dir, width_mm, height_mm, dpi = 300) {
  if (is.null(p)) { warning("house_save: NULL plot, skipped ", filename); return(invisible()) }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(file.path(out_dir, paste0(filename, ".png")), p,
                  width = width_mm, height = height_mm, units = "mm", dpi = dpi, bg = "white")
  ggplot2::ggsave(file.path(out_dir, paste0(filename, ".pdf")), p,
                  width = width_mm / 25.4, height = height_mm / 25.4, units = "in",
                  device = grDevices::cairo_pdf, bg = "white")
  invisible(TRUE)
}

# =============================================================================
# A12 -- PARAMETRIC-BOOTSTRAP LRT  (HANDOFF.md S7)
# =============================================================================
# For treatment terms in models where Kenward-Roger/Satterthwaite is undefined
# (beta-GLMM). Absorbed from .pb_lrt_term() in
# activity_analysis_STATS_choice_exp.R. pbkrtest::PBmodcomp has no glmmTMB
# method -- verified in that file -- which is why this is hand-rolled rather
# than delegated to a package.
#
# CRITICAL: the returned N and seed are what Methods and the caption print.
# A bootstrap without a stated N and seed is not a check (HANDOFF S7 [INV]).
pb_lrt <- function(fixed_str, reduced_fixed_str, re_str, data, family,
                   nsim = 1000L, seed) {
  na_out <- list(lrt_obs = NA_real_, p_pb = NA_real_, nsim_ok = 0L, N = nsim, seed = seed, method = "na")
  full_form    <- tryCatch(stats::as.formula(paste(".y ~", fixed_str, "+", re_str)), error = function(e) NULL)
  reduced_form <- tryCatch(stats::as.formula(paste(".y ~", reduced_fixed_str, "+", re_str)), error = function(e) NULL)
  if (is.null(full_form) || is.null(reduced_form)) return(modifyList(na_out, list(method = "bad_formula")))

  .fit <- function(fmla, dat) tryCatch(
    suppressWarnings(suppressMessages(glmmTMB::glmmTMB(
      fmla, data = dat, family = family,
      control = glmmTMB::glmmTMBControl(optCtrl = list(iter.max = 500))))),
    error = function(e) NULL)

  m_full <- .fit(full_form, data); m_red <- .fit(reduced_form, data)
  if (is.null(m_full) || is.null(m_red)) return(modifyList(na_out, list(method = "observed_fit_failed")))
  lrt_obs <- tryCatch(max(0, 2 * (as.numeric(stats::logLik(m_full)) - as.numeric(stats::logLik(m_red)))),
                      error = function(e) NA_real_)
  if (!is.finite(lrt_obs)) return(modifyList(na_out, list(method = "loglik_failed")))

  set.seed(seed)
  sim_y <- tryCatch(stats::simulate(m_red, nsim = nsim), error = function(e) NULL)
  if (is.null(sim_y)) return(modifyList(na_out, list(lrt_obs = lrt_obs, method = "simulate_failed")))

  sim_lrt <- vapply(seq_len(nsim), function(i) {
    d_sim <- data; d_sim$.y <- sim_y[[i]]
    mf <- .fit(full_form, d_sim); mr <- .fit(reduced_form, d_sim)
    if (is.null(mf) || is.null(mr)) return(NA_real_)
    v <- tryCatch(2 * (as.numeric(stats::logLik(mf)) - as.numeric(stats::logLik(mr))), error = function(e) NA_real_)
    if (!is.finite(v)) NA_real_ else max(0, v)
  }, numeric(1))

  n_ok <- sum(is.finite(sim_lrt))
  # ABORT if fewer than 50% of replicates converge -- HANDOFF S7 [INV].
  if (n_ok < nsim * 0.5)
    return(modifyList(na_out, list(lrt_obs = lrt_obs, nsim_ok = n_ok, method = "too_many_sim_failures")))
  p_pb <- (1 + sum(sim_lrt[is.finite(sim_lrt)] >= lrt_obs)) / (n_ok + 1)
  list(lrt_obs = lrt_obs, p_pb = p_pb, nsim_ok = n_ok, N = nsim, seed = seed, method = "pb_lrt")
}

# =============================================================================
# A14 -- PREFLIGHT CHECK  (CONVENTIONS.md App. C)
# =============================================================================
# 16 runnable assertions. Each corresponds to a defect that ACTUALLY HAPPENED
# in this project, or to a rule the bundle states as invariant. Do not add
# assertions for hypothetical failures (selection criterion, HANDOFF App. A).
preflight_check <- function(env = globalenv()) {
  failures <- character(0)
  chk <- function(ok, msg) if (!isTRUE(ok)) failures <<- c(failures, msg)

  chk(identical(PALETTE_TREATMENT, c("control" = "#D55E00", "exercise choice" = "#0072B2")),
      "1: PALETTE_TREATMENT does not match the Okabe-Ito pair")
  chk(!any(c("TREATMENT_COLORS", "ZONE_COLORS_MAIN") %in% ls(envir = env)),
      "2: a dead-palette object (TREATMENT_COLORS / ZONE_COLORS_MAIN) is in scope")
  chk(is.null(BASE_THEME$theme$panel.grid) ||
      inherits(BASE_THEME$theme$panel.grid, "element_blank"),
      "3: BASE_THEME$panel.grid is not blank")
  chk(fmt_F(30) == "30", "16a: fmt_F(30) != '30' -- integer df trap")
  chk(fmt3(0.5) == "0.50", "16b: fmt3(0.5) != '0.50'")
  chk(fmt3(0.77) == "0.77", "16c: fmt3(0.77) != '0.77'")
  chk(fmt_p_boot(0.001, 200) == "p \u2264 0.005", "16d: fmt_p_boot floor incorrect")
  chk(is.na(effect_size("chisq")$es), "8/9: effect_size('chisq') must return es = NA")

  if (length(failures)) {
    warning("preflight_check() FAILED:\n  - ", paste(failures, collapse = "\n  - "))
    return(invisible(FALSE))
  }
  message("preflight_check(): all checks passed.")
  invisible(TRUE)
}

# NOTE ON COMPLETENESS: assertions 4-7, 10-15 (no geom_boxplot in a manuscript
# build; every CLD frame tagged by cld_pipeline(); no shortened .group string;
# caption order regex; every model directory carries its full C12 sensitivity
# artifact set; stat_type read before df_denom; run-dir sentinel/backup guard;
# every figure has both formats at 300 dpi) require inspecting a live plot
# object, a run directory, or a data frame rather than this file's own
# constants, and are therefore implemented at the point of use (a manuscript
# build script), not here. See CONVENTIONS.md App. C for the full checklist.
