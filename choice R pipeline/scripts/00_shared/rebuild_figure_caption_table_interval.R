# =============================================================================
# rebuild_figure_caption_table_interval.R
# =============================================================================
# Regenerates easy_scripts/standalone/anova_statement_table_interval.csv from
# the CURRENT model outputs.
#
# WHY: the interval figures do not read their statistical annotations from the
# fitted models -- they look them up in that static CSV (see `.fig10_cap_csv()`).
# If the CSV is not rebuilt after a re-run, the figures silently display stale
# statistics from the previous fit, which is the single easiest way for a wrong
# number to reach the manuscript.
#
# HISTORY (2026-08-16). This script previously existed to keep the caption in
# step with a PER-OUTCOME choice between a continuous (1-df) and a categorical
# (2-df) parameterisation of interval, which AICc made differently for different
# indicators. That policy is withdrawn. Every interval model is now fitted with
# interval as a three-level factor (.ALLOW_CONTINUOUS_TP = FALSE in
# activity_analysis_STATS_choice_exp.R), so every treatment x interval term is a
# 2-df test and the parameterisation no longer varies between indicators. See
# METHODS_CHANGES.md section 2 for the reasoning and the AICc figures behind it.
#
# The winner cross-check below is retained: it now serves as an ASSERTION that
# the uniform rule actually held in the run being read, and will stop() if any
# model somehow reports a continuous-timepoint term.
#
# Run AFTER the behaviour engine has produced STEP5_stats/*/<model>/anova.csv
# and aicc_selection.csv, and BEFORE regenerating the manuscript figures.
# =============================================================================

# Auto-detect the most recent run unless the caller pre-sets STEP5. Hard-coding
# a timestamp here is how this file silently went stale before.
if (!exists("STEP5") || !nzchar(STEP5) || !dir.exists(STEP5)) {
  .p <- file.path(PROJECT_ROOT, "choice R pipeline/output/STEP5_stats")
  .r <- list.dirs(.p, full.names = TRUE, recursive = FALSE)
  .r <- .r[grepl("^STEP5_stats_\\d{8}_\\d{6}$", basename(.r))]
  if (!length(.r)) stop("No STEP5_stats_* run folder found under ", .p)
  STEP5 <- .r[which.max(file.mtime(.r))]
}
cat("Reading run:", STEP5, "\n")
OUT   <- file.path(PROJECT_ROOT, "choice R pipeline/easy_scripts/standalone/anova_statement_table_interval.csv")

# variable name used by the figures -> model output directory + pretty label
MAP <- list(
  mean_iid_cm            = c("iid_timepoint",            "Mean IID"),
  mean_nnd_cm            = c("nnd_timepoint",             "Mean NND"),
  mean_school_speed_cm_s = c("centroid_speed_timepoint",  "Mean school speed"),
  mean_school_area_cm2   = c("hull_area_timepoint",       "Mean school area"),
  # Polarisation removed 2026-08-16: excluded from the manuscript on
  # identity-invariance grounds (heading requires persistent per-fish identity,
  # which the tracker cannot guarantee; every other collective indicator is a
  # permutation-invariant function of the per-frame position set). No caption is
  # generated for it, so it cannot reach a figure.
  switches_per_session   = c("switches_timepoint",        "Switches per session"),
  main_cells_flow_tp     = c("zone_main_timepoint",       "Main cells flow"),
  sub_cells_high_clr_tp  = c("zone_sec_high_timepoint",   "Sub-cells high CLR")
)

# Unicode subscript digits, matching the existing file's formatting exactly.
SUB <- c("\u2080","\u2081","\u2082","\u2083","\u2084",
         "\u2085","\u2086","\u2087","\u2088","\u2089")
to_sub <- function(x) {
  s <- format(x, trim = TRUE)
  paste0(vapply(strsplit(s, "")[[1]], function(ch)
    if (ch %in% as.character(0:9)) SUB[as.integer(ch) + 1L]
    else if (ch == ".") "." else ch, character(1)), collapse = "")
}
# Unified number-display convention (author decision, 2026-08-07).
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
fmt3 <- function(x) {
  if (!is.finite(x)) return("NA")
  s <- sprintf("%.3f", x)
  if (grepl("0$", s)) s <- substr(s, 1, nchar(s) - 1)
  s
}
fmt_p <- function(p) if (p < 0.001) "p < 0.001" else sprintf("p = %s", fmt3(p))
# 2026-08-09: a test statistic or effect size below 0.001 renders as "< 0.001"
# rather than a long 4-sig-fig decimal ("0.0009258") or a rounded-away zero
# ("0.00"). Denominator df keep fmt_F (a df is never < 1).
fmt_Fstat <- function(x) {
  if (!is.finite(x)) return("NA")
  if (abs(x) < 0.001) return("< 0.001")
  fmt_F(x)
}
fmt_es3 <- function(x) {
  if (!is.finite(x)) return("NA")
  if (abs(x) < 0.001) return("< 0.001")
  fmt3(x)
}
num_fmt <- function(x) fmt_F(as.numeric(x))

rows <- list()
for (v in names(MAP)) {
  dir <- MAP[[v]][1]; lab <- MAP[[v]][2]
  f_anova <- file.path(STEP5, dir, "anova.csv")
  f_aicc  <- file.path(STEP5, dir, "aicc_selection.csv")
  if (!file.exists(f_anova)) { message("MISSING anova: ", f_anova); next }

  a <- read.csv(f_anova, stringsAsFactors = FALSE)
  # Interaction term is named "treatment:timepoint" (continuous) or
  # "treatment:timepoint_f" (categorical) depending on which fixed-effect
  # formula won AICc selection for this outcome.
  r <- a[a$term %in% c("treatment:timepoint", "treatment:timepoint_f"), , drop = FALSE]
  if (!nrow(r)) { message("no treatment:timepoint(_f) row: ", dir); next }
  r <- r[1, , drop = FALSE]

  # Under the uniform 2-df rule no interval model may report a continuous
  # timepoint term. Fail loudly rather than writing a 1-df caption.
  if (r$term[1] == "treatment:timepoint")
    stop("Continuous-timepoint term found in ", dir, " -- the uniform 2-df rule ",
         "did not hold in this run. Check .ALLOW_CONTINUOUS_TP.")

  # Cross-check against aicc_selection.csv that this term's df matches the
  # flagged winner, so a future re-run can't silently drift out of sync again.
  if (file.exists(f_aicc)) {
    sel <- read.csv(f_aicc, stringsAsFactors = FALSE)
    winner <- sel[sel$selected == TRUE, , drop = FALSE]
    if (nrow(winner)) {
      winner_is_cont <- grepl("treatment \\* timepoint$", winner$fixed_formula[1])
      row_is_cont    <- r$term[1] == "treatment:timepoint"
      if (!identical(winner_is_cont, row_is_cont)) {
        stop("Fixed-effect parameterisation mismatch for ", dir,
             ": anova.csv term = ", r$term[1],
             " but AICc winner fixed_formula = ", winner$fixed_formula[1])
      }
    }
  }

  F_val <- as.numeric(r$chisq[1]); df1 <- as.numeric(r$df[1])
  df2 <- as.numeric(r$df_denom[1]); p <- as.numeric(r$p_value[1])
  is_F <- grepl("^F", r$stat_type[1])

  if (is_F && is.finite(df2)) {
    # Denominator df printed as an INTEGER (author decision, 2026-08-17).
    # Kenward-Roger and Satterthwaite return fractional df; a caption reading
    # "F(2, 24.2)" invites the reader to wonder what a fifth of a degree of
    # freedom is. The fractional value is retained in anova.csv, which is the
    # analytic record; only the display is rounded.
    df2s <- format(round(df2))
    # 2026-08-09: partial eta-squared appended after p (author decision: every
    # Results-section figure reports the effect size in its ANOVA caption,
    # comma-separated, on one line). Same eta^2_p = F*df1/(F*df1+df2) identity
    # used by the trial-level table and by the engines themselves.
    eta <- (F_val * df1) / (F_val * df1 + df2)
    stmt <- sprintf("Treatment\u00d7Interval F%s,%s = %s, %s, \u03b7\u00b2%s = %s",
                    to_sub(df1), to_sub(df2s), fmt_Fstat(F_val), fmt_p(p),
                    "\u209a", fmt_es3(eta))
  } else {
    stmt <- sprintf("Treatment\u00d7Interval \u03c7\u00b2%s = %s, %s",
                    to_sub(df1), num_fmt(F_val), fmt_p(p))
  }
  rows[[length(rows) + 1L]] <- data.frame(variable = v, indicator = lab,
                                          statement = stmt,
                                          stringsAsFactors = FALSE)
}

tab <- do.call(rbind, rows)
if (file.exists(OUT) && !file.exists(paste0(OUT, ".bak_20260808")))
  file.copy(OUT, paste0(OUT, ".bak_20260808"))
con <- file(OUT, open = "w", encoding = "UTF-8")
utils::write.csv(tab, con, row.names = FALSE, fileEncoding = "UTF-8")
close(con)

cat("Rebuilt", OUT, "\n\n")
print(tab, right = FALSE)
