# =============================================================================
# rebuild_figure_caption_table.R
# =============================================================================
# Regenerates easy_scripts/standalone/anova_statement_table_trial.csv from the
# CURRENT model outputs.
#
# WHY: Figures 8 and 9 do not read their statistical annotations from the fitted
# models — they look them up in that static CSV (see `.fig9_cap_csv()` in
# activity_analysis_STATS_choice_exp.R). The file was last written 2026-06-02
# and still held pre-M1/M2 statistics, so the panels disagreed with the
# manuscript text after the models were refit (e.g. polarisation printed as
# chi-square 19.004 when the retained model gives F(1,7) = 16.738, and school
# speed printed as p = 0.015 when it is now p = 0.065). Regenerating the table
# from the anova.csv files keeps figure and text in step, and adds the partial
# eta-squared effect size in the project's X = Y form.
#
# Run AFTER the behaviour engine has produced STEP5_stats/*/<model>/anova.csv,
# and BEFORE regenerating the manuscript figures.
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
OUT   <- file.path(PROJECT_ROOT, "choice R pipeline/easy_scripts/standalone/anova_statement_table_trial.csv")

# variable name used by the figures  ->  model output directory + pretty label
MAP <- list(
  mean_iid_cm            = c("iid_aggregated",             "Mean IID"),
  mean_nnd_cm            = c("nnd_aggregated",             "Mean NND"),
  mean_school_speed_cm_s = c("centroid_speed_aggregated",  "Mean school speed"),
  mean_school_area_cm2   = c("hull_area_aggregated",       "Mean school area"),
  mean_polarisation      = c("polarisation_aggregated",    "Polarisation"),
  switches_per_session   = c("switches_aggregated",        "Switches per session"),
  logit_flow             = c("zone_flow_logit_aggregated", "Logit flow"),
  lr_high                = c("zone_sec_high_aggregated",   "High-flow log-ratio"),
  lr_medium              = c("zone_sec_medium_aggregated", "Medium-flow log-ratio"),
  lr_low                 = c("zone_sec_low_aggregated",    "Low-flow log-ratio")
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
# Unified number-display convention (author decision, 2026-08-07) -- see
# 00_shared/number_formatting.R. fmt_F: F/chi-square, 4 sig figs, trailing
# zeros trimmed. fmt3: everything else (p, eta2p), 3dp, single trailing-zero
# trim to 2 when the third decimal digit is exactly zero.
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
  f <- file.path(STEP5, dir, "anova.csv")
  if (!file.exists(f)) { message("MISSING: ", f); next }
  a <- read.csv(f, stringsAsFactors = FALSE)
  r <- a[a$term == "treatment", , drop = FALSE]
  if (!nrow(r)) { message("no treatment row: ", dir); next }
  F_val <- as.numeric(r$chisq[1]); df1 <- as.numeric(r$df[1])
  df2 <- as.numeric(r$df_denom[1]); p <- as.numeric(r$p_value[1])
  is_F <- grepl("^F", r$stat_type[1])

  if (is_F && is.finite(df2)) {
    eta <- (F_val * df1) / (F_val * df1 + df2)
    df2s <- fmt_F(df2)
    # 2026-08-09: single line, comma-separated (author decision). The panel
    # caption font size was reduced alongside this so F, p and eta^2_p fit
    # without clipping -- see the figure code's plot.caption sizing.
    stmt <- sprintf("Treatment F%s,%s = %s, %s, \u03b7\u00b2%s = %s",
                    to_sub(df1), to_sub(df2s), fmt_Fstat(F_val), fmt_p(p),
                    "\u209a", fmt_es3(eta))
  } else {
    # Wald chi-square: no partial eta-squared (undefined for this statistic).
    stmt <- sprintf("Treatment \u03c7\u00b2%s = %s, %s",
                    to_sub(df1), num_fmt(F_val), fmt_p(p))
  }
  rows[[length(rows) + 1L]] <- data.frame(variable = v, indicator = lab,
                                          statement = stmt,
                                          stringsAsFactors = FALSE)
}

tab <- do.call(rbind, rows)
if (file.exists(OUT) && !file.exists(paste0(OUT, ".bak_20260807")))
  file.copy(OUT, paste0(OUT, ".bak_20260807"))
con <- file(OUT, open = "w", encoding = "UTF-8")
utils::write.csv(tab, con, row.names = FALSE, fileEncoding = "UTF-8")
close(con)

cat("Rebuilt", OUT, "\n\n")
print(tab, right = FALSE)
