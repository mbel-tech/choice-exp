# =============================================================================
# rebuild_figure_caption_table_cells.R
# =============================================================================
# Regenerates easy_scripts/standalone/anova_statement_table_cells.csv from the
# CURRENT zone_main_cells_aggregated / zone_sec_cells_aggregated / (Figure 7bis)
# zone_main_timepoint anova.csv files.
#
# WHY: three panel captions across the manuscript are built directly from a
# live res$anova table with no CSV cache, unlike every sibling panel
# (.fig9_cap_csv() / .fig10_cap_csv()):
#   - Figure 8 panels C/D (cell-mean Treatment x Zone): .fig7_anova_cap()
#   - Figure 7bis (Treatment x Interval, ALR flow):     .mg_anova_cap()
# That made these three panels the remaining obstacle to a fast, model-fit-
# free figure-only regeneration path (figures_step5_module.R, 2026-08-09) --
# this script closes the gap using the SAME two caption formulas, applied to
# the already-saved anova.csv files instead of an in-memory model object.
# Note the two formulas render slightly differently on purpose (matching each
# function's own house style): .fig7_anova_cap -> "Treatment × Zone: ...",
# .mg_anova_cap -> "Treatment × Interval: ..." (colon, spaced ×).
#
# Run AFTER the STEP5 behaviour engine has produced
# STEP5_stats/*/zone_main_cells_aggregated/anova.csv,
# STEP5_stats/*/zone_sec_cells_aggregated/anova.csv, and
# STEP5_stats/*/zone_main_timepoint/anova.csv.
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
OUT   <- file.path(PROJECT_ROOT, "choice R pipeline/easy_scripts/standalone/anova_statement_table_cells.csv")

# Identical number-formatting convention used throughout the project (see
# 00_shared/number_formatting.R and the trial/interval caption-table siblings).
fmt_F <- function(x) {
  if (!is.finite(x)) return("NA")
  s <- format(signif(x, 4), scientific = FALSE, trim = TRUE)
  if (grepl(".", s, fixed = TRUE)) { s <- sub("0+$", "", s); s <- sub("\\.$", "", s) }
  s
}
fmt3 <- function(x) {
  if (!is.finite(x)) return("NA")
  s <- sprintf("%.3f", x)
  if (grepl("0$", s)) s <- substr(s, 1, nchar(s) - 1)
  s
}
fmt_p <- function(p) if (p < 0.001) "p < 0.001" else sprintf("p = %s", fmt3(p))
fmt_Fstat <- function(x) { if (!is.finite(x)) return("NA"); if (abs(x) < 0.001) return("< 0.001"); fmt_F(x) }
fmt_es3   <- function(x) { if (!is.finite(x)) return("NA"); if (abs(x) < 0.001) return("< 0.001"); fmt3(x) }
`%||%` <- function(a, b) if (is.null(a)) b else a

# Verbatim port of .fig7_anova_cap()'s formula, reading term/chisq/df/df_denom/
# p_value/stat_type from a saved anova.csv instead of a live res$anova table.
cap_from_anova_csv <- function(path) {
  if (!file.exists(path)) { message("MISSING: ", path); return(NA_character_) }
  av <- read.csv(path, stringsAsFactors = FALSE)
  row <- av[grepl("treatment.*zone$|^zone.*treatment$|treatment:zone|zone:treatment",
                  av$term, ignore.case = TRUE, perl = TRUE), ]
  if (nrow(row) == 0L) return(NA_character_)
  row <- row[1L, ]
  p <- suppressWarnings(as.numeric(row$p_value))
  if (!is.finite(p)) return(NA_character_)
  p_str <- fmt_p(p)
  chi   <- fmt_Fstat(as.numeric(row$chisq))
  st    <- row$stat_type[1] %||% "Wald-chisq"
  is_f  <- isTRUE(grepl("^F", st))
  if (is_f) {
    df1 <- as.integer(row$df)
    df2 <- suppressWarnings(as.numeric(if (!is.null(row$df_denom)) row$df_denom else NA_real_))
    if (is.finite(df2)) {
      eta <- fmt_es3((as.numeric(row$chisq) * df1) / (as.numeric(row$chisq) * df1 + df2))
      sprintf("Treatment \u00d7 Zone: F<sub>%d, %s</sub> = %s, %s, &eta;<sup>2</sup><sub>p</sub> = %s",
              df1, fmt_F(df2), chi, p_str, eta)
    } else {
      sprintf("Treatment \u00d7 Zone: F<sub>%d</sub> = %s, %s", df1, chi, p_str)
    }
  } else {
    sprintf("Treatment \u00d7 Zone: \u03c7\u00b2<sub>%d</sub> = %s, %s", as.integer(row$df), chi, p_str)
  }
}

# Verbatim port of .mg_anova_cap()'s formula (Figure 7bis: Treatment x Interval,
# ALR flow). Term pattern matches treatment:timepoint / treatment:timepoint_f,
# same convention as rebuild_figure_caption_table_interval.R.
cap_interval_from_anova_csv <- function(path) {
  if (!file.exists(path)) { message("MISSING: ", path); return(NA_character_) }
  av <- read.csv(path, stringsAsFactors = FALSE)
  row <- av[grepl("treatment.*timepoint|timepoint.*treatment", av$term,
                  ignore.case = TRUE, perl = TRUE), ]
  if (nrow(row) == 0L) return(NA_character_)
  row <- row[1L, ]
  p <- suppressWarnings(as.numeric(row$p_value))
  if (!is.finite(p)) return(NA_character_)
  p_str <- fmt_p(p)
  df1 <- if ("df" %in% names(row)) as.integer(round(row$df)) else NA_integer_
  df2 <- suppressWarnings(as.numeric(if ("df_denom" %in% names(row)) row$df_denom else NA_real_))
  chi <- fmt_F(as.numeric(row$chisq))
  if (is.finite(df2)) {
    eta <- fmt3((as.numeric(row$chisq) * df1) / (as.numeric(row$chisq) * df1 + df2))
    sprintf("Treatment × Interval: F<sub>%d, %s</sub> = %s, %s, &eta;<sup>2</sup><sub>p</sub> = %s",
            df1, fmt_F(df2), chi, p_str, eta)
  } else {
    sprintf("Treatment × Interval: F<sub>%d</sub> = %s, %s", df1, chi, p_str)
  }
}

tab <- data.frame(
  variable  = c("zone_main", "zone_sec", "zone_main_timepoint"),
  indicator = c("Zone occupancy (main)", "Area-normalized CLR (sub-zone)",
               "ALR flow x Interval (Figure 7bis)"),
  statement = c(
    cap_from_anova_csv(file.path(STEP5, "zone_main_cells_aggregated", "anova.csv")),
    cap_from_anova_csv(file.path(STEP5, "zone_sec_cells_aggregated",  "anova.csv")),
    cap_interval_from_anova_csv(file.path(STEP5, "zone_main_timepoint", "anova.csv"))
  ),
  stringsAsFactors = FALSE
)

con <- file(OUT, open = "w", encoding = "UTF-8")
utils::write.csv(tab, con, row.names = FALSE, fileEncoding = "UTF-8")
close(con)
cat("Rebuilt", OUT, "\n\n")
print(tab, right = FALSE)
