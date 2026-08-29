# =============================================================================
# export_supplementary_tables.R
# =============================================================================
# Builds the supplementary tables the revised Methods cites. Numbering continues
# the existing S-series, which runs to S9 in the MAV manuscript:
#
#   Table S10  Full interval-level ANOVA, EVERY analysed behavioural outcome
#   Table S11  Orthogonal linear / quadratic decomposition of the interval effect
#   Table S12  Treatment x interval post-hoc contrasts, final indicator set
#
# S10 and S11 exist to make the reporting selection-free: the main text shows a
# pre-specified subset, and these show everything that was analysed, so a reader
# can verify nothing was chosen on its result.
#
# Markdown is the source; the .docx is produced from it by Pandoc, so the two
# cannot drift apart.
# =============================================================================

PIPE   <- file.path(PROJECT_ROOT, "choice R pipeline")
OUT    <- file.path(PROJECT_ROOT, "FINAL behavioural figures and contrasts")
if (!exists("PANDOC")) source(file.path(PROJECT_ROOT, "config.R"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

.latest <- function(step, sentinel) {
  p <- file.path(PIPE, "output", step)
  d <- list.dirs(p, full.names = TRUE, recursive = FALSE)
  d <- d[grepl(paste0("^", step, "_\\d{8}_\\d{6}$"), basename(d))]
  d <- d[vapply(d, function(x) file.exists(file.path(x, sentinel)), logical(1))]
  if (!length(d)) stop("no complete ", step, " run")
  d[which.max(file.mtime(d))]
}
STEP5 <- .latest("STEP5_stats",  "centroid_speed_timepoint/anova.csv")
SEQ   <- .latest("SEQ_output",   "seq_anova_trial_x_timepoint.csv")
BOUT  <- .latest("BOUT_output",  "bout_anova_trial_x_timepoint.csv")
cat("STEP5:", basename(STEP5), "\nSEQ  :", basename(SEQ),
    "\nBOUT :", basename(BOUT), "\n")

fp <- function(p) if (!is.finite(p)) "—" else if (p < 0.001) "< 0.001" else sprintf("%.3f", p)
f2 <- function(x, d = 2) if (!is.finite(x)) "—" else formatC(x, format = "f", digits = d)
di <- function(x) if (!is.finite(x)) "—" else format(round(x))
esc <- function(x) gsub("|", "\\|", as.character(x), fixed = TRUE)

md_tab <- function(df, align) {
  sep <- vapply(align, function(a) switch(a, l = ":---", r = "---:", "---"), character(1))
  c(paste0("| ", paste(esc(names(df)), collapse = " | "), " |"),
    paste0("| ", paste(sep, collapse = " | "), " |"),
    apply(df, 1, function(r) paste0("| ", paste(esc(r), collapse = " | "), " |")))
}

L <- c("# Supplementary tables — behavioural analysis", "",
       "Interval is fitted as a three-level factor in every model, so interval and",
       "treatment × interval terms are 2-df tests throughout. Denominator degrees of",
       "freedom are shown as integers; the fractional Kenward-Roger values are retained",
       "in the analysis output. All p-values are raw: no adjustment for multiple",
       "comparisons is applied.", "")

# ---------------------------------------------------------------- S10 ------
rows <- list()
for (d in list.dirs(STEP5, full.names = TRUE, recursive = FALSE)) {
  f <- file.path(d, "anova.csv"); if (!file.exists(f)) next
  a <- read.csv(f, stringsAsFactors = FALSE)
  if (!any(grepl("timepoint", a$term))) next
  g <- function(tm) a[a$term == tm, , drop = FALSE]
  tr <- g("treatment"); it <- g("treatment:timepoint_f"); iv <- g("timepoint_f")
  if (!nrow(tr)) next
  rows[[length(rows) + 1L]] <- data.frame(
    Engine = "STEP5", Outcome = basename(d),
    `Treatment` = sprintf("F(%s, %s) = %s, p = %s", di(tr$df[1]), di(tr$df_denom[1]),
                          f2(tr$chisq[1]), fp(tr$p_value[1])),
    `Interval` = if (nrow(iv)) sprintf("F(%s, %s) = %s, p = %s", di(iv$df[1]),
                   di(iv$df_denom[1]), f2(iv$chisq[1]), fp(iv$p_value[1])) else "—",
    `Treatment × Interval` = if (nrow(it)) sprintf("F(%s, %s) = %s, p = %s",
                   di(it$df[1]), di(it$df_denom[1]), f2(it$chisq[1]), fp(it$p_value[1])) else "—",
    check.names = FALSE, stringsAsFactors = FALSE)
}
add_engine <- function(path, engine, has_alpha) {
  a <- read.csv(path, stringsAsFactors = FALSE)
  key <- if (has_alpha) paste(a$alphabet, a$metric) else a$metric
  for (k in unique(key)) {
    s <- a[key == k, , drop = FALSE]
    gg <- function(tm) s[s$term == tm, , drop = FALSE]
    tr <- gg("treatment"); it <- gg("treatment:timepoint_f"); iv <- gg("timepoint_f")
    if (!nrow(tr)) next
    rows[[length(rows) + 1L]] <<- data.frame(
      Engine = engine, Outcome = k,
      `Treatment` = sprintf("F(%s, %s) = %s, p = %s", di(tr$df1[1]), di(tr$df2[1]),
                            f2(tr$F[1]), fp(tr$p[1])),
      `Interval` = if (nrow(iv)) sprintf("F(%s, %s) = %s, p = %s", di(iv$df1[1]),
                     di(iv$df2[1]), f2(iv$F[1]), fp(iv$p[1])) else "—",
      `Treatment × Interval` = if (nrow(it)) sprintf("F(%s, %s) = %s, p = %s",
                     di(it$df1[1]), di(it$df2[1]), f2(it$F[1]), fp(it$p[1])) else "—",
      check.names = FALSE, stringsAsFactors = FALSE)
  }
}
add_engine(file.path(SEQ, "seq_anova_trial_x_timepoint.csv"), "SEQ", TRUE)
add_engine(file.path(BOUT, "bout_anova_trial_x_timepoint.csv"), "BOUT", FALSE)
S10 <- do.call(rbind, rows)
S10 <- S10[order(S10$Engine, S10$Outcome), ]

L <- c(L, "## Table S10. Interval-level ANOVA for every analysed behavioural outcome", "",
       paste("All", nrow(S10), "outcomes analysed at the trial × interval level",
             "(N = 48 sessions), reported irrespective of significance. The main text",
             "presents a pre-specified subset selected on construct coverage; this",
             "table exists so that selection can be checked."), "",
       md_tab(S10, c("l", "l", "l", "l", "l")), "")

# ---------------------------------------------------------------- S11 ------
prow <- list()
for (d in list.dirs(STEP5, full.names = TRUE, recursive = FALSE)) {
  f <- file.path(d, "poly_contrasts.csv"); if (!file.exists(f)) next
  p <- read.csv(f, stringsAsFactors = FALSE)
  p <- p[p$block == "treatment_x_interval", , drop = FALSE]
  for (comp in c("linear", "quadratic")) {
    r <- p[p$component == comp, , drop = FALSE]
    if (!nrow(r)) next
    prow[[length(prow) + 1L]] <- data.frame(
      Engine = "STEP5", Outcome = basename(d), Component = comp,
      Estimate = f2(r$estimate[1], 3), SE = f2(r$SE[1], 3),
      df = di(r$df[1]), t = f2(r$t.ratio[1]), p = fp(r$p.value[1]),
      check.names = FALSE, stringsAsFactors = FALSE)
  }
}
for (cfg in list(list(file.path(SEQ, "seq_poly_trial_x_timepoint.csv"), "SEQ", TRUE),
                 list(file.path(BOUT, "bout_poly_trial_x_timepoint.csv"), "BOUT", FALSE))) {
  if (!file.exists(cfg[[1]])) next
  a <- read.csv(cfg[[1]], stringsAsFactors = FALSE)
  a <- a[a$block == "treatment_x_interval", , drop = FALSE]
  key <- if (cfg[[3]]) paste(a$alphabet, a$metric) else a$metric
  for (k in unique(key)) for (comp in c("linear", "quadratic")) {
    r <- a[key == k & a$component == comp, , drop = FALSE]
    if (!nrow(r)) next
    prow[[length(prow) + 1L]] <- data.frame(
      Engine = cfg[[2]], Outcome = k, Component = comp,
      Estimate = f2(r$estimate[1], 3), SE = f2(r$SE[1], 3),
      df = di(r$df[1]), t = f2(r$t[1]), p = fp(r$p[1]),
      check.names = FALSE, stringsAsFactors = FALSE)
  }
}
S11 <- do.call(rbind, prow)
S11 <- S11[order(S11$Engine, S11$Outcome, S11$Component), ]

L <- c(L, "## Table S11. Orthogonal decomposition of the treatment × interval effect", "",
       paste("Because the three intervals are equally spaced (midpoints 15, 55 and 95 min),",
             "the two-degree-of-freedom interaction partitions orthogonally into a linear",
             "(monotone trend) and a quadratic (mid-trial deviation) component. These",
             "contrasts were pre-specified and are reported for every outcome, significant",
             "or not. The omnibus interaction in Table S10 remains the primary test; the",
             "decomposition characterises the shape of an effect and no outcome is declared",
             "significant on a component contrast alone. Entries shown as — are not",
             "estimable, which occurs where a treatment × interval cell is empty."), "",
       md_tab(S11, c("l", "l", "l", "r", "r", "r", "r", "r")), "")

# ---------------------------------------------------------------- S12 ------
tk <- file.path(OUT, "Interval_treatment_Tukey_contrasts.md")
L <- c(L, "## Table S12. Treatment × interval post-hoc contrasts, final indicator set", "",
       paste("All pairwise contrasts across the six treatment × interval cells, both",
             "simple-effect directions, and the compact letter display, for each of the",
             "outcomes reported in the main text. Reproduced in full in the companion",
             "file `Interval_treatment_Tukey_contrasts.docx`."), "")
if (file.exists(tk)) {
  body <- readLines(tk, warn = FALSE, encoding = "UTF-8")
  start <- grep("^## ", body)[1]
  L <- c(L, body[start:length(body)])
}

md <- file.path(OUT, "Supplementary_tables_S10_S12.md")
writeLines(L, md, useBytes = TRUE)
cat("markdown:", md, "\n")
dx <- file.path(OUT, "Supplementary_tables_S10_S12.docx")
if (file.exists(PANDOC)) {
  system2(PANDOC, c(shQuote(md), "-o", shQuote(dx)), stdout = TRUE, stderr = TRUE)
  cat("word    :", dx, " exists:", file.exists(dx), "\n")
}
cat("S10 rows:", nrow(S10), " S11 rows:", nrow(S11), "\n")
