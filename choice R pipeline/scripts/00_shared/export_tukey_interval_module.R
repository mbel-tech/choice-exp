# =============================================================================
# export_tukey_interval_module.R
# =============================================================================
# Exports the treatment x interval post-hoc contrasts for the FOURTEEN
# behavioural outcomes the manuscript reports, to Markdown and to Word.
#
# The set is Table 1's, grouped as the three figures:
#
#   Figure A  zone preference    ALR(flow), ALR(high), ALR(medium), ALR(low)
#   Figure B  engagement         crossings, flow bouts per minute, longest flow
#                                and calm bout, mean flow and calm bout duration
#   Figure C  collective         NND, IID, school area, school speed
#
# For each outcome the document carries:
#   1. the omnibus Type III ANOVA rows (treatment, interval, treatment x interval)
#   2. all pairwise contrasts across the six treatment x interval cells
#   3. simple effects: treatment within each interval, and interval within each
#      treatment
#   4. the compact letter display, under the canonical letter policy
#
# PROVENANCE. Outcomes fitted by STEP5 are read from that engine's own contrast
# CSVs -- the same files the figures and the manuscript cite, so the document
# cannot disagree with them. The SEQ and BOUT engines do not write Tukey
# contrasts (their post-hoc block is deliberately gated on a significant
# omnibus), so for those three outcomes the Level-B model recorded in the engine
# output is refitted here and the contrasts computed. The refit uses the exact
# model string the engine selected, so it is the same model, not a new one.
# Every table states its provenance.
#
# ADJUSTMENT. Tukey is requested throughout. emmeans substitutes Sidak on the
# six-cell family, because Tukey is exact only for a single set of pairwise
# comparisons; the substitution emmeans actually applied is recorded per table
# rather than assumed.
#
# CLD LETTERS come from cld_letter_policy.R, the same source the figures use, so
# a letter here is by construction the letter drawn in the panel.
#
# ONE SOURCE OF TRUTH. The .docx is produced from the .md via Pandoc rather than
# written independently, so the two documents cannot drift apart.
#
# USAGE:
#   "C:/Program Files/R/R-4.6.0/bin/Rscript.exe" --vanilla <this file>
#   (auto-detects the latest STEP5 / SEQ / BOUT runs; override by pre-setting
#    STEP5_OUT, SEQ_OUT, BOUT_OUT)
# =============================================================================

suppressPackageStartupMessages({
  library(emmeans)
  library(lme4)
})
has_lmerTest <- requireNamespace("lmerTest", quietly = TRUE)
if (has_lmerTest) suppressPackageStartupMessages(library(lmerTest))
has_multcomp <- requireNamespace("multcomp", quietly = TRUE)

PIPE_ROOT  <- file.path(PROJECT_ROOT, "choice R pipeline")
SHARED     <- file.path(PIPE_ROOT, "scripts", "00_shared")
OUT_DIR    <- file.path(PROJECT_ROOT, "interval_contrasts_export")
if (!exists("PANDOC")) source(file.path(PROJECT_ROOT, "config.R"))

source(file.path(SHARED, "cld_letter_policy.R"))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n")

# ---- locate the latest run of each engine ----------------------------------
# `sentinel` is a file written late in the run. Without it this picks up a run
# that is still being written and silently exports half-present results.
.latest_run <- function(step_name, sentinel = NULL) {
  parent <- file.path(PIPE_ROOT, "output", step_name)
  if (!dir.exists(parent)) return(NULL)
  d <- list.dirs(parent, full.names = TRUE, recursive = FALSE)
  d <- d[grepl(paste0("^", step_name, "_\\d{8}_\\d{6}$"), basename(d))]
  if (!is.null(sentinel)) {
    ok <- vapply(d, function(x) file.exists(file.path(x, sentinel)), logical(1))
    if (any(!ok)) message("skipping ", sum(!ok), " incomplete ", step_name, " run(s)")
    d <- d[ok]
  }
  if (!length(d)) return(NULL)
  d[which.max(file.mtime(d))]
}
if (!exists("STEP5_OUT") || !dir.exists(STEP5_OUT))
  STEP5_OUT <- .latest_run("STEP5_stats", "centroid_speed_timepoint/anova.csv")
if (!exists("SEQ_OUT")   || !dir.exists(SEQ_OUT))
  SEQ_OUT   <- .latest_run("SEQ_output",  "seq_anova_trial_x_timepoint.csv")
if (!exists("BOUT_OUT")  || !dir.exists(BOUT_OUT))
  BOUT_OUT  <- .latest_run("BOUT_output", "bout_anova_trial_x_timepoint.csv")
msg("STEP5: ", STEP5_OUT)
msg("SEQ  : ", SEQ_OUT)
msg("BOUT : ", BOUT_OUT)

# ---- the fourteen reported outcomes -----------------------------------------
# INDICATOR SET (2026-08-18), matching Table 1 and Figures A/B/C. Figure letters are
# placeholders; manuscript numbering is assigned when the figures are placed.
OUTCOMES <- list(
  # ---- Figure A: zone preference across the velocity gradient ----------------
  list(no =  1, family = "Figure A — Zone preference", panel = "A",
       label = "ALR(flow vs calm)", engine = "STEP5", dir = "zone_main_timepoint"),
  list(no =  2, family = "Figure A — Zone preference", panel = "B",
       label = "ALR(high vs calm)", engine = "STEP5", dir = "zone_sec_high_timepoint"),
  list(no =  3, family = "Figure A — Zone preference", panel = "C",
       label = "ALR(medium vs calm)", engine = "STEP5", dir = "zone_sec_medium_timepoint"),
  list(no =  4, family = "Figure A — Zone preference", panel = "D",
       label = "ALR(low vs calm)", engine = "STEP5", dir = "zone_sec_low_timepoint"),
  # ---- Figure B: pattern of engagement with the flow -------------------------
  list(no =  5, family = "Figure B — Engagement pattern", panel = "A",
       label = "Flow↔calm crossings of the school", engine = "STEP5", dir = "flux_timepoint"),
  list(no =  6, family = "Figure B — Engagement pattern", panel = "B",
       label = "Flow bouts per minute", engine = "BOUT", metric = "bout_rate_flow"),
  list(no =  7, family = "Figure B — Engagement pattern", panel = "C",
       label = "Longest flow bout (s)", engine = "BOUT", metric = "max_flow_bout_s"),
  list(no =  8, family = "Figure B — Engagement pattern", panel = "D",
       label = "Longest calm bout (s)", engine = "BOUT", metric = "max_calm_bout_s"),
  # Figure B rebuilt to 2x2 on 2026-08-17. The dwell times were dropped from the
  # figure but stay in Table 1 and the Results text -- they are the MEAN of the
  # same distribution whose MAXIMUM panels C/D plot, so their contrasts are still
  # needed. Labels follow Table 1, which names them mean bout duration rather
  # than dwell time.
  #
  # Entropy rate and the commitment index were dropped on 2026-08-18: both left
  # the manuscript entirely, and the Table S10 that used to justify keeping them
  # here no longer exists (the supplementary was restructured, and its full
  # outcome listing removed). Contrasts for outcomes the paper does not report
  # do not belong in a document that accompanies it.
  list(no =  9, family = "Figure B — Engagement pattern", panel = "not figured",
       label = "Mean flow bout duration (s)", engine = "SEQ",
       metric = "dwell_Flow", alphabet = "binary"),
  list(no = 10, family = "Figure B — Engagement pattern", panel = "not figured",
       label = "Mean calm bout duration (s)", engine = "SEQ",
       metric = "dwell_Calm", alphabet = "binary"),
  # ---- Figure C: collective movement ----------------------------------------
  list(no = 11, family = "Figure C — Collective movement", panel = "A",
       label = "Mean nearest-neighbour distance (NND, cm)", engine = "STEP5", dir = "nnd_timepoint"),
  list(no = 12, family = "Figure C — Collective movement", panel = "B",
       label = "Mean inter-individual distance (IID, cm)", engine = "STEP5", dir = "iid_timepoint"),
  list(no = 13, family = "Figure C — Collective movement", panel = "C",
       label = "Mean school area (cm²)", engine = "STEP5", dir = "hull_area_timepoint"),
  list(no = 14, family = "Figure C — Collective movement", panel = "D",
       label = "Mean school speed (cm s⁻¹)", engine = "STEP5", dir = "centroid_speed_timepoint")
)
stopifnot(!anyDuplicated(vapply(OUTCOMES, function(o) o$no, numeric(1))))

# ---- number formatting ------------------------------------------------------
f_num <- function(x, d = 3) {
  ifelse(is.na(x) | !is.finite(x), "—", formatC(as.numeric(x), format = "f", digits = d))
}
f_p <- function(p) {
  p <- suppressWarnings(as.numeric(p))
  ifelse(is.na(p) | !is.finite(p), "—",
         ifelse(p < 0.001, "< 0.001", formatC(p, format = "f", digits = 3)))
}
stars <- function(p) {
  p <- suppressWarnings(as.numeric(p))
  ifelse(is.na(p) | !is.finite(p), "",
  ifelse(p < 0.001, "\\*\\*\\*", ifelse(p < 0.01, "\\*\\*", ifelse(p < 0.05, "\\*", ""))))
}
# Escape the pipe so a value can never break a Markdown table row.
esc <- function(x) gsub("|", "\\|", as.character(x), fixed = TRUE)

md_table <- function(df, align = NULL) {
  if (is.null(df) || !nrow(df)) return("_No rows available._\n")
  hdr <- names(df)
  if (is.null(align)) align <- c("l", rep("r", length(hdr) - 1))
  sep <- vapply(align, function(a)
    switch(a, l = ":---", r = "---:", c = ":---:", "---"), character(1))
  body <- apply(df, 1, function(r) paste0("| ", paste(esc(r), collapse = " | "), " |"))
  paste0(paste0("| ", paste(esc(hdr), collapse = " | "), " |"), "\n",
         paste0("| ", paste(sep, collapse = " | "), " |"), "\n",
         paste(body, collapse = "\n"), "\n")
}

# ---- interval labels --------------------------------------------------------
TP_LAB <- c("1" = "5–25 min", "2" = "45–65 min", "3" = "85–105 min")
relab_tp <- function(x) {
  s <- as.character(x)
  s <- gsub("timepoint_f", "", s)
  for (k in names(TP_LAB)) s <- gsub(paste0("\\b", k, "\\b"), TP_LAB[[k]], s)
  s
}

# =============================================================================
# STEP5 outcomes: read the engine's own CSVs
# =============================================================================
read_step5 <- function(o) {
  d <- file.path(STEP5_OUT, o$dir)
  rd <- function(f) {
    p <- file.path(d, f)
    if (file.exists(p)) read.csv(p, stringsAsFactors = FALSE) else NULL
  }
  av <- rd("anova.csv")
  list(
    anova    = av,
    pairs    = rd("contrasts_treatmentxtimepoint_f.csv"),
    simple_t = rd("contrasts_treatmentxtimepoint_f__simple_treatment_within_timepoint_f.csv"),
    simple_i = rd("contrasts_treatmentxtimepoint_f__simple_timepoint_f_within_treatment.csv"),
    cld      = rd("cld_treatmentxtimepoint_f.csv"),
    source   = paste0("STEP5 engine output, `", basename(STEP5_OUT), "/", o$dir, "/`"),
    refit    = FALSE
  )
}

# =============================================================================
# SEQ / BOUT outcomes: refit the engine's own Level-B model, then contrast
# =============================================================================
read_seqbout <- function(o) {
  if (o$engine == "SEQ") {
    anova_fp <- file.path(SEQ_OUT, "seq_anova_trial_x_timepoint.csv")
    data_fp  <- file.path(SEQ_OUT, paste0("seq_metrics_per_trial_", o$alphabet, ".csv"))
    run_lbl  <- basename(SEQ_OUT)
  } else {
    anova_fp <- file.path(BOUT_OUT, "bout_anova_trial_x_timepoint.csv")
    data_fp  <- file.path(BOUT_OUT, "bout_metrics_per_session.csv")
    run_lbl  <- basename(BOUT_OUT)
  }
  if (!file.exists(anova_fp) || !file.exists(data_fp)) return(NULL)

  av_all <- read.csv(anova_fp, stringsAsFactors = FALSE)
  sel <- av_all$metric == o$metric
  if (!is.null(o$alphabet) && "alphabet" %in% names(av_all))
    sel <- sel & av_all$alphabet == o$alphabet
  av <- av_all[sel, , drop = FALSE]
  if (!nrow(av)) return(NULL)

  form <- av$model[1]
  if (is.na(form) || !nzchar(form)) return(NULL)

  raw <- read.csv(data_fp, stringsAsFactors = FALSE)
  if (!o$metric %in% names(raw)) return(NULL)
  d <- raw[is.finite(raw[[o$metric]]), , drop = FALSE]
  d$val         <- d[[o$metric]]
  d$timepoint_f <- factor(d$timepoint)
  d$treatment   <- factor(trimws(tolower(as.character(d$treatment))),
                          levels = .CLD_TREAT_ORDER)
  if (nlevels(droplevels(d$timepoint_f)) < 3) return(NULL)

  fit <- tryCatch(suppressMessages(suppressWarnings(
    if (has_lmerTest) lmerTest::lmer(stats::as.formula(form), data = d)
    else lme4::lmer(stats::as.formula(form), data = d))),
    error = function(e) NULL)
  if (is.null(fit)) return(NULL)

  em <- tryCatch(suppressMessages(emmeans::emmeans(fit, ~ treatment * timepoint_f)),
                 error = function(e) NULL)
  if (is.null(em)) return(NULL)

  pr   <- suppressMessages(pairs(em, adjust = "tukey"))
  smry <- summary(pr)
  adj_note <- {
    m <- attr(smry, "mesg")
    m <- if (is.null(m)) character(0) else grep("adjust", m, value = TRUE, ignore.case = TRUE)
    if (length(m)) paste(m, collapse = "; ") else "Tukey requested"
  }

  simple_t <- tryCatch(as.data.frame(suppressMessages(pairs(
    emmeans::emmeans(fit, ~ treatment | timepoint_f), adjust = "tukey"))),
    error = function(e) NULL)
  simple_i <- tryCatch(as.data.frame(suppressMessages(pairs(
    emmeans::emmeans(fit, ~ timepoint_f | treatment), adjust = "tukey"))),
    error = function(e) NULL)

  cld <- NULL
  if (has_multcomp) {
    cld <- tryCatch({
      x <- as.data.frame(multcomp::cld(em, Letters = letters, adjust = "tukey"))
      names(x)[names(x) == ".group"] <- ".group"
      x
    }, error = function(e) NULL)
  }

  # Engine ANOVA rows, mapped onto the STEP5 column names used downstream.
  av_std <- data.frame(
    term      = av$term,
    df        = av$df1,
    df_denom  = av$df2,
    chisq     = av$F,
    p_value   = av$p,
    stringsAsFactors = FALSE
  )

  list(anova = av_std, pairs = as.data.frame(smry),
       simple_t = simple_t, simple_i = simple_i, cld = cld,
       source = paste0("Level-B model refitted from `", run_lbl, "`: `", form, "`"),
       refit = TRUE, adj_note = adj_note, n_used = nrow(d))
}

# =============================================================================
# Rendering
# =============================================================================
TERM_LAB <- c(treatment = "Treatment", timepoint_f = "Interval",
              `treatment:timepoint_f` = "Treatment × Interval")

render_anova <- function(av) {
  if (is.null(av)) return("_ANOVA not available._\n")
  keep <- av[av$term %in% names(TERM_LAB), , drop = FALSE]
  if (!nrow(keep)) return("_ANOVA not available._\n")
  # df as integers (author decision 2026-08-17); the fractional Kenward-Roger /
  # Satterthwaite values remain in the engines' anova CSVs.
  out <- data.frame(
    Term    = unname(TERM_LAB[keep$term]),
    `F`     = f_num(keep$chisq, 3),
    df1     = f_num(keep$df, 0),
    df2     = f_num(keep$df_denom, 0),
    p       = f_p(keep$p_value),
    ` `     = stars(keep$p_value),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  md_table(out)
}

render_pairs <- function(pr, contrast_col = "contrast", extra = NULL) {
  if (is.null(pr) || !nrow(pr)) return("_Not available._\n")
  cn <- names(pr)
  est <- if ("estimate" %in% cn) pr$estimate else NA
  se  <- if ("SE" %in% cn) pr$SE else NA
  df  <- if ("df" %in% cn) pr$df else NA
  tt  <- if ("t.ratio" %in% cn) pr$t.ratio else if ("z.ratio" %in% cn) pr$z.ratio else NA
  pv  <- if ("p.value" %in% cn) pr$p.value else NA
  out <- data.frame(Contrast = relab_tp(pr[[contrast_col]]), stringsAsFactors = FALSE)
  if (!is.null(extra) && extra %in% cn)
    out[[if (extra == "timepoint_f") "Interval" else "Treatment"]] <-
      if (extra == "timepoint_f") relab_tp(pr[[extra]]) else as.character(pr[[extra]])
  out$Estimate <- f_num(est, 3)
  out$SE       <- f_num(se, 3)
  out$df       <- f_num(df, 0)   # integer, as in the ANOVA tables
  out$t        <- f_num(tt, 3)
  out$p        <- f_p(pv)
  out[[" "]]   <- stars(pv)
  md_table(out)
}

render_cld <- function(cld) {
  if (is.null(cld) || !nrow(cld)) return("_No compact letter display available._\n")
  tp_col <- if ("timepoint_f" %in% names(cld)) "timepoint_f" else
            if ("tp_num" %in% names(cld)) "tp_num" else return("_No CLD._\n")
  d <- cld
  d$tp_num <- suppressWarnings(as.integer(as.character(d[[tp_col]])))
  d$treatment <- trimws(as.character(d$treatment))
  d$.group <- trimws(as.character(d$.group))
  d <- .cld_remap_control_first(d)
  .cld_check_control_first(d, label = "export")
  d <- d[order(match(tolower(d$treatment), .CLD_TREAT_ORDER), d$tp_num), , drop = FALSE]
  out <- data.frame(
    Treatment = d$treatment,
    Interval  = unname(TP_LAB[as.character(d$tp_num)]),
    Mean      = f_num(if ("emmean" %in% names(d)) d$emmean else NA, 3),
    SE        = f_num(if ("SE" %in% names(d)) d$SE else NA, 3),
    Group     = d$.group,
    stringsAsFactors = FALSE, check.names = FALSE
  )
  md_table(out)
}

# =============================================================================
# Build the document
# =============================================================================
L <- c()
add <- function(...) L <<- c(L, paste0(...))

add("# Treatment × interval post-hoc contrasts")
add("")
add("Pairwise contrasts for the **final indicator set** (15 outcomes across ",
    "Figures A–C), at the trial × interval level (N = 48 sessions; 16 trials × ",
    "3 intervals).")
add("")
add("**The paper reports the trial × interval analysis only.** The aggregated / ",
    "trial-level module is not reported; it remains in the pipeline output as a ",
    "cross-check. Figure letters here are placeholders — manuscript numbering is ",
    "assigned when the figures are placed.")
add("")
add("Degrees of freedom are shown as integers throughout. Kenward-Roger and ",
    "Satterthwaite return fractional df; the fractional values are retained in ",
    "the engines' `anova.csv` files, which are the analytic record.")
add("")
add("Intervals are **5–25 min**, **45–65 min** and **85–105 min** after trial onset. ",
    "Interval is fitted as a three-level factor in every model, so all interval and ",
    "treatment × interval terms are 2-df tests (see `METHODS_CHANGES.md` §2).")
add("")
add("Compact letter display groups are relabelled so that **control at the first ",
    "interval always carries \"a\"**. This is a relabelling only — which cells share ",
    "a letter is unchanged — and it uses the same policy the figures use ",
    "(`cld_letter_policy.R`), so a letter here is the letter drawn in the panel.")
add("")
add("Significance marks: \\* p < 0.05, \\*\\* p < 0.01, \\*\\*\\* p < 0.001. ",
    "Contrasts are on the model scale and are not back-transformed.")
add("")
add("---")
add("")

cur_family <- ""
for (o in OUTCOMES) {
  if (!identical(o$family, cur_family)) {
    cur_family <- o$family
    add("## ", cur_family)
    add("")
  }
  msg("outcome ", o$no, ": ", o$label)
  r <- if (o$engine == "STEP5") read_step5(o) else read_seqbout(o)

  add("### ", o$no, ". ", o$label)
  add("")
  if (is.null(r)) {
    add("_Could not be assembled — engine output missing._")
    add(""); next
  }
  add("*Figure panel:* ", o$panel, "  ")
  add("*Source:* ", r$source, "  ")
  if (isTRUE(r$refit)) {
    add("*Sessions used:* ", r$n_used, " of 48  ")
    add("*p-value adjustment actually applied:* ", r$adj_note)
  } else {
    # Do not assert which adjustment emmeans finally applied: the engine
    # requested Tukey, but the CSV does not record whether emmeans substituted
    # Sidak (it does so for some family shapes, notably the compact letter
    # display, and not for others). The refitted outcomes above report the
    # applied method empirically because the model object is in hand here; for
    # these it is not, so the claim is limited to what is actually known.
    add("*p-value adjustment:* Tukey requested by the engine across the six ",
        "treatment × interval cells. emmeans substitutes Šidák for some family ",
        "shapes and the engine's CSV does not record which was finally applied; ",
        "compare the refitted outcomes above, where the applied method is stated.")
  }
  add("")
  add("**Omnibus (Type III)**"); add(""); add(render_anova(r$anova)); add("")
  add("**All pairwise contrasts across the six treatment × interval cells**")
  add(""); add(render_pairs(r$pairs)); add("")
  add("**Simple effect — treatment within each interval**")
  add(""); add(render_pairs(r$simple_t, extra = "timepoint_f")); add("")
  add("**Simple effect — interval within each treatment**")
  add(""); add(render_pairs(r$simple_i, extra = "treatment")); add("")
  add("**Compact letter display**")
  add(""); add(render_cld(r$cld)); add("")
}

add("---")
add("")
add("Generated by `scripts/00_shared/export_tukey_interval_module.R`.")

md_path <- file.path(OUT_DIR, "Interval_treatment_Tukey_contrasts.md")
writeLines(L, md_path, useBytes = TRUE)
msg("Markdown written: ", md_path)

# ---- Markdown -> Word, so the two documents cannot diverge -------------------
docx_path <- file.path(OUT_DIR, "Interval_treatment_Tukey_contrasts.docx")
if (file.exists(PANDOC)) {
  st <- system2(PANDOC, c(shQuote(md_path), "-o", shQuote(docx_path)),
                stdout = TRUE, stderr = TRUE)
  if (file.exists(docx_path)) msg("Word written: ", docx_path)
  else msg("Pandoc failed: ", paste(st, collapse = " | "))
} else {
  msg("Pandoc not found at ", PANDOC, " -- .docx not produced.")
}
msg("=== export complete ===")
