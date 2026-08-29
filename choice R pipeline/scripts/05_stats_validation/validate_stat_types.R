# =============================================================================
# validate_stat_types.R
# -----------------------------------------------------------------------------
# Audits the statistical test type actually used for every fitted indicator in
# both pipelines, producing TWO outputs:
#
#   1. CSV: <out_dir>/stat_type_validation.csv
#      One row per indicator with: label, family, model class, ANOVA mode
#      (F-KR / F-Satterthwaite / F-OLS / Wald-chisq), # terms tested, # of
#      terms with non-NA df_denom (the disambiguator), and per-term test
#      details for the focal interaction term.
#
#   2. Console summary block: pretty per-pipeline tally + flag list.
#
# INPUT
#   - Pipeline A: STEP5_OUT/<label>/anova.csv  + aicc_selection.csv
#                 (pre-existing artefacts of every run_lmm/glmm/betaglmm fit)
#   - Pipeline B: harvested LIVE from globalenv objects (fit_cort, mono_cell_models)
#                 — only if the b3 pipeline has been sourced this session.
#
# USAGE (single fresh R session AFTER pipeline runs):
#   source(file.path(PROJECT_ROOT, "choice R pipeline/scripts/05_stats_validation/validate_stat_types.R"))
#
# Or stand-alone scan of a specific STEP5_OUT directory:
#   STEP5_OUT_TO_AUDIT <- file.path(PROJECT_ROOT, "choice R pipeline/STEP5_stats/STEP5_stats_20260511_175522")
#   source(".../validate_stat_types.R")
# =============================================================================

suppressPackageStartupMessages({
  library(tibble); library(dplyr); library(readr); library(purrr)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# ---- 0. Locate Pipeline A STEP5 directory ----------------------------------
.find_step5_dir <- function() {
  if (exists("STEP5_OUT_TO_AUDIT", envir = globalenv(), inherits = FALSE) &&
      dir.exists(get("STEP5_OUT_TO_AUDIT", envir = globalenv())))
    return(get("STEP5_OUT_TO_AUDIT", envir = globalenv()))
  if (exists("STEP5_OUT", envir = globalenv(), inherits = FALSE) &&
      dir.exists(get("STEP5_OUT", envir = globalenv()))) {
    d <- get("STEP5_OUT", envir = globalenv())
    if (length(list.files(d, pattern = "anova\\.csv$", recursive = TRUE)) > 0)
      return(d)
  }
  # Fall back: most-recent STEP5_stats_* dir that actually has anova.csv files
  root <- file.path(PROJECT_ROOT, "choice R pipeline/STEP5_stats")
  if (!dir.exists(root)) return(NA_character_)
  cands <- list.files(root, pattern = "^STEP5_stats_", full.names = TRUE)
  cands <- cands[file.info(cands)$isdir %in% TRUE]
  if (length(cands) == 0) return(NA_character_)
  cands <- cands[order(file.info(cands)$mtime, decreasing = TRUE)]
  for (cd in cands) {
    if (length(list.files(cd, pattern = "anova\\.csv$", recursive = TRUE)) > 0)
      return(cd)
  }
  NA_character_
}

STEP5_DIR <- .find_step5_dir()
OUT_DIR   <- if (!is.na(STEP5_DIR)) STEP5_DIR else file.path(PROJECT_ROOT, "choice R pipeline")
message("[validate] Pipeline A STEP5 root: ", STEP5_DIR %||% "(none found)")

# ---- 1. Helper: classify a single anova.csv -------------------------------
# Returns a 1-row tibble describing the test type, inferred from columns.
# Also detects skipped.csv (response near-constant) and non-indicator dirs.
.classify_anova_csv <- function(label, anova_path, aicc_path) {
  ind_dir <- dirname(anova_path)
  skipped_path <- file.path(ind_dir, "skipped.csv")
  if (file.exists(skipped_path) && !file.exists(anova_path)) {
    note <- tryCatch(
      suppressMessages(readr::read_csv(skipped_path, show_col_types = FALSE))$note[1],
      error = function(e) "skipped (unknown reason)"
    )
    return(tibble(label = label, pipeline = "A", found = TRUE, parsed = FALSE,
                  anova_mode = "SKIPPED",
                  focal_stat_label = sprintf("[skipped] %s", note)))
  }
  if (!file.exists(anova_path)) {
    # Treat as non-indicator dir (no fit, no skip note) → drop silently downstream
    return(tibble(label = label, pipeline = "A", found = FALSE))
  }

  av <- tryCatch(suppressMessages(readr::read_csv(anova_path, show_col_types = FALSE)),
                 error = function(e) NULL)
  if (is.null(av) || nrow(av) == 0)
    return(tibble(label = label, pipeline = "A", found = TRUE, parsed = FALSE))

  # AICc → infer model class / family
  fit_method <- NA_character_
  re_winner  <- NA_character_
  if (!is.na(aicc_path) && file.exists(aicc_path)) {
    ai <- tryCatch(suppressMessages(readr::read_csv(aicc_path, show_col_types = FALSE)),
                   error = function(e) NULL)
    if (!is.null(ai) && nrow(ai) > 0) {
      win <- ai[ai$selected %in% TRUE, , drop = FALSE]
      if (nrow(win) > 0) {
        fit_method <- as.character(win$fit_method[1]) %||% NA_character_
        re_winner  <- as.character(win$re[1]) %||% NA_character_
      }
    }
  }

  # 2026-08-18: `df_denom` is NOT a reliable F-test disambiguator. The engine
  # deliberately repurposes it to carry nsim on parametric-bootstrap rows
  # (stat_type == "LRT-PB"; see the LRT-PB note in anova_caption_str). Every
  # beta-GLMM therefore has a df_denom of 1000 on its treatment row, which made
  # this auditor report 9 models as "F-KR" reading F(1, 1000) -- a bootstrap
  # replicate count rendered as denominator degrees of freedom. The engine's own
  # per-row `stat_type` column is authoritative; use it when present and keep the
  # heuristic only as a fallback for older runs that predate that column.
  stat_types <- if ("stat_type" %in% names(av)) unique(stats::na.omit(as.character(av$stat_type))) else character(0)
  pb_rows    <- if ("stat_type" %in% names(av)) which(as.character(av$stat_type) == "LRT-PB") else integer(0)
  # An F-test needs a df_denom on a row that is NOT a bootstrap row.
  .df_rows <- if ("df_denom" %in% names(av)) which(!is.na(av$df_denom) & av$df_denom > 0) else integer(0)
  has_df_denom <- length(setdiff(.df_rows, pb_rows)) > 0
  # Second disambiguator: presence of a "Residuals" row indicates lm/aov output
  # (F-OLS) where df_denom = residuals-row df. Used by zone_*_desc paths that
  # write F statistic into a column they mislabelled "chisq".
  has_residuals_row <- !is.null(av$term) && "Residuals" %in% av$term

  # Heuristic for family → inferred test family
  fam_glmm_dir <- file.exists(file.path(dirname(anova_path), "glmm_family.csv"))
  if (fam_glmm_dir) {
    family_used <- tryCatch(
      suppressMessages(readr::read_csv(
        file.path(dirname(anova_path), "glmm_family.csv"),
        show_col_types = FALSE
      ))$family[1],
      error = function(e) NA_character_)
  } else {
    family_used <- NA_character_
  }

  # ANOVA mode label:
  # MLE fit_method (with no df_denom) ⇒ glmmTMB beta-GLMM
  is_beta_glmm <- identical(fit_method, "MLE") ||
                  grepl("beta", label, ignore.case = TRUE) ||
                  file.exists(file.path(dirname(anova_path), "beta_glmm.note"))
  anova_mode <- if (length(pb_rows) > 0 && !has_df_denom) {
    # Focal inference came from the project's parametric-bootstrap LRT.
    "LRT-PB (parametric bootstrap)"
  } else if (has_df_denom) {
    # F-test: KR or OLS — distinguish by whether RE was fitted at all
    if (identical(fit_method, "OLS") || (is.na(re_winner) || identical(re_winner, "1")))
      "F-OLS"
    else
      "F-KR"
  } else if (has_residuals_row) {
    # lm-derived ANOVA written without df_denom (zone_*_desc paths). F-OLS.
    "F-OLS (label mismatch — CSV says chisq, actually F)"
  } else {
    # No df_denom column populated ⇒ Wald chi-sq (LMM fallback, GLMM, beta-GLMM)
    if (!is.na(family_used) && family_used %in% c("poisson", "negbin"))
      "Wald-chisq (GLMM)"
    else if (is_beta_glmm)
      "Wald-chisq (beta-GLMM)"
    else
      "Wald-chisq"
  }
  # Inject residuals df as df_denom for the focal term, when applicable
  if (has_residuals_row && !has_df_denom) {
    res_df <- suppressWarnings(as.numeric(av$df[av$term == "Residuals"][1]))
    if (is.finite(res_df) && "df_denom" %in% names(av)) av$df_denom <- res_df
    if (!("df_denom" %in% names(av))) av$df_denom <- res_df
    has_df_denom <- TRUE
  }

  # Identify focal interaction (highest-order term)
  terms_real <- av$term[av$term != "(Intercept)" & !is.na(av$term)]
  n_colons   <- nchar(terms_real) - nchar(gsub(":", "", terms_real))
  focal_term <- if (length(terms_real) > 0) terms_real[which.max(n_colons)] else NA_character_

  # Pull stats for focal term
  focal_row <- if (!is.na(focal_term)) av[av$term == focal_term, , drop = FALSE] else NULL
  focal_stat <- NA_real_; focal_df1 <- NA_real_; focal_df2 <- NA_real_; focal_p <- NA_real_
  if (!is.null(focal_row) && nrow(focal_row) > 0) {
    focal_stat <- suppressWarnings(as.numeric(focal_row$chisq[1])) %||% NA_real_
    focal_df1  <- suppressWarnings(as.numeric(focal_row$df[1]))    %||% NA_real_
    focal_df2  <- if ("df_denom" %in% names(focal_row))
                    suppressWarnings(as.numeric(focal_row$df_denom[1])) %||% NA_real_
                  else NA_real_
    focal_p    <- suppressWarnings(as.numeric(focal_row$p_value[1])) %||% NA_real_
  }

  tibble(
    label             = label,
    pipeline          = "A",
    found             = TRUE,
    parsed            = TRUE,
    aicc_fit_method   = fit_method,
    re_winner         = re_winner,
    family_used       = family_used %||% NA_character_,
    anova_mode        = anova_mode,
    has_df_denom      = has_df_denom,
    n_terms           = length(terms_real),
    focal_term        = focal_term,
    focal_statistic   = round(focal_stat, 4),
    focal_df1         = focal_df1,
    focal_df2         = round(focal_df2, 2),
    focal_p           = signif(focal_p, 4),
    focal_stat_label  = if (!is.na(focal_term)) {
      if (!is.na(focal_df2)) sprintf("F(%g, %.2f) = %.3f, p = %.4g",
                                      focal_df1, focal_df2, focal_stat, focal_p)
      else                   sprintf("χ²(%g) = %.3f, p = %.4g",
                                      focal_df1, focal_stat, focal_p)
    } else NA_character_
  )
}

# ---- 2. Walk Pipeline A STEP5 directory ------------------------------------
a_rows <- list()
if (!is.na(STEP5_DIR) && dir.exists(STEP5_DIR)) {
  ind_dirs <- list.dirs(STEP5_DIR, recursive = FALSE)
  for (d in ind_dirs) {
    label  <- basename(d)
    av_p   <- file.path(d, "anova.csv")
    ai_p   <- file.path(d, "aicc_selection.csv")
    a_rows[[length(a_rows) + 1]] <- .classify_anova_csv(label, av_p, ai_p)
  }
}
a_tbl <- if (length(a_rows) > 0) bind_rows(a_rows) else tibble()
# Drop non-indicator dirs (no anova.csv AND no skipped.csv)
if (nrow(a_tbl) > 0 && "found" %in% names(a_tbl))
  a_tbl <- a_tbl[a_tbl$found %in% TRUE, , drop = FALSE]

# ---- 3. Pipeline B: harvest live objects -----------------------------------
.classify_live_fit <- function(label, fit, family_hint = NA_character_) {
  if (is.null(fit)) return(tibble(label = label, pipeline = "B", found = FALSE))

  cls <- class(fit)[1]
  fam <- tryCatch(stats::family(fit)$family, error = function(e) family_hint)

  # Anova mode: for cortisol / monoamine cell LMMs we use KR
  av <- tryCatch({
    if (inherits(fit, "lmerMod") || inherits(fit, "lmerModLmerTest")) {
      car::Anova(fit, type = 3, test = "F", ddf = "Kenward-Roger")
    } else if (inherits(fit, "lm")) {
      car::Anova(fit, type = 3)
    } else if (inherits(fit, "glmmTMB")) {
      car::Anova(fit, type = "III")
    } else NULL
  }, error = function(e) NULL)

  if (is.null(av))
    return(tibble(label = label, pipeline = "B", found = TRUE, parsed = FALSE,
                  fit_class = cls, family_used = fam))

  d <- as.data.frame(av); d$term <- rownames(d)
  has_F   <- any(grepl("^F$|^F value", names(d)))
  has_chi <- any(grepl("^Chisq|^LR", names(d)))
  has_df2 <- any(grepl("^DenDF|^Df.den|^den\\.df", names(d)))

  pv_col   <- grep("^Pr|^p.value", names(d), ignore.case = TRUE, value = TRUE)[1]
  stat_col <- grep("^F$|^F value|^Chisq|^LR stat", names(d), ignore.case = TRUE, value = TRUE)[1]
  df_col   <- grep("^Df$|^NumDF|^num.df", names(d), ignore.case = TRUE, value = TRUE)[1]
  den_col  <- grep("^DenDF|^den.df|^df.den", names(d), ignore.case = TRUE, value = TRUE)[1]

  terms_real <- d$term[d$term != "(Intercept)"]
  n_colons   <- nchar(terms_real) - nchar(gsub(":", "", terms_real))
  focal_term <- if (length(terms_real) > 0) terms_real[which.max(n_colons)] else NA_character_
  fr <- if (!is.na(focal_term)) d[d$term == focal_term, , drop = FALSE] else NULL

  focal_stat <- NA_real_; focal_df1 <- NA_real_; focal_df2 <- NA_real_; focal_p <- NA_real_
  if (!is.null(fr) && nrow(fr) > 0) {
    if (!is.na(stat_col)) focal_stat <- suppressWarnings(as.numeric(fr[[stat_col]][1]))
    if (!is.na(df_col))   focal_df1  <- suppressWarnings(as.numeric(fr[[df_col]][1]))
    if (!is.na(den_col))  focal_df2  <- suppressWarnings(as.numeric(fr[[den_col]][1]))
    if (!is.na(pv_col))   focal_p    <- suppressWarnings(as.numeric(fr[[pv_col]][1]))
  }

  anova_mode <- if (has_F && has_df2) "F-KR" else
                if (has_F)            "F-OLS" else
                "Wald-chisq"

  tibble(
    label             = label,
    pipeline          = "B",
    found             = TRUE,
    parsed            = TRUE,
    aicc_fit_method   = NA_character_,
    re_winner         = NA_character_,
    family_used       = fam %||% NA_character_,
    anova_mode        = anova_mode,
    has_df_denom      = has_df2,
    n_terms           = length(terms_real),
    focal_term        = focal_term,
    focal_statistic   = round(focal_stat, 4),
    focal_df1         = focal_df1,
    focal_df2         = round(focal_df2, 2),
    focal_p           = signif(focal_p, 4),
    focal_stat_label  = if (!is.na(focal_term)) {
      if (!is.na(focal_df2)) sprintf("F(%g, %.2f) = %.3f, p = %.4g",
                                      focal_df1, focal_df2, focal_stat, focal_p)
      else                   sprintf("χ²(%g) = %.3f, p = %.4g",
                                      focal_df1, focal_stat, focal_p)
    } else NA_character_
  )
}

b_rows <- list()
if (exists("fit_cort", envir = globalenv(), inherits = FALSE)) {
  b_rows[[length(b_rows) + 1]] <-
    .classify_live_fit("cortisol__primary", get("fit_cort", envir = globalenv()))
}
if (exists("mono_cell_models", envir = globalenv(), inherits = FALSE)) {
  mcm <- get("mono_cell_models", envir = globalenv())
  if (is.list(mcm)) {
    for (cid in names(mcm)) {
      cell <- mcm[[cid]]
      if (is.list(cell) && !is.null(cell$fit)) {
        b_rows[[length(b_rows) + 1]] <-
          .classify_live_fit(paste0("monoamine_cell__", cid), cell$fit)
      }
    }
  }
}
b_tbl <- if (length(b_rows) > 0) bind_rows(b_rows) else tibble()

# ---- 4. Combine, save, and print summary ----------------------------------
all_tbl <- bind_rows(a_tbl, b_tbl)
out_csv <- file.path(OUT_DIR, "stat_type_validation.csv")
readr::write_csv(all_tbl, out_csv)
message("[validate] Wrote ", nrow(all_tbl), " rows → ", out_csv)

# Pretty per-pipeline summary
.summary_block <- function(tbl, pipeline_name) {
  if (nrow(tbl) == 0) {
    cat(sprintf("\n--- Pipeline %s: no rows ---\n", pipeline_name)); return(invisible())
  }
  cat(sprintf("\n=== Pipeline %s: %d indicators ===\n", pipeline_name, nrow(tbl)))
  tab <- table(tbl$anova_mode, useNA = "ifany")
  for (i in seq_along(tab)) {
    nm <- names(tab)[i]; if (is.na(nm) || !nzchar(nm)) nm <- "(NA / unparsed)"
    cat(sprintf("  %-25s : %d\n", nm, as.integer(tab[i])))
  }
  # Flag indicators where focal df2 is NA but mode says F (sanity check)
  flagged <- tbl[(grepl("^F", tbl$anova_mode) & is.na(tbl$focal_df2)) |
                 (tbl$anova_mode == "Wald-chisq" & !is.na(tbl$focal_df2)), , drop = FALSE]
  if (nrow(flagged) > 0) {
    cat(sprintf("\n  *** %d indicator(s) with mode/df mismatch: ***\n", nrow(flagged)))
    for (i in seq_len(nrow(flagged))) {
      cat(sprintf("    - %s : mode=%s, focal_df2=%s\n",
                  flagged$label[i], flagged$anova_mode[i],
                  as.character(flagged$focal_df2[i])))
    }
  }
}
.summary_block(a_tbl, "A (behaviour)")
.summary_block(b_tbl, "B (endocrine)")

cat("\n--- Per-indicator focal-term stat ---\n")
print(all_tbl[, c("pipeline", "label", "anova_mode", "focal_term", "focal_stat_label")],
      n = nrow(all_tbl))

invisible(all_tbl)
