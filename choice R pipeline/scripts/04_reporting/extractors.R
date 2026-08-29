# =============================================================================
# extractors.R
# -----------------------------------------------------------------------------
# Class-dispatch helpers that turn a fitted model + per-indicator manifest row
# into a single uniform tibble row for the Ruth_statement.docx master table.
#
# Supported classes:  lmerModLmerTest, lmerMod, glmmTMB, lm
# Stub class:         brmsfit  (warns if hit; none currently in either pipeline)
#
# Public entry point:
#   extract_row(obj_expr_text, meta_expr_text, manifest_row) -> 1-row tibble
#
# Helpers used downstream by build_ruth_statement.R:
#   prettify_terms(x)            character → manuscript-friendly factor names
#   pretty_fixed(formula_obj)    formula → "Treatment × Interval"
#   pretty_random(formula_obj)   formula → "Tank, Trial" or "None — OLS fallback"
#   selection_basis(tier, meta, fit_class, transform_text) → composed string
# =============================================================================

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
})

# ---- Label map --------------------------------------------------------------
.RUTH_TERM_MAP <- c(
  treatment       = "Treatment",
  timepoint_f     = "Interval",
  timepoint       = "Interval (numeric)",
  fish_density_f  = "Density",
  fish_density    = "Density (numeric)",
  condition       = "Condition",
  area            = "Brain region",
  zone            = "Zone",
  zone_lr         = "Zone (alr pair)",
  sex             = "Sex",
  tank            = "Tank",
  trial           = "Trial",
  trial_id        = "Trial",
  phys_trial_id   = "Physical trial",
  sample_id       = "Sample",
  plate           = "Plate",
  plate_date      = "Plate date",
  date            = "Date"
)

prettify_terms <- function(x) {
  if (length(x) == 0 || is.null(x)) return(character(0))
  out <- as.character(x)
  for (nm in names(.RUTH_TERM_MAP)) {
    out <- gsub(paste0("\\b", nm, "\\b"), .RUTH_TERM_MAP[[nm]], out, perl = TRUE)
  }
  # Tidy up operator glyphs and whitespace
  out <- gsub(":", " × ", out, fixed = TRUE)        # : → ×
  out <- gsub("\\*", " × (full) ", out, perl = TRUE) # * → × (full factorial)
  out <- gsub("\\s+", " ", out)
  out <- gsub("^\\s+|\\s+$", "", out)
  out
}

# ---- Formula parsing --------------------------------------------------------
pretty_fixed <- function(formula_obj) {
  if (is.null(formula_obj)) return(NA_character_)
  rhs_no_bars <- tryCatch(lme4::nobars(formula_obj), error = function(e) formula_obj)
  if (length(rhs_no_bars) < 3) return(NA_character_)
  fixed_chr <- paste(deparse(rhs_no_bars[[3]]), collapse = " ")
  fixed_chr <- gsub("\\s+", " ", fixed_chr)
  prettify_terms(fixed_chr)
}

pretty_random <- function(formula_obj) {
  if (is.null(formula_obj)) return("None — OLS fallback")
  bars <- tryCatch(lme4::findbars(formula_obj), error = function(e) NULL)
  if (is.null(bars) || length(bars) == 0) return("None — OLS fallback")
  parts <- vapply(bars, function(b) {
    s <- paste(deparse(b), collapse = " ")
    s <- gsub("\\s+", " ", s)
    s
  }, character(1))
  # Extract grouping variable on RHS of "|"
  pretty_groups <- vapply(parts, function(p) {
    if (grepl("\\|", p)) {
      grp <- sub(".*\\|\\s*", "", p)
      grp <- gsub("[()\\s]", "", grp, perl = TRUE)
      prettify_terms(grp)
    } else prettify_terms(p)
  }, character(1))
  pretty_groups <- unique(pretty_groups)
  paste(pretty_groups, collapse = ", ")
}

# ---- Response & transformation reconciliation -------------------------------
# Pull a transformation hint from the LHS of the model formula, if the LHS is
# a wrapped call like log(x), sqrt(x), atanh(x).
.lhs_transform_hint <- function(formula_obj) {
  if (is.null(formula_obj)) return(NA_character_)
  lhs <- formula_obj[[2]]
  if (is.call(lhs)) {
    fn <- as.character(lhs[[1]])
    if (fn %in% c("log", "sqrt", "atanh", "logit", "qlogis"))
      return(fn)
  }
  NA_character_
}

.response_name <- function(formula_obj) {
  if (is.null(formula_obj)) return(NA_character_)
  lhs <- formula_obj[[2]]
  if (is.call(lhs)) {
    # log(x) → "x" (the underlying response name)
    arg <- lhs[[2]]
    return(deparse(arg))
  }
  deparse(lhs)
}

# Build the transformation text shown in the "Response scale / transformation"
# column. Precedence: meta$transform (helper's own slot) > LHS wrapper > "raw".
.transform_text <- function(meta, formula_obj, b3_trans_label = NULL,
                            b3_trans_lambda = NULL) {
  # Pipeline B per-cell monoamine: trans_label + lambda win
  if (!is.null(b3_trans_label)) {
    if (!is.null(b3_trans_lambda) && is.finite(b3_trans_lambda)) {
      return(sprintf("%s (λ = %.2f)", b3_trans_label, b3_trans_lambda))
    }
    return(as.character(b3_trans_label))
  }
  # Pipeline A helpers: $transform
  if (!is.null(meta) && !is.null(meta$transform)) {
    t <- as.character(meta$transform)
    return(switch(t,
      "identity"   = "raw",
      "log"        = "log",
      "sqrt"       = "sqrt",
      "beta_logit" = "beta squeeze + logit link",
      t))
  }
  # Fallback: parse LHS
  hint <- .lhs_transform_hint(formula_obj)
  if (!is.na(hint)) return(hint)
  "raw"
}

# ---- Selection basis composer ----------------------------------------------
selection_basis <- function(tier, meta, fit_class, transform_text,
                            b3_re_label = NULL, b3_trans_rationale = NULL) {
  tier_tag <- sprintf("[%s]", tools::toTitleCase(tier))
  parts <- character(0)

  # AICc note
  if (!is.null(meta) && !is.null(meta$aicc_table)) {
    n_cand <- nrow(meta$aicc_table)
    sel    <- meta$aicc_table[meta$aicc_table$selected %in% TRUE |
                              meta$aicc_table$delta_AICc %in% 0, , drop = FALSE]
    winner <- if (nrow(sel) > 0) sel$re[1] else NA_character_
    if (!is.na(winner)) {
      parts <- c(parts, sprintf("AICc-selected from %d RE candidate(s): %s",
                                n_cand, winner))
    } else {
      parts <- c(parts, sprintf("AICc-selected from %d RE candidate(s)", n_cand))
    }
  } else if (!is.null(b3_re_label)) {
    parts <- c(parts, sprintf("AICc-selected RE: %s", b3_re_label))
  } else {
    parts <- c(parts, "Pre-specified RE")
  }

  # Transform note
  if (!is.null(b3_trans_rationale) && nzchar(b3_trans_rationale)) {
    parts <- c(parts, sprintf("transform: %s (%s)",
                              transform_text, b3_trans_rationale))
  } else {
    parts <- c(parts, sprintf("transform: %s", transform_text))
  }

  # Fallback note
  if (!is.null(meta)) {
    if (isTRUE(meta$simplified_to_lm))
      parts <- c(parts, "OLS fallback (all RE singular)")
    if (isTRUE(meta$re_fallback))
      parts <- c(parts, "intercept-only RE fallback")
    if (!is.null(meta$stat_type))
      parts <- c(parts, sprintf("inference: %s", meta$stat_type))
  } else if (identical(fit_class, "lm")) {
    parts <- c(parts, "OLS (no RE)")
  }

  paste(tier_tag, paste(parts, collapse = "; "))
}

# ---- Class dispatchers ------------------------------------------------------
.extract_lmm <- function(fit, meta, manifest_row) {
  fm <- stats::formula(fit)
  tibble(
    response_var    = .response_name(fm),
    fixed_effects   = pretty_fixed(fm),
    random_effects  = pretty_random(fm),
    family          = "Gaussian",
    link            = "identity",
    transform       = .transform_text(meta, fm),
    fit_class       = class(fit)[1]
  )
}

.extract_glmmTMB <- function(fit, meta, manifest_row) {
  fm   <- stats::formula(fit)
  fam  <- tryCatch(stats::family(fit), error = function(e) NULL)
  fam_name <- if (!is.null(fam)) fam$family else "glmmTMB (unknown family)"
  link_name <- if (!is.null(fam)) fam$link else NA_character_
  tibble(
    response_var    = .response_name(fm),
    fixed_effects   = pretty_fixed(fm),
    random_effects  = pretty_random(fm),
    family          = tools::toTitleCase(fam_name),
    link            = link_name,
    transform       = .transform_text(meta, fm),
    fit_class       = class(fit)[1]
  )
}

.extract_lm <- function(fit, meta, manifest_row) {
  fm <- stats::formula(fit)
  tibble(
    response_var    = .response_name(fm),
    fixed_effects   = pretty_fixed(fm),
    random_effects  = "None — OLS fallback",
    family          = "Gaussian",
    link            = "identity",
    transform       = .transform_text(meta, fm),
    fit_class       = "lm"
  )
}

.extract_brms <- function(fit, meta, manifest_row) {
  warning("brmsfit detected for indicator '", manifest_row$indicator,
          "' but extractor is a stub.")
  fm <- tryCatch(stats::formula(fit), error = function(e) NULL)
  fam_name <- tryCatch(fit$family$family, error = function(e) "brms (unknown)")
  tibble(
    response_var    = if (!is.null(fm)) .response_name(fm) else NA_character_,
    fixed_effects   = if (!is.null(fm)) pretty_fixed(fm) else NA_character_,
    random_effects  = if (!is.null(fm)) pretty_random(fm) else NA_character_,
    family          = fam_name,
    link            = NA_character_,
    transform       = "brms (see source)",
    fit_class       = class(fit)[1]
  )
}

# ---- Public entry point -----------------------------------------------------
# obj_expr_text : a character string holding the R expression that yields the
#                 fitted model (e.g. "res_pol_tp$model", "fit_cort",
#                 'mono_cell_models[["ht_5__DM"]]$fit')
# meta_expr_text: optional — character string for the parent res_* list (or NA)
# manifest_row  : the manifest tibble row (as a 1-row tibble)
extract_row <- function(obj_expr_text, meta_expr_text, manifest_row,
                        envir = globalenv()) {
  # Evaluate the fit
  fit <- tryCatch(eval(parse(text = obj_expr_text), envir = envir),
                  error = function(e) NULL)
  if (is.null(fit)) {
    return(tibble(
      response_var    = NA_character_,
      fixed_effects   = NA_character_,
      random_effects  = NA_character_,
      family          = NA_character_,
      link            = NA_character_,
      transform       = NA_character_,
      fit_class       = "NULL",
      selection_basis = sprintf("[%s] Object %s is NULL or does not exist",
                                 tools::toTitleCase(manifest_row$tier),
                                 obj_expr_text)
    ))
  }
  # Evaluate the meta object (parent res_* list), if provided
  meta <- NULL
  if (!is.null(meta_expr_text) && !is.na(meta_expr_text) &&
      nzchar(meta_expr_text)) {
    meta <- tryCatch(eval(parse(text = meta_expr_text), envir = envir),
                     error = function(e) NULL)
  }
  # Pipeline B per-cell extras
  b3_re_label        <- NULL
  b3_trans_label     <- NULL
  b3_trans_lambda    <- NULL
  b3_trans_rationale <- NULL
  if (!is.null(meta) && is.list(meta) &&
      all(c("trans_label", "re_label") %in% names(meta))) {
    b3_re_label        <- meta$re_label
    b3_trans_label     <- meta$trans_label
    b3_trans_lambda    <- meta$trans_lambda
    b3_trans_rationale <- meta$trans_rationale
  }

  cls <- class(fit)[1]
  base_row <- switch(cls,
    lmerModLmerTest = .extract_lmm(fit, meta, manifest_row),
    lmerMod         = .extract_lmm(fit, meta, manifest_row),
    glmmTMB         = .extract_glmmTMB(fit, meta, manifest_row),
    lm              = .extract_lm(fit, meta, manifest_row),
    brmsfit         = .extract_brms(fit, meta, manifest_row),
    {
      warning("Unsupported class '", cls, "' for indicator '",
              manifest_row$indicator, "'")
      tibble(
        response_var    = NA_character_,
        fixed_effects   = NA_character_,
        random_effects  = NA_character_,
        family          = NA_character_,
        link            = NA_character_,
        transform       = NA_character_,
        fit_class       = cls
      )
    }
  )

  # B3 per-cell transform override
  if (!is.null(b3_trans_label)) {
    fm <- tryCatch(stats::formula(fit), error = function(e) NULL)
    base_row$transform <- .transform_text(meta, fm,
                                          b3_trans_label  = b3_trans_label,
                                          b3_trans_lambda = b3_trans_lambda)
  }

  base_row$selection_basis <- selection_basis(
    tier            = manifest_row$tier,
    meta            = meta,
    fit_class       = base_row$fit_class,
    transform_text  = base_row$transform,
    b3_re_label     = b3_re_label,
    b3_trans_rationale = b3_trans_rationale
  )
  base_row
}

# ---- Convenience: harvest the whole manifest --------------------------------
harvest_manifest <- function(manifest, envir = globalenv()) {
  out <- purrr::map_dfr(seq_len(nrow(manifest)), function(i) {
    row <- manifest[i, , drop = FALSE]
    extracted <- extract_row(row$object_expr,
                             row$meta_expr,
                             row,
                             envir = envir)
    src <- sprintf("%s @ %s:%s", row$object_expr, row$source_path, row$source_line)
    tibble(
      block            = row$block,
      indicator        = row$indicator,
      response_pretty  = row$response_pretty,
      response_var     = extracted$response_var,
      unit             = row$unit,
      fixed_effects    = extracted$fixed_effects,
      random_effects   = extracted$random_effects,
      family           = extracted$family,
      link             = extracted$link,
      transform        = extracted$transform,
      selection_basis  = extracted$selection_basis,
      source           = src
    )
  })
  out
}
