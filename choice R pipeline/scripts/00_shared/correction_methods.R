# =============================================================================
# correction_methods.R -- multiple-comparison correction, as a STANDALONE module
# =============================================================================
# Extracted 2026-08-07 from analysis_b3_REVISED.R (endocrine engine) and
# activity_analysis_STATS_choice_exp_REVISED.R (behaviour engine), at the
# author's request: neither the manuscript nor the two analysis engines apply
# any multiple-comparison correction by default. Every treatment-effect p-value
# reported in the manuscript is now the RAW, unadjusted p from its model.
#
# This module is NOT sourced by either main pipeline script. It exists so the
# correction work already done (dopaminergic/serotonergic family BH+Holm,
# 7-procedure sensitivity comparison for the endocrine analysis; Jacobs' D and
# primary-non-zone-family BH for the behavioural analysis) is preserved,
# documented, and re-runnable if a future revision reinstates correction --
# without it being silently mixed into either engine's default output.
#
# To use: source("00_shared/correction_methods.R") and call the function you
# need on a data.frame with at minimum an `analyte`/`Analysis` column, an
# `area` column (endocrine only), and a `p_raw` column.
# =============================================================================


# ---- generic multiple-comparison adjusters ---------------------------------

#' Single-step Sidak adjustment (not covered by base::p.adjust).
sidak_adjust <- function(p) 1 - (1 - p)^length(p)

#' Apply one correction method, optionally within groups (families).
#'
#' @param p        numeric vector of raw p-values
#' @param group    optional grouping vector (e.g. neurotransmitter family, or
#'                 brain region); if NULL, correction pools ALL of `p` together
#' @param method   one of "BH", "BY", "holm", "hochberg", "hommel",
#'                 "bonferroni", "sidak"
adjust_p <- function(p, group = NULL, method = "BH") {
  one <- function(x) if (method == "sidak") sidak_adjust(x) else p.adjust(x, method = method)
  if (is.null(group)) return(one(p))
  out <- rep(NA_real_, length(p))
  for (lv in unique(group)) {
    idx <- group == lv
    out[idx] <- one(p[idx])
  }
  out
}


# =============================================================================
# ENDOCRINE (monoamine + cortisol) correction scheme, as used 2026-08-04 to
# 2026-08-07 before being removed from the default analysis.
# =============================================================================

# Declared neurotransmitter-system families (dopaminergic = primary/
# confirmatory, serotonergic = exploratory), flat across all four brain
# regions -- replaces an earlier per-region scheme (4 families of 6/7 tests).
FAMILY_DOPA <- c("DA", "DOPAC", "DOPAC/DA")
FAMILY_SERO <- c("5-HT", "5-HIAA", "5-HIAA/5-HT")

assign_endocrine_family <- function(analyte) {
  dplyr::case_when(
    analyte %in% FAMILY_DOPA ~ "dopaminergic (primary)",
    analyte %in% FAMILY_SERO ~ "serotonergic (exploratory)",
    TRUE                     ~ NA_character_
  )
}

#' Add BH (primary) and Holm (secondary, family-wise) columns to a per-cell
#' monoamine treatment-effect table, within the declared system families.
#'
#' @param df data.frame with columns `analyte`, `p_raw` (one row per
#'   analyte x brain-region cell)
add_endocrine_correction <- function(df) {
  df$family <- assign_endocrine_family(df$analyte)
  if (any(is.na(df$family)))
    stop("analyte with no declared correction family: ",
         paste(unique(df$analyte[is.na(df$family)]), collapse = ", "))
  df$p_BH   <- adjust_p(df$p_raw, df$family, "BH")
  df$p_Holm <- adjust_p(df$p_raw, df$family, "holm")
  df$p_BH_per_region <- adjust_p(df$p_raw, df$area, "BH")   # reference only
  df$confirmatory <- df$family == "dopaminergic (primary)"
  df$concordant_BH_Holm <- (df$p_BH < 0.05) == (df$p_Holm < 0.05)
  df
}

#' Every correction procedure, computed separately, for a robustness table.
#' `df` must have `analyte`, `area`, `p_raw`.
endocrine_correction_sensitivity <- function(df) {
  sys <- ifelse(df$analyte %in% FAMILY_DOPA, "dopaminergic", "serotonergic")
  methods <- c(BH = "BH", BY = "BY", Holm = "holm", Hochberg = "hochberg",
               Hommel = "hommel", Bonferroni = "bonferroni", Sidak = "sidak")
  types   <- c(BH = "FDR", BY = "FDR", Holm = "FWER", Hochberg = "FWER",
               Hommel = "FWER", Bonferroni = "FWER", Sidak = "FWER")
  grab <- function(q, an) q[df$area == "DM" & df$analyte == an]
  mk <- function(structure, method, errrate, q) data.frame(
    family_structure = structure, method = method, error_rate = errrate,
    Dm_DA = round(grab(q, "DA"), 4), Dm_DOPAC = round(grab(q, "DOPAC"), 4),
    Dm_5HT = round(grab(q, "5-HT"), 4), n_sig_of_24 = sum(q < 0.05, na.rm = TRUE))
  out <- list(mk("(no correction)", "none", "uncorrected", df$p_raw))
  for (nm in names(methods))
    out[[length(out) + 1]] <- mk("Declared: 2 system families (k=12 each)", nm,
                                 unname(types[nm]), adjust_p(df$p_raw, sys, methods[[nm]]))
  for (nm in names(methods))
    out[[length(out) + 1]] <- mk("Pooled: single family of 24", nm,
                                 unname(types[nm]), adjust_p(df$p_raw, NULL, methods[[nm]]))
  out[[length(out) + 1]] <- mk("Per brain region (4 families of 6)", "BH", "FDR",
                               adjust_p(df$p_raw, df$area, "BH"))
  out[[length(out) + 1]] <- mk("Dm-only system families (k=3 each)", "BH", "FDR",
                               ifelse(df$area == "DM", adjust_p(df$p_raw, paste(df$area, sys), "BH"), NA))
  do.call(rbind, out)
}


# =============================================================================
# BEHAVIOURAL correction scheme, as computed by
# activity_analysis_STATS_choice_exp.R (pre-dates this module; documented here
# for reference and reuse, not re-implemented from scratch).
# =============================================================================
# Families (as built by that script's `.bh_row()` calls, lines ~3520-3623):
#   bh_primary        : Jacobs' D beta-GLMM, zone preference (8 tests) -- HEADLINE
#   bh_suppl_primary  : non-zone primary indicators, trial-level (4 tests)
#   bh_suppl          : non-zone treatment:timepoint interactions (4 tests)
#   bh_jacobs_sens    : Jacobs' D atanh sensitivity check (8 tests, never BH-corrected)
#   bh_exploratory    : descriptive only, never BH-corrected
#
#' Apply BH within each of the three correctable behavioural families and
#' return them updated. `bh_jacobs_sens` and `bh_exploratory` are returned
#' unchanged (they were never BH-corrected even when this module was in use).
apply_behavioural_bh <- function(bh_primary, bh_suppl_primary, bh_suppl) {
  bh_primary$p_BH        <- adjust_p(bh_primary$p_raw,        method = "BH")
  bh_suppl_primary$p_BH  <- adjust_p(bh_suppl_primary$p_raw,  method = "BH")
  bh_suppl$p_BH          <- adjust_p(bh_suppl$p_raw,          method = "BH")
  list(bh_primary = bh_primary, bh_suppl_primary = bh_suppl_primary, bh_suppl = bh_suppl)
}
