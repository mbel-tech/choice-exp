# scripts/07_cross_validate/verify_easy_scripts_from_deposit.R
# -----------------------------------------------------------------------------
# Tier 3: the claim a stranger can check.
#
# Runs each mini-script that backs a reported behavioural outcome with the data
# source forced to the Zenodo deposit, and checks the treatment F and denominator
# df against the values published in the manuscript -- the same constants
# asserted by scripts/06_workbook/verify_zenodo_deposit.R.
#
# The two paths are genuinely independent, which is what makes agreement
# informative rather than circular:
#
#   verify_zenodo_deposit.R  fits a hand-written lmerTest::lmer() with the
#                            transform hardcoded per outcome.
#   the mini-scripts         route through run_lmm_analysis, which selects the
#                            random-effect structure by AICc over six candidates
#                            and searches for a transform.
#
# Two different machines, one deposited file, the same numbers.
#
# Requires only the deposit -- no pipeline output, no internal CSVs. That is the
# point: it is the reader's path, not the author's.
#
# Usage:
#   source("paths.R")
#   source(file.path(PROJECT_ROOT,
#     "choice R pipeline/scripts/07_cross_validate/verify_easy_scripts_from_deposit.R"))
# -----------------------------------------------------------------------------

if (!exists("PROJECT_ROOT")) stop("Run source(\"paths.R\") first.", call. = FALSE)

.ES <- file.path(PROJECT_ROOT, "choice R pipeline/easy_scripts/by_timepoint")

# outcome F and DenDF as published. Source: verify_zenodo_deposit.R:19-31.
.SPEC <- list(
  list(script = "logit_flow_by_tp.R",           label = "alr_flow (main zone)",   F = 17.503896, df2 = 14),
  list(script = "lr_high_by_tp.R",              label = "alr_high",               F =  7.529668, df2 = 14),
  list(script = "lr_medium_by_tp.R",            label = "alr_medium",             F =  5.275043, df2 = 14),
  list(script = "lr_low_by_tp.R",               label = "alr_low",                F =  2.178357, df2 = 14),
  list(script = "mean_nnd_cm_by_tp.R",          label = "mean_nnd_cm",            F = 10.697117, df2 = 14),
  list(script = "mean_iid_cm_by_tp.R",          label = "mean_iid_cm",            F =  3.192528, df2 = 14),
  list(script = "mean_hull_area_cm2_by_tp.R",   label = "mean_school_area_cm2",   F =  6.480284, df2 = 14),
  list(script = "mean_centroid_spd_cm_by_tp.R", label = "mean_school_speed_cm_s", F =  0.745289, df2 = 14)
)

.TOL <- 1e-3

# Known PRE-EXISTING divergence between a mini-script and the published result.
# Not caused by reading the deposit: verified to reproduce identically from the
# pipeline-side easy_scripts_dataset.csv, with byte-identical AICc tables.
#
# lr_medium_by_tp.R lets AICc choose freely over RE_WIDE and selects (1 | tank)
# -- AICc 177.10 against 181.46 for (1 | phys_trial_id), a clear preference and
# not a near-tie. The manuscript reports (1 | phys_trial_id), because
# DECISIONS_LOG.md D3 fixes that random effect BY DESIGN and lets information
# criteria decide only optional nuisance terms. The mini-script does not
# implement that constraint, so it answers a slightly different question and
# lands on F(1,39) = 8.93 rather than the reported F(1,14) = 5.28.
#
# Recorded rather than silently tolerated, and rather than "fixed" here:
# changing an RE selection rule is an analysis decision, not a data-plumbing one.
.PRE_EXISTING <- c("alr_medium")

cat("\n=== Tier 3: mini-scripts refit from the Zenodo deposit ===\n\n")
cat(sprintf("%-28s %12s %12s %9s %9s  %s\n",
            "outcome", "F_refit", "F_published", "df2", "transform", "verdict"))

.ok <- 0L
.rows <- list()

for (sp in .SPEC) {
  env <- new.env(parent = globalenv())
  out <- tryCatch({
    withCallingHandlers(
      suppressMessages(suppressWarnings(
        sys.source(file.path(.ES, sp$script), envir = env))),
      warning = function(w) invokeRestart("muffleWarning"))
    res <- get("res", envir = env)
    a <- res$anova
    row <- a[a$term == "treatment", ]
    list(F = row$chisq[1], df2 = row$df_denom[1], tf = res$transform,
         src = get("dat", envir = env))
  }, error = function(e) list(err = conditionMessage(e)))

  if (!is.null(out$err)) {
    cat(sprintf("%-28s %12s %12.6f %9s %9s  ERROR: %s\n",
                sp$label, "-", sp$F, "-", "-", substr(out$err, 1, 60)))
    .rows[[length(.rows)+1L]] <- data.frame(outcome = sp$label, status = "ERROR",
                                            detail = out$err, stringsAsFactors = FALSE)
    next
  }

  good <- abs(out$F - sp$F) < .TOL && abs(out$df2 - sp$df2) < .TOL
  known <- sp$label %in% .PRE_EXISTING
  if (good) .ok <- .ok + 1L
  status <- if (good) "MATCH" else if (known) "PRE-EXISTING" else "*** MISMATCH ***"
  cat(sprintf("%-28s %12.6f %12.6f %9.4f %9s  %s\n",
              sp$label, out$F, sp$F, out$df2, out$tf, status))
  .rows[[length(.rows)+1L]] <- data.frame(
    outcome = sp$label, status = status,
    detail = sprintf("F=%.6f df2=%.4f transform=%s", out$F, out$df2, out$tf),
    stringsAsFactors = FALSE)
}

.expected <- length(.SPEC) - length(.PRE_EXISTING)
cat("\n---\n")
cat(sprintf("%d / %d reported outcomes reproduced from the deposit alone\n",
            .ok, length(.SPEC)))
if (length(.PRE_EXISTING))
  cat(sprintf("%d known pre-existing divergence(s), source-independent: %s\n",
              length(.PRE_EXISTING), paste(.PRE_EXISTING, collapse = ", ")))
if (.ok < .expected) {
  stop("Tier 3 failed: ", .expected - .ok, " outcome(s) beyond the known ",
       "pre-existing divergence did not reproduce from the deposit.",
       call. = FALSE)
}
cat("Tier 3 GREEN (against the ", .expected, " outcomes the mini-scripts and the\n",
    "manuscript are expected to agree on).\n\n", sep = "")
invisible(do.call(rbind, .rows))
