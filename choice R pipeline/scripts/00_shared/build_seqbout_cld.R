# =============================================================================
# build_seqbout_cld.R
# =============================================================================
# The sequence and bout engines report ANOVAs and (conditionally) simple-effect
# contrasts, but they never write a compact letter display. Panels sourced from
# them therefore had no CLD, while STEP5-sourced panels in the same figure did.
# This builds the missing CLDs so every panel of Figures A-C can carry letters.
#
# HOW: for each metric, refit the SAME Level-B model the engine selected -- the
# formula is recorded in the engine's own anova table, so this is a refit, not a
# new model -- then take emmeans over treatment x timepoint_f and letter it.
#
# Letters go through .cld_remap_control_first() from cld_letter_policy.R, the
# same rule the figures and the Tukey export use, so a letter written here is
# by construction the letter drawn in the panel.
#
# Output: one cld_treatmentxtimepoint_f.csv per metric under
#   output/SEQBOUT_cld/<metric>/
# shaped exactly like a STEP5 model directory's CLD file, so the figures can
# consume it through the existing `cld_path` argument with no special-casing.
#
# NOTE ON ADJUSTMENT: emmeans substitutes Sidak for Tukey on a 6-cell family
# (Tukey is exact only for one set of pairwise comparisons). The substitution it
# actually applied is recorded per metric in the manifest, rather than assumed.
#
# USAGE: "C:/Program Files/R/R-4.6.0/bin/Rscript.exe" --vanilla <this file>
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(lme4); library(emmeans)
})
has_lmerTest <- requireNamespace("lmerTest", quietly = TRUE)
if (has_lmerTest) suppressPackageStartupMessages(library(lmerTest))
stopifnot(requireNamespace("multcomp", quietly = TRUE))

PIPE   <- file.path(PROJECT_ROOT, "choice R pipeline")
SHARED <- file.path(PIPE, "scripts", "00_shared")
OUTDIR <- file.path(PIPE, "output", "SEQBOUT_cld")
source(file.path(SHARED, "cld_letter_policy.R"))
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n")

.latest <- function(step, sentinel) {
  p <- file.path(PIPE, "output", step)
  d <- list.dirs(p, full.names = TRUE, recursive = FALSE)
  d <- d[grepl(paste0("^", step, "_\\d{8}_\\d{6}$"), basename(d))]
  d <- d[vapply(d, function(x) file.exists(file.path(x, sentinel)), logical(1))]
  if (!length(d)) stop("no complete ", step, " run")
  d[which.max(file.mtime(d))]
}
SEQ  <- .latest("SEQ_output",  "seq_anova_trial_x_timepoint.csv")
BOUT <- .latest("BOUT_output", "bout_anova_trial_x_timepoint.csv")
msg("SEQ : ", basename(SEQ))
msg("BOUT: ", basename(BOUT))

# metric, engine, alphabet (SEQ only)
JOBS <- list(
  list(metric = "entropy_rate",     engine = "SEQ",  alphabet = "binary"),
  list(metric = "dwell_Flow",       engine = "SEQ",  alphabet = "binary"),
  list(metric = "dwell_Calm",       engine = "SEQ",  alphabet = "binary"),
  # letter_order = "lone_first": control/interval 1 belongs to two
  # non-difference groups for this metric, so no labelling can give it a plain
  # "a" and the panel would carry no single "a" at all. Set here rather than
  # only in the figure so the CSV on disk holds the letters that are drawn --
  # cld_letter_policy.R's whole point is that there is one definition. The
  # remap is idempotent, so a consumer that applies it again is a no-op.
  list(metric = "bout_rate_flow",   engine = "BOUT", letter_order = "lone_first"),
  list(metric = "max_flow_bout_s",  engine = "BOUT"),
  # Added 2026-08-17 for Figure B panel D. Without it that panel would be the
  # only one of the four without letters. entropy_rate, dwell_* and
  # commitment_index are retained above even though they left the figure: they
  # still feed Table S10, and refitting them costs nothing.
  list(metric = "max_calm_bout_s",  engine = "BOUT"),
  list(metric = "commitment_index", engine = "BOUT")
)

manifest <- list()
for (j in JOBS) {
  if (j$engine == "SEQ") {
    av <- fread(file.path(SEQ, "seq_anova_trial_x_timepoint.csv"))
    av <- av[metric == j$metric & alphabet == j$alphabet]
    dat <- fread(file.path(SEQ, paste0("seq_metrics_per_trial_", j$alphabet, ".csv")))
  } else {
    av <- fread(file.path(BOUT, "bout_anova_trial_x_timepoint.csv"))
    av <- av[metric == j$metric]
    dat <- fread(file.path(BOUT, "bout_metrics_per_session.csv"))
  }
  if (!nrow(av) || !j$metric %in% names(dat)) { msg("SKIP ", j$metric, " (absent)"); next }
  form <- av$model[1]
  if (is.na(form) || !nzchar(form)) { msg("SKIP ", j$metric, " (no model)"); next }

  d <- dat[is.finite(get(j$metric))]
  d[, val := get(j$metric)]
  d[, timepoint_f := factor(timepoint)]
  d[, treatment := factor(trimws(tolower(as.character(treatment))),
                          levels = .CLD_TREAT_ORDER)]
  if (uniqueN(d$timepoint_f) < 3L) { msg("SKIP ", j$metric, " (<3 intervals)"); next }

  fit <- tryCatch(suppressMessages(suppressWarnings(
    if (has_lmerTest) lmerTest::lmer(stats::as.formula(form), data = d)
    else lme4::lmer(stats::as.formula(form), data = d))), error = function(e) NULL)
  if (is.null(fit)) { msg("SKIP ", j$metric, " (fit failed)"); next }

  em <- tryCatch(suppressMessages(emmeans::emmeans(fit, ~ treatment * timepoint_f)),
                 error = function(e) NULL)
  if (is.null(em)) { msg("SKIP ", j$metric, " (emmeans failed)"); next }

  cl <- tryCatch(as.data.frame(multcomp::cld(em, Letters = letters, adjust = "tukey")),
                 error = function(e) NULL)
  if (is.null(cl)) { msg("SKIP ", j$metric, " (cld failed)"); next }

  adj <- {
    m <- attr(summary(em), "mesg")
    m <- if (is.null(m)) character(0) else grep("adjust", m, value = TRUE, ignore.case = TRUE)
    if (length(m)) paste(m, collapse = "; ") else "tukey requested"
  }

  cl$.group    <- trimws(as.character(cl$.group))
  cl$treatment <- trimws(as.character(cl$treatment))

  # HOUSE RULE 3 (cld_letter_policy.R): drop letters that carry no claim.
  # multcomp's insert-and-absorb display is correct but not always minimal --
  # bout_rate_flow came out with "abd" on exercise-choice/interval-2 where the
  # "a" duplicated links already carried by "b" and "d". Reduced BEFORE the
  # remap, so the remap closes any gap the reduction opens in the a, b, c
  # sequence. The reduction refuses to act (with a warning) if the letters and
  # the p-value matrix disagree, so a mismatched adjustment cannot corrupt it.
  # contrast(method = "pairwise") rather than pairs(): emmeans exports the S3
  # METHOD pairs.emmGrid, not a `pairs` generic, so emmeans::pairs() is not a
  # thing that can be called and it silently fell through to the no-table branch.
  # Same adjustment as the cld() call above, so the p-values the reduction reads
  # are exactly the ones the letters were built from.
  prs <- tryCatch(emmeans::contrast(em, method = "pairwise", adjust = "tukey"),
                  error = function(e) NULL)
  if (is.null(prs)) {
    msg("  no pairwise table for ", j$metric, " -- letters left unreduced")
  } else {
    # coef() gives the contrast coefficients against the emmeans grid, so each
    # contrast's pair is read off as the +1 and -1 rows. Exact, and immune to
    # how emmeans chose to render the cell labels.
    grid <- as.data.frame(em)
    gid  <- paste(trimws(as.character(grid$treatment)),
                  trimws(as.character(grid$timepoint_f)))
    cid  <- paste(trimws(as.character(cl$treatment)),
                  trimws(as.character(cl$timepoint_f)))
    cf   <- stats::coef(prs)
    ccols <- names(cf)[grepl("^c[.][0-9]+$", names(cf))]
    pv   <- as.data.frame(prs)$p.value
    if (length(ccols) != length(pv) || nrow(cf) != nrow(grid)) {
      msg("  contrast coefficients do not line up for ", j$metric,
          " -- letters left unreduced")
    } else {
      pairs_df <- do.call(rbind, lapply(seq_along(ccols), function(k) {
        v  <- cf[[ccols[k]]]
        ia <- which(v > 0)[1]; ib <- which(v < 0)[1]
        data.frame(i = match(gid[ia], cid), j = match(gid[ib], cid),
                   p.value = pv[k])
      }))
      cl <- .cld_reduce_df(cl, cell_ids = cid, pairs_df = pairs_df,
                           label = j$metric)
    }
  }

  cl$tp_num    <- suppressWarnings(as.integer(as.character(cl$timepoint_f)))
  lo <- if (is.null(j$letter_order)) "control_first" else j$letter_order
  cl <- if (identical(lo, "lone_first")) .cld_remap_lone_first(cl)
        else .cld_remap_control_first(cl)
  .cld_check_control_first(cl, label = j$metric, mode = lo)
  cl$tp_num <- NULL

  dd <- file.path(OUTDIR, j$metric)
  dir.create(dd, recursive = TRUE, showWarnings = FALSE)
  fwrite(cl, file.path(dd, "cld_treatmentxtimepoint_f.csv"))
  manifest[[length(manifest) + 1L]] <- data.table(
    metric = j$metric, engine = j$engine,
    alphabet = if (is.null(j$alphabet)) NA_character_ else j$alphabet,
    n_sessions = nrow(d), model = form, adjustment = adj,
    letter_order = lo,
    groups = paste(cl$.group, collapse = " "))
  msg("wrote ", j$metric, "  groups: ", paste(cl$.group, collapse = " "))
}

if (length(manifest)) {
  mf <- rbindlist(manifest, fill = TRUE)
  fwrite(mf, file.path(OUTDIR, "seqbout_cld_manifest.csv"))
  msg("manifest: ", file.path(OUTDIR, "seqbout_cld_manifest.csv"))
}
msg("=== SEQ/BOUT CLD build complete ===")
