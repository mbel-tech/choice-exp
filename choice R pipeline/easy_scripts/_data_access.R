# easy_scripts/_data_access.R
# -----------------------------------------------------------------------------
# The IMPORT seam: read the Zenodo deposit and hand it back in the vocabulary
# the mini-scripts already speak.
#
# scripts/06_workbook/build_datasets_workbook.R is the EXPORT seam -- it renames
# internal columns to publication names on the way out. This file is its
# inverse, on the way back in. The two are asserted against each other at load
# (see .assert_no_drift below), so they cannot silently diverge.
#
# Why a mapping layer rather than rewriting the scripts: the mini-scripts share
# their vocabulary with the pipeline engines that produce the data
# (activity_analysis_STATS_choice_exp.R and analysis_b3_REVISED.R). Renaming the
# scripts to deposit vocabulary would make them disagree with the engines they
# are meant to reproduce. Only the arrival of the data changes; every filter,
# factor level, model call, post-hoc and plot stays exactly as it was.
#
# Public API
#   load_behaviour_dataset(source, require, quiet) -> 48-row data.frame
#   load_endocrine_dataset(source, require, quiet) -> 966-row data.frame
#   choice_exp_resolve_source(dataset, source)     -> resolution, no I/O
#   choice_exp_data_provenance()                   -> record of the last load
# -----------------------------------------------------------------------------

if (!exists("PROJECT_ROOT")) {
  stop("_data_access.R needs PROJECT_ROOT.\n",
       "  Run source(\"paths.R\") from the repository root first.", call. = FALSE)
}
if (!exists("DATA_ROOT")) source(file.path(PROJECT_ROOT, "config.R"))

suppressPackageStartupMessages({
  library(readr); library(dplyr)
})

.ES_DIR <- file.path(PROJECT_ROOT, "choice R pipeline", "easy_scripts")
.WB_BUILDER <- file.path(PROJECT_ROOT, "choice R pipeline", "scripts",
                         "06_workbook", "build_datasets_workbook.R")

# -----------------------------------------------------------------------------
# 1. The mapping, deposit -> internal
# -----------------------------------------------------------------------------
# Inverse of .BEH_RENAME_MAP (build_datasets_workbook.R) -- names are deposit
# columns, values are the internal names the mini-scripts select on.
#
# crossings_per_session -> zone_flux_per_session is the one that is easy to
# miss and expensive to get wrong: it is reported outcome 5, F(1,14) = 45.00.
# switches_per_session is a DIFFERENT column that happens to share its name
# across the two sources; the two are near-duplicates but not synonyms.
BEH_DEPOSIT_TO_INTERNAL <- c(
  alr_flow               = "logit_flow",
  alr_high               = "lr_high",
  alr_medium             = "lr_medium",
  alr_low                = "lr_low",
  crossings_per_session  = "zone_flux_per_session",
  bouts_per_min          = "bouts_per_min",
  longest_flow_bout_s    = "longest_flow_bout_s",
  longest_calm_bout_s    = "longest_calm_bout_s",
  mean_flow_bout_s       = "mean_flow_bout_s",
  mean_calm_bout_s       = "mean_calm_bout_s",
  mean_school_area_cm2   = "mean_hull_area_cm2",
  mean_school_speed_cm_s = "mean_centroid_spd_cm"
)

# Inverse of .BEH_CONTEXT_MAP (build_datasets_workbook.R).
BEH_CONTEXT_DEPOSIT_TO_INTERNAL <- c(
  obs_seconds    = "obs_seconds",
  n_frames       = "n_frames",
  prop_flow      = "prop_flow",
  prop_calm      = "prop_calm",
  prop_high_ac   = "prop_high",
  prop_medium_ac = "prop_medium",
  prop_low_ac    = "prop_low",
  prop_calm_ac   = "prop_calm_sec",
  clr_high       = "clr_sec_high_tp",
  clr_medium     = "clr_sec_medium_tp",
  clr_low        = "clr_sec_low_tp",
  clr_calm       = "clr_sec_calm_tp"
)

# Jacobs' D, per indicator_column_map.R.
BEH_JACOBS_DEPOSIT_TO_INTERNAL <- c(
  jacobs_d_flow   = "D_main_flow",
  jacobs_d_high   = "D_sec_high",
  jacobs_d_medium = "D_sec_medium",
  jacobs_d_low    = "D_sec_low"
)

# Endocrine recodes. The analyte map mirrors .analytes in
# _generate_easy_scripts.R; NE appears in neither source (excluded throughout).
ENDO_ANALYTE_DEPOSIT_TO_INTERNAL <- c(
  "5-HT"        = "ht_5",
  "5-HIAA"      = "hiaa_5",
  "5-HIAA/5-HT" = "hiaa_5_ratio",
  "DA"          = "da",
  "DOPAC"       = "dopac",
  "DOPAC/DA"    = "dopac_da_ratio"
)

# THE RECODE THAT FAILS SILENTLY IF FORGOTTEN. The deposit says
# "exercise choice"; the internal endocrine frame says "treat"; and all seven
# endocrine mini-scripts hardcode factor(levels = c("control","treat")). Without
# this map every exercise-choice row becomes NA and the models fit on half the
# data WITHOUT ERRORING. Asserted after application, below.
# (Behaviour needs no recode -- TREATMENT_LEVELS_g is already
#  c("control", "exercise choice"), which is what the deposit carries.)
ENDO_TREATMENT_DEPOSIT_TO_INTERNAL <- c(
  "control"         = "control",
  "exercise choice" = "treat"
)

# Columns a caller may ask for that no distributed dataset can supply.
CHOICE_EXP_UNAVAILABLE <- data.frame(
  column  = c("mean_polarisation", "trial_id"),
  dataset = c("behaviour", "behaviour"),
  ref     = c("D6", "-"),
  reason  = c(
    paste0(
      "Polarisation is excluded from reported results on identity-invariance ",
      "grounds (DECISIONS_LOG.md, D6). It is absent from BOTH distributed ",
      "datasets -- not a column of zenodo_dataset/pref_trials.csv, and not a ",
      "column of easy_scripts/easy_scripts_dataset.csv either. See ",
      "scripts/06_workbook/build_datasets_workbook.R. STEP2b still computes ",
      "it; run the tracking pipeline if you need it. The polarisation ",
      "mini-scripts are retained deliberately -- see README, 'Decisions that ",
      "govern interpretation'."),
    paste0(
      "trial_id is a row index in upstream STEP2 concatenation order. It is ",
      "not recoverable from either distributed dataset and is not used by any ",
      "model, grouping or plot. Supplied as NA_integer_ rather than ",
      "fabricated, so nobody mistakes a plausible 1..48 sequence for meaning.")
  ),
  stringsAsFactors = FALSE
)

# -----------------------------------------------------------------------------
# 2. Drift guard
# -----------------------------------------------------------------------------
# Parse a single top-level `symbol <- c(...)` out of a file and evaluate it in
# baseenv(). Side-effect free. Same idiom as _helpers.R's allowlisted-constant
# recovery. If the export seam is edited without editing this file, loading
# fails loudly instead of quietly mapping a column to the wrong name.
.extract_assignment <- function(file, symbol) {
  if (!file.exists(file)) return(NULL)
  exprs <- tryCatch(parse(file), error = function(e) NULL)
  if (is.null(exprs)) return(NULL)
  for (e in exprs) {
    if (is.call(e) && as.character(e[[1]])[1] %in% c("<-", "=") &&
        identical(as.character(e[[2]]), symbol)) {
      return(tryCatch(eval(e[[3]], envir = baseenv()), error = function(err) NULL))
    }
  }
  NULL
}

.assert_no_drift <- function() {
  checks <- list(
    list(sym = ".BEH_RENAME_MAP",  declared = BEH_DEPOSIT_TO_INTERNAL),
    list(sym = ".BEH_CONTEXT_MAP", declared = BEH_CONTEXT_DEPOSIT_TO_INTERNAL)
  )
  for (ck in checks) {
    found <- .extract_assignment(.WB_BUILDER, ck$sym)
    if (is.null(found)) next   # builder absent (e.g. partial checkout) -- skip
    if (!identical(found[order(names(found))],
                   ck$declared[order(names(ck$declared))])) {
      stop("[choice-exp] mapping drift: ", ck$sym, " in\n  ", .WB_BUILDER,
           "\nno longer matches the inverse declared in _data_access.R.\n",
           "  The export seam and the import seam must agree. Reconcile them ",
           "before running anything.", call. = FALSE)
    }
  }
  invisible(TRUE)
}

# -----------------------------------------------------------------------------
# 3. Source resolution
# -----------------------------------------------------------------------------
.DEPOSIT_FILES <- list(
  behaviour = "pref_trials.csv",
  endocrine = c("pref_fish.csv", "pref_cortisol.csv", "pref_monoamines_long.csv")
)
.INTERNAL_FILES <- list(
  behaviour = "easy_scripts_dataset.csv",
  endocrine = "easy_scripts_endo_dataset.csv"
)

choice_exp_resolve_source <- function(dataset = c("behaviour", "endocrine"),
                                      source  = "auto") {
  dataset <- match.arg(dataset)
  source  <- match.arg(source, c("auto", "deposit", "internal", "override"))

  dep <- file.path(DATA_ROOT, .DEPOSIT_FILES[[dataset]])
  int <- file.path(.ES_DIR,   .INTERNAL_FILES[[dataset]])
  dep_ok <- all(file.exists(dep))
  int_ok <- all(file.exists(int))

  pick <- switch(
    source,
    deposit  = if (dep_ok) "deposit" else NA_character_,
    internal = if (int_ok) "internal" else NA_character_,
    auto     = if (dep_ok) "deposit" else if (int_ok) "internal" else NA_character_
  )

  list(source  = pick,
       paths   = if (identical(pick, "deposit")) dep
                 else if (identical(pick, "internal")) int else character(0),
       deposit = list(paths = dep, present = dep_ok),
       internal = list(paths = int, present = int_ok))
}

.no_source_error <- function(dataset, res, source) {
  stop(
    "[choice-exp] no ", dataset, " data found (source = \"", source, "\").\n",
    "  Zenodo deposit : ", paste(res$deposit$paths, collapse = "\n                   "),
    "\n                   -> ", if (res$deposit$present) "present" else "NOT FOUND",
    "\n  pipeline CSV   : ", paste(res$internal$paths, collapse = "\n                   "),
    "\n                   -> ", if (res$internal$present) "present" else "NOT FOUND",
    "\n\n  Download the dataset from https://doi.org/10.5281/zenodo.22162227 and\n",
    "  unpack it to  ", DATA_ROOT, "\n",
    "  or set CHOICE_EXP_DATA_ROOT to wherever you put it.", call. = FALSE)
}

# -----------------------------------------------------------------------------
# 4. Provenance register
# -----------------------------------------------------------------------------
# Every load records what it read. scripts/07_cross_validate/refit_equivalence.R
# asserts against this to prove its data interception actually took effect -- a
# cross-validation that cannot fail is worse than none.
.CHOICE_EXP_PROV <- new.env(parent = emptyenv())
.CHOICE_EXP_PROV$last <- NULL

choice_exp_data_provenance <- function() .CHOICE_EXP_PROV$last

.record <- function(dataset, source, paths, n_row, n_col) {
  rec <- list(dataset = dataset, source = source, paths = paths,
              n_row = n_row, n_col = n_col, loaded_at = Sys.time())
  .CHOICE_EXP_PROV$last <- rec
  invisible(rec)
}

.check_available <- function(require, dataset) {
  if (!length(require)) return(invisible(TRUE))
  hit <- CHOICE_EXP_UNAVAILABLE[CHOICE_EXP_UNAVAILABLE$column %in% require &
                                CHOICE_EXP_UNAVAILABLE$dataset == dataset, ]
  hit <- hit[hit$column != "trial_id", ]   # trial_id is supplied as NA, not refused
  if (nrow(hit))
    stop("[choice-exp] '", hit$column[1], "' is not available.\n\n  ",
         hit$reason[1], call. = FALSE)
  invisible(TRUE)
}

# -----------------------------------------------------------------------------
# 5. Behaviour loader
# -----------------------------------------------------------------------------
load_behaviour_dataset <- function(source  = NULL,
                                   require = character(0),
                                   quiet   = FALSE) {
  .assert_no_drift()
  .check_available(require, "behaviour")

  if (is.null(source)) source <- getOption("choice_exp.data_source", DATA_SOURCE)
  override <- getOption("choice_exp.behaviour_csv", NULL)
  if (!is.null(override)) {
    dat <- readr::read_csv(override, show_col_types = FALSE, progress = FALSE)
    dat <- as.data.frame(dat)
    if (!quiet)
      message("[choice-exp] behaviour data: OVERRIDE  ", override,
              "  (", nrow(dat), " x ", ncol(dat), ")")
    .record("behaviour", "override", override, nrow(dat), ncol(dat))
    return(dat)
  }

  res <- choice_exp_resolve_source("behaviour", source)
  if (is.na(res$source)) .no_source_error("behaviour", res, source)

  if (identical(res$source, "internal")) {
    dat <- as.data.frame(readr::read_csv(res$paths, show_col_types = FALSE,
                                         progress = FALSE))
    if (!quiet)
      message("[choice-exp] behaviour data: pipeline CSV  ", res$paths,
              "  (", nrow(dat), " x ", ncol(dat), ")")
    .record("behaviour", "internal", res$paths, nrow(dat), ncol(dat))
    return(dat)
  }

  # ---- deposit ----
  # Explicit col_types, not readr's guesser: fish_density must come back numeric
  # so that as.character(DENSITY_LEVELS_g) matches the factor levels the
  # mini-scripts build at Step 3.
  raw <- readr::read_csv(
    res$paths, show_col_types = FALSE, progress = FALSE,
    col_types = readr::cols(
      phys_trial  = readr::col_integer(),
      trial_seq   = readr::col_integer(),
      interval    = readr::col_integer(),
      trial_date  = readr::col_character(),
      treatment   = readr::col_character(),
      tank        = readr::col_character(),
      motor_side  = readr::col_character(),
      .default    = readr::col_double()
    ))
  raw <- as.data.frame(raw)

  out <- raw
  .apply <- function(df, map) {
    for (dep in names(map)) {
      int <- map[[dep]]
      if (dep %in% names(df) && !identical(dep, int))
        names(df)[names(df) == dep] <- int
    }
    df
  }
  out <- .apply(out, BEH_DEPOSIT_TO_INTERNAL)
  out <- .apply(out, BEH_CONTEXT_DEPOSIT_TO_INTERNAL)
  out <- .apply(out, BEH_JACOBS_DEPOSIT_TO_INTERNAL)

  # Interval carries three internal spellings; build_datasets_workbook.R asserts
  # they always agree, so derive all three from the one deposited column.
  out$timepoint     <- as.integer(raw$interval)
  out$timepoint_f   <- as.integer(raw$interval)
  out$timepoint_num <- as.integer(raw$interval)
  out$interval      <- NULL

  out$phys_trial_id <- as.integer(raw$phys_trial)
  out$phys_trial    <- NULL

  out$fish_density_f <- as.numeric(raw$fish_density)

  # Deposit stores ISO dates; the internal CSV stores %d.%m.%Y. Reverse the
  # reformat done by build_zenodo_deposit.R.
  out$trial_date <- format(as.Date(raw$trial_date, format = "%Y-%m-%d"),
                           "%d.%m.%Y")

  # Not recoverable, not used by any model. NA rather than a fabricated index.
  out$trial_id <- NA_integer_

  stopifnot(nrow(out) == 48L,
            !anyNA(out$timepoint), !anyNA(out$phys_trial_id),
            !anyNA(out$trial_date),
            setequal(unique(out$treatment), c("control", "exercise choice")))

  missing_req <- setdiff(require, names(out))
  if (length(missing_req))
    stop("[choice-exp] requested column(s) absent from the deposit: ",
         paste(missing_req, collapse = ", "), call. = FALSE)

  if (!quiet)
    message("[choice-exp] behaviour data: zenodo deposit  ", res$paths,
            "  (", nrow(raw), " x ", ncol(raw), " -> ", nrow(out), " x ",
            ncol(out), " internal)")
  .record("behaviour", "deposit", res$paths, nrow(out), ncol(out))
  out
}

# -----------------------------------------------------------------------------
# 6. Endocrine loader
# -----------------------------------------------------------------------------
# The deposit splits what the internal CSV keeps in one long table. Rebuild it:
# cortisol_sample_id is the spine, carried by all three deposit tables;
# pref_fish.csv supplies tank / trial_key / sex, which the cortisol RE candidate
# set needs ((1|tank) and (1|trial) are both in RE_CORT_CAND).
load_endocrine_dataset <- function(source  = NULL,
                                   require = character(0),
                                   quiet   = FALSE) {
  if (is.null(source)) source <- getOption("choice_exp.data_source", DATA_SOURCE)
  override <- getOption("choice_exp.endocrine_csv", NULL)
  if (!is.null(override)) {
    dat <- as.data.frame(readr::read_csv(override, show_col_types = FALSE,
                                         progress = FALSE))
    if (!quiet)
      message("[choice-exp] endocrine data: OVERRIDE  ", override,
              "  (", nrow(dat), " x ", ncol(dat), ")")
    .record("endocrine", "override", override, nrow(dat), ncol(dat))
    return(dat)
  }

  res <- choice_exp_resolve_source("endocrine", source)
  if (is.na(res$source)) .no_source_error("endocrine", res, source)

  if (identical(res$source, "internal")) {
    dat <- as.data.frame(readr::read_csv(res$paths, show_col_types = FALSE,
                                         progress = FALSE))
    if (!quiet)
      message("[choice-exp] endocrine data: pipeline CSV  ", res$paths,
              "  (", nrow(dat), " x ", ncol(dat), ")")
    .record("endocrine", "internal", res$paths, nrow(dat), ncol(dat))
    return(dat)
  }

  # ---- deposit ----
  p <- setNames(res$paths, basename(res$paths))
  rd <- function(f, ...) as.data.frame(
    readr::read_csv(p[[f]], show_col_types = FALSE, progress = FALSE, ...))

  fish <- rd("pref_fish.csv")[, c("cortisol_sample_id", "tank",
                                  "trial_key", "sex")]
  # trial_key is read verbatim, not re-derived as paste0(tank, "_", density) --
  # build_zenodo_deposit.R already did that and shipped the result.

  cort <- rd("pref_cortisol.csv")
  cort <- cort[!is.na(cort$cortisol_plasma), ]
  cort <- merge(cort, fish, by = "cortisol_sample_id", all.x = TRUE)
  cort_out <- data.frame(
    sample_id  = as.character(cort$cortisol_sample_id),
    tank       = cort$tank,
    trial      = cort$trial_key,
    treatment  = unname(ENDO_TREATMENT_DEPOSIT_TO_INTERNAL[cort$treatment]),
    sex        = cort$sex,
    plate      = as.character(cort$plate),
    plate_date = as.character(cort$plate_date),
    area       = NA_character_,
    analyte    = "cort",
    value      = as.numeric(cort$cortisol_plasma),
    stringsAsFactors = FALSE)

  mono <- rd("pref_monoamines_long.csv")
  mono <- merge(mono, fish, by = "cortisol_sample_id", all.x = TRUE)
  mono_out <- data.frame(
    sample_id  = as.character(mono$cortisol_sample_id),
    tank       = mono$tank,
    trial      = mono$trial_key,
    treatment  = unname(ENDO_TREATMENT_DEPOSIT_TO_INTERNAL[mono$treatment]),
    sex        = mono$sex,
    plate      = NA_character_,
    plate_date = NA_character_,
    area       = as.character(mono$region),
    analyte    = unname(ENDO_ANALYTE_DEPOSIT_TO_INTERNAL[mono$analyte]),
    value      = as.numeric(mono$concentration),
    stringsAsFactors = FALSE)

  out <- rbind(cort_out, mono_out)
  out$log_value  <- ifelse(is.finite(out$value) & out$value > 0,
                           log(out$value), NA_real_)
  out$sqrt_value <- ifelse(is.finite(out$value) & out$value >= 0,
                           sqrt(out$value), NA_real_)

  # The recode is load-bearing; assert rather than hope.
  stopifnot(setequal(unique(out$treatment), c("control", "treat")),
            !anyNA(out$treatment),
            !anyNA(out$analyte))

  missing_req <- setdiff(require, names(out))
  if (length(missing_req))
    stop("[choice-exp] requested column(s) absent from the deposit: ",
         paste(missing_req, collapse = ", "), call. = FALSE)

  if (!quiet)
    message("[choice-exp] endocrine data: zenodo deposit  ", dirname(res$paths[1]),
            "  (", nrow(out), " x ", ncol(out), " internal)")
  .record("endocrine", "deposit", res$paths, nrow(out), ncol(out))
  out
}
