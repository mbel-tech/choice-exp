# easy_scripts/_helpers.R
# -----------------------------------------------------------------------------
# Bridge to the main pipeline (activity_analysis_STATS_choice_exp.R).
#
# Strategy:
#   1. Set EASY_SCRIPTS_HELPER_LOAD=1 so the STATS-file guard near §5 aborts
#      execution before any model fitting / figure assembly happens.
#   2. sys.source() the STATS file into globalenv() within a tryCatch that
#      catches the custom abort condition. This gives us all utility and
#      runner functions defined BEFORE the guard (run_lmm_analysis,
#      run_betaglmm_analysis, .build_jacobs, RE_WIDE/_AGG, etc.) plus the
#      data-prep tables (df, df_main_wide, etc.).
#   3. Parse the STATS file with parse() and evaluate ONLY the safe
#      assignments AFTER the guard line — function definitions and the
#      allowlisted theme / palette constants. This brings in the .mg_*
#      and .make_* graphing helpers (defined in §6) without re-running any
#      model fitting or figure rendering code.
#
# After sourcing this file, mini-scripts have direct access to:
#   - run_lmm_analysis(), run_betaglmm_analysis(), .build_jacobs()
#   - .mg_make_alr_scatter(), .mg_make_jacobs_panel(), .mg_make_cells_panel(),
#     .mg_mk_collective_panel(), .mg_make_line_with_tukey(),
#     .make_std_plot(), .make_indicator_line_plot()
#   - .mg_save(), .mg_frame(), .mg_treatment_xaxis(), .mg_wrap_caption()
#   - Theme + palette constants: BASE_THEME, TREATMENT_COLORS_PAL, color_map,
#     TIMEPOINT_LABELS, etc.
#
# Re-sourcing within the same RStudio session is fast (data is already in
# globalenv); first source per session takes ~5-10 s.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(stringr)
  library(lme4); library(lmerTest); library(glmmTMB); library(MuMIn)
  library(emmeans); library(multcomp); library(ggplot2)
})

STATS_FILE <- file.path(PROJECT_ROOT, "choice R pipeline/scripts/01_pipeline_analysis/activity_analysis_STATS_choice_exp.R")

if (!file.exists(STATS_FILE)) {
  stop("STATS file not found: ", STATS_FILE)
}

# ---- 1. Set guard signals + working directory -------------------------------
# Set the guard via THREE channels (env var, option, global var) — sys.source
# semantics around env vars are inconsistent across R versions, so we set all
# three to maximise the chance the guard at §5 of STATS picks it up.
Sys.setenv(EASY_SCRIPTS_HELPER_LOAD = "1")
options(easy_scripts_helper_load = TRUE)
assign(".EASY_HELPER_LOAD", TRUE, envir = globalenv())

.easy_orig_wd <- getwd()
.cleanup <- function() {
  Sys.unsetenv("EASY_SCRIPTS_HELPER_LOAD")
  options(easy_scripts_helper_load = NULL)
  if (exists(".EASY_HELPER_LOAD", envir = globalenv(), inherits = FALSE))
    rm(".EASY_HELPER_LOAD", envir = globalenv())
  setwd(.easy_orig_wd)
}

# Setwd to pipeline output root so .find_latest_csv() finds STEP2_output,
# STEP2b_output (the STATS file searches getwd() and dirname(getwd())).
PIPELINE_ROOT      <- file.path(PROJECT_ROOT, "choice R pipeline")
PIPELINE_OUTPUT    <- file.path(PIPELINE_ROOT, "output")
if (dir.exists(PIPELINE_OUTPUT)) setwd(PIPELINE_OUTPUT)

# Also pre-set STEP2_OUTPUT_DIR / STEP2b_OUTPUT_DIR globals to point at the
# latest dated output (the STATS file's data-loader checks these first).
.pick_latest <- function(parent_dir, prefix) {
  if (!dir.exists(parent_dir)) return(NULL)
  subs <- list.dirs(parent_dir, full.names = TRUE, recursive = FALSE)
  subs <- subs[grepl(paste0("^", prefix, "_\\d{8}_\\d{6}$"), basename(subs))]
  if (length(subs) == 0) return(NULL)
  subs[which.max(file.mtime(subs))]
}
.s2 <- .pick_latest(file.path(PIPELINE_OUTPUT, "STEP2_output"),  "STEP2_output")
.s2b <- .pick_latest(file.path(PIPELINE_OUTPUT, "STEP2b_output"), "STEP2b_output")
if (!is.null(.s2))  STEP2_OUTPUT_DIR  <- .s2
if (!is.null(.s2b)) STEP2b_OUTPUT_DIR <- .s2b

# ---- 2. sys.source STATS file (will abort at §5 guard) ----------------------
# Sourcing into globalenv lets data-prep tables and utility functions land
# directly where mini-scripts expect them. chdir=FALSE so getwd() stays at
# PIPELINE_ROOT for the duration of the source (.find_latest_csv depends on it).
tryCatch(
  sys.source(STATS_FILE, envir = globalenv(), chdir = FALSE),
  easy_scripts_helpers_only = function(c) {
    message("[easy_scripts] guard fired at §5 (model fitting skipped)")
    invisible(NULL)
  },
  error = function(e) {
    # If our classed condition somehow lands here, accept it as success
    if (inherits(e, "easy_scripts_helpers_only")) {
      message("[easy_scripts] guard fired at §5 (caught via error handler)")
      return(invisible(NULL))
    }
    # The STATS file loads STEP2 output in §1, roughly 2200 lines BEFORE the
    # easy-scripts guard in §5. On a clean clone there is no pipeline output --
    # choice R pipeline/output/ is gitignored -- so §1 stops and step 3 below,
    # which is what actually recovers run_lmm_analysis and the plot builders,
    # was never reached. That made every by_timepoint mini-script unrunnable
    # from a fresh checkout regardless of which dataset it read.
    #
    # Recognise that one failure precisely and fall through to step 3. Any
    # other error still re-raises: a genuine syntax error in the STATS file
    # must stay loud.
    if (grepl("^trial_(activity_summary|occupancy_long) not found",
              conditionMessage(e))) {
      message("[easy_scripts] no pipeline output present (", conditionMessage(e), ")")
      message("[easy_scripts] continuing: helper functions are recovered by ",
              "parsing, and the mini-scripts read their data directly.")
      return(invisible(NULL))
    }
    .cleanup()
    stop("[easy_scripts] sys.source failed: ", conditionMessage(e))
  }
)

# ---- 3. Parse STATS file for safe post-guard definitions --------------------
.easy_exprs <- parse(STATS_FILE, keep.source = TRUE)

# Allowlist of constant/theme symbols that mini-scripts may need.
# These are NOT function definitions; they're simple value assignments.
.EASY_CONST_ALLOWLIST <- c(
  # Theme + palette constants
  "BASE_THEME", "TREATMENT_COLORS", "TREATMENT_COLORS_PAL",
  "TREATMENT_SHAPES", "TIMEPOINT_LABELS", "LINE_TP_LABELS",
  "color_map", "zone_colors", "zone_shapes",
  ".MG_TAG_THEME", ".MG_GHOST_THEME", ".FIG10_CAP_OVR",
  ".MG_BODY_H", ".MG_LEG_H", ".MG_SPC_H",
  # Plot dimension/styling constants (3695-3705 of STATS)
  "PT_SIZE", "PT_ALPHA", "LW_MEAN", "LW_ERR", "ERR_W", "MEAN_W",
  "JITTER_W", "DODGE_W", "CLD_SIZE", "CLD_NUDGE_SEM_MULT",
  "CLD_NUDGE_RANGE_FRAC",
  # Position helpers used by .make_std_plot etc.
  "POS_JD", "POS_JITTER",
  # Spec lists used by the generator
  ".ind_tp_specs", ".gd_specs",
  # Model-fitting globals that run_lmm_analysis and check_and_transform read
  # but do not receive as arguments. Without these they fail with "object not
  # found" whenever §1 did not run -- which is every clean clone. All are
  # side-effect-free literals; the data frames and directory constants that §1
  # also defines are deliberately NOT here, because they are pipeline results,
  # not constants, and nothing on the mini-script path reads them.
  ".ALLOW_CONTINUOUS_TP", "N_PB_DEFAULT", "PB_SEED",
  ".UNBOUNDED_RESPONSE_PATTERNS", "ALL_TRIAL_DATES_BEH"
)

.is_safe_assignment <- function(e) {
  if (!is.call(e) || length(e) != 3L) return(FALSE)
  op <- tryCatch(as.character(e[[1]])[1], error = function(...) "")
  if (!op %in% c("<-", "=", "<<-")) return(FALSE)
  lhs <- e[[2]]; rhs <- e[[3]]

  # Pure function definition: `<- function(...) {...}` — always safe.
  if (is.call(rhs)) {
    rhs_head <- tryCatch(as.character(rhs[[1]])[1], error = function(...) "")
    if (rhs_head == "function") return(TRUE)
  }

  # Constant allowlist (theme + palette + spec lists)
  if (is.symbol(lhs)) {
    lhs_name <- as.character(lhs)
    if (lhs_name %in% .EASY_CONST_ALLOWLIST) return(TRUE)
  }
  FALSE
}

.n_fn   <- 0L
.n_const <- 0L
.n_skip <- 0L
for (e in .easy_exprs) {
  if (!.is_safe_assignment(e)) next
  # Skip if symbol already exists in globalenv (preserves the version loaded
  # by sys.source above when there's a redefinition order issue).
  lhs <- e[[2]]
  if (is.symbol(lhs) && exists(as.character(lhs), envir = globalenv(), inherits = FALSE)) {
    # Function definitions take precedence (overwrite); constants don't.
    is_fn <- is.call(e[[3]]) && tryCatch(as.character(e[[3]][[1]])[1], error = function(...) "") == "function"
    if (!is_fn) { .n_skip <- .n_skip + 1L; next }
  }
  res <- tryCatch({
    eval(e, envir = globalenv()); TRUE
  }, error = function(err) FALSE)
  if (isTRUE(res)) {
    if (is.call(e[[3]]) && tryCatch(as.character(e[[3]][[1]])[1], error = function(...) "") == "function") {
      .n_fn <- .n_fn + 1L
    } else {
      .n_const <- .n_const + 1L
    }
  } else {
    .n_skip <- .n_skip + 1L
  }
}

message(sprintf(
  "[easy_scripts/_helpers.R] loaded: %d functions, %d constants (skipped: %d)",
  .n_fn, .n_const, .n_skip
))

# STEP5_OUT is assigned inside an `if` block in the STATS file, so
# .is_safe_assignment correctly refuses it and step 3 cannot recover it. The
# model helpers write per-model CSVs there. When no pipeline output exists,
# point it at a temp directory rather than letting the first write fail.
#
# This also removes a standing side effect: sourcing the STATS file otherwise
# does dir.create(file.path(getwd(), "STEP5_stats")) unconditionally, littering
# whatever directory the user happened to be in.
if (!exists("STEP5_OUT", envir = globalenv(), inherits = FALSE)) {
  .step5_tmp <- file.path(tempdir(),
                          paste0("STEP5_stats_easy_",
                                 format(Sys.time(), "%Y%m%d_%H%M%S")))
  dir.create(.step5_tmp, recursive = TRUE, showWarnings = FALSE)
  assign("STEP5_OUT", .step5_tmp, envir = globalenv())
  message("[easy_scripts] STEP5_OUT -> ", .step5_tmp,
          " (no pipeline output; per-model CSVs go to a temp dir)")
}

# Clean up helper-load signals so subsequent calls behave normally
.cleanup()

invisible(NULL)
