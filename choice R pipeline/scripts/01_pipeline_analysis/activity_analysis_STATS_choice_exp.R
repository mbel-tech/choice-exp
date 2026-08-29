# =============================================================================
# STEP 5 — CHOICE EXPERIMENT  |  LMM STATISTICS + SKINNY GRAPHS
# Biological unit: trial × timepoint (school). Identity-free indicators only —
# fish_uid REs are removed because idtracker.ai cannot maintain persistent fish
# identities in choice arenas.
#
# PRIMARY ANALYSES (Module A — trial × timepoint):
#   A1  Main-zone occupancy (flow/calm)   × Treatment × Zone × Timepoint
#   A2  Sub-zone occupancy (h/m/l/calm)   × Treatment × Zone × Timepoint
#   A3  switches_per_session              × Treatment × Timepoint
#   A4  switches_per_session              × Treatment   (aggregated)
#   A5  prop_active (>1 BL/s, school)     × Treatment × Timepoint
#   A6  prop_active                       × Treatment   (aggregated)
#   A7-A16 group dynamics (NND, IID, hull area, centroid speed)
#
# WORKFLOW PER ANALYSIS:
#   1. Check normality of residuals (Shapiro-Wilk on lm residuals).
#   2. Check equality of variance (Levene's test).
#   3. If not normal: try log1p, sqrt; use best-transforming version.
#      (A3/A4 skip this — handled by GLMM family choice.)
#   4. Fit LMM candidates (lme4::lmer) with AICc model selection.
#      For A2/A4/A5: also fit a candidate with continuous timepoint predictor.
#      For A3/A4: use glmer (Poisson) or glmer.nb (if overdispersed).
#   5. Refit best model with REML = TRUE before inference.
#   6. Type III Wald chi-sq ANOVA on best model.
#   7. If any term significant (p < 0.05):
#        Tukey post-hoc + CLD (compact letter display) on each significant term.
#   8. Save residual plot + QQ plot.
#   9. Produce chi-sq(df) = value, p = value caption string per focal term.
#
# PSEUDOREPLICATION CHECK:
#   For A1 (zone) and A6 (moving), a parallel tank-level mean LM is run.
#   Treatment p-values are compared; concordance printed.
#
# MULTIPLE COMPARISONS:
#   No adjustment is applied anywhere in this study (author decision, 2026-08-07;
#   see Supplementary Table S7/2.7). All p-values are raw.
#
# DESIGN NOTE:
#   Fish density refers to housing-tank density before trials, not arena density
#   (arena always has 5 fish). It is therefore included as a random effect
#   candidate, not a fixed effect, so models test its contribution to variance
#   rather than treating it as an experimental manipulation.
#
# NOTE: idtracker.ai assigns fish IDs independently per segment.
# Tank is the unit repeated across timepoints; (1|tank) is the repeated-measures RE.
#
# OUTPUTS:
#   STEP5_stats/STEP5_stats_<ts>/
#     zone_aggregated/ zone_timepoint/
#     switches_aggregated/ switches_timepoint/
#     moving_timepoint/ moving_aggregated/
#     pseudoreplication_check.csv
#     bh_correction_summary.csv
#     skinny_graphs/
#       skinny_graphs_aggregated.pptx
#       skinny_graphs_timepoint.pptx
#       png_preview/
#   analysis_report_choice_exp.docx   (pipeline folder)
# =============================================================================


# =============================================================================
# ==== 0) SETUP ===============================================================
# =============================================================================

options(stringsAsFactors = FALSE, na.action = "na.omit")
set.seed(42)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggtext)
  library(lme4)
  library(MASS)
  library(MuMIn)
  library(car)
  library(emmeans)
  library(multcomp)
  library(multcompView)
  library(officer)
  library(rvg)
  library(flextable)
  library(readr)
  library(readxl)
  library(stringr)
  library(purrr)
  library(ggsignif)
})

# Optional packages: lmerTest + pbkrtest = Kenward-Roger df; glmmTMB = beta GLMM;
# performance = R²m/R²c. All gated behind capability flags so missing packages
# downgrade to current behaviour.
.HAS_KR     <- requireNamespace("pbkrtest",    quietly = TRUE) &&
               requireNamespace("lmerTest",   quietly = TRUE)
.HAS_TMB    <- requireNamespace("glmmTMB",    quietly = TRUE)
.HAS_PERF   <- requireNamespace("performance",quietly = TRUE)
.HAS_DHARMA <- requireNamespace("DHARMa",     quietly = TRUE)
# effectsize = noncentral-F CI on eta2_p/omega2_p (Nakagawa & Cuthill 2007).
.HAS_EFFECTSIZE <- requireNamespace("effectsize", quietly = TRUE)
# Dirichlet GLMM via brms (requires rstan/cmdstanr + Stan toolchain).
# Install: install.packages("brms")  then ensure rstan or cmdstanr is set up.
.HAS_BRMS   <- requireNamespace("brms",       quietly = TRUE)
if (.HAS_KR)     suppressPackageStartupMessages(library(lmerTest))
if (.HAS_TMB)    suppressPackageStartupMessages(library(glmmTMB))
if (.HAS_PERF)   suppressPackageStartupMessages(library(performance))
if (.HAS_EFFECTSIZE) suppressPackageStartupMessages(library(effectsize))
if (.HAS_BRMS)   suppressPackageStartupMessages(library(brms))

select    <- dplyr::select
filter    <- dplyr::filter
mutate    <- dplyr::mutate
summarise <- dplyr::summarise

ts_msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), " | ", ..., "\n", sep = "")
`%||%` <- function(a, b) if (!is.null(a)) a else b

# N-per-group caption text for graph footers.
# type = "agg" → unit is "Trial"; type = "tp" → unit is "Trial × Timepoint combination"
# df must contain columns: treatment, trial_id (and timepoint_f for "tp").
.n_cap <- function(df, type = c("agg", "tp")) {
  type <- match.arg(type)
  if (is.null(df) || !"treatment" %in% names(df)) return(NULL)
  tryCatch({
    if (type == "agg") {
      ns <- df %>%
        dplyr::group_by(treatment) %>%
        dplyr::summarise(
          n = if ("phys_trial_id" %in% names(df))
                dplyr::n_distinct(phys_trial_id)
              else
                dplyr::n(),
          .groups = "drop")
      unit <- "Trial"
    } else {
      ns <- df %>%
        dplyr::group_by(treatment) %>%
        dplyr::summarise(
          n = dplyr::n_distinct(interaction(trial_id, timepoint_f, drop = TRUE)),
          .groups = "drop")
      unit <- "Trial × Timepoint combination"
    }
    n_uniq <- dplyr::n_distinct(ns$n)
    n_str  <- if (n_uniq == 1) {
      paste0("N = ", ns$n[1], " per group")
    } else {
      paste(sprintf("%s: N = %d", ns$treatment, ns$n), collapse = "; ")
    }
    paste0(n_str, "; one datapoint = one ", unit, ".")
  }, error = function(e) NULL)
}

# Drop columns where every non-NA value is identical.
# Returns list(df = trimmed_df, notes = character vector of descriptive lines).
.drop_constant_cols <- function(df) {
  notes <- character(0)
  drop_cols <- character(0)
  first_col <- if (ncol(df) > 0) names(df)[1] else "group"
  for (cn in names(df)) {
    u <- unique(df[[cn]][!is.na(df[[cn]])])
    if (length(u) == 1) {
      notes <- c(notes, paste0("All ", first_col, " produced “", u, "” in ", cn, "."))
      drop_cols <- c(drop_cols, cn)
    }
  }
  list(df = df[, setdiff(names(df), drop_cols), drop = FALSE], notes = notes)
}

.get_global <- function(name, default = NULL) {
  if (exists(name, envir = .GlobalEnv, inherits = FALSE))
    get(name, envir = .GlobalEnv, inherits = FALSE)
  else default
}

.step5_dir <- tryCatch({
  frames     <- sys.frames()
  ofile_envs <- Filter(function(f) exists("ofile", envir = f, inherits = FALSE), frames)
  if (length(ofile_envs) > 0)
    normalizePath(dirname(get("ofile", envir = ofile_envs[[length(ofile_envs)]])),
                  winslash = "/", mustWork = FALSE)
  else getwd()
}, error = function(e) getwd())

TREATMENT_LEVELS_g  <- .get_global("TREATMENT_LEVELS",  c("control", "exercise choice"))
DENSITY_LEVELS_g    <- .get_global("DENSITY_LEVELS",    c(4L, 8L, 12L, 16L))
DENSITY_AS_FACTOR_g <- .get_global("DENSITY_AS_FACTOR", TRUE)
TIMEPOINT_LEVELS_g  <- .get_global("TIMEPOINT_LEVELS",  1:3)

# Parametric-bootstrap LRT replicates for beta-GLMM treatment tests (df_method_memo
# section 4: KR/Satterthwaite are undefined for GLMMs; pbkrtest::PBmodcomp has no
# glmmTMB method, confirmed 2026-08-08, so this project implements the same idea
# manually via simulate.glmmTMB).
#
# 2026-08-18: default lowered 1000 -> 200. Timing of the 2026-08-17 run showed
# the beta-GLMM bootstraps dominate the whole step: the 9 jacobs_beta_* models
# took 47.7% of a 53-minute run and zone_sec_cells_aggregated another 22.6%,
# while all 44 remaining model folders together took 2.8%. At 200 replicates
# those models run ~5x faster. The cost is p-value RESOLUTION: the smallest
# reportable PB p-value goes from ~1/1001 to ~1/201 (~0.005), which matters
# only for a p sitting right on a threshold. Raise it back per-run with
# N_PB <- 1000L in the calling environment before sourcing this script when
# publishing a number that depends on that precision.
N_PB_DEFAULT <- .get_global("N_PB", 200L)
PB_SEED      <- .get_global("PB_SEED", 20260808L)

.find_latest_csv <- function(step_name, csv_filename) {
  # Search getwd() first, then parent of getwd() (handles pipeline run from subdirectory)
  for (.base in c(getwd(), dirname(getwd()))) {
    parent <- file.path(.base, step_name)
    if (!dir.exists(parent)) next
    subdirs <- list.dirs(parent, full.names = TRUE, recursive = FALSE)
    subdirs <- subdirs[grepl(paste0("^", step_name, "_\\d{8}_\\d{6}$"), basename(subdirs))]
    if (length(subdirs) == 0) next
    latest    <- subdirs[which.max(file.mtime(subdirs))]
    candidate <- file.path(latest, csv_filename)
    if (file.exists(candidate)) return(candidate)
  }
  NULL
}

.step5_parent <- file.path(getwd(), "STEP5_stats")
if (!dir.exists(.step5_parent))
  dir.create(.step5_parent, recursive = TRUE, showWarnings = FALSE)
# Allow callers to pre-set STEP5_OUT (e.g. for manu_graphs-only regeneration).
if (!exists("STEP5_OUT") || !nzchar(STEP5_OUT) || !dir.exists(STEP5_OUT)) {
  STEP5_OUT <- file.path(.step5_parent,
                          paste0("STEP5_stats_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  dir.create(STEP5_OUT, recursive = TRUE, showWarnings = FALSE)
}
ts_msg("Output root: ", STEP5_OUT)


# =============================================================================
# ==== 1) LOAD DATA ===========================================================
# =============================================================================

# School-level (trial × timepoint) refactor: load trial_activity_summary +
# trial_occupancy_long produced by the new STEP2.
if (!exists("trial_activity_summary", envir = .GlobalEnv, inherits = FALSE)) {
  ts_msg("trial_activity_summary not in GlobalEnv; searching disk...")
  .step2_out <- .get_global("STEP2_OUTPUT_DIR")
  .csv_path  <- if (!is.null(.step2_out) &&
                    file.exists(file.path(.step2_out, "trial_activity_summary.csv")))
    file.path(.step2_out, "trial_activity_summary.csv") else NULL
  if (is.null(.csv_path))
    .csv_path <- .find_latest_csv("STEP2_output", "trial_activity_summary.csv")
  if (is.null(.csv_path))
    stop("trial_activity_summary not found. Run Steps 1-2 first.", call. = FALSE)
  trial_activity_summary <- readr::read_csv(.csv_path, show_col_types = FALSE)
  assign("trial_activity_summary", trial_activity_summary, envir = .GlobalEnv)
  ts_msg("Loaded: ", .csv_path, " (", nrow(trial_activity_summary), " rows)")
}
df <- get("trial_activity_summary", envir = .GlobalEnv, inherits = FALSE)

if (!exists("trial_occupancy_long", envir = .GlobalEnv, inherits = FALSE)) {
  .occ_path <- .find_latest_csv("STEP2_output", "trial_occupancy_long.csv")
  if (!is.null(.occ_path)) {
    trial_occupancy_long <- readr::read_csv(.occ_path, show_col_types = FALSE)
    assign("trial_occupancy_long", trial_occupancy_long, envir = .GlobalEnv)
  } else {
    stop("trial_occupancy_long not found. Re-run Step 2.", call. = FALSE)
  }
}
df_long <- get("trial_occupancy_long", envir = .GlobalEnv, inherits = FALSE)

# Timepoint guard
if (!"timepoint" %in% names(df)) {
  warning("'timepoint' column missing — timepoint analyses will be skipped.", call. = FALSE)
  df$timepoint <- NA_integer_
}

# Factor coercions (trial-level: no fish_uid)
df <- df %>%
  dplyr::mutate(
    treatment      = factor(trimws(tolower(as.character(treatment))),
                            levels = TREATMENT_LEVELS_g),
    fish_density_f = factor(fish_density, levels = DENSITY_LEVELS_g,
                             ordered = isTRUE(DENSITY_AS_FACTOR_g)),
    timepoint      = as.integer(timepoint),
    timepoint_f    = factor(timepoint, levels = TIMEPOINT_LEVELS_g),
    trial_id       = as.integer(trial_id),
    tank           = as.character(tank)
  )
df_long <- df_long %>%
  dplyr::mutate(
    treatment      = factor(trimws(tolower(as.character(treatment))),
                            levels = TREATMENT_LEVELS_g),
    timepoint      = as.integer(timepoint),
    timepoint_f    = factor(timepoint, levels = TIMEPOINT_LEVELS_g),
    trial_id       = as.integer(trial_id),
    tank           = as.character(tank),
    zone_level     = factor(zone_level, levels = c("main", "sec")),
    zone           = factor(zone, levels = c("flow", "calm", "high", "medium", "low"))
  )

# Physical trial ID: groups the 3 timepoint sessions from the same physical trial
# (tank × trial_date × treatment). Used as RE cluster in timepoint analyses.
if ("trial_date" %in% names(df)) {
  df <- df %>%
    dplyr::mutate(phys_trial_id = as.integer(
      factor(paste(tank, as.character(trial_date), treatment))))
} else {
  df$phys_trial_id <- df$trial_id  # fallback if date absent
}

ts_msg("Data: ", nrow(df), " rows | ",
       dplyr::n_distinct(df$tank), " tanks | ",
       dplyr::n_distinct(df$timepoint[!is.na(df$timepoint)]), " timepoints | ",
       dplyr::n_distinct(df$phys_trial_id), " physical trials")
readr::write_csv(df, file.path(STEP5_OUT, "analysis_ready.csv"))


# =============================================================================
# ==== 2) BUILD LONG-FORMAT DATASETS ==========================================
# =============================================================================
# trial_occupancy_long has columns: trial_id, tank, treatment, timepoint,
#   zone_level (main/sec), zone (flow/calm/high/medium/low), prop_time.
# Split into main- and sec-zone subsets for separate fits.

# Ensure metadata columns exist on df_long; backfill from df where missing.
.meta_cols <- intersect(c("fish_density", "fish_density_f", "trial_date"), names(df))
if (length(.meta_cols) > 0 &&
    !all(.meta_cols %in% names(df_long))) {
  .meta <- dplyr::distinct(
    df[, c("trial_id", "timepoint", .meta_cols)],
    trial_id, timepoint, .keep_all = TRUE
  )
  df_long <- dplyr::left_join(df_long, .meta, by = c("trial_id", "timepoint"))
}
if ("fish_density" %in% names(df_long)) {
  df_long$fish_density_f <- factor(df_long$fish_density, levels = DENSITY_LEVELS_g)
} else {
  df_long$fish_density_f <- NA
}

df_main_long <- dplyr::filter(df_long, zone_level == "main") %>%
  dplyr::mutate(zone = factor(as.character(zone), levels = c("flow", "calm")))

df_sec_long <- dplyr::filter(df_long, zone_level == "sec") %>%
  dplyr::mutate(zone = factor(as.character(zone),
                              levels = c("high", "medium", "low", "calm")))

# Propagate phys_trial_id to long format
if ("phys_trial_id" %in% names(df)) {
  .phys_map <- dplyr::distinct(df[, c("trial_id", "phys_trial_id")])
  df_long <- dplyr::left_join(df_long, .phys_map, by = "trial_id")
  df_main_long <- dplyr::left_join(df_main_long, .phys_map, by = "trial_id")
  df_sec_long  <- dplyr::left_join(df_sec_long,  .phys_map, by = "trial_id")
}

ts_msg("Long-format: main=", nrow(df_main_long),
       " | sec=", nrow(df_sec_long))

# Aggregated long-format: average prop_time across timepoints per physical trial × zone.
# Physical trial = tank × trial_date × treatment (N = 16: 4 tanks × 2 dates × 2 treatments).
# Used for A1_agg / A2_agg (treatment × zone, no timepoint term).
.agg_long_cols <- c("phys_trial_id", "zone", "treatment", "tank",
                    "fish_density_f", "trial_date")
df_main_long_agg <- df_main_long %>%
  dplyr::mutate(phys_trial_id = as.integer(factor(paste(tank, trial_date, treatment)))) %>%
  dplyr::group_by(dplyr::across(dplyr::any_of(.agg_long_cols))) %>%
  dplyr::summarise(prop_time = mean(prop_time, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(zone = factor(as.character(zone), levels = levels(df_main_long$zone)))

df_sec_long_agg <- df_sec_long %>%
  dplyr::mutate(phys_trial_id = as.integer(factor(paste(tank, trial_date, treatment)))) %>%
  dplyr::group_by(dplyr::across(dplyr::any_of(.agg_long_cols))) %>%
  dplyr::summarise(prop_time = mean(prop_time, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(zone = factor(as.character(zone), levels = levels(df_sec_long$zone)))

# =============================================================================
# COMPOSITIONALLY CORRECT ZONE DATASETS
# =============================================================================
# Why: prop_flow + prop_calm = 1 (main zones) and sum(sub-zone props) = 1.
# Fitting beta GLMMs or Gaussian LMMs on the long-format stacked data treats
# zone rows within a session as independent — they are not (algebraic dependency).
# Solution:
#   Main zones  → scalar logit(prop_flow) per session / physical trial.
#   Sub-zones   → three log-ratios vs calm (log(prop_z / prop_calm)) per session /
#                 physical trial. Calm is the reference simplex component and is
#                 implicitly modelled. ILR basis would be formally correct but
#                 harder to interpret; log-ratios (pairwise vs calm) are sufficient
#                 when calm is a natural reference and the principal research
#                 question is directional (treatment pushes fish away from calm).
#   Zeros: a small Haldane-type offset (.eps_lr = 1e-4, ~0.01% of session time)
#          is added to all four sub-zone proportions before log-ratio computation
#          to stabilise the ratio when any zone is never used. No offset is needed
#          for the main-zone logit because rows with prop_flow ∈ {0,1} are removed
#          (they arise from sessions where fish never used one zone at all; n very
#          small, removing them prevents infinite logit and is documented below).

.eps_lr <- 1e-4   # Haldane offset for sub-zone log-ratios

# — Session-level main zone (N = 48) ——————————————————————————————————————————
df_main_wide <- tryCatch({
  df_main_long %>%
    dplyr::filter(zone %in% c("flow", "calm")) %>%
    dplyr::select(dplyr::any_of(c("trial_id", "phys_trial_id", "tank", "trial_date",
                                   "treatment", "timepoint", "timepoint_f",
                                   "fish_density", "fish_density_f", "zone", "prop_time"))) %>%
    tidyr::pivot_wider(names_from = zone, values_from = prop_time, names_prefix = "prop_") %>%
    dplyr::filter(is.finite(prop_flow)) %>%
    dplyr::mutate(
      prop_flow_sq = pmin(pmax(prop_flow, 1e-4), 1 - 1e-4),
      logit_flow   = log(prop_flow_sq / (1 - prop_flow_sq))
    )
}, error = function(e) { warning("df_main_wide failed: ", e$message); NULL })

# — Physical-trial-level main zone (N = 16) ———————————————————————————————————
df_main_wide_agg <- tryCatch({
  df_main_long_agg %>%
    dplyr::filter(zone %in% c("flow", "calm")) %>%
    dplyr::select(dplyr::any_of(c("phys_trial_id", "tank", "trial_date", "treatment",
                                   "fish_density_f", "zone", "prop_time"))) %>%
    tidyr::pivot_wider(names_from = zone, values_from = prop_time, names_prefix = "prop_") %>%
    dplyr::filter(is.finite(prop_flow) & prop_flow > 0 & prop_flow < 1) %>%
    dplyr::mutate(logit_flow = log(prop_flow / (1 - prop_flow)))
}, error = function(e) { warning("df_main_wide_agg failed: ", e$message); NULL })

# — Session-level sub-zones (N = 48) ——————————————————————————————————————————
df_sec_wide <- tryCatch({
  df_sec_long %>%
    dplyr::filter(zone %in% c("high", "medium", "low", "calm")) %>%
    dplyr::select(dplyr::any_of(c("trial_id", "phys_trial_id", "tank", "trial_date",
                                   "treatment", "timepoint", "timepoint_f",
                                   "fish_density", "fish_density_f", "zone", "prop_time"))) %>%
    tidyr::pivot_wider(names_from = zone, values_from = prop_time, names_prefix = "prop_") %>%
    dplyr::mutate(
      # DO NOT area-normalize here. The sub-zone proportions arrive from STEP2
      # section 6b ALREADY divided by their area constants (high/10, medium/18,
      # low/13, calm/42) and rescaled to sum to 1 (prop_time_ac_*), which is what
      # trial_occupancy_long carries as prop_time for zone_level == "sec".
      #
      # THE 2026-08-18 DOUBLE-AREA-NORMALISATION FIX
      #   Until 2026-08-18 this block divided by the area constants a SECOND
      #   time, so every sub-zone log-ratio carried a zone-specific offset of
      #   log(A_calm / A_zone): high +1.435, medium +0.847, low +1.173.
      #   Consequences of the double correction:
      #     - treatment, interval and treatment x interval tests WITHIN a zone
      #       were unaffected (a constant per zone cancels in every contrast);
      #     - the zero point did NOT mean "used exactly in proportion to area",
      #       so the sign of an emmean was not interpretable as preference;
      #     - the stacked zone-vs-zone contrasts (A2_pairs) WERE biased, by the
      #       differences between those offsets (e.g. high vs medium by 0.588);
      #     - the Haldane floor did not mean what its comment says: applied to
      #       prop/A it clamped at A x 1e-4 of occupancy (0.42% of session time
      #       for calm), not at 1e-4.
      #   The area correction belongs in STEP2 and is applied there once.
      prop_high_an   = pmax(prop_high,   .eps_lr),
      prop_medium_an = pmax(prop_medium, .eps_lr),
      prop_low_an    = pmax(prop_low,    .eps_lr),
      prop_calm_an   = pmax(prop_calm,   .eps_lr),
      lr_high   = log(prop_high_an   / prop_calm_an),
      lr_medium = log(prop_medium_an / prop_calm_an),
      lr_low    = log(prop_low_an    / prop_calm_an)
    )
}, error = function(e) { warning("df_sec_wide failed: ", e$message); NULL })

# — Physical-trial-level sub-zones (N = 16) ———————————————————————————————————
df_sec_wide_agg <- tryCatch({
  df_sec_long_agg %>%
    dplyr::filter(zone %in% c("high", "medium", "low", "calm")) %>%
    dplyr::select(dplyr::any_of(c("phys_trial_id", "tank", "trial_date", "treatment",
                                   "fish_density_f", "zone", "prop_time"))) %>%
    tidyr::pivot_wider(names_from = zone, values_from = prop_time, names_prefix = "prop_") %>%
    dplyr::mutate(
      # DO NOT area-normalize here. The sub-zone proportions arrive from STEP2
      # section 6b ALREADY divided by their area constants (high/10, medium/18,
      # low/13, calm/42) and rescaled to sum to 1 (prop_time_ac_*), which is what
      # trial_occupancy_long carries as prop_time for zone_level == "sec".
      #
      # THE 2026-08-18 DOUBLE-AREA-NORMALISATION FIX
      #   Until 2026-08-18 this block divided by the area constants a SECOND
      #   time, so every sub-zone log-ratio carried a zone-specific offset of
      #   log(A_calm / A_zone): high +1.435, medium +0.847, low +1.173.
      #   Consequences of the double correction:
      #     - treatment, interval and treatment x interval tests WITHIN a zone
      #       were unaffected (a constant per zone cancels in every contrast);
      #     - the zero point did NOT mean "used exactly in proportion to area",
      #       so the sign of an emmean was not interpretable as preference;
      #     - the stacked zone-vs-zone contrasts (A2_pairs) WERE biased, by the
      #       differences between those offsets (e.g. high vs medium by 0.588);
      #     - the Haldane floor did not mean what its comment says: applied to
      #       prop/A it clamped at A x 1e-4 of occupancy (0.42% of session time
      #       for calm), not at 1e-4.
      #   The area correction belongs in STEP2 and is applied there once.
      prop_high_an   = pmax(prop_high,   .eps_lr),
      prop_medium_an = pmax(prop_medium, .eps_lr),
      prop_low_an    = pmax(prop_low,    .eps_lr),
      prop_calm_an   = pmax(prop_calm,   .eps_lr),
      lr_high   = log(prop_high_an   / prop_calm_an),
      lr_medium = log(prop_medium_an / prop_calm_an),
      lr_low    = log(prop_low_an    / prop_calm_an)
    )
}, error = function(e) { warning("df_sec_wide_agg failed: ", e$message); NULL })

ts_msg("Wide-format datasets: main_wide=", if (!is.null(df_main_wide)) nrow(df_main_wide) else "NULL",
       " | main_wide_agg=", if (!is.null(df_main_wide_agg)) nrow(df_main_wide_agg) else "NULL",
       " | sec_wide=", if (!is.null(df_sec_wide)) nrow(df_sec_wide) else "NULL",
       " | sec_wide_agg=", if (!is.null(df_sec_wide_agg)) nrow(df_sec_wide_agg) else "NULL")

# Wide-format aggregated: one row per physical trial (tank × trial_date × treatment).
# Averages rate/proportion columns across timepoints; sums count/time columns.
# N = 16 (4 tanks × 2 recording dates × 2 treatments).
# Used for all A*_agg analyses (no timepoint term).
{
  .agg_sum_cols  <- intersect(c("n_main_switches", "total_time_s"), names(df))
  .agg_mean_cols <- intersect(c("switches_per_session", "prop_active",
                                 "zone_flux_per_session"), names(df))
  df_agg <- df %>%
    dplyr::mutate(phys_trial_id = as.integer(factor(paste(tank, trial_date, treatment)))) %>%
    dplyr::group_by(phys_trial_id, tank, trial_date, treatment, fish_density, fish_density_f) %>%
    dplyr::summarise(
      dplyr::across(dplyr::all_of(.agg_sum_cols),  ~ sum(.x,  na.rm = TRUE)),
      dplyr::across(dplyr::all_of(.agg_mean_cols), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      treatment      = factor(treatment,      levels = levels(df$treatment)),
      fish_density_f = factor(fish_density_f, levels = levels(df$fish_density_f))
    )
  ts_msg("df_agg: ", nrow(df_agg), " physical trials | ",
         dplyr::n_distinct(df_agg$treatment), " treatments")
}

# =============================================================================
# ==== 2b) DESIGN BALANCE CHECK (F8) ==========================================
# =============================================================================
# Detect treatment × tank × density × date confounding and design structure.
# Design facts (choice exp):
#   4 tanks; each tank ran on 2 consecutive dates → 8 unique dates total,
#   nested within tank. Each (tank, date) carries 1 control + 1 exercise
#   trial → 16 physical trials. Each trial has 3 timepoint sessions → 48 rows.
{
  .bal <- tryCatch({
    .keep_cols <- intersect(c("treatment","tank","fish_density","trial_date"),
                            names(df))
    .d_meta <- df[!duplicated(df$trial_id), .keep_cols]
    .d_meta$tank_num <- as.numeric(as.factor(.d_meta$tank))

    .xt_td  <- as.data.frame(xtabs(~ treatment + tank,         data = .d_meta))
    .xt_dn  <- as.data.frame(xtabs(~ treatment + fish_density, data = .d_meta))
    .xt_all <- as.data.frame(xtabs(~ treatment + tank + fish_density, data = .d_meta))

    # Tank ↔ density: are they 1:1?
    .td_map <- unique(.d_meta[, c("tank","fish_density")])
    .tank_density_1to1 <- nrow(.td_map) == length(unique(.d_meta$tank))

    # Correlation between tank (numeric) and density
    .cor_td <- tryCatch(
      cor(.d_meta$tank_num, as.numeric(as.factor(.d_meta$fish_density))),
      error = function(e) NA_real_)

    # Date structure
    .n_dates       <- if ("trial_date" %in% .keep_cols)
                        dplyr::n_distinct(.d_meta$trial_date) else NA_integer_
    .date_per_tank <- if ("trial_date" %in% .keep_cols)
                        nrow(unique(.d_meta[, c("tank","trial_date")])) /
                        length(unique(.d_meta$tank)) else NA_real_
    # date nested in tank? TRUE if every date appears under exactly one tank
    .date_nested_in_tank <- if ("trial_date" %in% .keep_cols) {
      .dt <- unique(.d_meta[, c("tank","trial_date")])
      max(table(.dt$trial_date)) == 1L
    } else NA
    # Treatment x date balance (within-date treatment crossing)
    .xt_tr_date <- if ("trial_date" %in% .keep_cols)
      as.data.frame(xtabs(~ treatment + trial_date, data = .d_meta)) else NULL
    .treat_date_balanced <- if (!is.null(.xt_tr_date)) {
      all(.xt_tr_date$Freq == .xt_tr_date$Freq[1])
    } else NA

    list(treatment_tank = .xt_td, treatment_density = .xt_dn,
         treatment_date = .xt_tr_date,
         all = .xt_all, tank_density_1to1 = .tank_density_1to1,
         cor_tank_density = .cor_td,
         n_dates = .n_dates,
         date_per_tank = .date_per_tank,
         date_nested_in_tank = .date_nested_in_tank,
         treat_date_balanced = .treat_date_balanced)
  }, error = function(e) {
    warning("Balance check failed: ", conditionMessage(e)); NULL
  })

  if (!is.null(.bal)) {
    cat("\n=== DESIGN BALANCE CHECK ===\n")
    cat("Treatment x Tank:\n");    print(xtabs(~ treatment + tank, data = df[!duplicated(df$trial_id),]))
    cat("Treatment x Density:\n"); print(xtabs(~ treatment + fish_density, data = df[!duplicated(df$trial_id),]))
    if ("trial_date" %in% names(df)) {
      cat("Treatment x Date:\n");
      print(xtabs(~ treatment + trial_date, data = df[!duplicated(df$trial_id),]))
      cat("Tank x Date:\n");
      print(xtabs(~ tank + trial_date, data = df[!duplicated(df$trial_id),]))
      cat("N dates =", .bal$n_dates,
          "| dates per tank =", .bal$date_per_tank,
          "| date nested in tank:", .bal$date_nested_in_tank,
          "| treatment x date balanced:", .bal$treat_date_balanced, "\n")
    }
    if (isTRUE(.bal$tank_density_1to1))
      cat("WARNING: tank <-> density mapping is 1:1 (", paste(.bal$treatment_tank$treatment, collapse=","), "). These two REs are collinear.\n")
    if (!is.na(.bal$cor_tank_density) && abs(.bal$cor_tank_density) > 0.9)
      cat("WARNING: tank-density correlation =", round(.bal$cor_tank_density, 3),
          " > 0.9. Density RE absorbs tank RE; include only one.\n")
    if (isTRUE(.bal$date_nested_in_tank))
      cat("NOTE: trial_date is fully nested in tank (each date occurs under exactly one tank). ",
          "Use (1|tank/trial_date) rather than (1|tank) + (1|trial_date) to avoid collinear REs.\n",
          sep = "")
    cat("=============================\n\n")
    readr::write_csv(
      data.frame(
        check = c("tank_density_1to1", "cor_tank_density",
                  "n_dates", "dates_per_tank",
                  "date_nested_in_tank", "treatment_date_balanced",
                  "n_tanks", "n_phys_trials", "n_sessions"),
        value = c(as.character(.bal$tank_density_1to1),
                  round(.bal$cor_tank_density, 3),
                  as.character(.bal$n_dates),
                  as.character(.bal$date_per_tank),
                  as.character(.bal$date_nested_in_tank),
                  as.character(.bal$treat_date_balanced),
                  as.character(dplyr::n_distinct(df$tank)),
                  as.character(dplyr::n_distinct(
                    paste(df$tank, df$trial_date, df$treatment))),
                  as.character(nrow(df[!duplicated(df$trial_id), ])))
      ),
      file.path(STEP5_OUT, "design_balance.csv")
    )
  }
}


# =============================================================================
# ==== 3) RE CANDIDATES =======================================================
# =============================================================================
# fish_density is a housing-tank property (not an arena variable; arena always
# has 5 fish). It is included only as a RE candidate so models can account for
# variance attributable to pre-trial housing conditions.

# Trial × timepoint RE candidates (N = 48 sessions, 3 per physical trial).
#
# CHANGED 2026-08-08 (design-based random-effect policy; Reviewer 1 comments
# 6 & 7; Barr et al. 2013; Harrison et al. 2018): phys_trial_id is now FORCED
# into every candidate rather than competing against tank/date/density for
# AICc's favour. The 3 sessions per physical trial are genuine repeated
# measures on the same school -- that clustering is imposed by the design,
# not an empirical question AICc should be allowed to answer. Before this
# change, AICc dropped phys_trial_id entirely for several interval-level
# outcomes (e.g. some collective outcomes selected (1|trial_date) alone; zone-switch
# counts selected (1|tank) alone), leaving the within-trial clustering
# unmodelled. AICc is now reserved for the genuinely optional nuisance terms
# (tank, trial_date). "density" is removed as a random-effect candidate
# entirely: fish_density is a 4-level ORDERED DESIGN COVARIATE, not a
# sampling-unit grouping factor, and it is already fitted as a fixed term
# (see fixed_str at each call site) -- fitting the same variable as both a
# fixed covariate and a random grouping factor double-counts it.
# trial_id is session-level (row identifier) and must not be used as RE.
RE_WIDE <- c(
  phys_trial   = "(1 | phys_trial_id)",
  trial_tank   = "(1 | phys_trial_id) + (1 | tank)",
  trial_date_r = "(1 | phys_trial_id) + (1 | trial_date)"
)

# Long-format: zone rows from same session linked within physical trial.
# Unused by any run_lmm_analysis() call as of 2026-08-08 (referenced only in
# the Word-report methods-decisions text); kept trial-forced for consistency
# with RE_WIDE if it is ever wired up.
RE_LONG <- c(
  phys_trial   = "(1 | phys_trial_id)",
  trial_tank   = "(1 | phys_trial_id) + (1 | tank)"
)

# Aggregated-level RE candidates (N = 16 physical trials): INTENTIONALLY EMPTY.
#
# CHANGED 2026-08-08: at this level each physical trial contributes exactly
# one row, so a trial random intercept is not identifiable (it is confounded
# with the residual), and the previous (1|tank) candidate put a random
# intercept on a 4-level factor -- below the conventional 5-8-level floor for
# a stable variance estimate (Bolker et al. 2009; Harrison et al. 2018).
# Verified against the realised design (Table S1): fish_density is CROSSED
# with tank (a complete 4x4 grid, one trial per cell), not nested in it as an
# earlier code comment claimed, so it is not redundant with tank and both are
# now fitted as FIXED covariates at every "_aggregated" call site instead
# (factor(tank) + fish_density_f). An empty candidate set here routes through
# the existing re_fallback -> ordinary-least-squares path in
# run_lmm_analysis() unchanged; treatment remains orthogonal to both tank and
# density, so this does not alter the treatment estimate.
RE_WIDE_AGG <- character(0)
RE_LONG_AGG <- character(0)


# =============================================================================
# ==== 4) UTILITY FUNCTIONS ===================================================
# =============================================================================

# ---- Unified number-display convention (author decision, 2026-08-07) -------
# fmt_F: F/chi-square test statistics -> 4 significant figures, trailing
#        zeros trimmed (10.348 -> "10.35"; 9.724 stays "9.724").
# fmt3 : every other statistical value (p once below the "< 0.001" floor,
#        partial eta-squared, partial omega-squared, Cohen's d, R^2, ICC) ->
#        3 decimals, dropped to 2 ONLY when the 3rd decimal digit is exactly
#        zero (0.150 -> "0.15"; 0.019 stays "0.019"). See
#        00_shared/number_formatting.R for the full rationale; duplicated here
#        verbatim because this script is run standalone.
fmt_F <- function(x) {
  if (!is.finite(x)) return("NA")
  s <- format(signif(x, 4), scientific = FALSE, trim = TRUE)
  # Only trim trailing zeros when a decimal point is present -- otherwise an
  # integer denominator df like 30 is mangled to "3" (see project number-
  # formatting defect log).
  if (grepl(".", s, fixed = TRUE)) {
    s <- sub("0+$", "", s)
    s <- sub("\\.$", "", s)
  }
  s
}
fmt3 <- function(x) {
  if (!is.finite(x)) return("NA")
  s <- sprintf("%.3f", x)
  if (grepl("0$", s)) s <- substr(s, 1, nchar(s) - 1)
  s
}

# fmt_Fstat / fmt_es3 (2026-08-09): as fmt_F / fmt3, but a value below 0.001
# collapses to "< 0.001" instead of rendering either a long 4-sig-fig decimal
# ("0.0009258") or an uninformative rounded zero ("0.00"). Applied ONLY to the
# test statistic and to partial eta-squared -- denominator df keep fmt_F (a df
# is never < 1) and p-values keep fmt_p (which owns its own floor).
fmt_Fstat <- function(x) {
  if (!is.finite(x)) return("NA")
  if (abs(x) < 0.001) return("< 0.001")
  fmt_F(x)
}
fmt_es3 <- function(x) {
  if (!is.finite(x)) return("NA")
  if (abs(x) < 0.001) return("< 0.001")
  fmt3(x)
}

fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (round(p, 3) == 0) return("p < 0.001")
  paste0("p = ", fmt3(p))
}

fmt_num <- function(x, digits = 3) {
  if (is.null(x) || length(x) == 0 || is.na(x) || !is.finite(x)) return(NA_character_)
  trimws(sub("0+$", "", sub("\\.$", "",
    format(signif(x, digits), scientific = FALSE, trim = TRUE))))
}

# Build a single-term ANOVA caption string. Emits text for every term with a
# valid p-value, regardless of significance (Nakagawa & Cuthill 2007: effect
# size + CI is what makes a non-significant result "inconclusive" rather than
# "no effect"; gating on p < 0.05 discards exactly the information that makes
# null results interpretable). Statistic format follows the test type:
#   F-KR / F-SW \u2192 F(df1, df2) = X, p = Y     (df_denom required; falls back to F(df1) if absent)
#   Wald-chisq  \u2192 \u03c7\u00b2(df)      = X, p = Y
#   LRT-PB      \u2192 LRT = X, p = Y (parametric bootstrap, N = <df_denom>) -- used for
#                 beta-GLMM treatment terms, where KR/Satterthwaite are undefined
#                 (df_method_memo section 4); df_denom is repurposed to carry nsim_ok.
# Interaction terms render with "\u00d7" instead of ":". stat_type is read PER ROW
# (av$stat_type, set by .anova_to_legacy_df and overwritable per term -- e.g. a
# beta-GLMM's "treatment" row can carry "LRT-PB" while every other row in the
# same table keeps "Wald-chisq") if present, else from the explicit argument,
# else from attr(av, "stat_type"), defaulting to Wald-chisq.
anova_caption_str <- function(av, term_pattern, stat_type = NULL) {
  if (is.null(av)) return(NA_character_)
  row <- av[grepl(term_pattern, av$term, ignore.case = TRUE, perl = TRUE), ]
  if (nrow(row) == 0) return(NA_character_)
  row <- row[1, ]
  if (is.null(row$p_value) || length(row$p_value) == 0 || is.na(row$p_value)) return(NA_character_)

  row_stat_type <- if ("stat_type" %in% names(row) && length(row$stat_type) &&
                        !is.na(row$stat_type)) as.character(row$stat_type[1]) else NULL
  stat_type <- row_stat_type %||% stat_type %||% attr(av, "stat_type") %||% "Wald-chisq"
  num_fmt   <- fmt_F   # degrees of freedom use the same 4-sig-fig/trim rule
  # Test statistic at 4 significant figures (fmt_Fstat): 10.348 -> "10.35",
  # 9.724 stays "9.724" -- the project's single number-display convention --
  # with anything below 0.001 collapsing to "< 0.001" (2026-08-09).
  stat_val  <- fmt_Fstat(as.numeric(row$chisq))
  is_f_test <- isTRUE(grepl("^F[-_]", stat_type))
  is_pb_lrt <- isTRUE(identical(stat_type, "LRT-PB"))

  if (is_f_test) {
    df_num <- as.integer(row$df)
    df_den_raw <- if (!is.null(row$df_denom)) row$df_denom else NA_real_
    df_den <- if (length(df_den_raw) > 0 && !is.na(df_den_raw) &&
                  is.finite(as.numeric(df_den_raw))) num_fmt(df_den_raw) else NA_character_
    # Partial eta-squared for the ANOVA term, with a noncentral-F 95% CI via
    # effectsize::F_to_eta2 (Nakagawa & Cuthill 2007 \u00a7II.4). Falls back to the
    # point-estimate-only formula eta^2_p = (F*df1)/(F*df1+df2) -- no CI -- when
    # the effectsize package is unavailable. Reported only for F-tests (undefined
    # for the Wald chi-square branch below) and only when df_denom is available.
    eta_str <- ""
    if (!is.na(df_den)) {
      .f <- suppressWarnings(as.numeric(row$chisq))
      .d2 <- suppressWarnings(as.numeric(df_den_raw))
      if (is.finite(.f) && is.finite(.d2) && is.finite(df_num)) {
        if (isTRUE(.HAS_EFFECTSIZE)) {
          es <- tryCatch(effectsize::F_to_eta2(.f, df_num, .d2, ci = 0.95,
                                                alternative = "two.sided"),
                         error = function(e) NULL)
          if (!is.null(es) && nrow(es)) {
            # Point estimate only (no CI bracket) so the figure caption/
            # subtitle renders on a single line; the CI is still available
            # in the underlying anova.csv / Results text for readers who
            # want it.
            eta_str <- sprintf(", &eta;<sup>2</sup><sub>p</sub> = %s",
                                fmt_es3(es$Eta2_partial[1]))
          }
        }
        if (!nzchar(eta_str)) {
          .e <- (.f * df_num) / (.f * df_num + .d2)
          eta_str <- sprintf(", &eta;<sup>2</sup><sub>p</sub> = %s", fmt_es3(.e))
        }
      }
    }
    stat_str <- if (!is.na(df_den))
      sprintf("F<sub>%d, %s</sub> = %s, %s%s", df_num, df_den, stat_val,
              fmt_p(row$p_value), eta_str)
    else
      sprintf("F<sub>%d</sub> = %s, %s", df_num, stat_val, fmt_p(row$p_value))
  } else if (is_pb_lrt) {
    # KR/Satterthwaite are undefined for a beta GLMM (df_method_memo section 4);
    # pbkrtest::PBmodcomp has no glmmTMB method (confirmed 2026-08-08), so this
    # project runs its own parametric-bootstrap LRT (.pb_lrt_term()). df_denom
    # is repurposed here to carry the number of bootstrap replicates that
    # actually converged (nsim_ok), not a denominator df.
    .nsim <- suppressWarnings(as.integer(row$df_denom))
    stat_str <- sprintf("LRT = %s, %s (parametric bootstrap, N = %s)",
                        stat_val, fmt_p(row$p_value),
                        if (is.finite(.nsim)) as.character(.nsim) else "NA")
  } else {
    stat_str <- paste0("\u03c7\u00b2<sub>", as.integer(row$df), "</sub> = ",
                       stat_val, ", ", fmt_p(row$p_value))
  }

  term_str <- as.character(row$term[1])
  pretty_one <- function(x) {
    x <- sub("_f$", "", x)
    paste0(toupper(substring(x, 1, 1)), substring(x, 2))
  }
  pretty_term <- if (grepl(":", term_str, fixed = TRUE)) {
    # Interaction: replace ":" with " \u00d7 " in the rendered term
    paste(sapply(strsplit(term_str, ":")[[1]], pretty_one), collapse = " \u00d7 ")
  } else {
    pretty_one(term_str)
  }
  paste0(pretty_term, ": ", stat_str)
}

# Extract raw treatment p-value from an anova result data frame
.get_treatment_p <- function(res, treat_pat = "^treatment$") {
  if (is.null(res) || is.null(res$anova)) return(NA_real_)
  av  <- res$anova
  row <- av[grepl(treat_pat, av$term, ignore.case = TRUE, perl = TRUE), ]
  if (nrow(row) == 0) return(NA_real_)
  as.numeric(row$p_value[1])
}

# Responses that live on an unbounded scale by construction (logit, additive /
# centred log-ratio, Fisher-z of Jacobs' D). Methods 2.8.1.2 states these "were
# analysed without transformation"; before 2026-08-08 the code did not actually
# enforce that, and offered log1p/sqrt to them anyway.
.UNBOUNDED_RESPONSE_PATTERNS <- c("^logit_", "^lr$", "^lr_", "^clr$", "^clr_",
                                  "^D_z$", "^D_y$")
.is_unbounded_response <- function(nm) {
  if (is.null(nm) || !nzchar(nm)) return(FALSE)
  any(vapply(.UNBOUNDED_RESPONSE_PATTERNS, grepl, logical(1), x = nm))
}

# Try transforms; return list(y_trans, transform_name, sw_p_raw, sw_p_trans,
#                              levene_p_raw, levene_p_trans, normality_flag)
# normality_flag = TRUE when the chosen transform still fails SW or Levene (p < 0.05).
# Transform selection: prefers any transform that makes BOTH tests pass; if none
# does, falls back to the transform with the highest SW p (as before).
#
# FIXED 2026-08-08 -- two coupled defects:
#   (1) candidates were built as log1p(pmax(y, 0)) / sqrt(pmax(y, 0)). The
#       pmax() SILENTLY FLOORED every negative observation to zero rather than
#       excluding an out-of-domain transform. On signed responses this collapsed
#       many distinct observations onto an identical value:
#         jacobs_low_tp   (D_z)       21 of 48 obs (43.8%) -> 0
#         jacobs_high_tp  (D_z)        6 of 48 obs (12.5%) -> 0
#         zone_sec_medium_timepoint   (lr_medium)  6 of 47 -> 0
#         zone_main_timepoint         (logit_flow) 1 of 48 -> 0
#       Several of those values were <= -1, where log1p is not even defined.
#   (2) the transform search was applied to unbounded log-ratio / logit / Fisher-z
#       responses, which Methods 2.8.1.2 says are analysed untransformed.
# Both are fixed by: skipping the search entirely for unbounded responses, and
# offering a candidate only when the data lie inside its domain (never mutating
# the data to fit the transform).
check_and_transform <- function(y, fixed_str, data, response_name = NULL) {
  y_valid <- y[is.finite(y)]
  if (length(y_valid) < 8) return(list(y_trans=y, transform_name="none",
    sw_p_raw=NA, sw_p_trans=NA,
    levene_p_raw=NA, levene_p_trans=NA, normality_flag=FALSE))

  .fit_resid <- function(y_v, d) {
    d$.y <- y_v
    tryCatch({
      mod <- lm(as.formula(paste(".y ~", fixed_str)), data = d)
      r <- residuals(mod)[is.finite(residuals(mod))]
      if (length(r) > 5000) r <- sample(r, 5000)
      r
    }, error = function(e) y_v - mean(y_v, na.rm = TRUE))
  }

  .sw_p <- function(r) {
    tryCatch(shapiro.test(r)$p.value, error = function(e) NA_real_)
  }

  data_fit <- data[is.finite(y), ]
  .grp_raw <- tryCatch(
    interaction(data_fit[, intersect(c("treatment","zone"), names(data_fit))], drop = TRUE),
    error = function(e) factor(rep("all", nrow(data_fit)))
  )

  r_raw   <- .fit_resid(y_valid, data_fit)
  sw_raw  <- .sw_p(r_raw)
  lev_raw <- levene_p(y_valid, .grp_raw)

  best_name  <- "none"
  best_sw    <- if (is.na(sw_raw)) 1 else sw_raw
  best_lev   <- lev_raw
  y_trans    <- y

  # Unbounded responses are analysed on their own scale (Methods 2.8.1.2).
  .skip_search <- .is_unbounded_response(response_name)
  if (.skip_search && !is.na(sw_raw) && sw_raw < 0.05)
    ts_msg("    transform search skipped: '", response_name,
           "' is an unbounded (logit / log-ratio / Fisher-z) response")

  if (!.skip_search && !is.na(sw_raw) && sw_raw < 0.05) {
    # Domain guards: offer a transform only where it is defined for EVERY
    # observation. Never coerce data into the domain.
    candidates <- list()
    if (all(y_valid > -1)) candidates$log1p <- log1p(y_valid)
    if (all(y_valid >= 0)) candidates$sqrt  <- sqrt(y_valid)
    if (!length(candidates))
      ts_msg("    no transform is defined over the observed range of '",
             response_name %||% "response", "' (min = ",
             signif(min(y_valid), 4), "); analysed untransformed")
    for (nm in names(candidates)) {
      r_t   <- .fit_resid(candidates[[nm]], data_fit)
      sw_t  <- .sw_p(r_t)
      lev_t <- levene_p(candidates[[nm]], .grp_raw)

      sw_ok_new  <- !is.na(sw_t)  && sw_t  > 0.05
      lev_ok_new <- is.na(lev_t)  || lev_t > 0.05
      sw_ok_cur  <- !is.na(best_sw)  && best_sw  > 0.05
      lev_ok_cur <- is.na(best_lev) || best_lev > 0.05

      both_new <- sw_ok_new && lev_ok_new
      both_cur <- sw_ok_cur && lev_ok_cur

      improve <- both_new && !both_cur ||                         # new achieves both, cur doesn't
                 (!both_new && !both_cur &&                       # neither achieves both:
                  !is.na(sw_t) && sw_t > best_sw)                # pick best SW

      if (improve) {
        best_sw    <- if (is.na(sw_t)) best_sw else sw_t
        best_lev   <- lev_t
        best_name  <- nm
        y_trans[is.finite(y)] <- candidates[[nm]]
      }
    }
  }

  # Flag if chosen transform still fails either test
  sw_fail  <- !is.na(best_sw)  && best_sw  < 0.05
  lev_fail <- !is.na(best_lev) && best_lev < 0.05
  norm_flag <- sw_fail || lev_fail

  list(y_trans = y_trans, transform_name = best_name,
       sw_p_raw = sw_raw, sw_p_trans = best_sw,
       levene_p_raw = lev_raw, levene_p_trans = best_lev,
       normality_flag = norm_flag)
}

# Levene's test: returns p-value
levene_p <- function(y, group) {
  d <- data.frame(y = y, g = as.factor(group))[is.finite(y), ]
  if (nlevels(d$g) < 2) return(NA_real_)
  tryCatch(
    car::leveneTest(y ~ g, data = d)$`Pr(>F)`[1],
    error = function(e) NA_real_
  )
}

# ---- Hedges' g + 95% CI from Gaussian LMM (Nakagawa & Cuthill 2007 rebuild) -
# g = d * J, where J = 1 - 3/(4*df - 1) (Hedges 1981 small-sample correction,
# applied multiplicatively to the point estimate and both CI bounds), and d is
# SIGNED (direction retained -- do not report |d|; sign is part of what the
# effect size communicates) and derived from the contrast's own t-ratio and
# denominator df via noncentral-t inversion (effectsize::t_to_d), NOT a Wald
# normal-approximation CI. This propagates the model's actual df (previously
# n_obs - n_fixef, which ignored the random-effect structure entirely) into
# the interval, and is exact rather than asymptotic.
#
# Two standardizers are reported, explicitly labelled, per Nakagawa & Cuthill
# 2007 §III.5 (a residual-only "d" for clustered/repeated data is a WITHIN-
# cluster effect, not the total-variance effect a reader expects from "d"):
#   g        -- primary. Standardizer implicit in the contrast SE, which for a
#               between-cluster factor (e.g. treatment with (1|trial)) already
#               reflects total (residual + between-trial) variability.
#   g_resid  -- secondary/diagnostic. delta / sigma_residual only, i.e. the
#               within-cluster effect; retained for comparability with what
#               this project reported before the rebuild.
# sigma_total (sqrt of the sum of all VarCorr variance components) is reported
# alongside sigma_resid as an audit figure, though it is not itself the
# standardizer used for g (see above).
#
# NOTE: verified empirically (2026-08-08) that emmeans' `lmer.df` option
# defaults to "kenward-roger" whenever pbkrtest is installed (true in this
# environment), so this contrast's df already matches the study's primary KR
# ANOVA for the same term without needing to force it here.
#
# Returns NA for non-Gaussian families.
.calc_hedges_g <- function(res, label = "") {
  na_out <- list(g = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
                 g_resid = NA_real_, ci_lo_resid = NA_real_, ci_hi_resid = NA_real_,
                 delta = NA_real_, sigma_r = NA_real_, sigma_total = NA_real_,
                 J = NA_real_, method = "not_gaussian")
  if (is.null(res) || is.null(res$model)) return(na_out)
  fam <- res$family_used %||% res$family %||% "gaussian"
  if (fam != "gaussian") return(na_out)
  tryCatch({
    em  <- emmeans::emmeans(res$model, ~ treatment, type = "response")
    ctr <- as.data.frame(emmeans::contrast(em, method = "pairwise", infer = c(TRUE, TRUE)))
    if (nrow(ctr) == 0) return(modifyList(na_out, list(method = "no_contrast")))
    delta    <- ctr$estimate[1]
    se_delta <- ctr$SE[1]
    df_ctr   <- ctr$df[1]
    t_ctr    <- if ("t.ratio" %in% names(ctr)) ctr$t.ratio[1] else delta / se_delta
    # Try model sigma first; fall back to empirical pooled within-group SD
    sigma_r <- tryCatch(lme4::sigma(res$model), error = function(e) NA_real_)
    if (!is.finite(sigma_r) || sigma_r < 1e-10) {
      sigma_r <- tryCatch(sd(residuals(res$model)), error = function(e) NA_real_)
    }
    # Ultimate fallback: pooled within-treatment SD from raw data
    if (!is.finite(sigma_r) || sigma_r < 1e-10) {
      d_fallback <- res$data_fit
      if (!is.null(d_fallback) && "treatment" %in% names(d_fallback) &&
          ".y" %in% names(d_fallback)) {
        sigma_r <- tryCatch({
          pool_var <- d_fallback %>%
            dplyr::group_by(treatment) %>%
            dplyr::summarise(v = var(.y, na.rm = TRUE), n = dplyr::n(), .groups = "drop")
          sqrt(sum((pool_var$n - 1) * pool_var$v, na.rm = TRUE) /
               max(sum(pool_var$n - 1), 1))
        }, error = function(e) NA_real_)
      }
    }
    if (!is.finite(sigma_r) || sigma_r < 1e-10)
      return(modifyList(na_out, list(method = "sigma_zero")))
    sigma_total <- tryCatch({
      vc <- as.data.frame(lme4::VarCorr(res$model))
      sqrt(sum(vc$vcov, na.rm = TRUE))
    }, error = function(e) NA_real_)
    if (!is.finite(sigma_total) || sigma_total < 1e-10) sigma_total <- sigma_r
    if (!is.finite(df_ctr) || !is.finite(t_ctr) || df_ctr <= 1)
      return(modifyList(na_out, list(method = "no_df", sigma_r = round(sigma_r, 4),
                                      sigma_total = round(sigma_total, 4))))

    J <- 1 - 3 / (4 * df_ctr - 1)

    if (!isTRUE(.HAS_EFFECTSIZE))
      return(modifyList(na_out, list(method = "effectsize_missing",
                                      delta = round(delta, 4), sigma_r = round(sigma_r, 4),
                                      sigma_total = round(sigma_total, 4), J = round(J, 4))))
    d_res <- tryCatch(effectsize::t_to_d(t_ctr, df_ctr, ci = 0.95), error = function(e) NULL)
    if (is.null(d_res) || !nrow(d_res))
      return(modifyList(na_out, list(method = "t_to_d_failed")))

    g_resid    <- (delta / sigma_r) * J
    se_g_resid <- (se_delta / sigma_r) * J

    list(g           = round(d_res$d[1]      * J, 3),
         ci_lo       = round(d_res$CI_low[1] * J, 3),
         ci_hi       = round(d_res$CI_high[1]* J, 3),
         g_resid     = round(g_resid, 3),
         ci_lo_resid = round(g_resid - 1.96 * se_g_resid, 3),
         ci_hi_resid = round(g_resid + 1.96 * se_g_resid, 3),
         delta       = round(delta, 4),
         sigma_r     = round(sigma_r, 4),
         sigma_total = round(sigma_total, 4),
         J           = round(J, 4),
         method      = "t_to_d_noncentral")
  }, error = function(e) modifyList(na_out, list(method = paste0("error:", e$message))))
}

# Invert the two response transforms check_and_transform() can choose (Part 3f).
# For a Gaussian LMM fit on a manually pre-transformed .y, emmeans'
# type = "response" is a NO-OP: it only inverts a GLM/GLMM LINK function, and
# a Gaussian family has none, so a "response-scale" emmean/CI for a
# log1p/sqrt-transformed outcome was silently still on the transform scale.
.inv_transform <- function(x, transform_name) {
  switch(transform_name, log1p = expm1(x), sqrt = x^2, x)
}

# Run post-hoc + CLD for a single term spec (character like "treatment * zone")
.run_ph_cld <- function(model, term_spec, out_dir, term_tag, transform_name = "none") {
  tryCatch({
    em_fmla <- as.formula(paste("~", term_spec))
    # For numeric predictors (e.g. continuous timepoint), specify all observed
    # values so emmeans evaluates at each level, not just the mean.
    .mf  <- tryCatch(model.frame(model), error = function(e) NULL)
    .vars <- trimws(unlist(strsplit(term_spec, "\\s*[\\*:]\\s*")))
    .vars <- .vars[nzchar(.vars)]
    .at  <- list()
    if (!is.null(.mf)) {
      for (.v in .vars) {
        if (.v %in% names(.mf) && is.numeric(.mf[[.v]]) && !is.factor(.mf[[.v]])) {
          .at[[.v]] <- sort(unique(.mf[[.v]]))
        }
      }
    }
    em <- if (length(.at) > 0)
      emmeans::emmeans(model, em_fmla, type = "response", at = .at)
    else
      emmeans::emmeans(model, em_fmla, type = "response")
    # JOINT (all-pairs, Tukey across full cell grid). Conservative for interactions:
    # with k cells = K = choose(k,2) contrasts (e.g. 2×4 → 28). Reported for
    # completeness but typically masks within-stratum effects visible in graphs.
    cld_df  <- as.data.frame(multcomp::cld(em, Letters = letters, adjust = "tukey"))
    cld_df$.group <- trimws(cld_df$.group)
    cld_df$transform <- transform_name
    if (transform_name %in% c("log1p", "sqrt") && "emmean" %in% names(cld_df)) {
      cld_df$emmean_response   <- .inv_transform(cld_df$emmean,   transform_name)
      cld_df$lower.CL_response <- .inv_transform(cld_df$lower.CL, transform_name)
      cld_df$upper.CL_response <- .inv_transform(cld_df$upper.CL, transform_name)
    }
    ctr_df  <- as.data.frame(emmeans::contrast(em, method = "pairwise", adjust = "tukey",
                                                infer = c(TRUE, TRUE)))
    # Contrasts are NOT back-transformed: a difference of log1p/sqrt-transformed
    # values has no valid inverse-transform reading (unlike a single mean). The
    # `transform` column labels the scale honestly instead of leaving it implicit.
    ctr_df$transform <- transform_name
    readr::write_csv(cld_df,  file.path(out_dir, paste0("cld_",       term_tag, ".csv")))
    readr::write_csv(ctr_df,  file.path(out_dir, paste0("contrasts_", term_tag, ".csv")))

    # ---- SIMPLE EFFECTS (within-stratum) for interaction terms ----------------
    # When term_spec involves >= 2 factors, also compute Tukey within each
    # stratifier separately. e.g. for treatment*timepoint_f:
    #   ~ treatment | timepoint_f   → at each interval: control vs treat
    #   ~ timepoint_f | treatment   → within each treatment: T1 vs T2, ..., (Tukey-adjusted)
    # This is the cohort-level summary that matches what graphs typically show.
    simple <- list()
    if (length(.vars) >= 2) {
      for (stratifier in .vars) {
        focal_vars <- setdiff(.vars, stratifier)
        focal_spec <- paste(focal_vars, collapse = " * ")
        # Skip only if stratifier is truly continuous (>5 unique observed values).
        # Discrete numerics like timepoint (3 values) ARE valid stratifiers —
        # emmeans uses `at` to evaluate at their observed levels.
        .n_uniq <- if (!is.null(.mf) && stratifier %in% names(.mf))
                     length(unique(.mf[[stratifier]])) else Inf
        if (.n_uniq > 5) next
        em_simple <- tryCatch({
          em_f <- as.formula(paste("~", focal_spec, "|", stratifier))
          if (length(.at) > 0)
            emmeans::emmeans(model, em_f, type = "response", at = .at)
          else
            emmeans::emmeans(model, em_f, type = "response")
        }, error = function(e) NULL)
        if (is.null(em_simple)) next
        ctr_simple <- tryCatch(
          as.data.frame(emmeans::contrast(em_simple, method = "pairwise", adjust = "tukey",
                                           infer = c(TRUE, TRUE))),
          error = function(e) NULL)
        if (is.null(ctr_simple)) next
        ctr_simple$transform <- transform_name
        tag <- paste0("contrasts_", term_tag, "__simple_",
                      paste(focal_vars, collapse = "_"), "_within_", stratifier, ".csv")
        readr::write_csv(ctr_simple, file.path(out_dir, tag))
        slot_name <- paste0(paste(focal_vars, collapse = "_"), "_within_", stratifier)
        simple[[slot_name]] <- list(
          spec      = paste("~", focal_spec, "|", stratifier),
          contrasts = ctr_simple,
          emmeans   = as.data.frame(em_simple)
        )
      }
    }

    list(cld = cld_df, contrasts = ctr_df, emmeans = as.data.frame(em),
         contrasts_simple = if (length(simple) > 0) simple else NULL)
  }, error = function(e) {
    warning("Post-hoc/CLD failed for '", term_spec, "': ", conditionMessage(e))
    NULL
  })
}

# ---- Orthogonal polynomial decomposition of the interval effect --------------
# Companion to the uniform 2-df rule (.ALLOW_CONTINUOUS_TP). Forcing timepoint_f
# everywhere costs one degree of freedom on outcomes whose interval effect is
# purely linear; this recovers that power inside the SAME model, with no second
# fit and no parameterisation choice.
#
# The three tracked intervals are equally spaced (midpoints 15, 55 and 95 min --
# 40 min apart), so the 2-df interval effect decomposes ORTHOGONALLY into:
#   linear    -- a monotone trend across the trial. Numerically what the
#                withdrawn continuous parameterisation was testing, so it
#                returns essentially the power that parameterisation had.
#   quadratic -- deviation from linear, i.e. the mid-trial dip. This is the
#                formal test of the non-monotonic pattern; previously the dip
#                was described but never tested.
# Equal spacing is what makes contr.poly's weights the correct ones -- if the
# tracking schedule ever changes, these weights must be revisited.
#
# Reported EXHAUSTIVELY: this function runs for every interval model, including
# ones whose omnibus is null, and writes a file even when nothing is
# significant. That is deliberate and it is what keeps the decomposition
# selection-free -- a missing poly_contrasts.csv would itself be a silent
# selection. The omnibus interaction in anova.csv remains the PRIMARY test and
# alone governs whether an outcome qualifies for reporting; these contrasts
# characterise the SHAPE of an effect, they do not promote a null omnibus to a
# positive result.
#
# adjust = "none" is correct here and not an oversight: the two contrasts are
# mutually orthogonal, pre-specified, and reported in full rather than selected
# post hoc, so there is no multiplicity to correct for within the decomposition.
.write_poly_contrasts <- function(model, out_dir, label, transform_name = "none") {
  tryCatch({
    .mf <- tryCatch(model.frame(model), error = function(e) NULL)
    if (is.null(.mf)) return(invisible(NULL))
    if (!all(c("timepoint_f", "treatment") %in% names(.mf))) return(invisible(NULL))
    tp <- .mf[["timepoint_f"]]
    if (!is.factor(tp)) return(invisible(NULL))
    n_lev <- nlevels(droplevels(tp))
    if (n_lev < 3) return(invisible(NULL))

    poly_lbl <- c("linear", "quadratic", "cubic", "quartic")[seq_len(n_lev - 1)]

    # Rename emmeans' auto-generated "<factor>_poly" / "<factor>_pairwise"
    # columns to stable names, and map poly levels to readable labels.
    .norm <- function(df, block) {
      pc <- names(df)[grepl("_poly$", names(df))]
      if (length(pc) == 1) names(df)[names(df) == pc] <- "component"
      wc <- names(df)[grepl("_pairwise$", names(df))]
      if (length(wc) == 1) names(df)[names(df) == wc] <- "contrast"
      # The two emmeans calls label the polynomial term differently: the
      # interaction form emits "<factor>_poly", the within-stratum form emits a
      # plain "contrast" column holding "linear"/"quadratic". Normalise both
      # onto `component` so the column means the same thing in every row.
      if (!"component" %in% names(df) && "contrast" %in% names(df)) {
        cv <- tolower(trimws(as.character(df$contrast)))
        if (length(cv) && all(cv %in% c("linear", "quadratic", "cubic", "quartic"))) {
          df$component <- df$contrast
          df$contrast  <- NA_character_
        }
      }
      if ("component" %in% names(df)) {
        ci <- suppressWarnings(as.integer(factor(df$component,
                                                 levels = unique(df$component))))
        lab <- poly_lbl[ci]
        df$component <- ifelse(is.na(lab), as.character(df$component), lab)
      }
      df$block <- block
      df
    }

    out <- list()

    # (a) treatment x interval, decomposed: the 2-df omnibus split into its
    #     treatment x linear and treatment x quadratic components.
    out$inter <- tryCatch({
      em <- emmeans::emmeans(model, ~ timepoint_f * treatment)
      .norm(as.data.frame(emmeans::contrast(
        em, interaction = c("poly", "pairwise"), adjust = "none")),
        "treatment_x_interval")
    }, error = function(e) NULL)

    # (b) interval shape WITHIN each treatment arm -- describes each line in the
    #     figures separately, which is what a reader reads off the panel.
    out$within <- tryCatch({
      em <- emmeans::emmeans(model, ~ timepoint_f | treatment)
      .norm(as.data.frame(emmeans::contrast(em, "poly", adjust = "none")),
            "interval_within_treatment")
    }, error = function(e) NULL)

    out <- Filter(Negate(is.null), out)
    if (length(out) == 0) return(invisible(NULL))

    keep <- c("block", "component", "contrast", "treatment", "estimate",
              "SE", "df", "t.ratio", "z.ratio", "p.value")
    res <- do.call(rbind, lapply(out, function(df) {
      for (k in setdiff(keep, names(df))) df[[k]] <- NA
      df[, keep, drop = FALSE]
    }))
    rownames(res) <- NULL
    res$transform            <- transform_name
    res$equal_spacing        <- TRUE   # midpoints 15/55/95 min
    res$adjust               <- "none (orthogonal, pre-specified, exhaustive)"
    res$primary_test         <- "omnibus interaction in anova.csv"

    readr::write_csv(res, file.path(out_dir, "poly_contrasts.csv"))
    invisible(res)
  }, error = function(e) {
    warning("Polynomial decomposition failed for '", label, "': ",
            conditionMessage(e))
    invisible(NULL)
  })
}

# ---- ANOVA shim: converts any Anova output to legacy data.frame ----------------
# Guarantees columns: term, df, chisq, p_value regardless of whether
# KR/F or Wald chi-sq was used. Preserves df_denom (DenDF) when present so
# F-test captions can be rendered as F(df1, df2). Adds stat_type attribute.
.anova_to_legacy_df <- function(av_raw, stat_type = "Wald-chisq") {
  d <- as.data.frame(av_raw)
  d$term <- rownames(d); rownames(d) <- NULL

  # Unify column names: p_value
  pv_col <- grep("^Pr|^p.value|p_value", names(d), ignore.case = TRUE, value = TRUE)[1]
  if (!is.na(pv_col) && pv_col != "p_value") names(d)[names(d) == pv_col] <- "p_value"

  # Unify chi-sq / F → chisq column
  chisq_col <- grep("^Chisq|^F$|^F value|^LR stat", names(d), ignore.case = TRUE, value = TRUE)[1]
  if (!is.na(chisq_col) && chisq_col != "chisq") names(d)[names(d) == chisq_col] <- "chisq"
  if (!"chisq" %in% names(d)) d$chisq <- NA_real_

  # Unify df (numerator)
  df_col <- grep("^Df$|^NumDF|^num.df", names(d), ignore.case = TRUE, value = TRUE)[1]
  if (!is.na(df_col) && df_col != "df") names(d)[names(d) == df_col] <- "df"
  if (!"df" %in% names(d)) d$df <- NA_integer_

  # Preserve denominator df (KR / Satterthwaite F-tests)
  den_col <- grep("^DenDF|^den.df|^df.den", names(d), ignore.case = TRUE, value = TRUE)[1]
  if (!is.na(den_col) && den_col != "df_denom") names(d)[names(d) == den_col] <- "df_denom"
  if (!"df_denom" %in% names(d)) d$df_denom <- NA_real_

  out <- d[, intersect(c("term","df","df_denom","chisq","p_value"), names(d)), drop = FALSE]
  # Persist stat_type as a real column so it survives readr::write_csv (CSV
  # roundtrip drops R attributes). The attribute is kept too for back-compat.
  out$stat_type <- stat_type
  # Set attribute AFTER column subset (the [ operator drops user attrs).
  attr(out, "stat_type") <- stat_type
  out
}

# ---- KR-aware Type III ANOVA for lmer models ----------------------------------
# Returns list(kr = <legacy df>, sw = <legacy df or NULL>).
# kr: Kenward-Roger F (or Wald chi-sq fallback).
# sw: Satterthwaite F (NULL when KR unavailable or when sw computation fails).
.run_anova_lmm <- function(model, label) {
  tryCatch({
    # Plain lm (OLS fallback path): use car::Anova type-III F directly.
    # Inject residual df from the model into df_denom so .bh_row reports F(df1, df2).
    if (inherits(model, "lm") && !inherits(model, "lmerMod")) {
      a <- car::Anova(model, type = "III")
      out <- .anova_to_legacy_df(a, "F-OLS")
      .res_df <- tryCatch(model$df.residual, error = function(e) NA_real_)
      out$df_denom <- ifelse(is.na(out$df_denom) | is.null(out$df_denom),
                              .res_df, out$df_denom)
      return(list(kr = out, sw = NULL))
    }
    if (.HAS_KR) {
      m_lt <- lmerTest::as_lmerModLmerTest(model)
      a_kr <- anova(m_lt, type = 3, ddf = "Kenward-Roger")
      a_sw <- tryCatch(
        anova(m_lt, type = 3, ddf = "Satterthwaite"),
        error = function(e) NULL)
      list(
        kr = .anova_to_legacy_df(a_kr, "F-KR"),
        sw = if (!is.null(a_sw)) .anova_to_legacy_df(a_sw, "F-SW") else NULL)
    } else {
      a <- car::Anova(model, type = "III")
      list(kr = .anova_to_legacy_df(a, "Wald-chisq"), sw = NULL)
    }
  }, error = function(e) {
    warning("ANOVA (KR/Wald) failed for ", label, ": ", conditionMessage(e),
            " — retrying with Wald chi-sq.")
    tryCatch({
      a <- car::Anova(model, type = "III")
      list(kr = .anova_to_legacy_df(a, "Wald-chisq"), sw = NULL)
    }, error = function(e2) {
      warning("ANOVA fallback also failed for ", label, ": ", conditionMessage(e2))
      list(kr = NULL, sw = NULL)
    })
  })
}


# Attempt to simplify a singular LMM by falling back to simpler RE structures.
# Returns the updated best list entry; if no non-singular simplification exists,
# returns the original (with a logged message).
.simplify_singular <- function(best, d_fit, label) {
  if (!is.null(best$model) && inherits(best$model, "lm") &&
      !inherits(best$model, "lmerMod")) return(best)  # already plain lm
  if (is.null(best$model) || !lme4::isSingular(best$model)) return(best)

  fallback_res <- c("(1 | tank)", "(1 | phys_trial_id)", "(1 | trial_date)")
  fallback_res <- fallback_res[fallback_res != trimws(best$re)]

  for (re_s in fallback_res) {
    grp_var <- sub(".*\\|\\s*", "", gsub("[()\\s]", "", re_s))
    if (!grp_var %in% names(d_fit)) next
    fmla <- as.formula(paste(".y ~", best$fixed_formula, "+", re_s))
    m_s  <- tryCatch(lme4::lmer(fmla, data = d_fit, REML = FALSE),
                     error = function(e) NULL)
    if (is.null(m_s)) next
    if (!lme4::isSingular(m_s)) {
      ts_msg("    Singular → simplified RE to: ", re_s)
      best$model <- m_s
      best$re    <- re_s
      best$name  <- paste0(best$name, "_simplified")
      return(best)
    }
  }
  ts_msg("    Singular — no non-singular simplification found; keeping best-AICc model.")

  # Optional fallback to plain lm() when all RE candidates are singular
  if (isTRUE(.get_global("SIMPLIFY_SINGULAR_TO_LM", FALSE))) {
    fmla_lm <- as.formula(paste(".y ~", best$fixed_formula))
    m_lm <- tryCatch(lm(fmla_lm, data = d_fit), error = function(e) NULL)
    if (!is.null(m_lm)) {
      aic_lmm <- tryCatch(MuMIn::AICc(best$model), error = function(e) Inf)
      aic_lm  <- tryCatch(MuMIn::AICc(m_lm),       error = function(e) Inf)
      if (is.finite(aic_lm) && aic_lm < aic_lmm - 2) {
        ts_msg("    Singular LMM → lm() better by ΔAICc = ",
               round(aic_lmm - aic_lm, 2), "; simplifying to fixed-effects model.")
        best$model            <- m_lm
        best$re               <- "none (simplified from singular LMM)"
        best$name             <- paste0(best$name, "_lm_simplified")
        best$simplified_to_lm <- TRUE
        return(best)
      } else {
        ts_msg("    lm() not better (delta AICc = ",
               round(aic_lmm - aic_lm, 2), "); keeping singular LMM.")
      }
    }
  }
  best
}

# ---- Uniform interval parameterisation (2-df rule) ---------------------------
# Every trial x interval model uses timepoint as a 3-level FACTOR, so the
# interval and treatment:interval terms are 2-df tests in every outcome and are
# directly comparable across outcomes. The continuous (linear-trend, 1-df)
# alternative is withdrawn from the AICc candidate set entirely.
#
# Why the categorical form is the one to force, rather than the continuous one:
# the two parameterisations fail asymmetrically. Where the factor loses on AICc
# it loses only its parameter penalty (timepoint_f costs exactly two extra fixed
# parameters, ~4-5 AICc units; the observed gaps for nnd/iid/hull_area/lr_medium
# were 4.79-5.27, i.e. almost pure penalty, so the extra parameters bought
# essentially no log-likelihood and no signal is lost). Where the continuous
# form loses it loses the signal itself: on the primary endpoint (lr_flow) the
# factor model wins by 18.19 AICc DESPITE paying that penalty, because the
# treatment x interval pattern is strongly non-monotonic (a mid-trial dip). A
# linear trend cannot represent that and would report "no change across the
# trial", which is false. The factor is the general model; the linear model is
# a special case of it and the wrong one for the headline outcomes.
#
# The power the 2-df test costs on genuinely linear outcomes is recovered inside
# the SAME model by the pre-specified orthogonal polynomial decomposition -- see
# .write_poly_contrasts() below -- so nothing is given up by forcing the factor.
#
# Set to TRUE only to reproduce pre-revision output; leaving it FALSE is what
# enforces the 2-df rule in code rather than by convention. Call sites still
# pass cont_fixed_str = "treatment * timepoint" and are deliberately left
# untouched, so the intent stays visible and the switch is the single override.
.ALLOW_CONTINUOUS_TP <- FALSE

# ---- Master LMM analysis runner ---------------------------------------------
# cont_fixed_str: optional alternative fixed-effects formula using timepoint as
#   continuous integer instead of factor. When provided, candidates are fit with
#   both formulas; the globally best AICc wins. AICc table includes fixed_formula
#   column to make the winner transparent. Suppressed unless .ALLOW_CONTINUOUS_TP.
run_lmm_analysis <- function(label, data, response, fixed_str, re_cands,
                              focal_terms = NULL, cont_fixed_str = NULL) {
  if (!isTRUE(.ALLOW_CONTINUOUS_TP)) cont_fixed_str <- NULL  # uniform 2-df rule
  # Scope sum contrasts so Type III main effects are marginal (F1)
  .old_contrasts <- options(contrasts = c("contr.sum", "contr.poly"))
  on.exit(options(.old_contrasts), add = TRUE)

  out_dir <- file.path(STEP5_OUT, label)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  ts_msg("  [", label, "]")

  y_raw   <- data[[response]]
  valid   <- is.finite(y_raw) & !is.na(data$treatment)
  d_fit   <- data[valid, ]
  y_raw_v <- y_raw[valid]

  if (sum(valid) < 10) {
    warning("Too few valid rows (", sum(valid), ") for ", label)
    return(invisible(NULL))
  }

  # ---- 0) Near-constant response guard -----
  .y_var <- var(y_raw_v, na.rm = TRUE)
  if (!is.finite(.y_var) || .y_var < 1e-8) {
    ts_msg("    SKIPPED: response near-constant (var=", signif(.y_var, 3),
           "). LMM not fitted.")
    readr::write_csv(
      data.frame(label = label, response = response,
                 note  = "Response near-constant; LMM skipped."),
      file.path(out_dir, "skipped.csv")
    )
    return(invisible(NULL))
  }

  # ---- 1) Normality + transform -----
  trans <- check_and_transform(y_raw_v, fixed_str, d_fit, response_name = response)
  d_fit$.y <- trans$y_trans[seq_len(nrow(d_fit))]

  readr::write_csv(
    data.frame(response       = response,
               transform      = trans$transform_name,
               sw_p_raw       = round(trans$sw_p_raw,       4),
               sw_p_trans     = round(trans$sw_p_trans,     4),
               levene_p_raw   = round(trans$levene_p_raw,   4),
               levene_p_trans = round(trans$levene_p_trans, 4),
               normality_flag = trans$normality_flag),
    file.path(out_dir, "normality_check.csv")
  )
  .flag_str <- if (isTRUE(trans$normality_flag)) "  *** NORMALITY FLAG ***" else ""
  ts_msg("    transform: ", trans$transform_name,
         "  SW(raw)=",    round(trans$sw_p_raw,       3),
         "  SW(trans)=",  round(trans$sw_p_trans,     3),
         "  Lev(trans)=", round(trans$levene_p_trans, 3),
         .flag_str)

  # ---- 1b) Filter RE candidates (F7/F13/F17) -----
  # Drop tank RE when <5 unique tank levels (estimate unstable).
  # Strip phys_trial_id RE when data is aggregated (1 row per physical trial).
  re_cands_use <- re_cands
  if ("tank" %in% names(d_fit)) {
    n_tanks <- nlevels(factor(d_fit$tank))
    if (n_tanks < 4) {
      re_cands_use <- re_cands_use[!grepl("\\| tank\\b", re_cands_use)]
      ts_msg("    RE filter: dropped tank (only ", n_tanks, " levels)")
    }
  }
  # Strip phys_trial_id when data is aggregated (1 row per physical trial, no timepoints).
  if ("phys_trial_id" %in% names(d_fit)) {
    n_rows_per_phys <- max(table(d_fit$phys_trial_id))
    if (n_rows_per_phys == 1) {
      re_cands_use <- re_cands_use[!grepl("\\| phys_trial_id\\b", re_cands_use)]
      ts_msg("    RE filter: dropped phys_trial_id (aggregated: 1 row per physical trial)")
    }
  }
  if (length(re_cands_use) == 0) {
    ts_msg("    No valid RE candidates remain — fitting intercept-only (no RE)")
    re_cands_use <- c(none = "1")
    re_fallback  <- TRUE
  } else {
    re_fallback <- FALSE
  }

  # ---- 2) Fit LMM candidates by AICc (REML=TRUE for RE comparison, F6/F18) -----
  # When re_fallback=TRUE (no valid RE candidates remain), skip lme4 entirely
  # and fit ordinary OLS — lme4::lmer cannot fit a model without a RE term, so
  # the categorical-RE branch below would silently return zero candidates and
  # the analysis would be aborted. This path covers the "aggregated, 4-tank"
  # case after item-3 RE candidate cleanup.
  if (isTRUE(re_fallback)) {
    .fit_lm <- function(fs, tag) {
      m <- tryCatch(lm(as.formula(paste(".y ~", fs)), data = d_fit),
                    error = function(e) NULL)
      if (is.null(m)) return(NULL)
      list(name = tag, model = m, aicc = MuMIn::AICc(m),
           re = "1", fixed_formula = fs, fit_method = "OLS",
           simplified_to_lm = TRUE)
    }
    cand_results <- list()
    cr <- .fit_lm(fixed_str, "ols_cat"); if (!is.null(cr)) cand_results[[1]] <- cr
    if (!is.null(cont_fixed_str)) {
      cc <- .fit_lm(cont_fixed_str, "ols_cont"); if (!is.null(cc))
        cand_results[[length(cand_results) + 1]] <- cc
    }
  } else {
    # df_method_memo Part 2d: when cont_fixed_str is supplied, this AICc ranking
    # spans TWO different FIXED-effect structures (categorical vs continuous
    # timepoint), not just different random-effect structures -- REML
    # log-likelihoods are only comparable when the fixed effects are held
    # constant (Pinheiro & Bates 2000), so a REML-based AICc comparison across
    # fixed structures is invalid. Fit every candidate here (both
    # parameterizations, across every RE variant) with ML for this selection
    # step; the winner is refit REML before inference further down
    # ("REML refit before inference"), so the reported model is still REML.
    # Categorical timepoint candidates
    cand_results <- lapply(names(re_cands_use), function(rn) {
      fmla <- as.formula(paste(".y ~", fixed_str, "+", re_cands_use[[rn]]))
      m    <- tryCatch(lme4::lmer(fmla, data = d_fit, REML = FALSE), error = function(e) NULL)
      if (is.null(m)) return(NULL)
      list(name = rn, model = m, aicc = MuMIn::AICc(m),
           re = re_cands_use[[rn]], fixed_formula = fixed_str, fit_method = "ML",
           simplified_to_lm = FALSE)
    })
    cand_results <- Filter(Negate(is.null), cand_results)

    # Continuous timepoint candidates (when requested)
    if (!is.null(cont_fixed_str)) {
      cont_cands <- lapply(names(re_cands_use), function(rn) {
        fmla <- as.formula(paste(".y ~", cont_fixed_str, "+", re_cands_use[[rn]]))
        m    <- tryCatch(lme4::lmer(fmla, data = d_fit, REML = FALSE), error = function(e) NULL)
        if (is.null(m)) return(NULL)
        list(name = paste0(rn, "_cont"), model = m, aicc = MuMIn::AICc(m),
             re = re_cands_use[[rn]], fixed_formula = cont_fixed_str, fit_method = "ML",
             simplified_to_lm = FALSE)
      })
      cand_results <- c(cand_results, Filter(Negate(is.null), cont_cands))
    }
  }

  if (length(cand_results) == 0) {
    warning("No LMM converged for ", label); return(invisible(NULL))
  }

  # ---- Singularity guard (2026-08-18) ---------------------------------------
  # A singular fit means a variance component has collapsed to zero: that random
  # effect is not identifiable, so the candidate is not a model of the design and
  # must not win on AICc alone. Same principle already stated for phys_trial_id
  # above ("not an empirical question AICc should be allowed to answer") -- a
  # zero-variance trial term silently reintroduces the very pseudoreplication
  # that policy exists to prevent.
  #
  # Found via lr_medium (zone_sec_medium_timepoint): AICc picked
  # (1|phys_trial_id) + (1|tank) at 180.15 over (1|phys_trial_id) at 181.46, but
  # the winner was singular with phys_trial_id variance exactly 0. That dropped
  # the school term entirely, treating 48 sessions as independent for a
  # BETWEEN-school factor, and moved the treatment test from F(1,14) = 5.28,
  # p = 0.038 to F(1,39) = 8.93, p = 0.005.
  #
  # Singular candidates are excluded only when a non-singular one exists; if all
  # are singular the full set is kept and .simplify_singular() handles it below,
  # so no analysis loses its model. isSingular() is wrapped so non-lme4 fits
  # (OLS, glmmTMB) simply never test singular and behaviour there is unchanged.
  # Excluded candidates stay in aicc_selection.csv flagged `singular`, so the
  # selection remains auditable.
  .is_sing <- vapply(cand_results, function(cc)
    isTRUE(tryCatch(lme4::isSingular(cc$model), error = function(e) FALSE)), logical(1))
  .guard_on <- any(.is_sing) && !all(.is_sing)
  if (.guard_on)
    ts_msg("    Singularity guard: excluded ", sum(.is_sing), " singular candidate(s) (",
           paste(sapply(cand_results[.is_sing], `[[`, "name"), collapse = ", "), ")")
  .eligible <- if (.guard_on) !.is_sing else rep(TRUE, length(cand_results))

  aicc_vals <- sapply(cand_results, `[[`, "aicc")
  best_idx  <- which(.eligible)[which.min(aicc_vals[.eligible])]
  best      <- cand_results[[best_idx]]

  aicc_tbl <- data.frame(
    re            = sapply(cand_results, `[[`, "name"),
    re_formula    = sapply(cand_results, `[[`, "re"),
    fixed_formula = sapply(cand_results, `[[`, "fixed_formula"),
    fit_method    = sapply(cand_results, function(x) x$fit_method %||% "REML"),
    AICc          = round(aicc_vals, 2),
    delta_AICc    = round(aicc_vals - min(aicc_vals), 2),
    singular      = .is_sing,
    eligible      = .eligible,
    selected      = seq_along(cand_results) == best_idx
  )
  readr::write_csv(aicc_tbl, file.path(out_dir, "aicc_selection.csv"))
  ts_msg("    best RE: ", best$name,
         "  formula: ", if (grepl("_cont$", best$name)) "continuous-tp" else "categorical-tp",
         "  AICc=", round(best$aicc, 1))

  # ---- 2b) Singular simplification -----
  best <- .simplify_singular(best, d_fit, label)

  # ---- 3) Diagnostic plots -----
  r_diag <- residuals(best$model)
  f_diag <- fitted(best$model)
  .pngfile <- file.path(out_dir, "diagnostics.png")
  tryCatch({
    png(.pngfile, width = 1400, height = 600, res = 120)
    par(mfrow = c(1, 2))
    plot(f_diag, r_diag,
         xlab = "Fitted", ylab = "Residuals",
         main = paste(label, "— Residuals vs Fitted"),
         pch = 16, col = rgb(0, 0, 0, 0.4))
    abline(h = 0, lty = 2, col = "grey50")
    qqnorm(r_diag, main = "QQ — residuals", pch = 16, col = rgb(0, 0, 0, 0.4))
    qqline(r_diag, col = "steelblue", lwd = 2)
    dev.off()
  }, error = function(e) {
    if (dev.cur() != 1) dev.off()
    warning("Diagnostic plot failed: ", conditionMessage(e))
  })

  # ---- 3b) DHARMa simulation-based diagnostics (item #9) -----
  if (.HAS_DHARMA) {
    tryCatch({
      sim <- DHARMa::simulateResiduals(best$model, n = 500, plot = FALSE)
      png(file.path(out_dir, "dharma_diagnostics.png"), width = 1600, height = 700, res = 120)
      DHARMa::plot(sim, main = label)
      dev.off()
    }, error = function(e) {
      if (dev.cur() != 1) dev.off()
    })
  }

  # ---- 4) REML refit before inference (skip for plain lm) -----
  if (inherits(best$model, "lmerMod")) {
    best$model <- tryCatch(
      update(best$model, REML = TRUE),
      error = function(e) {
        warning("REML refit failed for ", label, ": ", conditionMessage(e),
                " — using ML fit for inference.")
        best$model
      }
    )
  }

  # ---- 5) Type III ANOVA (KR F if pbkrtest+lmerTest available, else Wald) -----
  .av_both <- .run_anova_lmm(best$model, label)
  av    <- .av_both$kr
  av_sw <- .av_both$sw
  if (!is.null(av)) {
    ts_msg("    ANOVA method: ", attr(av, "stat_type"),
           if (!is.null(av_sw)) " (+ Satterthwaite computed)" else "")
    readr::write_csv(av, file.path(out_dir, "anova.csv"))
    if (!is.null(av_sw))
      readr::write_csv(av_sw, file.path(out_dir, "anova_satterthwaite.csv"))
  }

  # Orthogonal linear/quadratic split of the interval effect. Runs for every
  # interval model, significant or not -- see .write_poly_contrasts().
  .write_poly_contrasts(best$model, out_dir, label,
                        transform_name = trans$transform_name)

  # ---- 6) ANOVA caption strings -----
  anova_caps <- list()
  if (!is.null(focal_terms) && !is.null(av)) {
    for (nm in names(focal_terms)) {
      anova_caps[[nm]] <- anova_caption_str(av, focal_terms[[nm]])
      ts_msg("    ", nm, ": ", anova_caps[[nm]])
    }
  }

  # ---- 7) Post-hoc + CLD: always on full interaction; also on significant terms -----
  posthoc <- list()
  if (!is.null(av)) {
    all_terms <- av$term[av$term != "(Intercept)" & !is.na(av$term)]
    # Full interaction = term(s) with the most colons (highest order)
    n_colons   <- nchar(all_terms) - nchar(gsub(":", "", all_terms))
    max_order  <- if (length(n_colons) > 0) max(n_colons) else 0
    full_inter <- all_terms[n_colons == max_order]
    # Also include any additionally significant terms
    sig_terms  <- av$term[!is.na(av$p_value) & av$p_value < 0.05 &
                           av$term != "(Intercept)"]
    run_terms  <- unique(c(full_inter, sig_terms))
    for (trm in run_terms) {
      vars      <- trimws(strsplit(trm, ":")[[1]])
      term_spec <- paste(vars, collapse = " * ")
      tag       <- gsub(":", "x", trm)
      ph        <- .run_ph_cld(best$model, term_spec, out_dir, tag,
                                transform_name = trans$transform_name)
      if (!is.null(ph)) posthoc[[trm]] <- ph
      ts_msg("    CLD for: ", trm)
    }
  }

  # ---- 8) R²m / R²c (F14) + ICC (Nakagawa & Schielzeth 2013 / Nakagawa,
  # Johnson & Schielzeth 2017 via 'performance') -----
  .r2 <- if (.HAS_PERF) {
    tryCatch(performance::r2(best$model), error = function(e) NULL)
  } else NULL
  if (!is.null(.r2)) {
    r2_val <- list(R2m = as.numeric(.r2$R2_marginal),
                   R2c = as.numeric(.r2$R2_conditional))
    ts_msg("    R²m=", round(r2_val$R2m, 3), "  R²c=", round(r2_val$R2c, 3))
  } else {
    r2_val <- NULL
  }
  .icc_val <- if (.HAS_PERF && inherits(best$model, "merMod")) {
    tryCatch(as.numeric(performance::icc(best$model)$ICC_adjusted), error = function(e) NA_real_)
  } else NA_real_
  if (!is.na(.icc_val)) ts_msg("    ICC(adjusted)=", round(.icc_val, 3))

  # ---- 9) Levene on residuals (F11/F3.2) -----
  .resid_lev_p <- tryCatch({
    .res_v <- residuals(best$model)
    .all_fe <- av$term[av$term != "(Intercept)" & !is.na(av$term)]
    .grp_vars <- unique(unlist(strsplit(.all_fe, ":")))
    .grp_vars <- intersect(.grp_vars, names(d_fit))
    .grp_fac  <- if (length(.grp_vars) > 0)
      interaction(d_fit[, .grp_vars, drop = FALSE], drop = TRUE)
    else
      factor(rep("all", length(.res_v)))
    levene_p(.res_v, .grp_fac)
  }, error = function(e) NA_real_)
  if (!is.na(.resid_lev_p))
    ts_msg("    Levene(residuals)=", round(.resid_lev_p, 3))

  invisible(list(
    model           = best$model,
    fixed_formula   = best$fixed_formula,
    transform       = trans$transform_name,
    normality_flag  = trans$normality_flag,
    re_fallback     = re_fallback,
    simplified_to_lm = isTRUE(best$simplified_to_lm),
    family_used     = "gaussian",
    r2              = r2_val,
    icc             = .icc_val,
    levene_resid_p  = .resid_lev_p,
    stat_type       = attr(av, "stat_type") %||% "Wald-chisq",
    aicc_table      = aicc_tbl,
    anova           = av,
    anova_sw        = av_sw,
    anova_caps      = anova_caps,
    posthoc         = posthoc,
    data_fit        = d_fit
  ))
}

# ---- GLMM analysis runner (for count outcomes: zone switches) ---------------
# Uses Poisson family; falls back to negative binomial if overdispersion ratio > 1.5.
# When force_family = "negbin", always uses glmmTMB::nbinom2 regardless of the
# overdispersion check — use this when the data type mandates NB a priori.
# REML does not apply to GLMMs; AICc selection proceeds on ML fits.
run_glmm_analysis <- function(label, data, response, fixed_str, re_cands,
                               focal_terms = NULL, cont_fixed_str = NULL,
                               offset_var = NULL, force_family = NULL) {
  if (!isTRUE(.ALLOW_CONTINUOUS_TP)) cont_fixed_str <- NULL  # uniform 2-df rule
  # Scope sum contrasts so Type III main effects are marginal (F1)
  .old_contrasts <- options(contrasts = c("contr.sum", "contr.poly"))
  on.exit(options(.old_contrasts), add = TRUE)

  out_dir <- file.path(STEP5_OUT, label)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  ts_msg("  [", label, "] (GLMM)")

  valid   <- is.finite(data[[response]]) & !is.na(data$treatment)
  d_fit   <- data[valid, ]
  y_raw_v <- data[[response]][valid]

  if (sum(valid) < 10) {
    warning("Too few valid rows (", sum(valid), ") for ", label)
    return(invisible(NULL))
  }

  d_fit$.y <- as.integer(round(pmax(y_raw_v, 0)))

  # Offset: log of exposure variable if supplied (P2.1)
  if (!is.null(offset_var) && offset_var %in% names(d_fit)) {
    d_fit$.offset <- log(pmax(d_fit[[offset_var]], 1))
  } else {
    d_fit$.offset <- NULL
  }
  .has_offset <- !is.null(d_fit$.offset)

  # Check overdispersion via simple Poisson GLM
  .od_ratio <- tryCatch({
    m_simp <- if (.has_offset)
      glm(.y ~ treatment + offset(.offset), data = d_fit, family = poisson)
    else
      glm(.y ~ treatment, data = d_fit, family = poisson)
    sum(residuals(m_simp, "pearson")^2) / m_simp$df.residual
  }, error = function(e) 1.0)

  if (!is.null(force_family) && force_family == "negbin") {
    overdisp    <- TRUE
    family_used <- "negbin"
    ts_msg("    family forced: negbin (glmmTMB::nbinom2)",
           "  (Pearson dispersion ratio = ", round(.od_ratio, 2), ")",
           if (.has_offset) "  [offset applied]" else "")
  } else {
    overdisp    <- isTRUE(.od_ratio > 1.5)
    family_used <- if (overdisp) "negbin" else "poisson"
    ts_msg("    family: ", family_used,
           "  (Pearson dispersion ratio = ", round(.od_ratio, 2), ")",
           if (.has_offset) "  [offset applied]" else "")
  }
  .use_tmb_nb <- isTRUE(!is.null(force_family) && force_family == "negbin")

  # Build fixed-formula set (with offset term if applicable)
  .fix_with_offset <- function(fs) {
    if (.has_offset) paste(fs, "+ offset(.offset)") else fs
  }
  fixed_strs <- list(cat = fixed_str)
  if (!is.null(cont_fixed_str)) fixed_strs[["cont"]] <- cont_fixed_str

  cand_results <- list()
  for (fs_nm in names(fixed_strs)) {
    fs   <- fixed_strs[[fs_nm]]
    foff <- .fix_with_offset(fs)
    for (rn in names(re_cands)) {
      fmla <- as.formula(paste(".y ~", foff, "+", re_cands[[rn]]))
      m <- if (.use_tmb_nb) {
        # Forced NB: use glmmTMB::nbinom2 (more stable than glmer.nb at small N)
        tryCatch(
          glmmTMB::glmmTMB(fmla, data = d_fit, family = glmmTMB::nbinom2()),
          warning = function(w) tryCatch(
            suppressWarnings(glmmTMB::glmmTMB(fmla, data = d_fit,
                                               family = glmmTMB::nbinom2())),
            error = function(e) NULL),
          error = function(e) NULL)
      } else {
        tryCatch({
          if (overdisp) lme4::glmer.nb(fmla, data = d_fit)
          else          lme4::glmer(fmla, data = d_fit, family = poisson)
        }, warning = function(w) {
          tryCatch({
            if (overdisp) suppressWarnings(lme4::glmer.nb(fmla, data = d_fit))
            else          suppressWarnings(lme4::glmer(fmla, data = d_fit, family = poisson))
          }, error = function(e) NULL)
        }, error = function(e) NULL)
      }
      if (is.null(m)) next
      aicc_val <- tryCatch(MuMIn::AICc(m), error = function(e) Inf)
      if (!is.finite(aicc_val)) next
      cand_results[[length(cand_results) + 1]] <- list(
        name          = if (fs_nm == "cat") rn else paste0(rn, "_cont"),
        model         = m,
        aicc          = aicc_val,
        re            = re_cands[[rn]],
        fixed_formula = fs
      )
    }
  }

  if (length(cand_results) == 0) {
    warning("No GLMM converged for ", label); return(invisible(NULL))
  }

  # ---- Singularity guard (2026-08-18) ---------------------------------------
  # A singular fit means a variance component has collapsed to zero: that random
  # effect is not identifiable, so the candidate is not a model of the design and
  # must not win on AICc alone. Same principle already stated for phys_trial_id
  # above ("not an empirical question AICc should be allowed to answer") -- a
  # zero-variance trial term silently reintroduces the very pseudoreplication
  # that policy exists to prevent.
  #
  # Found via lr_medium (zone_sec_medium_timepoint): AICc picked
  # (1|phys_trial_id) + (1|tank) at 180.15 over (1|phys_trial_id) at 181.46, but
  # the winner was singular with phys_trial_id variance exactly 0. That dropped
  # the school term entirely, treating 48 sessions as independent for a
  # BETWEEN-school factor, and moved the treatment test from F(1,14) = 5.28,
  # p = 0.038 to F(1,39) = 8.93, p = 0.005.
  #
  # Singular candidates are excluded only when a non-singular one exists; if all
  # are singular the full set is kept and .simplify_singular() handles it below,
  # so no analysis loses its model. isSingular() is wrapped so non-lme4 fits
  # (OLS, glmmTMB) simply never test singular and behaviour there is unchanged.
  # Excluded candidates stay in aicc_selection.csv flagged `singular`, so the
  # selection remains auditable.
  .is_sing <- vapply(cand_results, function(cc)
    isTRUE(tryCatch(lme4::isSingular(cc$model), error = function(e) FALSE)), logical(1))
  .guard_on <- any(.is_sing) && !all(.is_sing)
  if (.guard_on)
    ts_msg("    Singularity guard: excluded ", sum(.is_sing), " singular candidate(s) (",
           paste(sapply(cand_results[.is_sing], `[[`, "name"), collapse = ", "), ")")
  .eligible <- if (.guard_on) !.is_sing else rep(TRUE, length(cand_results))

  aicc_vals <- sapply(cand_results, `[[`, "aicc")
  best_idx  <- which(.eligible)[which.min(aicc_vals[.eligible])]
  best      <- cand_results[[best_idx]]

  aicc_tbl <- data.frame(
    re            = sapply(cand_results, `[[`, "name"),
    re_formula    = sapply(cand_results, `[[`, "re"),
    fixed_formula = sapply(cand_results, `[[`, "fixed_formula"),
    AICc          = round(aicc_vals, 2),
    delta_AICc    = round(aicc_vals - min(aicc_vals), 2),
    singular      = .is_sing,
    eligible      = .eligible,
    selected      = seq_along(cand_results) == best_idx
  )
  readr::write_csv(aicc_tbl, file.path(out_dir, "aicc_selection.csv"))
  readr::write_csv(
    data.frame(response = response, family = family_used,
               overdispersed = overdisp, dispersion_ratio = round(.od_ratio, 3)),
    file.path(out_dir, "glmm_family.csv")
  )
  ts_msg("    best: ", best$name, "  AICc=", round(best$aicc, 1))

  # Diagnostic plot
  tryCatch({
    r_diag <- residuals(best$model, type = "pearson")
    f_diag <- fitted(best$model)
    png(file.path(out_dir, "diagnostics.png"), width = 1400, height = 600, res = 120)
    par(mfrow = c(1, 2))
    plot(f_diag, r_diag,
         xlab = "Fitted", ylab = "Pearson residuals",
         main = paste(label, "— GLMM Residuals vs Fitted"),
         pch = 16, col = rgb(0, 0, 0, 0.4))
    abline(h = 0, lty = 2, col = "grey50")
    qqnorm(r_diag, main = "QQ — Pearson residuals", pch = 16, col = rgb(0, 0, 0, 0.4))
    qqline(r_diag, col = "steelblue", lwd = 2)
    dev.off()
  }, error = function(e) { if (dev.cur() != 1) dev.off() })

  # DHARMa diagnostics (#9)
  if (.HAS_DHARMA) {
    tryCatch({
      sim <- DHARMa::simulateResiduals(best$model, n = 500, plot = FALSE)
      png(file.path(out_dir, "dharma_diagnostics.png"), width = 1600, height = 700, res = 120)
      DHARMa::plot(sim, main = label)
      dev.off()
    }, error = function(e) { if (dev.cur() != 1) dev.off() })
  }

  # Type III ANOVA (Wald chi-sq; KR not applicable to GLMMs)
  av <- tryCatch({
    .anova_to_legacy_df(car::Anova(best$model, type = "III"), "Wald-chisq")
  }, error = function(e) {
    warning("ANOVA failed for ", label, ": ", conditionMessage(e)); NULL
  })
  if (!is.null(av)) readr::write_csv(av, file.path(out_dir, "anova.csv"))

  # Orthogonal linear/quadratic split of the interval effect (see LMM runner).
  .write_poly_contrasts(best$model, out_dir, label)

  anova_caps <- list()
  if (!is.null(focal_terms) && !is.null(av)) {
    for (nm in names(focal_terms)) {
      anova_caps[[nm]] <- anova_caption_str(av, focal_terms[[nm]])
      ts_msg("    ", nm, ": ", anova_caps[[nm]])
    }
  }

  posthoc <- list()
  if (!is.null(av)) {
    all_terms <- av$term[av$term != "(Intercept)" & !is.na(av$term)]
    n_colons   <- nchar(all_terms) - nchar(gsub(":", "", all_terms))
    max_order  <- if (length(n_colons) > 0) max(n_colons) else 0
    full_inter <- all_terms[n_colons == max_order]
    sig_terms  <- av$term[!is.na(av$p_value) & av$p_value < 0.05 &
                           av$term != "(Intercept)"]
    run_terms  <- unique(c(full_inter, sig_terms))
    for (trm in run_terms) {
      vars      <- trimws(strsplit(trm, ":")[[1]])
      term_spec <- paste(vars, collapse = " * ")
      tag       <- gsub(":", "x", trm)
      ph        <- .run_ph_cld(best$model, term_spec, out_dir, tag)
      if (!is.null(ph)) posthoc[[trm]] <- ph
      ts_msg("    CLD for: ", trm)
    }
  }

  .r2_glmm <- if (.HAS_PERF) tryCatch(performance::r2(best$model), error = function(e) NULL) else NULL
  r2_val_glmm <- if (!is.null(.r2_glmm)) {
    # performance::r2 may use R2_marginal OR R2m depending on version/model type
    .r2m <- .r2_glmm$R2_marginal %||% .r2_glmm$R2m %||% .r2_glmm[[grep("marginal", names(.r2_glmm), ignore.case=TRUE, value=TRUE)[1]]]
    .r2c <- .r2_glmm$R2_conditional %||% .r2_glmm$R2c %||% .r2_glmm[[grep("conditional", names(.r2_glmm), ignore.case=TRUE, value=TRUE)[1]]]
    if (!is.null(.r2m) && length(.r2m) == 1L && is.numeric(.r2m)) {
      ts_msg("    R²m=", round(as.numeric(.r2m), 3),
             "  R²c=", if (!is.null(.r2c) && length(.r2c)==1L) round(as.numeric(.r2c), 3) else "NA")
      list(R2m = as.numeric(.r2m), R2c = if (!is.null(.r2c)) as.numeric(.r2c) else NA_real_)
    } else {
      ts_msg("    R²m: performance::r2 returned no usable marginal R² for this model")
      NULL
    }
  } else NULL

  invisible(list(
    model          = best$model,
    fixed_formula  = best$fixed_formula,
    transform      = paste0("GLMM_", family_used),
    normality_flag = FALSE,
    re_fallback    = FALSE,
    family_used    = family_used,
    family         = family_used,
    r2              = r2_val_glmm,
    levene_resid_p  = NA_real_,
    stat_type       = "Wald-chisq",
    simplified_to_lm = FALSE,
    aicc_table      = aicc_tbl,
    anova           = av,
    anova_sw        = NULL,
    anova_caps      = anova_caps,
    posthoc         = posthoc,
    data_fit        = d_fit
  ))
}


# ---- Parametric-bootstrap LRT for glmmTMB (df_method_memo section 4) -------
# KR/Satterthwaite denominator-df approximations are defined only for LINEAR
# mixed models; neither applies to a GLMM (beta, Poisson, negbin). The memo's
# suggested remedy is pbkrtest::PBmodcomp, but PBmodcomp's methods are
# gls/lm/merMod ONLY -- confirmed empirically 2026-08-08 (no glmmTMB method
# exists) -- so it cannot be applied to the beta-GLMMs used throughout this
# engine. This reimplements the same idea by hand using simulate.glmmTMB
# (present, glmmTMB >= 1.1): simulate response data under the REDUCED model,
# refit both full and reduced models to each simulated dataset, and compare
# the observed LRT to that simulated null distribution.
#
# term_pattern matches the focal term (e.g. "^treatment$"); the reduced model
# drops that term AND any higher-order term containing it (e.g.
# "treatment:zone") via stats::drop.terms(), so marginality/hierarchy is
# preserved for any fixed formula, not just simple ones.
.pb_lrt_term <- function(fixed_str, re_str, data, term_pattern, family,
                          nsim = N_PB_DEFAULT, seed = PB_SEED) {
  na_out <- list(lrt_obs = NA_real_, p_pb = NA_real_, nsim_ok = 0L, method = "na")
  full_form <- tryCatch(as.formula(paste(".y ~", fixed_str, "+", re_str)),
                         error = function(e) NULL)
  if (is.null(full_form)) return(modifyList(na_out, list(method = "bad_formula")))
  tt   <- stats::terms(as.formula(paste(".y ~", fixed_str)))
  labs <- attr(tt, "term.labels")
  if (!length(labs)) return(modifyList(na_out, list(method = "no_terms")))
  drop_idx <- which(vapply(labs, function(l)
    any(grepl(term_pattern, strsplit(l, ":")[[1]], ignore.case = TRUE, perl = TRUE)),
    logical(1)))
  if (!length(drop_idx)) return(modifyList(na_out, list(method = "term_not_found")))
  reduced_tt <- tryCatch(stats::drop.terms(tt, drop_idx, keep.response = FALSE),
                          error = function(e) NULL)
  reduced_labs <- if (is.null(reduced_tt)) character(0) else attr(reduced_tt, "term.labels")
  reduced_fixed <- if (!length(reduced_labs)) "1" else paste(reduced_labs, collapse = " + ")
  reduced_form  <- tryCatch(as.formula(paste(".y ~", reduced_fixed, "+", re_str)),
                             error = function(e) NULL)
  if (is.null(reduced_form)) return(modifyList(na_out, list(method = "bad_reduced_formula")))

  .fit <- function(fmla, dat) tryCatch(
    suppressWarnings(suppressMessages(glmmTMB::glmmTMB(
      fmla, data = dat, family = family,
      control = glmmTMB::glmmTMBControl(optCtrl = list(iter.max = 500))))),
    error = function(e) NULL)

  m_full <- .fit(full_form, data)
  m_red  <- .fit(reduced_form, data)
  if (is.null(m_full) || is.null(m_red))
    return(modifyList(na_out, list(method = "observed_fit_failed")))
  lrt_obs <- tryCatch(max(0, 2 * (as.numeric(logLik(m_full)) - as.numeric(logLik(m_red)))),
                       error = function(e) NA_real_)
  if (!is.finite(lrt_obs)) return(modifyList(na_out, list(method = "loglik_failed")))

  set.seed(seed)
  sim_y <- tryCatch(stats::simulate(m_red, nsim = nsim), error = function(e) NULL)
  if (is.null(sim_y))
    return(modifyList(na_out, list(lrt_obs = round(lrt_obs, 4), method = "simulate_failed")))

  sim_lrt <- vapply(seq_len(nsim), function(i) {
    d_sim <- data
    d_sim$.y <- sim_y[[i]]
    mf <- .fit(full_form, d_sim)
    mr <- .fit(reduced_form, d_sim)
    if (is.null(mf) || is.null(mr)) return(NA_real_)
    v <- tryCatch(2 * (as.numeric(logLik(mf)) - as.numeric(logLik(mr))), error = function(e) NA_real_)
    if (!is.finite(v)) NA_real_ else max(0, v)
  }, numeric(1))

  n_ok <- sum(is.finite(sim_lrt))
  if (n_ok < nsim * 0.5)
    return(modifyList(na_out, list(lrt_obs = round(lrt_obs, 4), nsim_ok = n_ok,
                                    method = "too_many_sim_failures")))
  p_pb <- (1 + sum(sim_lrt[is.finite(sim_lrt)] >= lrt_obs)) / (n_ok + 1)
  list(lrt_obs = round(lrt_obs, 4), p_pb = p_pb, nsim_ok = n_ok, method = "pb_lrt")
}

# =============================================================================
# ==== 4b) BETA GLMM RUNNER (P2.2 / F3) =======================================
# =============================================================================
# For proportion responses bounded (0,1): zone occupancy + prop_active.
# Uses glmmTMB with beta_family(link="logit") + squeeze (y*(n-1)+0.5)/n.
# Falls back to Gaussian LMM when glmmTMB is unavailable (.HAS_TMB == FALSE).
run_betaglmm_analysis <- function(label, data, response, fixed_str, re_cands,
                                   focal_terms = NULL, cont_fixed_str = NULL) {
  if (!isTRUE(.ALLOW_CONTINUOUS_TP)) cont_fixed_str <- NULL  # uniform 2-df rule
  if (!.HAS_TMB) {
    ts_msg("  [", label, "] glmmTMB not available — falling back to Gaussian LMM")
    return(run_lmm_analysis(label = label, data = data, response = response,
                             fixed_str = fixed_str, re_cands = re_cands,
                             focal_terms = focal_terms,
                             cont_fixed_str = cont_fixed_str))
  }

  .old_contrasts <- options(contrasts = c("contr.sum", "contr.poly"))
  on.exit(options(.old_contrasts), add = TRUE)

  out_dir <- file.path(STEP5_OUT, label)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  ts_msg("  [", label, "] (beta GLMM)")

  y_raw   <- data[[response]]
  valid   <- is.finite(y_raw) & !is.na(data$treatment) & y_raw >= 0 & y_raw <= 1
  d_fit   <- data[valid, ]
  y_raw_v <- y_raw[valid]

  if (sum(valid) < 10) {
    warning("Too few valid rows (", sum(valid), ") for ", label)
    return(invisible(NULL))
  }

  # Beta squeeze: map 0/1 exactly to open interval
  n_obs   <- nrow(d_fit)
  d_fit$.y <- (y_raw_v * (n_obs - 1) + 0.5) / n_obs
  d_fit$.y <- pmin(pmax(d_fit$.y, 1e-6), 1 - 1e-6)

  # RE filter (same as run_lmm_analysis, F7/F13/F17)
  re_cands_use <- re_cands
  if ("tank" %in% names(d_fit)) {
    if (nlevels(factor(d_fit$tank)) < 4) {
      re_cands_use <- re_cands_use[!grepl("\\| tank\\b", re_cands_use)]
      ts_msg("    RE filter: dropped tank (<4 levels)")
    }
  }
  if ("phys_trial_id" %in% names(d_fit)) {
    if (max(table(d_fit$phys_trial_id)) == 1) {
      re_cands_use <- re_cands_use[!grepl("\\| phys_trial_id\\b", re_cands_use)]
      ts_msg("    RE filter: dropped phys_trial_id (aggregated: 1 row per physical trial)")
    }
  }
  if (length(re_cands_use) == 0) { re_cands_use <- c(none = "1"); re_fallback <- TRUE
  } else { re_fallback <- FALSE }

  # Fit candidates — AICc via glmmTMB
  .fit_one <- function(fs, rn, re_str) {
    fmla <- as.formula(paste(".y ~", fs, "+", re_str))
    tryCatch(
      glmmTMB::glmmTMB(fmla, data = d_fit, family = glmmTMB::beta_family(link = "logit"),
                        control = glmmTMB::glmmTMBControl(optCtrl = list(iter.max = 500))),
      error   = function(e) NULL,
      warning = function(w) {
        tryCatch(
          suppressWarnings(glmmTMB::glmmTMB(fmla, data = d_fit,
            family = glmmTMB::beta_family(link = "logit"),
            control = glmmTMB::glmmTMBControl(optCtrl = list(iter.max = 500)))),
          error = function(e2) NULL)
      })
  }

  fixed_strs <- list(cat = fixed_str)
  if (!is.null(cont_fixed_str)) fixed_strs[["cont"]] <- cont_fixed_str

  cand_results <- list()
  for (fs_nm in names(fixed_strs)) {
    fs <- fixed_strs[[fs_nm]]
    for (rn in names(re_cands_use)) {
      m <- .fit_one(fs, rn, re_cands_use[[rn]])
      if (is.null(m)) next
      aicc_val <- tryCatch(MuMIn::AICc(m), error = function(e) Inf)
      if (!is.finite(aicc_val)) next
      cand_results[[length(cand_results) + 1]] <- list(
        name = if (fs_nm == "cat") rn else paste0(rn, "_cont"),
        model = m, aicc = aicc_val,
        re = re_cands_use[[rn]], fixed_formula = fs)
    }
  }

  if (length(cand_results) == 0) {
    warning("No beta GLMM converged for ", label, " — falling back to LMM")
    out_lmm <- run_lmm_analysis(label = label, data = data, response = response,
                                 fixed_str = fixed_str, re_cands = re_cands,
                                 focal_terms = focal_terms,
                                 cont_fixed_str = cont_fixed_str)
    # Item 4 (prioritised corrections): if LMM also returns NULL (e.g. response
    # is near-constant after squeeze), fall back to a plain OLS fit on the
    # squeezed response so the BH/exploratory tables receive a valid p-value
    # instead of NA. The OLS fit is documented as `family = "OLS-fallback"`.
    if (is.null(out_lmm)) {
      ts_msg("    [", label, "] LMM also failed — running OLS fallback on .y")
      out_ols <- tryCatch({
        m <- lm(as.formula(paste(".y ~", fixed_str)), data = d_fit)
        av <- car::Anova(m, type = "III")
        av_legacy <- .anova_to_legacy_df(av, "Wald-F-OLS")
        readr::write_csv(
          data.frame(label = label, response = response,
                     note  = "Beta GLMM and LMM both failed; OLS on squeezed .y."),
          file.path(out_dir, "ols_fallback.csv"))
        list(model = m, anova = av_legacy, label = label,
             family_used = "OLS-fallback", stat_type = "Wald-F-OLS",
             re_fallback = TRUE, normality_flag = NA,
             r2 = list(R2m = NA_real_, R2c = NA_real_),
             data_fit = d_fit)
      }, error = function(e) {
        warning("  OLS fallback also failed for ", label, ": ", conditionMessage(e))
        NULL
      })
      return(out_ols)
    }
    return(out_lmm)
  }

  # ---- Singularity guard (2026-08-18) ---------------------------------------
  # A singular fit means a variance component has collapsed to zero: that random
  # effect is not identifiable, so the candidate is not a model of the design and
  # must not win on AICc alone. Same principle already stated for phys_trial_id
  # above ("not an empirical question AICc should be allowed to answer") -- a
  # zero-variance trial term silently reintroduces the very pseudoreplication
  # that policy exists to prevent.
  #
  # Found via lr_medium (zone_sec_medium_timepoint): AICc picked
  # (1|phys_trial_id) + (1|tank) at 180.15 over (1|phys_trial_id) at 181.46, but
  # the winner was singular with phys_trial_id variance exactly 0. That dropped
  # the school term entirely, treating 48 sessions as independent for a
  # BETWEEN-school factor, and moved the treatment test from F(1,14) = 5.28,
  # p = 0.038 to F(1,39) = 8.93, p = 0.005.
  #
  # Singular candidates are excluded only when a non-singular one exists; if all
  # are singular the full set is kept and .simplify_singular() handles it below,
  # so no analysis loses its model. isSingular() is wrapped so non-lme4 fits
  # (OLS, glmmTMB) simply never test singular and behaviour there is unchanged.
  # Excluded candidates stay in aicc_selection.csv flagged `singular`, so the
  # selection remains auditable.
  .is_sing <- vapply(cand_results, function(cc)
    isTRUE(tryCatch(lme4::isSingular(cc$model), error = function(e) FALSE)), logical(1))
  .guard_on <- any(.is_sing) && !all(.is_sing)
  if (.guard_on)
    ts_msg("    Singularity guard: excluded ", sum(.is_sing), " singular candidate(s) (",
           paste(sapply(cand_results[.is_sing], `[[`, "name"), collapse = ", "), ")")
  .eligible <- if (.guard_on) !.is_sing else rep(TRUE, length(cand_results))

  aicc_vals <- sapply(cand_results, `[[`, "aicc")
  best_idx  <- which(.eligible)[which.min(aicc_vals[.eligible])]
  best      <- cand_results[[best_idx]]

  aicc_tbl <- data.frame(
    re            = sapply(cand_results, `[[`, "name"),
    re_formula    = sapply(cand_results, `[[`, "re"),
    fixed_formula = sapply(cand_results, `[[`, "fixed_formula"),
    fit_method    = "MLE",
    AICc          = round(aicc_vals, 2),
    delta_AICc    = round(aicc_vals - min(aicc_vals), 2),
    singular      = .is_sing,
    eligible      = .eligible,
    selected      = seq_along(cand_results) == best_idx
  )
  readr::write_csv(aicc_tbl, file.path(out_dir, "aicc_selection.csv"))
  ts_msg("    best RE: ", best$name, "  AICc=", round(best$aicc, 1))

  # Diagnostics
  tryCatch({
    r_diag <- residuals(best$model, type = "pearson")
    f_diag <- fitted(best$model)
    png(file.path(out_dir, "diagnostics.png"), width = 1400, height = 600, res = 120)
    par(mfrow = c(1, 2))
    plot(f_diag, r_diag, xlab = "Fitted", ylab = "Pearson residuals",
         main = paste(label, "— beta GLMM"), pch = 16, col = rgb(0,0,0,0.4))
    abline(h = 0, lty = 2, col = "grey50")
    qqnorm(r_diag, main = "QQ — Pearson residuals", pch = 16, col = rgb(0,0,0,0.4))
    qqline(r_diag, col = "steelblue", lwd = 2)
    dev.off()
  }, error = function(e) { if (dev.cur() != 1) dev.off() })

  # DHARMa diagnostics (#9)
  if (.HAS_DHARMA) {
    tryCatch({
      sim <- DHARMa::simulateResiduals(best$model, n = 500, plot = FALSE)
      png(file.path(out_dir, "dharma_diagnostics.png"), width = 1600, height = 700, res = 120)
      DHARMa::plot(sim, main = label)
      dev.off()
    }, error = function(e) { if (dev.cur() != 1) dev.off() })
  }

  # ANOVA (Wald chi-sq; car::Anova supports glmmTMB)
  av <- tryCatch({
    .anova_to_legacy_df(car::Anova(best$model, type = "III"), "Wald-chisq")
  }, error = function(e) {
    warning("ANOVA failed for ", label, ": ", conditionMessage(e)); NULL
  })

  # Parametric-bootstrap LRT for the treatment term (df_method_memo section 4:
  # KR/Satterthwaite are undefined for a beta GLMM). Overwrites the treatment
  # row's chisq/p_value/stat_type in place when the bootstrap succeeds; every
  # OTHER term in `av` keeps its Wald-chisq test -- bootstrapping every term in
  # every one of ~18 beta-GLMM fits would multiply runtime well beyond the
  # already-substantial cost of nsim x 2 refits for the treatment term alone.
  pb_res <- .pb_lrt_term(best$fixed_formula, best$re, d_fit, "^treatment$",
                          family = glmmTMB::beta_family(link = "logit"))
  readr::write_csv(
    data.frame(label = label, term = "treatment", lrt_obs = pb_res$lrt_obs,
               p_pb = pb_res$p_pb, nsim_ok = pb_res$nsim_ok, method = pb_res$method,
               stringsAsFactors = FALSE),
    file.path(out_dir, "pb_lrt_treatment.csv"))
  if (!is.null(av) && identical(pb_res$method, "pb_lrt")) {
    trt_idx <- which(grepl("^treatment$", av$term, ignore.case = TRUE))
    if (length(trt_idx)) {
      av$chisq[trt_idx]     <- pb_res$lrt_obs
      av$p_value[trt_idx]   <- pb_res$p_pb
      av$df_denom[trt_idx]  <- pb_res$nsim_ok   # repurposed: nsim_ok, not a df
      av$stat_type[trt_idx] <- "LRT-PB"
      ts_msg("    [", label, "] treatment: PB-LRT = ", round(pb_res$lrt_obs, 3),
             ", p = ", signif(pb_res$p_pb, 4), " (N=", pb_res$nsim_ok, ")")
    }
  } else {
    ts_msg("    [", label, "] PB-LRT not applied to treatment (", pb_res$method,
           ") -- Wald-chisq row retained")
  }
  if (!is.null(av)) readr::write_csv(av, file.path(out_dir, "anova.csv"))

  # Orthogonal linear/quadratic split of the interval effect (see LMM runner).
  .write_poly_contrasts(best$model, out_dir, label)

  anova_caps <- list()
  if (!is.null(focal_terms) && !is.null(av)) {
    for (nm in names(focal_terms)) {
      anova_caps[[nm]] <- anova_caption_str(av, focal_terms[[nm]])
      ts_msg("    ", nm, ": ", anova_caps[[nm]])
    }
  }

  posthoc <- list()
  if (!is.null(av)) {
    all_terms <- av$term[av$term != "(Intercept)" & !is.na(av$term)]
    n_colons   <- nchar(all_terms) - nchar(gsub(":", "", all_terms))
    max_order  <- if (length(n_colons) > 0) max(n_colons) else 0
    full_inter <- all_terms[n_colons == max_order]
    sig_terms  <- av$term[!is.na(av$p_value) & av$p_value < 0.05 & av$term != "(Intercept)"]
    run_terms  <- unique(c(full_inter, sig_terms))
    for (trm in run_terms) {
      vars      <- trimws(strsplit(trm, ":")[[1]])
      term_spec <- paste(vars, collapse = " * ")
      tag       <- gsub(":", "x", trm)
      ph        <- .run_ph_cld(best$model, term_spec, out_dir, tag)
      if (!is.null(ph)) posthoc[[trm]] <- ph
      ts_msg("    CLD for: ", trm)
    }
  }

  # R²
  .r2 <- if (.HAS_PERF) tryCatch(performance::r2(best$model), error = function(e) NULL) else NULL
  r2_val <- if (!is.null(.r2)) list(R2m = as.numeric(.r2$R2_marginal),
                                     R2c = as.numeric(.r2$R2_conditional)) else NULL

  invisible(list(
    model          = best$model,
    fixed_formula  = best$fixed_formula,
    transform      = "beta_logit",
    normality_flag   = FALSE,
    re_fallback      = re_fallback,
    family_used      = "beta",
    family           = "beta",
    r2               = r2_val,
    levene_resid_p   = NA_real_,
    stat_type        = "Wald-chisq",
    simplified_to_lm = FALSE,
    aicc_table       = aicc_tbl,
    anova            = av,
    anova_sw         = NULL,
    anova_caps       = anova_caps,
    posthoc          = posthoc,
    data_fit         = d_fit
  ))
}


# =============================================================================
# ==== 5) RUN SIX ANALYSES ====================================================
# =============================================================================

# ---- easy_scripts helper-load guard -----------------------------------------
# When sourced by easy_scripts/_helpers.R (which sets EASY_SCRIPTS_HELPER_LOAD
# env var), abort the source here BEFORE any model fitting or figure assembly.
# easy_scripts/_helpers.R wraps the sys.source() call in tryCatch() that
# catches the custom condition below and continues silently. Effect on normal
# pipeline runs: zero (env var unset → guard is FALSE → execution continues).
# The easy_scripts bridge separately parses this same file for the .mg_*/
# .make_* function definitions in §6, so mini-scripts get all helpers without
# paying the cost of refitting every model on every load.
if (isTRUE(getOption("easy_scripts_helper_load", FALSE)) ||
    isTRUE(get0(".EASY_HELPER_LOAD", envir = globalenv(), inherits = FALSE,
                ifnotfound = FALSE)) ||
    nzchar(Sys.getenv("EASY_SCRIPTS_HELPER_LOAD"))) {
  if (exists("ts_msg")) ts_msg("[easy_scripts] aborting source at start of §5 (model fitting)")
  # Build a classed condition with .Data = list(message, call), class attribute set.
  .cond_data <- list(
    message = "easy_scripts helpers-only load completed at §5 boundary",
    call    = sys.call()
  )
  class(.cond_data) <- c("easy_scripts_helpers_only", "condition", "error")
  stop(.cond_data)
}

.has_tp <- any(!is.na(df$timepoint))

# =============================================================================
# A1: Main-zone occupancy — COMPOSITIONALLY CORRECT REPLACEMENT
# =============================================================================
# DESIGN CHOICE: The original beta GLMM on the stacked long-format (zone rows
# per session) doubles the observations (N=96 instead of 48) and introduces
# algebraic dependency (prop_flow + prop_calm = 1 within each session). A model
# with treatment × zone on that data has an inflated effective sample size and
# is formally rank-deficient once both zones are modelled jointly.
#
# REPLACEMENT: LMM on logit(prop_flow) from df_main_wide (N=48 sessions) with
# fixed effects treatment * timepoint_f and RE (1|phys_trial_id).
# — logit(prop_flow) is the natural 1-df scalar summary of the 2-zone split.
# — phys_trial_id as RE captures the within-physical-trial correlation across
#   the three repeated timepoint sessions.
# — Treatment × timepoint interaction tests whether the treatment effect changes
#   across the session (pre-/during-/post-exercise structure).
# — No zone factor: the zone contrast IS the response variable after logit.
ts_msg("=== A1: Main-zone logit(prop_flow) x Treatment x Timepoint [LMM, compositional] ===")
res_zone_main <- if (.has_tp && !is.null(df_main_wide) && nrow(df_main_wide) >= 12) {
  run_lmm_analysis(
    label          = "zone_main_timepoint",
    data           = df_main_wide,
    response       = "logit_flow",
    fixed_str      = "treatment * timepoint_f",
    cont_fixed_str = "treatment * timepoint",
    re_cands       = RE_WIDE,
    focal_terms    = list(
      "Treatment"           = "^treatment$",
      "Timepoint"           = "^timepoint",
      "Treatment:Timepoint" = "treatment.*timepoint|timepoint.*treatment"
    )
  )
} else {
  ts_msg("  Skipped: no timepoints or df_main_wide unavailable.")
  NULL
}

# =============================================================================
# A2: Sub-zone occupancy — COMPOSITIONALLY CORRECT REPLACEMENT
# =============================================================================
# DESIGN CHOICE: The original Gaussian LMM on stacked 4-zone long-format has the
# same compositional dependency problem as A1 (sum of all sub-zone props = 1).
# Additionally, the 4-zone LMM is rank-deficient: one zone is a linear function
# of the other three, so the ANOVA produces the numerically degenerate p≈1 /
# chi-sq≈1e-31 result seen in outputs. A formal ILR basis (3 orthonormal log-
# contrasts) is the full solution; log-ratios vs calm are used here as a simpler
# and biologically equivalent approximation when calm is a natural reference.
#
# REPLACEMENT: 3 separate LMMs, one per log-ratio (lr_high, lr_medium, lr_low)
# each modelled as logit_zone ~ treatment * timepoint_f + (1|phys_trial_id).
# These three models are tested independently; no combined BH family is applied
# here (each is exploratory).
ts_msg("=== A2: Sub-zone log-ratios x Treatment x Timepoint [LMM, 3 models] ===")
.run_lr_tp <- function(lr_var, label_suffix) {
  if (is.null(df_sec_wide) || !lr_var %in% names(df_sec_wide) ||
      nrow(df_sec_wide) < 12 || !.has_tp) {
    ts_msg("  ", label_suffix, " skipped.")
    return(NULL)
  }
  run_lmm_analysis(
    label          = paste0("zone_sec_", label_suffix, "_timepoint"),
    data           = df_sec_wide,
    response       = lr_var,
    fixed_str      = "treatment * timepoint_f",
    cont_fixed_str = "treatment * timepoint",
    re_cands       = RE_WIDE,
    focal_terms    = list(
      "Treatment"           = "^treatment$",
      "Timepoint"           = "^timepoint",
      "Treatment:Timepoint" = "treatment.*timepoint|timepoint.*treatment"
    )
  )
}
res_zone_sec_high_tp   <- .run_lr_tp("lr_high",   "high")
res_zone_sec_med_tp    <- .run_lr_tp("lr_medium",  "medium")
res_zone_sec_low_tp    <- .run_lr_tp("lr_low",     "low")

# Backward-compat alias for code that expects res_zone_sec (graph code, report).
# Points to the high-zone log-ratio model (primary interest: high current zone use).
res_zone_sec <- res_zone_sec_high_tp

# =============================================================================
# A1_agg / A1c: Aggregated main-zone — beta GLMM retained for graphs + descriptive
# =============================================================================
# A1_agg (beta GLMM on long-format) is kept for graph generation only. It is
# not used for inference (see 2.10); A1c (logit LMM) is the compositionally
# correct model used instead.
# The beta GLMM is technically misspecified (compositional dependency) but produces
# no numerical failures at the aggregated level; retaining it solely for its CLD
# annotation on the zone-proportion bar chart.
ts_msg("=== A1_agg: Main-zone occupancy (aggregated, DESCRIPTIVE ONLY) [beta GLMM] ===")
res_zone_main_agg <- run_betaglmm_analysis(
  label      = "zone_main_aggregated",
  data       = df_main_long_agg,
  response   = "prop_time",
  fixed_str  = "treatment * zone",
  re_cands   = RE_LONG_AGG,
  focal_terms = list(
    "Treatment"      = "^treatment$",
    "Zone"           = "^zone$",
    "Treatment:Zone" = "treatment.*zone|zone.*treatment"
  )
)

# A1c: Compositionally correct aggregated analysis — PRIMARY inference for main zones.
# Uses df_main_wide_agg (N=16 physical trials, one logit(prop_flow) per trial).
ts_msg("=== A1c: logit(prop_flow) x Treatment (compositional, aggregated) ===")
res_zone_flow_logit_agg <- if (!is.null(df_main_wide_agg) && nrow(df_main_wide_agg) >= 8) {
  run_lmm_analysis(
    label       = "zone_flow_logit_aggregated",
    data        = df_main_wide_agg,
    response    = "logit_flow",
    fixed_str   = "treatment + factor(tank) + fish_density_f",
    re_cands    = RE_WIDE_AGG,
    focal_terms = list("Treatment" = "^treatment$")
  )
} else {
  ts_msg("  A1c skipped: df_main_wide_agg unavailable or too few rows.")
  NULL
}

# =============================================================================
# A2_agg: Sub-zone aggregated — 3 log-ratio LMMs (REPLACES 4-zone Gaussian LMM)
# =============================================================================
# Same rationale as A2 TP. Three separate LMMs on df_sec_wide_agg (N=16).
ts_msg("=== A2_agg: Sub-zone log-ratios x Treatment (aggregated) [LMM, 3 models] ===")
.run_lr_agg <- function(lr_var, label_suffix) {
  if (is.null(df_sec_wide_agg) || !lr_var %in% names(df_sec_wide_agg) ||
      nrow(df_sec_wide_agg) < 8) {
    ts_msg("  ", label_suffix, " skipped.")
    return(NULL)
  }
  run_lmm_analysis(
    label       = paste0("zone_sec_", label_suffix, "_aggregated"),
    data        = df_sec_wide_agg,
    response    = lr_var,
    fixed_str   = "treatment + factor(tank) + fish_density_f",
    re_cands    = RE_WIDE_AGG,
    focal_terms = list("Treatment" = "^treatment$")
  )
}
res_zone_sec_high_agg   <- .run_lr_agg("lr_high",   "high")
res_zone_sec_med_agg    <- .run_lr_agg("lr_medium",  "medium")
res_zone_sec_low_agg    <- .run_lr_agg("lr_low",     "low")

# Backward-compat alias for downstream code expecting res_zone_sec_agg.
res_zone_sec_agg <- res_zone_sec_high_agg

# =============================================================================
# A2_pairs: Sub-zone pairwise comparisons via stacked alr LMM (EXPLORATORY)
# =============================================================================
# Tests sub-zone-vs-sub-zone preference (high vs medium, high vs low,
# medium vs low) on the alr scale. Within each (treatment x trial), the
# difference between two log-ratios vs calm equals the log-ratio between
# the two sub-zones (area constants cancel). Single LMM with treatment x
# zone_lr; emmeans Tukey contrasts on zone_lr give the pairwise sub-zone
# comparisons. Exploratory (no BH); Tukey within-model only.
ts_msg("=== A2_pairs: Stacked sub-zone alr LMM (pairwise zone contrasts) ===")

df_sec_lr_long_agg <- if (!is.null(df_sec_wide_agg)) {
  df_sec_wide_agg %>%
    dplyr::select(dplyr::any_of(c("phys_trial_id","tank","trial_date",
                                   "treatment","fish_density_f")),
                  lr_high, lr_medium, lr_low) %>%
    tidyr::pivot_longer(c(lr_high, lr_medium, lr_low),
                        names_to = "zone_lr", values_to = "lr") %>%
    dplyr::mutate(zone_lr = factor(sub("^lr_", "", zone_lr),
                                   levels = c("high","medium","low")))
} else NULL

df_sec_lr_long_tp <- if (!is.null(df_sec_wide)) {
  df_sec_wide %>%
    dplyr::select(dplyr::any_of(c("phys_trial_id","tank","trial_date",
                                   "treatment","timepoint_f",
                                   "fish_density_f")),
                  lr_high, lr_medium, lr_low) %>%
    tidyr::pivot_longer(c(lr_high, lr_medium, lr_low),
                        names_to = "zone_lr", values_to = "lr") %>%
    dplyr::mutate(zone_lr = factor(sub("^lr_", "", zone_lr),
                                   levels = c("high","medium","low")))
} else NULL

# Balance check
if (!is.null(df_sec_lr_long_agg)) {
  .bal_agg <- df_sec_lr_long_agg %>%
    dplyr::count(treatment, zone_lr)
  ts_msg("  A2_pairs agg balance: ",
         paste(sprintf("%s/%s=%d", .bal_agg$treatment, .bal_agg$zone_lr,
                       .bal_agg$n), collapse = "; "))
}

res_zone_sec_lr_pairs_agg <- if (!is.null(df_sec_lr_long_agg) &&
                                 nrow(df_sec_lr_long_agg) >= 24) {
  run_lmm_analysis(
    label       = "zone_sec_lr_pairs_aggregated",
    data        = df_sec_lr_long_agg,
    response    = "lr",
    fixed_str   = "treatment * zone_lr",
    re_cands    = c("phys_trial_id"      = "(1|phys_trial_id)",
                    "phys_trial_id+tank" = "(1|phys_trial_id) + (1|tank)"),
    focal_terms = list(
      "Treatment"      = "^treatment$",
      "Zone"           = "^zone_lr$",
      "Treatment:Zone" = "treatment.*zone_lr|zone_lr.*treatment"
    )
  )
} else {
  ts_msg("  A2_pairs agg skipped: df_sec_lr_long_agg unavailable.")
  NULL
}

res_zone_sec_lr_pairs_tp <- if (!is.null(df_sec_lr_long_tp) &&
                                nrow(df_sec_lr_long_tp) >= 24 && .has_tp) {
  run_lmm_analysis(
    label       = "zone_sec_lr_pairs_timepoint",
    data        = df_sec_lr_long_tp,
    response    = "lr",
    fixed_str   = "treatment * zone_lr * timepoint_f",
    re_cands    = c("phys_trial_id"        = "(1|phys_trial_id)",
                    "phys_trial_zone"      = "(1|phys_trial_id) + (1|zone_lr:phys_trial_id)",
                    "phys_trial_zone+tank" = "(1|phys_trial_id) + (1|zone_lr:phys_trial_id) + (1|tank)"),
    focal_terms = list(
      "Treatment"                = "^treatment$",
      "Zone"                     = "^zone_lr$",
      "Timepoint"                = "^timepoint_f$",
      "Treatment:Zone"           = "treatment.*zone_lr|zone_lr.*treatment",
      "Treatment:Timepoint"      = "treatment.*timepoint|timepoint.*treatment",
      "Zone:Timepoint"           = "zone_lr.*timepoint|timepoint.*zone_lr",
      "Treatment:Zone:Timepoint" = "treatment.*zone_lr.*timepoint|treatment.*timepoint.*zone_lr|zone_lr.*treatment.*timepoint|zone_lr.*timepoint.*treatment|timepoint.*treatment.*zone_lr|timepoint.*zone_lr.*treatment"
    )
  )
} else {
  ts_msg("  A2_pairs tp skipped: df_sec_lr_long_tp unavailable or no TP.")
  NULL
}

# =============================================================================
# CELL-MEAN COMPOSITIONAL MODELS — primary cross-cell inference
# =============================================================================
# Long-stacked models that compare every (treatment x zone) cell to every
# other cell. Two families, both EXPLORATORY (Tukey within model only; raw p,
# no adjustment for multiple comparisons):
#   (i)  MAIN-zone beta-GLMM on raw proportions p (flow + calm = 1 within
#        trial), logit link. Cell means on the proportion scale; 4 cells, 6
#        pairwise contrasts.
#   (ii) SUB-zone Gaussian LMM on the centred log-ratio (CLR) of the 4-part
#        area-normalised composition (high, medium, low, calm). Includes calm
#        as a first-class cell (alr cannot). 8 cells; emmeans Tukey gives the
#        full cell grid. RE structure (1|trial) + (1|zone:trial) absorbs the
#        sum-to-zero compositional constraint.
# =============================================================================
ts_msg("=== CELL-MEAN: zone x treatment cross-comparisons (main + sub) ===")

# --- Main-zone cell-mean long data (4 cells: flow|calm x 2 treatments) -------
df_main_cells_agg <- if (!is.null(df_main_wide_agg)) {
  df_main_wide_agg %>%
    dplyr::select(dplyr::any_of(c("phys_trial_id","tank","trial_date",
                                   "treatment","fish_density_f")),
                  prop_flow, prop_calm) %>%
    tidyr::pivot_longer(c(prop_flow, prop_calm),
                        names_to = "zone", values_to = "p",
                        names_prefix = "prop_") %>%
    dplyr::mutate(zone = factor(zone, levels = c("flow","calm")))
} else NULL

df_main_cells_tp <- if (!is.null(df_main_wide)) {
  df_main_wide %>%
    dplyr::select(dplyr::any_of(c("phys_trial_id","tank","trial_date",
                                   "treatment","timepoint_f","fish_density_f")),
                  prop_flow, prop_calm) %>%
    tidyr::pivot_longer(c(prop_flow, prop_calm),
                        names_to = "zone", values_to = "p",
                        names_prefix = "prop_") %>%
    dplyr::mutate(zone = factor(zone, levels = c("flow","calm")))
} else NULL

# --- Sub-zone CLR long data (8 cells: 4 zones x 2 treatments) ----------------
# CLR is the geometric-mean-centred log of the 4-part area-normalised
# composition. Computed within each composition (per trial for agg; per
# trial x timepoint for tp) so calm enters as a real cell with its own value.
.clr_long <- function(d, group_keys) {
  d %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_keys))) %>%
    dplyr::mutate(clr = log(p_an) - mean(log(p_an), na.rm = TRUE)) %>%
    dplyr::ungroup()
}

df_sec_cells_agg <- if (!is.null(df_sec_wide_agg)) {
  .d_agg <- df_sec_wide_agg %>%
    dplyr::select(dplyr::any_of(c("phys_trial_id","tank","trial_date",
                                   "treatment","fish_density_f")),
                  prop_high_an, prop_medium_an, prop_low_an, prop_calm_an) %>%
    tidyr::pivot_longer(c(prop_high_an, prop_medium_an, prop_low_an,
                          prop_calm_an),
                        names_to = "zone", values_to = "p_an") %>%
    dplyr::mutate(zone = factor(sub("^prop_(.+)_an$", "\\1", zone),
                                levels = c("high","medium","low","calm")))
  .clr_long(.d_agg, "phys_trial_id")
} else NULL

df_sec_cells_tp <- if (!is.null(df_sec_wide)) {
  .d_tp <- df_sec_wide %>%
    dplyr::select(dplyr::any_of(c("phys_trial_id","tank","trial_date",
                                   "treatment","timepoint_f","fish_density_f")),
                  prop_high_an, prop_medium_an, prop_low_an, prop_calm_an) %>%
    tidyr::pivot_longer(c(prop_high_an, prop_medium_an, prop_low_an,
                          prop_calm_an),
                        names_to = "zone", values_to = "p_an") %>%
    dplyr::mutate(zone = factor(sub("^prop_(.+)_an$", "\\1", zone),
                                levels = c("high","medium","low","calm")))
  .clr_long(.d_tp, c("phys_trial_id","timepoint_f"))
} else NULL

# Balance check
if (!is.null(df_sec_cells_agg)) {
  .bal_cells <- df_sec_cells_agg %>% dplyr::count(treatment, zone)
  ts_msg("  Cell-mean agg balance: ",
         paste(sprintf("%s/%s=%d", .bal_cells$treatment, .bal_cells$zone,
                       .bal_cells$n), collapse = "; "))
}

# --- Main-zone cell beta-GLMM (aggregated) -----------------------------------
res_zone_main_cells_agg <- if (!is.null(df_main_cells_agg) &&
                               nrow(df_main_cells_agg) >= 16) {
  run_betaglmm_analysis(
    label       = "zone_main_cells_aggregated",
    data        = df_main_cells_agg,
    response    = "p",
    fixed_str   = "treatment * zone",
    re_cands    = c("phys_trial_id"      = "(1|phys_trial_id)",
                    "phys_trial_id+tank" = "(1|phys_trial_id) + (1|tank)"),
    focal_terms = list(
      "Treatment"      = "^treatment$",
      "Zone"           = "^zone$",
      "Treatment:Zone" = "treatment.*zone|zone.*treatment"
    )
  )
} else {
  ts_msg("  Cell-mean main agg skipped: df_main_cells_agg unavailable.")
  NULL
}

# --- Main-zone cell beta-GLMM (timepoint) ------------------------------------
res_zone_main_cells_tp <- if (!is.null(df_main_cells_tp) &&
                              nrow(df_main_cells_tp) >= 16 && .has_tp) {
  run_betaglmm_analysis(
    label       = "zone_main_cells_timepoint",
    data        = df_main_cells_tp,
    response    = "p",
    fixed_str   = "treatment * zone * timepoint_f",
    re_cands    = c("phys_trial_id"        = "(1|phys_trial_id)",
                    "phys_trial_zone"      = "(1|phys_trial_id) + (1|zone:phys_trial_id)",
                    "phys_trial_zone+tank" = "(1|phys_trial_id) + (1|zone:phys_trial_id) + (1|tank)"),
    focal_terms = list(
      "Treatment"                = "^treatment$",
      "Zone"                     = "^zone$",
      "Timepoint"                = "^timepoint_f$",
      "Treatment:Zone"           = "treatment.*zone|zone.*treatment",
      "Treatment:Timepoint"      = "treatment.*timepoint|timepoint.*treatment",
      "Zone:Timepoint"           = "zone.*timepoint|timepoint.*zone",
      "Treatment:Zone:Timepoint" = "treatment.*zone.*timepoint|treatment.*timepoint.*zone|zone.*treatment.*timepoint|zone.*timepoint.*treatment|timepoint.*treatment.*zone|timepoint.*zone.*treatment"
    )
  )
} else {
  ts_msg("  Cell-mean main TP skipped: df_main_cells_tp unavailable or no TP.")
  NULL
}

# --- Sub-zone cell Dirichlet GLMM (aggregated) --------------------------------
# Compositionally rigorous: brms Dirichlet models the full 4-part area-
# normalised composition (high/medium/low/calm) per trial directly, with a
# trial-level random effect correlated across zones via (1|p|phys_trial_id).
# This eliminates the CLR sum-to-zero residual constraint that the CLR LMM
# absorbs only pragmatically.
# Falls back to CLR LMM when brms is not available.
# To enable: install.packages("brms") and set up rstan/cmdstanr + Stan toolchain.

# Helper: CLD letters from a boolean differ-matrix (absorption algorithm).
.dir_cld <- function(cell_ids, means_vec, differ_mat) {
  n    <- length(cell_ids)
  ord  <- order(means_vec[cell_ids], decreasing = TRUE)
  srt  <- cell_ids[ord]
  pool <- c(letters, paste0(letters, "2"))
  ltr_sets <- setNames(vector("list", n), cell_ids)
  idx  <- 1L
  for (anchor in srt) {
    non_diff <- cell_ids[!differ_mat[anchor, cell_ids]]
    already_share <- all(sapply(non_diff, function(cn)
      cn == anchor || length(intersect(ltr_sets[[anchor]], ltr_sets[[cn]])) > 0))
    if (!already_share) {
      ltr <- pool[idx]; idx <- idx + 1L
      for (cn in non_diff) ltr_sets[[cn]] <- c(ltr_sets[[cn]], ltr)
    }
  }
  for (cn in cell_ids)
    if (length(ltr_sets[[cn]]) == 0) { ltr_sets[[cn]] <- pool[idx]; idx <- idx + 1L }
  setNames(sapply(cell_ids, function(cn) paste(sort(ltr_sets[[cn]]), collapse = "")),
           cell_ids)
}

.run_zone_sec_dirichlet_agg <- function() {
  if (!.HAS_BRMS) return(NULL)
  if (is.null(df_sec_wide_agg)) return(NULL)

  out_dir <- file.path(STEP5_OUT, "zone_sec_cells_aggregated")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  ts_msg("  [zone_sec_cells_aggregated] (Dirichlet GLMM via brms)")

  # Prepare wide data: area-normalised proportions, normalised to sum = 1
  .an_cols <- c("prop_high_an","prop_medium_an","prop_low_an","prop_calm_an")
  .dir_df  <- df_sec_wide_agg[, c("phys_trial_id","tank","treatment", .an_cols),
                               drop = FALSE]
  .dir_df  <- .dir_df[complete.cases(.dir_df[, .an_cols]), ]
  if (nrow(.dir_df) < 8) { warning("Too few rows for Dirichlet"); return(NULL) }
  .rs <- rowSums(.dir_df[, .an_cols])
  for (.col in .an_cols) .dir_df[[.col]] <- .dir_df[[.col]] / .rs
  for (.col in .an_cols) .dir_df[[.col]] <- pmin(pmax(.dir_df[[.col]], 1e-6), 1-1e-6)
  .rs2 <- rowSums(.dir_df[, .an_cols])
  for (.col in .an_cols) .dir_df[[.col]] <- .dir_df[[.col]] / .rs2
  names(.dir_df)[names(.dir_df) %in% .an_cols] <-
    c("Y_high","Y_medium","Y_low","Y_calm")
  .resp   <- c("Y_high","Y_medium","Y_low","Y_calm")
  .zones  <- c("high","medium","low","calm")
  .trts   <- c("control","exercise choice")

  # Fit Dirichlet GLMM: calm as reference; correlated trial-level RE across zones
  .fmla <- as.formula(
    "cbind(Y_high, Y_medium, Y_low, Y_calm) ~ treatment + (1|p|phys_trial_id)")
  .fit <- tryCatch(
    brms::brm(
      formula = .fmla, data = .dir_df,
      family  = brms::dirichlet(refcat = "Y_calm"),
      chains  = 4L, iter = 4000L, warmup = 2000L,
      seed    = 42L,
      cores   = min(4L, max(1L, parallel::detectCores() - 1L)),
      silent  = 2L, refresh = 0L),
    error = function(e) {
      warning("[zone_sec_cells_aggregated] brms::brm failed: ", e$message); NULL })
  if (is.null(.fit)) return(NULL)

  # Cell-mean posterior predictions (RE = 0, fixed-effects only)
  .new_df <- data.frame(treatment = .trts)
  .pp <- tryCatch(
    brms::posterior_epred(.fit, newdata = .new_df, re_formula = NA, ndraws = 4000L),
    error = function(e) NULL)
  if (is.null(.pp)) return(NULL)
  # .pp: [ndraws, 2 (trt), 4 (zone)]

  # Cell summary (2 treatments × 4 zones = 8 cells)
  .cells <- expand.grid(treatment = .trts, zone = .zones, stringsAsFactors = FALSE)
  .ti <- c(control = 1L, "exercise choice" = 2L)
  .ki <- setNames(seq_along(.zones), .zones)
  .cell_smry <- do.call(rbind, lapply(seq_len(nrow(.cells)), function(r) {
    .dr <- .pp[, .ti[.cells$treatment[r]], .ki[.cells$zone[r]]]
    data.frame(treatment = .cells$treatment[r], zone = .cells$zone[r],
               emmean = mean(.dr), SE = sd(.dr),
               lower.CL = quantile(.dr, 0.025), upper.CL = quantile(.dr, 0.975))
  }))

  # Pairwise contrasts + differ matrix
  .nc    <- nrow(.cells)
  .cids  <- paste0(.cells$treatment, "/", .cells$zone)
  .dmat  <- matrix(FALSE, .nc, .nc, dimnames = list(.cids, .cids))
  .crows <- list()
  for (.i in seq_len(.nc - 1)) {
    for (.j in (.i + 1):.nc) {
      .d  <- .pp[, .ti[.cells$treatment[.i]], .ki[.cells$zone[.i]]] -
             .pp[, .ti[.cells$treatment[.j]], .ki[.cells$zone[.j]]]
      .pp_pos <- mean(.d > 0)
      .pval   <- 2 * min(.pp_pos, 1 - .pp_pos)
      .crows[[length(.crows)+1]] <- data.frame(
        contrast = paste0(.cids[.i], " - ", .cids[.j]),
        estimate = mean(.d), SE = sd(.d),
        lower.CL = quantile(.d, 0.025), upper.CL = quantile(.d, 0.975),
        p_value  = .pval, prob_pos = .pp_pos)
      .dmat[.cids[.i], .cids[.j]] <- .dmat[.cids[.j], .cids[.i]] <- .pval < 0.05
    }
  }
  .contrasts <- do.call(rbind, .crows)
  readr::write_csv(.contrasts, file.path(out_dir, "contrasts_treatmentxzone.csv"))

  # CLD
  .means <- setNames(.cell_smry$emmean, .cids)
  .cld_ltrs <- .dir_cld(.cids, .means, .dmat)
  .cld_df <- cbind(.cell_smry, .group = .cld_ltrs[.cids])
  readr::write_csv(.cld_df, file.path(out_dir, "cld_treatmentxzone.csv"))

  # ANOVA analogue: posterior treatment effect per zone
  .anov <- do.call(rbind, lapply(seq_along(.zones), function(k) {
    .d   <- .pp[, 2L, k] - .pp[, 1L, k]
    .ppp <- mean(.d > 0)
    data.frame(term = paste0("treatment [EC-Ctrl] zone=", .zones[k]),
               df = 1L, df_denom = NA_real_, chisq = NA_real_,
               p_value = 2 * min(.ppp, 1 - .ppp), prob_EC_gt_ctrl = .ppp,
               note = "Bayesian posterior probability (two-sided)")
  }))
  attr(.anov, "stat_type") <- "Bayesian-posterior"
  readr::write_csv(.anov, file.path(out_dir, "anova.csv"))

  saveRDS(.fit, file.path(out_dir, "brms_dirichlet_fit.rds"))
  ts_msg("    [zone_sec_cells_aggregated] Dirichlet GLMM done.")

  invisible(list(model = .fit, anova = .anov, label = "zone_sec_cells_aggregated",
                 family = "dirichlet", stat_type = "Bayesian-posterior",
                 data_fit = .dir_df))
}

res_zone_sec_cells_agg <- if (.HAS_BRMS && !is.null(df_sec_wide_agg)) {
  tryCatch(.run_zone_sec_dirichlet_agg(),
           error = function(e) {
             ts_msg("  Dirichlet failed (", e$message, ") — falling back to CLR LMM")
             NULL
           })
} else NULL

# Fall back to CLR LMM if Dirichlet not available or failed
if (is.null(res_zone_sec_cells_agg)) {
  if (.HAS_BRMS && is.null(df_sec_wide_agg))
    ts_msg("  Cell-mean sub agg: df_sec_wide_agg unavailable — using CLR LMM")
  else if (!.HAS_BRMS)
    ts_msg("  brms not installed — sub-zone cell model uses CLR LMM (Dirichlet fallback)")
  res_zone_sec_cells_agg <- if (!is.null(df_sec_cells_agg) &&
                                nrow(df_sec_cells_agg) >= 32) {
    run_lmm_analysis(
      label       = "zone_sec_cells_aggregated",
      data        = df_sec_cells_agg,
      response    = "clr",
      fixed_str   = "treatment * zone",
      re_cands    = c("phys_trial_id"        = "(1|phys_trial_id)",
                      "phys_trial_zone"      = "(1|phys_trial_id) + (1|zone:phys_trial_id)",
                      "phys_trial_zone+tank" = "(1|phys_trial_id) + (1|zone:phys_trial_id) + (1|tank)"),
      focal_terms = list(
        "Treatment"      = "^treatment$",
        "Zone"           = "^zone$",
        "Treatment:Zone" = "treatment.*zone|zone.*treatment"
      )
    )
  } else {
    ts_msg("  Cell-mean sub agg skipped: df_sec_cells_agg unavailable.")
    NULL
  }
}

# --- Sub-zone cell CLR LMM (timepoint) ---------------------------------------
res_zone_sec_cells_tp <- if (!is.null(df_sec_cells_tp) &&
                             nrow(df_sec_cells_tp) >= 32 && .has_tp) {
  run_lmm_analysis(
    label       = "zone_sec_cells_timepoint",
    data        = df_sec_cells_tp,
    response    = "clr",
    fixed_str   = "treatment * zone * timepoint_f",
    re_cands    = c("phys_trial_id"        = "(1|phys_trial_id)",
                    "phys_trial_zone"      = "(1|phys_trial_id) + (1|zone:phys_trial_id)",
                    "phys_trial_zone+tank" = "(1|phys_trial_id) + (1|zone:phys_trial_id) + (1|tank)"),
    focal_terms = list(
      "Treatment"                = "^treatment$",
      "Zone"                     = "^zone$",
      "Timepoint"                = "^timepoint_f$",
      "Treatment:Zone"           = "treatment.*zone|zone.*treatment",
      "Treatment:Timepoint"      = "treatment.*timepoint|timepoint.*treatment",
      "Zone:Timepoint"           = "zone.*timepoint|timepoint.*zone",
      "Treatment:Zone:Timepoint" = "treatment.*zone.*timepoint|treatment.*timepoint.*zone|zone.*treatment.*timepoint|zone.*timepoint.*treatment|timepoint.*treatment.*zone|timepoint.*zone.*treatment"
    )
  )
} else {
  ts_msg("  Cell-mean sub TP skipped: df_sec_cells_tp unavailable or no TP.")
  NULL
}

# =============================================================================
# JACOBS' PREFERENCE INDEX MODULE — Descriptive + sensitivity inference
# =============================================================================
# Jacobs' D = (r - p) / (r + p - 2*r*p) per (row x zone). Built from the
# pre-existing wide tables (df_main_wide{,_agg}, df_sec_wide{,_agg}) so the
# raw use proportions are used directly (no logit/log-ratio dependence).
#
# Availability:
#   main-zone domain   p_flow = p_calm = 0.5   (equal-area assumption)
#   secondary domain   p = ZONE_AREA_UNITS / sum(ZONE_AREA_UNITS)
#
# Inferential layer: two parallel paths per zone:
#   (a) beta-GLMM on (D+1)/2 with Smithson-Verkuilen squeeze, logit link
#       (`run_betaglmm_analysis`); honours bounded support natively. Reported
#       as the 8 zone x {agg, TP} tests in report section 6.1; raw p, no
#       adjustment for multiple comparisons.
#   (b) Gaussian LMM on atanh(D * (1 - 1/(2N))) (`.jac_lmm`); reported
#       alongside as a sensitivity / robustness check.
#
# Note on the alr/Jacobs relationship: both metrics respond monotonically to
# changes in zone use proportions and so generally agree in sign of the
# treatment effect, but they are not algebraically equivalent — the test
# statistics can differ in magnitude and significance. Within a session the
# four zone proportions sum to one, so within-session D values inherit an
# exact compositional constraint; per-zone modelling preserves inferential
# integrity. The legacy alr/logit tests (A1c, A2 sub-zone) are retained
# DESCRIPTIVELY only (no BH correction) as an interpretive bridge.
# =============================================================================
ts_msg("=== JACOBS: Preference index — descriptive + sensitivity LMMs ===")

.zarea  <- if (exists("ZONE_AREA_UNITS")) ZONE_AREA_UNITS else
             c(high = 10, medium = 18, low = 13, calm = 42)
.p_main <- c(flow = 0.5, calm = 0.5)
.p_sec  <- .zarea / sum(.zarea)

.jacobs_D <- function(r, p) {
  r <- pmin(pmax(r, 0), 1)
  ifelse(is.na(r) | is.na(p) | p == 0, NA_real_,
         (r - p) / (r + p - 2 * r * p))
}

# Fisher-style atanh squeeze: keeps |z| finite at boundaries.
.jacobs_z <- function(D, N) {
  shrink <- 1 - 1 / (2 * pmax(N, 2))
  atanh(pmin(pmax(D * shrink, -1 + 1e-12), 1 - 1e-12))
}

.build_jacobs <- function(df, prop_cols, p_vec, stratum_lbl) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  N <- nrow(df)
  df %>%
    dplyr::select(dplyr::any_of(c("phys_trial_id","tank","trial_date",
                                   "treatment","fish_density_f",
                                   "timepoint","timepoint_f")),
                  dplyr::all_of(prop_cols)) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(prop_cols),
                        names_to = "zone_raw", values_to = "r") %>%
    dplyr::mutate(
      zone    = sub("^prop_(time_in_)?", "", zone_raw),
      p       = p_vec[zone],
      D       = .jacobs_D(r, p),
      D_z     = .jacobs_z(D, N),
      D_y     = (D + 1) / 2,   # rescaled to [0,1] for beta-GLMM (SV squeeze applied inside runner)
      stratum = stratum_lbl
    ) %>%
    dplyr::select(-zone_raw)
}

.jac_main_agg <- .build_jacobs(df_main_wide_agg,
                                c("prop_flow","prop_calm"),
                                .p_main, "main_agg")
.jac_main_tp  <- .build_jacobs(df_main_wide,
                                c("prop_flow","prop_calm"),
                                .p_main, "main_tp")
.jac_sec_agg  <- .build_jacobs(df_sec_wide_agg,
                                c("prop_high","prop_medium","prop_low","prop_calm"),
                                .p_sec, "sec_agg")
.jac_sec_tp   <- .build_jacobs(df_sec_wide,
                                c("prop_high","prop_medium","prop_low","prop_calm"),
                                .p_sec, "sec_tp")

jacobs_long <- dplyr::bind_rows(.jac_main_agg, .jac_main_tp,
                                .jac_sec_agg,  .jac_sec_tp)

# Boundary-saturation diagnostic (% of |D| == 1 per stratum x zone).
.jac_sat <- jacobs_long %>%
  dplyr::group_by(stratum, zone) %>%
  dplyr::summarise(n     = sum(!is.na(D)),
                   n_sat = sum(abs(D) >= 1 - 1e-9, na.rm = TRUE),
                   pct_sat = round(100 * n_sat / pmax(n, 1), 1),
                   .groups = "drop")

# Descriptive summary: mean +/- SE +/- 95% CI per (treatment x zone x stratum).
.boot_ci <- function(x, R = 1000) {
  x <- x[!is.na(x)]
  if (length(x) < 3) return(c(NA_real_, NA_real_))
  bs <- replicate(R, mean(sample(x, length(x), replace = TRUE)))
  unname(quantile(bs, c(0.025, 0.975)))
}
.jac_summary <- jacobs_long %>%
  dplyr::group_by(stratum, zone, treatment) %>%
  dplyr::summarise(
    n        = sum(!is.na(D)),
    p_avail  = mean(p, na.rm = TRUE),
    D_mean   = mean(D, na.rm = TRUE),
    D_sd     = sd(D, na.rm = TRUE),
    D_se     = D_sd / sqrt(pmax(n, 1)),
    D_lo95   = D_mean - 1.96 * D_se,
    D_hi95   = D_mean + 1.96 * D_se,
    boot     = list(.boot_ci(D)),
    .groups  = "drop"
  ) %>%
  dplyr::mutate(D_boot_lo = sapply(boot, `[`, 1),
                D_boot_hi = sapply(boot, `[`, 2)) %>%
  dplyr::select(-boot) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 4)))

.jac_dir <- file.path(STEP5_OUT, "jacobs")
dir.create(.jac_dir, showWarnings = FALSE, recursive = TRUE)
readr::write_csv(jacobs_long,  file.path(.jac_dir, "jacobs_descriptive.csv"))
readr::write_csv(.jac_summary, file.path(.jac_dir, "jacobs_summary.csv"))
readr::write_csv(.jac_sat,     file.path(.jac_dir, "jacobs_saturation.csv"))
ts_msg("  Jacobs descriptive layer written: ", .jac_dir,
       "  | rows=", nrow(jacobs_long))

# -------- Sensitivity inferential layer --------
# Per-zone LMMs on atanh(D). Main-zone uses 'flow' only (flow + calm with
# equal areas carries 1 df). Sub-zone tests on high/medium/low; calm is the
# implicit baseline (descriptive only).
.jac_lmm <- function(jdf, zone_pick, fixed, re_cands, label) {
  if (is.null(jdf) || !nrow(jdf)) return(NULL)
  d <- jdf %>% dplyr::filter(zone == zone_pick, is.finite(D_z))
  if (nrow(d) < 8) {
    ts_msg("  jacobs_", label, " skipped (n = ", nrow(d), ").")
    return(NULL)
  }
  run_lmm_analysis(
    label       = paste0("jacobs_", label),
    data        = d,
    response    = "D_z",
    fixed_str   = fixed,
    re_cands    = re_cands,
    focal_terms = list(
      "Treatment"           = "^treatment$",
      "Treatment:Timepoint" = "treatment.*timepoint|timepoint.*treatment"
    )
  )
}

res_jac_flow_agg <- .jac_lmm(.jac_main_agg, "flow",   "treatment + factor(tank) + fish_density_f", RE_WIDE_AGG, "flow_agg")
res_jac_high_agg <- .jac_lmm(.jac_sec_agg,  "high",   "treatment + factor(tank) + fish_density_f", RE_WIDE_AGG, "high_agg")
res_jac_med_agg  <- .jac_lmm(.jac_sec_agg,  "medium", "treatment + factor(tank) + fish_density_f", RE_WIDE_AGG, "medium_agg")
res_jac_low_agg  <- .jac_lmm(.jac_sec_agg,  "low",    "treatment + factor(tank) + fish_density_f", RE_WIDE_AGG, "low_agg")

res_jac_flow_tp <- if (.has_tp)
  .jac_lmm(.jac_main_tp, "flow",   "treatment * timepoint_f", RE_WIDE, "flow_tp") else NULL
res_jac_high_tp <- if (.has_tp)
  .jac_lmm(.jac_sec_tp,  "high",   "treatment * timepoint_f", RE_WIDE, "high_tp") else NULL
res_jac_med_tp  <- if (.has_tp)
  .jac_lmm(.jac_sec_tp,  "medium", "treatment * timepoint_f", RE_WIDE, "medium_tp") else NULL
res_jac_low_tp  <- if (.has_tp)
  .jac_lmm(.jac_sec_tp,  "low",    "treatment * timepoint_f", RE_WIDE, "low_tp") else NULL

ts_msg("  Jacobs sensitivity LMMs fitted (8 tests, atanh path).")

# -------- Primary inferential layer: beta-GLMM on (D+1)/2 --------
# Bounded-support inference path. Uses the existing `run_betaglmm_analysis`
# helper which applies a Smithson-Verkuilen squeeze internally and runs an
# AICc-based RE-selection loop. 8 tests = 4 zones x {agg, TP}.
.jac_beta <- function(jdf, zone_pick, fixed, re_cands, label) {
  if (is.null(jdf) || !nrow(jdf)) return(NULL)
  d <- jdf %>% dplyr::filter(zone == zone_pick, is.finite(D_y))
  if (nrow(d) < 8) {
    ts_msg("  jacobs_beta_", label, " skipped (n = ", nrow(d), ").")
    return(NULL)
  }
  run_betaglmm_analysis(
    label       = paste0("jacobs_beta_", label),
    data        = d,
    response    = "D_y",
    fixed_str   = fixed,
    re_cands    = re_cands,
    focal_terms = list(
      "Treatment"           = "^treatment$",
      "Treatment:Timepoint" = "treatment.*timepoint|timepoint.*treatment"
    )
  )
}

res_jac_beta_flow_agg <- .jac_beta(.jac_main_agg, "flow",   "treatment + factor(tank) + fish_density_f", RE_WIDE_AGG, "flow_agg")
res_jac_beta_high_agg <- .jac_beta(.jac_sec_agg,  "high",   "treatment + factor(tank) + fish_density_f", RE_WIDE_AGG, "high_agg")
res_jac_beta_med_agg  <- .jac_beta(.jac_sec_agg,  "medium", "treatment + factor(tank) + fish_density_f", RE_WIDE_AGG, "medium_agg")
res_jac_beta_low_agg  <- .jac_beta(.jac_sec_agg,  "low",    "treatment + factor(tank) + fish_density_f", RE_WIDE_AGG, "low_agg")

res_jac_beta_flow_tp <- if (.has_tp)
  .jac_beta(.jac_main_tp, "flow",   "treatment * timepoint_f", RE_WIDE, "flow_tp") else NULL
res_jac_beta_high_tp <- if (.has_tp)
  .jac_beta(.jac_sec_tp,  "high",   "treatment * timepoint_f", RE_WIDE, "high_tp") else NULL
res_jac_beta_med_tp  <- if (.has_tp)
  .jac_beta(.jac_sec_tp,  "medium", "treatment * timepoint_f", RE_WIDE, "medium_tp") else NULL
res_jac_beta_low_tp  <- if (.has_tp)
  .jac_beta(.jac_sec_tp,  "low",    "treatment * timepoint_f", RE_WIDE, "low_tp") else NULL

ts_msg("  Jacobs PRIMARY beta-GLMMs fitted (8 tests, beta path).")

# ---- Joint Jacobs D beta-GLMMs: all zone × treatment cells in one model ------
# The per-zone models above test treatment within each zone separately and
# cannot produce a zone × treatment interaction or joint Tukey contrasts across
# the full zone × treatment grid. These joint models fit treatment * zone as
# crossed fixed effects on the stacked long data (one row per trial × zone),
# with trial-level random intercepts to account for within-trial correlation
# across zones. The treatment:zone interaction CLD (cld_treatmentxzone.csv)
# is the CLD source for Figure S7 Jacobs D panels.
ts_msg("=== Jacobs JOINT beta-GLMMs: treatment × zone (main + sub, aggregated) ===")

RE_JAC_JOINT_MAIN <- c(phys_trial      = "(1 | phys_trial_id)",
                        phys_trial_tank = "(1 | phys_trial_id) + (1 | tank)")
RE_JAC_JOINT_SEC  <- c(phys_trial      = "(1 | phys_trial_id)",
                        trial_zone      = "(1 | phys_trial_id) + (1 | phys_trial_id:zone)",
                        trial_zone_tank = "(1 | phys_trial_id) + (1 | phys_trial_id:zone) + (1 | tank)")

res_jac_beta_main_joint_agg <- tryCatch({
  d <- .jac_main_agg %>%
    dplyr::filter(zone %in% c("flow", "calm"), is.finite(D_y)) %>%
    dplyr::mutate(zone = factor(zone, levels = c("flow", "calm")))
  run_betaglmm_analysis(
    label       = "jacobs_beta_main_joint_agg",
    data        = d,
    response    = "D_y",
    fixed_str   = "treatment * zone",
    re_cands    = RE_JAC_JOINT_MAIN,
    focal_terms = list(
      "Treatment"       = "^treatment$",
      "Zone"            = "^zone$",
      "Treatment × Zone" = "treatment.*zone|zone.*treatment"
    )
  )
}, error = function(e) { warning("Joint main Jacobs beta-GLMM failed: ", e$message); NULL })

res_jac_beta_sec_joint_agg <- tryCatch({
  d <- .jac_sec_agg %>%
    dplyr::filter(zone %in% c("high", "medium", "low", "calm"), is.finite(D_y)) %>%
    dplyr::mutate(zone = factor(zone, levels = c("high", "medium", "low", "calm")))
  run_betaglmm_analysis(
    label       = "jacobs_beta_sec_joint_agg",
    data        = d,
    response    = "D_y",
    fixed_str   = "treatment * zone",
    re_cands    = RE_JAC_JOINT_SEC,
    focal_terms = list(
      "Treatment"       = "^treatment$",
      "Zone"            = "^zone$",
      "Treatment × Zone" = "treatment.*zone|zone.*treatment"
    )
  )
}, error = function(e) { warning("Joint sub Jacobs beta-GLMM failed: ", e$message); NULL })

ts_msg("  Jacobs JOINT beta-GLMMs fitted (main + sub, agg).")

# =============================================================================
# DESCRIPTIVE MODULE: Zone occupancies treated as independent observations
# [desc] — not used for inference; for graphical annotation only.
#
# These two OLS two-way ANOVAs (treatment × zone) ignore:
#   (a) the compositional constraint (zone proportions sum to 1 within trial), and
#   (b) within-trial clustering (multiple zones per trial are non-independent).
# They are retained because they produce compact-letter displays (CLDs) for the
# main-zone and sub-zone bar charts, where the CLD is clearly captioned as
# coming from a descriptive model. Results are saved with the _desc suffix and
# must not be cited as inferential evidence.
# =============================================================================

ts_msg("=== [desc] Main-zone occupancy — independent ANOVA (treatment x zone) ===")
res_zone_main_desc <- tryCatch({
  if (is.null(df_main_long_agg) || nrow(df_main_long_agg) < 4) stop("data unavailable")
  d_desc <- df_main_long_agg %>%
    dplyr::filter(is.finite(prop_time)) %>%
    dplyr::mutate(.y = prop_time,
                  treatment = factor(treatment),
                  zone      = factor(zone))
  m_desc <- lm(.y ~ treatment * zone, data = d_desc)
  av_desc <- car::Anova(m_desc, type = "III")
  # df_denom = residuals row df (lm-internal residual df, same for all rows)
  .res_df_main <- av_desc[["Df"]][rownames(av_desc) == "Residuals"]
  if (length(.res_df_main) == 0) .res_df_main <- m_desc$df.residual
  av_out  <- data.frame(
    term      = rownames(av_desc),
    df        = av_desc[["Df"]],
    df_denom  = ifelse(rownames(av_desc) == "Residuals", NA_real_, .res_df_main),
    chisq     = av_desc[["F value"]],     # F statistic (column name kept for legacy compat)
    p_value   = av_desc[["Pr(>F)"]],
    stat_type = "F-OLS",
    stringsAsFactors = FALSE
  )
  attr(av_out, "stat_type") <- "F-OLS"
  em_desc <- emmeans::emmeans(m_desc, ~ treatment | zone)
  cld_desc <- tryCatch(
    as.data.frame(multcomp::cld(em_desc, adjust = "tukey",
                                Letters = letters, reversed = FALSE)),
    error = function(e) NULL
  )
  if (!is.null(cld_desc)) names(cld_desc)[names(cld_desc) == ".group"] <- ".group"
  dir_desc <- file.path(STEP5_OUT, "zone_main_desc")
  dir.create(dir_desc, recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(av_out, file.path(dir_desc, "anova.csv"))
  if (!is.null(cld_desc)) readr::write_csv(cld_desc, file.path(dir_desc, "cld_treatment_zone.csv"))
  list(model = m_desc, anova = av_out, posthoc = list("treatment|zone" = list(cld = cld_desc)),
       stat_type = "F-OLS", family_used = "OLS-desc", label = "zone_main_desc",
       anova_caps = list(
         "Treatment:Zone" = anova_caption_str(av_out, "treatment.*zone|zone.*treatment", "F-OLS")
       ))
}, error = function(e) {
  warning("[desc] Main-zone independent ANOVA failed: ", e$message); NULL
})

ts_msg("=== [desc] Sub-zone occupancy — independent ANOVA (treatment x zone) ===")
res_zone_sec_desc <- tryCatch({
  if (is.null(df_sec_long_agg) || nrow(df_sec_long_agg) < 4) stop("data unavailable")
  d_desc_s <- df_sec_long_agg %>%
    dplyr::filter(is.finite(prop_time),
                  zone %in% c("high", "medium", "low", "calm")) %>%
    dplyr::mutate(.y = prop_time,
                  treatment = factor(treatment),
                  zone      = factor(zone, levels = c("high", "medium", "low", "calm")))
  m_desc_s <- lm(.y ~ treatment * zone, data = d_desc_s)
  av_desc_s <- car::Anova(m_desc_s, type = "III")
  .res_df_sec <- av_desc_s[["Df"]][rownames(av_desc_s) == "Residuals"]
  if (length(.res_df_sec) == 0) .res_df_sec <- m_desc_s$df.residual
  av_out_s  <- data.frame(
    term      = rownames(av_desc_s),
    df        = av_desc_s[["Df"]],
    df_denom  = ifelse(rownames(av_desc_s) == "Residuals", NA_real_, .res_df_sec),
    chisq     = av_desc_s[["F value"]],     # F statistic (legacy column name)
    p_value   = av_desc_s[["Pr(>F)"]],
    stat_type = "F-OLS",
    stringsAsFactors = FALSE
  )
  attr(av_out_s, "stat_type") <- "F-OLS"
  em_desc_s <- emmeans::emmeans(m_desc_s, ~ treatment | zone)
  cld_desc_s <- tryCatch(
    as.data.frame(multcomp::cld(em_desc_s, adjust = "tukey",
                                Letters = letters, reversed = FALSE)),
    error = function(e) NULL
  )
  dir_desc_s <- file.path(STEP5_OUT, "zone_sec_desc")
  dir.create(dir_desc_s, recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(av_out_s, file.path(dir_desc_s, "anova.csv"))
  if (!is.null(cld_desc_s)) readr::write_csv(cld_desc_s, file.path(dir_desc_s, "cld_treatment_zone.csv"))
  list(model = m_desc_s, anova = av_out_s, posthoc = list("treatment|zone" = list(cld = cld_desc_s)),
       stat_type = "F-OLS", family_used = "OLS-desc", label = "zone_sec_desc",
       anova_caps = list(
         "Treatment:Zone" = anova_caption_str(av_out_s, "treatment.*zone|zone.*treatment", "F-OLS")
       ))
}, error = function(e) {
  warning("[desc] Sub-zone independent ANOVA failed: ", e$message); NULL
})

# Wire descriptive results into the skinny-graph zone plot caption system.
# Existing res_zone_main_agg and res_zone_sec_agg already feed graph CLDs;
# the _desc objects are an independent source stored in their own directories.

# A3/A4: Switches modelled as Gaussian LMM on switches_per_session with F-KR inference.
ts_msg("=== A3: switches_per_session x Treatment x Timepoint [Gaussian LMM, F-KR] ===")
res_switch_tp <- if (.has_tp && "switches_per_session" %in% names(df)) run_lmm_analysis(
  label          = "switches_timepoint",
  data           = df,
  response       = "switches_per_session",
  fixed_str      = "treatment * timepoint_f",
  cont_fixed_str = "treatment * timepoint",
  re_cands       = RE_WIDE,
  focal_terms    = list(
    "Treatment"           = "^treatment$",
    "Timepoint"           = "^timepoint",
    "Treatment:Timepoint" = "treatment.*timepoint|timepoint.*treatment"
  )
) else { ts_msg("  Skipped (no timepoint or column absent)."); NULL }

ts_msg("=== A4: switches_per_session x Treatment (aggregated) [Gaussian LMM, F-KR] ===")
res_switch_agg <- if ("switches_per_session" %in% names(df_agg)) run_lmm_analysis(
  label      = "switches_aggregated",
  data       = df_agg,
  response   = "switches_per_session",
  fixed_str  = "treatment + factor(tank) + fish_density_f",
  re_cands   = RE_WIDE_AGG,
  focal_terms = list("Treatment" = "^treatment$")
) else { ts_msg("  Skipped (switches_per_session absent in df_agg)."); NULL }

ts_msg("=== A5: prop_active (>1 BL/s) x Treatment x Timepoint [beta GLMM] ===")
# P2.2: beta GLMM for bounded proportion (F3)
res_active_tp <- if (.has_tp) run_betaglmm_analysis(
  label           = "active_timepoint",
  data            = df,
  response        = "prop_active",
  fixed_str       = "treatment * timepoint_f",
  cont_fixed_str  = "treatment * timepoint",
  re_cands        = RE_WIDE,
  focal_terms = list(
    "Treatment"           = "^treatment$",
    "Timepoint"           = "^timepoint",
    "Treatment:Timepoint" = "treatment.*timepoint|timepoint.*treatment"
  )
) else { ts_msg("  Skipped (no timepoints)."); NULL }

ts_msg("=== A6: prop_active x Treatment (aggregated) [beta GLMM] ===")
res_active_agg <- run_betaglmm_analysis(
  label      = "active_aggregated",
  data       = df_agg,
  response   = "prop_active",
  fixed_str  = "treatment + factor(tank) + fish_density_f",
  re_cands   = RE_WIDE_AGG,
  focal_terms = list("Treatment" = "^treatment$")
)

ts_msg("=== A_flux_tp: zone_flux_per_session x Treatment x Timepoint ===")
res_flux_tp <- if (.has_tp && "zone_flux_per_session" %in% names(df)) run_lmm_analysis(
  label           = "flux_timepoint",
  data            = df,
  response        = "zone_flux_per_session",
  fixed_str       = "treatment * timepoint_f",
  cont_fixed_str  = "treatment * timepoint",
  re_cands        = RE_WIDE,
  focal_terms = list(
    "Treatment"           = "^treatment$",
    "Timepoint"           = "^timepoint",
    "Treatment:Timepoint" = "treatment.*timepoint|timepoint.*treatment"
  )
) else { ts_msg("  Skipped (no timepoint or column absent)."); NULL }

ts_msg("=== A_flux_agg: zone_flux_per_session x Treatment (aggregated) ===")
res_flux_agg <- if ("zone_flux_per_session" %in% names(df_agg)) run_lmm_analysis(
  label      = "flux_aggregated",
  data       = df_agg,
  response   = "zone_flux_per_session",
  fixed_str  = "treatment + factor(tank) + fish_density_f",
  re_cands   = RE_WIDE_AGG,
  focal_terms = list("Treatment" = "^treatment$")
) else { ts_msg("  Skipped (column absent)."); NULL }

# Backwards-compat aliases for downstream graph/report code that may still
# reference older variable names. Will be removed once STEP4/STEP5 graphs and
# the report builder are migrated.
res_zone_agg   <- res_zone_main_agg
res_zone_tp    <- res_zone_main
res_moving_tp  <- res_active_tp
res_moving_agg <- res_active_agg


# =============================================================================
# ==== 5b) GROUP DYNAMICS ANALYSES ============================================
# =============================================================================
# group_dynamics_summary: one row per trial_id x timepoint (tank / session level).
# RE_WIDE candidates apply — the unit is the trial/tank, not individual fish.

ts_msg("=== LOADING group_dynamics_summary ===")
if (!exists("group_dynamics_summary", envir = .GlobalEnv, inherits = FALSE)) {
  .gd_path <- .find_latest_csv("STEP2b_output", "group_dynamics_summary.csv")
  if (!is.null(.gd_path)) {
    group_dynamics_summary <- readr::read_csv(.gd_path, show_col_types = FALSE)
    assign("group_dynamics_summary", group_dynamics_summary, envir = .GlobalEnv)
    ts_msg("Loaded group_dynamics_summary from: ", .gd_path)
  } else {
    ts_msg("group_dynamics_summary not found. A7-A16 will be skipped.")
    assign("group_dynamics_summary", NULL, envir = .GlobalEnv)
  }
}

df_gd     <- .get_global("group_dynamics_summary")
df_gd_agg <- NULL   # initialised here; set inside block below if GD data is present

if (!is.null(df_gd) && nrow(df_gd) > 0) {
  # Attach metadata from fish_activity_summary — only columns df_gd doesn't already have
  .gd_want    <- c("trial_id", "trial", "tank", "treatment",
                   "fish_density", "fish_density_f", "timepoint",
                   "timepoint_f", "trial_date", "phys_trial_id")
  .gd_missing <- setdiff(.gd_want, names(df_gd))
  if (length(.gd_missing) > 0) {
    .gd_meta <- dplyr::distinct(
      df[, intersect(c("trial_id", "timepoint", .gd_missing), names(df))],
      trial_id, timepoint, .keep_all = TRUE
    )
    df_gd <- dplyr::left_join(df_gd, .gd_meta, by = c("trial_id", "timepoint"))
  }
  df_gd <- df_gd %>%
    dplyr::mutate(
      treatment      = factor(trimws(tolower(as.character(treatment))),
                              levels = TREATMENT_LEVELS_g),
      timepoint_f    = factor(timepoint, levels = TIMEPOINT_LEVELS_g),
      fish_density_f = factor(fish_density, levels = DENSITY_LEVELS_g)
      # NOTE (2026-08-08): a `mean_centroid_spd_cm * PIPELINE_FPS` conversion
      # used to sit here on the premise that STEP2b's `time` column could be
      # frame numbers rather than elapsed seconds. Verified against a raw
      # trajectories.csv (idtracker_STEP1_choice_exp.R:139, time column reads
      # 0.000, 0.040, 0.080, ... i.e. genuine elapsed seconds at 25 fps) and
      # against STEP2b's own dt_s_trial computation
      # (group_dynamics_STEP2b_choice_exp.R:343-378, centroid_spd = displacement
      # / dt_s_trial, dt_s_trial = median(time - lag(time))): STEP2b already
      # divides by real elapsed time, so mean_centroid_spd_cm is already true
      # cm/s. The multiplication here was double-applying the frame-rate
      # factor, inflating school speed 25x. Removed.
    )

  ts_msg("Group dynamics data: ", nrow(df_gd), " rows | ",
         dplyr::n_distinct(df_gd$trial_id), " sessions")

  # GD aggregated: average GD metrics across timepoints per physical trial.
  {
    .gd_resp_cols  <- intersect(
      c("mean_nnd_cm", "mean_iid_cm",
        "mean_hull_area_cm2", "mean_centroid_spd_cm"),
      names(df_gd))
    .gd_group_cols <- intersect(
      c("tank", "trial_date", "treatment", "fish_density", "fish_density_f"),
      names(df_gd))
    df_gd_agg <- df_gd %>%
      dplyr::mutate(phys_trial_id = as.integer(factor(paste(tank, trial_date, treatment)))) %>%
      dplyr::group_by(phys_trial_id, dplyr::across(dplyr::all_of(.gd_group_cols))) %>%
      dplyr::summarise(
        dplyr::across(dplyr::all_of(.gd_resp_cols), ~ mean(.x, na.rm = TRUE)),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        treatment      = factor(treatment,      levels = levels(df_gd$treatment)),
        fish_density_f = factor(fish_density_f, levels = levels(df_gd$fish_density_f))
      )
    ts_msg("df_gd_agg: ", nrow(df_gd_agg), " physical trials")
  }

  ts_msg("=== A7: Mean NND x Treatment x Timepoint ===")
  res_nnd_tp <- if (.has_tp) run_lmm_analysis(
    label           = "nnd_timepoint",
    data            = df_gd,
    response        = "mean_nnd_cm",
    fixed_str       = "treatment * timepoint_f",
    cont_fixed_str  = "treatment * timepoint",
    re_cands        = RE_WIDE,
    focal_terms = list(
      "Treatment"           = "^treatment$",
      "Timepoint"           = "^timepoint",
      "Treatment:Timepoint" = "treatment.*timepoint|timepoint.*treatment"
    )
  ) else { ts_msg("  Skipped (no timepoints)."); NULL }

  ts_msg("=== A8: Mean NND x Treatment (aggregated) ===")
  res_nnd_agg <- run_lmm_analysis(
    label      = "nnd_aggregated",
    data       = df_gd_agg,
    response   = "mean_nnd_cm",
    fixed_str  = "treatment + factor(tank) + fish_density_f",
    re_cands   = RE_WIDE_AGG,
    focal_terms = list("Treatment" = "^treatment$")
  )

  # A9 / A10 (alignment order parameter) REMOVED 2026-08-16. It is the only
  # collective indicator whose computation requires persistent per-fish
  # identity -- heading is a given fish's frame-to-frame displacement, so an
  # identity swap corrupts it -- while NND, IID, hull area and centroid speed
  # are all permutation-invariant functions of the per-frame position set. It
  # is therefore not fitted, not reported and not plotted anywhere.
  # See METHODS_CHANGES.md section 7. Analysis numbering A7-A16 is unchanged so
  # that archived output remains comparable.

  ts_msg("=== A11: Mean IID x Treatment x Timepoint ===")
  res_iid_tp <- if (.has_tp) run_lmm_analysis(
    label           = "iid_timepoint",
    data            = df_gd,
    response        = "mean_iid_cm",
    fixed_str       = "treatment * timepoint_f",
    cont_fixed_str  = "treatment * timepoint",
    re_cands        = RE_WIDE,
    focal_terms = list(
      "Treatment"           = "^treatment$",
      "Timepoint"           = "^timepoint",
      "Treatment:Timepoint" = "treatment.*timepoint|timepoint.*treatment"
    )
  ) else { ts_msg("  Skipped (no timepoints)."); NULL }

  ts_msg("=== A12: Mean IID x Treatment (aggregated) ===")
  res_iid_agg <- run_lmm_analysis(
    label      = "iid_aggregated",
    data       = df_gd_agg,
    response   = "mean_iid_cm",
    fixed_str  = "treatment + factor(tank) + fish_density_f",
    re_cands   = RE_WIDE_AGG,
    focal_terms = list("Treatment" = "^treatment$")
  )

  ts_msg("=== A13: Convex hull area x Treatment x Timepoint ===")
  res_hull_tp <- if (.has_tp) run_lmm_analysis(
    label           = "hull_area_timepoint",
    data            = df_gd,
    response        = "mean_hull_area_cm2",
    fixed_str       = "treatment * timepoint_f",
    cont_fixed_str  = "treatment * timepoint",
    re_cands        = RE_WIDE,
    focal_terms = list(
      "Treatment"           = "^treatment$",
      "Timepoint"           = "^timepoint",
      "Treatment:Timepoint" = "treatment.*timepoint|timepoint.*treatment"
    )
  ) else { ts_msg("  Skipped (no timepoints)."); NULL }

  ts_msg("=== A14: Convex hull area x Treatment (aggregated) ===")
  res_hull_agg <- run_lmm_analysis(
    label      = "hull_area_aggregated",
    data       = df_gd_agg,
    response   = "mean_hull_area_cm2",
    fixed_str  = "treatment + factor(tank) + fish_density_f",
    re_cands   = RE_WIDE_AGG,
    focal_terms = list("Treatment" = "^treatment$")
  )

  ts_msg("=== A15: Centroid speed x Treatment x Timepoint ===")
  res_cspd_tp <- if (.has_tp) run_lmm_analysis(
    label           = "centroid_speed_timepoint",
    data            = df_gd,
    response        = "mean_centroid_spd_cm",
    fixed_str       = "treatment * timepoint_f",
    cont_fixed_str  = "treatment * timepoint",
    re_cands        = RE_WIDE,
    focal_terms = list(
      "Treatment"           = "^treatment$",
      "Timepoint"           = "^timepoint",
      "Treatment:Timepoint" = "treatment.*timepoint|timepoint.*treatment"
    )
  ) else { ts_msg("  Skipped (no timepoints)."); NULL }

  ts_msg("=== A16: Centroid speed x Treatment (aggregated) ===")
  res_cspd_agg <- run_lmm_analysis(
    label      = "centroid_speed_aggregated",
    data       = df_gd_agg,
    response   = "mean_centroid_spd_cm",
    fixed_str  = "treatment + factor(tank) + fish_density_f",
    re_cands   = RE_WIDE_AGG,
    focal_terms = list("Treatment" = "^treatment$")
  )

} else {
  ts_msg("group_dynamics_summary unavailable — A7-A16 skipped.")
  res_nnd_tp  <- res_nnd_agg  <- NULL
  res_iid_tp  <- res_iid_agg  <- res_hull_tp <- res_hull_agg <- NULL
  res_cspd_tp <- res_cspd_agg <- NULL
}


# =============================================================================
# ==== 5b-bis) easy_scripts CSV export (behaviour) ============================
# =============================================================================
# Writes a frozen-snapshot wide CSV consumed by easy_scripts/ mini-scripts.
# One row per trial × timepoint (≈48 rows). Aggregated columns (Jacobs agg,
# CLR cells agg, etc.) repeat the per-physical-trial value across the 3
# timepoints of that trial, so aggregated mini-scripts can recover them by
# group_by(phys_trial_id) %>% summarise(mean) without loss.
# Skipped on helpers-only loads (guard at start of §5 returns first).
local({
  .easy_dir <- file.path(PROJECT_ROOT, "choice R pipeline/easy_scripts")
  dir.create(.easy_dir, showWarnings = FALSE, recursive = TRUE)
  ts_msg("[easy_scripts] building behaviour CSV → ", .easy_dir)

  # --- Base: df (trial × timepoint) -----------------------------------------
  .out <- df %>%
    dplyr::mutate(
      phys_trial_id = as.integer(factor(paste(tank, trial_date, treatment))),
      timepoint_num = as.integer(timepoint)
    )

  # --- Main-zone props + logit_flow (trial × timepoint grain) ----------------
  if (exists("df_main_wide") && !is.null(df_main_wide)) {
    .mw <- df_main_wide %>%
      dplyr::select(dplyr::any_of(c("trial_id","timepoint","prop_flow",
                                     "prop_calm","logit_flow")))
    .out <- dplyr::left_join(.out, .mw, by = c("trial_id","timepoint"))
  }

  # --- Sub-zone props + lr_* (trial × timepoint grain) -----------------------
  if (exists("df_sec_wide") && !is.null(df_sec_wide)) {
    .sw <- df_sec_wide %>%
      dplyr::select(dplyr::any_of(c("trial_id","timepoint",
        "prop_high","prop_medium","prop_low","prop_calm",
        "prop_high_an","prop_medium_an","prop_low_an","prop_calm_an",
        "lr_high","lr_medium","lr_low")))
    # Avoid prop_flow/prop_calm name clash by renaming sub-zone calm
    if ("prop_calm" %in% names(.sw)) {
      .sw <- dplyr::rename(.sw, prop_calm_sec = prop_calm)
    }
    .out <- dplyr::left_join(.out, .sw, by = c("trial_id","timepoint"))
  }

  # --- Group dynamics (FPS-corrected; trial × timepoint grain) ---------------
  if (exists("df_gd") && !is.null(df_gd)) {
    .gd <- df_gd %>%
      dplyr::select(dplyr::any_of(c("trial_id","timepoint",
        "mean_nnd_cm","mean_iid_cm",
        "mean_hull_area_cm2","mean_centroid_spd_cm")))
    .out <- dplyr::left_join(.out, .gd, by = c("trial_id","timepoint"))
  }

  # --- Jacobs D (trial × timepoint grain via main_tp / sec_tp) ---------------
  # .build_jacobs returns phys_trial_id + timepoint (not trial_id), so join
  # on those keys instead.
  .pivot_jac <- function(jac_df, zone_prefix) {
    if (is.null(jac_df) || nrow(jac_df) == 0) return(NULL)
    jac_df %>%
      dplyr::select(dplyr::any_of(c("phys_trial_id","timepoint","zone","D","D_z","D_y"))) %>%
      tidyr::pivot_wider(names_from = zone,
                          values_from = c(D, D_z, D_y),
                          names_glue  = paste0("{.value}_", zone_prefix, "_{zone}"))
  }
  if (exists(".jac_main_tp")) {
    .jt_main <- .pivot_jac(.jac_main_tp, "main")
    if (!is.null(.jt_main)) {
      .out <- dplyr::left_join(.out, .jt_main, by = c("phys_trial_id","timepoint"))
    }
  }
  if (exists(".jac_sec_tp")) {
    .jt_sec <- .pivot_jac(.jac_sec_tp, "sec")
    if (!is.null(.jt_sec)) {
      .out <- dplyr::left_join(.out, .jt_sec, by = c("phys_trial_id","timepoint"))
    }
  }

  # --- Jacobs D aggregated (phys_trial_id grain) → broadcast over timepoints -
  .pivot_jac_agg <- function(jac_df, zone_prefix) {
    if (is.null(jac_df) || nrow(jac_df) == 0) return(NULL)
    jac_df %>%
      dplyr::select(dplyr::any_of(c("phys_trial_id","zone","D","D_z","D_y"))) %>%
      tidyr::pivot_wider(names_from = zone,
                          values_from = c(D, D_z, D_y),
                          names_glue  = paste0("{.value}_", zone_prefix, "_{zone}_agg"))
  }
  if (exists(".jac_main_agg")) {
    .ja_main <- .pivot_jac_agg(.jac_main_agg, "main")
    if (!is.null(.ja_main)) {
      .out <- dplyr::left_join(.out, .ja_main, by = "phys_trial_id")
    }
  }
  if (exists(".jac_sec_agg")) {
    .ja_sec <- .pivot_jac_agg(.jac_sec_agg, "sec")
    if (!is.null(.ja_sec)) {
      .out <- dplyr::left_join(.out, .ja_sec, by = "phys_trial_id")
    }
  }

  # --- Cell-mean CLR (phys_trial_id grain, agg) → broadcast over timepoints --
  if (exists("df_main_cells_agg") && !is.null(df_main_cells_agg)) {
    .mc_agg <- df_main_cells_agg %>%
      dplyr::select(phys_trial_id, zone, p) %>%
      tidyr::pivot_wider(names_from = zone, values_from = p,
                          names_glue = "p_main_{zone}_agg")
    .out <- dplyr::left_join(.out, .mc_agg, by = "phys_trial_id")
  }
  if (exists("df_sec_cells_agg") && !is.null(df_sec_cells_agg)) {
    .sc_agg <- df_sec_cells_agg %>%
      dplyr::select(phys_trial_id, zone, p_an, clr) %>%
      tidyr::pivot_wider(names_from = zone, values_from = c(p_an, clr),
                          names_glue = "{.value}_sec_{zone}_agg")
    .out <- dplyr::left_join(.out, .sc_agg, by = "phys_trial_id")
  }

  # --- Cell-mean CLR (trial × timepoint grain) -------------------------------
  if (exists("df_sec_cells_tp") && !is.null(df_sec_cells_tp) &&
      "timepoint_f" %in% names(df_sec_cells_tp)) {
    .sc_tp <- df_sec_cells_tp %>%
      dplyr::select(phys_trial_id, timepoint_f, zone, p_an, clr) %>%
      tidyr::pivot_wider(names_from = zone, values_from = c(p_an, clr),
                          names_glue = "{.value}_sec_{zone}_tp")
    .out <- dplyr::left_join(.out, .sc_tp,
                              by = c("phys_trial_id","timepoint_f"))
  }
  if (exists("df_main_cells_tp") && !is.null(df_main_cells_tp) &&
      "timepoint_f" %in% names(df_main_cells_tp)) {
    .mc_tp <- df_main_cells_tp %>%
      dplyr::select(phys_trial_id, timepoint_f, zone, p) %>%
      tidyr::pivot_wider(names_from = zone, values_from = p,
                          names_glue = "p_main_{zone}_tp")
    .out <- dplyr::left_join(.out, .mc_tp,
                              by = c("phys_trial_id","timepoint_f"))
  }

  # --- trial_date format: DD.MM.YYYY (consumed by the workbook builder) ----
  if ("trial_date" %in% names(.out)) {
    .td <- tryCatch(as.Date(.out$trial_date), error = function(e) NULL)
    if (!is.null(.td) && !all(is.na(.td)))
      .out$trial_date <- format(.td, "%d.%m.%Y")
  }

  # --- Write -----------------------------------------------------------------
  .csv_path <- file.path(.easy_dir, "easy_scripts_dataset.csv")
  readr::write_csv(.out, .csv_path, na = "")
  ts_msg("[easy_scripts] wrote ", basename(.csv_path), ": ",
         nrow(.out), " rows × ", ncol(.out), " cols")
})


# =============================================================================
# ==== 5c) INDIVIDUAL SOCIAL CONTEXT ANALYSES (REMOVED) =======================
# =============================================================================
# Per-fish NND, centroid distance, and turning rate depend on persistent fish
# identity, which is unreliable in the choice arena. These analyses (formerly
# A17–A22) are dropped under the school-level refactor.

.social_avail <- character(0)

if (FALSE) {

  ts_msg("=== A17: Fish NND x Treatment x Timepoint ===")
  res_fnnd_tp <- if (.has_tp && "mean_nnd_fish_px" %in% .social_avail) run_lmm_analysis(
    label           = "fish_nnd_timepoint",
    data            = df,
    response        = "mean_nnd_fish_px",
    fixed_str       = "treatment * timepoint_f",
    cont_fixed_str  = "treatment * timepoint",
    re_cands        = RE_WIDE,
    focal_terms = list(
      "Treatment"           = "^treatment$",
      "Timepoint"           = "^timepoint",
      "Treatment:Timepoint" = "treatment.*timepoint|timepoint.*treatment"
    )
  ) else { ts_msg("  Skipped."); NULL }

  ts_msg("=== A18: Fish NND x Treatment (aggregated) ===")
  res_fnnd_agg <- if ("mean_nnd_fish_px" %in% .social_avail) run_lmm_analysis(
    label      = "fish_nnd_aggregated",
    data       = df,
    response   = "mean_nnd_fish_px",
    fixed_str  = "treatment",
    re_cands   = RE_WIDE,
    focal_terms = list("Treatment" = "^treatment$")
  ) else NULL

  ts_msg("=== A19: Centroid distance x Treatment x Timepoint ===")
  res_cdist_tp <- if (.has_tp && "mean_centdist_px" %in% .social_avail) run_lmm_analysis(
    label           = "centroid_dist_timepoint",
    data            = df,
    response        = "mean_centdist_px",
    fixed_str       = "treatment * timepoint_f",
    cont_fixed_str  = "treatment * timepoint",
    re_cands        = RE_WIDE,
    focal_terms = list(
      "Treatment"           = "^treatment$",
      "Timepoint"           = "^timepoint",
      "Treatment:Timepoint" = "treatment.*timepoint|timepoint.*treatment"
    )
  ) else { ts_msg("  Skipped."); NULL }

  ts_msg("=== A20: Centroid distance x Treatment (aggregated) ===")
  res_cdist_agg <- if ("mean_centdist_px" %in% .social_avail) run_lmm_analysis(
    label      = "centroid_dist_aggregated",
    data       = df,
    response   = "mean_centdist_px",
    fixed_str  = "treatment",
    re_cands   = RE_WIDE,
    focal_terms = list("Treatment" = "^treatment$")
  ) else NULL

  ts_msg("=== A21: Turning rate x Treatment x Timepoint ===")
  res_turn_tp <- if (.has_tp && "mean_turning_rate" %in% .social_avail) run_lmm_analysis(
    label           = "turning_rate_timepoint",
    data            = df,
    response        = "mean_turning_rate",
    fixed_str       = "treatment * timepoint_f",
    cont_fixed_str  = "treatment * timepoint",
    re_cands        = RE_WIDE,
    focal_terms = list(
      "Treatment"           = "^treatment$",
      "Timepoint"           = "^timepoint",
      "Treatment:Timepoint" = "treatment.*timepoint|timepoint.*treatment"
    )
  ) else { ts_msg("  Skipped."); NULL }

  ts_msg("=== A22: Turning rate x Treatment (aggregated) ===")
  res_turn_agg <- if ("mean_turning_rate" %in% .social_avail) run_lmm_analysis(
    label      = "turning_rate_aggregated",
    data       = df,
    response   = "mean_turning_rate",
    fixed_str  = "treatment",
    re_cands   = RE_WIDE,
    focal_terms = list("Treatment" = "^treatment$")
  ) else NULL

}
res_fnnd_tp  <- res_fnnd_agg  <- NULL
res_cdist_tp <- res_cdist_agg <- NULL
res_turn_tp  <- res_turn_agg  <- NULL


# =============================================================================
# ==== 5d) PSEUDOREPLICATION CHECK ============================================
# =============================================================================
# Under the school-level refactor the unit of analysis is already trial ×
# timepoint, so fish-level pseudoreplication is structurally absent. The check
# now compares trial × timepoint LMM vs tank-aggregated LM (timepoints averaged)
# to confirm repeated-measures within tank does not drive treatment p-values.

ts_msg("=== PSEUDOREPLICATION CHECK (trial vs tank) ===")

# Main-zone: tank-level check now uses logit(prop_flow) — compositionally correct.
# Trial-level model: res_zone_main (logit(prop_flow) ~ treatment * timepoint_f + (1|phys_trial_id)).
# Tank-level model: logit(mean_prop_flow) per tank (collapsed across all sessions) ~ treatment.
# CHANGE from prior version: was treatment × zone on long-format (compositionally flawed);
# now uses logit(prop_flow) on tank-averaged wide format, matching the new A1 model family.
.tank_zone_agg <- tryCatch({
  if (!is.null(df_main_wide)) {
    df_main_wide %>%
      dplyr::group_by(tank, treatment) %>%
      dplyr::summarise(mean_logit_flow = mean(logit_flow, na.rm = TRUE), .groups = "drop")
  } else {
    # Fallback: compute from long-format
    df_main_long %>%
      dplyr::filter(zone == "flow") %>%
      dplyr::group_by(tank, treatment) %>%
      dplyr::summarise(mean_p = mean(prop_time, na.rm = TRUE), .groups = "drop") %>%
      dplyr::mutate(mean_logit_flow = log(pmin(pmax(mean_p, 1e-4), 1 - 1e-4) /
                                          (1 - pmin(pmax(mean_p, 1e-4), 1 - 1e-4))))
  }
}, error = function(e) { warning("Tank zone agg failed: ", e$message); NULL })

.av_tank_zone <- tryCatch({
  m  <- lm(mean_logit_flow ~ treatment, data = .tank_zone_agg)
  a  <- car::Anova(m, type = "III")
  d  <- as.data.frame(a); d$term <- rownames(d); rownames(d) <- NULL
  names(d) <- c("ss", "df", "F_value", "p_value", "term"); d
}, error = function(e) { warning("Tank zone ANOVA failed: ", e$message); NULL })

.p_trial_zone <- .get_treatment_p(res_zone_main, "^treatment$")
.p_tank_zone  <- if (!is.null(.av_tank_zone)) {
  row <- .av_tank_zone[grepl("^treatment$", .av_tank_zone$term, perl = TRUE), ]
  if (nrow(row) > 0) as.numeric(row$p_value[1]) else NA_real_
} else NA_real_

# prop_active: tank-level (collapse timepoints)
.tank_active_agg <- df %>%
  dplyr::filter(is.finite(prop_active)) %>%
  dplyr::group_by(tank, treatment) %>%
  dplyr::summarise(mean_active = mean(prop_active, na.rm = TRUE), .groups = "drop")

.av_tank_active <- tryCatch({
  if (var(.tank_active_agg$mean_active, na.rm = TRUE) < 1e-8) {
    ts_msg("  Tank active: prop_active near-constant (",
           round(mean(.tank_active_agg$mean_active, na.rm = TRUE), 4),
           ") — ANOVA skipped.")
    NULL
  } else {
    m  <- lm(mean_active ~ treatment, data = .tank_active_agg)
    a  <- car::Anova(m, type = "III")
    d  <- as.data.frame(a); d$term <- rownames(d); rownames(d) <- NULL
    names(d) <- c("ss", "df", "F_value", "p_value", "term"); d
  }
}, error = function(e) {
  ts_msg("  Tank active ANOVA failed: ", conditionMessage(e))
  NULL
})

.p_trial_active <- .get_treatment_p(res_active_agg, "^treatment$")
.p_tank_active  <- if (!is.null(.av_tank_active)) {
  row <- .av_tank_active[grepl("^treatment$", .av_tank_active$term, perl = TRUE), ]
  if (nrow(row) > 0) as.numeric(row$p_value[1]) else NA_real_
} else NA_real_

# P3.4: replace binary Agree with effect-size table (F10)
# Extract treatment beta (log-odds or log-mean or raw coefficient) from each model
.get_trt_beta <- function(res, label) {
  tryCatch({
    m  <- res$model
    cf <- coef(summary(m))
    if (is.list(cf)) cf <- cf[[1]]   # glmmTMB returns list by component
    trt_rows <- grep("treatment", rownames(cf), ignore.case = TRUE)
    if (length(trt_rows) == 0) return(NA_real_)
    as.numeric(cf[trt_rows[1], 1])
  }, error = function(e) NA_real_)
}

.beta_trial_zone  <- .get_trt_beta(res_zone_main,  "A1")  # logit(prop_flow) LMM
.beta_trial_active <- .get_trt_beta(res_active_agg, "A6")

.beta_tank_zone <- tryCatch({
  m  <- lm(mean_logit_flow ~ treatment, data = .tank_zone_agg)
  cf <- coef(m); grep_trt <- grep("^treatment", names(cf)); if (length(grep_trt) > 0) cf[grep_trt[1]] else NA_real_
}, error = function(e) NA_real_)

.beta_tank_active <- tryCatch({
  if (!is.null(.av_tank_active)) {
    m  <- lm(mean_active ~ treatment, data = .tank_active_agg)
    cf <- coef(m); grep_trt <- grep("^treatment", names(cf)); if (length(grep_trt) > 0) cf[grep_trt[1]] else NA_real_
  } else NA_real_
}, error = function(e) NA_real_)

# Item 7 (prioritised corrections): extend pseudoreplication check from 2 to
# 6 rows so all 5 PRIMARY BH indicators are checked plus the prop_active
# diagnostic (A6 retained for transparency).
# Generic helper: collapse a session-level (or trial-level) data frame to a
# tank × treatment mean and fit lm(y ~ treatment). Returns p, beta, sign.
.tank_lm <- function(d, y_col, treat_col = "treatment") {
  if (is.null(d) || !y_col %in% names(d) || !"tank" %in% names(d))
    return(list(p = NA_real_, beta = NA_real_))
  tryCatch({
    d2 <- d %>%
          dplyr::group_by(tank, .data[[treat_col]]) %>%
          dplyr::summarise(.y = mean(.data[[y_col]], na.rm = TRUE), .groups = "drop")
    if (var(d2$.y, na.rm = TRUE) < 1e-10) return(list(p = NA_real_, beta = NA_real_))
    fmla <- as.formula(paste(".y ~", treat_col))
    m    <- lm(fmla, data = d2)
    a    <- car::Anova(m, type = "II")
    ad   <- as.data.frame(a); ad$term <- rownames(ad); rownames(ad) <- NULL
    pv_col <- grep("^Pr|^p.value", names(ad), ignore.case = TRUE, value = TRUE)[1]
    p_val  <- as.numeric(ad[grepl(treat_col, ad$term, ignore.case = TRUE), pv_col][1])
    cf     <- coef(m); ix <- grep(treat_col, names(cf), ignore.case = TRUE)
    list(p = p_val, beta = if (length(ix) > 0) as.numeric(cf[ix[1]]) else NA_real_)
  }, error = function(e) list(p = NA_real_, beta = NA_real_))
}

# Build the row set: name, response column, source data, trial-level model.
.pr_specs <- list(
  list(label = "A1c logit(prop_flow) x Treatment", y = "logit_flow",
       data = df_main_wide,        res = res_zone_flow_logit_agg),
  list(label = "A4  switches/session x Treatment", y = "switches_per_session",
       data = df,                  res = res_switch_agg),
  list(label = "A8  NND x Treatment",              y = "mean_nnd_cm",
       data = if (exists("df_gd"))     df_gd     else NULL, res = res_nnd_agg),
  list(label = "A16 School speed x Treatment",    y = "mean_centroid_spd_cm",
       data = if (exists("df_gd"))     df_gd     else NULL, res = res_cspd_agg),
  list(label = "A6  prop_active x Treatment",     y = "prop_active",
       data = df,                  res = res_active_agg)
)

.pr_rows <- lapply(.pr_specs, function(s) {
  trial_p    <- .get_treatment_p(s$res, "^treatment$")
  trial_beta <- .get_trt_beta(s$res, s$label)
  tank       <- .tank_lm(s$data, s$y)
  list(Analysis        = s$label,
       Trial_level_p   = round(trial_p,      4),
       Tank_level_p    = round(tank$p,       4),
       Trial_beta      = round(trial_beta,   4),
       Tank_beta       = round(tank$beta,    4),
       Beta_ratio      = round(abs(trial_beta / tank$beta), 3),
       Trial_sig       = isTRUE(trial_p < 0.05),
       Tank_sig        = isTRUE(tank$p  < 0.05),
       Direction_agree = isTRUE(sign(trial_beta) == sign(tank$beta)))
})
pseudorep_tbl <- do.call(rbind, lapply(.pr_rows, as.data.frame))
readr::write_csv(pseudorep_tbl, file.path(STEP5_OUT, "pseudoreplication_check.csv"))

.pseudorep_consistent <- all(pseudorep_tbl$Direction_agree, na.rm = TRUE)
if (.pseudorep_consistent) {
  cat("\nPseudoreplication check: trial-level vs tank-level effect directions CONSISTENT.\n\n")
} else {
  cat("\nPseudoreplication check: effect directions INCONSISTENT — interpret with caution.\n\n")
}

assign("pseudorep_tbl", pseudorep_tbl, envir = .GlobalEnv)


# =============================================================================
# ==== 5e) BENJAMINI-HOCHBERG CORRECTION ======================================
# =============================================================================
# BH FDR applied across all treatment p-values (A1-A22, NULL analyses excluded).

ts_msg("=== BH CORRECTION ===")

# No adjustment for multiple comparisons is applied anywhere in this study
# (author decision, 2026-08-07 — see Supplementary Table S7/2.7); every family
# below reports raw p only. The groupings retained here (zone preference,
# non-zone core indicators, Treatment:Timepoint interactions, sensitivity
# check, descriptive/exploratory) document the pre-specified analysis
# structure for the report — see the "6. Cross-Analysis Summary" section
# below for what each family contains.

BH_PRIMARY_ONLY <- TRUE  # kept for downstream compatibility

SIMPLIFY_SINGULAR_TO_LM <- FALSE

.bh_row <- function(label, res, pat = "^treatment$") {
  # Length-1 guard: every column must be exactly length 1 to avoid
  # data.frame() rejecting zero-length inputs.
  .one <- function(x, default = NA) {
    if (is.null(x) || length(x) == 0) default else x[1]
  }
  p_raw     <- .one(.get_treatment_p(res, pat), NA_real_)
  fam       <- .one(if (!is.null(res$family_used)) res$family_used else
                      if (!is.null(res$transform)) res$transform   else NA_character_,
                    NA_character_)
  re_fall   <- .one(if (!is.null(res$re_fallback)) res$re_fallback else NA, NA)
  stat      <- .one(if (!is.null(res$stat_type))  res$stat_type   else NA_character_,
                    NA_character_)
  r2m       <- .one(if (!is.null(res$r2$R2m))     round(res$r2$R2m, 3) else NA_real_,
                    NA_real_)
  norm_flag <- .one(if (!is.null(res$normality_flag)) res$normality_flag else NA, NA)

  # Extract df1, df2, and F/chi-sq for the matching ANOVA row (item 6).
  # Defensive: every extraction returns a length-1 numeric (or NA_real_) so
  # downstream data.frame() never sees a zero-length column.
  .scalar <- function(x) {
    if (is.null(x) || length(x) == 0) NA_real_ else as.numeric(x[1])
  }
  df1 <- df2 <- stat_val <- NA_real_
  if (!is.null(res$anova) && "term" %in% names(res$anova)) {
    av_row <- res$anova[grepl(pat, res$anova$term, ignore.case = TRUE, perl = TRUE), ]
    if (nrow(av_row) > 0) {
      df1      <- .scalar(av_row$df)
      df2      <- if ("df_denom" %in% names(av_row)) .scalar(av_row$df_denom) else NA_real_
      stat_val <- if ("chisq"    %in% names(av_row)) .scalar(av_row$chisq)    else NA_real_
    }
  }
  if (length(p_raw) != 1) p_raw <- NA_real_

  data.frame(Analysis = label, p_raw = p_raw,
             family = fam, re_fallback = re_fall,
             stat_type = stat,
             df1 = df1, df2 = round(df2, 2), stat_val = round(stat_val, 3),
             R2m = r2m,
             normality_flag = norm_flag,
             stringsAsFactors = FALSE)
}

# Zone preference: 8 Jacobs beta-GLMM tests (4 zones x {trial-aggregated,
# timepoint-resolved}). No adjustment for multiple comparisons is applied
# anywhere in this study (author decision, 2026-08-07 -- see Table S7/2.7);
# these families are retained purely as pre-specified organisational
# groupings for the report, not as a BH-correction scheme. The alr/logit
# zone-preference tests (A1c + A2 sub-zone) are reported descriptively
# (section 6.5); the non-zone indicators (switches, NND,
# centroid speed) are grouped separately in bh_suppl_primary.
bh_primary <- dplyr::bind_rows(
  .bh_row("J_beta_flow_agg   Jacobs D(flow)   x Treatment (agg)",   res_jac_beta_flow_agg),
  .bh_row("J_beta_high_agg   Jacobs D(high)   x Treatment (agg)",   res_jac_beta_high_agg),
  .bh_row("J_beta_med_agg    Jacobs D(medium) x Treatment (agg)",   res_jac_beta_med_agg),
  .bh_row("J_beta_low_agg    Jacobs D(low)    x Treatment (agg)",   res_jac_beta_low_agg),
  .bh_row("J_beta_flow_tp    Jacobs D(flow)   x Treatment:Tp",      res_jac_beta_flow_tp,
          "treatment.*timepoint|timepoint.*treatment"),
  .bh_row("J_beta_high_tp    Jacobs D(high)   x Treatment:Tp",      res_jac_beta_high_tp,
          "treatment.*timepoint|timepoint.*treatment"),
  .bh_row("J_beta_med_tp     Jacobs D(medium) x Treatment:Tp",      res_jac_beta_med_tp,
          "treatment.*timepoint|timepoint.*treatment"),
  .bh_row("J_beta_low_tp     Jacobs D(low)    x Treatment:Tp",      res_jac_beta_low_tp,
          "treatment.*timepoint|timepoint.*treatment")
)

# Joint Jacobs D models (treatment × zone): add interaction-term rows.
bh_primary <- dplyr::bind_rows(
  bh_primary,
  .bh_row("J_beta_main_joint_agg  Jacobs D main Treatment × Zone (agg)",
          res_jac_beta_main_joint_agg, "treatment.*zone|zone.*treatment"),
  .bh_row("J_beta_sec_joint_agg   Jacobs D sub  Treatment × Zone (agg)",
          res_jac_beta_sec_joint_agg,  "treatment.*zone|zone.*treatment")
)

# Non-zone core indicators (switches, NND, centroid speed),
# grouped separately from the zone-preference tests above. 4 tests; raw p,
# no adjustment for multiple comparisons (see note above bh_summary below).
bh_suppl_primary <- dplyr::bind_rows(
  .bh_row("A4  Switches/session x Treatment (agg)",        res_switch_agg),
  .bh_row("A8  NND x Treatment (agg)",                     res_nnd_agg),
  .bh_row("A16 School speed x Treatment (agg)",          res_cspd_agg)
)

# SUPPLEMENTARY-TP: Treatment:Timepoint interactions for the 4 non-zone-preference
# core metrics (A3/A7/A9/A15). A1 (logit prop_flow x Tp) is reported descriptively
# alongside the other alr/logit zone tests (section 6.5). Raw p throughout, no
# adjustment for multiple comparisons.
bh_suppl <- dplyr::bind_rows(
  .bh_row("A3  Switches x Treatment:Tp",
          res_switch_tp, "treatment.*timepoint|timepoint.*treatment"),
  .bh_row("A7  NND x Treatment:Tp",
          res_nnd_tp,    "treatment.*timepoint|timepoint.*treatment"),
  .bh_row("A15 School speed x Treatment:Tp",
          res_cspd_tp,   "treatment.*timepoint|timepoint.*treatment")
)

# JACOBS sensitivity family (atanh path, no BH).
# Reported alongside primary beta-GLMM as a robustness check.
bh_jacobs_sens <- dplyr::bind_rows(
  .bh_row("J_atanh_flow_agg  Jacobs atanh(D)(flow)   x Tx (agg)",   res_jac_flow_agg),
  .bh_row("J_atanh_high_agg  Jacobs atanh(D)(high)   x Tx (agg)",   res_jac_high_agg),
  .bh_row("J_atanh_med_agg   Jacobs atanh(D)(medium) x Tx (agg)",   res_jac_med_agg),
  .bh_row("J_atanh_low_agg   Jacobs atanh(D)(low)    x Tx (agg)",   res_jac_low_agg),
  .bh_row("J_atanh_flow_tp   Jacobs atanh(D)(flow)   x Tx:Tp",      res_jac_flow_tp,
          "treatment.*timepoint|timepoint.*treatment"),
  .bh_row("J_atanh_high_tp   Jacobs atanh(D)(high)   x Tx:Tp",      res_jac_high_tp,
          "treatment.*timepoint|timepoint.*treatment"),
  .bh_row("J_atanh_med_tp    Jacobs atanh(D)(medium) x Tx:Tp",      res_jac_med_tp,
          "treatment.*timepoint|timepoint.*treatment"),
  .bh_row("J_atanh_low_tp    Jacobs atanh(D)(low)    x Tx:Tp",      res_jac_low_tp,
          "treatment.*timepoint|timepoint.*treatment")
)

# DESCRIPTIVE / EXPLORATORY: raw p only, no adjustment for multiple comparisons.
# Includes the alr/logit zone-preference family (A1c, A1, A2 sub-zone
# log-ratios), marked with explicit "[desc]" tags, alongside the other
# exploratory indicators (IID, hull area, flux, prop_active).
bh_exploratory <- dplyr::bind_rows(
  # Demoted alr / logit zone-preference (descriptive only)
  .bh_row("A1c [desc] logit(prop_flow) x Treatment (agg)",       res_zone_flow_logit_agg),
  .bh_row("A1  [desc] logit(prop_flow) x Treatment:Tp",          res_zone_main,
          "treatment.*timepoint|timepoint.*treatment"),
  .bh_row("A2_lr_high   [desc] lr_high   x Treatment (agg)",     res_zone_sec_high_agg),
  .bh_row("A2_lr_med    [desc] lr_medium x Treatment (agg)",     res_zone_sec_med_agg),
  .bh_row("A2_lr_low    [desc] lr_low    x Treatment (agg)",     res_zone_sec_low_agg),
  .bh_row("A2_lr_high_tp [desc] lr_high   x Treatment:Tp",       res_zone_sec_high_tp,
          "treatment.*timepoint|timepoint.*treatment"),
  .bh_row("A2_lr_med_tp  [desc] lr_medium x Treatment:Tp",       res_zone_sec_med_tp,
          "treatment.*timepoint|timepoint.*treatment"),
  .bh_row("A2_lr_low_tp  [desc] lr_low    x Treatment:Tp",       res_zone_sec_low_tp,
          "treatment.*timepoint|timepoint.*treatment"),
  # Other exploratory (unchanged)
  .bh_row("A6  prop_active x Treatment (agg)",         res_active_agg),
  .bh_row("A_flux_agg zone_flux x Treatment (agg)",    res_flux_agg),
  .bh_row("A12 IID x Treatment (agg)",                 res_iid_agg),
  .bh_row("A14 Hull area x Treatment (agg)",           res_hull_agg),
  .bh_row("A5  prop_active x Treatment x Tp",          res_active_tp),
  .bh_row("A_flux_tp zone_flux x Treatment x Tp",      res_flux_tp),
  .bh_row("A11 IID x Treatment x Tp",                  res_iid_tp),
  .bh_row("A13 Hull area x Treatment x Tp",            res_hull_tp),
  # Legacy descriptive-only (compositionally misspecified — retained for reference)
  .bh_row("A1_agg [desc] Main-zone beta GLMM (agg)",  res_zone_main_agg)
)
bh_exploratory$p_BH    <- NA_real_
bh_exploratory$sig_raw <- !is.na(bh_exploratory$p_raw) & bh_exploratory$p_raw < 0.05
bh_exploratory$sig_BH  <- NA
bh_exploratory$set     <- "descriptive"

# -----------------------------------------------------------------------------
# BH CORRECTION DISABLED (author decision, 2026-08-09)
# -----------------------------------------------------------------------------
# Every family below now reports RAW p only; p_BH / sig_BH are NA throughout,
# exactly as bh_jacobs_sens and bh_exploratory already did.
#
# WHY: polarisation was removed from the manuscript (it is the only indicator
# whose computation requires persistent identity tracking -- heading is the
# frame-to-frame displacement of a given fish_id, so an identity swap corrupts
# it; see group_dynamics_STEP2b_choice_exp.R). Polarisation was a member of
# BOTH bh_suppl_primary (A10) and bh_suppl (A9), so dropping it re-ranks each
# family from 4 tests to 3 and mechanically INFLATES the neighbours' adjusted
# p-values: A8 NND x Treatment moved 0.0414 -> 0.0621 and would have lost BH
# significance purely because its BH rank changed (multiplier 4/2 -> 3/1), not
# because anything about the NND model changed. Making a retained result's
# significance contingent on the presence of a discarded, methodologically
# invalid test is not defensible, so the correction is withdrawn rather than
# recomputed on the reduced families.
#
# This is CONSISTENT with the rest of the study as already written: the
# endocrine engine (Table S8), the full monoamine grid (Table S10) and the
# Jacobs preference family (Table S12) all state "No adjustment for multiple
# comparisons is applied; p is raw." Nothing in the manuscript or the
# supplementary materials reports a BH-adjusted p for the behaviour families,
# so no text change is required -- verified 2026-08-09 by scanning both .docx
# files for Benjamini/Hochberg/FDR/p_BH (zero hits in the reported tables).
#
# The family groupings themselves are retained below (as `set` labels) because
# they still document which indicators were pre-specified as primary vs
# supplementary vs descriptive -- that structure is reported, the adjustment
# is not. To re-enable, restore the p.adjust() calls commented out below.
bh_summary           <- bh_primary
# bh_summary$p_BH    <- p.adjust(bh_summary$p_raw, method = "BH")   # disabled 2026-08-09
bh_summary$p_BH      <- NA_real_
bh_summary$sig_raw   <- !is.na(bh_summary$p_raw) & bh_summary$p_raw < 0.05
bh_summary$sig_BH    <- NA
bh_summary$set       <- "primary_jacobs_beta"

# bh_suppl_primary$p_BH    <- p.adjust(bh_suppl_primary$p_raw, method = "BH")  # disabled 2026-08-09
bh_suppl_primary$p_BH      <- NA_real_
bh_suppl_primary$sig_raw   <- !is.na(bh_suppl_primary$p_raw) & bh_suppl_primary$p_raw < 0.05
bh_suppl_primary$sig_BH    <- NA
bh_suppl_primary$set       <- "primary_nonzone"

# bh_suppl$p_BH      <- p.adjust(bh_suppl$p_raw, method = "BH")     # disabled 2026-08-09
bh_suppl$p_BH        <- NA_real_
bh_suppl$sig_raw     <- !is.na(bh_suppl$p_raw) & bh_suppl$p_raw < 0.05
bh_suppl$sig_BH      <- NA
bh_suppl$set         <- "supplementary_tp"

bh_jacobs_sens$p_BH    <- NA_real_
bh_jacobs_sens$sig_raw <- !is.na(bh_jacobs_sens$p_raw) & bh_jacobs_sens$p_raw < 0.05
bh_jacobs_sens$sig_BH  <- NA
bh_jacobs_sens$set     <- "jacobs_sensitivity_atanh"

bh_summary$p_raw         <- round(bh_summary$p_raw, 4)
bh_summary$p_BH          <- round(bh_summary$p_BH,  4)
bh_suppl_primary$p_raw   <- round(bh_suppl_primary$p_raw, 4)
bh_suppl_primary$p_BH    <- round(bh_suppl_primary$p_BH,  4)
bh_suppl$p_raw           <- round(bh_suppl$p_raw, 4)
bh_suppl$p_BH            <- round(bh_suppl$p_BH,  4)
bh_jacobs_sens$p_raw     <- round(bh_jacobs_sens$p_raw, 4)
bh_exploratory$p_raw     <- round(bh_exploratory$p_raw, 4)

readr::write_csv(bh_summary,        file.path(STEP5_OUT, "bh_correction_summary.csv"))
readr::write_csv(bh_suppl_primary,  file.path(STEP5_OUT, "bh_primary_nonzone.csv"))
readr::write_csv(bh_suppl,          file.path(STEP5_OUT, "bh_supplementary.csv"))
readr::write_csv(bh_jacobs_sens,    file.path(STEP5_OUT, "bh_jacobs_sensitivity.csv"))
readr::write_csv(bh_jacobs_sens,    file.path(STEP5_OUT, "jacobs", "bh_jacobs_sensitivity.csv"))
readr::write_csv(bh_exploratory,    file.path(STEP5_OUT, "bh_exploratory.csv"))
ts_msg("BH summary written (", nrow(bh_summary), " primary Jacobs beta | ",
       nrow(bh_suppl_primary), " primary non-zone | ",
       nrow(bh_suppl), " supplementary-TP | ",
       nrow(bh_jacobs_sens), " Jacobs atanh (no BH) | ",
       nrow(bh_exploratory), " descriptive)")

# =============================================================================
# ==== M4 (R2-5): TRIAL ORDER / DAY-OF-TESTING EFFECT ========================
# =============================================================================
# Each aggregated preference-experiment model is refit with order_idx
# (chronological rank of trial_date, 1..8) and its interaction with treatment
# added: y ~ treatment * order_idx + RE (same RE candidates as the primary fit).
#
# DESIGN LIMITATION (documented per M4's "test AND document" requirement, see
# design_balance.csv / the date-nested-in-tank note above): order_idx is not
# identifiable independently of tank or of fish_density in this design.
#   - Between tanks: each tank occupies its own exclusive, non-overlapping pair
#     of testing days, so the day sequence is a perfect function of tank
#     identity — order and tank cannot both be estimated when tank is also a
#     random effect (this is why RE_WIDE_AGG already excludes trial_date).
#   - Within a tank: the two testing days always run fish_density = {16,12}
#     then {8,4}, and the first trial each day always has the higher of that
#     day's two densities — order_idx and fish_density are related by an exact
#     affine transform in this design. An order test and a density test are
#     mathematically the same test here and cannot be told apart.
# Terms below should be read as "order-or-density", not as independent
# evidence of chronological drift.
ALL_TRIAL_DATES_BEH <- as.Date(sprintf("2024-06-%02d", 15:22))
.add_order_idx <- function(data) {
  if (is.null(data) || !"trial_date" %in% names(data)) return(data)
  data$order_idx <- as.numeric(factor(as.Date(substr(as.character(data$trial_date), 1, 10)),
                                      levels = ALL_TRIAL_DATES_BEH))
  data
}

run_order_sensitivity <- function(label, data, response, re_cands, is_beta = FALSE) {
  d2 <- .add_order_idx(data)
  if (is.null(d2) || !"order_idx" %in% names(d2) || all(is.na(d2$order_idx))) {
    ts_msg("  [", label, "_order_test] skipped: order_idx unavailable")
    return(NULL)
  }
  fn <- if (isTRUE(is_beta)) run_betaglmm_analysis else run_lmm_analysis
  tryCatch(
    fn(label = paste0(label, "_order_test"), data = d2, response = response,
       fixed_str = "treatment * order_idx", re_cands = re_cands,
       focal_terms = list("Treatment" = "^treatment$", "Order" = "^order_idx$",
                          "Treatment:Order" = "treatment.*order|order.*treatment")),
    error = function(e) { warning(label, " order test failed: ", e$message); NULL }
  )
}

.order_specs <- list(
  list(label = "zone_flow_logit_aggregated", data = df_main_wide_agg, response = "logit_flow",   re_cands = RE_WIDE_AGG),
  list(label = "switches_aggregated",        data = df_agg,          response = "switches_per_session", re_cands = RE_WIDE_AGG),
  list(label = "active_aggregated",          data = df_agg,          response = "prop_active",   re_cands = RE_WIDE_AGG, is_beta = TRUE),
  list(label = "nnd_aggregated",             data = df_gd_agg,       response = "mean_nnd_cm",         re_cands = RE_WIDE_AGG),
  list(label = "iid_aggregated",             data = df_gd_agg,       response = "mean_iid_cm",         re_cands = RE_WIDE_AGG),
  list(label = "hull_area_aggregated",       data = df_gd_agg,       response = "mean_hull_area_cm2",  re_cands = RE_WIDE_AGG),
  list(label = "centroid_speed_aggregated",  data = df_gd_agg,       response = "mean_centroid_spd_cm", re_cands = RE_WIDE_AGG)
)

ts_msg("=== M4: order/day-of-testing sensitivity (", length(.order_specs), " aggregated models) ===")
.order_results <- lapply(.order_specs, function(s) {
  if (is.null(s$data)) return(NULL)
  run_order_sensitivity(s$label, s$data, s$response, s$re_cands, isTRUE(s$is_beta))
})
names(.order_results) <- vapply(.order_specs, `[[`, character(1), "label")

.extract_order_terms <- function(res, label) {
  if (is.null(res) || is.null(res$anova)) return(NULL)
  a <- res$anova
  a[grepl("order", a$term, ignore.case = TRUE), , drop = FALSE] %>%
    dplyr::mutate(model = label, .before = 1)
}
order_effect_summary <- dplyr::bind_rows(
  Map(.extract_order_terms, .order_results, names(.order_results))
)
readr::write_csv(order_effect_summary, file.path(STEP5_OUT, "order_effect_summary.csv"))
ts_msg("Order/day effect (M4): ", nrow(order_effect_summary), " term-rows across ",
       sum(!vapply(.order_results, is.null, logical(1))), " models -> order_effect_summary.csv")

# Back-compat shim: some downstream report builders expect bh_subzone /
# bh_jacobs in the global env. Provide empty placeholders so they degrade
# gracefully rather than erroring.
bh_subzone <- bh_exploratory[grepl("^A2_lr", bh_exploratory$Analysis), , drop = FALSE]
bh_jacobs  <- bh_jacobs_sens
assign("bh_subzone",     bh_subzone,     envir = .GlobalEnv)
assign("bh_jacobs",      bh_jacobs,      envir = .GlobalEnv)
assign("bh_jacobs_sens", bh_jacobs_sens, envir = .GlobalEnv)
assign("bh_suppl_primary", bh_suppl_primary, envir = .GlobalEnv)
assign("bh_summary",     bh_summary,     envir = .GlobalEnv)

# ---- Hedges' g summary table (item #7) ------------------------------------
# Computed for all analyses that use gaussian LMM; NA for GLMM families.
ts_msg("=== HEDGES' g ===")
.hg_analyses <- list(
  # Primary inference models (compositionally correct)
  list(label = "A1c logit(prop_flow) x Treatment (agg)",    res = res_zone_flow_logit_agg),
  list(label = "A4  Switches/session x Treatment (agg)",    res = res_switch_agg),
  list(label = "A8  NND x Treatment (agg)",                 res = res_nnd_agg),
  list(label = "A16 School speed x Treatment (agg)",      res = res_cspd_agg),
  # Sub-zone log-ratio models (aggregated)
  list(label = "A2_lr_high  lr_high x Treatment (agg)",     res = res_zone_sec_high_agg),
  list(label = "A2_lr_med   lr_medium x Treatment (agg)",   res = res_zone_sec_med_agg),
  list(label = "A2_lr_low   lr_low x Treatment (agg)",      res = res_zone_sec_low_agg),
  # Other exploratory
  list(label = "A12 IID x Treatment (agg)",                 res = res_iid_agg),
  list(label = "A14 Hull area x Treatment (agg)",           res = res_hull_agg),
  list(label = "A_flux_agg zone_flux x Treatment (agg)",    res = res_flux_agg),
  # Legacy (descriptive only — beta GLMM, compositionally misspecified)
  list(label = "A1  [desc] Main-zone occ. x Treatment (agg)", res = res_zone_main_agg)
)
hedges_g_tbl <- do.call(rbind, lapply(.hg_analyses, function(x) {
  hg <- .calc_hedges_g(x$res, x$label)
  data.frame(Analysis = x$label, g = hg$g, ci_lo = hg$ci_lo, ci_hi = hg$ci_hi,
             g_resid = hg$g_resid, ci_lo_resid = hg$ci_lo_resid, ci_hi_resid = hg$ci_hi_resid,
             delta = hg$delta, sigma_r = hg$sigma_r, sigma_total = hg$sigma_total,
             J = hg$J, method = hg$method,
             stringsAsFactors = FALSE)
}))
readr::write_csv(hedges_g_tbl, file.path(STEP5_OUT, "hedges_g_summary.csv"))
ts_msg("Hedges' g written: ", nrow(hedges_g_tbl[!is.na(hedges_g_tbl$g), ]),
       " valid / ", nrow(hedges_g_tbl), " total")

# ---- Sub-zone alr pairwise contrasts summary -------------------------------
# Pulls the zone_lr Tukey contrasts from the stacked alr LMMs (A2_pairs)
# into a single descriptive table. Exploratory only (no BH).
ts_msg("=== Sub-zone alr pairwise contrasts summary ===")
.lr_pairs_table <- function(res, stratum_lbl) {
  if (is.null(res) || is.null(res$posthoc)) return(NULL)
  ph <- res$posthoc[grepl("^zone_lr$", names(res$posthoc), perl = TRUE)]
  if (length(ph) == 0) return(NULL)
  contr <- ph[[1]]$contrasts
  if (is.null(contr) || nrow(contr) == 0) return(NULL)
  contr$stratum <- stratum_lbl
  contr
}
.lr_pairs_summary <- dplyr::bind_rows(
  .lr_pairs_table(res_zone_sec_lr_pairs_agg, "agg"),
  .lr_pairs_table(res_zone_sec_lr_pairs_tp,  "tp")
)
if (!is.null(.lr_pairs_summary) && nrow(.lr_pairs_summary) > 0) {
  readr::write_csv(.lr_pairs_summary,
                   file.path(STEP5_OUT, "subzone_lr_pairwise_summary.csv"))
  ts_msg("  subzone_lr_pairwise_summary.csv written: ",
         nrow(.lr_pairs_summary), " rows")
} else {
  ts_msg("  subzone_lr_pairwise_summary: no contrasts available.")
}


# =============================================================================
# ==== 6) SKINNY GRAPHS =======================================================
# =============================================================================

ts_msg("=== SKINNY GRAPHS ===")

.sg_dir <- file.path(STEP5_OUT, "skinny_graphs")
.sg_png <- file.path(.sg_dir, "png_preview")
dir.create(.sg_png, recursive = TRUE, showWarnings = FALSE)

# ---- Plot settings -----------------------------------------------------------
TREATMENT_COLORS <- c("control" = "#2166AC", "exercise choice" = "#D6604D")
TREATMENT_SHAPES <- c("control" = 17L,       "exercise choice" = 16L)
ZONE_COLORS_MAIN <- c("flow" = "#4DAC26",    "calm" = "#8073AC")

JITTER_W             <- 0.12
DODGE_W              <- 0.70
PT_ALPHA             <- 0.80
PT_SIZE              <- 2.4
MEAN_W               <- 0.40
LW_MEAN              <- 0.85
LW_ERR               <- 0.65
ERR_W                <- 0.13
CLD_SIZE             <- 4.0
CLD_NUDGE_SEM_MULT   <- 1.6
CLD_NUDGE_RANGE_FRAC <- 0.05

BASE_THEME <- ggplot2::theme_minimal(base_size = 13) +
  ggplot2::theme(
    panel.grid        = ggplot2::element_blank(),
    panel.background  = ggplot2::element_blank(),
    axis.line.x       = ggplot2::element_line(color = "black", linewidth = 0.85),
    axis.line.y       = ggplot2::element_line(color = "black", linewidth = 0.85),
    axis.ticks        = ggplot2::element_line(color = "black", linewidth = 0.7),
    axis.text.x       = ggplot2::element_text(size = 13, face = "bold", color = "black"),
    axis.text.y       = ggplot2::element_text(size = 11, color = "black"),
    axis.title        = ggplot2::element_text(size = 13, face = "bold", color = "black"),
    axis.title.x      = ggplot2::element_blank(),
    legend.title      = ggplot2::element_text(size = 11, face = "bold"),
    legend.text       = ggplot2::element_text(size = 10),
    strip.text        = ggplot2::element_text(size = 12, face = "bold"),
    strip.background  = ggplot2::element_blank(),
    plot.margin       = ggplot2::margin(8, 10, 16, 10),
    plot.caption      = ggtext::element_markdown(size = 14, hjust = 0,
                                                 face = "italic", lineheight = 1.3,
                                                 margin = ggplot2::margin(t = 10)),
    plot.subtitle     = ggplot2::element_text(size = 9,  hjust = 0, face = "italic",
                                              colour = "grey40")
  )

TIMEPOINT_LABELS <- c("1" = "5\u201325 min",
                       "2" = "45\u201365 min",
                       "3" = "85\u2013105 min")

# Zone color maps (user-defined)
color_map       <- c(flow = "red4", calm = "lightskyblue3", high = "tomato3",
                     medium = "mediumturquoise", low = "goldenrod3")
color_map_broad <- c(flow = "red4", calm = "lightskyblue3")

# Okabe palette for treatment-colored indicator line plots
.pal_tmp <- tryCatch(as.character(cols4all::c4a("okabe")),
                     error = function(e) c("#E69F00","#56B4E9","#009E73","#F0E442",
                                           "#0072B2","#D55E00","#CC79A7","#000000"))
TREATMENT_COLORS_PAL <- c("control" = .pal_tmp[6], "exercise choice" = .pal_tmp[5])
rm(.pal_tmp)

LINE_TP_LABELS <- c("1" = "5\u201325 min", "2" = "45\u201365 min", "3" = "85\u2013105 min")

# ---- Summary helpers ---------------------------------------------------------
.sem <- function(x) { x <- x[is.finite(x)]; if (length(x) < 2) NA_real_ else sd(x)/sqrt(length(x)) }

.smry_zone <- function(d, y_col, groups = c("zone","treatment")) {
  d %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(groups))) %>%
    dplyr::summarise(mean_y = mean(.data[[y_col]], na.rm=TRUE),
                     sem_y  = .sem(.data[[y_col]]),
                     .groups = "drop")
}

.smry_std <- function(d, y_col, groups = "treatment") {
  d %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(groups))) %>%
    dplyr::summarise(mean_y = mean(.data[[y_col]], na.rm=TRUE),
                     sem_y  = .sem(.data[[y_col]]),
                     .groups = "drop")
}

# ---- CLD annotation helper ---------------------------------------------------
.get_cld_df <- function(analysis_res, term_key, y_nudge_frac = 0.08) {
  if (is.null(analysis_res) || is.null(analysis_res$posthoc)) return(NULL)
  ph_names <- names(analysis_res$posthoc)
  key_clean <- gsub("_", ":", gsub("x", ":", term_key))
  matched <- ph_names[grepl(gsub("\\.", ".*", key_clean), ph_names, ignore.case = TRUE)]
  if (length(matched) == 0) matched <- ph_names[1]
  if (length(matched) == 0 || is.null(analysis_res$posthoc[[matched[1]]])) return(NULL)
  cld <- analysis_res$posthoc[[matched[1]]]$cld
  if (is.null(cld) || !".group" %in% names(cld)) return(NULL)
  cld$.group <- trimws(cld$.group)
  cld
}

# ---- Zone plot builder -------------------------------------------------------
.make_zone_plot <- function(df_long, y_col, zone_colors, y_label, title,
                             caption_txt = NULL, cld_df = NULL, is_sig = FALSE,
                             facet_tp = FALSE, pct_scale = FALSE, n_txt = NULL) {
  groups_smry <- if (facet_tp) c("zone","treatment","timepoint_f") else c("zone","treatment")
  smry <- .smry_zone(df_long, y_col, groups_smry)

  POS_JD    <- ggplot2::position_jitterdodge(jitter.width=JITTER_W, jitter.height=0, dodge.width=DODGE_W)
  POS_DODGE <- ggplot2::position_dodge(width=DODGE_W)

  p <- ggplot2::ggplot(df_long, ggplot2::aes(
      x = zone, y = .data[[y_col]], colour = zone,
      shape = treatment, group = treatment)) +
    ggplot2::geom_point(position = POS_JD, alpha = PT_ALPHA, size = PT_SIZE) +
    ggplot2::geom_crossbar(
      data = smry,
      ggplot2::aes(y = mean_y, ymin = mean_y, ymax = mean_y, group = treatment),
      colour = "black",
      position = POS_DODGE, width = MEAN_W, linewidth = LW_MEAN, show.legend = FALSE) +
    ggplot2::geom_errorbar(
      data = smry,
      ggplot2::aes(y = mean_y, ymin = mean_y - sem_y, ymax = mean_y + sem_y,
                   group = treatment),
      colour = "black",
      position = POS_DODGE, width = ERR_W, linewidth = LW_ERR, show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = zone_colors, name = "Zone") +
    ggplot2::scale_shape_manual(values  = TREATMENT_SHAPES, name = "Treatment") +
    BASE_THEME +
    ggplot2::labs(title = title, x = "Zone", y = y_label,
                  caption = .compose_caption(n_txt, caption_txt))

  if (isTRUE(pct_scale)) {
    p <- p +
      ggplot2::scale_y_continuous(
        labels = function(x) round(x * 100),
        breaks = seq(0, 1, by = 0.2)) +
      ggplot2::expand_limits(y = c(0, 1))
  }

  # CLD letters: drawn whenever the cld table has structure and the letters
  # are not all identical (some pairwise contrasts were significant). This
  # holds even when the omnibus ANOVA term is n.s. — in that case the caller
  # also adds a disclosure note via .cap_zone(). If all letters are identical
  # the letters convey nothing and are omitted.
  if (!is.null(cld_df) && ".group" %in% names(cld_df) && nrow(cld_df) > 2) {
    .cld_unique <- unique(trimws(as.character(cld_df$.group[!is.na(cld_df$.group)])))
    if (length(.cld_unique) > 1) {
      .raw_y <- df_long[[y_col]][is.finite(df_long[[y_col]])]
      .y_max <- max(.raw_y, na.rm = TRUE)
      .y_min <- min(.raw_y, na.rm = TRUE)
      .y_rng <- max(.y_max - .y_min, 1e-10)
      # Spec: place CLD at 15% of data range above the highest scatter point;
      # all letters at the same vertical level.
      .lbl_y <- .y_max + 0.15 * .y_rng
      .y_top <- .lbl_y  + 0.05 * .y_rng
      .y_bot <- .y_min  - 0.02 * .y_rng
      cld_merge_cols <- intersect(c("zone","treatment","timepoint_f"), names(cld_df))
      cld_join <- dplyr::left_join(smry, cld_df[, c(cld_merge_cols, ".group")], by = cld_merge_cols)
      if (nrow(cld_join) > 0 && ".group" %in% names(cld_join)) {
        cld_join$label_y <- .lbl_y
        p <- p +
          ggplot2::expand_limits(y = c(.y_bot, .y_top)) +
          ggplot2::geom_text(
            data = cld_join,
            ggplot2::aes(x = zone, y = label_y, label = .group, group = treatment),
            position = POS_DODGE, size = CLD_SIZE, fontface = "bold",
            colour = "black", show.legend = FALSE)
      }
    }
  }

  if (facet_tp)
    p <- p + ggplot2::facet_wrap(
      ~timepoint_f, ncol = 1,
      labeller = ggplot2::labeller(timepoint_f = TIMEPOINT_LABELS))
  p
}

# ---- Standard scatter -------------------------------------------------------
.make_std_plot <- function(df, y_col, y_label, title, caption_txt = NULL,
                            cld_df = NULL, is_sig = FALSE, facet_tp = FALSE,
                            n_txt = NULL, sig_textsize = 4) {
  groups_smry <- if (facet_tp) c("treatment","timepoint_f") else "treatment"
  smry <- .smry_std(df, y_col, groups_smry)

  # Map treatment to numeric x so annotations can use x = 1.5 (rvg-safe).
  .trt_lvls <- levels(factor(df$treatment))
  .n_trt    <- length(.trt_lvls)
  .trt_map  <- stats::setNames(seq_len(.n_trt), .trt_lvls)
  df_plot   <- df[is.finite(df[[y_col]]), ]
  df_plot$.x_pos <- unname(.trt_map[as.character(df_plot$treatment)])
  smry$.x_pos    <- unname(.trt_map[as.character(smry$treatment)])

  POS_JITTER <- ggplot2::position_jitter(width = JITTER_W, height = 0)

  p <- ggplot2::ggplot(df_plot, ggplot2::aes(
      x = .x_pos, y = .data[[y_col]],
      colour = treatment, shape = treatment)) +
    ggplot2::geom_point(position = POS_JITTER, alpha = PT_ALPHA, size = PT_SIZE) +
    ggplot2::geom_crossbar(
      data = smry,
      ggplot2::aes(x = .x_pos, y = mean_y, ymin = mean_y, ymax = mean_y),
      colour = "black",
      width = MEAN_W, linewidth = LW_MEAN, show.legend = FALSE) +
    ggplot2::geom_errorbar(
      data = smry,
      ggplot2::aes(x = .x_pos, y = mean_y, ymin = mean_y - sem_y, ymax = mean_y + sem_y),
      colour = "black",
      width = ERR_W, linewidth = LW_ERR, show.legend = FALSE) +
    ggplot2::scale_x_continuous(
      breaks = seq_len(.n_trt),
      labels = stringr::str_to_title(.trt_lvls),
      limits = c(0.5, .n_trt + 0.5),
      expand = c(0, 0)) +
    ggplot2::scale_colour_manual(values = TREATMENT_COLORS_PAL, name = "Treatment") +
    ggplot2::scale_shape_manual(values  = TREATMENT_SHAPES, name = "Treatment") +
    BASE_THEME +
    ggplot2::labs(title = title, x = "Treatment", y = y_label,
                  caption = .compose_caption(n_txt, caption_txt))

  # Significance bracket for two-group comparisons via ggsignif.
  # Bracket sits 10% of the data range above the highest plotted point.
  if (isTRUE(is_sig) && .n_trt == 2L) {
    .raw_y  <- df[[y_col]][is.finite(df[[y_col]])]
    .y_max  <- max(.raw_y, na.rm = TRUE)
    .y_rng  <- max(.y_max - min(.raw_y, na.rm = TRUE), 1e-10)
    .bar_y  <- .y_max + 0.10 * .y_rng

    .sig_df <- data.frame(xmin = 1, xmax = 2,
                          y_position = .bar_y, annotations = "*")
    p <- p +
      ggplot2::expand_limits(y = .y_max + 0.22 * .y_rng) +
      ggsignif::geom_signif(
        data        = .sig_df,
        mapping     = ggplot2::aes(xmin = xmin, xmax = xmax,
                                   y_position = y_position,
                                   annotations = annotations),
        manual      = TRUE,
        tip_length  = 0.02,
        size        = 0.5,
        textsize    = sig_textsize,
        vjust       = 0.5,
        color       = "black",
        inherit.aes = FALSE
      )
  }

  # Suppress y-axis ticks below zero (values are always non-negative).
  p <- p + ggplot2::scale_y_continuous(
    breaks = function(lims) { b <- pretty(lims); b[b >= 0] }
  )

  if (facet_tp)
    p <- p + ggplot2::facet_wrap(
      ~timepoint_f, ncol = 1,
      labeller = ggplot2::labeller(timepoint_f = TIMEPOINT_LABELS))
  p
}

# ---- Zone line plot: x = timepoint, mean ± SEM per zone × treatment ---------
# color = zone (zone color map), shape = treatment, linetype = treatment
# error bars and mean cross bars black
.make_zone_line_plot <- function(df_long, y_col, zone_colors, y_label, title,
                                  pct_scale = FALSE, n_txt = NULL) {
  if (is.null(df_long) || !y_col %in% names(df_long)) return(NULL)
  smry <- df_long %>%
    dplyr::filter(is.finite(.data[[y_col]]), !is.na(timepoint_f)) %>%
    dplyr::group_by(zone, treatment, timepoint_f) %>%
    dplyr::summarise(mean_y = mean(.data[[y_col]], na.rm = TRUE),
                     sem_y  = .sem(.data[[y_col]]),
                     .groups = "drop") %>%
    dplyr::mutate(tp_num = as.integer(as.character(timepoint_f)))
  if (nrow(smry) == 0) return(NULL)

  p <- ggplot2::ggplot(smry, ggplot2::aes(
      x        = tp_num,
      y        = mean_y,
      colour   = zone,
      shape    = treatment,
      linetype = treatment,
      group    = interaction(zone, treatment))) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = PT_SIZE * 2) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = mean_y - sem_y, ymax = mean_y + sem_y),
      colour = "black", width = 0.10, linewidth = LW_ERR, show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = zone_colors, name = "Zone") +
    ggplot2::scale_shape_manual(values  = TREATMENT_SHAPES, name = "Treatment") +
    ggplot2::scale_linetype_manual(
      values = c("exercise choice" = "solid", "control" = "dashed"),
      name   = "Treatment") +
    ggplot2::scale_x_continuous(
      breaks = 1:3,
      labels = unname(LINE_TP_LABELS[as.character(1:3)])) +
    BASE_THEME +
    ggplot2::labs(title = title, x = "Timepoint", y = y_label, caption = n_txt)

  if (isTRUE(pct_scale))
    p <- p +
      ggplot2::scale_y_continuous(
        labels = function(x) round(x * 100),
        breaks = seq(0, 1, by = 0.2)) +
      ggplot2::expand_limits(y = c(0, 1))
  p
}

# ---- Indicator line plot: x = timepoint, mean ± SEM per treatment -----------
# color = treatment (pal), shape = treatment, linetype = treatment
# error bars and mean cross bars black
.make_indicator_line_plot <- function(df, y_col, y_label, title, n_txt = NULL) {
  if (is.null(df) || !y_col %in% names(df)) return(NULL)
  if (!"timepoint_f" %in% names(df)) return(NULL)
  smry <- df %>%
    dplyr::filter(is.finite(.data[[y_col]]), !is.na(timepoint_f)) %>%
    dplyr::group_by(treatment, timepoint_f) %>%
    dplyr::summarise(mean_y = mean(.data[[y_col]], na.rm = TRUE),
                     sem_y  = .sem(.data[[y_col]]),
                     .groups = "drop") %>%
    dplyr::mutate(tp_num = as.integer(as.character(timepoint_f)))
  if (nrow(smry) == 0) return(NULL)
  if (isTRUE(var(smry$mean_y, na.rm = TRUE) < 1e-10)) return(NULL)

  ggplot2::ggplot(smry, ggplot2::aes(
      x        = tp_num,
      y        = mean_y,
      colour   = treatment,
      shape    = treatment,
      linetype = treatment,
      group    = treatment)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = PT_SIZE * 2) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = mean_y - sem_y, ymax = mean_y + sem_y),
      colour = "black", width = 0.10, linewidth = LW_ERR, show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = TREATMENT_COLORS_PAL, name = "Treatment") +
    ggplot2::scale_shape_manual(values  = TREATMENT_SHAPES, name = "Treatment") +
    ggplot2::scale_linetype_manual(
      values = c("exercise choice" = "solid", "control" = "dashed"),
      name   = "Treatment") +
    ggplot2::scale_x_continuous(
      breaks = 1:3,
      labels = unname(LINE_TP_LABELS[as.character(1:3)])) +
    BASE_THEME +
    ggplot2::labs(title = title, x = "Timepoint", y = y_label, caption = n_txt)
}

# ---- PPTX writer -------------------------------------------------------------
.save_pptx <- function(plots, path, slide_w = 10, slide_h = 6) {
  prs <- officer::read_pptx()
  for (nm in names(plots)) {
    ts_msg("  Slide: ", nm)
    prs <- officer::add_slide(prs, layout = "Blank", master = "Office Theme")
    prs <- officer::ph_with(
      prs,
      value    = rvg::dml(ggobj = plots[[nm]]),
      location = officer::ph_location(left = 0, top = 0,
                                       width = slide_w, height = slide_h)
    )
    ts_msg("  Slide done: ", nm)
  }
  print(prs, target = path)
  ts_msg("Saved: ", basename(path))
}

# ---- Build ANOVA captions ----------------------------------------------------
.cap <- function(res, key) {
  val <- if (is.null(res)) NULL else res$anova_caps[[key]]
  if (is.null(val) || (length(val) == 1 && is.na(val))) return(NULL)
  val
}
# Combine all significant focal-term captions into a multi-line block.
# Returns NULL when no term is significant. Each significant term renders on
# its own line so plot captions can show every significant effect.
.cap_all_sig <- function(res) {
  if (is.null(res) || is.null(res$anova_caps)) return(NULL)
  caps <- vapply(res$anova_caps, function(v) {
    if (is.null(v) || (length(v) == 1L && is.na(v))) NA_character_ else as.character(v)
  }, character(1))
  caps <- caps[!is.na(caps) & nzchar(caps)]
  if (length(caps) == 0) return(NULL)
  paste(unname(caps), collapse = "\n")
}
# Compose the final plot footer caption: sample-size line first, ANOVA
# significant terms below (per spec: ANOVA caption directly below sample size).
.compose_caption <- function(n_txt, anova_txt) {
  parts <- c(n_txt, anova_txt)
  parts <- parts[!is.null(parts)]
  parts <- parts[!is.na(parts) & nzchar(parts)]
  if (length(parts) == 0) NULL else paste(parts, collapse = "\n")
}
.cld <- function(res, term_re) {
  if (is.null(res) || is.null(res$posthoc)) return(NULL)
  ph <- res$posthoc[grepl(term_re, names(res$posthoc), ignore.case=TRUE, perl=TRUE)]
  if (length(ph) == 0) return(NULL)
  ph[[1]]$cld
}
.is_sig <- function(res, key) {
  if (is.null(res) || is.null(res$anova_caps)) return(FALSE)
  val <- res$anova_caps[[key]]
  !is.null(val) && !(length(val) == 1L && is.na(val))
}

# CLD letters differ across groups (i.e. some pairwise contrasts were sig).
.cld_letters_differ <- function(cld_df) {
  if (is.null(cld_df) || !".group" %in% names(cld_df)) return(FALSE)
  u <- unique(trimws(as.character(cld_df$.group[!is.na(cld_df$.group)])))
  length(u) > 1L
}

# Build the formatted statistic string for a term regardless of p-value.
# Mirrors anova_caption_str's format selection (F-KR/F-SW vs Wald-chisq) but
# does not gate on significance — used to disclose a NS omnibus when CLD
# letters reveal post-hoc differences.
.term_str_anyp <- function(av, term_pattern, stat_type = NULL) {
  if (is.null(av)) return(NA_character_)
  row <- av[grepl(term_pattern, av$term, ignore.case = TRUE, perl = TRUE), ]
  if (nrow(row) == 0) return(NA_character_)
  row <- row[1, ]
  if (is.null(row$p_value) || length(row$p_value) == 0 || is.na(row$p_value))
    return(NA_character_)
  stat_type <- stat_type %||% attr(av, "stat_type") %||% "Wald-chisq"
  num_fmt   <- fmt_F
  stat_val  <- fmt_F(as.numeric(row$chisq))
  is_f_test <- isTRUE(grepl("^F[-_]", stat_type))

  if (is_f_test) {
    df_num <- as.integer(row$df)
    df_den_raw <- if (!is.null(row$df_denom)) row$df_denom else NA_real_
    df_den <- if (length(df_den_raw) > 0 && !is.na(df_den_raw) &&
                  is.finite(as.numeric(df_den_raw))) num_fmt(df_den_raw) else NA_character_
    eta_str <- ""
    if (!is.na(df_den)) {
      .f <- suppressWarnings(as.numeric(row$chisq))
      .d2 <- suppressWarnings(as.numeric(df_den_raw))
      if (is.finite(.f) && is.finite(.d2) && is.finite(df_num)) {
        .e <- (.f * df_num) / (.f * df_num + .d2)
        eta_str <- sprintf(", &eta;<sup>2</sup><sub>p</sub> = %s", fmt_es3(.e))
      }
    }
    if (!is.na(df_den))
      sprintf("F<sub>%d, %s</sub> = %s, %s%s", df_num, df_den, stat_val,
              fmt_p(row$p_value), eta_str)
    else
      sprintf("F<sub>%d</sub> = %s, %s", df_num, stat_val, fmt_p(row$p_value))
  } else {
    paste0("χ²<sub>", as.integer(row$df), "</sub> = ",
           stat_val, ", ", fmt_p(row$p_value))
  }
}

# Pretty-print a term name like "treatment:zone" → "Treatment × Zone".
.pretty_term <- function(term_str) {
  one <- function(x) {
    x <- sub("_f$", "", x)
    paste0(toupper(substring(x, 1, 1)), substring(x, 2))
  }
  if (grepl(":", term_str, fixed = TRUE))
    paste(sapply(strsplit(term_str, ":")[[1]], one), collapse = " × ")
  else one(term_str)
}

# Caption builder for occupancy/zone plots:
#   - includes every significant focal term (via .cap_all_sig), and
#   - appends a NS-but-CLD-differs note when the focal interaction is not
#     significant yet the CLD letters differ across groups (some pairwise
#     contrasts were significant). The note discloses the omnibus statistic
#     and points the reader to the letters.
.cap_zone <- function(res, term_key, cld_df) {
  base <- .cap_all_sig(res)
  if (is.null(res) || is.null(res$anova) ||
      is.null(term_key) || !.cld_letters_differ(cld_df)) return(base)
  is_term_sig <- .is_sig(res, term_key)
  if (isTRUE(is_term_sig)) return(base)
  # Map the focal_terms key (e.g. "Treatment:Zone:Timepoint") to a regex
  # against the anova table's term column. Underscores / colons separate parts.
  pattern <- paste0(strsplit(term_key, ":", fixed = TRUE)[[1]],
                    collapse = ".*")
  stat_str <- .term_str_anyp(res$anova, pattern, res$stat_type)
  if (is.na(stat_str)) return(base)
  note <- paste0(.pretty_term(term_key),
                 ": ANOVA n.s. (", stat_str,
                 "); post-hoc CLD letters differ — see panel.")
  .compose_caption(base, note)
}

# ---- Log-ratio plot builder --------------------------------------------------
# Like .make_std_plot but: (a) adds a dashed reference line at y = 0 to mark
# equal use of both zones, and (b) allows negative y-axis breaks (log-ratios
# span below zero, unlike all other indicators).
.make_logratio_plot <- function(df, y_col, y_label, title,
                                 caption_txt = NULL, is_sig = FALSE,
                                 facet_tp = FALSE, n_txt = NULL) {
  if (is.null(df) || !y_col %in% names(df)) return(NULL)
  .vals <- df[[y_col]][is.finite(df[[y_col]])]
  if (length(.vals) == 0) return(NULL)
  groups_smry <- if (facet_tp) c("treatment", "timepoint_f") else "treatment"
  smry <- .smry_std(df, y_col, groups_smry)
  .trt_lvls  <- levels(factor(df$treatment))
  .n_trt     <- length(.trt_lvls)
  .trt_map   <- stats::setNames(seq_len(.n_trt), .trt_lvls)
  df_plot    <- df[is.finite(df[[y_col]]), ]
  df_plot$.x_pos <- unname(.trt_map[as.character(df_plot$treatment)])
  smry$.x_pos    <- unname(.trt_map[as.character(smry$treatment)])
  POS_JITTER <- ggplot2::position_jitter(width = JITTER_W, height = 0)
  p <- ggplot2::ggplot(df_plot, ggplot2::aes(
      x = .x_pos, y = .data[[y_col]],
      colour = treatment, shape = treatment)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50",
                        linewidth = 0.6) +
    ggplot2::geom_point(position = POS_JITTER, alpha = PT_ALPHA, size = PT_SIZE) +
    ggplot2::geom_crossbar(
      data = smry,
      ggplot2::aes(x = .x_pos, y = mean_y, ymin = mean_y, ymax = mean_y),
      colour = "black", width = MEAN_W, linewidth = LW_MEAN, show.legend = FALSE) +
    ggplot2::geom_errorbar(
      data = smry,
      ggplot2::aes(x = .x_pos, y = mean_y,
                   ymin = mean_y - sem_y, ymax = mean_y + sem_y),
      colour = "black", width = ERR_W, linewidth = LW_ERR, show.legend = FALSE) +
    ggplot2::scale_x_continuous(
      breaks = seq_len(.n_trt),
      labels = stringr::str_to_title(.trt_lvls),
      limits = c(0.5, .n_trt + 0.5),
      expand = c(0, 0)) +
    ggplot2::scale_colour_manual(values = TREATMENT_COLORS_PAL, name = "Treatment") +
    ggplot2::scale_shape_manual(values  = TREATMENT_SHAPES,     name = "Treatment") +
    ggplot2::scale_y_continuous(breaks = pretty) +
    BASE_THEME +
    ggplot2::labs(title = title, x = "Treatment", y = y_label,
                  caption = .compose_caption(n_txt, caption_txt))
  if (isTRUE(is_sig) && .n_trt == 2L) {
    .y_max  <- max(.vals, na.rm = TRUE)
    .y_min  <- min(.vals, na.rm = TRUE)
    .y_rng  <- max(.y_max - .y_min, 1e-10)
    .bar_y  <- .y_max + 0.10 * .y_rng
    .cap_h  <- 0.025 * .y_rng
    .star_y <- .bar_y + 0.04 * .y_rng
    .y_top  <- .star_y + 0.06 * .y_rng
    .y_bot  <- .y_min  - 0.02 * .y_rng
    .x1 <- 1L; .x2 <- .n_trt; .xm <- (.x1 + .x2) / 2
    p <- p +
      ggplot2::expand_limits(y = c(.y_bot, .y_top)) +
      ggplot2::geom_segment(
        data = data.frame(x = .x1, xend = .x2, y = .bar_y, yend = .bar_y),
        ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
        colour = "black", linewidth = 0.7, inherit.aes = FALSE) +
      ggplot2::geom_segment(
        data = data.frame(x    = c(.x1, .x2), xend = c(.x1, .x2),
                          y    = rep(.bar_y - .cap_h, 2),
                          yend = rep(.bar_y, 2)),
        ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
        colour = "black", linewidth = 0.7, inherit.aes = FALSE) +
      ggplot2::geom_text(
        data = data.frame(x = .xm, y = .star_y),
        ggplot2::aes(x = x, y = y), label = "*",
        size = CLD_SIZE + 2, colour = "black", fontface = "bold",
        inherit.aes = FALSE, show.legend = FALSE)
  }
  if (facet_tp)
    p <- p + ggplot2::facet_wrap(
      ~timepoint_f, ncol = 1,
      labeller = ggplot2::labeller(timepoint_f = TIMEPOINT_LABELS))
  p
}

# ---- Build the 6 plots -------------------------------------------------------

p_zone_agg <- tryCatch(
  .make_zone_plot(
    df_long     = df_main_long_agg,
    y_col       = "prop_time",
    zone_colors = color_map_broad,
    y_label     = "Fish occupancy zone %",
    title       = "Main-zone occupancy by treatment",
    caption_txt = .cap_zone(res_zone_main_agg, "Treatment:Zone",
                            .cld(res_zone_main_agg, "treatment.*zone|zone.*treatment")),
    cld_df      = .cld(res_zone_main_agg, "treatment.*zone|zone.*treatment"),
    is_sig      = .is_sig(res_zone_main_agg, "Treatment:Zone"),
    facet_tp    = FALSE,
    pct_scale   = TRUE,
    n_txt       = .n_cap(df_main_long_agg, "agg")
  ), error = function(e) { warning("p_zone_agg failed: ", e$message); NULL }
)

p_zone_tp <- tryCatch(
  .make_zone_plot(
    df_long     = df_main_long,
    y_col       = "prop_time",
    zone_colors = color_map_broad,
    y_label     = "Fish occupancy zone %",
    title       = "Main-zone occupancy by treatment \u00d7 timepoint",
    caption_txt = .cap_zone(res_zone_main, "Treatment:Zone:Timepoint",
                            .cld(res_zone_main, "treatment.*zone.*timepoint")),
    cld_df      = .cld(res_zone_main, "treatment.*zone.*timepoint"),
    is_sig      = .is_sig(res_zone_main, "Treatment:Zone:Timepoint"),
    facet_tp    = TRUE,
    pct_scale   = TRUE,
    n_txt       = .n_cap(df_main_long, "tp")
  ), error = function(e) { warning("p_zone_tp failed: ", e$message); NULL }
)

p_subzone_tp <- tryCatch(
  .make_zone_plot(
    df_long     = df_sec_long,
    y_col       = "prop_time",
    zone_colors = color_map[c("high", "medium", "low", "calm")],
    y_label     = "Fish occupancy zone %",
    title       = "Sub-zone occupancy by treatment \u00d7 timepoint",
    caption_txt = .cap_zone(res_zone_sec, "Treatment:Zone:Timepoint",
                            .cld(res_zone_sec, "treatment.*zone.*timepoint")),
    cld_df      = .cld(res_zone_sec, "treatment.*zone.*timepoint"),
    is_sig      = .is_sig(res_zone_sec, "Treatment:Zone:Timepoint"),
    facet_tp    = TRUE,
    pct_scale   = TRUE,
    n_txt       = .n_cap(df_sec_long, "tp")
  ), error = function(e) { warning("p_subzone_tp failed: ", e$message); NULL }
)

p_subzone_agg <- tryCatch(
  .make_zone_plot(
    df_long     = df_sec_long_agg,
    y_col       = "prop_time",
    zone_colors = color_map[c("high", "medium", "low", "calm")],
    y_label     = "Fish occupancy zone %",
    title       = "Sub-zone occupancy by treatment",
    caption_txt = .cap_zone(res_zone_sec_agg, "Treatment:Zone",
                            .cld(res_zone_sec_agg, "treatment.*zone|zone.*treatment")),
    cld_df      = .cld(res_zone_sec_agg, "treatment.*zone|zone.*treatment"),
    is_sig      = .is_sig(res_zone_sec_agg, "Treatment:Zone"),
    facet_tp    = FALSE,
    pct_scale   = TRUE,
    n_txt       = .n_cap(df_sec_long_agg, "agg")
  ), error = function(e) { warning("p_subzone_agg failed: ", e$message); NULL }
)

p_switch_agg <- tryCatch(
  .make_std_plot(
    df          = df_agg,
    y_col       = "switches_per_session",
    y_label     = "Transitions between main zones",
    title       = "",
    caption_txt = .cap_all_sig(res_switch_agg),
    cld_df      = .cld(res_switch_agg, "^treatment$"),
    is_sig      = .is_sig(res_switch_agg, "Treatment"),
    facet_tp    = FALSE,
    n_txt       = .n_cap(df_agg, "agg")
  ), error = function(e) { warning("p_switch_agg failed: ", e$message); NULL }
)

p_switch_tp <- tryCatch(
  .make_std_plot(
    df          = df,
    y_col       = "switches_per_session",
    y_label     = "Transitions between main zones",
    title       = "Switches/session by treatment \u00d7 timepoint",
    caption_txt = .cap_all_sig(res_switch_tp),
    cld_df      = .cld(res_switch_tp, "treatment.*timepoint|timepoint.*treatment"),
    is_sig      = .is_sig(res_switch_tp, "Treatment:Timepoint"),
    facet_tp    = TRUE,
    n_txt       = .n_cap(df, "tp")
  ), error = function(e) { warning("p_switch_tp failed: ", e$message); NULL }
)

p_moving_tp <- tryCatch(
  .make_std_plot(
    df          = df,
    y_col       = "prop_active",
    y_label     = "Proportion active (> 1 BL s\u207b\u00b9)",
    title       = "School active fraction by treatment \u00d7 timepoint",
    caption_txt = .cap_all_sig(res_active_tp),
    cld_df      = .cld(res_active_tp, "treatment.*timepoint|timepoint.*treatment"),
    is_sig      = .is_sig(res_active_tp, "Treatment:Timepoint"),
    facet_tp    = TRUE,
    n_txt       = .n_cap(df, "tp")
  ), error = function(e) { warning("p_moving_tp failed: ", e$message); NULL }
)

p_moving_agg <- tryCatch(
  .make_std_plot(
    df          = df_agg,
    y_col       = "prop_active",
    y_label     = "Proportion active (> 1 BL s\u207b\u00b9)",
    title       = "School active fraction by treatment",
    caption_txt = .cap_all_sig(res_active_agg),
    cld_df      = .cld(res_active_agg, "^treatment$"),
    is_sig      = .is_sig(res_active_agg, "Treatment"),
    facet_tp    = FALSE,
    n_txt       = .n_cap(df_agg, "agg")
  ), error = function(e) { warning("p_moving_agg failed: ", e$message); NULL }
)

p_flux_tp <- tryCatch(
  .make_std_plot(
    df          = df,
    y_col       = "zone_flux_per_session",
    y_label     = "Zone flux per session",
    title       = paste0("Zone flux by treatment × timepoint"),
    caption_txt = .cap_all_sig(res_flux_tp),
    cld_df      = .cld(res_flux_tp, "treatment.*timepoint|timepoint.*treatment"),
    is_sig      = .is_sig(res_flux_tp, "Treatment:Timepoint"),
    facet_tp    = TRUE,
    n_txt       = .n_cap(df, "tp")
  ), error = function(e) { warning("p_flux_tp failed: ", e$message); NULL }
)

p_flux_agg <- tryCatch(
  .make_std_plot(
    df          = df_agg,
    y_col       = "zone_flux_per_session",
    y_label     = "Zone flux per session",
    title       = "Zone flux by treatment",
    caption_txt = .cap_all_sig(res_flux_agg),
    cld_df      = .cld(res_flux_agg, "^treatment$"),
    is_sig      = .is_sig(res_flux_agg, "Treatment"),
    facet_tp    = FALSE,
    n_txt       = .n_cap(df_agg, "agg")
  ), error = function(e) { warning("p_flux_agg failed: ", e$message); NULL }
)

# ---- Log-ratio preference plots (A1c + A2) -----------------------------------
# Main-zone preference score: logit(prop_flow/prop_calm) — aggregated
p_logitflow_agg <- tryCatch(
  .make_logratio_plot(
    df          = df_main_wide_agg,
    y_col       = "logit_flow",
    y_label     = "Flow preference score (log-ratio)\n↑ flow   |   0 = equal   |   ↓ calm",
    title       = "Main-zone choice by treatment",
    caption_txt = .cap_all_sig(res_zone_flow_logit_agg),
    is_sig      = .is_sig(res_zone_flow_logit_agg, "Treatment"),
    facet_tp    = FALSE,
    n_txt       = .n_cap(df_main_wide_agg, "agg")
  ), error = function(e) { warning("p_logitflow_agg failed: ", e$message); NULL }
)

# Main-zone preference score — timepoint-resolved
p_logitflow_tp <- tryCatch(
  .make_logratio_plot(
    df          = df_main_wide,
    y_col       = "logit_flow",
    y_label     = "Flow preference score (log-ratio)\n↑ flow   |   0 = equal   |   ↓ calm",
    title       = "Main-zone choice by treatment × timepoint",
    caption_txt = .cap_all_sig(res_zone_main),
    is_sig      = .is_sig(res_zone_main, "Treatment:Timepoint"),
    facet_tp    = TRUE,
    n_txt       = .n_cap(df_main_wide, "tp")
  ), error = function(e) { warning("p_logitflow_tp failed: ", e$message); NULL }
)

# Sub-zone log-ratio plots: log(subzone / calm) — aggregated (3 zones)
.lr_agg_specs <- list(
  list(y = "lr_high",   zone = "High-flow",   res = res_zone_sec_high_agg),
  list(y = "lr_medium", zone = "Medium-flow", res = res_zone_sec_med_agg),
  list(y = "lr_low",    zone = "Low-flow",    res = res_zone_sec_low_agg)
)
.lr_plots_agg <- lapply(.lr_agg_specs, function(s) {
  tryCatch(
    .make_logratio_plot(
      df          = df_sec_wide_agg,
      y_col       = s$y,
      y_label     = paste0(s$zone, " preference (log-ratio vs calm)\n",
                           "↑ subzone   |   0 = equal   |   ↓ calm"),
      title       = paste0(s$zone, " sub-zone choice by treatment"),
      caption_txt = .cap_all_sig(s$res),
      is_sig      = .is_sig(s$res, "Treatment"),
      facet_tp    = FALSE,
      n_txt       = .n_cap(df_sec_wide_agg, "agg")
    ), error = function(e) { warning(s$y, "_agg lr plot failed: ", e$message); NULL }
  )
})
names(.lr_plots_agg) <- paste0(c("lr_high", "lr_medium", "lr_low"), "_agg")

# Sub-zone log-ratio plots — timepoint-resolved
.lr_tp_specs <- list(
  list(y = "lr_high",   zone = "High-flow",   res = res_zone_sec_high_tp),
  list(y = "lr_medium", zone = "Medium-flow", res = res_zone_sec_med_tp),
  list(y = "lr_low",    zone = "Low-flow",    res = res_zone_sec_low_tp)
)
.lr_plots_tp <- lapply(.lr_tp_specs, function(s) {
  tryCatch(
    .make_logratio_plot(
      df          = df_sec_wide,
      y_col       = s$y,
      y_label     = paste0(s$zone, " preference (log-ratio vs calm)\n",
                           "↑ subzone   |   0 = equal   |   ↓ calm"),
      title       = paste0(s$zone, " sub-zone choice by treatment × timepoint"),
      caption_txt = .cap_all_sig(s$res),
      is_sig      = .is_sig(s$res, "Treatment:Timepoint"),
      facet_tp    = TRUE,
      n_txt       = .n_cap(df_sec_wide, "tp")
    ), error = function(e) { warning(s$y, "_tp lr plot failed: ", e$message); NULL }
  )
})
names(.lr_plots_tp) <- paste0(c("lr_high", "lr_medium", "lr_low"), "_tp")

# ---- Line plots: occupancy × timepoint --------------------------------------
p_mainzone_line <- tryCatch(
  .make_zone_line_plot(
    df_long     = df_main_long,
    y_col       = "prop_time",
    zone_colors = color_map_broad,
    y_label     = "Fish occupancy zone %",
    title       = "Main-zone occupancy across timepoints",
    pct_scale   = TRUE,
    n_txt       = .n_cap(df_main_long, "tp")
  ), error = function(e) { warning("p_mainzone_line failed: ", e$message); NULL }
)

p_seczone_line <- tryCatch(
  .make_zone_line_plot(
    df_long     = df_sec_long,
    y_col       = "prop_time",
    zone_colors = color_map[c("high","medium","low","calm")],
    y_label     = "Fish occupancy zone %",
    title       = "Sub-zone occupancy across timepoints",
    pct_scale   = TRUE,
    n_txt       = .n_cap(df_sec_long, "tp")
  ), error = function(e) { warning("p_seczone_line failed: ", e$message); NULL }
)

# ---- Indicator line plots × timepoint (pal colors) --------------------------
.ind_tp_specs <- list(
  list(y = "switches_per_session",  lbl = "Transitions between main zones",
       ttl = "Zone switches across timepoints", src = "df"),
  list(y = "zone_flux_per_session", lbl = "Zone flux per session",
       ttl = "Zone flux across timepoints", src = "df"),
  list(y = "mean_nnd_cm",     lbl = "Mean NND (cm)",
       ttl = "Nearest-neighbour distance across timepoints", src = "gd"),
  list(y = "mean_iid_cm",    lbl = "Mean IID (cm)",
       ttl = "Inter-individual distance across timepoints", src = "gd"),
  list(y = "mean_hull_area_cm2",  lbl = "School area (cm²)",
       ttl = "Group spread across timepoints", src = "gd"),
  list(y = "mean_centroid_spd_cm", lbl = "School speed (cm/s)",
       ttl = "School speed across timepoints", src = "gd")
)

.ind_line_plots <- lapply(.ind_tp_specs, function(s) {
  .dat <- if (s$src == "gd") {
    if (!is.null(df_gd) && s$y %in% names(df_gd)) df_gd else NULL
  } else {
    if (s$y %in% names(df)) df else NULL
  }
  tryCatch(
    .make_indicator_line_plot(.dat, s$y, s$lbl, s$ttl,
                               n_txt = .n_cap(.dat, "tp")),
    error = function(e) { warning("indicator line plot '", s$y, "' failed: ", e$message); NULL }
  )
})
names(.ind_line_plots) <- paste0(sapply(.ind_tp_specs, `[[`, "y"), "_line_tp")

# ---- Group-dynamics indicator scatter plots ---------------------------------
.gd_specs <- list(
  list(y = "mean_nnd_cm",          lbl = "Mean NND (cm)",
       ttl_agg = "NND by treatment",              ttl_tp = "NND by treatment × timepoint",
       res_agg = "res_nnd_agg",  res_tp = "res_nnd_tp"),
  list(y = "mean_iid_cm",          lbl = "Mean IID (cm)",
       ttl_agg = "IID by treatment",              ttl_tp = "IID by treatment × timepoint",
       res_agg = "res_iid_agg",  res_tp = "res_iid_tp"),
  list(y = "mean_hull_area_cm2",   lbl = "School area (cm²)",
       ttl_agg = "School area by treatment",       ttl_tp = "School area by treatment × timepoint",
       res_agg = "res_hull_agg", res_tp = "res_hull_tp"),
  list(y = "mean_centroid_spd_cm", lbl = "School speed (cm/s)",
       ttl_agg = "School speed by treatment",      ttl_tp = "School speed by treatment × timepoint",
       res_agg = "res_cspd_agg", res_tp = "res_cspd_tp")
)

.gd_scatter_agg <- lapply(.gd_specs, function(s) {
  r   <- tryCatch(get(s$res_agg), error = function(e) NULL)
  dat <- if (!is.null(df_gd_agg) && s$y %in% names(df_gd_agg)) df_gd_agg else NULL
  if (is.null(dat)) return(NULL)
  tryCatch(.make_std_plot(
    df = dat, y_col = s$y, y_label = s$lbl, title = s$ttl_agg,
    caption_txt = .cap_all_sig(r),
    cld_df      = .cld(r, "^treatment$"),
    is_sig      = .is_sig(r, "Treatment"),
    facet_tp    = FALSE,
    n_txt       = .n_cap(dat, "agg")
  ), error = function(e) { warning(s$y, "_agg plot failed: ", e$message); NULL })
})
names(.gd_scatter_agg) <- paste0(sapply(.gd_specs, `[[`, "y"), "_agg_scatter")

.gd_scatter_tp <- lapply(.gd_specs, function(s) {
  r   <- tryCatch(get(s$res_tp), error = function(e) NULL)
  dat <- if (!is.null(df_gd) && s$y %in% names(df_gd)) df_gd else NULL
  if (is.null(dat)) return(NULL)
  tryCatch(.make_std_plot(
    df = dat, y_col = s$y, y_label = s$lbl, title = s$ttl_tp,
    caption_txt = .cap_all_sig(r),
    cld_df      = .cld(r, "treatment.*timepoint"),
    is_sig      = .is_sig(r, "Treatment:Timepoint"),
    facet_tp    = TRUE,
    n_txt       = .n_cap(dat, "tp")
  ), error = function(e) { warning(s$y, "_tp plot failed: ", e$message); NULL })
})
names(.gd_scatter_tp) <- paste0(sapply(.gd_specs, `[[`, "y"), "_tp_scatter")

# ---- Save PPTX ---------------------------------------------------------------
plots_agg <- Filter(Negate(is.null), c(
  list(zone_main_aggregated    = p_zone_agg,
       zone_sec_aggregated     = p_subzone_agg,
       logratio_flow_agg       = p_logitflow_agg,
       switches_aggregated     = p_switch_agg,
       flux_aggregated         = p_flux_agg,
       active_aggregated       = p_moving_agg),
  .lr_plots_agg,
  .gd_scatter_agg
))

plots_tp  <- Filter(Negate(is.null), c(
  list(zone_main_timepoint     = p_zone_tp,
       zone_sec_timepoint      = p_subzone_tp,
       logratio_flow_tp        = p_logitflow_tp,
       switches_timepoint      = p_switch_tp,
       flux_timepoint          = p_flux_tp,
       active_timepoint        = p_moving_tp,
       zone_main_line_tp       = p_mainzone_line,
       zone_sec_line_tp        = p_seczone_line),
  .lr_plots_tp,
  .gd_scatter_tp,
  .ind_line_plots
))

if (length(plots_agg) > 0)
  tryCatch(
    .save_pptx(plots_agg,
               file.path(.sg_dir, "skinny_graphs_aggregated.pptx"),
               slide_w = 10, slide_h = 6),
    error = function(e)
      warning("skinny_graphs_aggregated.pptx failed: ", e$message, call. = FALSE)
  )

if (length(plots_tp) > 0)
  tryCatch(
    .save_pptx(plots_tp,
               file.path(.sg_dir, "skinny_graphs_timepoint.pptx"),
               slide_w = 10, slide_h = 9),
    error = function(e)
      warning("skinny_graphs_timepoint.pptx failed: ", e$message, call. = FALSE)
  )

# ---- PNG previews ------------------------------------------------------------
all_sg <- c(plots_agg, plots_tp)
invisible(lapply(names(all_sg), function(nm) {
  h <- if (grepl("timepoint|tp", nm)) 9 else 6
  tryCatch(
    ggplot2::ggsave(file.path(.sg_png, paste0(nm, ".png")),
                    all_sg[[nm]], width = 10, height = h, dpi = 150),
    error = function(e) NULL
  )
}))
ts_msg("PNG previews: ", .sg_png)

# ---- sign_graph: export only plots with significant results ------------------
.sg_sign_dir <- file.path(STEP5_OUT, "sign_graph")
.sg_sign_png <- file.path(.sg_sign_dir, "png")
dir.create(.sg_sign_png, showWarnings = FALSE, recursive = TRUE)

.sig_check <- function(res, cld_pattern) {
  isTRUE(!is.null(.cap_all_sig(res))) ||
    isTRUE(.cld_letters_differ(.cld(res, cld_pattern)))
}
.has_sig_caption <- function(p) {
  if (is.null(p)) return(FALSE)
  cap <- tryCatch(p$labels$caption, error = function(e) NULL)
  if (is.null(cap)) return(FALSE)
  grepl("F[(]|[xX]2[(]|CLD|ANOVA|p\\s*[<=]\\s*0[.]", as.character(cap))
}

.sig_flags <- c(
  # A1_agg: kept for graph annotation even though not used for inference
  zone_main_aggregated = .sig_check(res_zone_main_agg, "treatment.*zone|zone.*treatment"),
  # A1 (logit): treatment main effect or treatment:timepoint interaction
  zone_main_timepoint  = .sig_check(res_zone_main,     "treatment.*timepoint|timepoint.*treatment"),
  # A2 log-ratio models: treatment main effect (high-zone as representative)
  zone_sec_aggregated  = .sig_check(res_zone_sec_high_agg, "^treatment$"),
  zone_sec_timepoint   = .sig_check(res_zone_sec_high_tp,  "treatment.*timepoint|timepoint.*treatment"),
  switches_aggregated  = .sig_check(res_switch_agg,    "^treatment$"),
  switches_timepoint   = .sig_check(res_switch_tp,     "treatment.*timepoint|timepoint.*treatment"),
  flux_aggregated      = .sig_check(res_flux_agg,      "^treatment$"),
  flux_timepoint       = .sig_check(res_flux_tp,       "treatment.*timepoint|timepoint.*treatment"),
  active_aggregated    = .sig_check(res_active_agg,    "^treatment$"),
  active_timepoint     = .sig_check(res_active_tp,     "treatment.*timepoint|timepoint.*treatment")
)

for (.s in .gd_specs) {
  .r_agg <- tryCatch(get(.s$res_agg), error = function(e) NULL)
  .r_tp  <- tryCatch(get(.s$res_tp),  error = function(e) NULL)
  .sig_flags[paste0(.s$y, "_agg_scatter")] <- .sig_check(.r_agg, "^treatment$")
  .sig_flags[paste0(.s$y, "_tp_scatter")]  <- .sig_check(.r_tp,
    "treatment.*timepoint|timepoint.*treatment")
}

.sig_names_known  <- names(.sig_flags)[as.logical(.sig_flags)]
.sig_names_extra  <- Filter(function(nm) .has_sig_caption(all_sg[[nm]]),
                            setdiff(names(all_sg), names(.sig_flags)))
.all_sig_names    <- unique(c(.sig_names_known, .sig_names_extra))

plots_sig_agg <- Filter(Negate(is.null),
  all_sg[intersect(.all_sig_names, names(plots_agg))])
plots_sig_tp  <- Filter(Negate(is.null),
  all_sg[intersect(.all_sig_names, names(plots_tp))])

ts_msg("sign_graph: ", length(plots_sig_agg), " agg + ",
       length(plots_sig_tp), " tp significant plot(s)")

if (length(plots_sig_agg) > 0)
  tryCatch(
    .save_pptx(plots_sig_agg,
               file.path(.sg_sign_dir, "sign_graphs_aggregated.pptx"),
               slide_w = 10, slide_h = 6),
    error = function(e)
      warning("sign_graphs_aggregated.pptx failed: ", e$message, call. = FALSE)
  )

if (length(plots_sig_tp) > 0)
  tryCatch(
    .save_pptx(plots_sig_tp,
               file.path(.sg_sign_dir, "sign_graphs_timepoint.pptx"),
               slide_w = 10, slide_h = 9),
    error = function(e)
      warning("sign_graphs_timepoint.pptx failed: ", e$message, call. = FALSE)
  )

invisible(lapply(.all_sig_names, function(nm) {
  p <- all_sg[[nm]]
  if (is.null(p)) return(invisible(NULL))
  h <- if (nm %in% names(plots_agg)) 6 else 9
  tryCatch(
    ggplot2::ggsave(file.path(.sg_sign_png, paste0(nm, ".png")),
                    p, width = 10, height = h, dpi = 150),
    error = function(e) NULL
  )
}))
ts_msg("sign_graph written: ", .sg_sign_dir)

# ---- sign_graph/stacked: combined layout plots -------------------------------
suppressPackageStartupMessages(library(patchwork))

.sg_stacked_dir <- file.path(.sg_sign_dir, "stacked")
dir.create(.sg_stacked_dir, showWarnings = FALSE, recursive = TRUE)

# Helper: drop title, add minimal white frame with margin
.no_title_frame <- function(p) {
  if (is.null(p)) return(NULL)
  p + ggplot2::theme(
    plot.title      = ggplot2::element_blank(),
    plot.margin     = ggplot2::margin(10, 10, 10, 10),
    plot.background = ggplot2::element_rect(fill = "white", colour = "grey80",
                                             linewidth = 0.5)
  )
}

# Helper: truncate .group labels to n chars in all geom_text layers of a plot
.trunc_cld_in_plot <- function(p, n = 2L) {
  if (is.null(p)) return(NULL)
  for (i in seq_along(p$layers)) {
    lyr <- p$layers[[i]]
    if (inherits(lyr$geom, "GeomText") && is.data.frame(lyr$data) &&
        ".group" %in% names(lyr$data)) {
      p$layers[[i]]$data$.group <-
        substr(trimws(as.character(lyr$data$.group)), 1L, n)
    }
  }
  p
}

# ---- Combo 1: zone_main / zone_sec stacked (A above B) ----------------------
.c1_A <- .trunc_cld_in_plot(.no_title_frame(p_zone_agg))
.c1_B <- .trunc_cld_in_plot(.no_title_frame(p_subzone_agg))

if (!is.null(.c1_A) && !is.null(.c1_B)) {
  .combo1 <- (.c1_A / .c1_B) +
    patchwork::plot_annotation(tag_levels = list(c("A", "B"))) &
    ggplot2::theme(
      plot.tag          = ggplot2::element_text(face = "bold", size = 18),
      plot.tag.position = c(0.02, 0.98)
    )
  tryCatch(
    ggplot2::ggsave(file.path(.sg_stacked_dir, "combo_AB_zone_stacked.png"),
                    .combo1, width = 10, height = 11, dpi = 150),
    error = function(e) warning("combo1 ggsave failed: ", e$message)
  )
  ts_msg("stacked combo1 (A/B zone) saved")
}

# ---- Combo 2: NND/IID/Hull row with 3× bigger CLD + truncated labels -------
.trunc2 <- function(cld_df) {
  if (is.null(cld_df) || !".group" %in% names(cld_df)) return(cld_df)
  cld_df$.group <- substr(trimws(as.character(cld_df$.group)), 1, 2)
  cld_df
}

.OLD_CLD_SIZE <- CLD_SIZE
CLD_SIZE <- .OLD_CLD_SIZE * 3  # 3× bigger CLD letters and asterisks

.mk_gd_big <- function(res, y_col, y_label) {
  dat <- if (!is.null(df_gd_agg) && y_col %in% names(df_gd_agg)) df_gd_agg else NULL
  if (is.null(dat)) return(NULL)
  tryCatch(
    .no_title_frame(.make_std_plot(
      df          = dat,
      y_col       = y_col,
      y_label     = y_label,
      title       = "",
      caption_txt = .cap_all_sig(res),
      cld_df      = .trunc2(.cld(res, "^treatment$")),
      is_sig      = .is_sig(res, "Treatment"),
      facet_tp    = FALSE,
      n_txt       = .n_cap(dat, "agg")
    )),
    error = function(e) { warning("stacked gd '", y_col, "' failed: ", e$message); NULL }
  )
}

.p2_nnd  <- .trunc_cld_in_plot(.mk_gd_big(res_nnd_agg,  "mean_nnd_cm",        "Mean NND (cm)"))
.p2_iid  <- .trunc_cld_in_plot(.mk_gd_big(res_iid_agg,  "mean_iid_cm",        "Mean IID (cm)"))
.p2_hull <- .trunc_cld_in_plot(.mk_gd_big(res_hull_agg, "mean_hull_area_cm2", "School area (cm²)"))

CLD_SIZE <- .OLD_CLD_SIZE  # restore global

# Was a 2x2 (nnd|pol)/(hull|iid). The alignment panel was removed 2026-08-16,
# leaving three collective indicators, so this is now a single row: A=NND,
# B=IID, C=Hull.
if (!any(sapply(list(.p2_nnd, .p2_iid, .p2_hull), is.null))) {
  .combo2 <- (.p2_nnd | .p2_iid | .p2_hull) +
    patchwork::plot_annotation(tag_levels = list(c("A", "B", "C"))) &
    ggplot2::theme(
      plot.tag          = ggplot2::element_text(face = "bold", size = 18),
      plot.tag.position = c(0.02, 0.98)
    )
  tryCatch(
    ggplot2::ggsave(file.path(.sg_stacked_dir, "combo_ABCD_gd_grid.png"),
                    .combo2, width = 14, height = 10, dpi = 150),
    error = function(e) warning("combo2 ggsave failed: ", e$message)
  )
  ts_msg("stacked combo2 (A-D 2x2 gd grid) saved")
} else {
  ts_msg("stacked combo2 skipped — one or more GD plots NULL")
}

# ---- Combo 3: IID/NND/Hull row — hull legend only, pts 50% bigger ----------
# All other changes (3× CLD, truncated labels, no-title frame) kept from combo2.
.OLD_CLD_SIZE3 <- CLD_SIZE
.OLD_PT_SIZE3  <- PT_SIZE
CLD_SIZE <- .OLD_CLD_SIZE3 * 3   # 3× bigger CLD as in combo2
PT_SIZE  <- .OLD_PT_SIZE3  * 1.5 # 50% bigger scatter points

.mk_gd_big3 <- function(res, y_col, y_label, keep_legend = FALSE) {
  dat <- if (!is.null(df_gd_agg) && y_col %in% names(df_gd_agg)) df_gd_agg else NULL
  if (is.null(dat)) return(NULL)
  p <- tryCatch(
    .no_title_frame(.make_std_plot(
      df          = dat,
      y_col       = y_col,
      y_label     = y_label,
      title       = "",
      caption_txt = .cap_all_sig(res),
      cld_df      = .trunc2(.cld(res, "^treatment$")),
      is_sig      = .is_sig(res, "Treatment"),
      facet_tp    = FALSE,
      n_txt       = .n_cap(dat, "agg")
    )),
    error = function(e) { warning("combo3 gd '", y_col, "' failed: ", e$message); NULL }
  )
  if (is.null(p)) return(NULL)
  if (isTRUE(keep_legend)) {
    p + ggplot2::theme(
      legend.position = "bottom",
      legend.text     = ggplot2::element_text(size = 16),
      legend.title    = ggplot2::element_text(size = 16, face = "bold"),
      legend.key.size = ggplot2::unit(1.2, "cm")
    )
  } else {
    p + ggplot2::theme(legend.position = "none")
  }
}

.p3_iid  <- .trunc_cld_in_plot(.mk_gd_big3(res_iid_agg,  "mean_iid_cm",        "Mean IID (cm)"))
.p3_nnd  <- .trunc_cld_in_plot(.mk_gd_big3(res_nnd_agg,  "mean_nnd_cm",        "Mean NND (cm)"))
.p3_hull <- .trunc_cld_in_plot(.mk_gd_big3(res_hull_agg, "mean_hull_area_cm2", "School area (cm²)",
                                             keep_legend = TRUE))

CLD_SIZE <- .OLD_CLD_SIZE3  # restore globals
PT_SIZE  <- .OLD_PT_SIZE3

# Was a 2x2 (iid|nnd)/(hull|pol); now a single row after the alignment panel was
# removed: A=IID, B=NND, C=Hull (hull carries the legend).
if (!any(sapply(list(.p3_iid, .p3_nnd, .p3_hull), is.null))) {
  .combo3 <- (.p3_iid | .p3_nnd | .p3_hull) +
    patchwork::plot_annotation(tag_levels = list(c("A", "B", "C"))) &
    ggplot2::theme(
      plot.tag          = ggplot2::element_text(face = "bold", size = 18),
      plot.tag.position = c(0.02, 0.98)
    )
  tryCatch(
    ggplot2::ggsave(file.path(.sg_stacked_dir, "combo_IID_NND_Hull.png"),
                    .combo3, width = 14, height = 10, dpi = 150),
    error = function(e) warning("combo3 ggsave failed: ", e$message)
  )
  ts_msg("stacked combo3 (IID/NND/Hull row, hull legend) saved")
} else {
  ts_msg("stacked combo3 skipped — one or more GD plots NULL")
}

ts_msg("stacked folder written: ", .sg_stacked_dir)


# =============================================================================
# ==== 6b) MANUSCRIPT FIGURES (Figures 7-10 + S1) ============================
# Composite manuscript figures written to STEP5_OUT/manu_graphs/.
# All inherit the existing pipeline aesthetics: theme_minimal(base_size=13),
# treatment colours (TREATMENT_COLORS), zone palette (color_map / color_map_broad),
# jittered raw points (PT_SIZE / PT_ALPHA), mean crossbar + SE errorbar.
# Panel tags A/B/C/... rendered in bold 18 pt at top-left outside the data
# region (plot.tag.position = c(0.02, 0.98)).
# =============================================================================
ts_msg("=== MANU GRAPHS ===")
.mg_dir <- file.path(STEP5_OUT, "manu_graphs")
dir.create(.mg_dir, recursive = TRUE, showWarnings = FALSE)
.mg_dir_new     <- file.path(PROJECT_ROOT, "manu_graphs_26.05.2026")
.mg_dir_new_pdf <- file.path(.mg_dir_new, "PDF")
dir.create(.mg_dir_new,     recursive = TRUE, showWarnings = FALSE)
dir.create(.mg_dir_new_pdf, recursive = TRUE, showWarnings = FALSE)
ts_msg("  Dated output folder (PNG): ", .mg_dir_new)
ts_msg("  Dated output folder (PDF): ", .mg_dir_new_pdf)

# Common patchwork tag theme. Tag positioned in the top-right margin of each
# panel so that it sits inside the 10 mm margin without colliding with data.
.MG_TAG_THEME <- ggplot2::theme(
  plot.tag          = ggplot2::element_text(face = "bold", size = 18),
  plot.tag.position = "topleft"
)

# Helper: drop title, white background, generous top margin so the bold
# panel tag has room without colliding with axis labels.
.mg_frame <- function(p) {
  if (is.null(p)) return(NULL)
  p + ggplot2::theme(
    plot.title       = ggplot2::element_blank(),
    plot.margin      = ggplot2::margin(10, 10, 10, 10, unit = "mm"),
    plot.background  = ggplot2::element_rect(fill = "white", colour = NA)
  )
}

# Helper: tighten x-axis labels for treatment scatter panels. The underlying
# .make_std_plot uses scale_x_continuous with str_to_title labels ("Control",
# "Exercise Choice"); when panels are narrow the two-word "Exercise Choice"
# label crowds against "Control". This helper replaces the scale labels with
# a stacked two-line version so they fit at any sensible panel width.
.mg_treatment_xaxis <- function(p) {
  if (is.null(p)) return(NULL)
  p + ggplot2::scale_x_continuous(
        breaks = c(1, 2),
        labels = c("Control", "Exercise\nchoice"),
        limits = c(0.5, 2.5), expand = c(0, 0)
      ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 12, face = "bold",
                                          lineheight = 0.9, hjust = 0.5)
    )
}

# Save helper:
#   (1) PNG + PDF  → STEP5_OUT/manu_graphs/   (standard pipeline location)
#   (2) PNG        → D:/CHOICE R SCRIPTS/manu_graphs_26.05.2026/
#   (3) PDF (font-embedded, cairo_pdf device)
#              → D:/CHOICE R SCRIPTS/manu_graphs_26.05.2026/PDF/
.mg_save <- function(p, filename, width_mm, height_mm, dpi = 300) {
  if (is.null(p)) { ts_msg("  manu_graphs: skipped ", filename, " (NULL)"); return(invisible()) }

  # --- (1) Standard pipeline location: PNG + PDF ---
  for (.ext in c("png", "pdf")) {
    .out <- file.path(.mg_dir, paste0(filename, ".", .ext))
    tryCatch(
      ggplot2::ggsave(.out, p, width = width_mm, height = height_mm,
                      units = "mm", dpi = dpi, bg = "white"),
      error = function(e) warning(filename, " ", .ext, " save failed: ", e$message)
    )
  }

  # --- (2) Dated folder: PNG ---
  if (exists(".mg_dir_new") && dir.exists(.mg_dir_new)) {
    .out_png <- file.path(.mg_dir_new, paste0(filename, ".png"))
    tryCatch(
      ggplot2::ggsave(.out_png, p, width = width_mm, height = height_mm,
                      units = "mm", dpi = dpi, bg = "white"),
      error = function(e) warning(filename, " dated PNG save failed: ", e$message)
    )
  }

  # --- (3) Dated PDF subfolder: font-embedded PDF via cairo_pdf device ---
  if (exists(".mg_dir_new_pdf") && dir.exists(.mg_dir_new_pdf)) {
    .out_pdf <- file.path(.mg_dir_new_pdf, paste0(filename, ".pdf"))
    tryCatch(
      ggplot2::ggsave(.out_pdf, p,
                      width  = width_mm / 25.4,   # ggsave with cairo_pdf expects inches
                      height = height_mm / 25.4,
                      units  = "in",
                      device = cairo_pdf,
                      bg     = "white"),
      error = function(e) warning(filename, " dated PDF save failed: ", e$message)
    )
  }

  ts_msg("  manu_graphs: saved ", filename, " (.png + .pdf | dated .png + cairo .pdf)")
}

# Helper: add per-group (N = X) labels just above each x-axis tick.
# Works for 2-group scatter panels that use numeric x positions 1 and 2.
# Expands the lower y-limit to create visible padding below the data region.
.mg_add_n_labels <- function(p, dat, y_col, tmap) {
  if (is.null(p) || is.null(dat) || !y_col %in% names(dat)) return(p)
  tryCatch({
    .n_df <- dat %>%
      dplyr::group_by(treatment) %>%
      dplyr::summarise(
        n = if ("phys_trial_id" %in% names(dat))
              dplyr::n_distinct(phys_trial_id) else dplyr::n(),
        .groups = "drop") %>%
      dplyr::mutate(
        .x    = unname(tmap[as.character(treatment)]),
        label = paste0("(N = ", n, ")"))
    .raw   <- dat[[y_col]][is.finite(dat[[y_col]])]
    .y_min <- min(.raw, na.rm = TRUE)
    .y_rng <- max(diff(range(.raw, na.rm = TRUE)), 1e-10)
    # Place label between x-axis line and y = 0 (or below y_min for panels
    # that already span negative values).
    if (.y_min >= 0) {
      .y_n      <- -0.08 * .y_rng
      .y_expand <- -0.18 * .y_rng
    } else {
      .y_n      <- .y_min - 0.10 * .y_rng
      .y_expand <- .y_min - 0.20 * .y_rng
    }
    .n_df$lbl_y <- .y_n
    p +
      ggplot2::expand_limits(y = .y_expand) +
      ggplot2::geom_text(data = .n_df,
        ggplot2::aes(x = .x, y = lbl_y, label = label),
        inherit.aes = FALSE, size = 3.0, colour = "grey35", vjust = 0.5,
        show.legend = FALSE)
  }, error = function(e) p)
}

# Helper: add per-zone-treatment (N = X) labels for grouped zone scatter panels
# (Jacobs D and alr-sec), where zones sit on the x-axis.
.mg_add_n_labels_zone <- function(p, dat, y_col, y_fixed_min = NULL) {
  if (is.null(p) || is.null(dat) || !y_col %in% names(dat)) return(p)
  tryCatch({
    .raw   <- dat[[y_col]][is.finite(dat[[y_col]])]
    .y_rng <- max(diff(range(.raw, na.rm = TRUE)), 1e-10)
    .y_data_min <- min(.raw, na.rm = TRUE)
    # For panels with a fixed y scale (e.g. Jacobs D [-1,1]), use fixed_min.
    .y_base <- if (!is.null(y_fixed_min)) y_fixed_min else .y_data_min
    .y_n      <- .y_base - 0.10 * .y_rng
    .y_expand <- .y_base - 0.20 * .y_rng
    .n_df <- dat %>%
      dplyr::group_by(treatment) %>%
      dplyr::summarise(
        n = if ("phys_trial_id" %in% names(dat))
              dplyr::n_distinct(phys_trial_id) else dplyr::n(),
        .groups = "drop") %>%
      dplyr::mutate(label = paste0("N = ", n))
    .n_unique <- dplyr::n_distinct(.n_df$n)
    # For a single-N-per-group display: show one centred label
    .lbl <- if (.n_unique == 1L) paste0("(N = ", .n_df$n[1], " per group)")
            else paste(.n_df$label, collapse = " / ")
    p +
      ggplot2::expand_limits(y = .y_expand) +
      ggplot2::annotate("text", x = Inf, y = .y_n, label = .lbl,
        hjust = 1.05, vjust = 0.5, size = 2.8, colour = "grey35")
  }, error = function(e) p)
}

# ---- Figure 7 — Trial-level zone preference (2x2: A/B=log-ratio, C/D=Jacobs D) ----
# Panel A: logit_flow = log(r_flow/(1-r_flow)) ~ treatment (main-zone alr)
# Panel B: lr_high/medium/low ~ zone x treatment (sub-zone alr; calm excluded as alr ref)
# Panel C: Jacobs' D for main zones (flow, calm)
# Panel D: Jacobs' D for sub-zones (high, medium, low, calm)

# Helper: scatter for a single log-ratio variable (allows negative y; dashed ref at 0).
.mg_make_alr_scatter <- function(df_in, y_col, y_label, res_obj, n_cap_df = NULL,
                                  show_n = TRUE, cap_csv = NULL, sig_from_csv = NULL) {
  if (is.null(df_in) || !y_col %in% names(df_in)) return(NULL)
  tryCatch({
    .d    <- df_in[is.finite(df_in[[y_col]]), ]
    if (nrow(.d) == 0) return(NULL)
    .lvls <- levels(factor(.d$treatment))
    .tmap <- stats::setNames(seq_along(.lvls), .lvls)
    .d$.x <- unname(.tmap[as.character(.d$treatment)])
    .smry <- .d %>%
      dplyr::group_by(treatment) %>%
      dplyr::summarise(mean_y = mean(.data[[y_col]], na.rm = TRUE),
                       se_y   = .sem(.data[[y_col]]),
                       .groups = "drop") %>%
      dplyr::mutate(.x = unname(.tmap[as.character(treatment)]))
    .n_df <- if (!show_n) NULL else if (is.null(n_cap_df)) df_in else n_cap_df
    .cap  <- if (!is.null(cap_csv)) .compose_caption(.n_cap(.n_df, "agg"), cap_csv) else
             .compose_caption(.n_cap(.n_df, "agg"), .cap_all_sig(res_obj))

    p <- ggplot2::ggplot(.d, ggplot2::aes(x = .x, y = .data[[y_col]],
                                           colour = treatment, shape = treatment)) +
      ggplot2::geom_point(
        position = ggplot2::position_jitter(width = JITTER_W, height = 0),
        alpha = PT_ALPHA, size = PT_SIZE) +
      ggplot2::geom_crossbar(data = .smry,
        ggplot2::aes(x = .x, y = mean_y, ymin = mean_y, ymax = mean_y),
        colour = "black", width = MEAN_W, linewidth = LW_MEAN, show.legend = FALSE) +
      ggplot2::geom_errorbar(data = .smry,
        ggplot2::aes(x = .x, y = mean_y, ymin = mean_y - se_y, ymax = mean_y + se_y),
        colour = "black", width = ERR_W, linewidth = LW_ERR, show.legend = FALSE) +
      ggplot2::scale_x_continuous(
        breaks = seq_along(.lvls), labels = stringr::str_to_title(.lvls),
        limits = c(0.5, length(.lvls) + 0.5), expand = c(0, 0)) +
      ggplot2::scale_colour_manual(values = TREATMENT_COLORS_PAL, name = "Treatment") +
      ggplot2::scale_shape_manual(values  = TREATMENT_SHAPES,     name = "Treatment") +
      BASE_THEME +
      ggplot2::labs(x = "Treatment", y = y_label, caption = .cap) +
      ggplot2::theme(legend.position = "bottom")

    # Asterisk significance bracket for 2-group comparisons via ggsignif.
    # Use CSV-parsed significance when provided, otherwise fall back to LMM result.
    .is_sig_val <- if (!is.null(sig_from_csv)) sig_from_csv else .is_sig(res_obj, "Treatment")
    if (.is_sig_val && length(.lvls) == 2L) {
      .raw  <- .d[[y_col]]
      .ym   <- max(.raw, na.rm = TRUE)
      .rng  <- max(diff(range(.raw, na.rm = TRUE)), 1e-10)
      .bary <- .ym + 0.10 * .rng
      .sig_df <- data.frame(xmin = 1, xmax = 2,
                            y_position = .bary, annotations = "*")
      p <- p +
        ggplot2::expand_limits(y = .ym + 0.22 * .rng) +
        ggsignif::geom_signif(
          data        = .sig_df,
          mapping     = ggplot2::aes(xmin = xmin, xmax = xmax,
                                     y_position = y_position,
                                     annotations = annotations),
          manual      = TRUE,
          tip_length  = 0.02,
          size        = 0.5,
          textsize    = 8,
          vjust       = 0.5,
          color       = "black",
          inherit.aes = FALSE
        )
    }
    # Per-group (N = X) labels above x-axis ticks — only when show_n = TRUE.
    if (show_n) p <- .mg_add_n_labels(p, .d, y_col, .tmap)
    p
  }, error = function(e) { warning("alr scatter '", y_col, "' failed: ", e$message); NULL })
}

# Helper: grouped scatter for secondary-zone log-ratios (lr_high/medium/low; calm excluded).
.mg_make_alr_sec_panel <- function() {
  if (is.null(df_sec_wide_agg)) return(NULL)
  .need <- intersect(c("lr_high", "lr_medium", "lr_low"), names(df_sec_wide_agg))
  if (length(.need) == 0) return(NULL)
  tryCatch({
    .lr_long <- df_sec_wide_agg %>%
      dplyr::select(dplyr::any_of(
        c("treatment", "phys_trial_id", "fish_density", .need))) %>%
      tidyr::pivot_longer(cols = dplyr::all_of(.need),
                          names_to = "zone_raw", values_to = "lr") %>%
      dplyr::mutate(zone = factor(sub("^lr_", "", zone_raw),
                                  levels = c("high", "medium", "low"))) %>%
      dplyr::filter(is.finite(lr))
    if (nrow(.lr_long) == 0) return(NULL)
    .smry <- .lr_long %>%
      dplyr::group_by(zone, treatment) %>%
      dplyr::summarise(mean_y = mean(lr, na.rm = TRUE),
                       se_y   = .sem(lr), .groups = "drop")
    .zone_pal <- color_map[c("high", "medium", "low")]
    .n_str <- .n_cap(df_sec_wide_agg, "agg")
    .p_str <- paste(Filter(nzchar, Filter(Negate(is.null),
                   lapply(list(res_zone_sec_high_agg,
                               res_zone_sec_med_agg,
                               res_zone_sec_low_agg),
                          .cap_all_sig))),
                   collapse = "\n")
    if (!nzchar(.p_str)) .p_str <- NULL
    .cap <- .compose_caption(.n_str, .p_str)

    ggplot2::ggplot(.lr_long,
                    ggplot2::aes(x = zone, y = lr, colour = zone, shape = treatment)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                         colour = "grey50", linewidth = 0.4) +
      ggplot2::geom_point(
        position = ggplot2::position_jitterdodge(jitter.width = JITTER_W,
                                                 dodge.width  = DODGE_W),
        alpha = PT_ALPHA, size = PT_SIZE) +
      ggplot2::geom_crossbar(data = .smry,
        ggplot2::aes(x = zone, y = mean_y, ymin = mean_y, ymax = mean_y,
                     group = interaction(zone, treatment)),
        position  = ggplot2::position_dodge(width = DODGE_W),
        width = MEAN_W, fatten = 0, linewidth = LW_MEAN,
        colour = "black", inherit.aes = FALSE) +
      ggplot2::geom_errorbar(data = .smry,
        ggplot2::aes(x = zone, ymin = mean_y - se_y, ymax = mean_y + se_y,
                     group = interaction(zone, treatment)),
        position  = ggplot2::position_dodge(width = DODGE_W),
        width = ERR_W, linewidth = LW_ERR, colour = "black", inherit.aes = FALSE) +
      ggplot2::scale_colour_manual(values = .zone_pal, name = "Zone", guide = "none") +
      ggplot2::scale_shape_manual(values  = TREATMENT_SHAPES, name = "Treatment") +
      BASE_THEME +
      ggplot2::labs(x = NULL, y = "log-ratio vs calm", caption = .cap) +
      ggplot2::theme(legend.position = "bottom") ->
      .p_sec
    .mg_add_n_labels_zone(.p_sec, .lr_long, "lr")
  }, error = function(e) { warning("alr sec panel failed: ", e$message); NULL })
}

# Helper: single-panel Jacobs D scatter for a subset of zones.
# jac_cld_suffix: reads per-zone Jacobs D CLDs from
#   jacobs_beta_{zone}_{jac_cld_suffix}/cld_treatment.csv for each zone in
#   zones_keep. Suppresses labels within a zone if both treatments share the
#   same CLD letter. Calm is silently skipped if no CLD file exists.
.mg_make_jacobs_panel <- function(jac_df, zones_keep, res_list = NULL,
                                   n_cap_df = NULL, label_dir = NULL) {
  if (is.null(jac_df)) return(NULL)
  tryCatch({
    .dat <- jac_df %>%
      dplyr::filter(zone %in% zones_keep, is.finite(D)) %>%
      dplyr::mutate(zone = factor(zone, levels = zones_keep))
    if (nrow(.dat) == 0) return(NULL)
    .smry <- .dat %>%
      dplyr::group_by(zone, treatment) %>%
      dplyr::summarise(D_mean = mean(D, na.rm = TRUE),
                       D_se   = .sem(D), .groups = "drop")
    .pal <- color_map[zones_keep]; .pal <- .pal[!is.na(.pal)]
    .n_df  <- if (is.null(n_cap_df)) jac_df else n_cap_df
    .n_str <- .n_cap(.n_df, "agg")
    .p_str <- if (!is.null(res_list)) {
      .caps <- Filter(nzchar, Filter(Negate(is.null),
                                     lapply(res_list, .cap_all_sig)))
      .caps <- lapply(.caps, function(x) gsub("\n", "<br>", x, fixed = TRUE))
      if (length(.caps)) paste(unlist(.caps), collapse = "<br>") else NULL
    } else NULL
    # Each ANOVA term on its own line; blank line between n_str and ANOVA block.
    .cap_parts <- c(.n_str, .p_str)
    .cap_parts <- .cap_parts[!is.null(.cap_parts) & !is.na(.cap_parts) &
                              nzchar(.cap_parts)]
    .cap <- if (length(.cap_parts) == 0) NULL else paste(.cap_parts, collapse = "<br><br>")

    ggplot2::ggplot(.dat, ggplot2::aes(x = zone, y = D,
                                       colour = zone, shape = treatment)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                         colour = "grey50", linewidth = 0.6) +
      ggplot2::geom_point(
        position = ggplot2::position_jitterdodge(jitter.width = JITTER_W,
                                                 dodge.width  = DODGE_W),
        alpha = PT_ALPHA, size = PT_SIZE) +
      ggplot2::geom_crossbar(data = .smry,
        ggplot2::aes(x = zone, y = D_mean, ymin = D_mean, ymax = D_mean,
                     group = interaction(zone, treatment)),
        position  = ggplot2::position_dodge(width = DODGE_W),
        width = MEAN_W, fatten = 0, linewidth = LW_MEAN,
        colour = "black", inherit.aes = FALSE) +
      ggplot2::geom_errorbar(data = .smry,
        ggplot2::aes(x = zone, ymin = D_mean - D_se, ymax = D_mean + D_se,
                     group = interaction(zone, treatment)),
        position  = ggplot2::position_dodge(width = DODGE_W),
        width = ERR_W, linewidth = LW_ERR, colour = "black", inherit.aes = FALSE) +
      ggplot2::scale_colour_manual(values = .pal, name = "Zone", guide = "none") +
      ggplot2::scale_shape_manual(values  = TREATMENT_SHAPES, name = "Treatment") +
      ggplot2::scale_y_continuous(limits = c(-1.3, 1),
                                  breaks = c(-1, -0.5, 0, 0.5, 1)) +
      BASE_THEME +
      ggplot2::labs(x = NULL, y = "Jacob's D index", caption = .cap) +
      ggplot2::theme(legend.position = "bottom") ->
      .p_jac
    .p_jac <- .mg_add_n_labels_zone(.p_jac, .dat, "D", y_fixed_min = -1)

    # CLD letters per (zone x treatment) cell from joint Jacobs D beta-GLMM
    # (treatment * zone fixed effects). Suppress within a zone if both
    # treatments share the same CLD letter.
    .cld_jac <- if (!is.null(label_dir)) .mg_read_cld_treatxzone(label_dir) else NULL

    if (!is.null(.cld_jac) && nrow(.cld_jac) > 0) {
      # Keep only zones where the two treatments have different CLD letters
      .cld_show <- do.call(rbind, lapply(zones_keep, function(z) {
        rows <- .cld_jac[.cld_jac$zone == z, , drop = FALSE]
        if (nrow(rows) < 2 || length(unique(rows$.group)) == 1L) return(NULL)
        rows
      }))
      if (!is.null(.cld_show) && nrow(.cld_show) > 0) {
        .smry2 <- merge(.smry, .cld_show, by = c("treatment", "zone"), all.x = FALSE)
        .smry2 <- .smry2[!is.na(.smry2$.group), , drop = FALSE]
        if (nrow(.smry2) > 0) {
          .trt_lvls <- c("control", "exercise choice")
          .tx_off   <- stats::setNames(
            seq(-DODGE_W / 2, DODGE_W / 2, length.out = length(.trt_lvls)),
            .trt_lvls)
          .smry2$cell_x <- as.integer(factor(.smry2$zone, levels = zones_keep)) +
                            unname(.tx_off[.smry2$treatment])
          # Uniform CLD height: max raw datapoint across whole panel + 10% of
          # the full data range. All labels sit at the same y level.
          .D_all_max   <- max(.dat$D, na.rm = TRUE)
          .D_all_min   <- min(.dat$D, na.rm = TRUE)
          .D_full_range <- if (is.finite(.D_all_max - .D_all_min) && .D_all_max > .D_all_min)
                             .D_all_max - .D_all_min else 1
          .cld_y_unif  <- .D_all_max + 0.10 * .D_full_range
          .smry2$cld_y <- .cld_y_unif
          .p_jac <- .p_jac +
            ggplot2::geom_text(
              data = .smry2,
              ggplot2::aes(x = cell_x, y = cld_y, label = .group),
              inherit.aes = FALSE, fontface = "bold", size = 6,
              colour = "black", show.legend = FALSE) +
            ggplot2::coord_cartesian(ylim = c(-1.3, max(1, .cld_y_unif + 0.05 * .D_full_range)))
        }
      }
    }
    .p_jac
  }, error = function(e) { warning("Jacobs panel failed: ", e$message); NULL })
}

# ---- Cell-mean panel helpers (used by new Figure 7) -------------------------
# Helper: wrap caption text to a per-panel character width (forward-declared
# here because .mg_make_cells_panel below uses it; the canonical definition
# at line ~5234 is identical).
if (!exists(".mg_wrap_caption", inherits = FALSE)) {
  .mg_wrap_caption <- function(text, panel_width_chars = 26) {
    # No longer wraps (2026-08-09): the ANOVA caption must render on one
    # line. panel_width_chars is kept as a no-op parameter so call sites
    # don't need editing; text is returned unchanged.
    text
  }
}

# Reads the (treatment x zone) CLD table written by run_lmm_analysis (or the
# beta-GLMM analogue). Trims whitespace, truncates letters to max 2 chars.
.mg_read_cld_treatxzone <- function(label_dir) {
  fp <- file.path(STEP5_OUT, label_dir, "cld_treatmentxzone.csv")
  if (!file.exists(fp)) return(NULL)
  d <- tryCatch(read.csv(fp, stringsAsFactors = FALSE),
                error = function(e) NULL)
  if (is.null(d) || !all(c("treatment","zone",".group") %in% names(d)))
    return(NULL)
  d$.group    <- substr(trimws(as.character(d$.group)), 1, 2)
  d$treatment <- trimws(as.character(d$treatment))
  d$zone      <- trimws(as.character(d$zone))
  d <- d[nzchar(d$.group), , drop = FALSE]
  if (nrow(d) == 0) return(NULL)
  d
}

# ANOVA caption for cell-mean panels: emits the treatment x zone interaction
# F/chisq stat only when p < 0.05 (matches .fig8_anova_cap formatting).
.fig7_anova_cap <- function(res) {
  if (is.null(res) || is.null(res$anova)) return(NULL)
  av  <- res$anova
  row <- av[grepl("treatment.*zone$|^zone.*treatment$|treatment:zone|zone:treatment",
                   av$term, ignore.case = TRUE, perl = TRUE), ]
  if (nrow(row) == 0L) return(NULL)
  row <- row[1L, ]
  p   <- suppressWarnings(as.numeric(row$p_value))
  if (!is.finite(p)) return(NULL)
  p_str <- fmt_p(p)
  chi   <- fmt_Fstat(as.numeric(row$chisq))
  st    <- attr(av, "stat_type") %||% "Wald-chisq"
  is_f  <- isTRUE(grepl("^F", st))
  if (is_f) {
    df1 <- as.integer(row$df)
    df2 <- suppressWarnings(as.numeric(
             if (!is.null(row$df_denom)) row$df_denom else NA_real_))
    if (is.finite(df2)) {
      eta <- fmt_es3((as.numeric(row$chisq) * df1) / (as.numeric(row$chisq) * df1 + df2))
      sprintf("Treatment × Zone: F<sub>%d, %s</sub> = %s, %s, &eta;<sup>2</sup><sub>p</sub> = %s",
              df1, fmt_F(df2), chi, p_str, eta)
    } else
      sprintf("Treatment × Zone: F<sub>%d</sub> = %s, %s", df1, chi, p_str)
  } else {
    sprintf("Treatment × Zone: χ²<sub>%d</sub> = %s, %s",
            as.integer(row$df), chi, p_str)
  }
}

# Remap CLD letters for zone × treatment cell panels. Precedence: zones_order
# first (e.g. flow > calm or high > medium > low > calm), then treatments_order
# within each zone (control > EC). First new letter seen → "a", etc.
.mg_remap_cld_zone <- function(cld_df, zones_order,
                                treatments_order = c("control","exercise choice")) {
  if (is.null(cld_df) || !".group" %in% names(cld_df)) return(cld_df)
  prec_idx <- order(
    match(cld_df$zone,      zones_order),
    match(cld_df$treatment, treatments_order)
  )
  seen <- character(0)
  for (i in prec_idx) {
    lets <- strsplit(cld_df$.group[i], "")[[1]]
    for (l in lets) if (!l %in% seen) seen <- c(seen, l)
  }
  if (!length(seen)) return(cld_df)
  letter_map <- stats::setNames(letters[seq_along(seen)], seen)
  cld_df$.group <- vapply(cld_df$.group, function(g) {
    lets     <- strsplit(g, "")[[1]]
    remapped <- letter_map[lets]
    remapped <- remapped[!is.na(remapped)]
    if (!length(remapped)) return(g)
    paste(sort(remapped), collapse = "")
  }, character(1L), USE.NAMES = FALSE)
  cld_df
}

# Cell-mean scatter panel: x = zone (categorical), treatment dodged within
# each zone group. Cell mean ± SE crossbar/errorbar; trial-level points
# jittered around each cell. CLD letters per cell sit just above each
# cell's +SE tip (7% of yrange offset, black, size 6 — same as Figure 9).
.mg_make_cells_panel <- function(df_long, y_col, y_label, zones_order,
                                  label_dir, res_obj,
                                  hline_at = NULL,
                                  panel_width_chars = 30,
                                  treatments_order = c("control","exercise choice"),
                                  remap_cld = FALSE,
                                  always_show_cld = FALSE,
                                  show_zone_legend = TRUE,
                                  show_treatment_legend = TRUE,
                                  cap_csv = NULL) {
  if (is.null(df_long) || !y_col %in% names(df_long) ||
      !all(c("treatment","zone") %in% names(df_long))) return(NULL)
  tryCatch({
    .d <- df_long
    .d$treatment <- factor(as.character(.d$treatment), levels = treatments_order)
    .d$zone      <- factor(as.character(.d$zone),      levels = zones_order)
    .d <- .d[is.finite(.d[[y_col]]) & !is.na(.d$treatment) & !is.na(.d$zone), ]
    if (nrow(.d) == 0) return(NULL)

    # Build cell coordinates: zones on integer x; treatments dodged ±0.18.
    .n_trts <- length(treatments_order)
    .dodge  <- 0.36 / .n_trts
    .tx_off <- stats::setNames(
      ((seq_along(treatments_order) - mean(seq_along(treatments_order))) * 2) * .dodge,
      treatments_order)
    .d$cell_x <- as.integer(.d$zone) + unname(.tx_off[as.character(.d$treatment)])

    .smry <- .d %>%
      dplyr::group_by(treatment, zone) %>%
      dplyr::summarise(mean_y = mean(.data[[y_col]], na.rm = TRUE),
                       sem_y  = .sem(.data[[y_col]]),
                       max_y  = max(.data[[y_col]], na.rm = TRUE),
                       .groups = "drop") %>%
      dplyr::mutate(cell_x = as.integer(factor(zone, levels = zones_order)) +
                              unname(.tx_off[as.character(treatment)]))

    ymin <- min(c(.smry$mean_y - .smry$sem_y, .d[[y_col]]), na.rm = TRUE)
    ymax <- max(c(.smry$mean_y + .smry$sem_y, .d[[y_col]]), na.rm = TRUE)
    yrange <- ymax - ymin
    if (!is.finite(yrange) || yrange == 0) yrange <- 1

    p <- ggplot2::ggplot(.d, ggplot2::aes(x = cell_x, y = .data[[y_col]],
                                          colour = zone, shape = treatment))
    if (!is.null(hline_at))
      p <- p + ggplot2::geom_hline(yintercept = hline_at, linetype = "dashed",
                                   colour = "grey50", linewidth = 0.4)
    p <- p +
      ggplot2::geom_point(
        position = ggplot2::position_jitter(width = 0.05, height = 0),
        alpha = PT_ALPHA, size = PT_SIZE) +
      ggplot2::geom_crossbar(
        data = .smry,
        ggplot2::aes(x = cell_x, y = mean_y, ymin = mean_y, ymax = mean_y),
        colour = "black", width = 0.20, linewidth = LW_MEAN, show.legend = FALSE,
        inherit.aes = FALSE) +
      ggplot2::geom_errorbar(
        data = .smry,
        ggplot2::aes(x = cell_x, y = mean_y, ymin = mean_y - sem_y,
                     ymax = mean_y + sem_y),
        colour = "black", width = 0.10, linewidth = LW_ERR, show.legend = FALSE,
        inherit.aes = FALSE) +
      ggplot2::scale_x_continuous(
        breaks = seq_along(zones_order),
        labels = stringr::str_to_title(zones_order),
        limits = c(0.5, length(zones_order) + 0.5),
        expand = c(0, 0)) +
      ggplot2::scale_colour_manual(values = color_map[zones_order], name = "Zone",
        guide = if (show_zone_legend) "legend" else "none") +
      ggplot2::scale_shape_manual(values  = TREATMENT_SHAPES,    name = "Treatment",
        guide = if (show_treatment_legend) "legend" else "none") +
      BASE_THEME +
      ggplot2::labs(x = "Zone", y = y_label,
                    caption = .mg_wrap_caption(
                      if (!is.null(cap_csv)) cap_csv else .fig7_anova_cap(res_obj),
                      panel_width_chars)) +
      ggplot2::theme(
        # 14 -> 11 on 2026-08-09: these cell-mean panels (Figure 8C/8D) carry
        # the treatment x zone caption, and once eta^2_p joined it on a single
        # line the sub-zone version ("Treatment x Zone: F3, 42 = 8.131,
        # p < 0.001, eta^2_p = 0.367", ~55 chars) ran past the panel edge and
        # lost the effect size. Confirmed by inspecting the rendered PNG.
        plot.caption = ggtext::element_markdown(size = 11, hjust = 0,
                                                face = "italic", lineheight = 1.3,
                                                margin = ggplot2::margin(t = 10))
      )

    # CLD letters per cell (one per treatment x zone), black, bold, size 6,
    # centred on the cell and placed 15% of the panel's data range above that
    # cell's highest plotted point (author decision, 2026-08-09) — the previous
    # mean+SE anchor put letters inside the point cloud wherever a cell's
    # scatter extended well above its own error bar.
    cld_df <- .mg_read_cld_treatxzone(label_dir)
    if (remap_cld && !is.null(cld_df))
      cld_df <- .mg_remap_cld_zone(cld_df, zones_order, treatments_order)
    if (!is.null(cld_df)) {
      cld_join <- merge(cld_df,
                        .smry[, c("treatment","zone","mean_y","sem_y","max_y","cell_x")],
                        by = c("treatment","zone"), all.x = TRUE)
      cld_join <- cld_join[!is.na(cld_join$mean_y), , drop = FALSE]
      if (nrow(cld_join) > 0 && (always_show_cld || length(unique(cld_join$.group)) > 1L)) {
        cld_join$y <- cld_join$max_y + 0.15 * yrange
        p <- p + ggplot2::geom_text(
          data = cld_join,
          ggplot2::aes(x = cell_x, y = y, label = .group),
          inherit.aes = FALSE, fontface = "bold", size = 6,
          colour = "black", show.legend = FALSE)
        p <- p + ggplot2::coord_cartesian(
          ylim = c(ymin - 0.05 * yrange, ymax + 0.26 * yrange))
      }
    }
    p
  }, error = function(e) {
    warning("cells panel '", y_col, "' failed: ", e$message); NULL
  })
}

# ---- Layout helpers: panel-body / caption / legend-ghost stacking -----------
# Relative height constants used in plot_layout(heights = ...) for the
# body / legend-ghost and body / spacer patchwork stacks below.
# .MG_SPC_H = 0.5 × .MG_LEG_H gives 0.5× caption-to-bottom padding for
# panels that carry no legend.
.MG_BODY_H <- 10
.MG_LEG_H  <- 1.5
.MG_SPC_H  <- 0.75

# Ghost theme: applied to a copy of a legend-carrying panel to produce a
# legend-only strip (all plot elements blanked, caption suppressed).
# Stacking body / ghost gives the order: plot → ANOVA caption → legend.
.MG_GHOST_THEME <- ggplot2::theme(
  panel.background = ggplot2::element_blank(),
  panel.border     = ggplot2::element_blank(),
  panel.grid.major = ggplot2::element_blank(),
  panel.grid.minor = ggplot2::element_blank(),
  axis.line        = ggplot2::element_blank(),
  axis.title       = ggplot2::element_blank(),
  axis.text        = ggplot2::element_blank(),
  axis.ticks       = ggplot2::element_blank(),
  plot.background  = ggplot2::element_blank(),
  plot.caption     = ggplot2::element_blank(),
  plot.margin      = ggplot2::margin(0, 0, 0, 0),
  legend.margin    = ggplot2::margin(2, 0, 0, 0)
)

# Caption theme override for Figure_10 (interval line plots): black text.
# Was 13 pt to match legend text; dropped to 11 on 2026-08-09 when partial
# eta-squared was appended to these captions (longest statement grew 43 -> 56
# characters, which at 13 pt would have run to the panel edge in this 2-column
# layout). Measured against the regenerated caption table, not guessed.
.FIG10_CAP_OVR <- ggplot2::theme(
  plot.caption = ggtext::element_markdown(
    size   = 11, hjust  = 0, face   = "italic",
    colour = "black", margin = ggplot2::margin(t = 10)
  )
)

# ---------------------------------------------------------------------------
# Parse p-value from an ANOVA statement string (e.g. "Treatment F₁,₁₁ = 10.348, p = 0.008"
# or "Treatment χ²₁ = 19.004, p < 0.001") and return TRUE if p < 0.05.
# Returns FALSE if the string is NULL, NA, or the p-value cannot be parsed.
# Convert Unicode subscript digits (U+2080–U+2089) to <sub>...</sub> HTML so
# that ggtext::element_markdown renders them correctly in both PNG and PDF.
# Standard PDF fonts (Helvetica, Arial) lack these codepoints, causing blank
# squares; <sub> tags are handled by ggtext's own renderer and are device-agnostic.
.sub_to_html <- function(s) {
  if (is.null(s) || is.na(s) || !nzchar(s)) return(s)
  sub_chars <- "₀₁₂₃₄₅₆₇₈₉"
  # Match runs of subscript digits (with optional embedded commas for df1,df2).
  pat <- paste0("[", sub_chars, "][", sub_chars, ",]*[", sub_chars, "]|[", sub_chars, "]")
  stringr::str_replace_all(s, pat, function(m) {
    paste0("<sub>", chartr(sub_chars, "0123456789", m), "</sub>")
  })
}
# ---------------------------------------------------------------------------

.csv_is_sig <- function(statement) {
  if (is.null(statement) || is.na(statement) || !nzchar(statement)) return(FALSE)
  # "p < 0.001" → always significant
  if (grepl("p\\s*<\\s*0\\.0", statement)) return(TRUE)
  # "p = 0.NNN" → extract numeric value
  m <- regmatches(statement, regexpr("p\\s*=\\s*([0-9]+\\.[0-9]+)", statement))
  if (length(m) == 0 || !nzchar(m)) return(FALSE)
  pval <- suppressWarnings(as.numeric(sub("p\\s*=\\s*", "", m)))
  is.finite(pval) && pval < 0.05
}
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Trial-level ANOVA captions — sourced from anova_statement_table_trial.csv.
# Lookup by variable name (column "variable"); returns the "statement" string,
# or NULL if the file or row is missing. Used for Figure_8 (zone/ALR) and
# Figure_9 (collective trial-level) panels.
.fig9_cap_csv <- local({
  .csv_path <- file.path(PROJECT_ROOT, "choice R pipeline/easy_scripts/standalone/anova_statement_table_trial.csv")
  .tbl <- tryCatch(read.csv(.csv_path, stringsAsFactors = FALSE, encoding = "UTF-8"),
                   error = function(e) NULL)
  function(variable) {
    if (is.null(.tbl)) return(NULL)
    row <- .tbl[trimws(.tbl$variable) == variable, , drop = FALSE]
    if (nrow(row) == 0) return(NULL)
    .sub_to_html(trimws(row$statement[1]))
  }
})
# ---------------------------------------------------------------------------

# ---- Figure 7: 2 x 2 — alr panels (top) + cell-mean CLD panels (bottom) ----
.fig7_A <- .mg_frame(.mg_make_alr_scatter(
  df_in        = if (!is.null(df_main_wide_agg) && "logit_flow" %in% names(df_main_wide_agg))
                   df_main_wide_agg else NULL,
  y_col        = "logit_flow",
  y_label      = "ALR(flow vs calm)",
  res_obj      = if (exists("res_zone_flow_logit_agg")) res_zone_flow_logit_agg else NULL,
  n_cap_df     = df_main_wide_agg,
  show_n       = FALSE,
  cap_csv      = .fig9_cap_csv("logit_flow"),
  sig_from_csv = .csv_is_sig(.fig9_cap_csv("logit_flow"))))

.fig7_B <- .mg_frame(.mg_make_alr_scatter(
  df_in        = if (!is.null(df_sec_wide_agg) && "lr_high" %in% names(df_sec_wide_agg))
                   df_sec_wide_agg else NULL,
  y_col        = "lr_high",
  y_label      = "ALR(high vs calm)",
  res_obj      = if (exists("res_zone_sec_high_agg")) res_zone_sec_high_agg else NULL,
  n_cap_df     = df_sec_wide_agg,
  show_n       = FALSE,
  cap_csv      = .fig9_cap_csv("lr_high"),
  sig_from_csv = .csv_is_sig(.fig9_cap_csv("lr_high"))))

# Trial-level cell-mean ANOVA captions -- sourced from
# anova_statement_table_cells.csv (rebuild_figure_caption_table_cells.R).
# Added 2026-08-09 alongside figures_step5_module.R: this is the one panel
# family whose caption previously came ONLY from a live res_obj$anova table
# (.fig7_anova_cap()), with no CSV cache like every sibling panel -- the last
# obstacle to a model-fit-free figure-only regeneration path. res_obj is kept
# as a fallback so the live pipeline is unaffected if the cache is stale/missing.
.fig7_cells_cap_csv <- local({
  .csv_path <- file.path(PROJECT_ROOT, "choice R pipeline/easy_scripts/standalone/anova_statement_table_cells.csv")
  .tbl <- tryCatch(read.csv(.csv_path, stringsAsFactors = FALSE, encoding = "UTF-8"),
                   error = function(e) NULL)
  function(variable) {
    if (is.null(.tbl)) return(NULL)
    row <- .tbl[trimws(.tbl$variable) == variable, , drop = FALSE]
    if (nrow(row) == 0 || is.na(row$statement[1])) return(NULL)
    trimws(row$statement[1])
  }
})

.fig7_C <- .mg_frame(.mg_make_cells_panel(
  df_long               = df_main_cells_agg,
  y_col                 = "p",
  y_label               = "Zone occupancy",
  zones_order           = c("flow","calm"),
  label_dir             = "zone_main_cells_aggregated",
  res_obj               = res_zone_main_cells_agg,
  cap_csv               = .fig7_cells_cap_csv("zone_main"),
  hline_at              = NULL,
  panel_width_chars     = 32,
  remap_cld             = TRUE,
  show_zone_legend      = FALSE,
  show_treatment_legend = TRUE))

.fig7_D <- .mg_frame(.mg_make_cells_panel(
  df_long               = df_sec_cells_agg,
  y_col                 = "clr",
  y_label               = "area normalized CLR",
  zones_order           = c("high","medium","low","calm"),
  label_dir             = "zone_sec_cells_aggregated",
  res_obj               = res_zone_sec_cells_agg,
  cap_csv               = .fig7_cells_cap_csv("zone_sec"),
  hline_at              = NULL,
  panel_width_chars     = 32,
  remap_cld             = TRUE,
  show_zone_legend      = FALSE,
  show_treatment_legend = FALSE))

# ============ Line-with-Tukey panel helpers (moved up from the Figure 9
# section, 2026-08-09) ============
# .mg_make_line_with_tukey() and its dependents were originally defined only
# where Figure 9 (interval) first needed them, further down. Figure 7bis
# (below) now also needs it, moved earlier so it can double as Figure_8 panel
# F -- so the whole self-contained helper bundle (pure function definitions,
# no side effects) moved up here, ahead of every place that calls it. Note:
# .mg_wrap_caption is redefined identically here to the forward-declared
# no-op version above (both just return text unchanged) -- harmless.
# Helper: build raw cell means/SE table from a long df by treatment x timepoint.
.mg_cell_summary <- function(df_in, y_col) {
  if (is.null(df_in) || !y_col %in% names(df_in) ||
      !"timepoint_f" %in% names(df_in)) return(NULL)
  df_in %>%
    dplyr::filter(is.finite(.data[[y_col]]), !is.na(timepoint_f)) %>%
    dplyr::group_by(treatment, timepoint_f) %>%
    dplyr::summarise(mean_y = mean(.data[[y_col]], na.rm = TRUE),
                     sem_y  = .sem(.data[[y_col]]),
                     .groups = "drop") %>%
    dplyr::mutate(
      tp_num = as.integer(as.character(timepoint_f)),
      cell   = paste(as.character(treatment), as.character(timepoint_f))
    )
}

# Helper: parse Tukey contrast label "control timepoint1 - exercise choice timepoint3"
# into a row of (cell1_treatment, cell1_tp, cell2_treatment, cell2_tp, p).
.mg_parse_contrast <- function(contr_row) {
  txt <- as.character(contr_row$contrast)
  parts <- strsplit(txt, " - ", fixed = TRUE)[[1]]
  if (length(parts) != 2) return(NULL)
  pat_tp <- "^(.+?)\\s+timepoint_?f?(\\d)$"
  m1 <- regmatches(parts[1], regexec(pat_tp, parts[1]))[[1]]
  m2 <- regmatches(parts[2], regexec(pat_tp, parts[2]))[[1]]
  if (length(m1) != 3 || length(m2) != 3) return(NULL)
  data.frame(
    tx1 = m1[2], tp1 = as.integer(m1[3]),
    tx2 = m2[2], tp2 = as.integer(m2[3]),
    p   = contr_row$p.value,
    stringsAsFactors = FALSE
  )
}

.mg_sig_stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  ""
}

# Helper: wrap subtitle text to a per-panel character width.
.mg_wrap_caption <- function(text, panel_width_chars = 26) {
  # No longer wraps (2026-08-09): the ANOVA caption must render on one
  # line. panel_width_chars is kept as a no-op parameter so call sites
  # don't need editing; text is returned unchanged.
  text
}

# Helper: read the (treatment x timepoint) CLD table written by
# run_lmm_analysis. Trims whitespace, truncates letters to max 2 chars,
# coerces timepoint to integer.
.mg_read_cld_tp <- function(label_dir, cld_path = NULL) {
  # cld_path: explicit CSV path, for panels whose CLD is not produced by a
  # STEP5 model directory (commitment_index comes from the bout engine).
  if (!is.null(cld_path) && nzchar(cld_path) && file.exists(cld_path)) {
    use_fp <- cld_path
  } else {
  if (is.null(label_dir) || !nzchar(label_dir)) return(NULL)
  fp   <- file.path(STEP5_OUT, label_dir, "cld_treatmentxtimepoint.csv")
  fp_f <- file.path(STEP5_OUT, label_dir, "cld_treatmentxtimepoint_f.csv")
  # Prefer continuous-time CLD; fall back to factor-time CLD.
  # hull_area_timepoint (and any other indicator where factor timepoint wins AICc)
  # only writes the _f variant, so the fallback is necessary.
  use_fp <- if (file.exists(fp)) fp else if (file.exists(fp_f)) fp_f else return(NULL)
  }
  d <- tryCatch(read.csv(use_fp, stringsAsFactors = FALSE),
                error = function(e) NULL)
  if (is.null(d) || !all(c("treatment",".group") %in% names(d))) return(NULL)
  # Support "timepoint" (continuous model) or "timepoint_f" (factor model) as tp source
  tp_col <- if ("timepoint" %in% names(d)) "timepoint" else
            if ("timepoint_f" %in% names(d)) "timepoint_f" else return(NULL)
  d$.group    <- substr(trimws(as.character(d$.group)), 1, 2)
  d$treatment <- trimws(as.character(d$treatment))
  d$tp_num    <- suppressWarnings(as.integer(as.character(d[[tp_col]])))
  d <- d[!is.na(d$tp_num) & nzchar(d$.group), , drop = FALSE]
  if (nrow(d) == 0) return(NULL)
  d
}

# Build the per-cell label table for plotting. Places a label just above the
# +SE error bar tip (7% of the data range gap) per (tp, treatment); when both
# treatments at the same tp share the same letter string, emits a single
# merged label above the higher of the two SE tips. Each row carries x_adj = 0
# (callers may override specific cells via cld_nudges).
# collapse_tp: integer vector of tp values where a single label is shown above
# max(both tops)+off, using the first letter of the Control treatment's group.
.mg_cld_label_rows <- function(cld_df, smry, yrange, force_split_tp = integer(0),
                                collapse_tp = integer(0)) {
  off <- 0.11 * yrange
  out <- list()
  for (tp in sort(unique(cld_df$tp_num))) {
    rows <- cld_df[cld_df$tp_num == tp, , drop = FALSE]
    sm   <- smry[smry$tp_num == tp, , drop = FALSE]
    if (!nrow(rows) || !nrow(sm)) next
    sm$top <- sm$mean_y + sm$sem_y
    rows <- merge(rows, sm[, c("treatment","top")],
                  by = "treatment", all.x = TRUE)
    rows <- rows[!is.na(rows$top), , drop = FALSE]
    if (!nrow(rows)) next
    if (tp %in% collapse_tp) {
      prim <- rows[rows$treatment == "control", , drop = FALSE]
      if (!nrow(prim)) prim <- rows[1L, , drop = FALSE]
      out[[length(out) + 1L]] <- data.frame(
        tp_num    = tp,
        treatment = NA_character_,
        y         = max(rows$top, na.rm = TRUE) + off,
        label     = substr(prim$.group[1], 1L, 1L),
        merged    = TRUE,
        x_adj     = 0,
        stringsAsFactors = FALSE
      )
    } else if (length(unique(rows$.group)) == 1L && !tp %in% force_split_tp) {
      out[[length(out) + 1L]] <- data.frame(
        tp_num    = tp,
        treatment = NA_character_,
        y         = max(rows$top, na.rm = TRUE) + off,
        label     = rows$.group[1],
        merged    = TRUE,
        x_adj     = 0,
        stringsAsFactors = FALSE
      )
    } else {
      for (i in seq_len(nrow(rows))) {
        out[[length(out) + 1L]] <- data.frame(
          tp_num    = tp,
          treatment = rows$treatment[i],
          y         = rows$top[i] + off,
          label     = rows$.group[i],
          merged    = FALSE,
          x_adj     = 0,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (!length(out)) return(NULL)
  do.call(rbind, out)
}

# Remap CLD letter assignments so the first group encountered under
# precedence order (control → EC; interval 1→2→3) is labelled "a",
# the second new letter "b", etc.  Letters within each .group string
# are re-sorted alphabetically after remapping.
.mg_remap_cld <- function(cld_df) {
  if (is.null(cld_df) || !".group" %in% names(cld_df)) return(cld_df)
  prec_idx <- order(
    match(cld_df$treatment, c("control", "exercise choice")),
    cld_df$tp_num
  )
  seen <- character(0)
  for (i in prec_idx) {
    lets <- strsplit(cld_df$.group[i], "")[[1]]
    for (l in lets) if (!l %in% seen) seen <- c(seen, l)
  }
  if (!length(seen)) return(cld_df)
  letter_map <- stats::setNames(letters[seq_along(seen)], seen)
  cld_df$.group <- vapply(cld_df$.group, function(g) {
    lets     <- strsplit(g, "")[[1]]
    remapped <- letter_map[lets]
    remapped <- remapped[!is.na(remapped)]
    if (!length(remapped)) return(g)
    paste(sort(remapped), collapse = "")
  }, character(1L), USE.NAMES = FALSE)
  cld_df
}

.mg_make_line_with_tukey <- function(df_src, y_col, y_label, label_dir,
                                      anova_caption,
                                      panel_width_chars = 20,
                                      cld_nudges = NULL,
                                      remap_cld = FALSE,
                                      collapse_tp = integer(0),
                                      cld_path = NULL) {
  smry <- .mg_cell_summary(df_src, y_col)
  if (is.null(smry) || nrow(smry) == 0) return(NULL)

  cld_df <- .mg_read_cld_tp(label_dir, cld_path = cld_path)
  if (remap_cld && !is.null(cld_df)) cld_df <- .mg_remap_cld(cld_df)

  ymin <- min(smry$mean_y - smry$sem_y, na.rm = TRUE)
  ymax <- max(smry$mean_y + smry$sem_y, na.rm = TRUE)
  yrange <- ymax - ymin
  if (!is.finite(yrange) || yrange == 0) yrange <- 1

  p <- ggplot2::ggplot(
      smry,
      ggplot2::aes(x = tp_num, y = mean_y,
                   colour   = treatment,
                   shape    = treatment,
                   linetype = treatment,
                   group    = treatment)
    ) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = PT_SIZE * 1.6) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = mean_y - sem_y, ymax = mean_y + sem_y),
      colour = "black", width = 0.10, linewidth = LW_ERR,
      show.legend = FALSE
    ) +
    ggplot2::scale_colour_manual(values = TREATMENT_COLORS_PAL, name = "Treatment") +
    ggplot2::scale_shape_manual(values  = TREATMENT_SHAPES,    name = "Treatment") +
    ggplot2::scale_linetype_manual(
      values = c("exercise choice" = "solid", "control" = "dashed"),
      name   = "Treatment"
    ) +
    ggplot2::scale_x_continuous(
      breaks = 1:3,
      labels = c("5–25", "45–65", "85–105"),
      limits = c(0.6, 3.4)
    ) +
    BASE_THEME +
    ggplot2::labs(x = "Interval (min)", y = y_label,
                  caption = .mg_wrap_caption(anova_caption,
                                             panel_width_chars)) +
    ggplot2::theme(
      plot.caption = ggtext::element_markdown(size = 11, hjust = 0,
                                              face = "italic", colour = "grey25",
                                              margin = ggplot2::margin(t = 10)),
      axis.title.x = ggplot2::element_text(size = 13, face = "bold", color = "black")
    )

  # CLD letters per (treatment x timepoint) cell, all black, size 6, positioned
  # just above the +SE error bar tip (7% of yrange gap). Shared-letter cells at
  # the same interval render as a single label above the higher SE tip.
  # Per-cell x nudges are applied via cld_nudges (data frame: treatment, tp_num, x_adj,
  # optionally y_extra_frac, y_at_mean, suppress, force_separate).
  if (!is.null(cld_df)) {
    # force_separate: tp_num values where merging is prevented even when letters match
    force_split_tp <- if (!is.null(cld_nudges) && "force_separate" %in% names(cld_nudges))
      cld_nudges$tp_num[!is.na(cld_nudges$force_separate) & cld_nudges$force_separate]
    else integer(0)
    lab <- .mg_cld_label_rows(cld_df, smry, yrange, force_split_tp = force_split_tp,
                               collapse_tp = collapse_tp)
    if (!is.null(lab) && nrow(lab) > 0) {
      if (!is.null(cld_nudges) && nrow(cld_nudges) > 0 &&
          all(c("treatment","tp_num") %in% names(cld_nudges))) {
        lab <- merge(lab, cld_nudges, by = c("treatment","tp_num"), all.x = TRUE,
                     suffixes = c("", ".nudge"))
        # x_adj: replace default 0 with nudge value where provided
        if ("x_adj.nudge" %in% names(lab)) {
          lab$x_adj <- ifelse(is.na(lab$x_adj.nudge), lab$x_adj, lab$x_adj.nudge)
          lab$x_adj.nudge <- NULL
        }
        # y_extra_frac: additional vertical offset (fraction of yrange) for nudged letters
        if ("y_extra_frac" %in% names(lab)) {
          lab$y <- lab$y + ifelse(is.na(lab$y_extra_frac), 0, lab$y_extra_frac * yrange)
          lab$y_extra_frac <- NULL
        }
        # y_at_mean: override y to the treatment x timepoint cell mean (not above SE bar)
        if ("y_at_mean" %in% names(lab)) {
          need_y <- !is.na(lab$y_at_mean) & lab$y_at_mean & !is.na(lab$treatment)
          if (any(need_y)) {
            mean_lkp <- stats::setNames(smry$mean_y,
                                        paste(smry$treatment, smry$tp_num))
            for (i in which(need_y)) {
              mv <- mean_lkp[paste(lab$treatment[i], lab$tp_num[i])]
              if (!is.na(mv) && is.finite(mv)) lab$y[i] <- mv
            }
          }
          lab$y_at_mean <- NULL
        }
        # y_at_treat_avg: set y to grand mean of ALL treatment cell means at that
        # tp_num (from smry) plus the standard 0.11*yrange offset.  Used for a
        # single shared label positioned vertically at the midpoint between groups.
        if ("y_at_treat_avg" %in% names(lab)) {
          need_avg <- !is.na(lab$y_at_treat_avg) & lab$y_at_treat_avg
          if (any(need_avg)) {
            for (i in which(need_avg)) {
              tp_means <- smry$mean_y[smry$tp_num == lab$tp_num[i]]
              tp_means <- tp_means[is.finite(tp_means)]
              if (length(tp_means) >= 1)
                lab$y[i] <- mean(tp_means) + 0.11 * yrange
            }
          }
          lab$y_at_treat_avg <- NULL
        }
        # y_centered_gap: centre the letter vertically in the empty band between
        # the two treatments' error bars at that interval (own SE tip facing the
        # other group, to the other group's facing SE tip). Used where the
        # default "above my own error bar" position would collide with the other
        # treatment's marker.
        if ("y_centered_gap" %in% names(lab)) {
          need_gap <- !is.na(lab$y_centered_gap) & lab$y_centered_gap &
                      !is.na(lab$treatment)
          if (any(need_gap)) {
            for (i in which(need_gap)) {
              sm_tp <- smry[smry$tp_num == lab$tp_num[i], , drop = FALSE]
              own   <- sm_tp[sm_tp$treatment == lab$treatment[i], , drop = FALSE]
              oth   <- sm_tp[sm_tp$treatment != lab$treatment[i], , drop = FALSE]
              if (nrow(own) == 1L && nrow(oth) == 1L) {
                if (oth$mean_y > own$mean_y) {
                  lab$y[i] <- mean(c(own$mean_y + own$sem_y,
                                     oth$mean_y - oth$sem_y))
                } else {
                  lab$y[i] <- mean(c(own$mean_y - own$sem_y,
                                     oth$mean_y + oth$sem_y))
                }
              }
            }
          }
          lab$y_centered_gap <- NULL
        }
        # label_override: replace computed label text for specific cells
        if ("label_override" %in% names(lab)) {
          has_ovr <- !is.na(lab$label_override) & nzchar(lab$label_override)
          if (any(has_ovr)) lab$label[has_ovr] <- lab$label_override[has_ovr]
          lab$label_override <- NULL
        }
        # suppress: drop specific cells from the label set
        if ("suppress" %in% names(lab)) {
          lab <- lab[is.na(lab$suppress) | !lab$suppress, , drop = FALSE]
          if ("suppress" %in% names(lab)) lab$suppress <- NULL
        }
        # force_separate: already consumed above; remove from lab if present
        if ("force_separate" %in% names(lab)) lab$force_separate <- NULL
      }
      lab_split <- !lab$merged
      if (any(lab_split)) {
        p <- p + ggplot2::geom_text(
          data = lab[lab_split, , drop = FALSE],
          ggplot2::aes(x = tp_num + x_adj, y = y, label = label),
          inherit.aes = FALSE, fontface = "bold", size = 6,
          colour = "black", show.legend = FALSE
        )
      }
      if (any(lab$merged)) {
        p <- p + ggplot2::geom_text(
          data = lab[lab$merged, , drop = FALSE],
          ggplot2::aes(x = tp_num + x_adj, y = y, label = label),
          inherit.aes = FALSE, fontface = "bold", size = 6,
          colour = "black"
        )
      }
      # Expand top margin: base offset 11% + max nudge extra 4% + buffer → 24%.
      p <- p + ggplot2::coord_cartesian(
        ylim = c(ymin - 0.05 * yrange, ymax + 0.24 * yrange))
    }
  }
  p
}

# Build a short ANOVA subtitle for each panel: Treatment × Interval interaction
# term only, with full-word labels to avoid abbreviation in the figure.
# Wrapping is applied at render time via .mg_wrap_caption(panel_width_chars).
.mg_anova_cap <- function(res_tp_obj, label_main_term) {
  if (is.null(res_tp_obj) || is.null(res_tp_obj$anova)) return(NULL)
  av <- res_tp_obj$anova
  row <- av[grepl("treatment.*timepoint|timepoint.*treatment", av$term,
                   ignore.case = TRUE, perl = TRUE), ]
  if (nrow(row) == 0) return(NULL)
  p <- suppressWarnings(as.numeric(row$p_value[1]))
  if (!is.finite(p)) return(NULL)
  p_str <- fmt_p(p)
  df1 <- if ("df" %in% names(row)) as.integer(round(row$df[1])) else NA_integer_
  df2_raw <- if ("df_denom" %in% names(row)) as.numeric(row$df_denom[1]) else NA_real_
  chi <- fmt_F(as.numeric(row$chisq[1]))
  if (is.finite(df2_raw)) {
    eta <- fmt3((as.numeric(row$chisq[1]) * df1) / (as.numeric(row$chisq[1]) * df1 + df2_raw))
    sprintf("Treatment × Interval: F<sub>%d, %s</sub> = %s, %s, &eta;<sup>2</sup><sub>p</sub> = %s",
            df1, fmt_F(df2_raw), chi, p_str, eta)
  } else {
    sprintf("Treatment × Interval: F<sub>%d</sub> = %s, %s", df1, chi, p_str)
  }
}

# ---------------------------------------------------------------------------
# Figure_10 ANOVA captions — sourced from anova_statement_table_interval.csv.
# The helper returns the "statement" string for a given indicator label,
# falling back to NULL (which suppresses the caption) if the file or row
# is missing.
.fig10_cap_csv <- local({
  .csv_path <- file.path(PROJECT_ROOT, "choice R pipeline/easy_scripts/standalone/anova_statement_table_interval.csv")
  .tbl <- tryCatch(read.csv(.csv_path, stringsAsFactors = FALSE, encoding = "UTF-8"),
                   error = function(e) NULL)
  function(indicator) {
    if (is.null(.tbl)) return(NULL)
    row <- .tbl[trimws(.tbl$indicator) == indicator, , drop = FALSE]
    if (nrow(row) == 0) return(NULL)
    .sub_to_html(trimws(row$statement[1]))
  }
})
# ---------------------------------------------------------------------------

# ---- Figure 7bis panel: alr(flow) by treatment x interval (single panel) ----
# Built here (moved up from its original standalone-only position, 2026-08-09)
# so the same native ggplot object can double as Figure_8 panel F below, in
# addition to still being saved on its own further down. Already
# .mg_frame()-wrapped inside .mg_make_line_with_tukey()'s caller just like
# .fig7_A-.fig7_D, so it needs no raster embedding and no scale compensation
# to match its neighbours -- same theme, same wrapper, same native panel path.
.fig7bis_cap <- .fig7_cells_cap_csv("zone_main_timepoint") %||%
                (if (exists("res_zone_main")) .mg_anova_cap(res_zone_main) else NULL)
.fig7bis <- if (!is.null(df_main_wide) && "logit_flow" %in% names(df_main_wide)) {
  .mg_frame(.mg_make_line_with_tukey(
    df_src        = df_main_wide,
    y_col         = "logit_flow",
    y_label       = "ALR(flow vs calm)",
    label_dir     = "zone_main_timepoint",
    anova_caption = .fig7bis_cap,
    remap_cld     = TRUE,
    # Control at interval 2 sits well below exercise choice; its default label
    # position (own SE tip + 11% of range) landed on the exercise-choice marker.
    # Centre it in the gap between the two groups' error bars instead.
    cld_nudges    = data.frame(
      treatment      = "control",
      tp_num         = 2L,
      y_centered_gap = TRUE,
      stringsAsFactors = FALSE)
  ))
} else NULL

# ---- Legend images for Figures 8–10 (external PNG files) ----------------
# Loaded once here; re-used by Figure 9 and Figure 10 assembly below.
# Hard-stop if either file is missing — do not silently substitute.
{
  .leg_dir  <- file.path(PROJECT_ROOT, "manu_graphs_26.05.2026/Legends")
  .ab_path  <- file.path(.leg_dir, "Figure_10_legend_panelAB.PNG")
  .trt_path <- file.path(.leg_dir, "treatment_legend.PNG")
  for (.p in c(.ab_path, .trt_path)) {
    if (!file.exists(.p))
      stop("[Figure assembly] Legend file not found: ", .p,
           "\n  Place the file there and re-run.")
  }
  # Wrap a PNG as a ggplot panel that preserves the original pixel aspect ratio.
  # coord_equal() enforces asp = px_height/px_width; any extra cell space is
  # filled by the white plot background — no stretching of legend text/symbols.
  .gg_legend <- function(img) {
    asp <- dim(img)[1] / dim(img)[2]        # px_height / px_width
    ggplot2::ggplot() +
      ggplot2::annotation_custom(
        grid::rasterGrob(img, interpolate = TRUE),
        xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
      ) +
      ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
      ggplot2::scale_y_continuous(limits = c(0, asp), expand = c(0, 0)) +
      ggplot2::coord_equal() +
      ggplot2::theme_void() +
      ggplot2::theme(
        plot.background  = ggplot2::element_rect(fill = "white", colour = NA),
        plot.margin      = ggplot2::margin(0, 0, 0, 0)
      )
  }
  .legend_ab_elem  <- .gg_legend(png::readPNG(.ab_path))
  .legend_trt_elem <- .gg_legend(png::readPNG(.trt_path))
  ts_msg("[legends] loaded: Figure_10_legend_panelAB.PNG + treatment_legend.PNG")
}

# Panel E: p_stay_Flow (raster-embedded, sequence-module engine). DESCRIPTIVE
# by the ALR collinearity screen (r = 0.87 vs ALR(flow); collinearity_vs_ALR.R)
# -- reported here, not as independent confirmation, and captioned accordingly
# (author decision, 2026-08-09; see the manuscript's second-level collinearity
# screen, Methods).
# Panel F: Figure 7bis (alr(flow) by treatment x interval), reusing the native
# ggplot object built above -- same theme/wrapper as A-D, so no raster
# embedding or scale compensation is needed for this panel. (Author decision,
# 2026-08-09: keep panel F as Figure 7bis, not p_stay_High -- reconciled here
# and in figures_step5_module.R to match.)
.seq_fig_dir     <- file.path(PROJECT_ROOT, "all_manu_graphs/sequence_and_patterns_graphs")
.pstay_flow_path <- file.path(.seq_fig_dir, "seq_metric_binary_p_stay_Flow.png")
if (!file.exists(.pstay_flow_path))
  stop("[Figure assembly] Sequence-module figure not found: ", .pstay_flow_path,
       "\n  Run figures_seq_bout_module.R first.")
.fig7_E <- .mg_frame(.gg_legend(png::readPNG(.pstay_flow_path)))
.fig7_F <- .fig7bis

if (!any(sapply(list(.fig7_A, .fig7_B, .fig7_C, .fig7_D, .fig7_E, .fig7_F), is.null))) {
  # Legend-free panels + panel tags
  .f7_body_A <- .fig7_A + ggplot2::theme(legend.position = "none") + .MG_TAG_THEME + ggplot2::labs(tag = "A")
  .f7_body_B <- .fig7_B + ggplot2::theme(legend.position = "none") + .MG_TAG_THEME + ggplot2::labs(tag = "B")
  .f7_body_C <- .fig7_C + ggplot2::theme(legend.position = "none") + .MG_TAG_THEME + ggplot2::labs(tag = "C")
  .f7_body_D <- .fig7_D + ggplot2::theme(legend.position = "none") + .MG_TAG_THEME + ggplot2::labs(tag = "D")
  .f7_body_E <- .fig7_E + .MG_TAG_THEME + ggplot2::labs(tag = "E")
  .f7_body_F <- .fig7_F + ggplot2::theme(legend.position = "none") + .MG_TAG_THEME + ggplot2::labs(tag = "F")
  # Single shared legend centred below the lowest panel row (E/F), matching
  # Figure 9's disposition. Shape encodes treatment in A/B/E/F; colour maps to
  # treatment in A/B/E/F but to zone in C/D, so only the shape (treatment)
  # legend is shared across the whole grid.
  .f7_top <- (patchwork::wrap_elements(full = .f7_body_A) |
               patchwork::wrap_elements(full = .f7_body_B))
  .f7_mid <- (patchwork::wrap_elements(full = .f7_body_C) |
               patchwork::wrap_elements(full = .f7_body_D))
  .f7_bot <- (patchwork::wrap_elements(full = .f7_body_E) |
               patchwork::wrap_elements(full = .f7_body_F))
  .fig7 <- (patchwork::wrap_elements(full = .f7_top) /
             patchwork::wrap_elements(full = .f7_mid) /
             patchwork::wrap_elements(full = .f7_bot) /
             .legend_trt_elem) +
             patchwork::plot_layout(heights = c(6, 6, 6, 0.5)) +
             patchwork::plot_annotation(
               theme = ggplot2::theme(
                 plot.background = ggplot2::element_rect(fill = "white", colour = NA),
                 plot.margin     = ggplot2::margin(5, 5, 5, 5, unit = "mm")
               )
             )
  .mg_save(.fig7, "Figure_8", width_mm = 280, height_mm = 320)
} else {
  ts_msg("  Figure 8 skipped — missing component plot(s)")
}
writeLines(c(
  "Figure 7. Zone preference by treatment (N = 8 trials per treatment).",
  "",
  "(A) Additive log-ratio of flow vs calm occupancy (alr(flow)) per trial.",
  "    Points jittered around treatment group. Crossbar = mean; error bar = +/- 1 SE.",
  "    Asterisk (*) indicates significant treatment effect (LMM, Kenward-Roger).",
  "",
  "(B) Additive log-ratio of high sub-zone vs calm occupancy (alr(high)) per trial.",
  "    As panel A.",
  "",
  "(C) Cell-mean proportion of observation time in each main zone (flow, calm)",
  "    by treatment. Cells from a beta generalised linear mixed model with logit",
  "    link (p ~ treatment * zone + (1|trial)). Horizontal dashed line at 0.5",
  "    indicates equal flow:calm occupancy. CLD letters denote Tukey-adjusted",
  "    pairwise cell groupings (all cells group a: no pairwise differences).",
  "",
  "(D) Cell-mean centred log-ratio (CLR) of the area-normalised 4-part composition",
  "    (high, medium, low, calm) by treatment. Cells from a Gaussian linear mixed",
  "    model (clr ~ treatment * zone + (1|trial) + (1|zone:trial)). Horizontal",
  "    dashed line at 0 indicates use proportional to the geometric mean of the",
  "    composition. CLD letters denote Tukey-adjusted pairwise cell groupings.",
  "",
  "(E) Probability the school remains in the flow zone from one 1-s bin to the",
  "    next (binary alphabet), per trial. Gaussian LMM, treatment fixed effect.",
  "    Sourced from the sequence (bin-state) analysis, not refitted here.",
  "    DESCRIPTIVE: shares 76% of its school-level variance with ALR(flow)",
  "    (r = 0.87; Table S19) and is reported as a re-expression of that result,",
  "    not independent confirmation.",
  "",
  "(F) Additive log-ratio of flow vs calm occupancy (alr(flow)) by treatment",
  "    across observation intervals -- the interval-resolved view of panel A",
  "    (N = 8 trials per treatment per interval; intervals: 5-25 min, 45-65 min,",
  "    85-105 min). Lines connect treatment group means; error bars = +/- 1 SE;",
  "    points are trial-level values. CLD letters denote Tukey-adjusted pairwise",
  "    groupings across treatment x interval cells.",
  "",
  "Bars show cell mean +/- 1 SE; points are trial-level values (jittered).",
  "ANOVA captions emitted per panel only when the treatment, treatment x zone,",
  "or treatment x interval term is significant (p < 0.05; 3 d.p.). Cell-mean",
  "analyses (C, D) are exploratory — Tukey-adjusted within-model only; no BH",
  "correction applied."
), file.path(.mg_dir, "figure7_caption.txt"))
ts_msg("  Figure 7 caption written.")

# ---- Figure S7: Jacobs' D supplementary panels (1 x 2: main + sub) ---------
.figS7_A <- .mg_frame(.mg_make_jacobs_panel(
  jac_df     = if (exists(".jac_main_agg")) .jac_main_agg else NULL,
  zones_keep = c("flow", "calm"),
  res_list   = list(res_jac_beta_flow_agg),
  n_cap_df   = df_main_wide_agg,
  label_dir  = "jacobs_beta_main_joint_agg"))
.figS7_B <- .mg_frame(.mg_make_jacobs_panel(
  jac_df     = if (exists(".jac_sec_agg")) .jac_sec_agg else NULL,
  zones_keep = c("high", "medium", "low", "calm"),
  res_list   = list(res_jac_beta_high_agg, res_jac_beta_med_agg,
                        res_jac_beta_low_agg),
  n_cap_df  = df_sec_wide_agg,
  label_dir = "jacobs_beta_sec_joint_agg"))

if (!any(sapply(list(.figS7_A, .figS7_B), is.null))) {
  .figS7 <- (.figS7_A | .figS7_B) +
    patchwork::plot_annotation(tag_levels = list(c("A", "B"))) +
    patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(legend.position = "bottom", legend.box = "horizontal") &
    .MG_TAG_THEME
  .mg_save(.figS7, "figureS7_jacobs_D", width_mm = 280, height_mm = 130)
} else {
  ts_msg("  Figure S7 skipped — missing component plot(s)")
}
writeLines(c(
  "Figure S7. Jacobs' preference index D as an interpretive companion to",
  "Figure 7 (N = 8 trials per treatment).",
  "",
  "(A) Jacobs' D for main zones (flow, calm) under each treatment. Per-zone",
  "    beta-GLMMs on (D + 1)/2 with Smithson-Verkuilen squeeze. No adjustment",
  "    for multiple comparisons; p is raw.",
  "",
  "(B) Jacobs' D for sub-zones (high, medium, low, calm) under each treatment.",
  "    Per-zone beta-GLMMs. No adjustment for multiple comparisons; p is raw.",
  "",
  "Bars show cell mean +/- 1 SE. Points are trial-level values (jittered).",
  "Horizontal dashed line at D = 0 indicates use proportional to availability",
  "(no preference). Jacobs' D is a familiar ethology preference index that",
  "folds availability into the metric; treated here as an interpretive",
  "companion to the cell-mean log-ratio inference shown in Figure 7."
), file.path(.mg_dir, "figureS7_caption.txt"))
ts_msg("  Figure S7 caption written.")

# ---- Figure 8 — Trial-level collective + transitions (6 panels A-F, 3x2) ----
# Panel A: flow-calm transitions (moved here from standalone Figure 9)
# Panels B-F: NND, IID, Convex hull, Centroid speed
# CLD letters instead of asterisks; ANOVA statement + N label per panel.

# ANOVA caption formatter for Figure 8: emits F / χ² stat only when p < 0.05;
# p formatted to exactly 3 decimal places (e.g., p = 0.019, p < 0.001).
.fig8_anova_cap <- function(res) {
  if (is.null(res) || is.null(res$anova)) return(NULL)
  av  <- res$anova
  row <- av[grepl("^treatment$", av$term, ignore.case = TRUE, perl = TRUE), ]
  if (nrow(row) == 0L) return(NULL)
  row <- row[1L, ]
  p   <- suppressWarnings(as.numeric(row$p_value))
  if (!is.finite(p)) return(NULL)
  p_str  <- fmt_p(p)
  chi    <- fmt_F(as.numeric(row$chisq))
  st     <- attr(av, "stat_type") %||% "Wald-chisq"
  is_f   <- isTRUE(grepl("^F", st))
  if (is_f) {
    df1 <- as.integer(row$df)
    df2 <- suppressWarnings(as.numeric(
             if (!is.null(row$df_denom)) row$df_denom else NA_real_))
    if (is.finite(df2)) {
      eta <- fmt_es3((as.numeric(row$chisq) * df1) / (as.numeric(row$chisq) * df1 + df2))
      return(sprintf("F<sub>%d, %s</sub> = %s, %s, &eta;<sup>2</sup><sub>p</sub> = %s",
                     df1, fmt_F(df2), chi, p_str, eta))
    } else
      return(sprintf("F<sub>%d</sub> = %s, %s", df1, chi, p_str))
  } else {
    return(sprintf("χ²<sub>%d</sub> = %s, %s", as.integer(row$df), chi, p_str))
  }
}

# Figure 9 only: split the one-line ANOVA statement so the F/df/p part fills
# the first line and the effect size drops to a second. The six panels of this
# figure are the narrowest in the manuscript; on one line the caption had to
# run at 7.5 pt to fit, which is below the journal's minimum readable size.
# Every other figure keeps the single-line form (author decision, 2026-08-09).
# Two spellings of the effect size reach here: the caption CSV carries a literal
# U+03B7 ("η²ₚ"), the in-script fallback .fig8_anova_cap() emits the HTML entity
# "&eta;". Match either so the break lands regardless of the caption's source.
.mg_break_before_eta <- function(txt) {
  if (is.null(txt) || is.na(txt) || !nzchar(txt)) return(txt)
  sub(",\\s*(&eta;|η)", ",<br>\\1", txt)
}

.mg_mk_collective_panel <- function(res, y_col, y_label, dat = NULL,
                                     y_limits = NULL,
                                     show_colour_legend = TRUE,
                                     show_n = TRUE,
                                     cap_csv = NULL,
                                     sig_from_csv = NULL) {
  .dat <- if (!is.null(dat)) dat else {
    if (!is.null(df_gd_agg) && y_col %in% names(df_gd_agg)) df_gd_agg else NULL
  }
  if (is.null(.dat)) return(NULL)
  tryCatch({
    .lvls <- levels(factor(.dat$treatment))
    .tmap <- stats::setNames(seq_along(.lvls), .lvls)
    # Determine significance for bracket: use CSV-parsed p-value when provided,
    # falling back to the LMM result object. This ensures the bracket matches
    # the caption statement from anova_statement_table_trial.csv.
    .is_sig_val <- if (!is.null(sig_from_csv)) sig_from_csv else .is_sig(res, "Treatment")
    p <- .make_std_plot(
      df          = .dat,
      y_col       = y_col,
      y_label     = y_label,
      title       = "",
      caption_txt = .mg_break_before_eta(
                      if (!is.null(cap_csv)) cap_csv else .fig8_anova_cap(res)),
      cld_df      = NULL,
      is_sig      = .is_sig_val,
      facet_tp    = FALSE,
      n_txt       = NULL,
      sig_textsize = 8
    )
    if (show_n) p <- .mg_add_n_labels(p, .dat, y_col, .tmap)
    # When show_colour_legend = FALSE, suppress colour from the legend so only
    # shape remains (colours are still rendered on the points themselves).
    if (!show_colour_legend)
      p <- p + ggplot2::scale_colour_manual(
                 values = TREATMENT_COLORS_PAL, name = "Treatment", guide = "none")
    if (!is.null(y_limits))
      p <- p + ggplot2::coord_cartesian(ylim = y_limits)
    .mg_frame(p)
  }, error = function(e) {
    warning("Figure 8 panel '", y_col, "' failed: ", e$message); NULL
  })
}

# PANEL A REPLACED 2026-08-09: switches_per_session (the STEP5 collective-
# metrics engine's own raw flow<->calm transition count, F(1,8) = 0.105,
# p = 0.754, ns) is swapped for switch_rate (the sequence engine's per-minute
# transition RATE from the 1-s bin state series, F(1,10) = 59.67, p < 0.001,
# eta2p = 0.856 -- see seq_anova_trial_level.csv). Same underlying behaviour
# (flow<->calm switching), a materially stronger and better-normalised
# measurement of it. Raster-embedded from figures_seq_bout_module.R's
# fig_seq_manuscript_embeds() for the same reason as Figure 7 panels E/F: the
# sequence engine's own model is the source of truth, never re-fitted here.
# Panel F is the BINARY alphabet entropy rate (F1,10 = 68.11, p < 0.001,
# eta2p = 0.872), matching panel A's binary switch_rate -- NOT the graded
# version (an intermediate edit had this pointed at the graded alphabet;
# corrected to match the approved plan).
.seq_fig_dir2      <- file.path(PROJECT_ROOT, "all_manu_graphs/sequence_and_patterns_graphs")
.switch_rate_path  <- file.path(.seq_fig_dir2, "seq_switch_rate_manuscript.png")
.entropy_bin_path  <- file.path(.seq_fig_dir2, "seq_metric_binary_entropy_rate.png")
for (.p in c(.switch_rate_path, .entropy_bin_path)) {
  if (!file.exists(.p))
    stop("[Figure assembly] Sequence-module figure not found: ", .p,
         "\n  Run figures_seq_bout_module.R first.")
}
.fig8_A <- .mg_frame(.gg_legend(png::readPNG(.switch_rate_path)))
.fig8_F <- .mg_frame(.gg_legend(png::readPNG(.entropy_bin_path)))
.fig8_B <- .mg_treatment_xaxis(.mg_mk_collective_panel(res_nnd_agg,
  "mean_nnd_cm",          "Mean NND (cm)",                                  show_n = FALSE,
  cap_csv = .fig9_cap_csv("mean_nnd_cm"),
  sig_from_csv = .csv_is_sig(.fig9_cap_csv("mean_nnd_cm"))))
# Every indicator in this figure (transitions, NND, IID, school area, school
# speed) is a permutation-invariant function of the per-frame position set and
# is therefore immune to identity error. Alignment/heading-based measures are
# not; none is fitted or plotted anywhere in this pipeline (removed 2026-08-16,
# see METHODS_CHANGES.md section 7).
# Panels renumbered 2026-08-09: old D/E/F -> new C/D/E, giving a 5-panel figure.
.fig8_C <- .mg_treatment_xaxis(.mg_mk_collective_panel(res_iid_agg,
  "mean_iid_cm",          "Mean IID (cm)",          show_colour_legend = FALSE, show_n = FALSE,
  cap_csv = .fig9_cap_csv("mean_iid_cm"),
  sig_from_csv = .csv_is_sig(.fig9_cap_csv("mean_iid_cm"))))
.fig8_D <- .mg_treatment_xaxis(.mg_mk_collective_panel(res_hull_agg,
  "mean_hull_area_cm2",   "School area (cm²)",       show_colour_legend = FALSE, show_n = FALSE,
  cap_csv = .fig9_cap_csv("mean_school_area_cm2"),
  sig_from_csv = .csv_is_sig(.fig9_cap_csv("mean_school_area_cm2"))))
.fig8_E <- .mg_treatment_xaxis(.mg_mk_collective_panel(res_cspd_agg,
  "mean_centroid_spd_cm", "School speed (cm/s)",
  show_colour_legend = FALSE, show_n = FALSE,
  cap_csv = .fig9_cap_csv("mean_school_speed_cm_s"),
  sig_from_csv = .csv_is_sig(.fig9_cap_csv("mean_school_speed_cm_s"))))

if (!any(sapply(list(.fig8_A, .fig8_B, .fig8_C, .fig8_D, .fig8_E, .fig8_F), is.null))) {
  # Panel tag theme for Figure 9: tags placed in the left margin (further left
  # than "topleft" which anchors to the panel boundary).
  .f8_tag <- .MG_TAG_THEME + ggplot2::theme(plot.tag.position = c(0, 1.0))
  # Caption size reduced from BASE_THEME (14 pt) so all 5 captions fit within
  # each narrow panel without truncation. Raised 7.5 -> 9 on 2026-08-09 once
  # .mg_break_before_eta() moved the effect size onto a second line: line 1 is
  # then at most 35 characters ("Treatment F1,8 = 0.1048, p = 0.754,") vs the
  # 47-character single-line string that was CONFIRMED to fit at 7.5 pt in
  # this exact layout, so 7.5 * (47/35) * 0.9 safety margin = 9.06 -> 9 pt.
  # (11 pt was tried first and clipped panel A/B on the rendered PNG -- do
  # not re-attempt a value above ~9.5 without re-verifying on real output;
  # a full STEP5 re-run costs ~50 min, so this is deliberately conservative.)
  .f8_cap_theme <- ggplot2::theme(
    plot.caption = ggtext::element_markdown(size = 9, hjust = 0, face = "italic",
                                            lineheight = 1.3,
                                            margin = ggplot2::margin(t = 10))
  )
  # Legend-free panels + tags
  .f8_body_A <- .fig8_A + ggplot2::theme(legend.position = "none") + .f8_tag + .f8_cap_theme + ggplot2::labs(tag = "A")
  .f8_body_B <- .fig8_B + ggplot2::theme(legend.position = "none") + .f8_tag + .f8_cap_theme + ggplot2::labs(tag = "B")
  .f8_body_C <- .fig8_C + ggplot2::theme(legend.position = "none") + .f8_tag + .f8_cap_theme + ggplot2::labs(tag = "C")
  # Bottom-row panels: halve bottom margin (10→5 mm) to reduce gap to legend
  .f8_body_D <- .fig8_D + ggplot2::theme(legend.position = "none",
    plot.margin = ggplot2::margin(10, 10, 5, 10, unit = "mm")) + .f8_tag + .f8_cap_theme + ggplot2::labs(tag = "D")
  .f8_body_E <- .fig8_E + ggplot2::theme(legend.position = "none",
    plot.margin = ggplot2::margin(10, 10, 5, 10, unit = "mm")) + .f8_tag + .f8_cap_theme + ggplot2::labs(tag = "E")
  # Panel F (2026-08-09): entropy_rate, graded alphabet, raster-embedded from
  # the sequence module -- fills the slot left empty when the fifth collective
  # old panel C) was removed, restoring the 3x2 grid.
  .f8_body_F <- .fig8_F + ggplot2::theme(legend.position = "none",
    plot.margin = ggplot2::margin(10, 10, 5, 10, unit = "mm")) + .f8_tag + .f8_cap_theme + ggplot2::labs(tag = "F")
  # 3x2 panel grid (restored 2026-08-09 -- see panel F above).
  .f8_grid <- (patchwork::wrap_elements(full = .f8_body_A) |
               patchwork::wrap_elements(full = .f8_body_B) |
               patchwork::wrap_elements(full = .f8_body_C)) /
              (patchwork::wrap_elements(full = .f8_body_D) |
               patchwork::wrap_elements(full = .f8_body_E) |
               patchwork::wrap_elements(full = .f8_body_F))
  # Figure_10_legend_panelAB centred below grid (0.5x reduced height); 0.5 cm white frame
  .fig8 <- (patchwork::wrap_elements(full = .f8_grid) /
             .legend_ab_elem) +
             patchwork::plot_layout(heights = c(6, 0.5)) +
             patchwork::plot_annotation(
               theme = ggplot2::theme(
                 plot.background = ggplot2::element_rect(fill = "white", colour = NA),
                 plot.margin     = ggplot2::margin(5, 5, 5, 5, unit = "mm")
               )
             )
  .mg_save(.fig8, "Figure_9", width_mm = 320, height_mm = 220)
} else {
  ts_msg("  Figure 9 skipped — missing component plot(s)")
}
writeLines(c(
  "Figure 9. Trial-level collective movement and transition indicators",
  "(N = 8 trials per treatment).",
  "",
  "(A) School transition rate between main zones (switch_rate; transitions",
  "    per minute of the flow/calm bin-state sequence), per trial. Gaussian",
  "    LMM, treatment fixed effect. Sourced from the sequence (bin-state)",
  "    analysis (behavioural_sequence_analysis_choice_exp.R), not refitted",
  "    here; replaces the earlier switches_per_session count (F(1,8) = 0.105,",
  "    p = 0.754, ns), a cruder, non-rate-normalised measure of the same",
  "    behaviour. DESCRIPTIVE: shares 74% of its school-level variance with",
  "    ALR(flow) (r = -0.86; Table S19).",
  "(B) Mean nearest-neighbour distance (NND, cm). Gaussian LMM.",
  "(C) Mean inter-individual distance (IID, cm). Gaussian LMM.",
  "(D) School area (cm2). Gaussian LMM.",
  "(E) School speed (cm/s). Gaussian LMM.",
  "(F) Normalised entropy rate of the binary (Flow/Calm) engagement-state",
  "    sequence (0-1), per trial. Gaussian LMM, treatment fixed effect.",
  "    Sourced from the sequence analysis, as panel A. DESCRIPTIVE: shares",
  "    78% of its school-level variance with ALR(flow) (r = -0.88; Table S19).",
  "",
  "Bars show group mean +/- 1 SE. Points are raw trial-level values (jittered).",
  "CLD letters denote Tukey-adjusted pairwise groupings (max 2 characters).",
  "p-values from LMMs with treatment as fixed effect and fish density as",
  "random effect (RE-selection by AICc); Kenward-Roger df correction applied.",
  "* p < 0.05; ** p < 0.01; *** p < 0.001."
), file.path(.mg_dir, "figure8_caption.txt"))
ts_msg("  Figure 9 caption written.")

# ---- Figure 10 — Interval-resolved indicators with Tukey brackets ----------
# (internal object names below still say "fig9"; manuscript calls this Figure 10)
# Panels: A NND, B IID, C Hull. Centroid speed dropped
# (no significant ANOVA term and no significant Tukey contrast).
# Each panel is a line plot (mean ± SE per timepoint × treatment); significant
# Tukey contrasts (p < 0.05) are annotated with brackets and asterisks.

# (helper function bundle moved up above -- see 'Line-with-Tukey panel helpers' near FIGURE 7/8 section)
.fig9_A <- if (exists("res_nnd_tp")) .mg_make_line_with_tukey(
  df_src        = if (!is.null(df_gd) && "mean_nnd_cm" %in% names(df_gd)) df_gd else NULL,
  y_col         = "mean_nnd_cm",
  y_label       = "Mean NND (cm)",
  label_dir     = "nnd_timepoint",
  anova_caption = .fig10_cap_csv("Mean NND"),
  remap_cld     = TRUE,
  # Display-label overrides for Panel A only.
  # Positions unchanged; x_adj/y_at_mean kept for ex.choice/interval-1 only.
  cld_nudges    = data.frame(
    treatment      = c("control",  "control",  "control",
                       "exercise choice", "exercise choice", "exercise choice"),
    tp_num         = c(1L, 2L, 3L,
                       1L, 2L, 3L),
    x_adj          = c(NA_real_, NA_real_, NA_real_,
                       -0.22,    NA_real_, NA_real_),
    y_at_mean      = c(NA, NA, NA,
                       TRUE, NA, NA),
    label_override = c("a",  "a",  "a",
                       "ab", "ab", "b"),
    stringsAsFactors = FALSE)
) else NULL

# Panel removed 2026-08-09 on the same methodological ground as the Figure 9
# removal above: an alignment/heading indicator requires persistent per-fish
# identity, which the tracker cannot guarantee. The model itself was removed
# from the pipeline on 2026-08-16, so nothing is fitted or written for it.
# Panels renumbered: old C/D (IID/school area) -> new B/C.

# label_dir was NULL until 2026-08-09, which suppressed this panel's CLD letters
# even though iid_timepoint/contrasts_treatmentxtimepoint.csv holds three
# significant within-exercise-choice Tukey contrasts (intervals 1-2, 1-3 and
# 2-3, all p = 0.046) — the very contrasts the figure legend already describes.
# The letters are now drawn so the panel matches its own caption and the
# contrast table.
.fig9_B <- if (exists("res_iid_tp")) .mg_make_line_with_tukey(
  df_src        = if (!is.null(df_gd) && "mean_iid_cm" %in% names(df_gd)) df_gd else NULL,
  y_col         = "mean_iid_cm",
  y_label       = "Mean IID (cm)",
  label_dir     = "iid_timepoint",
  anova_caption = .fig10_cap_csv("Mean IID"),
  remap_cld     = TRUE
) else NULL

.fig9_C <- if (exists("res_hull_tp")) .mg_make_line_with_tukey(
  df_src        = if (!is.null(df_gd) && "mean_hull_area_cm2" %in% names(df_gd)) df_gd else NULL,
  y_col         = "mean_hull_area_cm2",
  y_label       = "School area (cm²)",
  label_dir     = "hull_area_timepoint",
  anova_caption = .fig10_cap_csv("Mean school area"),
  remap_cld     = TRUE,
  # Panel D CLD overrides:
  #   Interval 1: single shared "a" at mean(ctrl, ec) mean height + 0.11*yrange,
  #               x offset left. Control label suppressed; ec label carries text.
  #   Intervals 2-3: per-treatment labels at default error-bar-top positions.
  cld_nudges    = data.frame(
    treatment      = c("control",  "exercise choice",
                       "control",  "exercise choice",
                       "control",  "exercise choice"),
    tp_num         = c(1L, 1L,
                       2L, 2L,
                       3L, 3L),
    x_adj          = c(NA_real_, -0.22,
                       NA_real_, NA_real_,
                       NA_real_, NA_real_),
    force_separate = c(TRUE, TRUE,
                       TRUE, TRUE,
                       TRUE, TRUE),
    suppress       = c(TRUE, NA,
                       NA,   NA,
                       NA,   NA),
    y_at_treat_avg = c(NA, TRUE,
                       NA, NA,
                       NA, NA),
    label_override = c(NA_character_, "a",
                       "a",           "ab",
                       "a",           "b"),
    stringsAsFactors = FALSE)
) else NULL


# ---- Panel D: commitment_index by treatment x interval (bout-structure
# engine, 2026-08-09). Sourced from that engine's own session-level export
# (df_gd only holds NND/IID/area/speed) but rendered with the SAME
# .mg_make_line_with_tukey() function as A-C, so the aesthetic matches
# exactly -- same theme, point/line sizes, colour/shape scales, panel frame.
# Panel D carries CLD letters and its own Treatment x Interval ANOVA statement,
# matching panels A-C (author decision, 2026-08-10). commitment_index is flagged
# `tested = FALSE` by the bout engine's collinearity gate (school-level
# |rho| ~ 0.96/0.95 vs switch_rate/entropy_rate), but that caveat is carried in
# Methods 2.8.1 + Tables S13/S19 for EVERY collinear metric in the paper rather
# than per panel -- an earlier per-panel "DESCRIPTIVE (collinear...)" annotation
# was deliberately removed from the Figure 9 caption and Results 3.2.4, so
# annotating only this panel would be inconsistent.
# The letters come from a dedicated emmeans refit of the SAME Level-B model the
# bout engine reports (val ~ treatment * timepoint_f + fish_density + (1|trial)),
# written to easy_scripts/standalone/cld_commitment_index_treatmentxtimepoint.csv
# with letters remapped so control/interval 1 carries "a". NOTE: emmeans
# substitutes Sidak for Tukey on this 6-cell family (Tukey is only exact for a
# single set of pairwise comparisons), so these letters are Sidak-adjusted --
# disclosed in the Figure 10 caption.
# The statement is derived LIVE from the engine's own export (never hand-typed),
# so it cannot go stale if the data are refitted.
.bout_session_path <- .find_latest_csv("BOUT_output", "bout_metrics_per_session.csv")
.df_bout_session <- if (!is.null(.bout_session_path)) {
  d <- readr::read_csv(.bout_session_path, show_col_types = FALSE)
  d$treatment   <- factor(trimws(tolower(as.character(d$treatment))), levels = TREATMENT_LEVELS_g)
  d$timepoint_f <- factor(d$timepoint, levels = TIMEPOINT_LEVELS_g)
  d
} else NULL

.cld_commit_path <- file.path(file.path(PROJECT_ROOT, "choice R pipeline/easy_scripts/standalone"),
                              "cld_commitment_index_treatmentxtimepoint.csv")
.bout_tp_path <- .find_latest_csv("BOUT_output", "bout_anova_trial_x_timepoint.csv")
.commit_cap <- local({
  if (is.null(.bout_tp_path)) return("")
  d <- readr::read_csv(.bout_tp_path, show_col_types = FALSE)
  r <- d[d$metric == "commitment_index" & d$term == "treatment:timepoint_f", ]
  if (!nrow(r)) return("")
  p_txt <- if (r$p[1] < 0.001) "p < 0.001" else sprintf("p = %.3f", r$p[1])
  sprintf("Treatment\u00d7Interval *F*<sub>%s,%s</sub> = %s, %s, \u03b7\u00b2<sub>p</sub> = %s",
          format(r$df1[1]), format(round(r$df2[1], 2)),
          format(round(r$F[1], 3)), p_txt, format(round(r$eta2_p[1], 3)))
})

.fig9_D <- if (!is.null(.df_bout_session)) .mg_make_line_with_tukey(
  df_src        = .df_bout_session,
  y_col         = "commitment_index",
  y_label       = "Commitment index",
  label_dir     = NULL,
  cld_path      = .cld_commit_path,
  anova_caption = .commit_cap,
  remap_cld     = FALSE
) else NULL

.fig9_panels <- list(.fig9_A, .fig9_B, .fig9_C, .fig9_D)
.fig9_panels <- lapply(.fig9_panels, function(p) if (is.null(p)) NULL else .mg_frame(p))
if (!any(sapply(.fig9_panels, is.null))) {
  # Legend-free panels + tags; panel B (IID) retains scale overrides (suppress duplicate legend)
  .f9_body_A <- .fig9_panels[[1]] + ggplot2::theme(legend.position = "none") + .FIG10_CAP_OVR + .MG_TAG_THEME + ggplot2::labs(tag = "A")
  .f9_body_B <- .fig9_panels[[2]] + ggplot2::theme(legend.position = "none") + .FIG10_CAP_OVR +
    ggplot2::scale_colour_manual(values = TREATMENT_COLORS_PAL, name = "Treatment", guide = "none") +
    ggplot2::scale_linetype_manual(values = c("exercise choice" = "solid", "control" = "dashed"),
      name = "Treatment", guide = "none") +
    .MG_TAG_THEME + ggplot2::labs(tag = "B")
  .f9_body_C <- .fig9_panels[[3]] + ggplot2::theme(legend.position = "none") + .FIG10_CAP_OVR + .MG_TAG_THEME + ggplot2::labs(tag = "C")
  .f9_body_D <- .fig9_panels[[4]] + ggplot2::theme(legend.position = "none") + .FIG10_CAP_OVR + .MG_TAG_THEME + ggplot2::labs(tag = "D")
  # True 2x2 grid (2026-08-09): the trailing plot_spacer() that used to hold
  # this cell (kept from when the fifth collective panel was removed, to preserve
  # each panel's ~130mm width and the already-verified 11pt caption size) is
  # now filled by panel D -- same cell, same width, no re-verification needed.
  .f9_grid <- (patchwork::wrap_elements(full = .f9_body_A) |
               patchwork::wrap_elements(full = .f9_body_B)) /
              (patchwork::wrap_elements(full = .f9_body_C) |
               patchwork::wrap_elements(full = .f9_body_D))
  # Figure_10_legend_panelAB centred below grid (height decreased by 70% from 0.7); 0.5 cm white frame
  .fig9_iv <- (patchwork::wrap_elements(full = .f9_grid) /
                .legend_ab_elem) +
                patchwork::plot_layout(heights = c(6, 0.5)) +
                patchwork::plot_annotation(
                  theme = ggplot2::theme(
                    plot.background = ggplot2::element_rect(fill = "white", colour = NA),
                    plot.margin     = ggplot2::margin(5, 5, 5, 5, unit = "mm")
                  )
                )
  .mg_save(.fig9_iv, "Figure_10", width_mm = 260, height_mm = 220)
  writeLines(c(
    "Figure 10. Interval-resolved collective movement indicators with significant",
    "treatment-related terms or Tukey contrasts (N = 8 trials per treatment",
    "per interval; intervals: 5-25 min, 45-65 min, 85-105 min).",
    "",
    "(A) Mean nearest-neighbour distance (NND, cm).",
    "    Significant Treatment x Interval interaction (LMM with F-KR).",
    "(B) Mean inter-individual distance (IID, cm).",
    "    Significant within-exercise-choice Tukey contrasts across intervals.",
    "(C) School area (cm2).",
    "    Significant Tukey contrast (control Interval 2 vs exercise choice",
    "    Interval 3).",
    "(D) Commitment index (longest flow bout / interval duration), from the",
    "    bout-structure analysis. Collinear with switch_rate and entropy_rate",
    "    at school level (|r| ~ 0.98, already reported in Figure 9 panels A",
    "    and F), so it is read as a re-expression of those outcomes rather",
    "    than as independent evidence; its Treatment x Interval term and CLD",
    "    letters are shown for completeness. Letters here are Sidak-adjusted",
    "    across the six treatment x interval cells (emmeans substitutes Sidak",
    "    for Tukey on a 6-cell family); panels A-C are Tukey-adjusted.",
    "",
    "Lines show group mean +/- 1 SE per interval. Compact letter display",
    "(CLD) labels above each cell denote Tukey-adjusted pairwise groupings",
    "(treatment x interval); cells sharing a letter do not differ at p < 0.05.",
    "When both treatments at the same interval share the same letter(s),",
    "a single label is shown. Letters are truncated to a maximum of two",
    "characters. Subtitle shows the Treatment x Interval interaction",
    "F-statistic from the LMM for every panel, including panel D."
  ), file.path(.mg_dir, "figure9_caption.txt"))
  ts_msg("  Figure 10 (interval) caption written.")
} else {
  ts_msg("  Figure 10 (interval) skipped — missing component plot(s)")
}

# ---- Figure 7bis: alr(flow) by treatment × interval (single panel) ----------
# .fig7bis itself is built earlier now (with the Figure 7/8 panels, so it can
# double as Figure_8 panel F) -- this just saves it standalone and writes its
# own caption, unchanged from before.
if (!is.null(.fig7bis)) {
  .mg_save(.fig7bis, "figure7bis_alr_flow_timepoint", width_mm = 140, height_mm = 130)
  writeLines(c(
    "Figure 7bis. Additive log-ratio of flow vs calm occupancy (alr(flow)) by",
    "treatment across observation intervals (N = 8 trials per treatment per",
    "interval; intervals: 5-25 min, 45-65 min, 85-105 min).",
    "",
    "Lines connect treatment group means; error bars = +/- 1 SE. Points are",
    "trial-level values. CLD letters denote Tukey-adjusted pairwise groupings",
    "across all treatment x interval cells (precedence: control > exercise choice;",
    "interval 1 > 2 > 3). ANOVA caption emitted only when the treatment x",
    "interval interaction term is significant (p < 0.05; 3 d.p.)."
  ), file.path(.mg_dir, "figure7bis_caption.txt"))
  ts_msg("  Figure 7bis saved and caption written.")
} else {
  ts_msg("  Figure 7bis skipped — df_main_wide not available")
}

# ---- Figure 11: independent temporal structure (raster-embedded) ----------
# The 4 metrics that survive the school-level collinearity screen against the
# ALR primary endpoint (|r| < 0.70 vs every ALR outcome; collinearity_vs_ALR.R,
# Table S19): burst_flow, dwell_Flow, dwell_Calm, memory_flow -- all null at
# trial level (Level A), two significant at trial x interval (Level B):
# burst_flow's treatment main effect (p = 0.026) and dwell_Calm's treatment x
# interval interaction (p = 0.031). Built by fig11_temporal_structure() in
# figures_seq_bout_module.R (6 native ggplot panels, not re-fitted here) and
# embedded the same way as Figure 7 panels E/F. Saved under the filename
# "figure11_temporal_structure", NOT "Figure_11" -- the copy-rename step below
# already claims "Figure_11.png" for the (manuscript-numbered) Figure 13
# serotonin panel; the manuscript text and caption cite this as "Figure 11".
.fig11_path <- file.path(.seq_fig_dir, "figure11_temporal_structure.png")
if (!file.exists(.fig11_path)) {
  ts_msg("  Figure 11 skipped — ", .fig11_path,
         " not found. Run figures_seq_bout_module.R first.")
} else {
  .fig11 <- .gg_legend(png::readPNG(.fig11_path))
  .mg_save(.fig11, "figure11_temporal_structure", width_mm = 320, height_mm = 220)
  writeLines(c(
    "Figure 11. Independent temporal structure of exercise-choice behaviour",
    "(N = 8 trials/schools per treatment).",
    "",
    "(A) Burstiness of flow bouts (bout-structure analysis), per school.",
    "(B) Mean flow dwell time per visit (binary alphabet), per trial.",
    "(C) Mean calm dwell time per visit (binary alphabet), per trial.",
    "(D) Lag-1 autocorrelation of successive flow-bout durations within a",
    "    school (memory_flow).",
    "(E) Burstiness of flow bouts by interval. Significant treatment main",
    "    effect (F1,22.84 = 5.70, p = 0.026, eta2p = 0.20) despite panel A",
    "    being null at trial level -- resolvable only at session level.",
    "(F) Mean calm dwell time by interval. Significant treatment x interval",
    "    interaction (F2,24.2 = 4.01, p = 0.031, eta2p = 0.25).",
    "",
    "All four metrics clear the school-level collinearity screen against the",
    "ALR primary endpoint (|r| < 0.70 vs every ALR outcome; max pairwise |r|",
    "among the four = 0.56; Table S19) and are the independent test of",
    "temporal structure, as distinct from Figures 8 and 9's descriptive",
    "(collinear) sequence metrics. Points (A-D) are trial/school-level values",
    "(jittered); black crossbar = mean, black bar = +/- 1 SE. Lines (E-F)",
    "show interval mean +/- 1 SE. Panels A-D: one-line caption gives the",
    "treatment F-statistic. Panels E-F: two-line caption gives the treatment",
    "and treatment x interval F-statistics."
  ), file.path(.mg_dir, "figure11_temporal_structure_caption.txt"))
  ts_msg("  Figure 11 saved and caption written.")
}

ts_msg("manu_graphs folder written: ", .mg_dir)
ts_msg("Dated manu_graphs folder: ", .mg_dir_new)

# Mirror every figure from manu_graphs to D:/CHOICE R SCRIPTS/all_manu_graphs
.all_manu_dir <- file.path(PROJECT_ROOT, "all_manu_graphs")
dir.create(.all_manu_dir, recursive = TRUE, showWarnings = FALSE)
.mg_files <- list.files(.mg_dir, pattern = "\\.(png|pdf|txt)$",
                        full.names = TRUE, recursive = FALSE)
invisible(file.copy(.mg_files, file.path(.all_manu_dir, basename(.mg_files)),
                    overwrite = TRUE))
ts_msg("Mirrored ", length(.mg_files), " manu_graphs file(s) → ", .all_manu_dir)
# Copy TXT caption files to the dated PNG folder and PDF subfolder.
.mg_txt_files <- list.files(.mg_dir, pattern = "\\.txt$", full.names = TRUE, recursive = FALSE)
if (length(.mg_txt_files) && dir.exists(.mg_dir_new)) {
  invisible(file.copy(.mg_txt_files, file.path(.mg_dir_new, basename(.mg_txt_files)),
                      overwrite = TRUE))
}
if (length(.mg_txt_files) && exists(".mg_dir_new_pdf") && dir.exists(.mg_dir_new_pdf)) {
  invisible(file.copy(.mg_txt_files, file.path(.mg_dir_new_pdf, basename(.mg_txt_files)),
                      overwrite = TRUE))
}
ts_msg("Captions copied to dated folders (", length(.mg_txt_files), " txt files)")

# ---- Copy-rename Figures 11–13 from shared manu_graphs archive --------------
# Sources are written by analysis_b3.R to D:/CHOICE R SCRIPTS/all_manu_graphs/.
# We copy-rename them into the dated output folders (PNG + PDF if available).
.fig11_13_map <- list(
  list(src = "DM_serotonin_group_SEM",     dst = "Figure_11"),
  list(src = "DM_catecholamine_group_SEM", dst = "Figure_12"),
  list(src = "Cort_jitterplot_legacy",      dst = "Figure_13")
)
for (.m in .fig11_13_map) {
  # PNG → dated PNG folder
  .src_png <- file.path(.all_manu_dir, paste0(.m$src, ".png"))
  if (file.exists(.src_png) && exists(".mg_dir_new") && dir.exists(.mg_dir_new)) {
    invisible(file.copy(.src_png,
                        file.path(.mg_dir_new, paste0(.m$dst, ".png")),
                        overwrite = TRUE))
    ts_msg("  Copied ", .m$src, ".png  ->  ", .m$dst, ".png")
  } else {
    ts_msg("  WARN: source not found or output folder missing: ", .src_png)
  }
  # PDF → dated PDF subfolder: copy source PDF if available, otherwise render
  # a cairo PDF from the source PNG (rasterised at native PNG resolution).
  .src_pdf <- file.path(.all_manu_dir, paste0(.m$src, ".pdf"))
  .dst_pdf <- file.path(.mg_dir_new_pdf, paste0(.m$dst, ".pdf"))
  if (exists(".mg_dir_new_pdf") && dir.exists(.mg_dir_new_pdf)) {
    if (file.exists(.src_pdf)) {
      invisible(file.copy(.src_pdf, .dst_pdf, overwrite = TRUE))
      ts_msg("  Copied ", .m$src, ".pdf  ->  ", .m$dst, ".pdf")
    } else if (file.exists(.src_png)) {
      tryCatch({
        .img   <- png::readPNG(.src_png)
        .w_in  <- dim(.img)[2L] / 300
        .h_in  <- dim(.img)[1L] / 300
        grDevices::cairo_pdf(.dst_pdf, width = .w_in, height = .h_in)
        grid::grid.raster(.img)
        grDevices::dev.off()
        ts_msg("  Rendered ", .m$src, ".png  ->  ", .m$dst, ".pdf (cairo)")
      }, error = function(e) {
        ts_msg("  WARN: cairo PDF render failed for ", .m$dst, ": ", e$message)
      })
    }
  }
}
ts_msg("Figures 11-13 copy-rename complete.")


# =============================================================================
# ==== 7) DESCRIPTIVE SUMMARY =================================================
# =============================================================================

.sum_cols <- intersect(c("prop_active", "switches_per_session",
                          "n_main_switches", "zone_flux_per_session",
                          "prop_time_in_flow", "prop_time_in_calm",
                          "prop_time_in_high", "prop_time_in_medium",
                          "prop_time_in_low",
                          "mean_nnd_cm",
                          "mean_hull_area_cm2", "mean_centroid_spd_cm"),
                        names(df))

.grp_cols_tp   <- intersect(c("treatment","fish_density","timepoint"), names(df))
.grp_cols_marg <- intersect(c("treatment","fish_density"), names(df))

desc_tp <- df %>%
  dplyr::group_by(dplyr::across(dplyr::all_of(.grp_cols_tp))) %>%
  dplyr::summarise(n = dplyr::n(),
    dplyr::across(dplyr::all_of(.sum_cols),
      list(mean=~mean(.x,na.rm=TRUE), sd=~sd(.x,na.rm=TRUE), sem=~.sem(.x)),
      .names = "{.col}__{.fn}"),
    .groups = "drop")
readr::write_csv(desc_tp, file.path(STEP5_OUT, "descriptive_summary_by_timepoint.csv"))

desc_marg <- df %>%
  dplyr::group_by(dplyr::across(dplyr::all_of(.grp_cols_marg))) %>%
  dplyr::summarise(n = dplyr::n(),
    dplyr::across(dplyr::all_of(.sum_cols),
      list(mean=~mean(.x,na.rm=TRUE), sd=~sd(.x,na.rm=TRUE), sem=~.sem(.x)),
      .names = "{.col}__{.fn}"),
    .groups = "drop")
readr::write_csv(desc_marg, file.path(STEP5_OUT, "descriptive_summary_marginal.csv"))
ts_msg("Descriptive summaries written")


# =============================================================================
# ==== 8) WORD REPORT (officer / flextable) ===================================
# =============================================================================

ts_msg("=== WORD REPORT ===")

.rpt_path_step5 <- file.path(STEP5_OUT, "analysis_report_choice_exp.docx")
.rpt_path_main  <- file.path(
  if (exists(".pipeline_dir_choice", envir = .GlobalEnv, inherits = FALSE))
    get(".pipeline_dir_choice", envir = .GlobalEnv)
  else getwd(),
  "analysis_report_choice_exp.docx"
)

# ---- Pretty column names for flextable display --------------------------------
.pretty_names <- function(df) {
  col_map <- c(
    term             = "Term",
    df               = "df",
    chisq            = "χ² / F",
    p_value          = "p",
    sig              = "Sig.",
    re               = "RE Name",
    re_formula       = "RE Formula",
    fixed_formula    = "Fixed Formula",
    fit_method       = "Fit",
    AICc             = "AICc",
    delta_AICc       = "ΔAICc",
    selected         = "Selected",
    response         = "Response",
    transform        = "Transform",
    family_used      = "Family",
    family           = "Family",
    stat_type        = "Test",
    R2m              = "R²m",
    sw_p_raw         = "SW p (raw)",
    sw_p_trans       = "SW p (trans)",
    levene_p_raw     = "Levene p (raw)",
    levene_p_trans   = "Levene p (trans)",
    levene_resid_p   = "Levene p (resid)",
    normality_flag   = "Flag",
    re_fallback      = "RE fallback",
    set              = "BH set",
    treatment        = "Treatment",
    zone             = "Zone",
    timepoint_f      = "Timepoint",
    emmean           = "EMmean",
    SE               = "SE",
    lower.CL         = "Lower 95% CI",
    upper.CL         = "Upper 95% CI",
    `.group`         = "CLD",
    Analysis         = "Analysis",
    p_raw            = "p (raw)",
    p_BH             = "p (BH)",
    sig_raw          = "Sig. (raw)",
    sig_BH           = "Sig. (BH)",
    best_RE          = "Best RE",
    Trial_level_p    = "Trial p",
    Tank_level_p     = "Tank p",
    Trial_beta       = "Trial β",
    Tank_beta        = "Tank β",
    Beta_ratio       = "β ratio",
    Direction_agree  = "Direction agree",
    Trial_sig        = "Trial sig",
    Tank_sig         = "Tank sig",
    timepoint        = "Timepoint",
    contrast         = "Contrast",
    estimate         = "Estimate",
    Item             = "Item",
    Value            = "Value",
    n                = "n",
    mean             = "Mean",
    SEM              = "SEM"
  )
  rename_vec <- col_map[names(col_map) %in% names(df)]
  if (length(rename_vec) > 0)
    names(df)[match(names(rename_vec), names(df))] <- unname(rename_vec)
  df
}

# ---- Helper: standardise ANOVA df for flextable ----------------------------
.anova_ft <- function(res) {
  if (is.null(res) || is.null(res$anova)) return(NULL)
  av <- res$anova
  av <- av[av$term != "(Intercept)", ]
  av$chisq <- vapply(av$chisq, fmt_num, character(1))
  av$sig   <- dplyr::case_when(
    is.na(av$p_value)   ~ "",
    av$p_value < 0.001  ~ "***",
    av$p_value < 0.01   ~ "**",
    av$p_value < 0.05   ~ "*",
    TRUE                ~ ""
  )
  av$p_value <- vapply(av$p_value,
    function(p) if (is.na(p)) NA_character_ else fmt_p(p), character(1))
  av <- av[, intersect(c("term","df","chisq","p_value","sig"), names(av))]
  # Rename chisq column to match actual test used
  st <- if (!is.null(res$stat_type)) res$stat_type else "Wald-chisq"
  if ("χ² / F" %in% names(.pretty_names(av))) {
    # .pretty_names already maps chisq → "χ² / F"
  }
  .pretty_names(av)
}

.norm_ft <- function(dir_name) {
  p <- file.path(STEP5_OUT, dir_name, "normality_check.csv")
  if (!file.exists(p)) return(NULL)
  d <- readr::read_csv(p, show_col_types = FALSE)
  .pretty_names(as.data.frame(d))
}

.aicc_ft <- function(res) {
  if (is.null(res) || is.null(res$aicc_table)) return(NULL)
  tbl <- res$aicc_table
  tbl$AICc       <- round(tbl$AICc, 2)
  tbl$delta_AICc <- round(tbl$delta_AICc, 2)
  tbl <- tbl[, intersect(c("re","re_formula","fixed_formula","AICc","delta_AICc","selected"), names(tbl))]
  .pretty_names(tbl)
}

.ft_styled <- function(df, bold_col = NULL, bold_val = TRUE) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  ft <- flextable::flextable(df) %>%
    flextable::font(fontname = "Arial", part = "all") %>%
    flextable::fontsize(size = 9, part = "all") %>%
    flextable::fontsize(size = 10, part = "header") %>%
    flextable::bold(part = "header") %>%
    flextable::bg(bg = "#DEEAF1", part = "header") %>%
    flextable::border_outer(part = "all",
      border = officer::fp_border(color = "#AAAAAA", width = 0.5)) %>%
    flextable::border_inner_h(
      border = officer::fp_border(color = "#CCCCCC", width = 0.5)) %>%
    flextable::set_table_properties(layout = "autofit", width = 1)
  if (!is.null(bold_col) && bold_col %in% names(df)) {
    sel_rows <- which(df[[bold_col]] == bold_val)
    if (length(sel_rows) > 0)
      ft <- flextable::bold(ft, i = sel_rows, part = "body")
  }
  ft
}

# ---- Officer helpers --------------------------------------------------------
.add_h1  <- function(doc, txt) officer::body_add_par(doc, txt, style = "heading 1")
.add_h2  <- function(doc, txt) officer::body_add_par(doc, txt, style = "heading 2")
.add_h3  <- function(doc, txt) officer::body_add_par(doc, txt, style = "heading 3")
.add_par <- function(doc, txt) officer::body_add_par(doc, txt, style = "Normal")
.add_ft  <- function(doc, ft)  {
  if (!is.null(ft)) flextable::body_add_flextable(doc, ft)
  else officer::body_add_par(doc, "(no data)", style = "Normal")
}

# Strip Word auto-numbering from heading styles so manual numbers in text are
# not duplicated by the theme's multilevel list numbering.
# Modifies styles.xml in officer's temp package_dir before the file is written.
.strip_heading_autonumber <- function(doc) {
  tryCatch({
    pkg <- doc[["package_dir"]]
    if (!is.character(pkg) || length(pkg) != 1L || !dir.exists(pkg)) {
      ts_msg("  [headings] package_dir not found — skipping")
      return(doc)
    }
    styles_path <- file.path(pkg, "word", "styles.xml")
    if (!file.exists(styles_path)) {
      ts_msg("  [headings] styles.xml not found")
      return(doc)
    }
    xml_s <- xml2::read_xml(styles_path)
    ns <- c(w = "http://schemas.openxmlformats.org/wordprocessingml/2006/main")
    n_removed <- 0L
    for (sel in c(
      "//w:style[@w:styleId='Heading1']//w:pPr/w:numPr",
      "//w:style[@w:styleId='Heading2']//w:pPr/w:numPr",
      "//w:style[@w:styleId='Heading3']//w:pPr/w:numPr",
      "//w:style[w:name[@w:val='heading 1']]//w:pPr/w:numPr",
      "//w:style[w:name[@w:val='heading 2']]//w:pPr/w:numPr",
      "//w:style[w:name[@w:val='heading 3']]//w:pPr/w:numPr"
    )) {
      nodes <- xml2::xml_find_all(xml_s, sel, ns)
      if (length(nodes) > 0L) { xml2::xml_remove(nodes); n_removed <- n_removed + length(nodes) }
    }
    xml2::write_xml(xml_s, styles_path)
    ts_msg("  [headings] numPr nodes removed: ", n_removed)
    doc
  }, error = function(e) {
    warning("Could not strip heading autonumber: ", conditionMessage(e))
    doc
  })
}

# ---- Effect size: Cohen's d (LMM), OR (beta GLMM), IRR (Poisson/NB GLMM) ---
.calc_effect_size <- function(res) {
  if (is.null(res) || is.null(res$model)) return(NULL)
  fam <- res$family_used %||% res$family %||% "gaussian"
  tryCatch({
    cf <- tryCatch(
      lme4::fixef(res$model),
      error = function(e) tryCatch(
        glmmTMB::fixef(res$model)$cond,
        error = function(e2) NULL))
    if (is.null(cf) || length(cf) == 0) return(NULL)
    trt_idx <- grep("^treatment", names(cf), ignore.case = TRUE, perl = TRUE)[1]
    if (is.na(trt_idx)) return(NULL)
    if (fam %in% c("beta", "negbin", "poisson")) {
      coef_val <- as.numeric(cf[trt_idx])
      list(value = exp(coef_val),
           label = if (fam == "beta") "OR" else "IRR")
    } else {
      # Use the emmeans pairwise contrast (the full between-level difference),
      # not the raw contr.sum model coefficient directly. Under the sum-to-zero
      # contrasts set project-wide for Type III marginality, a 2-level factor's
      # coefficient is HALF the between-level difference (it is the deviation of
      # one level from the grand mean), so coef_val / sigma previously understated
      # Cohen's d by a factor of 2.
      sigma_val <- tryCatch(lme4::sigma(res$model), error = function(e) NA_real_)
      if (is.na(sigma_val) || sigma_val < 1e-10) return(NULL)
      em  <- tryCatch(emmeans::emmeans(res$model, ~ treatment), error = function(e) NULL)
      ctr <- if (!is.null(em)) tryCatch(as.data.frame(emmeans::contrast(em, method = "pairwise")),
                                         error = function(e) NULL) else NULL
      if (is.null(ctr) || !nrow(ctr)) return(NULL)
      list(value = ctr$estimate[1] / sigma_val, label = "Cohen's d")
    }
  }, error = function(e) NULL)
}

# ---- Per-analysis sub-report block ------------------------------------------
.analysis_block <- function(doc, section_num, label, display_name, dir_name, res,
                              is_glmm = FALSE) {
  doc <- .add_h2(doc, paste0(section_num, "  ", display_name))

  if (is.null(res)) {
    doc <- .add_par(doc, "(Analysis skipped or model did not converge.)")
    return(doc)
  }

  # Detect actual model type from res slots — is_glmm arg is kept for compat only
  .fam   <- res$family_used %||% res$family %||% "gaussian"
  .is_lmm <- .fam == "gaussian"
  .stat  <- res$stat_type %||% if (.is_lmm) "Wald-chisq" else "Wald-chisq"

  # Summary line
  .lev_str <- if (!is.na(res$levene_resid_p %||% NA_real_))
    paste0("  Levene(resid)=", fmt_num(res$levene_resid_p)) else ""
  .re_str  <- if (isTRUE(res$re_fallback)) "  [RE fallback: no valid RE, intercept-only]" else ""
  .norm_warn <- if (isTRUE(res$normality_flag))
    "\nWARNING: residual normality / homoscedasticity not fully met. Interpret with caution." else ""

  doc <- .add_par(doc, paste0(
    "Family: ", .fam,
    "  Transform: ", if (!is.null(res$transform)) res$transform else "—",
    "  Test: ", .stat,
    "  Fixed: ", if (!is.null(res$fixed_formula)) res$fixed_formula else "—",
    .lev_str, .re_str, .norm_warn))

  # Biological effect size — reported for every term regardless of significance
  # (Nakagawa & Cuthill 2007: gating on p < 0.05 discards exactly the
  # information that makes a non-significant result "inconclusive" rather than
  # "no effect").
  .trt_row_es <- if (!is.null(res$anova))
    res$anova[grepl("^treatment$", res$anova$term, perl = TRUE, ignore.case = TRUE), ]
  else NULL
  .trt_p_es <- if (!is.null(.trt_row_es) && nrow(.trt_row_es) > 0)
    suppressWarnings(as.numeric(.trt_row_es$p_value[1])) else NA_real_
  if (!is.na(.trt_p_es)) {
    .es <- .calc_effect_size(res)
    if (!is.null(.es))
      doc <- .add_par(doc, paste0(
        "Biological effect size: ", .es$label, " = ", fmt3(.es$value),
        " (", fmt_p(.trt_p_es), ")."))
    # Hedges' g for Gaussian LMMs
    .hg <- .calc_hedges_g(res, label)
    if (!is.na(.hg$g))
      doc <- .add_par(doc, paste0(
        "Hedges' g = ", .hg$g, " [95% CI: ", .hg$ci_lo, ", ", .hg$ci_hi, "] ",
        "(within-cluster variant g_resid = ", .hg$g_resid, " [", .hg$ci_lo_resid, ", ",
        .hg$ci_hi_resid, "]). ",
        "Interpretation: |g| < 0.2 trivial, 0.2–0.5 small, 0.5–0.8 medium, > 0.8 large."))
  }

  if (isTRUE(res$simplified_to_lm))
    doc <- .add_par(doc, paste0(
      "Note: All LMM random-effects candidates were singular (RE variance ≈ 0). ",
      "A fixed-effects-only model (lm) had lower AICc and was used for inference. ",
      "Random-effects structure is not present in this model."))

  # Sub-section counter: LMM starts at 1 (has normality section), GLMM at 0
  .sub <- if (.is_lmm) 1L else 0L

  # Sub-section 1: Normality & variance check (Gaussian LMM only)
  if (.is_lmm) {
    nc <- .norm_ft(dir_name)
    doc <- .add_h3(doc, paste0(section_num, ".1  Normality & Variance Check"))
    if (!is.null(nc)) {
      doc <- .add_ft(doc, .ft_styled(nc, bold_col = "Flag", bold_val = TRUE))
    } else {
      doc <- .add_par(doc, "(normality_check.csv not found)")
    }
  }

  # Sub-section 2 (LMM) / 1 (GLMM): AICc / RE selection
  doc <- .add_h3(doc, paste0(section_num, ".", .sub + 1L,
                              "  Random Effects Selection (AICc)"))
  at <- .aicc_ft(res)
  if (!is.null(at) && "fit_method" %in% names(res$aicc_table))
    at <- .pretty_names(res$aicc_table[, intersect(
      c("re","re_formula","fixed_formula","fit_method","AICc","delta_AICc","selected"),
      names(res$aicc_table))])
  doc <- .add_ft(doc, .ft_styled(at, bold_col = "Selected"))

  # Sub-section 3 (LMM) / 2 (GLMM): R²m / R²c
  doc <- .add_h3(doc, paste0(section_num, ".", .sub + 2L,
                              "  Model R² (Marginal & Conditional)"))
  if (!is.null(res$r2) && !is.null(res$r2$R2m)) {
    .r2m_val <- if (length(res$r2$R2m) >= 1) res$r2$R2m[1] else NA_real_
    .r2c_val <- if (length(res$r2$R2c) >= 1) res$r2$R2c[1] else NA_real_
    .icc_val <- res$icc %||% NA_real_
    r2_tbl <- data.frame(
      Metric = c("R²m (marginal — fixed effects only)",
                 "R²c (conditional — fixed + random effects)",
                 "ICC (adjusted, random-effect structure)"),
      Value  = c(fmt3(.r2m_val), fmt3(.r2c_val),
                 if (is.finite(.icc_val)) fmt3(.icc_val) else "—"),
      stringsAsFactors = FALSE)
    doc <- .add_ft(doc, .ft_styled(r2_tbl))
    if (isTRUE(is.finite(.r2m_val)) && isTRUE(is.finite(.r2c_val)) &&
        abs(.r2m_val - .r2c_val) < 0.005) {
      doc <- .add_par(doc, paste0(
        "Note: R²m ≈ R²c (", fmt3(.r2m_val), " vs ",
        fmt3(.r2c_val), "). Random effects explain negligible additional ",
        "variance. The mixed-effects structure was retained for consistency across ",
        "analyses; consider a fixed-effects-only model for this indicator."))
    }
  } else {
    doc <- .add_par(doc, "(R² not computed — GLMM family or performance package unavailable.)")
  }

  # Sub-section 4 (LMM) / 3 (GLMM): ANOVA
  .anova_title <- paste0(section_num, ".", .sub + 3L,
    "  ANOVA (Type III — ",
    if (.stat == "F-KR") "Kenward-Roger F"
    else if (.fam %in% c("beta","poisson","negbin")) "Wald χ²"
    else "Wald χ²", ")")
  doc <- .add_h3(doc, .anova_title)
  doc <- .add_ft(doc, .ft_styled(.anova_ft(res)))

  # KR vs Satterthwaite comparison (LMM only, when both available)
  if (isTRUE(.is_lmm) && !is.null(res$anova) && !is.null(res$anova_sw)) {
    tryCatch({
      .kr <- res$anova[, intersect(c("term","df","chisq","p_value"), names(res$anova)), drop=FALSE]
      .sw <- res$anova_sw[, intersect(c("term","df","chisq","p_value"), names(res$anova_sw)), drop=FALSE]
      names(.kr)[names(.kr) == "df"]      <- "df_KR"
      names(.kr)[names(.kr) == "chisq"]   <- "F_KR"
      names(.kr)[names(.kr) == "p_value"] <- "p_KR"
      names(.sw)[names(.sw) == "df"]      <- "df_SW"
      names(.sw)[names(.sw) == "chisq"]   <- "F_SW"
      names(.sw)[names(.sw) == "p_value"] <- "p_SW"
      .cmp <- merge(.kr, .sw, by = "term", sort = FALSE)
      .cmp$p_KR <- vapply(.cmp$p_KR, function(p)
        if (is.na(p)) NA_character_ else fmt_p(as.numeric(p)), character(1))
      .cmp$p_SW <- vapply(.cmp$p_SW, function(p)
        if (is.na(p)) NA_character_ else fmt_p(as.numeric(p)), character(1))
      .cmp[, c("F_KR","F_SW")] <- lapply(
        .cmp[, c("F_KR","F_SW")], function(x) round(as.numeric(x), 3))
      doc <- .add_par(doc, "Kenward-Roger vs Satterthwaite df comparison:")
      doc <- .add_ft(doc, .ft_styled(.cmp))
    }, error = function(e)
      doc <<- .add_par(doc, paste0("(KR/SW comparison failed: ", e$message, ")")))
  }

  # Sub-sections 5 and 6 (LMM) / 4 and 5 (GLMM): Pairwise Contrasts + CLD (separated)
  if (!is.null(res$posthoc) && length(res$posthoc) > 0) {
    # Omnibus note
    if (!is.null(res$anova)) {
      .trt_p <- tryCatch({
        av_row <- res$anova[grepl("^treatment$", res$anova$term, perl=TRUE, ignore.case=TRUE), ]
        if (nrow(av_row) > 0) as.numeric(av_row$p_value[1]) else NA_real_
      }, error = function(e) NA_real_)
    } else { .trt_p <- NA_real_ }

    # Sub-section 5 (LMM) / 4 (GLMM): Pairwise Contrasts
    doc <- .add_h3(doc, paste0(section_num, ".", .sub + 4L,
                               "  Pairwise Contrasts (Tukey)"))
    if (!is.na(.trt_p) && .trt_p >= 0.05)
      doc <- .add_par(doc, paste0("(Omnibus treatment ", fmt_p(.trt_p),
                                  "; contrasts shown for description only.)"))

    .any_contrasts <- FALSE
    for (ph_nm in names(res$posthoc)) {
      ph_entry <- res$posthoc[[ph_nm]]
      ctr_d    <- ph_entry$contrasts
      if (!is.null(ctr_d) && is.data.frame(ctr_d) && nrow(ctr_d) > 1) {
        .any_contrasts <- TRUE
        doc <- .add_par(doc, paste0("Term: ", ph_nm))
        ctr_show <- ctr_d[, intersect(
          c("contrast","estimate","SE","df","t.ratio","z.ratio","p.value"),
          names(ctr_d)), drop = FALSE]
        if ("p.value" %in% names(ctr_show))
          ctr_show$p.value <- vapply(ctr_show$p.value,
            function(p) if (is.na(p)) NA_character_ else fmt_p(p), character(1))
        ctr_show <- dplyr::mutate(ctr_show,
          dplyr::across(where(is.numeric), ~as.numeric(signif(.x, 3))))
        names(ctr_show)[names(ctr_show) == "t.ratio"] <- "t"
        names(ctr_show)[names(ctr_show) == "z.ratio"] <- "z"
        names(ctr_show)[names(ctr_show) == "p.value"] <- "p"
        .ctr_clean <- .drop_constant_cols(as.data.frame(ctr_show))
        for (.nt in .ctr_clean$notes) doc <- .add_par(doc, .nt)
        if (ncol(.ctr_clean$df) > 0)
          doc <- .add_ft(doc, .ft_styled(.ctr_clean$df))
      }
    }
    if (!.any_contrasts)
      doc <- .add_par(doc, "(Only two groups — no pairwise contrasts table; see CLD below.)")

    # Sub-section 5 (LMM) / 4+1 (GLMM): Simple-effects (within-stratum Tukey) ----------
    # For interaction terms: Tukey adjusted only within each level of the stratifier
    # (e.g. control vs exercise at each interval separately) — far fewer contrasts
    # and less conservative than the joint 28-contrast Tukey above.
    doc <- .add_h3(doc, paste0(section_num, ".", .sub + 5L,
                               "  Simple-Effects Contrasts (within-stratum Tukey)"))
    doc <- .add_par(doc, paste0(
      "Within-stratum Tukey: adjustment applied within each stratum level separately. ",
      "Typically 1–4 contrasts per stratum vs 28 in the joint table — ",
      "less conservative and matches patterns visible in graphs."))
    .any_simple <- FALSE
    for (ph_nm in names(res$posthoc)) {
      ph_entry <- res$posthoc[[ph_nm]]
      if (!is.null(ph_entry$contrasts_simple)) {
        for (.slot in names(ph_entry$contrasts_simple)) {
          cs_d <- ph_entry$contrasts_simple[[.slot]]$contrasts
          if (is.null(cs_d) || !is.data.frame(cs_d) || nrow(cs_d) == 0) next
          .any_simple <- TRUE
          .view_lbl <- gsub("_within_", " within ", .slot)
          doc <- .add_par(doc, paste0("Term: ", ph_nm, " | ", .view_lbl))
          tryCatch({
            cs_show <- cs_d[, intersect(
              c("contrast","treatment","timepoint","timepoint_f","zone",
                "estimate","SE","df","t.ratio","p.value"),
              names(cs_d)), drop = FALSE]
            if ("p.value" %in% names(cs_show))
              cs_show$p.value <- vapply(cs_show$p.value,
                function(p) if (is.na(p)) NA_character_ else fmt_p(p), character(1))
            cs_show <- dplyr::mutate(cs_show,
              dplyr::across(where(is.numeric), ~as.numeric(signif(.x, 3))))
            names(cs_show)[names(cs_show) == "t.ratio"] <- "t"
            names(cs_show)[names(cs_show) == "p.value"] <- "p"
            .cs_clean <- .drop_constant_cols(as.data.frame(cs_show))
            if (ncol(.cs_clean$df) > 0)
              doc <- .add_ft(doc, .ft_styled(.cs_clean$df))
          }, error = function(e)
            doc <<- .add_par(doc, paste0("(simple-effects table error: ", e$message, ")")))
        }
      }
    }
    if (!.any_simple)
      doc <- .add_par(doc, "(No interaction terms — simple-effects not applicable.)")

    # Sub-section 6 (LMM) / 6 (GLMM): CLD
    doc <- .add_h3(doc, paste0(section_num, ".", .sub + 6L,
                               "  Compact Letter Display (CLD)"))
    .any_cld <- FALSE
    for (ph_nm in names(res$posthoc)) {
      ph_entry <- res$posthoc[[ph_nm]]
      cld_d    <- ph_entry$cld
      if (!is.null(cld_d) && ".group" %in% names(cld_d)) {
        .any_cld <- TRUE
        doc <- .add_par(doc, paste0("Term: ", ph_nm))
        cld_show <- cld_d[, intersect(
          c("treatment","zone","timepoint_f","timepoint",
            "emmean","response","SE","lower.CL","upper.CL",".group"),
          names(cld_d))]
        cld_show <- dplyr::mutate(cld_show,
          dplyr::across(where(is.numeric), ~as.numeric(signif(.x, 3))))
        .cld_clean <- .drop_constant_cols(.pretty_names(as.data.frame(cld_show)))
        for (.nt in .cld_clean$notes) doc <- .add_par(doc, .nt)
        if (ncol(.cld_clean$df) > 0)
          doc <- .add_ft(doc, .ft_styled(.cld_clean$df))
      }
    }
    if (!.any_cld)
      doc <- .add_par(doc, "(CLD not available.)")

  } else {
    doc <- .add_h3(doc, paste0(section_num, ".", .sub + 4L,
                               "  Post-hoc Comparisons (Tukey)"))
    doc <- .add_par(doc,
      "(Post-hoc not available — model may not have converged or no valid terms.)")
  }

  doc
}

# ---- Safe per-section wrapper -----------------------------------------------
.doc_section <- function(doc, label, expr_fn) {
  tryCatch(expr_fn(doc),
    error = function(e) {
      cat("!! Word report section FAILED: ", label, "\n   ", conditionMessage(e), "\n", sep = "")
      tryCatch(officer::body_add_par(doc,
        paste0("[Section '", label, "' failed: ", conditionMessage(e), "]"),
        style = "Normal"), error = function(e2) doc)
    })
}

# ---- Build document ---------------------------------------------------------
doc <- tryCatch(officer::read_docx(), error = function(e) {
  warning("officer::read_docx() failed: ", conditionMessage(e))
  NULL
})
if (!is.null(doc)) doc <- .strip_heading_autonumber(doc)
if (is.null(doc)) {
  warning("Word report skipped — could not create document object.")
} else {

  # Title
  doc <- tryCatch({
    doc <- officer::body_add_par(doc, "Choice Experiment: Statistical Analysis Report",
                                  style = "heading 1")
    doc <- .add_par(doc, paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M")))
    doc <- .add_par(doc, "")
    doc
  }, error = function(e) {
    cat("!! Word report title section failed: ", conditionMessage(e), "\n")
    doc
  })

  # =========================================================
  # 1. Study Design and Data Summary
  # =========================================================
  doc <- .add_h1(doc, "1.  Study Design and Data Summary")
  doc <- .add_par(doc, paste0(
    "This experiment tested the effect of exercise choice on zebrafish activity and ",
    "zone preference in a split-flow arena. Two treatment groups were compared: ",
    "control and exercise choice. Fish were always housed five per arena; the fish ",
    "density variable records housing-tank density before the trial (four levels: ",
    paste(DENSITY_LEVELS_g, collapse = ", "), " fish per tank) and is treated as a ",
    "random effect, not a fixed experimental manipulation."
  ))

  .n_trial <- dplyr::n_distinct(df$trial_id)
  .n_tank  <- dplyr::n_distinct(df$tank)
  .n_tp    <- dplyr::n_distinct(df$timepoint[!is.na(df$timepoint)])
  .n_rows  <- nrow(df)
  .trt_tab <- table(df$treatment)
  .den_levels <- if ("fish_density" %in% names(df))
    paste(sort(unique(df$fish_density)), collapse = " / ") else "n/a"

  design_tbl <- data.frame(
    Item  = c("Total trial \u00d7 timepoint rows", "Unique trials", "Unique tanks",
              "Timepoints per trial", "Treatment levels", "Housing density levels"),
    Value = c(.n_rows, .n_trial, .n_tank, .n_tp,
              paste(names(.trt_tab), collapse = " / "),
              .den_levels),
    stringsAsFactors = FALSE
  )
  doc <- .add_ft(doc, .ft_styled(design_tbl))
  doc <- .add_par(doc, "")

  doc <- .add_h2(doc, "1.1  Random Effects Reference")
  doc <- .add_par(doc, "Each RE candidate reflects a grouping structure in the experimental design. Candidates compete by AICc; invalid or near-singular candidates are dropped.")
  re_ref_tbl <- data.frame(
    `RE formula` = c(
      "(1 | trial_id)",
      "(1 | tank)",
      "(1 | trial_id) + (1 | tank)",
      "(1 | tank / trial_id)",
      "(1 | trial_date)",
      "(1 | trial_id) + (1 | trial_date)",
      "(1 | fish_density_f)",
      "intercept-only (no RE)"
    ),
    Explanation = c(
      "Trial-level random intercept. Absorbs session-to-session variability; accounts for repeated timepoints from the same trial.",
      "Tank-level random intercept. Absorbs systematic differences between the four housing tanks (dropped when fewer than 5 tank levels).",
      "Independent trial and tank random intercepts. Partitions within-trial and between-tank variance simultaneously.",
      "Trials nested within tanks. Separates tank-level variance from within-tank trial-level variance.",
      "Date random intercept. Absorbs day-to-day environmental or procedural variation across recording sessions.",
      "Trial and date random intercepts. Used when both session identity and recording date contribute unexplained variance.",
      "Housing-density random intercept. Absorbs pre-trial group-size effects (4, 8, 12, or 16 fish per housing tank).",
      "No random grouping applied. Used when all RE candidates yield too few levels (<5) or singular fits with no valid simplification."
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  doc <- .add_ft(doc, .ft_styled(re_ref_tbl))
  doc <- .add_par(doc, "")

  doc <- .add_h2(doc, "1.2  Table Column Glossary")
  doc <- .add_par(doc, "All column headings used in tables throughout this report.")
  col_gloss_tbl <- data.frame(
    Name = c(
      "Response","Transform","SW p (raw)","SW p (trans)","Levene p (raw)","Levene p (trans)","Flag",
      "RE Formula","Fixed Formula","Fit","AICc","ΔAICc","Selected",
      "R²m","R²c",
      "Term","df","χ² / F","p","Sig.",
      "Contrast","Estimate","SE","t / z",
      "Treatment","Zone","Timepoint","EMmean","Response","Lower 95% CI","Upper 95% CI","CLD",
      "Analysis","Family","Test","p (raw)","Sig. (raw)","R²m","RE fallback",
      "Analysis","Trial p","Tank p","Trial β","Tank β","β ratio","Trial sig","Tank sig","Direction agree",
      "N","Satterthwaite df"
    ),
    Section = c(
      rep("Normality & Variance", 7),
      rep("Random Effects Selection", 6),
      rep("Model R²", 2),
      rep("ANOVA", 5),
      rep("Pairwise Contrasts", 4),
      rep("Post-hoc CLD", 8),
      rep("Cross-Analysis Summary", 7),
      rep("Pseudoreplication", 9),
      rep("General", 2)
    ),
    Explanation = c(
      # Normality
      "Outcome variable being analysed.",
      "Transformation applied to meet normality (none / log1p / sqrt).",
      "Shapiro-Wilk p-value on raw residuals.",
      "Shapiro-Wilk p-value on transformed residuals.",
      "Levene test p-value on raw data (variance homogeneity).",
      "Levene test p-value on transformed data.",
      "TRUE = normality or homoscedasticity assumption not met; interpret with caution.",
      # RE selection
      "Random effects formula evaluated as a candidate model.",
      "Fixed-effects structure label (categorical-tp = timepoint as factor; continuous-tp = timepoint as integer).",
      "Model fitting method used for AICc comparison (REML for LMM).",
      "Akaike Information Criterion corrected for small samples; lower = better fit.",
      "AICc difference from best-fit model; ΔAICc > 2 indicates substantially worse fit.",
      "TRUE for the model selected as best by AICc.",
      # R²
      "Marginal R²: proportion of variance explained by fixed effects only (Nakagawa & Schielzeth method via 'performance' package).",
      "Conditional R²: proportion of variance explained by fixed + random effects combined. If R²m ≈ R²c, random effects contribute negligible additional variance.",
      # ANOVA
      "Fixed effect or interaction term being tested.",
      paste0("Degrees of freedom. For LMM: Kenward-Roger (KR) df — primary inference. ",
             "Satterthwaite df also computed and shown in the comparison table for reference. ",
             "For GLMM: large-sample Wald df (typically 1 or #levels-1)."),
      "Test statistic: Kenward-Roger F (primary, LMM with lmerTest) or Wald χ² (GLMM / fallback). Satterthwaite F shown as comparison only.",
      "p-value for the term; exact to 3 sig. figs. when > 0.0001, otherwise p < 0.0001.",
      "Significance code: *** p < 0.001, ** p < 0.01, * p < 0.05, (blank) p ≥ 0.05.",
      # Contrasts
      "Label describing the pairwise comparison (e.g. 'exercise choice − control').",
      "Difference in marginal means on the model scale (log-odds for beta/Poisson; raw for Gaussian).",
      "Standard error of the contrast estimate.",
      "t or z test statistic for the contrast.",
      # CLD
      "Treatment group (control / exercise choice).",
      "Arena zone (flow / calm; or high / medium / low / calm for sub-zones).",
      "Recording timepoint (1 = pre-training, 2 = mid-training, 3 = post-training).",
      "Estimated marginal mean on the model (link) scale.",
      "Estimated marginal mean back-transformed to the response scale (shown when type='response').",
      "Lower bound of 95% confidence interval.",
      "Upper bound of 95% confidence interval.",
      "Compact letter display. Groups sharing a letter are not significantly different (Tukey p > 0.05).",
      # Cross-Analysis Summary
      "Analysis identifier (e.g. A1 Main-zone occupancy x Treatment x Tp).",
      "Model distribution family used: gaussian / beta / poisson / negbin.",
      "Inference method: F-KR = Kenward-Roger F-test; Wald-chisq = Wald chi-squared. Satterthwaite shown as comparison only in individual analysis sections.",
      "Treatment p-value; raw, no adjustment for multiple comparisons is applied anywhere in this study.",
      "Significant at raw p < 0.05.",
      "Marginal R² (fixed effects only).",
      "TRUE = optimal RE structure was invalid; intercept-only model fitted instead.",
      # Pseudoreplication
      "Analysis compared between trial-level and tank-level models.",
      "Treatment p-value from trial-level mixed model (primary analysis).",
      "Treatment p-value from tank-level ordinary LM on tank means.",
      "Treatment coefficient (beta) from trial-level model.",
      "Treatment coefficient (beta) from tank-level model.",
      "Ratio |Trial β / Tank β|; values near 1 indicate consistent effect magnitude across levels.",
      "TRUE = treatment significant at p < 0.05 in trial-level model.",
      "TRUE = treatment significant at p < 0.05 in tank-level model.",
      "TRUE = both models agree on the direction of the treatment effect.",
      # General
      "Sample size (number of datapoints per group). For aggregated analyses: number of Trials. For timepoint analyses: number of Trial × Timepoint combinations.",
      "Satterthwaite denominator df are computed alongside KR for comparison only (see ANOVA subsections). Primary inference uses Kenward-Roger. Wald df for GLMMs."
    ),
    stringsAsFactors = FALSE
  )
  doc <- .add_ft(doc, .ft_styled(col_gloss_tbl))
  doc <- .add_par(doc, "")

  # =========================================================
  # 2. Statistical Decision Log
  # =========================================================
  doc <- .add_h1(doc, "2.  Statistical Decision Log")

  doc <- .add_h2(doc, "2.1  Observation unit")
  doc <- .add_par(doc, paste0(
    "The biological unit of analysis is the trial \u00d7 timepoint (school). ",
    "Persistent within-trial fish identity is unreliable in the choice arena, so ",
    "all indicators are identity-free school-level summaries. Tank-level random ",
    "effects absorb repeated-measures dependence between timepoints of the same tank."
  ))

  doc <- .add_h2(doc, "2.2  Housing density as random effect")
  doc <- .add_par(doc, paste0(
    "Fish density records the number of fish in the housing tank prior to the trial. ",
    "Because all trials use five fish in the arena, density does not vary within the ",
    "experimental unit and is therefore not a fixed effect. It is included as a ",
    "random-effects candidate (1 | fish_density_f) so that variance attributable to ",
    "pre-trial housing conditions can be partitioned away from treatment effects."
  ))

  doc <- .add_h2(doc, "2.3  Random effects candidate set")
  doc <- .add_par(doc, paste0(
    "Long-format zone occupancy analyses (A1, A2) use RE candidates: ",
    paste(names(RE_LONG), collapse = ", "), ". ",
    "Wide-format trial-level analyses (A3-A16) use: ",
    paste(names(RE_WIDE), collapse = ", "), ". ",
    "(1 | tank / trial_id) absorbs both tank-level and trial-level dependence ",
    "between zone rows. RE candidates are compared by AICc on REML fits. ",
    "(1 | tank) is dropped automatically when fewer than five tank levels are present; ",
    "(1 | trial_id) is dropped for aggregated analyses with one row per trial."
  ))

  doc <- .add_h2(doc, "2.4  Model family and transformation protocol")
  doc <- .add_par(doc, paste0(
    "Proportion outcomes (zone occupancy A1-A2; proportion active A5-A6) are modelled ",
    "with beta GLMM via glmmTMB (family = beta_family(link = 'logit')); response ",
    "values are squeezed to the open interval (0, 1) using (y*(n-1)+0.5)/n. ",
    "Count outcomes (zone switches A3-A4) are modelled with Poisson or negative-",
    "binomial GLMM with an offset for observation time. ",
    "All remaining continuous outcomes (NND, IID, school area, school ",
    "speed, zone flux) use Gaussian LMM with optional log1p or sqrt transform selected ",
    "by Shapiro-Wilk p on fixed-effects-only LM residuals. ",
    "When glmmTMB is unavailable a Gaussian fallback with logit transform is used for ",
    "proportions; the 'Family' column in Section 6 records the family actually used."
  ))

  doc <- .add_h2(doc, "2.5  Inference method")
  doc <- .add_par(doc, paste0(
    "Random-effects structure is selected by AICc on REML fits. ",
    "The winning structure is then refit with REML = TRUE (LMM) or full likelihood ",
    "(GLMM) before inference. ",
    "For LMM, Type III ANOVA uses Kenward-Roger (KR) F-tests (lmerTest/pbkrtest) as the ",
    "primary inference method; otherwise Wald chi-squared is used. ",
    "Satterthwaite df are also computed and displayed alongside KR in each analysis section ",
    "for comparison; they are NOT used for primary inference. ",
    "The 'Test' column in Section 6 records the primary method: F-KR = Kenward-Roger F; ",
    "Wald-chisq = Wald chi-squared. ",
    "Post-hoc pairwise comparisons use emmeans with Tukey adjustment."
  ))

  doc <- .add_h2(doc, "2.6  Multiple testing correction")
  doc <- .add_par(doc, paste0(
    "No adjustment for multiple comparisons is applied anywhere in this study ",
    "(author decision, 2026-08-07; see Supplementary Table S7/2.7). The pre-specified ",
    "analysis groupings are retained below for organisational reference only: ",
    "(1) non-zone core indicators (5 tests, aggregated treatment main effects): ",
    "logit(prop_flow) [A1c, compositionally correct], switches/session [A4], NND [A8], ",
    "school speed [A16]. ",
    "NOTE: A1_agg (beta GLMM on long-format zones) is excluded from this group ",
    "because it is compositionally misspecified (prop_flow + prop_calm = 1 within each ",
    "session; the long-format stacking doubles observations and introduces algebraic ",
    "dependency). It is replaced by A1c. ",
    "(2) Treatment:Timepoint interactions (5 tests): ",
    "logit(prop_flow):Timepoint [A1], switches:Timepoint [A3], NND:Timepoint [A7], ",
    "school speed:Timepoint [A15]. ",
    "All remaining analyses are reported descriptively: ",
    "sub-zone log-ratios (A2a-c agg + TP), IID, school area, flux, prop_active, and the ",
    "legacy A1_agg beta GLMM (retained for graphs). ",
    "All p-values throughout are raw."
  ))

  doc <- .add_h2(doc, "2.7  Family-specific modelling notes")
  doc <- .add_par(doc, paste0(
    "Zone switches (A3/A4): switches_per_session modelled as Gaussian LMM with F-KR inference. ",
    "Proportion active (A5/A6): beta GLMM. ",
    "Zone occupancy (A1/A1c): logit(prop_flow) Gaussian LMM (compositionally correct scalar). ",
    "Zone occupancy (A2a-c): log-ratio LMMs using calm as reference zone. ",
    "Zone flux (A_flux): zone_flux_per_session is a continuous rate; Gaussian LMM ",
    "with optional transformation unless a raw integer count column is available."
  ))

  doc <- .add_h2(doc, "2.8  Continuous timepoint model variant")
  doc <- .add_par(doc, paste0(
    "For analyses involving timepoint (A2, A4, A5), an additional set of candidates ",
    "is fit treating timepoint as a continuous integer predictor rather than a ",
    "factor. Both categorical and continuous timepoint candidates compete in the ",
    "same AICc selection; the best across both sets is reported. This allows the ",
    "data to indicate whether the temporal effect is better described as a linear ",
    "trend or as segment-specific means."
  ))

  doc <- .add_h2(doc, "2.9  Contrast coding")
  doc <- .add_par(doc, paste0(
    "All unordered factors use sum-to-zero contrasts (contr.sum) scoped within each ",
    "model-fitting call. This ensures that Type III main effects are marginal (averaged ",
    "over other factor levels) rather than effects at the reference level, which is ",
    "required for correct interpretation of main effects in interaction models."
  ))

  doc <- .add_h2(doc, "2.10  Aggregated vs timepoint analyses — rationale and unit")
  doc <- .add_par(doc, paste0(
    "Physical trial: one physical trial = one unique (tank × date × treatment) combination. ",
    "With 4 tanks, 2 recording dates per tank, and 2 treatments per date, there are ",
    "N = 16 physical trials. Each physical trial comprises 3 successive 20-minute ",
    "recording sessions (timepoints 1–3), giving 48 total session rows in the data. "
  ))
  doc <- .add_par(doc, paste0(
    "Two-level analysis approach: ",
    "(i) Aggregated (primary, N = 16): metrics averaged across the 3 timepoints per ",
    "physical trial before modelling. This eliminates the repeated-measures dependence ",
    "between timepoints and gives the cleanest, most conservative test of the overall ",
    "treatment effect. The RE set at the aggregated level is limited to (1|tank), ",
    "(1|trial_date), and (1|fish_density_f) because phys_trial_id is the observational ",
    "unit (one row per trial) and cannot also be a grouping factor. ",
    "(ii) Timepoint (supplementary, N = 48): all three sessions modelled jointly with ",
    "treatment × timepoint interaction and (1|phys_trial_id) as the repeated-measures RE. ",
    "This tests whether the treatment effect changes over the recording period. "
  ))
  doc <- .add_par(doc, paste0(
    "Physical trial identifier: phys_trial_id is derived as as.integer(factor(paste(tank, ",
    "trial_date, treatment))) and is constructed identically in df, df_long, df_gd, and all ",
    "aggregated datasets, ensuring consistent clustering across analyses. ",
    "Date nesting: each trial_date appears under exactly one tank ",
    "(date is fully nested within tank), so (1|tank) and (1|trial_date) are correlated. ",
    "When both appear as RE candidates, AICc selects the better-fitting one; the winner ",
    "is reported in the aicc_selection.csv for each analysis."
  ))
  doc <- .add_par(doc, paste0(
    "Compositional correction (A1 / A1c / A2): zone occupancy proportions sum to 1 ",
    "within each session by construction. Fitting a Gaussian LMM or beta GLMM on ",
    "stacked long-format zone rows treats two (or four) algebraically dependent observations ",
    "as independent, inflating effective sample size and producing biased standard errors. ",
    "Two separate corrections were implemented: ",
    "(1) Main zones (flow/calm): the logit(prop_flow) scalar is the natural 1-df summary. ",
    "A1c (logit(prop_flow) ~ treatment, N = 16) replaces A1_agg as the compositionally correct ",
    "aggregated main-zone test (section 6.3). A1_agg (beta GLMM on long-format) is ",
    "retained solely for graph annotation; its p-value is NOT used for inference. ",
    "The TP version (A1) is refactored to logit(prop_flow) ~ treatment * timepoint_f + ",
    "(1|phys_trial_id) (N = 48) and is reported in section 6.4. ",
    "(2) Sub-zones (high/medium/low/calm): the 4-zone Gaussian LMM is rank-deficient ",
    "(sum-to-1 constraint makes the design matrix singular), producing numerically ",
    "degenerate p ≈ 1 results in earlier runs. Replacement: 3 separate LMMs on log-ratio ",
    "responses — log(prop_high/prop_calm), log(prop_medium/prop_calm), ",
    "log(prop_low/prop_calm) — with calm as the reference simplex component. A Haldane ",
    "offset of 0.0001 is added to all proportions before log-ratio computation to handle ",
    "sessions where a zone is never entered. All 6 sub-zone log-ratio models (3 agg + 3 TP) ",
    "are classified as exploratory (section 6.5); raw p, no adjustment for multiple comparisons."
  ))

  # =========================================================
  # =========================================================
  # 4. Part 1: Aggregated Results (Treatment Main Effect)
  # =========================================================
  doc <- .add_h1(doc, "4.  Part 1: Aggregated Analyses (Treatment Main Effect)")
  doc <- .add_par(doc, paste0(
    "Analyses A1c (logit(prop_flow), compositionally correct), A2a-c (sub-zone log-ratios vs calm), ",
    "A2_pairs (sub-zone pairwise comparisons, exploratory), A2_cells (cell-mean ",
    "cross-comparisons across treatment x zone, exploratory), A4, A6, ",
    "A_flux_agg, A8, A10, A12, A14, A16 aggregate across timepoints and test the ",
    "overall treatment main effect. A1_agg (beta GLMM) is retained for graph annotation ",
    "only; it is compositionally misspecified (see 2.10) and its p-value is not used ",
    "for inference. All p-values throughout are raw. ",
    "Post-hoc CLD tables are included for all analyses regardless of p-value."
  ))

  # A1c is the compositionally correct main-zone model; A1_agg retained for graph
  # annotation only (compositionally misspecified -- see 2.10)
  doc <- .analysis_block(doc, "4.1c", "zone_flow_logit_aggregated",
    "A1c: logit(prop_flow) x Treatment (compositional, aggregated)",
    "zone_flow_logit_aggregated", res_zone_flow_logit_agg, FALSE)
  doc <- .analysis_block(doc, "4.1", "zone_main_aggregated",
    "A1_agg: Main-zone Occupancy beta GLMM (aggregated, compositionally misspecified, retained for graphs only)",
    "zone_main_aggregated", res_zone_main_agg, FALSE)
  # A2_agg: three log-ratio LMMs replacing the compositionally flawed 4-zone Gaussian LMM
  doc <- .analysis_block(doc, "4.2a", "zone_sec_high_aggregated",
    "A2a: log(high/calm) x Treatment (aggregated)",
    "zone_sec_high_aggregated", res_zone_sec_high_agg, FALSE)
  doc <- .analysis_block(doc, "4.2b", "zone_sec_medium_aggregated",
    "A2b: log(medium/calm) x Treatment (aggregated)",
    "zone_sec_medium_aggregated", res_zone_sec_med_agg, FALSE)
  doc <- .analysis_block(doc, "4.2c", "zone_sec_low_aggregated",
    "A2c: log(low/calm) x Treatment (aggregated)",
    "zone_sec_low_aggregated", res_zone_sec_low_agg, FALSE)
  # A2_pairs_agg: stacked alr LMM — sub-zone pairwise contrasts (high/medium/low vs each other)
  doc <- .add_h2(doc, "4.2d  A2_pairs [EXPLORATORY]: Sub-zone alr Pairwise Comparisons (aggregated)")
  doc <- .add_par(doc, paste0(
    "Stacked LMM fitting all three sub-zone log-ratios vs calm simultaneously ",
    "(lr ~ treatment * zone_lr, where zone_lr ∈ {high, medium, low}). ",
    "The zone_lr Tukey contrasts yield pairwise sub-zone comparisons (high vs medium, ",
    "high vs low, medium vs low) on the alr scale. Exploratory only — no BH correction. ",
    "The treatment * zone_lr interaction term tests whether the treatment effect differs ",
    "across sub-zones."
  ))
  doc <- .analysis_block(doc, "4.2d", "zone_sec_lr_pairs_aggregated",
    "A2_pairs: Stacked sub-zone alr LMM — pairwise sub-zone contrasts (aggregated)",
    "zone_sec_lr_pairs_aggregated", res_zone_sec_lr_pairs_agg, FALSE)
  # A2_cells main: cell-mean beta-GLMM on (flow + calm) x treatment grid
  doc <- .add_h2(doc, "4.2e  A2_cells [EXPLORATORY]: Main-zone cell-mean beta-GLMM (aggregated)")
  doc <- .add_par(doc, paste0(
    "Beta GLMM on raw proportions (p ~ treatment * zone + (1|trial)) over the ",
    "long-stacked main-zone composition (flow + calm). Yields cell means on the ",
    "proportion scale and Tukey-adjusted pairwise contrasts across the 4 ",
    "(treatment x zone) cells, including diagonal contrasts such as ",
    "calm x Control vs flow x Exercise Choice. Exploratory only — no BH ",
    "correction; the per-trial logit_flow LMM (4.1c) carries the headline test."
  ))
  doc <- .analysis_block(doc, "4.2e", "zone_main_cells_aggregated",
    "A2_cells: Main-zone cell-mean beta-GLMM (aggregated)",
    "zone_main_cells_aggregated", res_zone_main_cells_agg, FALSE)
  # A2_cells sub: cell-mean CLR LMM on the 4-part sub-zone composition
  doc <- .add_h2(doc, "4.2f  A2_cells [EXPLORATORY]: Sub-zone cell-mean CLR LMM (aggregated)")
  doc <- .add_par(doc, paste0(
    "Gaussian LMM on the centred log-ratio (CLR) of the 4-part area-normalised ",
    "sub-zone composition (high, medium, low, calm). CLR = log(p_an,zone) - mean(log(p_an,.)) ",
    "computed within each trial's composition. clr ~ treatment * zone + ",
    "(1|trial) + (1|zone:trial); the trial RE absorbs the sum-to-zero ",
    "constraint. Yields all 8 (treatment x zone) cell means including calm, ",
    "and Tukey-adjusted pairwise contrasts (28 pairs). Exploratory only."
  ))
  doc <- .analysis_block(doc, "4.2f", "zone_sec_cells_aggregated",
    "A2_cells: Sub-zone cell-mean CLR LMM (aggregated)",
    "zone_sec_cells_aggregated", res_zone_sec_cells_agg, FALSE)
  # Jacobs' D joint beta-GLMMs: treatment × zone (aggregated)
  doc <- .add_h2(doc, "4.2g  Jacobs' D: Main-zone Treatment × Zone (joint beta-GLMM, aggregated)")
  doc <- .add_par(doc, paste0(
    "Joint beta-GLMM on Jacobs' D for main zones (flow, calm): ",
    "treatment * zone as crossed fixed effects on stacked long data ",
    "(one row per trial × zone, area-availability adjusted). ",
    "Yields a treatment:zone interaction test and joint Tukey CLDs across all ",
    "4 (treatment × zone) cells. These CLDs are the source for Figure S7 panel A. ",
    "Exploratory in the same sense as 4.2e-4.2f; the per-zone Jacobs' D beta-GLMMs ",
    "(section 6.1) report the treatment main effect per zone, while this joint ",
    "model tests the treatment x zone interaction directly."
  ))
  doc <- .analysis_block(doc, "4.2g", "jacobs_beta_main_joint_agg",
    "Jacobs' D main zones: Treatment × Zone (joint beta-GLMM, aggregated)",
    "jacobs_beta_main_joint_agg",
    if (exists("res_jac_beta_main_joint_agg")) res_jac_beta_main_joint_agg else NULL,
    is_glmm = TRUE)
  doc <- .add_h2(doc, "4.2h  Jacobs' D: Sub-zone Treatment × Zone (joint beta-GLMM, aggregated)")
  doc <- .add_par(doc, paste0(
    "Joint beta-GLMM on Jacobs' D for sub-zones (high, medium, low, calm): ",
    "treatment * zone as crossed fixed effects on stacked long data ",
    "(one row per trial × zone). ",
    "Yields a treatment:zone interaction test and joint Tukey CLDs across all ",
    "8 (treatment × zone) cells. These CLDs are the source for Figure S7 panel B."
  ))
  doc <- .analysis_block(doc, "4.2h", "jacobs_beta_sec_joint_agg",
    "Jacobs' D sub-zones: Treatment × Zone (joint beta-GLMM, aggregated)",
    "jacobs_beta_sec_joint_agg",
    if (exists("res_jac_beta_sec_joint_agg")) res_jac_beta_sec_joint_agg else NULL,
    is_glmm = TRUE)
  doc <- .analysis_block(doc, "4.3", "switches_aggregated",
    "A4: Zone Switches x Treatment (aggregated)",
    "switches_aggregated", res_switch_agg, FALSE)
  doc <- .analysis_block(doc, "4.4", "active_aggregated",
    "A6: Proportion Active x Treatment (aggregated)",
    "active_aggregated", res_active_agg, FALSE)
  doc <- .analysis_block(doc, "4.5", "flux_aggregated",
    "A_flux: Zone Flux x Treatment (aggregated)",
    "flux_aggregated", res_flux_agg, FALSE)
  if (!is.null(df_gd)) {
    doc <- .analysis_block(doc, "4.6", "nnd_aggregated",
      "A8: Mean NND x Treatment (aggregated)",
      "nnd_aggregated", res_nnd_agg, FALSE)
    doc <- .analysis_block(doc, "4.8", "iid_aggregated",
      "A12: Mean IID x Treatment (aggregated)",
      "iid_aggregated", res_iid_agg, FALSE)
    doc <- .analysis_block(doc, "4.9", "hull_area_aggregated",
      "A14: School Area x Treatment (aggregated)",
      "hull_area_aggregated", res_hull_agg, FALSE)
    doc <- .analysis_block(doc, "4.10", "centroid_speed_aggregated",
      "A16: School Speed x Treatment (aggregated)",
      "centroid_speed_aggregated", res_cspd_agg, FALSE)
  }

  # =========================================================
  # 5. Part 2: Timepoint Results
  # =========================================================
  doc <- .add_h1(doc, "5.  Part 2: Timepoint Analyses (Treatment x Timepoint Interaction)")
  doc <- .add_par(doc, paste0(
    "Analyses A1 (logit(prop_flow) × timepoint, compositional LMM), A2a-c (sub-zone ",
    "log-ratios × timepoint), A2_pairs (sub-zone pairwise comparisons × timepoint, exploratory), ",
    "A2_cells (cell-mean cross-comparisons × timepoint, exploratory), ",
    "A3, A5, A_flux_tp, A7, A9, A11, A13, A15 model the ",
    "Treatment × Timepoint interaction. Post-hoc CLD tables are included for all ",
    "analyses regardless of p-value."
  ))

  # A1 TP: now logit(prop_flow) ~ treatment * timepoint_f LMM (compositionally correct)
  doc <- .analysis_block(doc, "5.1", "zone_main_timepoint",
    "A1: logit(prop_flow) x Treatment x Timepoint (compositional LMM)",
    "zone_main_timepoint", res_zone_main, FALSE)
  # A2 TP: three log-ratio LMMs
  doc <- .analysis_block(doc, "5.2a", "zone_sec_high_timepoint",
    "A2a: log(high/calm) x Treatment x Timepoint",
    "zone_sec_high_timepoint", res_zone_sec_high_tp, FALSE)
  doc <- .analysis_block(doc, "5.2b", "zone_sec_medium_timepoint",
    "A2b: log(medium/calm) x Treatment x Timepoint",
    "zone_sec_medium_timepoint", res_zone_sec_med_tp, FALSE)
  doc <- .analysis_block(doc, "5.2c", "zone_sec_low_timepoint",
    "A2c: log(low/calm) x Treatment x Timepoint",
    "zone_sec_low_timepoint", res_zone_sec_low_tp, FALSE)
  # A2_pairs_tp: stacked alr LMM — sub-zone pairwise contrasts x timepoint
  doc <- .add_h2(doc, "5.2d  A2_pairs [EXPLORATORY]: Sub-zone alr Pairwise Comparisons (timepoint-resolved)")
  doc <- .add_par(doc, paste0(
    "Same stacked alr LMM as 4.2d but with a three-way treatment × zone_lr × timepoint_f ",
    "fixed structure. The zone_lr Tukey contrasts yield interval-resolved pairwise sub-zone ",
    "comparisons. Exploratory only — no BH correction."
  ))
  doc <- .analysis_block(doc, "5.2d", "zone_sec_lr_pairs_timepoint",
    "A2_pairs: Stacked sub-zone alr LMM — pairwise sub-zone contrasts (timepoint-resolved)",
    "zone_sec_lr_pairs_timepoint", res_zone_sec_lr_pairs_tp, FALSE)
  # A2_cells main TP: cell-mean beta-GLMM x timepoint
  doc <- .add_h2(doc, "5.2e  A2_cells [EXPLORATORY]: Main-zone cell-mean beta-GLMM (timepoint-resolved)")
  doc <- .add_par(doc, paste0(
    "Same beta-GLMM as 4.2e but with a three-way treatment x zone x ",
    "timepoint_f fixed structure. Cell-grid contrasts can be requested at ",
    "specific intervals via emmeans subsetting. Exploratory only."
  ))
  doc <- .analysis_block(doc, "5.2e", "zone_main_cells_timepoint",
    "A2_cells: Main-zone cell-mean beta-GLMM (timepoint-resolved)",
    "zone_main_cells_timepoint", res_zone_main_cells_tp, FALSE)
  # A2_cells sub TP: cell-mean CLR LMM x timepoint
  doc <- .add_h2(doc, "5.2f  A2_cells [EXPLORATORY]: Sub-zone cell-mean CLR LMM (timepoint-resolved)")
  doc <- .add_par(doc, paste0(
    "Same CLR LMM as 4.2f but with a three-way treatment x zone x timepoint_f ",
    "fixed structure. Provides cell-by-cell Tukey contrasts within and across ",
    "intervals. Exploratory only."
  ))
  doc <- .analysis_block(doc, "5.2f", "zone_sec_cells_timepoint",
    "A2_cells: Sub-zone cell-mean CLR LMM (timepoint-resolved)",
    "zone_sec_cells_timepoint", res_zone_sec_cells_tp, FALSE)
  doc <- .analysis_block(doc, "5.3", "switches_timepoint",
    "A3: Zone Switches x Treatment x Timepoint",
    "switches_timepoint", res_switch_tp, FALSE)
  doc <- .analysis_block(doc, "5.4", "active_timepoint",
    "A5: Proportion Active x Treatment x Timepoint",
    "active_timepoint", res_active_tp, FALSE)
  doc <- .analysis_block(doc, "5.5", "flux_timepoint",
    "A_flux_tp: Zone Flux x Treatment x Timepoint",
    "flux_timepoint", res_flux_tp, FALSE)
  if (!is.null(df_gd)) {
    doc <- .analysis_block(doc, "5.6", "nnd_timepoint",
      "A7: Mean NND x Treatment x Timepoint",
      "nnd_timepoint", res_nnd_tp, FALSE)
    doc <- .analysis_block(doc, "5.8", "iid_timepoint",
      "A11: Mean IID x Treatment x Timepoint",
      "iid_timepoint", res_iid_tp, FALSE)
    doc <- .analysis_block(doc, "5.9", "hull_area_timepoint",
      "A13: School Area x Treatment x Timepoint",
      "hull_area_timepoint", res_hull_tp, FALSE)
    doc <- .analysis_block(doc, "5.10", "centroid_speed_timepoint",
      "A15: School Speed x Treatment x Timepoint",
      "centroid_speed_timepoint", res_cspd_tp, FALSE)
  }

  # =========================================================
  # 6. Cross-Analysis BH-Corrected Summary
  # =========================================================
  doc <- .add_h1(doc, "6.  Cross-Analysis Summary (raw p; no multiple-comparisons adjustment)")
  doc <- .add_par(doc, paste0(
    "No adjustment for multiple comparisons is applied anywhere in this study ",
    "(consistent with Methods 2.8 and Supplementary Table S7/S8); all p-values below ",
    "are raw. The groupings in 6.1-6.5 (zone preference, non-zone core indicators, ",
    "Treatment:Timepoint interactions, sensitivity check, descriptive) reflect the ",
    "pre-specified analysis families and are retained for organisational reference only. ",
    "Sig. (raw): p_raw < 0.05. ",
    "Family: model distribution used. Test: F-KR = Kenward-Roger F; Wald-chisq = Wald χ². ",
    "R²m = marginal R² (fixed effects only). RE fallback = TRUE when optimal RE was dropped. ",
    "The alr/logit zone-preference tests (A1c, A1, A2 sub-zone) are reported ",
    "in section 6.5 alongside the other descriptive/exploratory indicators."
  ))
  .bh_cols_avail <- function(df) {
    want <- c("Analysis","family","stat_type","p_raw","sig_raw",
              "R2m","re_fallback","normality_flag")
    out <- df[, intersect(want, names(df)), drop = FALSE]
    .fp <- function(p) vapply(p, function(x) if (is.na(x)) NA_character_ else fmt_p(x), character(1))
    if ("p_raw" %in% names(out)) out$p_raw <- .fp(out$p_raw)
    if ("R2m"   %in% names(out)) out$R2m   <- vapply(out$R2m, fmt_num, character(1))
    out
  }
  doc <- .add_h2(doc, "6.1  Zone preference: Jacobs' D per-zone beta-GLMM")
  doc <- .add_par(doc, paste0(
    "Per-zone Jacobs' D beta-GLMMs: 8 tests (4 zones × {aggregated, timepoint-resolved}); ",
    "raw p, no adjustment for multiple comparisons. Jacobs' D = (r − p)/(r + p − 2rp); ",
    "rescaled to (D+1)/2 for beta family, Smithson-Verkuilen squeeze."
  ))
  doc <- .add_ft(doc, .ft_styled(.pretty_names(.bh_cols_avail(bh_primary))))

  # 6.1b: Jacobs' D treatment × zone — placed immediately after 6.1
  if (exists("res_jac_beta_main_joint_agg") || exists("res_jac_beta_sec_joint_agg")) {
    doc <- .add_h2(doc, "6.1b  Jacobs' D: Treatment × Zone interaction (joint beta-GLMM)")
    doc <- .add_par(doc, paste0(
      "Joint beta-GLMMs testing the treatment × zone interaction on Jacobs' D: ",
      "a single model per zone set (main: flow + calm; sub: high + medium + low + calm) ",
      "with treatment * zone crossed fixed effects on stacked long data. ",
      "Yields the treatment:zone interaction term and joint Tukey CLDs across all ",
      "(treatment × zone) cells. CLD letters are used in Figure S7. ",
      "Kept separate from the per-zone table in 6.1; raw p, no adjustment for multiple comparisons."
    ))
    bh_jac_joint <- dplyr::bind_rows(
      if (exists("res_jac_beta_main_joint_agg") && !is.null(res_jac_beta_main_joint_agg))
        .bh_row("Jacobs' D main zones — Treatment × Zone (agg)",
                res_jac_beta_main_joint_agg, "treatment.*zone|zone.*treatment")
      else NULL,
      if (exists("res_jac_beta_sec_joint_agg") && !is.null(res_jac_beta_sec_joint_agg))
        .bh_row("Jacobs' D sub zones  — Treatment × Zone (agg)",
                res_jac_beta_sec_joint_agg,  "treatment.*zone|zone.*treatment")
      else NULL
    )
    if (!is.null(bh_jac_joint) && nrow(bh_jac_joint) > 0)
      doc <- .add_ft(doc, .ft_styled(.pretty_names(.bh_cols_avail(bh_jac_joint))))
  }

  if (exists("bh_jacobs_sens") && is.data.frame(bh_jacobs_sens) && nrow(bh_jacobs_sens) > 0) {
    doc <- .add_h2(doc, "6.2  Jacobs' D sensitivity (atanh-LMM, no BH)")
    doc <- .add_par(doc, paste0(
      "Robustness check: same 8 zone x stratum cells as 6.1 but fitted as ",
      "Gaussian LMMs on atanh(D * (1 - 1/(2N))) with F-KR inference. ",
      "Reported alongside the beta-GLMM results in 6.1; raw p only. ",
      "Disagreement in significance between 6.1 and 6.2 typically reflects ",
      "boundary-saturation effects (compression by atanh vs squeeze under beta)."
    ))
    doc <- .add_ft(doc, .ft_styled(.pretty_names(.bh_cols_avail(bh_jacobs_sens))))
  }

  if (exists("bh_suppl_primary") && is.data.frame(bh_suppl_primary) && nrow(bh_suppl_primary) > 0) {
    doc <- .add_h2(doc, "6.3  Non-zone core indicators: Switches, NND, School speed")
    doc <- .add_par(doc, paste0(
      nrow(bh_suppl_primary),
      " non-zone core indicators (trial-aggregated, treatment main effect); ",
      "raw p, no adjustment for multiple comparisons."
    ))
    doc <- .add_ft(doc, .ft_styled(.pretty_names(.bh_cols_avail(bh_suppl_primary))))
  }

  if (exists("bh_suppl") && is.data.frame(bh_suppl) && nrow(bh_suppl) > 0) {
    doc <- .add_h2(doc, "6.4  Supplementary-TP: Treatment x Timepoint Interactions")
    doc <- .add_par(doc, paste0(
      nrow(bh_suppl),
      " Treatment:Timepoint interaction tests (one per non-zone core metric); ",
      "raw p, no adjustment for multiple comparisons."
    ))
    doc <- .add_ft(doc, .ft_styled(.pretty_names(.bh_cols_avail(bh_suppl))))
  }

  if (exists("bh_exploratory") && is.data.frame(bh_exploratory) && nrow(bh_exploratory) > 0) {
    doc <- .add_h2(doc, "6.5  Descriptive / Exploratory")
    doc <- .add_par(doc, paste0(
      "Raw p-values only. Includes the alr/logit zone-preference tests ",
      "(A1c, A1, A2 sub-zone) alongside other exploratory indicators. These ", nrow(bh_exploratory),
      " analyses are presented for completeness; interpret with caution."
    ))
    doc <- .add_ft(doc, .ft_styled(.pretty_names(.bh_cols_avail(bh_exploratory))))
  }

  # 6.6 Sub-zone alr pairwise contrasts (exploratory)
  doc <- .add_h2(doc, "6.6  Sub-zone alr Pairwise Contrasts (exploratory, no BH)")
  doc <- .add_par(doc, paste0(
    "Tukey-adjusted pairwise contrasts from the stacked sub-zone alr LMM (A2_pairs). ",
    "Each row compares two sub-zones (high, medium, or low) on the area-normalised log-ratio ",
    "scale. The estimate equals the log-ratio difference between the two sub-zones relative to ",
    "calm (i.e. log(r_zone1/r_zone2) plus a constant area correction that cancels). ",
    "Stratum 'agg' = trial-aggregated (N = 16); 'tp' = timepoint-resolved (N = 48). ",
    "Exploratory only — no BH correction."
  ))
  if (exists(".lr_pairs_summary") && is.data.frame(.lr_pairs_summary) &&
      nrow(.lr_pairs_summary) > 0) {
    .lps_show <- .lr_pairs_summary
    .lps_cols <- intersect(c("contrast","estimate","SE","df","t.ratio","p.value","stratum"),
                           names(.lps_show))
    .lps_show <- .lps_show[, .lps_cols, drop = FALSE]
    if ("p.value" %in% names(.lps_show))
      .lps_show$p.value <- vapply(.lps_show$p.value,
                                  function(p) if (is.na(p)) NA_character_ else fmt_p(p),
                                  character(1))
    doc <- .add_ft(doc, .ft_styled(.lps_show))
  } else {
    doc <- .add_par(doc, "(Sub-zone pairwise contrasts not available.)")
  }

  # 6.7 Hedges' g effect size summary (renumbered from 6.4)
  doc <- .add_h2(doc, "6.7  Effect Sizes (Hedges' g, Gaussian LMM only)")
  doc <- .add_par(doc, paste0(
    "Hedges' g = |ΔM| / σ_residual × J(df_resid), where J = 1 − 3/(4×df−1) ",
    "is the small-sample correction factor (Hedges 1981). ",
    "95% CI derived by propagating the emmeans contrast SE through sigma. ",
    "g ≈ 0.2 small, 0.5 medium, 0.8 large (Cohen 1988). ",
    "Reported only for Gaussian LMM (aggregated analyses); NA for beta/Poisson/NB families."
  ))
  if (exists("hedges_g_tbl") && is.data.frame(hedges_g_tbl)) {
    .hg_show <- hedges_g_tbl[, c("Analysis","g","ci_lo","ci_hi","method"), drop = FALSE]
    names(.hg_show) <- c("Analysis", "Hedges' g", "95% CI lower", "95% CI upper", "Method")
    doc <- .add_ft(doc, .ft_styled(.hg_show))
  } else {
    doc <- .add_par(doc, "(Hedges' g not computed.)")
  }

  # =========================================================
  # 7. Pseudoreplication Assessment
  # =========================================================
  doc <- .add_h1(doc, "7.  Pseudoreplication Assessment")
  doc <- .add_h2(doc, "7.1  Rationale")
  doc <- .add_par(doc, paste0(
    "The analysis unit is already trial \u00d7 timepoint, but multiple timepoints ",
    "from the same tank are still non-independent. To check that conclusions are ",
    "not driven by repeated-measures structure, A1 and A6 were re-run as simple ",
    "tank-level LMs on tank-aggregated means. Effect directions and beta coefficients ",
    "are compared between the trial-level model (Trial \u03b2) and the tank-level model ",
    "(Tank \u03b2). Direction_agree = TRUE indicates both models find the treatment effect ",
    "in the same direction. Beta_ratio = Trial \u03b2 / Tank \u03b2 (values near 1 indicate ",
    "consistent effect magnitude)."
  ))
  doc <- .add_h2(doc, "7.2  Comparison: trial-level model vs tank-level model")
  doc <- .add_ft(doc, .ft_styled(.pretty_names(pseudorep_tbl)))
  doc <- .add_h2(doc, "7.3  Conclusion")
  pr_txt <- if (.pseudorep_consistent)
    "Trial-level vs tank-level treatment directions are CONSISTENT: pseudoreplication is unlikely to drive conclusions."
  else
    "INCONSISTENT directions: interpret trial-level results with caution; tank-level model contradicts trial-level model for at least one analysis."
  doc <- .add_par(doc, pr_txt)

  # =========================================================
  # 8. Descriptive Statistics
  # =========================================================
  doc <- tryCatch({
    doc <- .add_h1(doc, "8.  Descriptive Statistics")
    doc <- .add_par(doc, paste0(
      "Mean ± SEM by treatment and timepoint for all response variables. ",
      "Zone occupancy is expressed as a proportion (0–1); multiply by 100 for %."))

    .sec8_idx <- 0L

    # 8.A Wide-format indicators from df (one row per trial × timepoint)
    .desc_wide <- intersect(c("prop_active","switches_per_session",
                               "zone_flux_per_session","n_main_switches",
                               "mean_nnd_cm","mean_iid_cm",
                               "mean_hull_area_cm2","mean_centroid_spd_cm"),
                             names(df))
    .desc_wide_src <- setNames(rep("df", length(.desc_wide)), .desc_wide)
    if (!is.null(df_gd)) {
      .gd_vars <- intersect(c("mean_nnd_cm","mean_iid_cm",
                               "mean_hull_area_cm2","mean_centroid_spd_cm"),
                             names(df_gd))
      for (.gv in .gd_vars)
        if (!.gv %in% names(df)) .desc_wide_src[[.gv]] <- "gd"
    }

    for (vn in names(.desc_wide_src)) {
      .sec8_idx <- .sec8_idx + 1L
      .src_df <- if (.desc_wide_src[[vn]] == "gd") df_gd else df
      if (!vn %in% names(.src_df)) next
      doc <- .add_h2(doc, paste0("8.", .sec8_idx, "  ", vn))
      tryCatch({
        grp_d <- .src_df %>%
          dplyr::filter(is.finite(.data[[vn]])) %>%
          dplyr::group_by(treatment, timepoint_f) %>%
          dplyr::summarise(
            n    = dplyr::n(),
            Mean = round(mean(.data[[vn]], na.rm = TRUE), 3),
            SEM  = round(.sem(.data[[vn]]), 4),
            .groups = "drop")
        doc <- .add_ft(doc, .ft_styled(.pretty_names(as.data.frame(grp_d))))
      }, error = function(e)
        doc <<- .add_par(doc, paste0("(Error computing descriptive stats for ", vn, ": ", e$message, ")")))
    }

    # 8.B Zone occupancy (long format: one row per trial × timepoint × zone)
    for (.zl in list(
      list(d = df_main_long, lbl = "Main-zone occupancy (prop_time)", grps = c("treatment","timepoint_f","zone")),
      list(d = df_sec_long,  lbl = "Sub-zone occupancy (prop_time)",  grps = c("treatment","timepoint_f","zone"))
    )) {
      if (is.null(.zl$d)) next
      .sec8_idx <- .sec8_idx + 1L
      doc <- .add_h2(doc, paste0("8.", .sec8_idx, "  ", .zl$lbl))
      tryCatch({
        grp_z <- .zl$d %>%
          dplyr::filter(is.finite(prop_time)) %>%
          dplyr::group_by(dplyr::across(dplyr::all_of(.zl$grps))) %>%
          dplyr::summarise(
            n    = dplyr::n(),
            Mean = round(mean(prop_time, na.rm = TRUE), 3),
            SEM  = round(.sem(prop_time), 4),
            .groups = "drop")
        doc <- .add_ft(doc, .ft_styled(.pretty_names(as.data.frame(grp_z))))
      }, error = function(e)
        doc <<- .add_par(doc, paste0("(Error: ", e$message, ")")))
    }
    doc
  }, error = function(e) {
    cat("!! Section 8 (Descriptive Statistics) failed: ", conditionMessage(e), "\n")
    tryCatch(.add_par(doc, paste0("[Section 8 failed: ", conditionMessage(e), "]")),
             error = function(e2) doc)
  })

  # Save
  tryCatch({
    print(doc, target = .rpt_path_step5)
    ts_msg("Report saved: ", .rpt_path_step5)
    tryCatch({
      file.copy(.rpt_path_step5, .rpt_path_main, overwrite = TRUE)
      ts_msg("Report copied: ", .rpt_path_main)
    }, error = function(e)
      warning("Could not copy report to pipeline folder: ", e$message))
    assign("STEP5_REPORT_PATH", .rpt_path_main, envir = .GlobalEnv)
  }, error = function(e) {
    cat("!! Word report print/save FAILED: ", conditionMessage(e), "\n", sep = "")
    warning("Word report print/save failed: ", conditionMessage(e),
            "\nAll CSV outputs still available in: ", STEP5_OUT, call. = FALSE)
  })

} # end if (!is.null(doc))


# =============================================================================
# ==== 9) PUBLISH TO GLOBALENV ================================================
# =============================================================================

assign("res_zone_main",     res_zone_main,     envir = .GlobalEnv)
assign("res_zone_sec",      res_zone_sec,      envir = .GlobalEnv)
assign("res_zone_main_agg", res_zone_main_agg, envir = .GlobalEnv)
assign("res_zone_sec_agg",  res_zone_sec_agg,  envir = .GlobalEnv)
assign("res_switch_agg",  res_switch_agg, envir = .GlobalEnv)
assign("res_switch_tp",   res_switch_tp,  envir = .GlobalEnv)
assign("res_active_tp",   res_active_tp,  envir = .GlobalEnv)
assign("res_active_agg",  res_active_agg, envir = .GlobalEnv)
assign("res_flux_tp",     res_flux_tp,    envir = .GlobalEnv)
assign("res_flux_agg",    res_flux_agg,   envir = .GlobalEnv)
assign("bh_summary",      bh_summary,     envir = .GlobalEnv)
assign("pseudorep_tbl",   pseudorep_tbl,  envir = .GlobalEnv)
assign("STEP5_OUTPUT_DIR",  STEP5_OUT,    envir = .GlobalEnv)

ts_msg("Step 5 (Module A) complete — results in: ", STEP5_OUT)
ts_msg("  School-level analyses | pseudoreplication check | BH correction | Word report | skinny_graphs/")