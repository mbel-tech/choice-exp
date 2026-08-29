# =============================================================================
# STEP 5 (Module B) — CHOICE EXPERIMENT  |  TRIAL-LEVEL STATISTICS
# Companion to activity_analysis_STATS_choice_exp.R (Module A).
#
# Unit of analysis: TRIAL (one row per trial_id), timepoints collapsed.
# Fixed effect: treatment.  Random effect: (1 | tank) when multiple trials per
# tank are available; otherwise none (a plain LM is used as fallback).
#
# AGGREGATION RULES:
#   - prop_active, prop_time_*, NND, polarisation, hull area, centroid speed:
#       mean over timepoints within trial.
#   - switches_per_session and zone_flux_per_session: mean over timepoints
#       (each timepoint is already normalised to a 20-min session, so the mean
#       is itself a 20-min-session-comparable quantity).
#   - n_main_switches (raw count): summed over timepoints; written for
#       diagnostics only (not used as a statistical response).
#   - long-format prop_time: mean per (trial_id, zone_level, zone).
#
# Reuses: run_lmm_analysis, fmt_p, .get_treatment_p, .find_latest_csv from
# Module A — Module A must therefore have been sourced first in the same
# session (the master pipeline guarantees this order).
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lme4)
  library(MASS)
  library(MuMIn)
  library(car)
  library(emmeans)
  library(multcomp)
  library(multcompView)
  library(readr)
})

# Guard: helpers from Module A must be present
if (!exists("run_lmm_analysis", envir = .GlobalEnv, inherits = FALSE) ||
    !exists("fmt_p",             envir = .GlobalEnv, inherits = FALSE)) {
  stop("Module B requires Module A helpers. Source Module A first.", call. = FALSE)
}

ts_msg_b <- function(...) cat(format(Sys.time(), "%H:%M:%S"), " | [Mod B] ", ..., "\n", sep = "")

# ---- 1) Output dir ----------------------------------------------------------
.parent_b <- file.path(getwd(), "STEP5_stats_trial")
if (!dir.exists(.parent_b)) dir.create(.parent_b, recursive = TRUE, showWarnings = FALSE)
STEP5B_OUT <- file.path(.parent_b,
                       paste0("STEP5_stats_trial_",
                              format(Sys.time(), "%Y%m%d_%H%M%S")))
dir.create(STEP5B_OUT, recursive = TRUE, showWarnings = FALSE)

# run_lmm_analysis (defined in Module A) writes under STEP5_OUT — temporarily
# override the global so Module B output lands in STEP5B_OUT, then restore.
.prev_STEP5_OUT <- if (exists("STEP5_OUT", envir = .GlobalEnv, inherits = FALSE))
  get("STEP5_OUT", envir = .GlobalEnv, inherits = FALSE) else NULL
STEP5_OUT <- STEP5B_OUT
assign("STEP5_OUT", STEP5_OUT, envir = .GlobalEnv)
ts_msg_b("Output root: ", STEP5_OUT)

# ---- 2) Collapse trial x timepoint -> trial ---------------------------------
df_tp <- get("trial_activity_summary", envir = .GlobalEnv, inherits = FALSE)

.mean_cols <- intersect(c("prop_active", "switches_per_session",
                          "zone_flux_per_session",
                          "prop_time_in_flow", "prop_time_in_calm",
                          "prop_time_in_high", "prop_time_in_medium",
                          "prop_time_in_low",
                          "mean_nnd_cm", "mean_polarisation", "mean_iid_cm",
                          "mean_hull_area_cm2", "mean_centroid_spd_cm",
                          "bl_cm_trial"),
                        names(df_tp))
.sum_cols <- intersect(c("n_main_switches"), names(df_tp))

.meta_keep <- intersect(c("tank", "treatment", "fish_density",
                          "fish_density_f", "trial_date"), names(df_tp))

df_b <- df_tp %>%
  dplyr::group_by(trial_id) %>%
  dplyr::summarise(
    dplyr::across(dplyr::all_of(.meta_keep), ~ dplyr::first(.x)),
    dplyr::across(dplyr::all_of(.mean_cols), ~ mean(.x, na.rm = TRUE)),
    dplyr::across(dplyr::all_of(.sum_cols),  ~ sum(.x,  na.rm = TRUE)),
    n_timepoints = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    treatment = factor(trimws(tolower(as.character(treatment))),
                       levels = .get_global("TREATMENT_LEVELS",
                                            c("control", "exercise choice"))),
    tank      = as.character(tank)
  )

readr::write_csv(df_b, file.path(STEP5_OUT, "analysis_ready_trial.csv"))

# Long-format occupancy collapsed to trial level
df_long_tp <- get("trial_occupancy_long", envir = .GlobalEnv, inherits = FALSE)
df_long_b  <- df_long_tp %>%
  dplyr::group_by(trial_id, tank, treatment, zone_level, zone) %>%
  dplyr::summarise(prop_time = mean(prop_time, na.rm = TRUE),
                   .groups = "drop") %>%
  dplyr::mutate(
    treatment  = factor(trimws(tolower(as.character(treatment))),
                        levels = .get_global("TREATMENT_LEVELS",
                                             c("control", "exercise choice"))),
    zone_level = factor(zone_level, levels = c("main", "sec")),
    zone       = factor(zone, levels = c("flow", "calm", "high", "medium", "low")),
    tank       = as.character(tank)
  )

df_main_b <- dplyr::filter(df_long_b, zone_level == "main") %>%
  dplyr::mutate(zone = factor(as.character(zone), levels = c("flow", "calm")))
df_sec_b  <- dplyr::filter(df_long_b, zone_level == "sec") %>%
  dplyr::mutate(zone = factor(as.character(zone),
                              levels = c("high", "medium", "low", "calm")))

ts_msg_b("Trial-level data: ", nrow(df_b), " trials | ",
         dplyr::n_distinct(df_b$tank), " tanks")

# ---- 3) RE candidates -------------------------------------------------------
# Use (1|tank) only when multiple trials share a tank.
.tanks_with_multi <- sum(table(df_b$tank) > 1)
RE_TRIAL <- if (.tanks_with_multi >= 2) c(tank = "(1 | tank)") else c()

if (length(RE_TRIAL) == 0) {
  ts_msg_b("Each tank has \u2264 1 trial — fitting plain LMs (no RE).")
  # Provide a dummy RE that lme4 will collapse harmlessly; or fall back to lm.
  # Simpler: run LMs via lm() wrapped in a small adapter.
}

# ---- 4) Helper: trial-level LM fallback when no RE candidate ----------------
.run_lm_trial <- function(label, data, response, fixed_str, focal_terms = NULL) {
  out_dir <- file.path(STEP5_OUT, label)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  ts_msg_b("[", label, "] (LM)")

  d <- data[is.finite(data[[response]]) & !is.na(data$treatment), , drop = FALSE]
  if (nrow(d) < 6) { warning("Too few rows for ", label); return(invisible(NULL)) }
  d$.y <- d[[response]]

  m <- tryCatch(lm(as.formula(paste(".y ~", fixed_str)), data = d),
                error = function(e) NULL)
  if (is.null(m)) { warning("LM fit failed for ", label); return(invisible(NULL))}

  av <- tryCatch({
    a <- car::Anova(m, type = "III")
    dav <- as.data.frame(a); dav$term <- rownames(dav); rownames(dav) <- NULL
    names(dav) <- c("ss", "df", "F_value", "p_value", "term")
    dav$chisq <- dav$F_value  # alias so downstream caption code works
    dav
  }, error = function(e) NULL)
  if (!is.null(av)) readr::write_csv(av, file.path(out_dir, "anova.csv"))

  caps <- list()
  if (!is.null(focal_terms) && !is.null(av)) {
    for (nm in names(focal_terms)) {
      row <- av[grepl(focal_terms[[nm]], av$term, ignore.case = TRUE, perl = TRUE), ]
      if (nrow(row) > 0)
        caps[[nm]] <- paste0("F(", row$df[1], ",", m$df.residual, ") = ",
                              signif(row$F_value[1], 4), ", ", fmt_p(row$p_value[1]))
    }
  }

  invisible(list(
    model         = m,
    fixed_formula = fixed_str,
    transform     = "LM",
    aicc_table    = NULL,
    anova         = av,
    anova_caps    = caps,
    posthoc       = list(),
    data_fit      = d
  ))
}

.fit <- function(label, data, response, fixed_str, focal_terms = NULL) {
  if (length(RE_TRIAL) > 0) {
    run_lmm_analysis(label = label, data = data, response = response,
                     fixed_str = fixed_str, re_cands = RE_TRIAL,
                     focal_terms = focal_terms)
  } else {
    .run_lm_trial(label, data, response, fixed_str, focal_terms)
  }
}

# ---- 5) Analyses (trial-level) ----------------------------------------------
ts_msg_b("=== B1: Main-zone occupancy x Treatment x Zone (trial) ===")
res_b_zone_main <- .fit("zone_main", df_main_b, "prop_time",
                        "treatment * zone",
                        list("Treatment"      = "^treatment$",
                             "Zone"           = "^zone$",
                             "Treatment:Zone" = "treatment.*zone|zone.*treatment"))

ts_msg_b("=== B2: Sub-zone occupancy x Treatment x Zone (trial) ===")
res_b_zone_sec <- .fit("zone_sec", df_sec_b, "prop_time",
                       "treatment * zone",
                       list("Treatment"      = "^treatment$",
                            "Zone"           = "^zone$",
                            "Treatment:Zone" = "treatment.*zone|zone.*treatment"))

ts_msg_b("=== B3: switches_per_session x Treatment (trial) ===")
res_b_switch <- .fit("switches", df_b, "switches_per_session",
                     "treatment", list("Treatment" = "^treatment$"))

ts_msg_b("=== B4: prop_active x Treatment (trial) ===")
res_b_active <- .fit("active", df_b, "prop_active",
                     "treatment", list("Treatment" = "^treatment$"))

# Group dynamics (when columns available)
.has <- function(x) x %in% names(df_b)

res_b_nnd  <- if (.has("mean_nnd_cm"))
  .fit("nnd", df_b, "mean_nnd_cm", "treatment",
       list("Treatment" = "^treatment$")) else NULL
res_b_pol  <- if (.has("mean_polarisation"))
  .fit("polarisation", df_b, "mean_polarisation", "treatment",
       list("Treatment" = "^treatment$")) else NULL
res_b_hull <- if (.has("mean_hull_area_cm2"))
  .fit("hull_area", df_b, "mean_hull_area_cm2", "treatment",
       list("Treatment" = "^treatment$")) else NULL
res_b_cspd <- if (.has("mean_centroid_spd_cm"))
  .fit("centroid_speed", df_b, "mean_centroid_spd_cm", "treatment",
       list("Treatment" = "^treatment$")) else NULL

# ---- 6) BH correction across trial-level analyses ---------------------------
.bh_row_b <- function(label, res, pat = "^treatment$") {
  data.frame(Analysis = label,
             p_raw    = .get_treatment_p(res, pat),
             stringsAsFactors = FALSE)
}

bh_b <- dplyr::bind_rows(
  .bh_row_b("B1 Main-zone x Treatment",  res_b_zone_main),
  .bh_row_b("B2 Sub-zone x Treatment",   res_b_zone_sec),
  .bh_row_b("B3 Switches x Treatment",   res_b_switch),
  .bh_row_b("B4 prop_active x Treatment",res_b_active),
  .bh_row_b("B5 NND x Treatment",        res_b_nnd),
  .bh_row_b("B6 Polarisation x Treatment", res_b_pol),
  .bh_row_b("B7 School area x Treatment",  res_b_hull),
  .bh_row_b("B8 School speed x Treatment", res_b_cspd)
)
bh_b$p_BH    <- p.adjust(bh_b$p_raw, method = "BH")
bh_b$sig_raw <- !is.na(bh_b$p_raw) & bh_b$p_raw < 0.05
bh_b$sig_BH  <- !is.na(bh_b$p_BH)  & bh_b$p_BH  < 0.05
bh_b$p_raw   <- round(bh_b$p_raw, 4)
bh_b$p_BH    <- round(bh_b$p_BH,  4)

readr::write_csv(bh_b, file.path(STEP5_OUT, "bh_correction_summary_trial.csv"))

# ---- 7) Publish to GlobalEnv ------------------------------------------------
assign("res_b_zone_main", res_b_zone_main, envir = .GlobalEnv)
assign("res_b_zone_sec",  res_b_zone_sec,  envir = .GlobalEnv)
assign("res_b_switch",    res_b_switch,    envir = .GlobalEnv)
assign("res_b_active",    res_b_active,    envir = .GlobalEnv)
assign("res_b_nnd",       res_b_nnd,       envir = .GlobalEnv)
assign("res_b_pol",       res_b_pol,       envir = .GlobalEnv)
assign("res_b_hull",      res_b_hull,      envir = .GlobalEnv)
assign("res_b_cspd",      res_b_cspd,      envir = .GlobalEnv)
assign("bh_summary_trial", bh_b,           envir = .GlobalEnv)
assign("STEP5B_OUTPUT_DIR", STEP5B_OUT,    envir = .GlobalEnv)

# Restore Module A's STEP5_OUT so a re-run of Module A in the same session
# does not write into the Module B directory.
if (is.null(.prev_STEP5_OUT)) {
  rm("STEP5_OUT", envir = .GlobalEnv)
} else {
  assign("STEP5_OUT", .prev_STEP5_OUT, envir = .GlobalEnv)
}

ts_msg_b("Module B complete \u2014 results in: ", STEP5B_OUT)
