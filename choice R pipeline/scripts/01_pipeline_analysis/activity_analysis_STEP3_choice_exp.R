# =============================================================================
# STEP 3 — CHOICE EXPERIMENT
# Linear mixed-effects models for activity and zone occupancy.
#
# ---- TWO ANALYTICAL APPROACHES ----
#
# NOTE: idtracker.ai assigns fish identities independently per segment.
# The same fish_id across segments is NOT the same individual.
# The TANK is the unit repeated across timepoints (same group, not same fish).
# Repeated-measures RE is therefore (1|tank), not a per-fish identifier.
#
#   TIMEPOINT ANALYSIS  — treatment × timepoint interaction
#     Data: df_wide (one row per fish per timepoint); RE candidates include (1|tank)
#     Standard Module A : treatment * timepoint_f + fish_density_f
#     Standard Module B : treatment * timepoint_f
#     Zone    Module A  : treatment * zone * timepoint_f + fish_density_f
#     Zone    Module B  : treatment * zone * timepoint_f
#
#   AGGREGATED ANALYSIS — treatment only (timepoint dropped from fixed effects)
#     Same data (df_wide); timepoint_f omitted from fixed effects.
#     (1|tank) RE accounts for the non-independence of the same tank across timepoints.
#     Standard Module A : treatment * fish_density_f
#     Standard Module B : treatment
#     Zone    Module A  : treatment * zone + fish_density_f
#     Zone    Module B  : treatment * zone
#
# ---- ZONE RESPONSE ----
#     Main zones     : seconds in zone (time_in_flow_s, time_in_calm_s)
#     Secondary zones: area-corrected seconds (area_corr_s_*)
#
# ---- PRIMARY POST-HOC ----
#     Zone models: formal treatment × zone interaction contrasts via
#     emmeans::contrast(em, interaction = c("pairwise","pairwise"))
#
# OUTPUTS:
#   STEP3_stats/STEP3_stats_<timestamp>/
#     timepoint/
#       standard/ main_zones/ sec_zones/
#         module_A(B)/ anova/ posthoc/ model_selection/ diagnostics/ reports/
#     aggregated/
#       standard/ main_zones/ sec_zones/
#         module_A(B)/ ...
# =============================================================================


# =============================================================================
# ==== 0) SETUP ===============================================================
# =============================================================================

options(stringsAsFactors = FALSE, na.action = "na.fail")
set.seed(123)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(lme4)
  library(MuMIn)
  library(car)
  library(emmeans)
  library(multcomp)
  library(multcompView)
  library(officer)
  library(flextable)
  library(readr)
  library(readxl)
  library(stringr)
})

select    <- dplyr::select
filter    <- dplyr::filter
mutate    <- dplyr::mutate
summarise <- dplyr::summarise

ts_msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), " | ", ..., "\n", sep = "")

.get_global <- function(name, default = NULL) {
  if (exists(name, envir = .GlobalEnv, inherits = FALSE))
    get(name, envir = .GlobalEnv, inherits = FALSE)
  else default
}

.step3_dir <- tryCatch({
  frames     <- sys.frames()
  ofile_envs <- Filter(function(f) exists("ofile", envir = f, inherits = FALSE), frames)
  if (length(ofile_envs) > 0)
    normalizePath(dirname(get("ofile", envir = ofile_envs[[length(ofile_envs)]])),
                  winslash = "/", mustWork = FALSE)
  else getwd()
}, error = function(e) getwd())

.find_latest_csv <- function(step_name, csv_filename) {
  parent <- file.path(getwd(), step_name)
  if (!dir.exists(parent)) return(NULL)
  subdirs <- list.dirs(parent, full.names = TRUE, recursive = FALSE)
  subdirs <- subdirs[grepl(paste0("^", step_name, "_\\d{8}_\\d{6}$"), basename(subdirs))]
  if (length(subdirs) == 0) return(NULL)
  latest <- subdirs[which.max(file.mtime(subdirs))]
  candidate <- file.path(latest, csv_filename)
  if (file.exists(candidate)) candidate else NULL
}

.step3_parent <- file.path(getwd(), "STEP3_stats")
if (!dir.exists(.step3_parent))
  dir.create(.step3_parent, recursive = TRUE, showWarnings = FALSE)
STEP3_OUT <- file.path(.step3_parent,
                        paste0("STEP3_stats_", format(Sys.time(), "%Y%m%d_%H%M%S")))
dir.create(STEP3_OUT, recursive = TRUE, showWarnings = FALSE)
ts_msg("Output root: ", STEP3_OUT)

DENSITY_AS_FACTOR_g <- .get_global("DENSITY_AS_FACTOR", TRUE)
DENSITY_LEVELS_g    <- .get_global("DENSITY_LEVELS",    c(4L, 8L, 12L, 16L))
TREATMENT_LEVELS_g  <- .get_global("TREATMENT_LEVELS",  c("control", "exercise choice"))
TIMEPOINT_LEVELS_g  <- .get_global("TIMEPOINT_LEVELS",  1:3)
.ZONE_AREA          <- .get_global("ZONE_AREA_UNITS",   c(high = 10, medium = 18,
                                                           low  = 13, calm   = 42))


# =============================================================================
# ==== 1) LOAD DATA ===========================================================
# =============================================================================

if (!exists("master_fish_by_frame", envir = .GlobalEnv, inherits = FALSE)) {
  ts_msg("master_fish_by_frame not in GlobalEnv; searching disk...")
  .csv_path <- .find_latest_csv("STEP2_output", "master_fish_by_frame_step2.csv")
  if (is.null(.csv_path))
    .csv_path <- .find_latest_csv("STEP1_output", "master_fish_by_frame.csv")
  if (is.null(.csv_path))
    stop("Cannot find master_fish_by_frame. Run Steps 1-2 first.", call. = FALSE)
  master_fish_by_frame <- readr::read_csv(.csv_path, show_col_types = FALSE)
  assign("master_fish_by_frame", master_fish_by_frame, envir = .GlobalEnv)
  ts_msg("Loaded from: ", .csv_path)
}

df_raw <- get("master_fish_by_frame", envir = .GlobalEnv, inherits = FALSE)

if (!exists("TRIAL_META", envir = .GlobalEnv, inherits = FALSE)) {
  .xl <- file.path(.step3_dir, "trial_summary_choice_exp.xlsx")
  if (file.exists(.xl)) {
    .tm <- readxl::read_excel(.xl)
    if ("condition" %in% names(.tm) && !("treatment" %in% names(.tm)))
      .tm <- dplyr::rename(.tm, treatment = condition)
    assign("TRIAL_META", .tm, envir = .GlobalEnv)
  }
}
TRIAL_META <- .get_global("TRIAL_META")


# =============================================================================
# ==== 2) BUILD WIDE MODELLING DATASET ========================================
# =============================================================================
# One row per fish × timepoint.  timepoint is preserved from master_fish_by_frame
# (parsed from session folder name in Step 1).

if (!"timepoint" %in% names(df_raw)) {
  warning("'timepoint' column not found — setting to NA. Re-run Step 1.", call. = FALSE)
  df_raw$timepoint <- NA_integer_
}

df_raw <- df_raw %>%
  dplyr::arrange(trial_id, fish_id, time) %>%
  dplyr::group_by(trial_id, fish_id) %>%
  dplyr::mutate(
    dt_row = time - lag(time),
    dt_row = dplyr::if_else(is.na(dt_row) | dt_row < 0, 0, dt_row)
  ) %>%
  dplyr::ungroup()

to_moving01 <- function(x) {
  x <- tolower(trimws(as.character(x)))
  ifelse(x == "moving", 1L, ifelse(x == "immobile", 0L, NA_integer_))
}

df_wide <- df_raw %>%
  dplyr::mutate(
    moving_TS1 = to_moving01(locomotor_status_TS1),
    is_flow    = as.integer(main_zone == "flow"),
    is_calm    = as.integer(main_zone == "calm"),
    is_high    = as.integer(sec_zone  == "high"),
    is_medium  = as.integer(sec_zone  == "medium"),
    is_low     = as.integer(sec_zone  == "low")
  ) %>%
  dplyr::group_by(
    trial_id, trial_date, fish_id, fish_in_session,
    treatment, fish_density, tank, motor_side, timepoint
  ) %>%
  dplyr::summarise(
    total_time_s          = sum(dt_row, na.rm = TRUE),
    mean_speed_cm_s       = weighted.mean(speed_cm_s, w = dt_row, na.rm = TRUE),
    max_speed_cm_s        = max(speed_cm_s, na.rm = TRUE),
    total_distance_m      = sum(speed_cm_s * dt_row, na.rm = TRUE) / 100,
    moving_time_TS1_s     = sum(moving_TS1  * dt_row, na.rm = TRUE),
    prop_moving_TS1       = moving_time_TS1_s / total_time_s,
    time_in_flow_s        = sum(is_flow   * dt_row, na.rm = TRUE),
    time_in_calm_s        = sum(is_calm   * dt_row, na.rm = TRUE),
    prop_time_in_flow     = time_in_flow_s / total_time_s,
    t_high                = sum(is_high   * dt_row, na.rm = TRUE),
    t_medium              = sum(is_medium * dt_row, na.rm = TRUE),
    t_low                 = sum(is_low    * dt_row, na.rm = TRUE),
    t_calm                = sum(is_calm   * dt_row, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    .nh = t_high   / .ZONE_AREA[["high"]],
    .nm = t_medium / .ZONE_AREA[["medium"]],
    .nl = t_low    / .ZONE_AREA[["low"]],
    .nc = t_calm   / .ZONE_AREA[["calm"]],
    .ns = .nh + .nm + .nl + .nc,
    area_corr_s_flow_high   = .nh,
    area_corr_s_flow_medium = .nm,
    area_corr_s_flow_low    = .nl,
    area_corr_s_calm        = .nc
  ) %>%
  dplyr::select(-.nh, -.nm, -.nl, -.nc, -.ns,
                -t_high, -t_medium, -t_low, -t_calm) %>%
  dplyr::mutate(
    fish_uid    = paste(trial_id, fish_id, sep = "_"),
    treatment   = factor(trimws(tolower(as.character(treatment))),
                         levels = TREATMENT_LEVELS_g),
    fish_density_f = if (isTRUE(DENSITY_AS_FACTOR_g))
      factor(fish_density, levels = DENSITY_LEVELS_g, ordered = TRUE)
    else
      as.numeric(fish_density),
    timepoint_f = factor(timepoint, levels = TIMEPOINT_LEVELS_g),
    trial_id    = as.integer(trial_id),
    tank        = as.character(tank)
  )

if (!"trial_date" %in% names(df_wide) && !is.null(TRIAL_META)) {
  # TRIAL_META now uses 'trial' as the key (not 'trial_id')
  .dm_key <- if ("trial_id" %in% names(TRIAL_META)) "trial_id" else "trial"
  .dm <- dplyr::distinct(TRIAL_META, !!rlang::sym(.dm_key), trial_date) %>%
    dplyr::rename(trial_id = !!rlang::sym(.dm_key)) %>%
    dplyr::mutate(trial_id = as.integer(trial_id))
  df_wide <- dplyr::left_join(df_wide, .dm, by = "trial_id")
}

ts_msg("Wide data (per fish × timepoint): ", nrow(df_wide), " rows; ",
       dplyr::n_distinct(df_wide$tank), " tanks; ",
       dplyr::n_distinct(df_wide$timepoint), " timepoints")
readr::write_csv(df_wide, file.path(STEP3_OUT, "wide_data.csv"))


# =============================================================================
# ==== 3) OUTCOMES ============================================================
# =============================================================================

STANDARD_OUTCOMES <- intersect(
  c("mean_speed_cm_s", "max_speed_cm_s", "total_distance_m",
    "moving_time_TS1_s", "prop_moving_TS1"),
  names(df_wide)
)


# =============================================================================
# ==== 4) LONG-FORMAT ZONE DATASETS ===========================================
# =============================================================================
# Both analyses (timepoint and aggregated) use the same df_wide and df_*_long
# datasets.  The aggregated analysis simply omits timepoint_f from fixed effects;
# the (1|tank) RE accounts for the same tank appearing at multiple timepoints.

# ---- 4a) Main zones -------
df_main_long <- df_wide %>%
  dplyr::select(fish_uid, trial_id, trial_date,
                fish_id, fish_in_session, treatment, fish_density,
                fish_density_f, tank, timepoint, timepoint_f,
                time_in_flow_s, time_in_calm_s) %>%
  tidyr::pivot_longer(
    cols      = c(time_in_flow_s, time_in_calm_s),
    names_to  = "zone", values_to = "occupancy_s"
  ) %>%
  dplyr::mutate(
    zone = dplyr::recode(zone,
      time_in_flow_s = "flow", time_in_calm_s = "calm"),
    zone = factor(zone, levels = c("calm", "flow"))
  )

# ---- 4b) Secondary zones -------
df_sec_long <- df_wide %>%
  dplyr::select(fish_uid, trial_id, trial_date,
                fish_id, fish_in_session, treatment, fish_density,
                fish_density_f, tank, timepoint, timepoint_f,
                area_corr_s_flow_high, area_corr_s_flow_medium,
                area_corr_s_flow_low,  area_corr_s_calm) %>%
  tidyr::pivot_longer(
    cols      = c(area_corr_s_flow_high, area_corr_s_flow_medium,
                  area_corr_s_flow_low,  area_corr_s_calm),
    names_to  = "zone", values_to = "occupancy_s"
  ) %>%
  dplyr::mutate(
    zone = dplyr::recode(zone,
      area_corr_s_flow_high   = "high",  area_corr_s_flow_medium = "medium",
      area_corr_s_flow_low    = "low",   area_corr_s_calm        = "calm"),
    zone = factor(zone, levels = c("calm", "low", "medium", "high"))
  )

ts_msg("Long data: main=", nrow(df_main_long), " sec=", nrow(df_sec_long))
readr::write_csv(df_main_long, file.path(STEP3_OUT, "main_zone_long_data.csv"))
readr::write_csv(df_sec_long,  file.path(STEP3_OUT, "sec_zone_long_data.csv"))


# =============================================================================
# ==== 5) RANDOM EFFECTS CANDIDATES ===========================================
# =============================================================================

# Wide-format RE candidates — used for both timepoint and aggregated analyses.
# (1|tank) is the key repeated-measures RE: the same tank of fish appears at
# each of the 3 timepoints.  For the aggregated analysis (no timepoint in
# fixed effects) (1|tank) still accounts for that non-independence.
RE_WIDE <- c(
  trial        = "(1 | trial_id)",
  tank         = "(1 | tank)",
  date         = "(1 | trial_date)",
  trial_tank   = "(1 | trial_id) + (1 | tank)",
  trial_date_r = "(1 | trial_id) + (1 | trial_date)"
)

# Long-format RE candidates — used for both timepoint and aggregated zone analyses.
# (1|fish_uid) captures within-session correlation (same fish in multiple zones).
# (1|tank) captures between-session correlation (same tank across timepoints).
RE_LONG <- c(
  fish           = "(1 | fish_uid)",
  fish_trial     = "(1 | fish_uid) + (1 | trial_id)",
  fish_tank      = "(1 | fish_uid) + (1 | tank)",
  nested         = "(1 | trial_id / fish_uid)",
  fish_date      = "(1 | fish_uid) + (1 | trial_date)"
)


# =============================================================================
# ==== 6) UTILITY FUNCTIONS ===================================================
# =============================================================================

fmt_p <- function(p) {
  # Report exact value to 3 significant figures when p > 0.0001.
  # Use threshold notation only when p <= 0.0001.
  if (is.na(p)) return("NA")
  if (p <= 0.0001) return("< 0.0001")
  trimws(sub("0+$", "", sub("\\.$", "",
    format(signif(p, 3), scientific = FALSE, trim = TRUE))))
}

select_best_re <- function(response, fixed_str, df, re_cands) {
  results <- lapply(names(re_cands), function(rn) {
    fmla <- as.formula(paste(response, "~", fixed_str, "+", re_cands[[rn]]))
    m    <- tryCatch(lme4::lmer(fmla, data = df, REML = FALSE), error = function(e) NULL)
    if (is.null(m)) return(NULL)
    list(name = rn, re = re_cands[[rn]], model = m,
         aicc = MuMIn::AICc(m), formula = deparse(fmla))
  })
  results <- Filter(Negate(is.null), results)
  if (length(results) == 0) return(NULL)
  aicc_vals <- sapply(results, `[[`, "aicc")
  best_idx  <- which.min(aicc_vals)
  list(
    best_model = results[[best_idx]]$model,
    aicc_table = data.frame(
      re_label     = sapply(results, `[[`, "name"),
      re_formula   = sapply(results, `[[`, "re"),
      full_formula = sapply(results, `[[`, "formula"),
      AICc         = round(aicc_vals, 2),
      delta_AICc   = round(aicc_vals - min(aicc_vals), 2),
      selected     = seq_along(results) == best_idx
    ),
    best_re = results[[best_idx]]$name
  )
}

run_anova_type3 <- function(model, label) {
  tryCatch({
    av  <- car::Anova(model, type = "III")
    out <- as.data.frame(av)
    out$term  <- rownames(out); out$label <- label; rownames(out) <- NULL
    names(out) <- gsub("Chisq", "chisq", names(out))
    names(out) <- gsub("Df",    "df",    names(out))
    names(out) <- gsub("Pr\\(>Chisq\\)", "p_value", names(out))
    out
  }, error = function(e) {
    warning("ANOVA failed for ", label, ": ", conditionMessage(e)); NULL
  })
}

run_posthoc_standard <- function(model, label, with_timepoint = FALSE) {
  out <- list(treatment = NULL, treatment_by_tp = NULL, tp_by_treatment = NULL)
  # Marginal treatment comparison
  tryCatch({
    em  <- emmeans::emmeans(model, ~ treatment, type = "response")
    ctr <- as.data.frame(emmeans::contrast(em, method = "pairwise", adjust = "tukey"))
    ctr$label <- label; out$treatment <- ctr
  }, error = function(e) warning("Post-hoc treatment failed for ", label))
  if (with_timepoint) {
    # Treatment within each timepoint
    tryCatch({
      em_t <- emmeans::emmeans(model, ~ treatment | timepoint_f, type = "response")
      ctr_t <- as.data.frame(emmeans::contrast(em_t, method = "pairwise", adjust = "tukey"))
      ctr_t$label <- label; out$treatment_by_tp <- ctr_t
    }, error = function(e) warning("Post-hoc treatment|timepoint failed for ", label))
    # Timepoint within each treatment
    tryCatch({
      em_tp <- emmeans::emmeans(model, ~ timepoint_f | treatment, type = "response")
      ctr_tp <- as.data.frame(emmeans::contrast(em_tp, method = "pairwise", adjust = "tukey"))
      ctr_tp$label <- label; out$tp_by_treatment <- ctr_tp
    }, error = function(e) warning("Post-hoc timepoint|treatment failed for ", label))
  }
  out
}

run_posthoc_zone <- function(model, label) {
  out <- list(interaction_contrasts = NULL,
              treatment_within_zone = NULL,
              zone_within_treatment = NULL)
  tryCatch({
    em <- emmeans::emmeans(model, ~ treatment * zone, type = "response")
    int_ctr <- as.data.frame(
      emmeans::contrast(em, interaction = c("pairwise", "pairwise"), adjust = "tukey"))
    int_ctr$label <- label; out$interaction_contrasts <- int_ctr
  }, error = function(e)
    warning("Post-hoc interaction contrasts failed for ", label, ": ", conditionMessage(e)))
  tryCatch({
    em_t <- emmeans::emmeans(model, ~ treatment | zone, type = "response")
    ctr_t <- as.data.frame(emmeans::contrast(em_t, method = "pairwise", adjust = "tukey"))
    ctr_t$label <- label; out$treatment_within_zone <- ctr_t
  }, error = function(e)
    warning("Post-hoc treatment|zone failed for ", label, ": ", conditionMessage(e)))
  tryCatch({
    em_z <- emmeans::emmeans(model, ~ zone | treatment, type = "response")
    ctr_z <- as.data.frame(emmeans::contrast(em_z, method = "pairwise", adjust = "tukey"))
    ctr_z$label <- label; out$zone_within_treatment <- ctr_z
  }, error = function(e)
    warning("Post-hoc zone|treatment failed for ", label, ": ", conditionMessage(e)))
  out
}

save_diag_plots <- function(model, label, dir_diag) {
  r <- residuals(model); f <- fitted(model)
  d <- data.frame(fitted = f, resid = r)
  p1 <- ggplot2::ggplot(d, ggplot2::aes(fitted, resid)) +
    ggplot2::geom_point(alpha = 0.4) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = paste(label, "— Residuals vs Fitted"))
  p2 <- ggplot2::ggplot(d, ggplot2::aes(sample = resid)) +
    ggplot2::stat_qq(alpha = 0.4) + ggplot2::stat_qq_line() +
    ggplot2::theme_minimal() + ggplot2::labs(title = paste(label, "— QQ plot"))
  slug <- gsub("[^A-Za-z0-9_]", "_", label)
  ggplot2::ggsave(file.path(dir_diag, paste0(slug, "_resid.png")), p1, width=6, height=4, dpi=150)
  ggplot2::ggsave(file.path(dir_diag, paste0(slug, "_qq.png")),    p2, width=6, height=4, dpi=150)
}

add_ft <- function(doc, df_tbl, heading) {
  doc <- officer::body_add_par(doc, heading, style = "heading 2")
  ft  <- flextable::flextable(df_tbl) %>%
    flextable::fontsize(size = 9, part = "all") %>%
    flextable::autofit() %>%
    flextable::set_table_properties(layout = "autofit", width = 1)
  doc <- flextable::body_add_flextable(doc, ft)
  officer::body_add_par(doc, "", style = "Normal")
}


# =============================================================================
# ==== 7) STANDARD OUTCOMES MODULE ============================================
# =============================================================================

run_standard_module <- function(module_label, fixed_str, df, outcomes,
                                 re_cands, out_root, with_timepoint = FALSE) {
  ts_msg("  [Standard] Module ", module_label, " | ", fixed_str)
  mod_dir  <- file.path(out_root, "standard", paste0("module_", module_label))
  dir_aicc <- file.path(mod_dir, "model_selection"); dir.create(dir_aicc,  recursive = TRUE)
  dir_anov <- file.path(mod_dir, "anova");           dir.create(dir_anov,  recursive = TRUE)
  dir_ph   <- file.path(mod_dir, "posthoc");         dir.create(dir_ph,    recursive = TRUE)
  dir_diag <- file.path(mod_dir, "diagnostics");     dir.create(dir_diag,  recursive = TRUE)
  dir_rep  <- file.path(mod_dir, "reports");         dir.create(dir_rep,   recursive = TRUE)

  aicc_list <- list(); anova_list <- list(); ph_list <- list()

  for (outcome in outcomes) {
    ts_msg("    ", outcome)
    df_fit <- df[is.finite(df[[outcome]]) & !is.na(df$treatment), ]
    if (nrow(df_fit) < 20) { warning("Too few rows for ", outcome); next }
    ms <- select_best_re(outcome, fixed_str, df_fit, re_cands)
    if (is.null(ms)) next
    aicc_list[[outcome]] <- ms$aicc_table
    readr::write_csv(ms$aicc_table, file.path(dir_aicc, paste0(outcome, "_AICc.csv")))
    av <- run_anova_type3(ms$best_model, outcome)
    if (!is.null(av)) {
      anova_list[[outcome]] <- av
      readr::write_csv(av, file.path(dir_anov, paste0(outcome, "_ANOVA.csv")))
    }
    ph <- run_posthoc_standard(ms$best_model, outcome, with_timepoint)
    ph_list[[outcome]] <- ph
    if (!is.null(ph$treatment))
      readr::write_csv(ph$treatment, file.path(dir_ph, paste0(outcome, "_posthoc_treatment.csv")))
    if (!is.null(ph$treatment_by_tp))
      readr::write_csv(ph$treatment_by_tp, file.path(dir_ph, paste0(outcome, "_posthoc_treatment_by_tp.csv")))
    if (!is.null(ph$tp_by_treatment))
      readr::write_csv(ph$tp_by_treatment, file.path(dir_ph, paste0(outcome, "_posthoc_tp_by_treatment.csv")))
    tryCatch(save_diag_plots(ms$best_model, outcome, dir_diag), error = function(e) NULL)
  }

  doc <- officer::read_docx()
  doc <- officer::body_add_par(doc, paste("Standard Outcomes —", module_label),
                                style = "heading 1")
  for (nm in names(anova_list)) {
    doc <- officer::body_add_par(doc, nm, style = "heading 1")
    if (!is.null(aicc_list[[nm]]))
      doc <- add_ft(doc, aicc_list[[nm]], "Model selection (AICc)")
    av <- anova_list[[nm]]; av$p_value <- sapply(av$p_value, fmt_p)
    doc <- add_ft(doc, av, "Type III ANOVA")
    if (!is.null(ph_list[[nm]]$treatment))
      doc <- add_ft(doc, ph_list[[nm]]$treatment, "Tukey — marginal treatment")
    if (!is.null(ph_list[[nm]]$treatment_by_tp))
      doc <- add_ft(doc, ph_list[[nm]]$treatment_by_tp, "Tukey — treatment within timepoint")
    if (!is.null(ph_list[[nm]]$tp_by_treatment))
      doc <- add_ft(doc, ph_list[[nm]]$tp_by_treatment, "Tukey — timepoint within treatment")
  }
  print(doc, target = file.path(dir_rep, paste0("Standard_", module_label, ".docx")))
  ts_msg("    Report written: module ", module_label)
  invisible(list(aicc = aicc_list, anova = anova_list, posthoc = ph_list))
}


# =============================================================================
# ==== 8) ZONE OCCUPANCY MODULE ===============================================
# =============================================================================

run_zone_module <- function(module_label, fixed_str, df_long, re_cands,
                             zone_label, out_root) {
  ts_msg("  [Zone: ", zone_label, "] Module ", module_label, " | ", fixed_str)
  mod_dir  <- file.path(out_root, zone_label, paste0("module_", module_label))
  dir_aicc <- file.path(mod_dir, "model_selection"); dir.create(dir_aicc,  recursive = TRUE)
  dir_anov <- file.path(mod_dir, "anova");           dir.create(dir_anov,  recursive = TRUE)
  dir_ph   <- file.path(mod_dir, "posthoc");         dir.create(dir_ph,    recursive = TRUE)
  dir_diag <- file.path(mod_dir, "diagnostics");     dir.create(dir_diag,  recursive = TRUE)
  dir_rep  <- file.path(mod_dir, "reports");         dir.create(dir_rep,   recursive = TRUE)

  label  <- paste(zone_label, module_label)
  df_fit <- df_long[is.finite(df_long$occupancy_s) &
                      !is.na(df_long$treatment) &
                      !is.na(df_long$zone), ]
  if (nrow(df_fit) < 40) {
    warning("Too few rows for zone model: ", label); return(invisible(NULL))
  }

  ms <- select_best_re("occupancy_s", fixed_str, df_fit, re_cands)
  if (is.null(ms)) return(invisible(NULL))
  readr::write_csv(ms$aicc_table,
                   file.path(dir_aicc, paste0("zone_AICc_", module_label, ".csv")))

  av <- run_anova_type3(ms$best_model, label)
  if (!is.null(av)) {
    av$key_term <- grepl("treatment.*zone|zone.*treatment", av$term, ignore.case = TRUE)
    readr::write_csv(av, file.path(dir_anov, paste0("zone_ANOVA_", module_label, ".csv")))
  }

  ph <- run_posthoc_zone(ms$best_model, label)
  if (!is.null(ph$interaction_contrasts))
    readr::write_csv(ph$interaction_contrasts,
                     file.path(dir_ph, paste0("interaction_contrasts_", module_label, ".csv")))
  if (!is.null(ph$treatment_within_zone))
    readr::write_csv(ph$treatment_within_zone,
                     file.path(dir_ph, paste0("treatment_within_zone_", module_label, ".csv")))
  if (!is.null(ph$zone_within_treatment))
    readr::write_csv(ph$zone_within_treatment,
                     file.path(dir_ph, paste0("zone_within_treatment_", module_label, ".csv")))

  tryCatch(save_diag_plots(ms$best_model, label, dir_diag), error = function(e) NULL)

  tryCatch({
    em_int <- emmeans::emmeans(ms$best_model, ~ treatment * zone, type = "response")
    readr::write_csv(as.data.frame(em_int),
                     file.path(mod_dir, "emmeans_treatment_x_zone.csv"))
  }, error = function(e) NULL)

  doc <- officer::read_docx()
  doc <- officer::body_add_par(doc,
    paste("Zone Occupancy —", zone_label, "—", module_label), style = "heading 1")
  doc <- officer::body_add_par(doc, paste("Fixed effects:", fixed_str), style = "Normal")
  doc <- officer::body_add_par(doc,
    paste("Best RE:", ms$best_re, "|",
          ms$aicc_table$re_formula[ms$aicc_table$selected]), style = "Normal")
  doc <- officer::body_add_par(doc, "", style = "Normal")
  doc <- add_ft(doc, ms$aicc_table, "Model selection (AICc)")

  if (!is.null(av)) {
    av_rpt <- av; av_rpt$p_value <- sapply(av_rpt$p_value, fmt_p)
    doc <- officer::body_add_par(doc,
      "Type III ANOVA — key term: treatment:zone interaction", style = "heading 2")
    ft_av <- flextable::flextable(av_rpt) %>%
      flextable::fontsize(size = 9, part = "all") %>%
      flextable::autofit() %>%
      flextable::set_table_properties(layout = "autofit", width = 1) %>%
      flextable::bold(
        i = ~ grepl("treatment.*zone|zone.*treatment", term, ignore.case = TRUE),
        bold = TRUE)
    doc <- flextable::body_add_flextable(doc, ft_av)
    doc <- officer::body_add_par(doc, "", style = "Normal")
  }

  if (!is.null(ph$interaction_contrasts))
    doc <- add_ft(doc, ph$interaction_contrasts,
                  "Post-hoc: treatment \u00d7 zone interaction contrasts (primary)")
  if (!is.null(ph$treatment_within_zone))
    doc <- add_ft(doc, ph$treatment_within_zone,
                  "Post-hoc: treatment within each zone (supplementary)")
  if (!is.null(ph$zone_within_treatment))
    doc <- add_ft(doc, ph$zone_within_treatment,
                  "Post-hoc: zone within each treatment (supplementary)")

  print(doc, target = file.path(dir_rep,
    paste0("ZoneOcc_", zone_label, "_", module_label, ".docx")))
  ts_msg("    Report written: ", zone_label, " ", module_label)

  invisible(list(aicc = ms$aicc_table, anova = av, posthoc = ph,
                 best_model = ms$best_model))
}


# =============================================================================
# ==== 9) RUN — TIMEPOINT ANALYSIS ============================================
# =============================================================================

.out_tp <- file.path(STEP3_OUT, "timepoint")
dir.create(.out_tp, recursive = TRUE)
ts_msg("=== TIMEPOINT ANALYSIS (treatment * timepoint) ===")

# ---- Standard outcomes -------------------------------------------------------
fixed_std_tp_A <- if (isTRUE(DENSITY_AS_FACTOR_g))
  "treatment * timepoint_f + fish_density_f" else
  "treatment * timepoint_f + fish_density"

std_tp_A <- run_standard_module("A_treatment_density", fixed_std_tp_A,
                                 df_wide, STANDARD_OUTCOMES, RE_WIDE,
                                 .out_tp, with_timepoint = TRUE)
std_tp_B <- run_standard_module("B_treatment_only", "treatment * timepoint_f",
                                 df_wide, STANDARD_OUTCOMES, RE_WIDE,
                                 .out_tp, with_timepoint = TRUE)

# ---- Main zone occupancy ------------------------------------------------------
fixed_main_tp_A <- if (isTRUE(DENSITY_AS_FACTOR_g))
  "treatment * zone * timepoint_f + fish_density_f" else
  "treatment * zone * timepoint_f + fish_density"

ts_msg("=== ZONE OCCUPANCY: MAIN ZONES (timepoint) ===")
main_tp_A <- run_zone_module("A_treatment_density", fixed_main_tp_A,
                              df_main_long, RE_LONG, "main_zones", .out_tp)
main_tp_B <- run_zone_module("B_treatment_only", "treatment * zone * timepoint_f",
                              df_main_long, RE_LONG, "main_zones", .out_tp)

# ---- Secondary zone occupancy -------------------------------------------------
fixed_sec_tp_A <- if (isTRUE(DENSITY_AS_FACTOR_g))
  "treatment * zone * timepoint_f + fish_density_f" else
  "treatment * zone * timepoint_f + fish_density"

ts_msg("=== ZONE OCCUPANCY: SECONDARY ZONES (timepoint) ===")
sec_tp_A <- run_zone_module("A_treatment_density", fixed_sec_tp_A,
                             df_sec_long, RE_LONG, "sec_zones", .out_tp)
sec_tp_B <- run_zone_module("B_treatment_only", "treatment * zone * timepoint_f",
                             df_sec_long, RE_LONG, "sec_zones", .out_tp)


# =============================================================================
# ==== 10) RUN — AGGREGATED ANALYSIS ==========================================
# =============================================================================

.out_agg <- file.path(STEP3_OUT, "aggregated")
dir.create(.out_agg, recursive = TRUE)
ts_msg("=== AGGREGATED ANALYSIS (timepoint dropped from fixed effects) ===")
# Same df_wide and df_*_long datasets; timepoint_f omitted from fixed effects.
# (1|tank) from RE_WIDE accounts for the same tank appearing at multiple timepoints.

# ---- Standard outcomes -------------------------------------------------------
fixed_std_agg_A <- if (isTRUE(DENSITY_AS_FACTOR_g))
  "treatment * fish_density_f" else "treatment * fish_density"

std_agg_A <- run_standard_module("A_treatment_density", fixed_std_agg_A,
                                  df_wide, STANDARD_OUTCOMES, RE_WIDE,
                                  .out_agg, with_timepoint = FALSE)
std_agg_B <- run_standard_module("B_treatment_only", "treatment",
                                  df_wide, STANDARD_OUTCOMES, RE_WIDE,
                                  .out_agg, with_timepoint = FALSE)

# ---- Main zone occupancy ------------------------------------------------------
fixed_main_agg_A <- if (isTRUE(DENSITY_AS_FACTOR_g))
  "treatment * zone + fish_density_f" else "treatment * zone + fish_density"

ts_msg("=== ZONE OCCUPANCY: MAIN ZONES (aggregated) ===")
main_agg_A <- run_zone_module("A_treatment_density", fixed_main_agg_A,
                               df_main_long, RE_LONG, "main_zones", .out_agg)
main_agg_B <- run_zone_module("B_treatment_only", "treatment * zone",
                               df_main_long, RE_LONG, "main_zones", .out_agg)

# ---- Secondary zone occupancy -------------------------------------------------
fixed_sec_agg_A <- if (isTRUE(DENSITY_AS_FACTOR_g))
  "treatment * zone + fish_density_f" else "treatment * zone + fish_density"

ts_msg("=== ZONE OCCUPANCY: SECONDARY ZONES (aggregated) ===")
sec_agg_A <- run_zone_module("A_treatment_density", fixed_sec_agg_A,
                              df_sec_long, RE_LONG, "sec_zones", .out_agg)
sec_agg_B <- run_zone_module("B_treatment_only", "treatment * zone",
                              df_sec_long, RE_LONG, "sec_zones", .out_agg)


# =============================================================================
# ==== 11) PUBLISH TO GLOBALENV ===============================================
# =============================================================================

assign("step3_tp_standard_A",   std_tp_A,  envir = .GlobalEnv)
assign("step3_tp_standard_B",   std_tp_B,  envir = .GlobalEnv)
assign("step3_tp_main_zones_A", main_tp_A, envir = .GlobalEnv)
assign("step3_tp_main_zones_B", main_tp_B, envir = .GlobalEnv)
assign("step3_tp_sec_zones_A",  sec_tp_A,  envir = .GlobalEnv)
assign("step3_tp_sec_zones_B",  sec_tp_B,  envir = .GlobalEnv)

assign("step3_agg_standard_A",   std_agg_A,  envir = .GlobalEnv)
assign("step3_agg_standard_B",   std_agg_B,  envir = .GlobalEnv)
assign("step3_agg_main_zones_A", main_agg_A, envir = .GlobalEnv)
assign("step3_agg_main_zones_B", main_agg_B, envir = .GlobalEnv)
assign("step3_agg_sec_zones_A",  sec_agg_A,  envir = .GlobalEnv)
assign("step3_agg_sec_zones_B",  sec_agg_B,  envir = .GlobalEnv)

assign("df_wide_step3",      df_wide,      envir = .GlobalEnv)
assign("df_main_long_step3", df_main_long, envir = .GlobalEnv)
assign("df_sec_long_step3",  df_sec_long,  envir = .GlobalEnv)
assign("STEP3_OUTPUT_DIR",      STEP3_OUT,       envir = .GlobalEnv)

ts_msg("Step 3 complete — timepoint + aggregated results in: ", STEP3_OUT)
