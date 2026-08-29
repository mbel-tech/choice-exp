################################################################################
# FINAL_easy_script_JUN2026_AICc.R
# -----------------------------------------------------------------------------
# Copy of FINAL_easy_script_JUN2026.R with AICc-based random-effect (RE)
# selection replacing the hardcoded (1|trial_date) / (1|phys_trial) used in
# the original manual script.
#
# Key change: for every indicator, a set of RE candidate structures is fitted
# and the one with the lowest AICc is selected.  The same logic is used in the
# main pipeline (activity_analysis_STATS_choice_exp.R §3).
#
# Why this matters: the original manual script used (1|trial_date) at the
# trial-aggregated level; trial_date is perfectly confounded with tank
# (1 control + 1 EC per date) and absorbed huge between-tank variance, inflating
# SE ~5× and suppressing treatment significance for lr_high.
#
# Diagnostic checks kept:
#   - Shapiro-Wilk on model residuals (reported for Gaussian LMMs)
#   - Levene's test for variance homogeneity
#   - Normality flag when chosen transform still fails SW or Levene (p < 0.05)
#   - Singularity flag on selected LMM
#   - AICc table printed and stored for every indicator
#   - DHARMa-style: Shapiro on raw residuals post-fit
#
# Requires: MuMIn, car  (in addition to packages in the original script)
################################################################################

################################################################################
################# LIBRARIES ####################################################
################################################################################

library(readxl)
library(dplyr)
library(tidyr)
library(purrr)
library(car)        # Levene's test, Type-III ANOVA
library(tibble)
library(lme4)
library(lmerTest)
library(glmmTMB)
library(MuMIn)      # AICc()
library(emmeans)
library(MASS)
library(officer)    # Word document creation
library(flextable)  # Pretty tables in Word

options(contrasts = c("contr.sum", "contr.poly"))


################################################################################
################# IMPORT DATA ##################################################
################################################################################

file_path <- file.path(
  file.path(PROJECT_ROOT, "choice R pipeline/easy_scripts/standalone"),
  "Behavioural and neuroendocrine correlates datasets.xlsx"
)

Preference_exp_behavior <- read_excel(
  path  = file_path,
  sheet = "Preference_exp_behavior"
)

str(Preference_exp_behavior)
head(Preference_exp_behavior)


################################################################################
################# SELECT VARIABLES #############################################
################################################################################

df <- Preference_exp_behavior %>%
  dplyr::select(
    phys_trial,
    interval,
    trial_date,
    treatment,
    tank,
    fish_density,
    motor_side,
    mean_iid_cm,
    mean_nnd_cm,
    mean_school_speed_cm_s,
    mean_school_area_cm2,
    mean_polarisation,
    switches_per_session,
    main_cells_flow_tp,
    sub_cells_high_clr_tp,
    logit_flow,
    lr_high,
    lr_medium,
    lr_low
  ) %>%
  mutate(
    treatment  = as.factor(treatment),
    interval   = as.factor(interval),
    phys_trial = as.factor(phys_trial),
    trial_date = as.factor(trial_date),
    tank       = as.factor(tank),
    motor_side = as.factor(motor_side)
  )

str(df)
head(df)


################################################################################
################# JOIN RAW SUB-ZONE PROPORTIONS (for log-of-means alr) #########
################################################################################
# The Excel sheet carries only the per-interval log-ratios (lr_*) and logit_flow.
# To reproduce the MAIN PIPELINE's trial-level (aggregated) alr exactly, we must
# aggregate the RAW sub-zone proportions across intervals FIRST, then area-
# normalise and take the log-ratio (log-of-means), rather than averaging the
# per-interval log-ratios (mean-of-logs). The two differ by Jensen's inequality
# and the main pipeline uses log-of-means (see activity_analysis_STATS §3,
# df_sec_long_agg -> df_sec_wide_agg).
#
# Raw proportions live only in the canonical easy_scripts dataset. They are joined
# here by (phys_trial, interval) == (phys_trial_id, timepoint_num); this key map
# is verified identical (lr_high matches to < 1e-14 across all 48 rows).
#
# Sub-zone composition: prop_high + prop_medium + prop_low + prop_calm_sec = 1.
#   NOTE: prop_calm_sec is the SUB-ZONE calm (4-way split), NOT prop_calm
#         (= 1 - prop_flow, the main-zone calm). The pipeline uses prop_calm_sec.
.eps_lr <- 1e-4   # Haldane-type floor, identical to main pipeline (.eps_lr)

.canonical_csv <- file.path(PROJECT_ROOT, "choice R pipeline/easy_scripts/easy_scripts_dataset.csv")
.raw_props <- readr::read_csv(.canonical_csv, show_col_types = FALSE) %>%
  dplyr::transmute(
    phys_trial = as.factor(as.character(phys_trial_id)),
    interval   = as.factor(as.character(timepoint_num)),
    prop_high, prop_medium, prop_low, prop_calm_sec, prop_flow
  )

df <- df %>% dplyr::left_join(.raw_props, by = c("phys_trial", "interval"))

cat("\nRaw-proportion join: ",
    sum(!is.na(df$prop_high)), " of ", nrow(df),
    " rows matched (prop_high non-NA)\n", sep = "")
stopifnot(all(!is.na(df$prop_high)))


################################################################################
################# RANDOM-EFFECT CANDIDATE SETS #################################
################################################################################
# Mirrors activity_analysis_STATS_choice_exp.R §3.
#
# Interval-level (N = 48, 3 intervals per physical trial):
#   phys_trial: trial-level random intercept (within-trial correlation)
#   tank:       tank-level random intercept  (shared physical space)
#   trial_date: recording-date random intercept
#   phys_trial_tank: crossed trial + tank REs
#
# Trial-aggregated level (N = 16 physical trials):
#   tank:       main pipeline choice (4 tanks, treatment varies within tank)
#   trial_date: original manual script choice — included here so AICc can
#               reject it empirically if it fits worse

RE_INTERVAL <- c(
  phys_trial       = "(1 | phys_trial)",
  tank             = "(1 | tank)",
  trial_date       = "(1 | trial_date)",
  phys_trial_tank  = "(1 | phys_trial) + (1 | tank)"
)

RE_TRIAL <- c(
  tank       = "(1 | tank)",
  trial_date = "(1 | trial_date)"
)


################################################################################
################# SHARED HELPER FUNCTIONS ######################################
################################################################################

adjust_beta_01 <- function(x) {
  n <- sum(!is.na(x))
  (x * (n - 1) + 0.5) / n
}

# Levene's test — returns p-value (NA if < 2 groups)
levene_p <- function(y, group) {
  d <- data.frame(y = y, g = as.factor(group))[is.finite(y), ]
  if (nlevels(d$g) < 2) return(NA_real_)
  tryCatch(
    car::leveneTest(y ~ g, data = d)$`Pr(>F)`[1],
    error = function(e) NA_real_
  )
}

# Shapiro-Wilk on model residuals
shapiro_residual_p <- function(model) {
  res <- residuals(model)
  res <- res[is.finite(res)]
  if (length(res) < 3 || length(res) > 5000) res <- sample(res, min(length(res), 5000))
  tryCatch(shapiro.test(res)$p.value, error = function(e) NA_real_)
}

# Normality check + optional log1p/sqrt transform (mirrors main pipeline logic)
# Returns list: y_trans, transform_name, sw_p_raw, sw_p_trans,
#               levene_p_raw, levene_p_trans, normality_flag
check_and_transform <- function(y, fixed_str, data) {
  y_valid <- y[is.finite(y)]
  if (length(y_valid) < 8)
    return(list(y_trans = y, transform_name = "none",
                sw_p_raw = NA, sw_p_trans = NA,
                levene_p_raw = NA, levene_p_trans = NA,
                normality_flag = FALSE))

  .fit_resid <- function(y_v, d) {
    d$.y <- y_v
    tryCatch({
      mod <- lm(as.formula(paste(".y ~", fixed_str)), data = d)
      r <- residuals(mod)[is.finite(residuals(mod))]
      if (length(r) > 5000) r <- sample(r, 5000)
      r
    }, error = function(e) y_v - mean(y_v, na.rm = TRUE))
  }

  .sw_p <- function(r) tryCatch(shapiro.test(r)$p.value, error = function(e) NA_real_)

  data_fit <- data[is.finite(y), ]
  .grp <- tryCatch(
    interaction(data_fit[, intersect(c("treatment", "zone"), names(data_fit))], drop = TRUE),
    error = function(e) factor(rep("all", nrow(data_fit)))
  )

  r_raw   <- .fit_resid(y_valid, data_fit)
  sw_raw  <- .sw_p(r_raw)
  lev_raw <- levene_p(y_valid, .grp)

  best_name <- "none"
  best_sw   <- if (is.na(sw_raw)) 1 else sw_raw
  best_lev  <- lev_raw
  y_trans   <- y

  if (!is.na(sw_raw) && sw_raw < 0.05) {
    candidates <- list(
      log1p = log1p(pmax(y_valid, 0)),
      sqrt  = sqrt(pmax(y_valid, 0))
    )
    for (nm in names(candidates)) {
      r_t   <- .fit_resid(candidates[[nm]], data_fit)
      sw_t  <- .sw_p(r_t)
      lev_t <- levene_p(candidates[[nm]], .grp)

      both_new <- !is.na(sw_t) && sw_t > 0.05 && (is.na(lev_t) || lev_t > 0.05)
      both_cur <- !is.na(best_sw) && best_sw > 0.05 && (is.na(best_lev) || best_lev > 0.05)
      improve  <- (both_new && !both_cur) ||
                  (!both_new && !both_cur && !is.na(sw_t) && sw_t > best_sw)
      if (improve) {
        best_sw   <- if (is.na(sw_t)) best_sw else sw_t
        best_lev  <- lev_t
        best_name <- nm
        y_trans[is.finite(y)] <- candidates[[nm]]
      }
    }
  }

  norm_flag <- (!is.na(best_sw) && best_sw < 0.05) ||
               (!is.na(best_lev) && best_lev < 0.05)

  list(y_trans = y_trans, transform_name = best_name,
       sw_p_raw = sw_raw, sw_p_trans = best_sw,
       levene_p_raw = lev_raw, levene_p_trans = best_lev,
       normality_flag = norm_flag)
}

# Gaussian LMM: AICc-based RE selection + normality check
# Returns: list(model, aicc_table, transform_name, sw_p_raw, sw_p_trans,
#               levene_p_trans, normality_flag, re_chosen, singular)
fit_gaussian_aicc <- function(data, response_var, fixed_str, re_cands, label = "",
                              skip_transform = FALSE) {
  y_raw <- data[[response_var]]
  valid <- is.finite(y_raw) & !is.na(data$treatment)
  d_fit <- data[valid, , drop = FALSE]
  y_v   <- y_raw[valid]

  # Normality + transform.
  # skip_transform = TRUE for already-transformed responses (log-ratios, logit):
  #   re-applying log1p/sqrt(pmax()) would clamp legitimate negatives to 0 and
  #   destroy the signal. We still report Shapiro-Wilk on the raw values.
  if (skip_transform) {
    sw_raw <- tryCatch(shapiro.test(y_v)$p.value, error = function(e) NA_real_)
    trans <- list(y_trans = y_v, transform_name = "none",
                  sw_p_raw = sw_raw, sw_p_trans = sw_raw,
                  levene_p_raw = NA_real_, levene_p_trans = NA_real_,
                  normality_flag = (!is.na(sw_raw) && sw_raw < 0.05))
  } else {
    trans <- check_and_transform(y_v, fixed_str, d_fit)
  }
  d_fit$.y <- trans$y_trans

  cat("\n  [", label, "] transform:", trans$transform_name,
      " | SW(raw) p=", round(trans$sw_p_raw, 3),
      " | SW(trans) p=", round(trans$sw_p_trans, 3),
      " | Levene(trans) p=", round(trans$levene_p_trans, 3),
      if (trans$normality_flag) "  *** NORMALITY FLAG ***" else "", "\n")

  # Fit candidates with REML = TRUE (correct for RE comparison, Verbeke & Molenberghs)
  cand_results <- lapply(names(re_cands), function(rn) {
    fmla <- as.formula(paste(".y ~", fixed_str, "+", re_cands[[rn]]))
    m    <- tryCatch(lmerTest::lmer(fmla, data = d_fit, REML = TRUE), error = function(e) NULL)
    if (is.null(m)) return(NULL)
    list(name = rn, re = re_cands[[rn]], model = m,
         aicc = tryCatch(MuMIn::AICc(m), error = function(e) Inf))
  })
  cand_results <- Filter(Negate(is.null), cand_results)

  if (length(cand_results) == 0) {
    warning("No LMM converged for ", label)
    return(NULL)
  }

  aicc_vals <- sapply(cand_results, `[[`, "aicc")
  best_idx  <- which.min(aicc_vals)
  best      <- cand_results[[best_idx]]

  aicc_tbl <- tibble(
    RE            = sapply(cand_results, `[[`, "name"),
    re_formula    = sapply(cand_results, `[[`, "re"),
    fixed_formula = fixed_str,
    AICc          = round(aicc_vals, 2),
    delta_AICc    = round(aicc_vals - min(aicc_vals), 2),
    selected      = seq_along(cand_results) == best_idx
  )

  cat("  AICc table:\n"); print(aicc_tbl, n = Inf)
  cat("  Selected RE:", best$name, "(", best$re, ")\n")

  # Singularity check — try fallback to simpler RE if singular
  singular_orig <- lme4::isSingular(best$model)
  if (singular_orig) {
    cat("  *** SINGULAR FIT — attempting simplification ***\n")
    fallbacks <- setdiff(names(re_cands), best$name)
    for (fb in fallbacks) {
      fmla_fb <- as.formula(paste(".y ~", fixed_str, "+", re_cands[[fb]]))
      m_fb    <- tryCatch(lmerTest::lmer(fmla_fb, data = d_fit, REML = TRUE),
                          error = function(e) NULL)
      if (!is.null(m_fb) && !lme4::isSingular(m_fb)) {
        cat("  Simplification: using", fb, "instead (non-singular)\n")
        best <- list(name = fb, re = re_cands[[fb]], model = m_fb,
                     aicc = tryCatch(MuMIn::AICc(m_fb), error = function(e) Inf))
        break
      }
    }
  }

  singular_final <- lme4::isSingular(best$model)
  if (singular_final)
    cat("  *** SINGULAR FIT (retained — no non-singular simplification found) ***\n")

  list(model          = best$model,
       re_chosen      = best$name,
       re_formula     = best$re,
       transform_name = trans$transform_name,
       sw_p_raw       = trans$sw_p_raw,
       sw_p_trans     = trans$sw_p_trans,
       levene_p_trans = trans$levene_p_trans,
       normality_flag = trans$normality_flag,
       singular       = singular_final,
       aicc_table     = aicc_tbl)
}

# Beta GLMM: AICc-based RE selection
fit_beta_aicc <- function(data, response_var, re_cands, label = "") {
  cand_results <- lapply(names(re_cands), function(rn) {
    fmla <- as.formula(paste(response_var, "~ treatment +", re_cands[[rn]]))
    m    <- tryCatch(glmmTMB(fmla, data = data, family = beta_family(link = "logit")),
                     error = function(e) NULL)
    if (is.null(m)) return(NULL)
    list(name = rn, re = re_cands[[rn]], model = m,
         aicc = tryCatch(MuMIn::AICc(m), error = function(e) Inf))
  })
  cand_results <- Filter(Negate(is.null), cand_results)
  if (length(cand_results) == 0) { warning("No beta GLMM converged for ", label); return(NULL) }

  aicc_vals <- sapply(cand_results, `[[`, "aicc")
  best_idx  <- which.min(aicc_vals)
  best      <- cand_results[[best_idx]]

  aicc_tbl <- tibble(
    RE = sapply(cand_results, `[[`, "name"), re_formula = sapply(cand_results, `[[`, "re"),
    AICc = round(aicc_vals, 2), delta_AICc = round(aicc_vals - min(aicc_vals), 2),
    selected = seq_along(cand_results) == best_idx
  )
  cat("  AICc table (beta):\n"); print(aicc_tbl, n = Inf)
  cat("  Selected RE:", best$name, "\n")
  list(model = best$model, re_chosen = best$name, re_formula = best$re, aicc_table = aicc_tbl)
}

# Beta GLMM — interval-level (needs treatment * interval fixed effects)
fit_beta_aicc_interval <- function(data, response_var, re_cands, label = "") {
  cand_results <- lapply(names(re_cands), function(rn) {
    fmla <- as.formula(paste(response_var, "~ treatment * interval +", re_cands[[rn]]))
    m    <- tryCatch(glmmTMB(fmla, data = data, family = beta_family(link = "logit")),
                     error = function(e) NULL)
    if (is.null(m)) return(NULL)
    list(name = rn, re = re_cands[[rn]], model = m,
         aicc = tryCatch(MuMIn::AICc(m), error = function(e) Inf))
  })
  cand_results <- Filter(Negate(is.null), cand_results)
  if (length(cand_results) == 0) { warning("No beta GLMM converged for ", label); return(NULL) }

  aicc_vals <- sapply(cand_results, `[[`, "aicc")
  best_idx  <- which.min(aicc_vals)
  best      <- cand_results[[best_idx]]

  aicc_tbl <- tibble(
    RE = sapply(cand_results, `[[`, "name"), re_formula = sapply(cand_results, `[[`, "re"),
    AICc = round(aicc_vals, 2), delta_AICc = round(aicc_vals - min(aicc_vals), 2),
    selected = seq_along(cand_results) == best_idx
  )
  cat("  AICc table (beta interval):\n"); print(aicc_tbl, n = Inf)
  cat("  Selected RE:", best$name, "\n")
  list(model = best$model, re_chosen = best$name, re_formula = best$re, aicc_table = aicc_tbl)
}

# Negative-binomial GLMM: AICc RE selection
fit_negbin_aicc <- function(data, response_var, fixed_str, re_cands, label = "") {
  cand_results <- lapply(names(re_cands), function(rn) {
    fmla <- as.formula(paste(response_var, "~", fixed_str, "+", re_cands[[rn]]))
    m    <- tryCatch(glmmTMB(fmla, data = data, family = nbinom2(link = "log")),
                     error = function(e) NULL)
    if (is.null(m)) return(NULL)
    list(name = rn, re = re_cands[[rn]], model = m,
         aicc = tryCatch(MuMIn::AICc(m), error = function(e) Inf))
  })
  cand_results <- Filter(Negate(is.null), cand_results)
  if (length(cand_results) == 0) { warning("No negbin GLMM converged for ", label); return(NULL) }

  aicc_vals <- sapply(cand_results, `[[`, "aicc")
  best_idx  <- which.min(aicc_vals)
  best      <- cand_results[[best_idx]]

  aicc_tbl <- tibble(
    RE = sapply(cand_results, `[[`, "name"), re_formula = sapply(cand_results, `[[`, "re"),
    AICc = round(aicc_vals, 2), delta_AICc = round(aicc_vals - min(aicc_vals), 2),
    selected = seq_along(cand_results) == best_idx
  )
  cat("  AICc table (negbin):\n"); print(aicc_tbl, n = Inf)
  cat("  Selected RE:", best$name, "\n")
  list(model = best$model, re_chosen = best$name, re_formula = best$re, aicc_table = aicc_tbl)
}


################################################################################
################################################################################
################# INTERVAL-LEVEL PIPELINE ######################################
################# MAIN EFFECTS: TREATMENT * INTERVAL ###########################
################################################################################
################################################################################

vars_to_test <- c(
  "mean_iid_cm",
  "mean_nnd_cm",
  "mean_school_speed_cm_s",
  "mean_school_area_cm2",
  "mean_polarisation",
  "switches_per_session",
  "main_cells_flow_tp",
  "sub_cells_high_clr_tp"
)

count_vars     <- c("switches_per_session")
force_beta_vars <- c("mean_polarisation", "main_cells_flow_tp")
force_gaussian_vars <- c("sub_cells_high_clr_tp")   # already transformed; Gaussian

model_list         <- list()
model_choice_table <- list()
anova_list         <- list()
emmeans_list       <- list()
posthoc_list       <- list()
aicc_list          <- list()

for (v in vars_to_test) {

  cat("\n\n============================================================\n")
  cat("Interval-level variable:", v, "\n")
  cat("============================================================\n")

  df_var <- df %>%
    dplyr::select(phys_trial, interval, treatment, tank, trial_date, all_of(v)) %>%
    filter(!is.na(.data[[v]]))

  var_min <- min(df_var[[v]], na.rm = TRUE)
  var_max <- max(df_var[[v]], na.rm = TRUE)
  cat("Range:", round(var_min, 4), "to", round(var_max, 4), "\n")

  ############################################################
  #### CASE 0: COUNT VARIABLE ################################
  ############################################################

  if (v %in% count_vars) {
    cat("Decision: count variable — negative binomial GLMM with AICc RE selection.\n")

    fit_res <- fit_negbin_aicc(df_var, v, "treatment * interval", RE_INTERVAL, label = v)
    if (is.null(fit_res)) next

    model_final  <- fit_res$model
    model_type   <- "Negative binomial GLMM"
    transform    <- "None"
    shapiro_p    <- NA_real_
    norm_flag    <- NA
    re_chosen      <- fit_res$re_chosen
    re_formula_str <- fit_res$re_formula

    anova_final  <- car::Anova(model_final, type = 3)
    emm_final    <- emmeans(model_final, ~ treatment * interval, type = "response")
    posthoc_final <- pairs(emm_final, adjust = "tukey")
    aicc_list[[v]] <- fit_res$aicc_table
  }

  ############################################################
  #### CASE 1: FORCED BETA ###################################
  ############################################################

  else if (v %in% force_beta_vars) {
    cat("Decision: 0-1 proportion — beta GLMM with AICc RE selection.\n")

    response_beta <- paste0(v, "_beta")
    df_var <- df_var %>%
      mutate("{response_beta}" := adjust_beta_01(.data[[v]]))

    fit_res <- fit_beta_aicc_interval(df_var, response_beta, RE_INTERVAL, label = v)
    if (is.null(fit_res)) next

    model_final    <- fit_res$model
    model_type     <- "Beta GLMM"
    transform      <- "Adjusted 0-1 beta"
    shapiro_p      <- NA_real_
    norm_flag      <- NA
    re_chosen      <- fit_res$re_chosen
    re_formula_str <- fit_res$re_formula

    anova_final   <- car::Anova(model_final, type = 3)
    emm_final     <- emmeans(model_final, ~ treatment * interval, type = "response")
    posthoc_final <- pairs(emm_final, adjust = "tukey")
    aicc_list[[v]] <- fit_res$aicc_table
  }

  ############################################################
  #### CASE 2: FORCED GAUSSIAN (already on linear scale) ####
  ############################################################

  else if (v %in% force_gaussian_vars) {
    cat("Decision: forced Gaussian LME (pre-transformed) with AICc RE selection.\n")
    cat("      Response already on log/logit scale -> transform check skipped.\n")

    fit_res <- fit_gaussian_aicc(df_var, v, "treatment * interval",
                                 RE_INTERVAL, label = v, skip_transform = TRUE)
    if (is.null(fit_res)) next

    model_final    <- fit_res$model
    model_type     <- "Gaussian LME"
    transform      <- fit_res$transform_name
    shapiro_p      <- shapiro_residual_p(model_final)
    norm_flag      <- fit_res$normality_flag
    re_chosen      <- fit_res$re_chosen
    re_formula_str <- fit_res$re_formula

    anova_final   <- anova(model_final, type = 3)
    emm_final     <- emmeans(model_final, ~ treatment * interval)
    posthoc_final <- pairs(emm_final, adjust = "tukey")
    aicc_list[[v]] <- fit_res$aicc_table
  }

  ############################################################
  #### CASE 3: MIXED SIGN — GAUSSIAN, NO TRANSFORM ##########
  ############################################################

  else if (var_min < 0 && var_max > 0) {
    cat("Decision: mixed-sign variable — Gaussian LME, no transform, AICc RE selection.\n")

    fit_res <- fit_gaussian_aicc(df_var, v, "treatment * interval",
                                 RE_INTERVAL, label = v)
    if (is.null(fit_res)) next

    model_final    <- fit_res$model
    model_type     <- "Gaussian LME"
    transform      <- fit_res$transform_name
    shapiro_p      <- shapiro_residual_p(model_final)
    norm_flag      <- fit_res$normality_flag
    re_chosen      <- fit_res$re_chosen
    re_formula_str <- fit_res$re_formula

    anova_final   <- anova(model_final, type = 3)
    emm_final     <- emmeans(model_final, ~ treatment * interval)
    posthoc_final <- pairs(emm_final, adjust = "tukey")
    aicc_list[[v]] <- fit_res$aicc_table
  }

  ############################################################
  #### CASE 4: POSITIVE CONTINUOUS — Gaussian + transform ####
  ############################################################

  else if (var_min >= 0 && var_max > 1) {
    cat("Decision: positive continuous — Gaussian LME with transform check, AICc RE selection.\n")

    fit_res <- fit_gaussian_aicc(df_var, v, "treatment * interval",
                                 RE_INTERVAL, label = v)
    if (is.null(fit_res)) next

    model_final    <- fit_res$model
    model_type     <- "Gaussian LME"
    transform      <- fit_res$transform_name
    shapiro_p      <- shapiro_residual_p(model_final)
    norm_flag      <- fit_res$normality_flag
    re_chosen      <- fit_res$re_chosen
    re_formula_str <- fit_res$re_formula

    anova_final   <- anova(model_final, type = 3)
    emm_final     <- emmeans(model_final, ~ treatment * interval)
    posthoc_final <- pairs(emm_final, adjust = "tukey")
    aicc_list[[v]] <- fit_res$aicc_table
  }

  ############################################################
  #### CASE 5: FALLBACK GAUSSIAN #############################
  ############################################################

  else {
    cat("Decision: fallback — Gaussian LME with AICc RE selection.\n")

    fit_res <- fit_gaussian_aicc(df_var, v, "treatment * interval",
                                 RE_INTERVAL, label = v)
    if (is.null(fit_res)) next

    model_final    <- fit_res$model
    model_type     <- "Gaussian LME"
    transform      <- fit_res$transform_name
    shapiro_p      <- shapiro_residual_p(model_final)
    norm_flag      <- fit_res$normality_flag
    re_chosen      <- fit_res$re_chosen
    re_formula_str <- fit_res$re_formula

    anova_final   <- anova(model_final, type = 3)
    emm_final     <- emmeans(model_final, ~ treatment * interval)
    posthoc_final <- pairs(emm_final, adjust = "tukey")
    aicc_list[[v]] <- fit_res$aicc_table
  }

  ############################################################
  #### STORE OUTPUTS #########################################
  ############################################################

  model_list[[v]]     <- model_final
  anova_list[[v]]     <- anova_final
  emmeans_list[[v]]   <- emm_final
  posthoc_list[[v]]   <- posthoc_final

  model_choice_table[[v]] <- tibble(
    variable           = v,
    min_value          = var_min,
    max_value          = var_max,
    model_type         = model_type,
    transformation     = transform,
    re_selected        = re_chosen,
    re_formula         = re_formula_str,
    residual_shapiro_p = shapiro_p,
    normality_flag     = if (is.na(norm_flag)) NA else as.character(norm_flag)
  )

  cat("\nChosen interval-level model for", v, ":\n")
  print(model_choice_table[[v]])
  cat("\nANOVA for", v, ":\n"); print(anova_final)
  cat("\nTukey post hoc for", v, ":\n"); print(posthoc_final)
}

model_decision_table <- bind_rows(model_choice_table)
cat("\n\n=== INTERVAL-LEVEL MODEL DECISION SUMMARY ===\n")
print(model_decision_table)


################################################################################
################################################################################
################# TRIAL-LEVEL DATASET ##########################################
################################################################################
################################################################################

first_non_na <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA) else return(x[1])
}

df_trial <- df %>%
  group_by(phys_trial) %>%
  summarise(
    trial_date             = first_non_na(trial_date),
    treatment              = first_non_na(treatment),
    tank                   = first_non_na(tank),           # needed for RE selection
    fish_density           = first_non_na(fish_density),
    motor_side             = first_non_na(motor_side),
    intervals_included     = paste(sort(unique(as.character(interval))), collapse = ", "),
    n_intervals            = n_distinct(interval),
    n_rows                 = n(),
    mean_iid_cm            = mean(mean_iid_cm,            na.rm = TRUE),
    mean_nnd_cm            = mean(mean_nnd_cm,            na.rm = TRUE),
    mean_school_speed_cm_s = mean(mean_school_speed_cm_s, na.rm = TRUE),
    mean_school_area_cm2   = mean(mean_school_area_cm2,   na.rm = TRUE),
    mean_polarisation      = mean(mean_polarisation,      na.rm = TRUE),
    switches_per_session   = sum(switches_per_session,    na.rm = TRUE),
    # --- LOG-OF-MEANS for alr / logit (matches main pipeline) -----------------
    # Aggregate RAW sub-zone proportions across intervals first ...
    prop_high      = mean(prop_high,     na.rm = TRUE),
    prop_medium    = mean(prop_medium,   na.rm = TRUE),
    prop_low       = mean(prop_low,      na.rm = TRUE),
    prop_calm_sec  = mean(prop_calm_sec, na.rm = TRUE),
    prop_flow      = mean(prop_flow,     na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    # ... then area-normalise (floor at .eps_lr) and take the log-ratio vs calm.
    prop_high_an   = pmax(prop_high     / 10, .eps_lr),
    prop_medium_an = pmax(prop_medium   / 18, .eps_lr),
    prop_low_an    = pmax(prop_low      / 13, .eps_lr),
    prop_calm_an   = pmax(prop_calm_sec / 42, .eps_lr),
    lr_high    = log(prop_high_an   / prop_calm_an),
    lr_medium  = log(prop_medium_an / prop_calm_an),
    lr_low     = log(prop_low_an    / prop_calm_an),
    # logit(prop_flow) on the trial-aggregated flow proportion (no offset; main
    # pipeline drops prop_flow in {0,1}, which does not occur at trial level here)
    logit_flow = log(prop_flow / (1 - prop_flow)),
    treatment  = as.factor(treatment),
    phys_trial = as.factor(phys_trial),
    trial_date = as.factor(trial_date),
    tank       = as.factor(tank)
  )

cat("\nTrial-level dataset: N =", nrow(df_trial), "\n")
cat("tank levels:", nlevels(df_trial$tank),
    "  trial_date levels:", nlevels(df_trial$trial_date), "\n")
print(table(df_trial$treatment))
print(table(df_trial$tank, df_trial$treatment))
print(table(df_trial$trial_date, df_trial$treatment))


################################################################################
################################################################################
################# TRIAL-LEVEL PIPELINE #########################################
################# MAIN EFFECT: TREATMENT #######################################
################################################################################
################################################################################

vars_to_test_trial <- c(
  "mean_iid_cm",
  "mean_nnd_cm",
  "mean_school_speed_cm_s",
  "mean_school_area_cm2",
  "mean_polarisation",
  "switches_per_session",
  "logit_flow",
  "lr_high",
  "lr_medium",
  "lr_low"
)

count_vars_trial        <- c("switches_per_session")
force_beta_vars_trial   <- c("mean_polarisation")
force_gaussian_vars_trial <- c("logit_flow", "lr_high", "lr_medium", "lr_low")

model_list_trial         <- list()
model_choice_table_trial <- list()
anova_list_trial         <- list()
emmeans_list_trial       <- list()
posthoc_list_trial       <- list()
aicc_list_trial          <- list()

for (v in vars_to_test_trial) {

  cat("\n\n============================================================\n")
  cat("Trial-level variable:", v, "\n")
  cat("============================================================\n")

  df_var <- df_trial %>%
    dplyr::select(phys_trial, trial_date, treatment, tank, all_of(v)) %>%
    filter(!is.na(.data[[v]]))

  var_min <- min(df_var[[v]], na.rm = TRUE)
  var_max <- max(df_var[[v]], na.rm = TRUE)
  cat("Range:", round(var_min, 4), "to", round(var_max, 4), "\n")

  ############################################################
  #### CASE 0: COUNT ########################################
  ############################################################

  if (v %in% count_vars_trial) {
    cat("Decision: count — negative binomial GLMM with AICc RE selection.\n")

    fit_res <- fit_negbin_aicc(df_var, v, "treatment", RE_TRIAL, label = v)
    if (is.null(fit_res)) next

    model_final    <- fit_res$model
    model_type     <- "Negative binomial GLMM"
    transform      <- "None"
    shapiro_p      <- NA_real_
    norm_flag      <- NA
    re_chosen      <- fit_res$re_chosen
    re_formula_str <- fit_res$re_formula

    anova_final   <- car::Anova(model_final, type = 3)
    emm_final     <- emmeans(model_final, ~ treatment, type = "response")
    posthoc_final <- pairs(emm_final, adjust = "tukey")
    aicc_list_trial[[v]] <- fit_res$aicc_table
  }

  ############################################################
  #### CASE 1: FORCED BETA ###################################
  ############################################################

  else if (v %in% force_beta_vars_trial) {
    cat("Decision: 0-1 proportion — beta GLMM with AICc RE selection.\n")

    response_beta <- paste0(v, "_beta")
    df_var <- df_var %>%
      mutate("{response_beta}" := adjust_beta_01(.data[[v]]))

    fit_res <- fit_beta_aicc(df_var, response_beta, RE_TRIAL, label = v)
    if (is.null(fit_res)) next

    model_final    <- fit_res$model
    model_type     <- "Beta GLMM"
    transform      <- "Adjusted 0-1 beta"
    shapiro_p      <- NA_real_
    norm_flag      <- NA
    re_chosen      <- fit_res$re_chosen
    re_formula_str <- fit_res$re_formula

    anova_final   <- car::Anova(model_final, type = 3)
    emm_final     <- emmeans(model_final, ~ treatment, type = "response")
    posthoc_final <- pairs(emm_final, adjust = "tukey")
    aicc_list_trial[[v]] <- fit_res$aicc_table
  }

  ############################################################
  #### CASE 2: FORCED GAUSSIAN (alr, logit — pre-transformed)
  ############################################################

  else if (v %in% force_gaussian_vars_trial) {
    cat("Decision: forced Gaussian LME (pre-transformed) with AICc RE selection.\n")
    cat("NOTE: AICc selects between tank and trial_date as RE.\n")
    cat("      Main pipeline uses (1|tank) for aggregated level (AICc-supported).\n")
    cat("      Response already on log/logit scale -> transform check skipped.\n")

    fit_res <- fit_gaussian_aicc(df_var, v, "treatment", RE_TRIAL, label = v,
                                 skip_transform = TRUE)
    if (is.null(fit_res)) next

    model_final    <- fit_res$model
    model_type     <- "Gaussian LME"
    transform      <- fit_res$transform_name
    shapiro_p      <- shapiro_residual_p(model_final)
    norm_flag      <- fit_res$normality_flag
    re_chosen      <- fit_res$re_chosen
    re_formula_str <- fit_res$re_formula

    anova_final   <- anova(model_final, type = 3)
    emm_final     <- emmeans(model_final, ~ treatment)
    posthoc_final <- pairs(emm_final, adjust = "tukey")
    aicc_list_trial[[v]] <- fit_res$aicc_table
  }

  ############################################################
  #### CASE 3: MIXED SIGN ####################################
  ############################################################

  else if (var_min < 0 && var_max > 0) {
    cat("Decision: mixed-sign — Gaussian LME, AICc RE selection.\n")

    fit_res <- fit_gaussian_aicc(df_var, v, "treatment", RE_TRIAL, label = v)
    if (is.null(fit_res)) next

    model_final    <- fit_res$model
    model_type     <- "Gaussian LME"
    transform      <- fit_res$transform_name
    shapiro_p      <- shapiro_residual_p(model_final)
    norm_flag      <- fit_res$normality_flag
    re_chosen      <- fit_res$re_chosen
    re_formula_str <- fit_res$re_formula

    anova_final   <- anova(model_final, type = 3)
    emm_final     <- emmeans(model_final, ~ treatment)
    posthoc_final <- pairs(emm_final, adjust = "tukey")
    aicc_list_trial[[v]] <- fit_res$aicc_table
  }

  ############################################################
  #### CASE 4: POSITIVE CONTINUOUS ###########################
  ############################################################

  else if (var_min >= 0 && var_max > 1) {
    cat("Decision: positive continuous — Gaussian LME, AICc RE selection.\n")

    fit_res <- fit_gaussian_aicc(df_var, v, "treatment", RE_TRIAL, label = v)
    if (is.null(fit_res)) next

    model_final    <- fit_res$model
    model_type     <- "Gaussian LME"
    transform      <- fit_res$transform_name
    shapiro_p      <- shapiro_residual_p(model_final)
    norm_flag      <- fit_res$normality_flag
    re_chosen      <- fit_res$re_chosen
    re_formula_str <- fit_res$re_formula

    anova_final   <- anova(model_final, type = 3)
    emm_final     <- emmeans(model_final, ~ treatment)
    posthoc_final <- pairs(emm_final, adjust = "tukey")
    aicc_list_trial[[v]] <- fit_res$aicc_table
  }

  ############################################################
  #### CASE 5: FALLBACK ######################################
  ############################################################

  else {
    cat("Decision: fallback — Gaussian LME, AICc RE selection.\n")

    fit_res <- fit_gaussian_aicc(df_var, v, "treatment", RE_TRIAL, label = v)
    if (is.null(fit_res)) next

    model_final    <- fit_res$model
    model_type     <- "Gaussian LME"
    transform      <- fit_res$transform_name
    shapiro_p      <- shapiro_residual_p(model_final)
    norm_flag      <- fit_res$normality_flag
    re_chosen      <- fit_res$re_chosen
    re_formula_str <- fit_res$re_formula

    anova_final   <- anova(model_final, type = 3)
    emm_final     <- emmeans(model_final, ~ treatment)
    posthoc_final <- pairs(emm_final, adjust = "tukey")
    aicc_list_trial[[v]] <- fit_res$aicc_table
  }

  ############################################################
  #### STORE OUTPUTS #########################################
  ############################################################

  model_list_trial[[v]]     <- model_final
  anova_list_trial[[v]]     <- anova_final
  emmeans_list_trial[[v]]   <- emm_final
  posthoc_list_trial[[v]]   <- posthoc_final

  model_choice_table_trial[[v]] <- tibble(
    variable           = v,
    min_value          = var_min,
    max_value          = var_max,
    model_type         = model_type,
    transformation     = transform,
    re_selected        = re_chosen,
    re_formula         = re_formula_str,
    residual_shapiro_p = shapiro_p,
    normality_flag     = if (is.na(norm_flag)) NA else as.character(norm_flag)
  )

  cat("\nChosen trial-level model for", v, ":\n")
  print(model_choice_table_trial[[v]])
  cat("\nANOVA for", v, ":\n"); print(anova_final)
  cat("\nTukey post hoc for", v, ":\n"); print(posthoc_final)
}

model_decision_table_trial <- bind_rows(model_choice_table_trial)
cat("\n\n=== TRIAL-LEVEL MODEL DECISION SUMMARY ===\n")
print(model_decision_table_trial)


################################################################################
################# OUTPUT HELPERS ###############################################
################################################################################

output_dir <- file.path(
  file.path(PROJECT_ROOT, "choice R pipeline/easy_scripts/standalone")
)

pretty_indicator_labels <- c(
  mean_iid_cm            = "Mean IID",
  mean_nnd_cm            = "Mean NND",
  mean_school_speed_cm_s = "Mean school speed",
  mean_school_area_cm2   = "Mean school area",
  mean_polarisation      = "Polarisation",
  switches_per_session   = "Switches per session",
  main_cells_flow_tp     = "Main cells flow",
  sub_cells_high_clr_tp  = "Sub-cells high CLR",
  logit_flow             = "Logit flow",
  lr_high                = "High-flow log-ratio",
  lr_medium              = "Medium-flow log-ratio",
  lr_low                 = "Low-flow log-ratio"
)

pretty_variable_name <- function(x)
  ifelse(x %in% names(pretty_indicator_labels), pretty_indicator_labels[x], x)

format_p_statement <- function(p)
  ifelse(is.na(p), "NA",
         ifelse(p < 0.001, "< 0.001", paste0("= ", sprintf("%.3f", p))))

to_subscript_statement <- function(x) {
  s <- as.character(x)
  for (i in 0:9)
    s <- gsub(as.character(i), intToUtf8(0x2080L + i), s, fixed = TRUE)
  s
}

format_numeric_columns <- function(tab, digits = 4)
  tab %>% mutate(across(where(is.numeric), ~ round(.x, digits)))

format_p_columns <- function(tab) {
  p_cols <- names(tab)[grepl("p|Pr\\(", names(tab), ignore.case = TRUE)]
  for (pc in p_cols) {
    if (is.numeric(tab[[pc]])) {
      tab[[pc]] <- ifelse(is.na(tab[[pc]]), NA,
                          ifelse(tab[[pc]] < 0.001, "< 0.001",
                                 sprintf("%.4f", tab[[pc]])))
    }
  }
  tab
}

clean_table_for_markdown <- function(tab) {
  tab <- as.data.frame(tab)
  if ("variable" %in% names(tab)) tab$variable <- pretty_variable_name(tab$variable)
  tab <- format_numeric_columns(tab)
  tab <- format_p_columns(tab)
  tab
}

make_markdown_table <- function(tab) {
  tab <- clean_table_for_markdown(tab)
  if (nrow(tab) == 0) return("_No rows._")
  tab[] <- lapply(tab, function(x) { x <- as.character(x); x[is.na(x)] <- ""; x })
  hdr <- paste0("| ", paste(names(tab), collapse = " | "), " |")
  sep <- paste0("| ", paste(rep("---", ncol(tab)), collapse = " | "), " |")
  body <- apply(tab, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
  paste(c(hdr, sep, body), collapse = "\n")
}

anova_to_table <- function(anova_object) {
  tab <- as.data.frame(anova_object)
  tab$term <- rownames(tab)
  rownames(tab) <- NULL
  tab %>% relocate(term)
}

posthoc_to_table <- function(posthoc_object, significant_only = FALSE) {
  tab <- as.data.frame(posthoc_object)
  if (significant_only && "p.value" %in% names(tab))
    tab <- tab %>% filter(p.value < 0.05)
  tab
}

write_model_markdown_report <- function(
    file_path, title, analysis_description,
    model_decision_table_object, anova_list_object,
    emmeans_list_object, posthoc_list_object, aicc_list_object = NULL,
    other_decision_table = NULL, other_decision_label = NULL) {

  md <- c(paste0("# ", title), "", analysis_description, "",
          paste0("Generated on: ", Sys.time()), "", "---", "",
          "## Model-selection summary", "",
          make_markdown_table(model_decision_table_object), "")

  # Optionally include the companion-level decision table (e.g. trial table
  # shown inside the interval report, and vice versa)
  if (!is.null(other_decision_table)) {
    lbl <- if (!is.null(other_decision_label)) other_decision_label else "Companion analysis model-selection"
    md <- c(md, paste0("### ", lbl), "",
            make_markdown_table(other_decision_table), "")
  }

  for (v in names(anova_list_object)) {
    pretty_v <- pretty_variable_name(v)
    md <- c(md, "", "---", "", paste0("## ", pretty_v), "",
            paste0("**Variable:** `", v, "`"), "")

    # AICc table
    if (!is.null(aicc_list_object) && !is.null(aicc_list_object[[v]])) {
      md <- c(md, "### AICc RE selection", "",
              make_markdown_table(aicc_list_object[[v]]), "")
    }

    # Model decision row
    if ("variable" %in% names(model_decision_table_object)) {
      dr <- model_decision_table_object %>% filter(variable == v)
      md <- c(md, "### Selected model", "", make_markdown_table(dr), "")
    }

    # ANOVA
    if (!is.null(anova_list_object[[v]]))
      md <- c(md, "### Type III ANOVA", "",
              make_markdown_table(anova_to_table(anova_list_object[[v]])), "")

    # EMMs
    if (!is.null(emmeans_list_object[[v]]))
      md <- c(md, "### Estimated marginal means", "",
              make_markdown_table(as.data.frame(emmeans_list_object[[v]])), "")

    # Significant Tukey contrasts
    if (!is.null(posthoc_list_object[[v]])) {
      sig_ph <- posthoc_to_table(posthoc_list_object[[v]], significant_only = TRUE)
      md <- c(md, "### Significant Tukey contrasts (p < 0.05)", "",
              make_markdown_table(sig_ph), "")
    }

    # All Tukey contrasts
    if (!is.null(posthoc_list_object[[v]])) {
      all_ph <- posthoc_to_table(posthoc_list_object[[v]], significant_only = FALSE)
      md <- c(md, "### All Tukey contrasts", "",
              make_markdown_table(all_ph), "")
    }
  }

  writeLines(md, con = file_path, useBytes = TRUE)
  normalizePath(file_path)
}


################################################################################
################# EXPORT REPORTS ###############################################
################################################################################

interval_path <- file.path(output_dir, "interval_treatment_x_interval_model_results.md")
write_model_markdown_report(
  file_path                = interval_path,
  title                    = "Interval-level model results (AICc RE selection)",
  analysis_description     = paste(
    "Fixed: treatment × interval. RE: AICc-selected from:",
    paste(paste0(names(RE_INTERVAL), " = ", RE_INTERVAL), collapse = "; "), ".",
    "Normality: Shapiro-Wilk + Levene on OLS residuals; log1p/sqrt if SW p < 0.05."),
  model_decision_table_object = model_decision_table,
  anova_list_object           = anova_list,
  emmeans_list_object         = emmeans_list,
  posthoc_list_object         = posthoc_list,
  aicc_list_object            = aicc_list,
  other_decision_table        = model_decision_table_trial,
  other_decision_label        = "Trial-aggregated model-selection summary (for reference)"
)

trial_path <- file.path(output_dir, "trial_aggregated_treatment_model_results.md")
write_model_markdown_report(
  file_path                = trial_path,
  title                    = "Trial-aggregated model results (AICc RE selection)",
  analysis_description     = paste(
    "Fixed: treatment. RE: AICc-selected from:",
    paste(paste0(names(RE_TRIAL), " = ", RE_TRIAL), collapse = "; "), ".",
    "NOTE: original manual script used (1|trial_date); trial_date is confounded",
    "with tank (1 control + 1 EC per date), so AICc prefers (1|tank)."),
  model_decision_table_object = model_decision_table_trial,
  anova_list_object           = anova_list_trial,
  emmeans_list_object         = emmeans_list_trial,
  posthoc_list_object         = posthoc_list_trial,
  aicc_list_object            = aicc_list_trial,
  other_decision_table        = model_decision_table,
  other_decision_label        = "Interval-level model-selection summary (for reference)"
)

cat("\nInterval-level report:", interval_path, "\n")
cat("Trial-level report:   ", trial_path, "\n")


################################################################################
################# EXPORT ANOVA STATEMENT CSVs ##################################
################################################################################

format_anova_statement <- function(anova_obj, model_obj = NULL,
                                   indicator, term, term_label) {
  tab <- as.data.frame(anova_obj)
  tab$term <- rownames(tab)

  row <- tab[tab$term == term, , drop = FALSE]
  if (nrow(row) == 0) return(NA_character_)

  if (all(c("NumDF", "DenDF", "F value", "Pr(>F)") %in% names(row))) {
    # round() (not ceiling) so the reported denominator df matches the printed
    # lmerTest ANOVA table and the main pipeline (Satterthwaite DenDF ~ 11.x -> 11)
    df1 <- round(row$NumDF); df2 <- round(row$DenDF)
    paste0(term_label, " F", to_subscript_statement(df1), ",",
           to_subscript_statement(df2), " = ",
           sprintf("%.3f", row$`F value`), ", p ",
           format_p_statement(row$`Pr(>F)`))
  } else if (all(c("Df", "F value", "Pr(>F)") %in% names(row))) {
    df1 <- ceiling(row$Df)
    df2 <- if (!is.null(model_obj)) ceiling(df.residual(model_obj)) else NA
    paste0(term_label, " F", to_subscript_statement(df1), ",",
           to_subscript_statement(df2), " = ",
           sprintf("%.3f", row$`F value`), ", p ",
           format_p_statement(row$`Pr(>F)`))
  } else if (all(c("Chisq", "Df", "Pr(>Chisq)") %in% names(row))) {
    paste0(term_label, " χ²", to_subscript_statement(ceiling(row$Df)),
           " = ", sprintf("%.3f", row$Chisq), ", p ",
           format_p_statement(row$`Pr(>Chisq)`))
  } else NA_character_
}

make_statement_table <- function(anova_list_obj, model_list_obj, term, term_label) {
  indicators <- names(anova_list_obj)
  tibble(
    variable  = indicators,
    indicator = pretty_variable_name(indicators),
    statement = map_chr(indicators, function(v)
      tryCatch(
        format_anova_statement(anova_list_obj[[v]], model_list_obj[[v]],
                               v, term, term_label),
        error = function(e) NA_character_
      )
    )
  )
}

anova_stmt_interval <- make_statement_table(
  anova_list, model_list, "treatment:interval", "Treatment×Interval")
anova_stmt_trial    <- make_statement_table(
  anova_list_trial, model_list_trial, "treatment", "Treatment")

print(anova_stmt_interval)
print(anova_stmt_trial)

csv_interval <- file.path(output_dir, "anova_statement_table_interval.csv")
csv_trial    <- file.path(output_dir, "anova_statement_table_trial.csv")

write_csv_utf8bom <- function(df, path) {
  tmp <- tempfile(fileext = ".csv")
  write.csv(df, tmp, row.names = FALSE, fileEncoding = "UTF-8")
  raw_content <- readBin(tmp, raw(), file.info(tmp)$size)
  writeBin(c(as.raw(c(0xEF, 0xBB, 0xBF)), raw_content), path)
  unlink(tmp)
}

write_csv_utf8bom(anova_stmt_interval, csv_interval)
write_csv_utf8bom(anova_stmt_trial,    csv_trial)

cat("\nInterval ANOVA statements CSV:", csv_interval, "\n")
cat("Trial ANOVA statements CSV:   ", csv_trial,    "\n")


################################################################################
################# EXPORT MODEL DECISION TABLES – WORD ##########################
################################################################################

# Pretty column names for Word table headers
pretty_decision_cols <- c(
  variable           = "Variable",
  min_value          = "Min",
  max_value          = "Max",
  model_type         = "Model type",
  transformation     = "Transformation",
  re_selected        = "RE selected",
  re_formula         = "RE formula",
  residual_shapiro_p = "Shapiro-Wilk p",
  normality_flag     = "Normal?"
)

format_decision_table_for_word <- function(dt) {
  dt2 <- as.data.frame(dt)
  # Pretty indicator labels
  dt2$variable <- pretty_variable_name(dt2$variable)
  # Round numeric columns
  if ("min_value" %in% names(dt2))
    dt2$min_value <- round(dt2$min_value, 3)
  if ("max_value" %in% names(dt2))
    dt2$max_value <- round(dt2$max_value, 3)
  # Format Shapiro p
  if ("residual_shapiro_p" %in% names(dt2)) {
    pv <- suppressWarnings(as.numeric(dt2$residual_shapiro_p))
    dt2$residual_shapiro_p <- ifelse(
      is.na(pv), "",
      ifelse(pv < 0.001, "< 0.001", sprintf("%.4f", pv))
    )
  }
  # Rename columns to pretty headings
  for (old_nm in names(pretty_decision_cols)) {
    new_nm <- pretty_decision_cols[[old_nm]]
    idx <- which(names(dt2) == old_nm)
    if (length(idx) == 1L) names(dt2)[idx] <- new_nm
  }
  dt2
}

make_decision_flextable <- function(dt) {
  formatted <- format_decision_table_for_word(dt)
  ft <- flextable(formatted)
  ft <- bold(ft, part = "header")
  ft <- bg(ft, part = "header", bg = "#2E75B6")
  ft <- color(ft, part = "header", color = "white")
  ft <- fontsize(ft, part = "all", size = 10)
  ft <- font(ft, part = "all", fontname = "Calibri")
  ft <- align(ft, part = "header", align = "center")
  ft <- align(ft, part = "body",   align = "left")
  ft <- set_table_properties(ft, width = 1, layout = "autofit")
  ft
}

docx_path <- file.path(
  file.path(PROJECT_ROOT, "choice R pipeline/easy_scripts/standalone/manual"),
  "FINAL_easy_script_JUN2026_AICc.docx"
)

doc <- read_docx()
doc <- body_add_par(doc, "Model Decision Tables", style = "heading 1")
doc <- body_add_par(doc,
                    paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M")),
                    style = "Normal")
doc <- body_add_par(doc, "", style = "Normal")
doc <- body_add_par(doc, "Interval-level Model Selection", style = "heading 2")
doc <- body_add_flextable(doc, make_decision_flextable(model_decision_table))
doc <- body_add_par(doc, "", style = "Normal")
doc <- body_add_par(doc, "Trial-aggregated Model Selection", style = "heading 2")
doc <- body_add_flextable(doc, make_decision_flextable(model_decision_table_trial))
print(doc, target = docx_path)

cat("\nModel decision tables Word document:", docx_path, "\n")
