# =============================================================================
# sensitivity_forced_trial_re.R
# =============================================================================
# Reviewer 1, comment 7 (robustness_analysis_memo.md, 2026-08-08): AICc-based
# selection of the random-effects structure is not label-blind in its
# consequences -- when real between-trial variance exists, AICc sometimes
# drops the trial random intercept, and every time it does, the analysis
# reverts to treating replicates from the same trial as independent. The
# memo's simulation (matched to this design) shows this inflates the
# false-positive rate from a nominal 5% up to 9% at realistic between-trial
# SD, while always retaining trial is calibrated (slightly conservative).
#
# NOTE on scope: this check is only meaningful where the model's rows are NOT
# already one-per-trial. The behavioural trial-level (Level A, N=16) models
# ARE already one row per trial by construction (aggregated), so there is no
# within-trial replication left for a trial random intercept to absorb --
# lme4 cannot even fit it (grouping factor levels == N). The memo's own
# framing ("matters most for monoamine and cortisol models, where fish
# within a trial are the replicates") is therefore applied here to:
#   (a) behavioural INTERVAL-level (Level B, N=48, 3 sessions/trial) models,
#   (b) cortisol (fish nested in trial, n=78),
#   (c) monoamine headline cells (fish nested in trial, n=29-37/cell).
#
# Output: sensitivity_forced_trial_re_<domain>.csv per domain.
# =============================================================================

suppressMessages({
  library(data.table)
  library(lme4)
  library(lmerTest)
  library(performance)
})

STEP5 <- file.path(PROJECT_ROOT, "choice R pipeline/output/STEP5_stats/STEP5_stats_20260808_123821")
set.seed(20260706)

.eta2 <- function(F, df1, df2) (F * df1) / (F * df1 + df2)

.forced_trial_check <- function(data, response, fixed_extra, label, trial_col = "trial") {
  fx <- if (nchar(fixed_extra)) paste("treatment +", fixed_extra) else "treatment"
  f_forced <- as.formula(paste(response, "~", fx, "+ (1|", trial_col, ")"))
  m_forced <- tryCatch(lmer(f_forced, data = data, REML = TRUE), error = function(e) NULL)
  if (is.null(m_forced)) return(NULL)
  an <- tryCatch(anova(m_forced, ddf = "Kenward-Roger"), error = function(e) NULL)
  if (is.null(an) || !("treatment" %in% rownames(an))) return(NULL)
  F_f <- an["treatment", "F value"]; df1_f <- an["treatment", "NumDF"]; df2_f <- an["treatment", "DenDF"]
  p_f <- an["treatment", "Pr(>F)"]
  r2 <- tryCatch(performance::r2_nakagawa(m_forced), error = function(e) list(R2_marginal = NA, R2_conditional = NA))
  icc_val <- tryCatch(performance::icc(m_forced)$ICC_adjusted, error = function(e) NA)

  trials <- unique(data[[trial_col]])
  loo <- vapply(trials, function(tr) {
    d2 <- data[data[[trial_col]] != tr, ]
    m2 <- tryCatch(lmer(formula(m_forced), data = d2, REML = TRUE), error = function(e) NULL)
    if (is.null(m2)) return(NA_real_)
    a2 <- tryCatch(anova(m2, ddf = "Kenward-Roger"), error = function(e) NULL)
    if (is.null(a2) || !("treatment" %in% rownames(a2))) return(NA_real_)
    a2["treatment", "Pr(>F)"]
  }, numeric(1))

  data.frame(
    outcome = label, singular = isSingular(m_forced),
    F_forced = F_f, df1_forced = df1_f, df2_forced = df2_f, p_forced = p_f,
    eta2p_forced = .eta2(F_f, df1_f, df2_f),
    R2m = r2$R2_marginal, R2c = r2$R2_conditional, ICC_trial = icc_val,
    loo_n_sig = sum(loo < 0.05, na.rm = TRUE), loo_n_total = sum(!is.na(loo)),
    loo_p_min = min(loo, na.rm = TRUE), loo_p_max = max(loo, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}
.pull_existing <- function(dir) {
  f <- file.path(STEP5, dir, "anova.csv")
  if (!file.exists(f)) return(c(F = NA, p = NA, eta2p = NA))
  a <- read.csv(f)
  r <- a[grepl("^treatment$", a$term), ]
  if (!nrow(r)) return(c(F = NA, p = NA, eta2p = NA))
  Fv <- r$chisq[1]; d1 <- r$df[1]; d2 <- r$df_denom[1]
  c(F = Fv, p = r$p_value[1], eta2p = .eta2(Fv, d1, d2))
}
.attach_existing <- function(out, map_dir) {
  existing <- lapply(out$outcome, function(o) .pull_existing(map_dir[[o]]))
  out$F_aicc_selected <- vapply(existing, function(x) unname(x["F"]), numeric(1))
  out$p_aicc_selected <- vapply(existing, function(x) unname(x["p"]), numeric(1))
  out$eta2p_aicc_selected <- vapply(existing, function(x) unname(x["eta2p"]), numeric(1))
  out$sig_flip <- (out$p_aicc_selected < 0.05) != (out$p_forced < 0.05)
  out
}

# =============================================================================
# (a) Behavioural interval-level (Level B, N=48, 3 sessions/trial)
# =============================================================================
dt <- fread(file.path(STEP5, "analysis_ready.csv"))
dt[, treatment := factor(treatment)]
dt[, fish_density_f := factor(fish_density)]
dt[, timepoint_f := factor(timepoint)]
dt[, trial := factor(phys_trial_id)]
dt[, logit_flow := qlogis(pmin(pmax(prop_time_in_flow, 1e-6), 1 - 1e-6))]

gd_path <- file.path(PROJECT_ROOT, "choice R pipeline/output/STEP2b_output/STEP2b_output_20260807_214103/group_dynamics_summary.csv")
gd <- fread(gd_path)
gd[, treatment := factor(treatment)]
gd[, fish_density_f := factor(fish_density)]
gd[, trial := factor(paste(tank, trial_date, treatment))]

beh_results <- list(
  .forced_trial_check(dt, "logit_flow", "fish_density_f + timepoint_f", "logit_flow (main-zone flow, interval)"),
  .forced_trial_check(dt, "switches_per_session", "fish_density_f + timepoint_f", "switches_per_session (interval)"),
  .forced_trial_check(gd, "mean_nnd_cm", "fish_density_f + timepoint_f", "mean_nnd_cm (interval)"),
  .forced_trial_check(gd, "mean_polarisation", "fish_density_f + timepoint_f", "mean_polarisation (interval)"),
  .forced_trial_check(gd, "mean_hull_area_cm2", "fish_density_f + timepoint_f", "mean_hull_area_cm2 (interval)"),
  .forced_trial_check(gd, "mean_iid_cm", "fish_density_f + timepoint_f", "mean_iid_cm (interval)"),
  .forced_trial_check(gd, "mean_centroid_spd_cm", "fish_density_f + timepoint_f", "mean_centroid_spd_cm (interval)")
)
beh_out <- do.call(rbind, beh_results[!vapply(beh_results, is.null, logical(1))])
map_dir_tp <- c(
  "logit_flow (main-zone flow, interval)"  = "zone_main_timepoint",
  "switches_per_session (interval)"        = "switches_timepoint",
  "mean_nnd_cm (interval)"                 = "nnd_timepoint",
  "mean_polarisation (interval)"           = "polarisation_timepoint",
  "mean_hull_area_cm2 (interval)"          = "hull_area_timepoint",
  "mean_iid_cm (interval)"                 = "iid_timepoint",
  "mean_centroid_spd_cm (interval)"        = "centroid_speed_timepoint"
)
beh_out <- .attach_existing(beh_out, map_dir_tp)
fwrite(beh_out, file.path(STEP5, "sensitivity_forced_trial_re_behavioural_interval.csv"))
cat("=== (a) Behavioural interval-level (N=48) forced-trial-RE check ===\n")
print(beh_out[, c("outcome","F_aicc_selected","p_aicc_selected","F_forced","p_forced","sig_flip","ICC_trial","loo_n_sig","loo_n_total")], digits = 4)

# =============================================================================
# (b) Cortisol (fish nested in trial)
# =============================================================================
cort <- fread(file.path(PROJECT_ROOT, "choice exp for claude mono and cortisol/claude output REVISED/b3_cortisol_clean.csv"))
cort[, condition := factor(condition)]
names(cort)[names(cort) == "condition"] <- "treatment"
cort[, trial := factor(trial)]
cort[, log_cort := log(plasma_cortisol)]
cort_res <- .forced_trial_check(cort[is.finite(log_cort) & !is.na(trial)], "log_cort", "", "plasma_cortisol (log)")
cort_out <- do.call(rbind, list(cort_res))
cort_out$F_aicc_selected <- 9.45613567858656   # from anova_kr_cortisol.csv (Data/)
cort_out$p_aicc_selected <- 0.00295857427243374
cort_out$eta2p_aicc_selected <- .eta2(9.45613567858656, 1, 73.1274511183651)
cort_out$sig_flip <- (cort_out$p_aicc_selected < 0.05) != (cort_out$p_forced < 0.05)
fwrite(cort_out, file.path(STEP5, "sensitivity_forced_trial_re_cortisol.csv"))
cat("\n=== (b) Cortisol forced-trial-RE check ===\n")
print(cort_out[, c("outcome","F_aicc_selected","p_aicc_selected","F_forced","p_forced","sig_flip","ICC_trial","loo_n_sig","loo_n_total")], digits = 4)

# =============================================================================
# (c) Monoamine headline cells (Dm 5-HT, DA, DOPAC)
# =============================================================================
mono <- fread(file.path(PROJECT_ROOT, "choice exp for claude mono and cortisol/claude output REVISED/Data/monoamine_long.csv"))
mono[, treatment := factor(treatment)]
mono[, trial := factor(trial)]
headline <- list(
  c(area = "Dm", analyte = "5-HT"),
  c(area = "Dm", analyte = "DA"),
  c(area = "Dm", analyte = "DOPAC")
)
mono_results <- lapply(headline, function(h) {
  sub <- mono[area == h["area"] & analyte == h["analyte"] & is.finite(value)]
  sub[, logval := log(value)]
  .forced_trial_check(sub, "logval", "", paste0(h["area"], " ", h["analyte"]))
})
mono_out <- do.call(rbind, mono_results[!vapply(mono_results, is.null, logical(1))])
# existing (AICc-selected) KR results, from anova_kr_monoamines.csv, treatment row
mono_existing <- data.frame(
  outcome = c("Dm 5-HT", "Dm DA", "Dm DOPAC"),
  F_aicc_selected = c(7.41041254845869, 2.51564506294918, 0.765778770769319),
  df1_aicc = c(1,1,1), df2_aicc = c(33.7974527103854, 34.6099918418756, 34.6928927364204),
  p_aicc_selected = c(0.0101784791170107, 0.121817125792, 0.387544434716305)
)
mono_out <- merge(mono_out, mono_existing, by = "outcome")
mono_out$eta2p_aicc_selected <- .eta2(mono_out$F_aicc_selected, mono_out$df1_aicc, mono_out$df2_aicc)
mono_out$sig_flip <- (mono_out$p_aicc_selected < 0.05) != (mono_out$p_forced < 0.05)
fwrite(mono_out, file.path(STEP5, "sensitivity_forced_trial_re_monoamines.csv"))
cat("\n=== (c) Monoamine headline cells (Dm) forced-trial-RE check ===\n")
print(mono_out[, c("outcome","F_aicc_selected","p_aicc_selected","F_forced","p_forced","sig_flip","ICC_trial","loo_n_sig","loo_n_total")], digits = 4)

cat("\nAll sensitivity CSVs written to:", STEP5, "\n")
