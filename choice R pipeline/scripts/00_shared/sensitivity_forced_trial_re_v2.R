# =============================================================================
# sensitivity_forced_trial_re_v2.R
# =============================================================================
# Corrected version: for every outcome where AICc dropped trial entirely from
# the winning random-effects structure, ADD (1|trial) to that SAME winning
# structure (never replace it, never change the fixed-effects formula), and
# compare. This directly operationalises the memo's recommendation ("retain
# trial unconditionally") as an apples-to-apples sensitivity check rather than
# a different model.
#
# Confirmed AICc winners that DROPPED trial (from re_selection_aicc_*.csv /
# aicc_selection.csv, read directly, 2026-08-08):
#   - plasma_cortisol:        (1|plate)                  [fixed: treatment]
#   - 5-HT, 5-HIAA, DA, DOPAC,
#     5-HIAA/5-HT, DOPAC/DA:  (1|sample_id)               [fixed: treatment*area]
#   - polarisation_timepoint: (1|trial_date)              [fixed: treatment*timepoint + fish_density]
#   - switches_timepoint:     (1|tank)                    [fixed: treatment*timepoint_f]
# =============================================================================

suppressMessages({
  library(data.table); library(lme4); library(lmerTest); library(performance)
})
# Sum-to-zero contrasts so a main-effect Type III test in the presence of an
# interaction (treatment:area, treatment:timepoint) is the marginal effect
# the rest of the manuscript reports, matching run_lmm_analysis()'s
# convention. The interaction terms themselves are contrast-invariant.
options(contrasts = c("contr.sum", "contr.poly"))
set.seed(20260706)
STEP5 <- file.path(PROJECT_ROOT, "choice R pipeline/output/STEP5_stats/STEP5_stats_20260808_123821")
.eta2 <- function(F, df1, df2) (F * df1) / (F * df1 + df2)

.compare <- function(m_orig, m_added, term, label) {
  a1 <- tryCatch(anova(m_orig, ddf = "Kenward-Roger"), error = function(e) NULL)
  a2 <- tryCatch(anova(m_added, ddf = "Kenward-Roger"), error = function(e) NULL)
  r1 <- if (!is.null(a1) && term %in% rownames(a1)) a1[term, ] else NULL
  r2 <- if (!is.null(a2) && term %in% rownames(a2)) a2[term, ] else NULL
  data.frame(
    outcome = label, term = term,
    F_orig = if (!is.null(r1)) r1[["F value"]] else NA, p_orig = if (!is.null(r1)) r1[["Pr(>F)"]] else NA,
    F_trial_added = if (!is.null(r2)) r2[["F value"]] else NA, p_trial_added = if (!is.null(r2)) r2[["Pr(>F)"]] else NA,
    sig_flip = (if(!is.null(r1)) r1[["Pr(>F)"]]<0.05 else NA) != (if(!is.null(r2)) r2[["Pr(>F)"]]<0.05 else NA),
    ICC_trial_added = tryCatch(performance::icc(m_added)$ICC_adjusted, error = function(e) NA),
    singular_added = isSingular(m_added),
    stringsAsFactors = FALSE
  )
}

out <- list()

# ---- (1) Cortisol: (1|plate) -> (1|plate) + (1|trial) ----------------------
cort <- fread(file.path(PROJECT_ROOT, "choice exp for claude mono and cortisol/claude output REVISED/b3_cortisol_clean.csv"))
names(cort)[names(cort) == "condition"] <- "treatment"
cort[, treatment := factor(treatment)][, trial := factor(trial)][, plate := factor(plate)]
cort[, log_cort := log(plasma_cortisol)]
cort <- cort[is.finite(log_cort) & !is.na(trial) & !is.na(plate)]
m1 <- lmer(log_cort ~ treatment + (1|plate), data = cort, REML = TRUE)
m2 <- lmer(log_cort ~ treatment + (1|plate) + (1|trial), data = cort, REML = TRUE)
out$cort <- .compare(m1, m2, "treatment", "plasma_cortisol (log)")

# ---- (2) Monoamines: (1|sample_id) -> (1|sample_id) + (1|trial) -----------
mono <- fread(file.path(PROJECT_ROOT, "choice exp for claude mono and cortisol/claude output REVISED/Data/monoamine_long.csv"))
mono[, treatment := factor(treatment)][, trial := factor(trial)][, sample_id := factor(sample_id)]
mono[, area := factor(area)]
for (an in c("5-HT", "5-HIAA", "DA", "DOPAC", "5-HIAA/5-HT", "DOPAC/DA")) {
  sub <- mono[analyte == an & is.finite(value) & !is.na(trial)]
  sub[, logval := log(value)]
  m1 <- tryCatch(lmer(logval ~ treatment * area + (1|sample_id), data = sub, REML = TRUE), error = function(e) NULL)
  m2 <- tryCatch(lmer(logval ~ treatment * area + (1|sample_id) + (1|trial), data = sub, REML = TRUE), error = function(e) NULL)
  if (!is.null(m1) && !is.null(m2)) {
    out[[paste0("mono_", an, "_treat")]] <- .compare(m1, m2, "treatment", paste0(an, " (treatment)"))
    out[[paste0("mono_", an, "_int")]]   <- .compare(m1, m2, "treatment:area", paste0(an, " (treatment:area)"))
  }
}

# ---- (3) Behavioural interval-level: polarisation, switches ----------------
gd_path <- file.path(PROJECT_ROOT, "choice R pipeline/output/STEP2b_output/STEP2b_output_20260807_214103/group_dynamics_summary.csv")
gd <- fread(gd_path)
gd[, treatment := factor(treatment)][, fish_density_f := factor(fish_density)]
gd[, trial := factor(paste(tank, trial_date, treatment))][, trial_date := factor(trial_date)]
m1 <- lmer(mean_polarisation ~ treatment * timepoint + fish_density_f + (1|trial_date), data = gd, REML = TRUE)
m2 <- lmer(mean_polarisation ~ treatment * timepoint + fish_density_f + (1|trial_date) + (1|trial), data = gd, REML = TRUE)
out$polar_treat <- .compare(m1, m2, "treatment", "mean_polarisation (treatment, interval)")
out$polar_int   <- .compare(m1, m2, "treatment:timepoint", "mean_polarisation (treatment:timepoint, interval)")

dt <- fread(file.path(STEP5, "analysis_ready.csv"))
dt[, treatment := factor(treatment)][, timepoint_f := factor(timepoint)][, tank := factor(tank)]
m1 <- lmer(switches_per_session ~ treatment * timepoint_f + (1|tank), data = dt, REML = TRUE)
m2 <- lmer(switches_per_session ~ treatment * timepoint_f + (1|tank) + (1|phys_trial_id), data = dt, REML = TRUE)
out$switch_treat <- .compare(m1, m2, "treatment", "switches_per_session (treatment, interval)")
out$switch_int   <- .compare(m1, m2, "treatment:timepoint_f", "switches_per_session (treatment:timepoint_f, interval)")

final <- do.call(rbind, out[!vapply(out, is.null, logical(1))])
fwrite(final, file.path(STEP5, "sensitivity_trial_added_v2.csv"))
cat("=== Trial-added sensitivity check (AICc-dropped-trial cases only) ===\n")
print(final, digits = 4)
