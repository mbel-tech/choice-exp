# =============================================================================
# sensitivity_part2_loo_icc_equiv.R
# =============================================================================
# Completes three memo items on the headline outcomes:
#  (1) fixes the 5-HT/5-HIAA trial-added check (excludes the single non-finite
#      log(0) row per analyte; log() is otherwise the pipeline's own choice
#      for these two analytes' pooled-region model -- see decision log);
#  (2) leave-one-trial-out influence analysis (memo 3.5);
#  (3) R2m/R2c/ICC(trial) reporting (memo 3.3);
#  (4) TOST equivalence testing for the non-significant brain regions
#      (Vv, Vd, POA) against a pre-specified SESOI (memo 3.6). SESOI is set
#      to d = 0.70, the same threshold the manuscript's own sample-size
#      justification (2.8.5) already uses as its planning effect size --
#      re-using a value the manuscript has already committed to rather than
#      inventing a new one.
# =============================================================================

suppressMessages({
  library(data.table); library(lme4); library(lmerTest); library(performance)
})
options(contrasts = c("contr.sum", "contr.poly"))
# TOSTER is not installed in this environment; TOST (two one-sided tests) is
# implemented by hand below on log-transformed values with a Welch
# (unequal-variance) df, following Lakens et al. (2018) eq. 1-2 directly --
# no functional difference from the package for the simple two-sample case
# used here.
.manual_tost <- function(x1, x2, low_d, high_d) {
  n1 <- length(x1); n2 <- length(x2)
  m1 <- mean(x1); m2 <- mean(x2); v1 <- var(x1); v2 <- var(x2)
  sp <- sqrt(((n1 - 1) * v1 + (n2 - 1) * v2) / (n1 + n2 - 2))  # pooled SD for d
  diff <- m1 - m2
  se <- sqrt(v1 / n1 + v2 / n2)
  df <- (v1 / n1 + v2 / n2)^2 / ((v1 / n1)^2 / (n1 - 1) + (v2 / n2)^2 / (n2 - 1))
  low_raw <- low_d * sp; high_raw <- high_d * sp
  t1 <- (diff - low_raw) / se;  p1 <- 1 - pt(t1, df)   # H0: diff <= low_raw
  t2 <- (diff - high_raw) / se; p2 <- pt(t2, df)        # H0: diff >= high_raw
  list(d_obs = diff / sp, p_tost = max(p1, p2), df = df)
}
set.seed(20260706)
STEP5 <- file.path(PROJECT_ROOT, "choice R pipeline/output/STEP5_stats/STEP5_stats_20260808_123821")
OUTDIR <- STEP5
.eta2 <- function(F, df1, df2) (F * df1) / (F * df1 + df2)

# =============================================================================
# (1) 5-HT / 5-HIAA trial-added check, excluding the single non-finite row
# =============================================================================
mono <- fread(file.path(PROJECT_ROOT, "choice exp for claude mono and cortisol/claude output REVISED/Data/monoamine_long.csv"))
mono[, treatment := factor(treatment)][, trial := factor(trial)][, sample_id := factor(sample_id)]
mono[, area := factor(area)]

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
    stringsAsFactors = FALSE
  )
}

fixed_out <- list()
for (an in c("5-HT", "5-HIAA")) {
  sub <- mono[analyte == an & is.finite(value) & value > 0 & !is.na(trial)]
  sub[, logval := log(value)]
  m1 <- lmer(logval ~ treatment * area + (1|sample_id), data = sub, REML = TRUE)
  m2 <- lmer(logval ~ treatment * area + (1|sample_id) + (1|trial), data = sub, REML = TRUE)
  fixed_out[[paste0(an, "_treat")]] <- .compare(m1, m2, "treatment", paste0(an, " (treatment)"))
  fixed_out[[paste0(an, "_int")]]   <- .compare(m1, m2, "treatment:area", paste0(an, " (treatment:area)"))
}
fixed_tab <- do.call(rbind, fixed_out)
cat("=== (1) 5-HT / 5-HIAA trial-added check (fixed) ===\n"); print(fixed_tab, digits = 4)
fwrite(fixed_tab, file.path(OUTDIR, "sensitivity_5ht_5hiaa_trial_added.csv"))

# =============================================================================
# (2) + (3) Leave-one-trial-out and R2/ICC for headline outcomes
# =============================================================================
# Headline behavioural (trial-level, N=16) -- LOO drops one trial at a time.
dt <- fread(file.path(STEP5, "analysis_ready.csv"))
dt[, treatment := factor(treatment)][, fish_density_f := factor(fish_density)][, tank := factor(tank)]
dt[, trial := factor(phys_trial_id)]
dt[, logit_flow := qlogis(pmin(pmax(prop_time_in_flow, 1e-6), 1 - 1e-6))]
agg <- dt[, .(logit_flow = mean(logit_flow), switches_per_session = mean(switches_per_session)),
          by = .(trial, treatment, tank)]

gd_path <- file.path(PROJECT_ROOT, "choice R pipeline/output/STEP2b_output/STEP2b_output_20260807_214103/group_dynamics_summary.csv")
gd <- fread(gd_path)
gd[, treatment := factor(treatment)][, tank := factor(tank)]
gd[, trial := factor(paste(tank, trial_date, treatment))]
agg_gd <- gd[, .(mean_nnd_cm = mean(mean_nnd_cm), mean_polarisation = mean(mean_polarisation),
                  mean_hull_area_cm2 = mean(mean_hull_area_cm2), mean_iid_cm = mean(mean_iid_cm),
                  mean_centroid_spd_cm = mean(mean_centroid_spd_cm)),
              by = .(trial, treatment, tank)]

.loo_trial_level <- function(data, response, label) {
  f <- as.formula(paste(response, "~ treatment + (1|tank)"))
  m_full <- tryCatch(lmer(f, data = data, REML = TRUE), error = function(e) NULL)
  if (is.null(m_full)) return(NULL)
  units <- unique(data$trial)
  loo_p <- vapply(units, function(u) {
    d2 <- data[data$trial != u, ]
    m2 <- tryCatch(lmer(f, data = d2, REML = TRUE), error = function(e) NULL)
    if (is.null(m2)) return(NA_real_)
    a2 <- tryCatch(anova(m2, ddf = "Kenward-Roger"), error = function(e) NULL)
    if (is.null(a2) || !("treatment" %in% rownames(a2))) return(NA_real_)
    a2["treatment", "Pr(>F)"]
  }, numeric(1))
  a_full <- tryCatch(anova(m_full, ddf = "Kenward-Roger"), error = function(e) NULL)
  p_full <- if (!is.null(a_full) && "treatment" %in% rownames(a_full)) a_full["treatment", "Pr(>F)"] else NA
  r2 <- tryCatch(performance::r2_nakagawa(m_full), error = function(e) list(R2_marginal = NA, R2_conditional = NA))
  icc_val <- tryCatch(performance::icc(m_full)$ICC_adjusted, error = function(e) NA)
  data.frame(outcome = label, p_full = p_full, R2m = r2$R2_marginal, R2c = r2$R2_conditional,
             ICC_tank = icc_val, loo_n = length(units), loo_n_sig = sum(loo_p < 0.05, na.rm = TRUE),
             loo_p_min = min(loo_p, na.rm = TRUE), loo_p_max = max(loo_p, na.rm = TRUE),
             loo_flip = (p_full < 0.05) != all(loo_p < 0.05, na.rm = TRUE) && (p_full<0.05 || any(loo_p<0.05,na.rm=TRUE)),
             stringsAsFactors = FALSE)
}
loo_list <- list(
  .loo_trial_level(agg, "logit_flow", "logit_flow (ALR flow)"),
  .loo_trial_level(agg, "switches_per_session", "switches_per_session"),
  .loo_trial_level(agg_gd, "mean_nnd_cm", "mean_nnd_cm"),
  .loo_trial_level(agg_gd, "mean_polarisation", "mean_polarisation"),
  .loo_trial_level(agg_gd, "mean_hull_area_cm2", "mean_hull_area_cm2"),
  .loo_trial_level(agg_gd, "mean_iid_cm", "mean_iid_cm"),
  .loo_trial_level(agg_gd, "mean_centroid_spd_cm", "mean_centroid_spd_cm")
)
loo_tab <- do.call(rbind, loo_list[!vapply(loo_list, is.null, logical(1))])
cat("\n=== (2)+(3) Leave-one-trial-out + R2/ICC, behavioural trial-level (N=16) ===\n")
print(loo_tab, digits = 4)
fwrite(loo_tab, file.path(OUTDIR, "sensitivity_loo_r2_icc_behavioural.csv"))

# Cortisol LOO (fish-level, drop one trial's fish at a time)
cort <- fread(file.path(PROJECT_ROOT, "choice exp for claude mono and cortisol/claude output REVISED/b3_cortisol_clean.csv"))
names(cort)[names(cort) == "condition"] <- "treatment"
cort[, treatment := factor(treatment)][, trial := factor(trial)][, plate := factor(plate)]
cort[, log_cort := log(plasma_cortisol)]
cort <- cort[is.finite(log_cort) & !is.na(trial) & !is.na(plate)]
f_cort <- log_cort ~ treatment + (1|plate) + (1|trial)
m_cort <- lmer(f_cort, data = cort, REML = TRUE)
trials_c <- unique(cort$trial)
loo_cort_p <- vapply(trials_c, function(tr) {
  d2 <- cort[cort$trial != tr, ]
  m2 <- tryCatch(lmer(f_cort, data = d2, REML = TRUE), error = function(e) NULL)
  if (is.null(m2)) return(NA_real_)
  a2 <- tryCatch(anova(m2, ddf = "Kenward-Roger"), error = function(e) NULL)
  if (is.null(a2) || !("treatment" %in% rownames(a2))) return(NA_real_)
  a2["treatment", "Pr(>F)"]
}, numeric(1))
r2_cort <- tryCatch(performance::r2_nakagawa(m_cort), error = function(e) list(R2_marginal = NA, R2_conditional = NA))
icc_cort <- tryCatch(performance::icc(m_cort)$ICC_adjusted, error = function(e) NA)
cort_loo_tab <- data.frame(
  outcome = "plasma_cortisol (log, trial+plate RE)", R2m = r2_cort$R2_marginal, R2c = r2_cort$R2_conditional,
  ICC_trial = icc_cort, loo_n = length(trials_c), loo_n_sig = sum(loo_cort_p < 0.05, na.rm = TRUE),
  loo_p_min = min(loo_cort_p, na.rm = TRUE), loo_p_max = max(loo_cort_p, na.rm = TRUE)
)
cat("\n=== Cortisol LOO + R2/ICC (trial+plate RE model) ===\n")
print(cort_loo_tab, digits = 4)
fwrite(cort_loo_tab, file.path(OUTDIR, "sensitivity_loo_r2_icc_cortisol.csv"))

# =============================================================================
# (4) TOST equivalence testing, null brain regions (Vv, Vd, POA), SESOI d=0.70
# =============================================================================
SESOI <- 0.70
equiv_out <- list()
for (an in c("5-HT", "5-HIAA", "DA", "DOPAC", "5-HIAA/5-HT", "DOPAC/DA")) {
  for (ar in c("VV", "VD", "POA")) {
    sub <- mono[analyte == an & area == ar & is.finite(value) & value > 0]
    ctrl <- sub[treatment == "control", value]
    exer <- sub[treatment == "treat", value]
    if (length(ctrl) < 3 || length(exer) < 3) next
    tt <- tryCatch(.manual_tost(log(exer), log(ctrl), -SESOI, SESOI), error = function(e) NULL)
    if (is.null(tt)) next
    equiv_out[[paste(an, ar)]] <- data.frame(
      analyte = an, area = ar, n_ctrl = length(ctrl), n_exer = length(exer),
      d_obs = tt$d_obs, p_tost = tt$p_tost, equivalent_at_d070 = tt$p_tost < 0.05,
      stringsAsFactors = FALSE
    )
  }
}
equiv_tab <- do.call(rbind, equiv_out)
cat("\n=== (4) TOST equivalence testing, null brain regions (SESOI d=0.70) ===\n")
print(equiv_tab, digits = 3)
fwrite(equiv_tab, file.path(OUTDIR, "sensitivity_equivalence_null_regions.csv"))

cat("\nAll part-2 sensitivity CSVs written to:", OUTDIR, "\n")
