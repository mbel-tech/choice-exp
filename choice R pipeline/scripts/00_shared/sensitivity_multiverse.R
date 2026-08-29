# =============================================================================
# sensitivity_multiverse.R -- specification-curve analysis (memo 3.1)
# =============================================================================
# Primary flow-preference result (ALR(flow), trial-level) and four headline
# neurochemical results (Dm 5-HT, DA, DOPAC; plasma cortisol). For each, fit
# every defensible combination of the analytical forks identified in the
# memo and report the distribution of the treatment-effect p-value/estimate
# across specifications, plus the proportion significant.
#
# Grid for ALR(flow) (curated, not a literal full factorial -- some
# combinations are not meaningful together):
#   response scale : ALR(flow) logit  |  raw proportion (beta GLMM, logit link)
#   aggregation    : trial-level (mean of session proportions)  |
#                    interval-level with trial random intercept
#   random effects : (1|trial)  |  (1|trial)+(1|tank)  |  AICc-selected (existing)
#   transformation : none  |  log1p  |  sqrt   [applied to the trial-level
#                    continuous ALR response only; beta GLMM is already on
#                    its native bounded scale and logit-linked, so "none" only]
#
# Grid for each neurochemical headline outcome (Dm 5-HT/DA/DOPAC) and cortisol:
#   random effects : AICc-selected (existing, sample_id/plate only) |
#                    + (1|trial) added | (1|trial) alone (no sample_id/plate)
#   transformation : log  |  sqrt  |  Box-Cox (grid search, same as pipeline)
#   estimation     : REML/Kenward-Roger  |  parametric bootstrap (PBmodcomp)
# =============================================================================

suppressMessages({
  library(data.table); library(lme4); library(lmerTest); library(MASS); library(pbkrtest)
})
# Sum-to-zero contrasts so Type III tests of a main effect are marginal in
# the presence of an interaction, matching the project-wide convention
# (activity_analysis_STATS_choice_exp.R run_lmm_analysis(), line ~1095).
# Without this, a main-effect F-test in a model that also contains an
# interaction is NOT the same quantity the rest of the manuscript reports
# (the highest-order interaction term itself is contrast-invariant and is
# unaffected either way).
options(contrasts = c("contr.sum", "contr.poly"))
set.seed(20260706)
STEP5 <- file.path(PROJECT_ROOT, "choice R pipeline/output/STEP5_stats/STEP5_stats_20260808_123821")
N_BOOT <- 1000

.boxcox_lambda <- function(x) {
  x <- x[is.finite(x) & x > 0]
  bc <- MASS::boxcox(x ~ 1, plotit = FALSE, lambda = seq(-2, 2, 0.1))
  bc$x[which.max(bc$y)]
}
.bc_transform <- function(x, lam) if (abs(lam) < 1e-6) log(x) else (x^lam - 1) / lam

# =============================================================================
# PART A: ALR(flow) specification curve
# =============================================================================
dt <- fread(file.path(STEP5, "analysis_ready.csv"))
dt[, treatment := factor(treatment)][, fish_density_f := factor(fish_density)]
dt[, tank := factor(tank)][, trial := factor(phys_trial_id)][, timepoint_f := factor(timepoint)]
dt[, logit_flow := qlogis(pmin(pmax(prop_time_in_flow, 1e-6), 1 - 1e-6))]

agg <- dt[, .(prop_time_in_flow = mean(prop_time_in_flow), logit_flow = mean(logit_flow)),
          by = .(trial, treatment, fish_density_f, tank)]

specs <- list()
.add_spec <- function(desc, F_val, df1, df2, p_val) {
  specs[[length(specs) + 1]] <<- data.frame(spec = desc, F = F_val, df1 = df1, df2 = df2, p = p_val,
                                             sig = p_val < 0.05, stringsAsFactors = FALSE)
}
.try_lmm_p <- function(formula_str, data) {
  m <- tryCatch(lmer(as.formula(formula_str), data = data, REML = TRUE), error = function(e) NULL)
  if (is.null(m)) return(NULL)
  a <- tryCatch(anova(m, ddf = "Kenward-Roger"), error = function(e) NULL)
  if (is.null(a) || !("treatment" %in% rownames(a))) return(NULL)
  list(F = a["treatment", "F value"], df1 = a["treatment", "NumDF"], df2 = a["treatment", "DenDF"],
       p = a["treatment", "Pr(>F)"])
}

# -- trial-level, response=logit(ALR flow), RE x transform grid --
re_variants <- list("trial" = "(1|trial)", "trial_tank" = "(1|trial) + (1|tank)", "tank_only" = "(1|tank)")
trans_variants <- list(none = function(x) x, log1p = function(x) log1p(x - min(x) + 0.01),
                        sqrt = function(x) sqrt(x - min(x) + 0.01))
for (re_name in names(re_variants)) {
  for (tr_name in names(trans_variants)) {
    d2 <- copy(agg); d2[, y := trans_variants[[tr_name]](logit_flow)]
    r <- .try_lmm_p(paste("y ~ treatment + fish_density_f +", re_variants[[re_name]]), d2)
    if (!is.null(r)) .add_spec(paste0("trial-level, logit(ALR), RE=", re_name, ", trans=", tr_name),
                                r$F, r$df1, r$df2, r$p)
  }
}
# -- trial-level, response = raw proportion via beta-like transform (Smithson-Verkuilen squeeze) --
n <- nrow(agg); sv <- (agg$prop_time_in_flow * (n - 1) + 0.5) / n
agg[, prop_sv := sv]
for (re_name in names(re_variants)) {
  d2 <- copy(agg); d2[, y := qlogis(prop_sv)]
  r <- .try_lmm_p(paste("y ~ treatment + fish_density_f +", re_variants[[re_name]]), d2)
  if (!is.null(r)) .add_spec(paste0("trial-level, logit(raw prop, SV-squeezed), RE=", re_name), r$F, r$df1, r$df2, r$p)
}
# -- interval-level (N=48), with trial RE, categorical/continuous timepoint --
interval_re <- c(trial = "(1|trial)", trial_tank = "(1|trial) + (1|tank)")
for (fixed in c("treatment * timepoint_f + fish_density_f", "treatment * timepoint + fish_density_f")) {
  for (re_key in names(interval_re)) {
    r <- .try_lmm_p(paste("logit_flow ~", fixed, "+", interval_re[[re_key]]), dt)
    if (!is.null(r)) .add_spec(paste0("interval-level, ", ifelse(grepl("timepoint_f", fixed), "categorical", "continuous"),
                                       " interval, RE=", re_key), r$F, r$df1, r$df2, r$p)
  }
}
# -- existing AICc-selected primary (for reference row) --
a_exist <- read.csv(file.path(STEP5, "zone_flow_logit_aggregated/anova.csv"))
r_exist <- a_exist[a_exist$term == "treatment", ]
.add_spec("AICc-selected primary (as reported)", r_exist$chisq[1], r_exist$df[1], r_exist$df_denom[1], r_exist$p_value[1])

alr_curve <- do.call(rbind, specs)
alr_curve <- alr_curve[order(alr_curve$p), ]
fwrite(alr_curve, file.path(STEP5, "multiverse_alr_flow.csv"))
cat("=== ALR(flow) specification curve: N =", nrow(alr_curve), "specifications ===\n")
cat("Proportion significant (p<0.05):", round(mean(alr_curve$sig, na.rm = TRUE), 3), "\n")
cat("Median p:", round(median(alr_curve$p, na.rm = TRUE), 5), "  Range:", round(min(alr_curve$p, na.rm=TRUE),6), "-", round(max(alr_curve$p, na.rm=TRUE),4), "\n")
print(alr_curve, digits = 4)

# =============================================================================
# PART B: Neurochemical headline outcomes -- RE x transform x estimation grid
# =============================================================================
mono <- fread(file.path(PROJECT_ROOT, "choice exp for claude mono and cortisol/claude output REVISED/Data/monoamine_long.csv"))
mono[, treatment := factor(treatment)][, trial := factor(trial)][, sample_id := factor(sample_id)][, area := factor(area)]

.mono_spec_curve <- function(an, ar, boot = FALSE) {
  sub <- mono[analyte == an & area == ar & is.finite(value) & value > 0]
  out <- list()
  lam <- tryCatch(.boxcox_lambda(sub$value), error = function(e) NA)
  transforms <- list(log = log(sub$value), sqrt = sqrt(sub$value))
  if (is.finite(lam)) transforms$boxcox <- .bc_transform(sub$value, lam)
  re_cands <- c(sample_id = "(1|sample_id)", sample_id_trial = "(1|sample_id) + (1|trial)", trial_only = "(1|trial)")
  for (tr_name in names(transforms)) {
    sub2 <- copy(sub); sub2[, y := transforms[[tr_name]]]
    for (re_key in names(re_cands)) {
      re_form <- re_cands[[re_key]]
      m <- tryCatch(lmer(as.formula(paste("y ~ treatment +", re_form)), data = sub2, REML = TRUE), error = function(e) NULL)
      if (is.null(m)) next
      a <- tryCatch(anova(m, ddf = "Kenward-Roger"), error = function(e) NULL)
      if (is.null(a) || !("treatment" %in% rownames(a))) next
      p_kr <- a["treatment", "Pr(>F)"]
      p_boot <- NA
      if (boot) {
        m0 <- tryCatch(lme4::lmer(as.formula(paste("y ~ 1 +", re_form)), data = sub2, REML = FALSE), error = function(e) NULL)
        m1 <- tryCatch(lme4::lmer(as.formula(paste("y ~ treatment +", re_form)), data = sub2, REML = FALSE), error = function(e) NULL)
        if (!is.null(m0) && !is.null(m1)) {
          pb <- tryCatch(pbkrtest::PBmodcomp(m1, m0, nsim = N_BOOT, seed = 20260706), error = function(e) NULL)
          if (!is.null(pb)) p_boot <- pb$test["PBtest", "p.value"]
        }
      }
      out[[length(out) + 1]] <- data.frame(
        outcome = paste(an, ar), transform = tr_name, re = re_key,
        F = a["treatment", "F value"], p_KR = p_kr, p_bootstrap = p_boot, stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, out)
}

headline_cells <- list(c("5-HT", "DM"), c("DA", "DM"), c("DOPAC", "DM"))
mono_curves <- lapply(headline_cells, function(h) .mono_spec_curve(h[1], h[2], boot = TRUE))
mono_curve_tab <- do.call(rbind, mono_curves)
fwrite(mono_curve_tab, file.path(STEP5, "multiverse_neurochemical.csv"))
cat("\n=== Neurochemical headline specification curve (Dm 5-HT/DA/DOPAC), with parametric bootstrap ===\n")
for (h in headline_cells) {
  sub <- mono_curve_tab[mono_curve_tab$outcome == paste(h[1], h[2]), ]
  cat(sprintf("%s: %d specs, %.1f%% sig (KR), median p_KR=%.4g, median p_boot=%.4g\n",
              paste(h[1], h[2]), nrow(sub), 100 * mean(sub$p_KR < 0.05, na.rm = TRUE),
              median(sub$p_KR, na.rm = TRUE), median(sub$p_bootstrap, na.rm = TRUE)))
}
print(mono_curve_tab, digits = 4)

# =============================================================================
# PART C: Cortisol specification curve + bootstrap
# =============================================================================
cort <- fread(file.path(PROJECT_ROOT, "choice exp for claude mono and cortisol/claude output REVISED/b3_cortisol_clean.csv"))
names(cort)[names(cort) == "condition"] <- "treatment"
cort[, treatment := factor(treatment)][, trial := factor(trial)][, plate := factor(plate)][, tank := factor(tank)]
cort <- cort[is.finite(plasma_cortisol) & plasma_cortisol > 0]
lam_cort <- tryCatch(.boxcox_lambda(cort$plasma_cortisol), error = function(e) NA)
cort_transforms <- list(log = log(cort$plasma_cortisol), sqrt = sqrt(cort$plasma_cortisol))
if (is.finite(lam_cort)) cort_transforms$boxcox <- .bc_transform(cort$plasma_cortisol, lam_cort)
cort_re <- c("plate" = "(1|plate)", "plate_trial" = "(1|plate) + (1|trial)", "trial_only" = "(1|trial)", "tank_only" = "(1|tank)")
cort_out <- list()
for (tr_name in names(cort_transforms)) {
  cort2 <- copy(cort); cort2[, y := cort_transforms[[tr_name]]]
  for (re_name in names(cort_re)) {
    m <- tryCatch(lmer(as.formula(paste("y ~ treatment +", cort_re[[re_name]])), data = cort2, REML = TRUE), error = function(e) NULL)
    if (is.null(m)) next
    a <- tryCatch(anova(m, ddf = "Kenward-Roger"), error = function(e) NULL)
    if (is.null(a) || !("treatment" %in% rownames(a))) next
    m0 <- tryCatch(lme4::lmer(as.formula(paste("y ~ 1 +", cort_re[[re_name]])), data = cort2, REML = FALSE), error = function(e) NULL)
    m1 <- tryCatch(lme4::lmer(as.formula(paste("y ~ treatment +", cort_re[[re_name]])), data = cort2, REML = FALSE), error = function(e) NULL)
    p_boot <- NA
    if (!is.null(m0) && !is.null(m1)) {
      pb <- tryCatch(pbkrtest::PBmodcomp(m1, m0, nsim = N_BOOT, seed = 20260706), error = function(e) NULL)
      if (!is.null(pb)) p_boot <- pb$test["PBtest", "p.value"]
    }
    cort_out[[length(cort_out) + 1]] <- data.frame(
      transform = tr_name, re = re_name, F = a["treatment", "F value"],
      p_KR = a["treatment", "Pr(>F)"], p_bootstrap = p_boot, stringsAsFactors = FALSE
    )
  }
}
cort_curve <- do.call(rbind, cort_out)
fwrite(cort_curve, file.path(STEP5, "multiverse_cortisol.csv"))
cat("\n=== Cortisol specification curve, with parametric bootstrap ===\n")
cat(sprintf("%d specs, %.1f%% sig (KR), median p_KR=%.4g, median p_boot=%.4g\n",
            nrow(cort_curve), 100 * mean(cort_curve$p_KR < 0.05, na.rm = TRUE),
            median(cort_curve$p_KR, na.rm = TRUE), median(cort_curve$p_bootstrap, na.rm = TRUE)))
print(cort_curve, digits = 4)

cat("\n=== ALL MULTIVERSE + BOOTSTRAP ANALYSES COMPLETE ===\n")
