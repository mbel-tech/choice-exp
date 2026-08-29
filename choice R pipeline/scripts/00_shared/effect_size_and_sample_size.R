# =============================================================================
# effect_size_and_sample_size.R
# =============================================================================
# Implements `sample_size_recommendation.md` against the CURRENT analysis
# outputs, and computes the project-wide ANOVA effect sizes.
#
# WHY THIS FILE EXISTS
#   The recommendation document was written against an earlier analysis pass;
#   its worked table (its section 3) quotes F values that the M1/M2 revision has
#   since superseded (e.g. ALR(flow) F = 6.27 -> 7.53; Dm DA 8.32 -> 9.72; Dm
#   DOPAC 5.81 -> 7.06; Dm 5-HT 5.51 -> 6.11). This script therefore applies the
#   document's METHOD to the CURRENT numbers rather than copying its arithmetic.
#
# WHAT IT PRODUCES  ->  choice R pipeline/output/POWER_output/
#   realised_effect_sizes.csv      per-term eta^2_p, omega^2_p, Cohen's d + CI
#   design_as_run_MDE.csv          minimum detectable effect of the design as run
#   sample_size_planning_table.csv trials-per-arm needed, by d x ICC
#   icc_estimates.csv              ICCs extracted from the retained mixed models
#   sample_size_design.png         the planning curves
#
# ---------------------------------------------------------------------------
# EFFECT-SIZE ESTIMATORS (and why these)
#
#   partial eta-squared      eta^2_p = (F * df1) / (F * df1 + df2)
#       Proportion of variance attributable to a term once the model's other
#       terms are held constant. This is the estimator already used elsewhere in
#       this project, it is defined for every F-test regardless of whether the
#       fit is `lm`, `lmer` (Satterthwaite) or `lmer` (Kenward-Roger), and it
#       needs only the term's own F and degrees of freedom -- so it can be
#       recomputed from an ANOVA table without refitting. PARTIAL, not classical:
#       the denominator excludes variance explained by the other fixed terms.
#
#   partial omega-squared    omega^2_p = df1(F - 1) / (df1(F - 1) + df2 + 1)
#       The small-sample BIAS-CORRECTED counterpart. eta^2_p is upward-biased,
#       badly so at the small df2 that a 16-trial design produces, so omega^2_p
#       is reported alongside it and is the value that should be carried into
#       any future power calculation. Clamped at 0 (it can go negative for
#       near-null effects, which is not interpretable as variance).
#
#   Cohen's d                d = 2*sqrt(F) / sqrt(df2)      [two-group contrast]
#       Used ONLY for the sample-size work, because the cluster-randomised
#       formulation in the recommendation document is parameterised in d. Note
#       this project's F-tests come from mixed models whose df2 is already the
#       CLUSTER-level (trial) denominator for the behavioural outcomes, so the
#       resulting d is a cluster-level standardised effect for those rows and an
#       individual-level one for the fish-level endocrine rows; the two are
#       flagged in the `unit` column and must not be pooled.
#
#   REPEATED-MEASURES FACTORS
#       For the trial x interval sequence models the interval and treatment x
#       interval terms are within-school. eta^2_p / omega^2_p as defined above
#       are computed per term from that term's own F and its own Satterthwaite
#       df2, so the within-school error term is used automatically and no
#       generalised eta-squared conversion is applied. Generalised eta-squared
#       would be preferable for comparison ACROSS studies with different designs;
#       it is not used here because the project reports partial eta-squared
#       throughout and mixing the two would make the manuscript inconsistent.
#
#   LIMITATIONS
#       eta^2_p from a mixed model is not a variance decomposition of the total
#       outcome variance (the random-effect variance is excluded), so it is not
#       comparable with an eta^2_p from a fixed-effects-only ANOVA on the same
#       data. Confidence intervals on d use the noncentral t and assume the
#       two-group contrast is the only fixed effect of interest.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  has_effectsize <- requireNamespace("effectsize", quietly = TRUE)
})
ts_msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), " | ", ..., "\n", sep = "")

PIPE  <- file.path(PROJECT_ROOT, "choice R pipeline")
STEP5 <- file.path(PIPE, "STEP5_stats", "STEP5_stats_20260511_175522")
ENDO  <- file.path(PROJECT_ROOT, "choice exp for claude mono and cortisol/claude output REVISED/Data")
OUT   <- file.path(PIPE, "output", "POWER_output")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ---- estimators -------------------------------------------------------------
eta2_p   <- function(F, df1, df2) ifelse(is.finite(F) & is.finite(df1) & is.finite(df2),
                                         (F * df1) / (F * df1 + df2), NA_real_)
omega2_p <- function(F, df1, df2) ifelse(is.finite(F) & is.finite(df1) & is.finite(df2),
                                         pmax(0, (df1 * (F - 1)) / (df1 * (F - 1) + df2 + 1)),
                                         NA_real_)
# Two-group contrast: d = 2*sqrt(F)/sqrt(df2). Equivalent to effectsize::F_to_d.
F_to_d   <- function(F, df2) ifelse(is.finite(F) & is.finite(df2) & df1_is_one(F),
                                    2 * sqrt(F) / sqrt(df2), NA_real_)
df1_is_one <- function(x) TRUE  # guarded at call sites instead

# noncentral-t CI on d for a two-group contrast (df2 = denominator df)
d_ci <- function(F, df2, conf = 0.95) {
  if (!is.finite(F) || !is.finite(df2) || df2 <= 0) return(c(NA_real_, NA_real_))
  tval <- sqrt(F)
  lo_p <- (1 - conf) / 2; hi_p <- 1 - lo_p
  f <- function(ncp, target) stats::pt(tval, df2, ncp) - target
  safe <- function(target) {
    out <- tryCatch(stats::uniroot(f, interval = c(-50, 50), target = target,
                                   extendInt = "yes")$root, error = function(e) NA_real_)
    out
  }
  ncp_lo <- safe(hi_p); ncp_hi <- safe(lo_p)
  c(2 * ncp_lo / sqrt(df2), 2 * ncp_hi / sqrt(df2))
}

# =============================================================================
# 1) REALISED EFFECT SIZES from the current analysis outputs
# =============================================================================
rows <- list()
.add <- function(domain, outcome, term, F, df1, df2, p, unit) {
  ci <- if (isTRUE(df1 == 1)) d_ci(F, df2) else c(NA_real_, NA_real_)
  rows[[length(rows) + 1L]] <<- data.table(
    domain = domain, outcome = outcome, term = term,
    F = F, df1 = df1, df2 = df2, p = p,
    eta2_p = eta2_p(F, df1, df2), omega2_p = omega2_p(F, df1, df2),
    cohens_d = if (isTRUE(df1 == 1)) 2 * sqrt(F) / sqrt(df2) else NA_real_,
    d_lo = ci[1], d_hi = ci[2], unit = unit)
}

# --- behaviour: trial-level aggregated models (df2 is already cluster-level) --
beh <- c(zone_flow_logit_aggregated = "ALR(flow) zone preference",
         nnd_aggregated             = "Nearest-neighbour distance",
         polarisation_aggregated    = "Polarisation",
         centroid_speed_aggregated  = "School speed",
         switches_aggregated        = "Flow-calm transitions",
         hull_area_aggregated       = "School area",
         iid_aggregated             = "Inter-individual distance")
for (d in names(beh)) {
  f <- file.path(STEP5, d, "anova.csv")
  if (!file.exists(f)) next
  a <- data.table::fread(f)
  r <- a[term == "treatment"]
  if (!nrow(r)) next
  .add("behaviour", beh[[d]], "treatment",
       r$chisq[1], r$df[1], r$df_denom[1], r$p_value[1], "trial (cluster)")
}

# --- endocrine: cortisol + the 24 monoamine cells (fish-level) ----------------
fc <- file.path(ENDO, "anova_kr_cortisol.csv")
if (file.exists(fc)) {
  a <- data.table::fread(fc)
  .add("endocrine", "Plasma cortisol (log)", "condition",
       a$F[1], a$Df[1], a$Df.res[1], a[["Pr(>F)"]][1], "fish")
}
fm <- file.path(ENDO, "cell_cross_summary_bh.csv")
if (file.exists(fm)) {
  a <- data.table::fread(fm)
  for (i in seq_len(nrow(a)))
    .add("endocrine", sprintf("%s %s", a$area[i], a$analyte[i]), "treatment",
         a$F_stat[i], a$df1[i], a$df2[i], a$p_raw[i], "fish")
}

# --- behavioural sequence: both analysis levels -------------------------------
seqdirs <- list.dirs(file.path(PIPE, "output", "SEQ_output"),
                     recursive = FALSE, full.names = TRUE)
if (length(seqdirs)) {
  sq <- seqdirs[which.max(file.mtime(seqdirs))]
  for (fn in c(trial = "seq_anova_trial_level.csv",
               tp    = "seq_anova_trial_x_timepoint.csv")) {
    f <- file.path(sq, fn)
    if (!file.exists(f)) next
    a <- data.table::fread(f)
    for (i in seq_len(nrow(a)))
      .add(paste0("sequence (", a$level[i], ")"),
           sprintf("%s %s", a$alphabet[i], a$metric[i]), a$term[i],
           a$F[i], a$df1[i], a$df2[i], a$p[i],
           if (a$level[i] == "trial") "trial (cluster)" else "trial x interval")
  }
}

realised <- data.table::rbindlist(rows, fill = TRUE)
data.table::fwrite(realised, file.path(OUT, "realised_effect_sizes.csv"))
ts_msg("realised_effect_sizes.csv: ", nrow(realised), " terms")

# =============================================================================
# 2) ICC from the retained mixed models  (recommendation section 2.2)
# =============================================================================
# The document asks for the empirical ICC prior to be taken from the fitted
# models rather than assumed. Cortisol is the model whose clustering matters
# most (the stressor is delivered per trial), so it is the one refit here.
icc_rows <- list()
if (requireNamespace("performance", quietly = TRUE) &&
    requireNamespace("lme4", quietly = TRUE)) {
  # b3_cortisol_clean.csv sits one level ABOVE Data/ in the engine's output tree.
  fcort <- Filter(file.exists, c(file.path(ENDO, "b3_cortisol_clean.csv"),
                                 file.path(dirname(ENDO), "b3_cortisol_clean.csv")))
  fcort <- if (length(fcort)) fcort[1] else NA_character_
  if (!is.na(fcort)) {
    cd <- data.table::fread(fcort)
    nm <- names(cd)
    ycol <- if ("plasma_cortisol" %in% nm) "plasma_cortisol" else NA
    if (!is.na(ycol) && "trial" %in% nm) {
      cd <- cd[is.finite(get(ycol)) & get(ycol) > 0]
      cd[, y := log(get(ycol))]
      fit <- tryCatch(suppressMessages(suppressWarnings(
        lme4::lmer(y ~ 1 + (1 | trial), data = cd))), error = function(e) NULL)
      if (!is.null(fit)) {
        ic <- tryCatch(performance::icc(fit), error = function(e) NULL)
        if (!is.null(ic))
          icc_rows[[length(icc_rows) + 1L]] <- data.table(
            model = "log(plasma cortisol) ~ 1 + (1|trial)",
            icc_adjusted = as.numeric(ic$ICC_adjusted),
            icc_unadjusted = as.numeric(ic$ICC_unadjusted),
            n_obs = nrow(cd), n_trial = data.table::uniqueN(cd$trial))
      }
    }
  }
}
icc_tab <- if (length(icc_rows)) data.table::rbindlist(icc_rows, fill = TRUE) else
  data.table(model = character(), icc_adjusted = numeric(),
             icc_unadjusted = numeric(), n_obs = integer(), n_trial = integer())
data.table::fwrite(icc_tab, file.path(OUT, "icc_estimates.csv"))
if (nrow(icc_tab)) ts_msg("Empirical ICC (cortisol, fish within trial): ",
                          sprintf("%.3f", icc_tab$icc_adjusted[1]))

# =============================================================================
# 3) CLUSTER-RANDOMISED POWER  (recommendation section 2.1)
# =============================================================================
#   DE     = 1 + (m - 1) * ICC
#   lambda = d / sqrt(2 * DE / (m * k))
#   df     = 2 * (k - 1)
#   power  = P(|t_nc(df, lambda)| > t_crit)
crt_power <- function(d, m, k, icc, alpha = 0.05) {
  DE  <- 1 + (m - 1) * icc
  df  <- 2 * (k - 1)
  if (df < 1) return(NA_real_)
  lam <- d / sqrt(2 * DE / (m * k))
  tc  <- stats::qt(1 - alpha / 2, df)
  stats::pt(-tc, df, lam) + (1 - stats::pt(tc, df, lam))
}
# smallest k (trials/arm) reaching the target power
crt_k_needed <- function(d, m, icc, power = 0.80, alpha = 0.05, kmax = 200) {
  for (k in 2:kmax) if (isTRUE(crt_power(d, m, k, icc, alpha) >= power)) return(k)
  NA_integer_
}
# minimum detectable effect of a design as run
crt_mde <- function(m, k, icc, power = 0.80, alpha = 0.05) {
  f <- function(d) crt_power(d, m, k, icc, alpha) - power
  tryCatch(stats::uniroot(f, c(0.01, 10))$root, error = function(e) NA_real_)
}

# --- design as run: 8 trials/arm x 5 fish, across the plausible ICC range -----
mde <- data.table::CJ(m = 5L, k = 8L, icc = c(0, 0.10, 0.20, 0.30))
mde[, design_effect := 1 + (m - 1) * icc]
mde[, n_eff_per_arm := k * m / design_effect]
mde[, MDE_cohens_d := mapply(crt_mde, m, k, icc)]
# the monoamine subsample (~18 fish/arm, i.e. ~3.6 fish/trial over 5 trials)
mono <- data.table::CJ(m = 5L, k = 4L, icc = c(0.10, 0.30))
mono[, design_effect := 1 + (m - 1) * icc]
mono[, n_eff_per_arm := k * m / design_effect]
mono[, MDE_cohens_d := mapply(crt_mde, m, k, icc)]
mde_all <- rbind(cbind(scenario = "Behaviour/cortisol: 8 trials/arm x 5 fish", mde),
                 cbind(scenario = "Monoamine subsample: ~4 trials/arm x 5 fish", mono))
data.table::fwrite(mde_all, file.path(OUT, "design_as_run_MDE.csv"))

# --- planning table: trials per arm needed, by target d x ICC x group size ----
plan <- data.table::CJ(d = c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0),
                       icc = c(0.10, 0.20, 0.30),
                       m = c(5L, 8L, 12L))
plan[, k_per_arm := mapply(crt_k_needed, d, m, icc)]
plan[, fish_per_arm := k_per_arm * m]
plan[, power_at_k := mapply(function(d, m, k, i)
  if (is.na(k)) NA_real_ else crt_power(d, m, k, i), d, m, k_per_arm, icc)]
data.table::fwrite(plan, file.path(OUT, "sample_size_planning_table.csv"))

# --- headline recommendation (document section 5.2): d = 0.7, ICC = 0.20 ------
rec_k5  <- crt_k_needed(0.7, 5, 0.20)
rec_k8  <- crt_k_needed(0.7, 8, 0.20)
ts_msg(sprintf("Recommendation @ d=0.7, ICC=0.20: %d trials/arm at 5 fish/trial (%d fish/arm); %d at 8 fish/trial (%d fish/arm)",
               rec_k5, rec_k5 * 5, rec_k8, rec_k8 * 8))

# =============================================================================
# 4) FIGURE — power curves
# =============================================================================
curve <- data.table::CJ(k = 4:24, m = c(5L, 8L, 12L), icc = c(0.10, 0.20, 0.30))
curve[, power := mapply(function(k, m, i) crt_power(0.7, m, k, i), k, m, icc)]
curve[, `Fish per trial` := factor(m)]
.fmt3 <- function(x) {
  s <- sprintf("%.3f", x)
  ifelse(grepl("0$", s), substr(s, 1, nchar(s) - 1), s)
}
curve[, icc_lab := factor(sprintf("ICC = %s", .fmt3(icc)))]

p_pow <- ggplot(curve, aes(k, power, colour = `Fish per trial`)) +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey40") +
  geom_line(linewidth = 1) +
  geom_point(size = 1.6) +
  facet_wrap(~ icc_lab) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  scale_colour_manual(values = c("5" = "#E69F00", "8" = "#0072B2", "12" = "#009E73")) +
  labs(title = "Cluster-randomised power for the exercise-choice design",
       subtitle = paste0("Target effect d = 0.70 (lower-bound-informed); treatment ",
                         "allocated at the level of the arena trial.\nDashed line = 80% power. ",
                         "Trials, not fish per trial, are the effective lever."),
       x = "Trials per treatment arm (k)", y = "Power") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        plot.title = element_text(face = "bold"))
ggsave(file.path(OUT, "sample_size_design.png"), p_pow,
       width = 250, height = 120, units = "mm", dpi = 300, bg = "white")

# =============================================================================
ts_msg("==================== POWER / EFFECT SIZE DONE ====================")
cat("\n--- Design as run: minimum detectable effect (80% power, alpha .05) ---\n")
print(mde_all[, .(scenario, m, k, icc, design_effect,
                  n_eff_per_arm = round(n_eff_per_arm, 1),
                  MDE_d = round(MDE_cohens_d, 2))])
cat("\n--- Trials per arm needed (target power 80%) ---\n")
print(plan[m == 5L, .(d, icc, m, k_per_arm, fish_per_arm,
                      power = round(power_at_k, 3))])
cat("\n--- Largest realised effects (by omega^2_p) ---\n")
print(head(realised[order(-omega2_p),
                    .(domain, outcome, term, F = round(F, 2), df1,
                      df2 = round(df2, 1), p = round(p, 4),
                      eta2_p = round(eta2_p, 3), omega2_p = round(omega2_p, 3),
                      d = round(cohens_d, 2))], 12))
cat("\nOutputs written to: ", OUT, "\n")
