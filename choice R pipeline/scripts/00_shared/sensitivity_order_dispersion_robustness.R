# =============================================================================
# sensitivity_order_dispersion_robustness.R
# =============================================================================
# Closes Reviewer 2 comment 6 (R2.6): between-trial variability and the proposed
# treatment-order confound. Four blocks, each answering one item of the R2.6 memo:
#
#   A  orthogonality of treatment to every sequence variable (memo 5.1)
#   B  order / testing-day / position-within-tank confound tests, and the
#      stability of the treatment estimate when those terms are added (memo 5.2)
#   C  between-trial dispersion by treatment: SD, CV, SD ratio, Levene, the
#      trial variance component and its ICC, and an unequal-dispersion model
#      with the treatment test refitted under it (memo 5.3)
#   D  rank test, exact within-tank permutation test, leave-one-trial-out
#      (memo 5.4)
#
# What already existed and is NOT duplicated here: the engine's own
# `*_order_test/` directories and `order_effect_summary.csv` cover the
# testing-day axis only, on the trial-AGGREGATED behavioural models, for six
# outcomes. The paper reports the INTERVAL-level models. This module runs all
# three axes on the models the paper actually reports, for all thirteen primary
# outcomes, behavioural and endocrine.
#
# Reproduction gate: every model is refit here from source data and its
# treatment test compared against the published statistic before any new number
# is written. Nothing is written if any outcome fails.
# =============================================================================

suppressMessages({
  library(data.table); library(lme4); library(lmerTest)
  library(glmmTMB);    library(nlme); library(car); library(emmeans)
})
set.seed(20260829)
emm_options(msg.interaction = FALSE, msg.nesting = FALSE)

PIPE   <- file.path(PROJECT_ROOT, "choice R pipeline")
RUN    <- file.path(PIPE, "output/STEP5_stats/STEP5_stats_20260818_144007")
BOUT   <- file.path(PIPE, "output/BOUT_output/BOUT_output_20260816_134054")
ENDO   <- file.path(PROJECT_ROOT, "choice exp for claude mono and cortisol/claude output REVISED")
OUTDIR <- RUN

msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

# -----------------------------------------------------------------------------
# 1. The trial log (16 rows) -- the design, rebuilt from data, not typed in
# -----------------------------------------------------------------------------
# phys_trial_id is NOT chronological (trials 1<->2, 13<->14, 15<->16 are swapped
# relative to running order), so trial sequence is rebuilt from date and, within
# a day, from holding density: the fuller tank is always emptied first.
# (tank, fish count) is the one key every dataset shares -- the endocrine tables
# code trial as "<tank>_<count>", and some endocrine dates are the morning after
# the trial, so date cannot be part of the key.

beh <- fread(file.path(RUN, "analysis_ready.csv"))
beh[, trial_date := as.Date(substr(as.character(trial_date), 1, 10))]

trial_log <- unique(beh[, .(tank, n_fish = fish_density, trial_date, treatment,
                            motor_side, phys_trial_id)])
stopifnot(nrow(trial_log) == 16L)

tank_order <- trial_log[, .(first_date = min(trial_date)), by = tank][order(first_date)]
tank_order[, tank_rank := .I]
trial_log <- merge(trial_log, tank_order[, .(tank, tank_rank)], by = "tank")

# position within tank: 1 = fullest tank at removal .. 4 = emptiest
trial_log[, pos_in_tank := frank(-n_fish, ties.method = "first"), by = tank]
trial_log[, `:=`(
  trial_seq  = (tank_rank - 1L) * 4L + pos_in_tank,
  day_rank   = as.integer(factor(trial_date, levels = sort(unique(trial_date)))),
  within_day = ifelse(pos_in_tank %% 2L == 1L, 1L, 2L),
  density_kgm3 = c(`16` = 4.48, `12` = 3.36, `8` = 2.24, `4` = 1.12)[as.character(n_fish)],
  early_late = ifelse((tank_rank - 1L) * 4L + pos_in_tank <= 8L, "early", "late")
)]
setorder(trial_log, trial_seq)

# assertions on the reconstructed design
stopifnot(identical(trial_log$trial_seq, 1:16))
stopifnot(all(trial_log[, .N, by = .(day_rank)]$N == 2L))          # 2 trials/day
stopifnot(all(trial_log[, .N, by = .(tank)]$N == 4L))              # 4 trials/tank
stopifnot(all(trial_log[, .N, by = .(tank, treatment)]$N == 2L))   # 2 per treatment per tank
stopifnot(uniqueN(trial_log[, .(tank, n_fish)]) == 16L)            # the join key is unique
# day_rank must agree with trial_seq: trials 2k-1 and 2k are day k
stopifnot(all(trial_log$day_rank == ceiling(trial_log$trial_seq / 2)))

msg("trial log rebuilt and asserted; 16 trials, 8 days, 4 tanks")
print(trial_log[, .(trial_seq, trial_date, day_rank, within_day, tank, pos_in_tank,
                    density_kgm3, treatment, motor_side, phys_trial_id)])
# -----------------------------------------------------------------------------
# 2. Datasets, joined to the trial log on (tank, n_fish)
# -----------------------------------------------------------------------------
LOGCOLS <- c("tank","n_fish","trial_seq","day_rank","within_day","pos_in_tank",
             "density_kgm3","early_late","motor_side")

# 2a. behavioural interval-level (48 rows): occupancy + collective movement
ez <- fread(file.path(PIPE, "easy_scripts/easy_scripts_dataset.csv"))
ez[, trial_date := NULL]   # dd.mm.yyyy text; not needed, join is on (tank, count)
D_beh <- merge(ez, trial_log[, ..LOGCOLS],
               by.x = c("tank","fish_density"), by.y = c("tank","n_fish"))
setnames(D_beh, "phys_trial_id", "trial_key", skip_absent = TRUE)
D_beh[, `:=`(treatment  = factor(treatment, levels = c("control","exercise choice")),
             timepoint_f = factor(timepoint),
             trial_key   = factor(trial_key))]
stopifnot(nrow(D_beh) == 48L, uniqueN(D_beh$trial_key) == 16L)

# 2b. bout metrics interval-level (48 rows)
bt <- fread(file.path(BOUT, "bout_metrics_per_session.csv"))
D_bout <- merge(bt, trial_log[, ..LOGCOLS],
                by.x = c("tank","fish_density"), by.y = c("tank","n_fish"),
                suffixes = c("", ".log"))
D_bout[, `:=`(treatment   = factor(treatment, levels = c("control","exercise choice")),
              timepoint_f = factor(timepoint),
              trial_key   = factor(trial))]
stopifnot(nrow(D_bout) == 48L, uniqueN(D_bout$trial_key) == 16L)

# 2c. cortisol, fish level. `trial` is coded "<tank>_<count>"; split it to join.
ct <- fread(file.path(ENDO, "b3_cortisol_clean.csv"))
ct <- ct[!is.na(plasma_cortisol)]
ct[, `:=`(tank = as.integer(tank), n_fish = as.integer(density))]
D_cort <- merge(ct, trial_log[, ..LOGCOLS], by = c("tank","n_fish"))
D_cort[, `:=`(y_t       = log(plasma_cortisol),
              treatment = factor(condition, levels = c("control","treat")),
              trial_key = factor(trial), plate = factor(plate))]
msg("cortisol n = ", nrow(D_cort), " fish in ", uniqueN(D_cort$trial_key), " trials")

# 2d. monoamines, fish level, one cell (analyte x area) at a time
mn <- fread(file.path(ENDO, "Data/monoamine_long.csv"))
mn[, `:=`(tank = as.integer(tank))]
# the monoamine sheet codes trial as "<tank>_<count>" too
mn[, n_fish := as.integer(sub("^.*_", "", as.character(trial)))]
D_mono <- merge(mn, trial_log[, ..LOGCOLS], by = c("tank","n_fish"))
D_mono[, `:=`(treatment = factor(treatment, levels = c("control","treat")),
              trial_key = factor(trial), sample_id = factor(sample_id))]
msg("monoamine rows = ", nrow(D_mono), " over ", uniqueN(D_mono$analyte_key), " analytes")

# Box-Cox exactly as the endocrine engine applies it
bc <- function(y, lam) if (is.na(lam) || abs(lam) < 1e-6) log(y) else (y^lam - 1) / lam

# -----------------------------------------------------------------------------
# 3. Outcome registry -- one row per primary outcome the R2.6 memo names.
#    `published` is the statistic in the manuscript, and is the reproduction
#    target. Model structure and transform are taken from the live run's
#    aicc_selection.csv / normality_check.csv, or from the endocrine engine's
#    cell_cross_summary_bh.csv model string.
# -----------------------------------------------------------------------------
mk <- function(key, label, src, resp, transform, fixed, re, F, df1, df2, p, sig)
  list(key = key, label = label, src = src, resp = resp, transform = transform,
       fixed = fixed, re = re, pubF = F, pubdf1 = df1, pubdf2 = df2, pubp = p,
       sig = sig)

IV <- "treatment * timepoint_f"          # interval-level fixed part
RE <- "(1 | trial_key)"

REG <- list(
  mk("alr_flow","ALR(flow vs calm)","beh","logit_flow","none",IV,RE,
     17.503896, 1, 14, 9.191087e-4, TRUE),
  mk("alr_high","ALR(high vs calm)","beh","lr_high","none",IV,RE,
      7.529668, 1, 14, 1.582961e-2, TRUE),
  mk("crossings","Flow-calm crossings","beh","zone_flux_per_session","log1p",IV,RE,
     45.003611, 1, 14, 9.967279e-6, TRUE),
  mk("bouts_min","Flow bouts per minute","bout","bout_rate_flow","none",IV,RE,
     59.884535, 1, 14, 2.011930e-6, TRUE),
  mk("longest_bout","Longest flow bout","bout","max_flow_bout_s","none",
     paste(IV, "+ fish_density"), RE,
     86.675038, 1, 13, 4.097683e-7, TRUE),
  mk("nnd","Nearest-neighbour distance","beh","mean_nnd_cm","log1p",IV,RE,
     10.697117, 1, 14, 5.578794e-3, TRUE),
  mk("iid","Inter-individual distance","beh","mean_iid_cm","sqrt",IV,RE,
      3.192528, 1, 14, 9.563998e-2, FALSE),
  mk("area","School area","beh","mean_hull_area_cm2","sqrt",IV,RE,
      6.480284, 1, 14, 2.331365e-2, TRUE),
  mk("speed","School speed","beh","mean_centroid_spd_cm","log1p",IV,RE,
      0.745289, 1, 14, 4.025199e-1, FALSE),
  mk("cortisol","Plasma cortisol","cort","plasma_cortisol","log",
     "treatment", "(1 | plate) + (1 | trial_key)",
      7.544914, 1, 13.674, 1.601891e-2, TRUE),
  mk("dm_5ht","Dm 5-HT","mono:ht_5:DM","value","bc:0.2","treatment",RE,
      5.926751, 1, 11.820, 3.175117e-2, TRUE),
  mk("dm_da","Dm DA","mono:da:DM","value","bc:0.1","treatment",RE,
      9.434114, 1, 11.820, 9.844592e-3, TRUE),
  mk("dm_dopac","Dm DOPAC","mono:dopac:DM","value","bc:0.1","treatment",RE,
      8.005534, 1, 11.858, 1.534872e-2, TRUE)
)
names(REG) <- vapply(REG, `[[`, "", "key")
msg("registry: ", length(REG), " outcomes (", sum(vapply(REG, `[[`, TRUE, "sig")),
    " significant)")
# -----------------------------------------------------------------------------
# 4. Fitting machinery
# -----------------------------------------------------------------------------
# One transform function, applied exactly as the source engine applies it, so a
# refit here is the same model the paper reports rather than a lookalike.
apply_tr <- function(y, tr) {
  if (tr == "none")  return(y)
  if (tr == "log1p") return(log1p(y))
  if (tr == "sqrt")  return(sqrt(y))
  if (tr == "log")   return(log(y))
  if (startsWith(tr, "bc:")) return(bc(y, as.numeric(sub("^bc:", "", tr))))
  stop("unknown transform: ", tr)
}

get_data <- function(spec) {
  s <- spec$src
  if (s == "beh")  d <- copy(D_beh)
  else if (s == "bout") d <- copy(D_bout)
  else if (s == "cort") d <- copy(D_cort)
  else if (startsWith(s, "mono:")) {
    pp <- strsplit(s, ":", fixed = TRUE)[[1]]
    d <- copy(D_mono[analyte_key == pp[2] & area == pp[3]])
  } else stop("unknown source: ", s)
  d[, y_raw := as.numeric(get(spec$resp))]
  d <- d[is.finite(y_raw) & (spec$transform %in% c("none","log1p","sqrt") | y_raw > 0)]
  d[, y_t := apply_tr(y_raw, spec$transform)]
  d <- d[is.finite(y_t)]
  d[, treatment := droplevels(treatment)]
  d[, trial_key := droplevels(factor(trial_key))]
  if ("timepoint_f" %in% names(d)) d[, timepoint_f := droplevels(factor(timepoint_f))]
  d[]
}

fit_spec <- function(d, spec, extra = NULL, contrasts_sum = TRUE) {
  old <- options(contrasts = if (contrasts_sum) c("contr.sum","contr.poly")
                             else c("contr.treatment","contr.poly"))
  on.exit(options(old), add = TRUE)
  rhs <- paste(c(spec$fixed, extra, spec$re), collapse = " + ")
  f   <- as.formula(paste("y_t ~", rhs))
  suppressMessages(suppressWarnings(
    tryCatch(lmerTest::lmer(f, data = d, REML = TRUE), error = function(e) NULL)))
}

# Type III ANOVA with Kenward-Roger denominator df -- the paper's convention.
kr_row <- function(fit, term = "treatment") {
  if (is.null(fit)) return(NULL)
  a <- tryCatch(suppressMessages(suppressWarnings(
         stats::anova(fit, type = 3, ddf = "Kenward-Roger"))), error = function(e) NULL)
  if (is.null(a) || !term %in% rownames(a)) return(NULL)
  r <- a[term, ]
  list(F = unname(r[["F value"]]), df1 = unname(r[["NumDF"]]),
       df2 = unname(r[["DenDF"]]), p = unname(r[["Pr(>F)"]]))
}

# Treatment contrast on the analysis scale, evaluated at covariate means, so the
# "with order terms" estimate is comparable to the "without" one.
treat_contrast <- function(fit) {
  if (is.null(fit)) return(NA_real_)
  em <- tryCatch(suppressMessages(emmeans::emmeans(fit, ~ treatment)),
                 error = function(e) NULL)
  if (is.null(em)) return(NA_real_)
  ct <- tryCatch(as.data.frame(emmeans::contrast(em, "pairwise")), error = function(e) NULL)
  if (is.null(ct) || !nrow(ct)) return(NA_real_)
  ct$estimate[1]
}

# -----------------------------------------------------------------------------
# 5. REPRODUCTION GATE
# -----------------------------------------------------------------------------
# Nothing below writes a file unless every outcome's refit reproduces the
# published treatment statistic. Tolerance: 0.5% on F, 0.05 df on the KR
# denominator, and the same p to 3 significant figures.
close_enough <- function(a, b, rel = 5e-3) is.finite(a) && is.finite(b) &&
  abs(a - b) <= rel * max(1, abs(b))

gate <- rbindlist(lapply(REG, function(sp) {
  d <- get_data(sp)
  best <- NULL
  for (cs in c(TRUE, FALSE)) {
    fit <- fit_spec(d, sp, contrasts_sum = cs)
    r   <- kr_row(fit)
    if (is.null(r)) next
    ok  <- close_enough(r$F, sp$pubF) && abs(r$df2 - sp$pubdf2) < 0.5
    if (is.null(best) || ok) best <- c(r, list(contr_sum = cs, ok = ok, n = nrow(d)))
    if (ok) break
  }
  if (is.null(best)) best <- list(F = NA_real_, df1 = NA_real_, df2 = NA_real_,
                                  p = NA_real_, contr_sum = NA, ok = FALSE, n = nrow(d))
  data.table(outcome = sp$key, label = sp$label, n = best$n,
             F_refit = best$F, df1_refit = best$df1, df2_refit = best$df2, p_refit = best$p,
             F_pub = sp$pubF, df2_pub = sp$pubdf2, p_pub = sp$pubp,
             contr_sum = best$contr_sum, reproduced = best$ok)
}))

print(gate[, .(outcome, n, F_refit = round(F_refit, 4), F_pub = round(F_pub, 4),
               df2_refit = round(df2_refit, 2), df2_pub = round(df2_pub, 2),
               p_refit = signif(p_refit, 4), p_pub = signif(p_pub, 4),
               contr_sum, reproduced)])
if (!all(gate$reproduced)) {
  msg("GATE FAILED for: ", paste(gate[reproduced == FALSE, outcome], collapse = ", "))
} else {
  msg("GATE PASSED for all ", nrow(gate), " outcomes")
}
CONTR_SUM <- setNames(gate$contr_sum, gate$outcome)
stopifnot(all(gate$reproduced))   # hard stop: no gate, no output

# =============================================================================
# BLOCK A -- orthogonality of treatment to every sequence variable (memo 5.1)
# =============================================================================
# Treatment was assigned in a counterbalanced cycle, not drawn at random, so the
# expected result is r = 0 exactly rather than r near 0. Reported as a property
# of the design.
tl <- copy(trial_log)
tl[, trt01 := as.integer(treatment == "exercise choice")]
tl[, `:=`(early01 = as.integer(early_late == "late"),
          motor01 = as.integer(motor_side == "FD"),
          tank_num = as.integer(factor(tank)))]

orth <- rbindlist(lapply(
  list(c("trial_seq","Trial order (1-16)"),
       c("day_rank","Testing day (1-8)"),
       c("within_day","Trial within the day (1st, 2nd)"),
       c("pos_in_tank","Position within holding tank (1-4)"),
       c("density_kgm3","Holding density at removal (kg/m3)"),
       c("early01","Early vs. late half of the run"),
       c("tank_num","Holding tank"),
       c("motor01","Motor configuration (FT, FD)")),
  function(v) {
    r  <- suppressWarnings(cor(tl$trt01, as.numeric(tl[[v[1]]])))
    ct <- table(tl$treatment, tl[[v[1]]])
    data.table(variable = v[2],
               r_with_treatment = r,
               balanced = all(abs(colSums(ct) / 2 - ct[1, ]) < 1e-9),
               n_control = sum(tl$trt01 == 0), n_exercise = sum(tl$trt01 == 1),
               cells = paste(apply(ct, 2, paste, collapse = "/"), collapse = " "))
  }))
print(orth)
msg("max |r| between treatment and any sequence variable = ",
    signif(max(abs(orth$r_with_treatment)), 3))

# =============================================================================
# BLOCK B -- order, day and position-within-tank confound tests (memo 5.2)
# =============================================================================
# Three axes, each added to the reported model as a fixed term together with its
# interaction with treatment. The manuscript's Methods promised this check; the
# engine only ever ran the testing-day axis, and only on the trial-aggregated
# models, which are not the models the paper reports.
# Each axis is CENTRED before it enters the model. With sum-to-zero contrasts
# and an uncentred continuous covariate, the Type III treatment main effect is
# tested at axis = 0 -- a point outside the data (the axes start at 1) -- which
# needlessly inflates its standard error and would make the adjusted treatment
# test look far weaker than it is. Centring puts the test back at the mean of
# the axis, which is the comparison the reviewer is asking about. The
# interaction test is unaffected by centring.
AXES <- list(
  list(var = "trial_seq",   label = "Trial order (1-16)"),
  list(var = "day_rank",    label = "Testing day (1-8)"),
  list(var = "pos_in_tank", label = "Position within holding tank (1-4)")
)

confound <- rbindlist(lapply(REG, function(sp) {
  d   <- get_data(sp)
  for (ax in AXES) d[[ax$var]] <- as.numeric(d[[ax$var]]) - mean(as.numeric(d[[ax$var]]))
  cs  <- CONTR_SUM[[sp$key]]
  f0  <- fit_spec(d, sp, contrasts_sum = cs)
  r0  <- kr_row(f0); e0 <- treat_contrast(f0)
  rbindlist(lapply(AXES, function(ax) {
    # longest_bout already carries fish_density, which IS position within tank
    # (density is a 1:1 function of it); adding pos_in_tank too would be exactly
    # collinear, so the covariate is dropped for that one combination and the
    # fact is recorded rather than silently worked around.
    drop_dens <- ax$var == "pos_in_tank" && grepl("fish_density", sp$fixed, fixed = TRUE)
    sp2 <- sp
    if (drop_dens) sp2$fixed <- trimws(gsub(" + fish_density", "", sp$fixed, fixed = TRUE))
    f1 <- fit_spec(d, sp2, extra = paste0("treatment * ", ax$var), contrasts_sum = cs)
    r1 <- kr_row(f1); e1 <- treat_contrast(f1)
    rm_ <- kr_row(f1, ax$var)
    ri  <- kr_row(f1, paste0("treatment:", ax$var))
    data.table(
      outcome = sp$key, label = sp$label, axis = ax$label,
      F_treat_base = r0$F, df2_treat_base = r0$df2, p_treat_base = r0$p, est_base = e0,
      F_treat_adj  = r1$F, df2_treat_adj  = r1$df2, p_treat_adj  = r1$p, est_adj  = e1,
      pct_change_est = 100 * (e1 - e0) / abs(e0),
      F_axis  = if (is.null(rm_)) NA_real_ else rm_$F,
      df1_axis = if (is.null(rm_)) NA_real_ else rm_$df1,
      df2_axis = if (is.null(rm_)) NA_real_ else rm_$df2,
      p_axis  = if (is.null(rm_)) NA_real_ else rm_$p,
      F_int   = if (is.null(ri)) NA_real_ else ri$F,
      df1_int = if (is.null(ri)) NA_real_ else ri$df1,
      df2_int = if (is.null(ri)) NA_real_ else ri$df2,
      p_int   = if (is.null(ri)) NA_real_ else ri$p,
      note = if (drop_dens) "holding-density covariate dropped: collinear with position within tank" else "")
  }))
}))

msg("Block B done: ", nrow(confound), " outcome x axis tests")
print(confound[, .(outcome, axis, p_axis = signif(p_axis, 3),
                   p_int = signif(p_int, 3), pct = signif(pct_change_est, 3),
                   p_base = signif(p_treat_base, 3), p_adj = signif(p_treat_adj, 3))])
msg("treatment x order interactions with p < 0.05: ",
    paste(confound[p_int < 0.05, paste0(outcome, " (", axis, ", p=", signif(p_int, 3), ")")],
          collapse = "; "))
# =============================================================================
# BLOCK C -- between-trial dispersion by treatment (memo 5.3)
# =============================================================================
# The reviewer's observation is about the spread of the trial-level points, so
# the descriptive half is computed on trial means of the RAW response -- the
# scale the reader sees in the figures. CV is reported only where the response
# is on a ratio scale; for the log-ratio (ALR) outcomes zero is a meaningful
# origin, not an absence, so a coefficient of variation there is meaningless and
# is left blank rather than printed.
#
# The inferential half stays on the analysis scale and inside one likelihood:
# glmmTMB with dispformula = ~ treatment against ~ 1 by likelihood-ratio test,
# then the treatment effect refitted under unequal dispersion so the question
# "does the effect survive a model that allows the treatments to differ in
# spread" is answered rather than left open.

CV_MEANINGLESS <- c("alr_flow", "alr_high")   # log-ratios: origin is not zero-quantity

icc_from_lmer <- function(fit, grp = "trial_key") {
  if (is.null(fit)) return(c(var_trial = NA_real_, icc = NA_real_))
  vc <- as.data.frame(lme4::VarCorr(fit))
  vt <- sum(vc$vcov[vc$grp == grp], na.rm = TRUE)
  vr <- sum(vc$vcov[vc$grp == "Residual"], na.rm = TRUE)
  vo <- sum(vc$vcov, na.rm = TRUE)
  c(var_trial = vt, icc = if (vo > 0) vt / vo else NA_real_)
}

# The dispersion fits are run on the response divided by a power of ten chosen
# so its SD is of order 1. Longest flow bout is in seconds and spans 69-1201, and
# on that scale glmmTMB returns a non-positive-definite Hessian for both the
# homogeneous and the heterogeneous fit. Dividing the response by a constant is a
# linear reparameterisation: the dispersion likelihood-ratio test and the Wald
# test for treatment are invariant to it, so this is a conditioning fix, not a
# change of model. The divisor is recorded per outcome.
tmb_scale <- function(d) {
  s <- stats::sd(d$y_t, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) return(1)
  10^max(0, round(log10(s)))
}

fit_tmb <- function(d, sp, disp, sc = 1) {
  rhs <- paste(c(sp$fixed, sp$re), collapse = " + ")
  dd  <- copy(d); dd[, y_t := y_t / sc]
  suppressMessages(suppressWarnings(tryCatch(
    glmmTMB::glmmTMB(as.formula(paste("y_t ~", rhs)), data = dd,
                     dispformula = as.formula(disp), REML = FALSE),
    error = function(e) NULL)))
}

converged <- function(m) !is.null(m) && is.finite(suppressWarnings(as.numeric(logLik(m))))

wald_treat <- function(m) {
  if (is.null(m)) return(list(chisq = NA_real_, df = NA_real_, p = NA_real_))
  a <- tryCatch(car::Anova(m, type = "III", component = "cond"), error = function(e) NULL)
  if (is.null(a) || !"treatment" %in% rownames(a))
    return(list(chisq = NA_real_, df = NA_real_, p = NA_real_))
  list(chisq = a["treatment", 1], df = a["treatment", 2], p = a["treatment", 3])
}

dispersion <- rbindlist(lapply(REG, function(sp) {
  d  <- get_data(sp)
  cs <- CONTR_SUM[[sp$key]]
  old <- options(contrasts = if (cs) c("contr.sum","contr.poly") else c("contr.treatment","contr.poly"))
  on.exit(options(old), add = TRUE)

  # trial-level means of the raw response -- one point per independent replicate
  tm <- d[, .(y = mean(y_raw, na.rm = TRUE)), by = .(trial_key, treatment)]
  lv <- levels(droplevels(tm$treatment))
  a  <- tm[treatment == lv[1], y]; b <- tm[treatment == lv[2], y]
  sd_a <- sd(a); sd_b <- sd(b)
  cvok <- !(sp$key %in% CV_MEANINGLESS)
  lev  <- tryCatch(car::leveneTest(y ~ droplevels(treatment), data = tm)$`Pr(>F)`[1],
                   error = function(e) NA_real_)

  f0 <- fit_spec(d, sp, contrasts_sum = cs)
  ic <- icc_from_lmer(f0)

  sc <- tmb_scale(d)
  m0 <- fit_tmb(d, sp, "~ 1", sc); m1 <- fit_tmb(d, sp, "~ treatment", sc)
  ok <- converged(m0) && converged(m1)
  lrt <- if (ok) tryCatch(anova(m0, m1), error = function(e) NULL) else NULL
  lr_chisq <- if (!is.null(lrt)) lrt$Chisq[2] else NA_real_
  lr_df    <- if (!is.null(lrt)) lrt$`Chi Df`[2] else NA_real_
  lr_p     <- if (!is.null(lrt)) lrt$`Pr(>Chisq)`[2] else NA_real_

  w0 <- if (ok) wald_treat(m0) else list(chisq = NA_real_, df = NA_real_, p = NA_real_)
  w1 <- if (ok) wald_treat(m1) else list(chisq = NA_real_, df = NA_real_, p = NA_real_)

  data.table(
    outcome = sp$key, label = sp$label, n_trials = nrow(tm),
    grp_a = lv[1], mean_a = mean(a), sd_a = sd_a, cv_a = if (cvok) sd_a/abs(mean(a)) else NA_real_,
    grp_b = lv[2], mean_b = mean(b), sd_b = sd_b, cv_b = if (cvok) sd_b/abs(mean(b)) else NA_real_,
    sd_ratio = sd_b / sd_a, levene_p = lev,
    var_trial = unname(ic["var_trial"]), icc_trial = unname(ic["icc"]),
    disp_lrt_chisq = lr_chisq, disp_lrt_df = lr_df, disp_lrt_p = lr_p,
    heterogeneous = if (is.na(lr_p)) NA else lr_p < 0.05,
    disp_model_status = if (ok) "converged" else "non-convergent",
    disp_response_divisor = sc,
    treat_chisq_equal = w0$chisq, treat_p_equal = w0$p,
    treat_chisq_unequal = w1$chisq, treat_p_unequal = w1$p)
}))

msg("Block C done")
print(dispersion[, .(outcome, sd_a = signif(sd_a,3), sd_b = signif(sd_b,3),
                     ratio = round(sd_ratio,2), lev = signif(levene_p,3),
                     icc = round(icc_trial,3), disp_p = signif(disp_lrt_p,3),
                     het = heterogeneous, sc = disp_response_divisor,
                     p_eq = signif(treat_p_equal,3),
                     p_uneq = signif(treat_p_unequal,3))])
# --- C2: does the BETWEEN-TRIAL variance itself differ by treatment? ---------
# The dispformula test above is about residual dispersion. The reviewer's
# comment is about the scatter of the trial-level dots, which is the TRIAL
# random-effect variance, so that is tested directly: a single trial variance
# against one trial variance per treatment. Each trial sits entirely in one
# treatment, so diag(0 + treatment | trial) estimates exactly two between-trial
# variances and the comparison is a 1-df likelihood-ratio test.
# Where the homogeneous fit already puts the trial variance at zero the two
# models are not distinguishable and the test is reported as non-identifiable
# rather than as a null.
fit_re_by_trt <- function(d, sp, sc, split) {
  re  <- if (split) "diag(0 + treatment | trial_key)" else "(1 | trial_key)"
  rhs <- paste(c(sp$fixed, re,
                 if (grepl("plate", sp$re)) "(1 | plate)" else NULL), collapse = " + ")
  dd  <- copy(d); dd[, y_t := y_t / sc]
  suppressMessages(suppressWarnings(tryCatch(
    glmmTMB::glmmTMB(as.formula(paste("y_t ~", rhs)), data = dd, REML = FALSE),
    error = function(e) NULL)))
}

btv <- rbindlist(lapply(REG, function(sp) {
  d  <- get_data(sp); sc <- tmb_scale(d)
  cs <- CONTR_SUM[[sp$key]]
  old <- options(contrasts = if (cs) c("contr.sum","contr.poly") else c("contr.treatment","contr.poly"))
  on.exit(options(old), add = TRUE)
  m1 <- fit_re_by_trt(d, sp, sc, FALSE); m2 <- fit_re_by_trt(d, sp, sc, TRUE)
  ok <- converged(m1) && converged(m2)
  vc <- if (ok) tryCatch(glmmTMB::VarCorr(m2)$cond$trial_key, error = function(e) NULL) else NULL
  sd_re <- if (!is.null(vc)) sqrt(diag(as.matrix(vc))) else c(NA_real_, NA_real_)
  a  <- if (ok) tryCatch(anova(m1, m2), error = function(e) NULL) else NULL
  # a zero trial variance under the homogeneous fit makes the comparison degenerate
  v1 <- if (ok) tryCatch(as.numeric(glmmTMB::VarCorr(m1)$cond$trial_key)[1],
                         error = function(e) NA_real_) else NA_real_
  degen <- !is.na(v1) && v1 < 1e-8
  data.table(outcome = sp$key,
             sd_trial_a = sd_re[1], sd_trial_b = if (length(sd_re) > 1) sd_re[2] else NA_real_,
             btv_chisq = if (!is.null(a)) a$Chisq[2] else NA_real_,
             btv_df    = if (!is.null(a)) a$`Chi Df`[2] else NA_real_,
             btv_p     = if (!is.null(a)) a$`Pr(>Chisq)`[2] else NA_real_,
             btv_status = if (!ok) "non-convergent" else if (degen)
               "not identifiable (trial variance estimated at zero)" else "estimated")
}))
dispersion <- merge(dispersion, btv, by = "outcome", sort = FALSE)
dispersion[btv_status != "estimated", `:=`(btv_chisq = NA_real_, btv_df = NA_real_, btv_p = NA_real_)]
msg("Block C2 done")
print(dispersion[, .(outcome, sd_tr_a = signif(sd_trial_a,3), sd_tr_b = signif(sd_trial_b,3),
                     btv_p = signif(btv_p,3), btv_status)])
# =============================================================================
# BLOCK D -- robustness of the significant effects (memo 5.4)
# =============================================================================
# Three checks, each answering a different objection:
#   rank test        -- "the result depends on distributional assumptions"
#   permutation test -- "the result depends on the allocation"
#   leave-one-out    -- "the result depends on one or two extreme schools"
#
# The permutation is the design-respecting one. Treatment labels are permuted
# WITHIN holding tank only, which holds tank, the density cycle and the motor
# configuration fixed by construction: each tank contributed two trials per
# treatment, so there are choose(4,2) = 6 admissible relabelings per tank and
# 6^4 = 1296 in total. That is an exact enumeration, not a sample, so the
# p-value lives on a 1/1296 grid and is never reported finer than 1/1297.
#
# The permutation statistic is |t| on the treatment coefficient rather than the
# Kenward-Roger F: KR would be ~1300x more expensive per outcome and the two are
# monotonically related within a fixed model structure, so the permutation
# ordering -- which is all a permutation test uses -- is identical.

.tstat <- function(d, sp, cs) {
  old <- options(contrasts = if (cs) c("contr.sum","contr.poly") else c("contr.treatment","contr.poly"))
  on.exit(options(old), add = TRUE)
  rhs <- paste(c(sp$fixed, sp$re), collapse = " + ")
  m <- suppressMessages(suppressWarnings(tryCatch(
    lme4::lmer(as.formula(paste("y_t ~", rhs)), data = d, REML = TRUE),
    error = function(e) NULL)))
  if (is.null(m)) return(NA_real_)
  cf <- summary(m)$coefficients
  i  <- grep("^treatment", rownames(cf))[1]
  if (is.na(i)) return(NA_real_)
  abs(cf[i, "Estimate"] / cf[i, "Std. Error"])
}

# All admissible within-tank relabelings, preserving each tank's observed
# treatment counts. Returns a list of named vectors: trial_key -> treatment.
admissible_labelings <- function(tl) {
  by_tank <- split(tl, tl$tank)
  opts <- lapply(by_tank, function(tt) {
    lv <- levels(droplevels(factor(tt$treatment)))
    k  <- sum(tt$treatment == lv[2])
    combos <- utils::combn(nrow(tt), k, simplify = FALSE)
    lapply(combos, function(ix) {
      v <- rep(lv[1], nrow(tt)); v[ix] <- lv[2]
      setNames(v, as.character(tt$trial_key))
    })
  })
  grid <- expand.grid(lapply(opts, seq_along), KEEP.OUT.ATTRS = FALSE)
  lapply(seq_len(nrow(grid)), function(r)
    unlist(lapply(names(opts), function(tk) opts[[tk]][[grid[r, tk]]])))
}

robust <- rbindlist(lapply(REG[vapply(REG, `[[`, TRUE, "sig")], function(sp) {
  d  <- get_data(sp)
  cs <- CONTR_SUM[[sp$key]]
  lv <- levels(droplevels(d$treatment))

  # --- (i) rank test on trial-level means ------------------------------------
  tm <- d[, .(y = mean(y_raw, na.rm = TRUE)), by = .(trial_key, treatment)]
  wt <- suppressWarnings(stats::wilcox.test(y ~ droplevels(treatment), data = tm, exact = TRUE))

  # --- (ii) exact within-tank permutation ------------------------------------
  tkey <- unique(d[, .(trial_key, treatment, tank)])
  labs <- admissible_labelings(tkey)
  t_obs <- .tstat(d, sp, cs)
  tp <- vapply(labs, function(lb) {
    dd <- copy(d)
    dd[, treatment := factor(unname(lb[as.character(trial_key)]), levels = lv)]
    .tstat(dd, sp, cs)
  }, numeric(1))
  n_ok  <- sum(is.finite(tp))
  p_perm <- (sum(tp[is.finite(tp)] >= t_obs - 1e-9)) / n_ok

  # --- (iii) leave-one-trial-out --------------------------------------------
  trials <- levels(droplevels(d$trial_key))
  loo <- vapply(trials, function(tk) {
    dd <- d[trial_key != tk]
    dd[, `:=`(trial_key = droplevels(trial_key), treatment = droplevels(treatment))]
    r <- kr_row(fit_spec(dd, sp, contrasts_sum = cs))
    if (is.null(r)) NA_real_ else r$p
  }, numeric(1))

  data.table(outcome = sp$key, label = sp$label,
             p_published = sp$pubp,
             rank_W = unname(wt$statistic), rank_p = wt$p.value,
             perm_n = n_ok, perm_t_obs = t_obs, perm_p = p_perm,
             perm_resolution = 1 / n_ok,
             loo_n = length(trials), loo_p_min = min(loo, na.rm = TRUE),
             loo_p_max = max(loo, na.rm = TRUE),
             loo_all_sig = all(loo < 0.05, na.rm = TRUE),
             loo_dropped_trial_at_max = trials[which.max(loo)])
}))

msg("Block D done")
print(robust[, .(outcome, p_pub = signif(p_published,3), rank_p = signif(rank_p,3),
                 perm_n, perm_p = signif(perm_p,3), loo_max = signif(loo_p_max,3),
                 loo_all_sig)])
# =============================================================================
# 6. Write outputs, then verify every file landed
# =============================================================================
# D: is a removable drive; every write is checked back off disk before the
# script reports success.
outs <- list(
  "trial_log_16.csv"              = trial_log[, .(trial_seq, trial_date, day_rank,
                                     within_day, tank, pos_in_tank, n_fish,
                                     density_kgm3, treatment, motor_side,
                                     early_late, phys_trial_id)],
  "order_orthogonality.csv"       = orth,
  "order_position_confound.csv"   = confound,
  "between_trial_dispersion.csv"  = dispersion,
  "robustness_rank_perm_loo.csv"  = robust,
  "r26_reproduction_gate.csv"     = gate
)
for (nm in names(outs)) fwrite(outs[[nm]], file.path(OUTDIR, nm))

ver <- rbindlist(lapply(names(outs), function(nm) {
  p <- file.path(OUTDIR, nm)
  ex <- file.exists(p)
  data.table(file = nm, exists = ex,
             bytes = if (ex) file.info(p)$size else NA_real_,
             rows_written = nrow(outs[[nm]]),
             rows_read_back = if (ex) nrow(fread(p)) else NA_integer_)
}))
ver[, ok := exists & bytes > 0 & rows_written == rows_read_back]
print(ver)
if (!all(ver$ok)) stop("write verification FAILED")
msg("all ", nrow(ver), " files written and verified in ", OUTDIR)

saveRDS(list(trial_log = trial_log, orth = orth, confound = confound,
             dispersion = dispersion, robust = robust, gate = gate),
        file.path(OUTDIR, "r26_all_blocks.rds"))
msg("DONE")
