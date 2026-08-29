# =============================================================================
# BEHAVIOURAL SEQUENCE ANALYSIS — CHOICE EXPERIMENT (identity-free, group-level)
#
# WHY THIS SCRIPT
#   idtracker.ai cannot guarantee persistent fish identities within a session
#   (positions are correct, identities get reassigned; see STEP2 header). That
#   makes classic per-individual sequence analysis (one fish's behaviour states
#   over time) impossible. Instead we analyse the SCHOOL as the unit: the
#   temporal sequence of *collective* exercise-choice states within each
#   trial x timepoint.
#
#   The "behaviour" is spatial: what fraction of the school is in the flow
#   (exercise) zone vs the calm zone at each moment. We discretise that into
#   states and study the sequence's DYNAMICS — transitions, dwell times,
#   predictability — and whether they differ between the "exercise choice"
#   and "control" treatments.
#
#   TWO STATE ALPHABETS ARE RUN IN PARALLEL:
#     (1) GRADED  : prop_flow -> Low / Med / High  (data-driven terciles)
#     (2) BINARY  : prop_flow -> Flow-dominant / Calm-dominant  (majority > 0.5)
#
# INPUT
#   Latest output/STEP1_output/STEP1_output_*/master_fish_by_frame.csv
#   Columns used: trial_id, frame, time, timepoint, treatment, fish_density,
#                 tank, main_zone   (main_zone in {flow, calm})
#
# OUTPUT  ->  output/SEQ_output/SEQ_output_<timestamp>/
#   seq_metrics_per_trial_graded.csv     per (trial,timepoint) metrics, alphabet 1
#   seq_metrics_per_trial_binary.csv     per (trial,timepoint) metrics, alphabet 2
#   transition_matrix_pooled_graded.csv  treatment-pooled transition probs
#   transition_matrix_pooled_binary.csv  treatment-pooled transition probs
#   treatment_contrasts.csv              Wilcoxon + LMM per metric, both alphabets
#   figures/  (state ribbons, transition heatmaps, dwell boxplots, metric compare)
#
# STANDALONE: run with
#   "C:/Program Files/R/R-4.6.0/bin/Rscript.exe" --vanilla <this file>
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(stringr)
  library(ggtext)
  has_lme4 <- requireNamespace("lme4", quietly = TRUE)
  has_lmerTest <- requireNamespace("lmerTest", quietly = TRUE)
  # pbkrtest backs lmerTest's ddf = "Kenward-Roger"; without it KR is unavailable
  # and inference silently degrades to Satterthwaite.
  has_pbkrtest <- requireNamespace("pbkrtest", quietly = TRUE)
})

ts_msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), " | ", ..., "\n", sep = "")

# ---- CONFIG ----------------------------------------------------------------
PIPE_ROOT <- file.path(PROJECT_ROOT, "choice R pipeline")
BIN_S     <- 1.0        # sequence bin width in seconds (~25 frames @ 25 fps)
MIN_BINS  <- 10L        # drop a (trial,timepoint) sequence shorter than this
GRADED_LEVELS <- c("Low", "Med", "High")
BINARY_LEVELS <- c("Calm", "Flow")   # Calm-dominant, Flow-dominant

TREAT_COLORS <- c("control" = "#E69F00", "exercise choice" = "#0072B2")
STATE_FILL   <- c("Low" = "#FDE0A9", "Med" = "#F0A860", "High" = "#B8410E",
                  "Calm" = "#8FC7E8", "Flow" = "#0072B2")

# =============================================================================
# 1) LOCATE + LOAD FRAME-LEVEL DATA
# =============================================================================
.find_latest_master <- function() {
  parent <- file.path(PIPE_ROOT, "output", "STEP1_output")
  subs <- list.dirs(parent, full.names = TRUE, recursive = FALSE)
  subs <- subs[grepl("STEP1_output_\\d{8}_\\d{6}$", basename(subs))]
  if (!length(subs)) stop("No STEP1_output_* folder found under ", parent)
  latest <- subs[which.max(file.mtime(subs))]
  f <- file.path(latest, "master_fish_by_frame.csv")
  if (!file.exists(f)) stop("master_fish_by_frame.csv not found in ", latest)
  f
}

master_csv <- .find_latest_master()
ts_msg("Reading (selected cols): ", master_csv)
dt <- data.table::fread(
  master_csv,
  select = c("trial_id", "trial", "frame", "time", "timepoint", "treatment",
             "fish_density", "tank", "main_zone"),
  showProgress = FALSE
)
ts_msg("Loaded ", format(nrow(dt), big.mark = ","), " detection-rows")

# Keep only detections assigned to one of the two competing main zones. The
# exercise "choice" is flow vs calm, so that is the natural denominator.
dt <- dt[main_zone %in% c("flow", "calm")]
dt[, treatment := tolower(trimws(treatment))]
dt <- dt[treatment %in% names(TREAT_COLORS)]
ts_msg("After zone/treatment filter: ", format(nrow(dt), big.mark = ","), " rows")

# =============================================================================
# 2) PER-FRAME OCCUPANCY  ->  PER-BIN prop_flow
# =============================================================================
# Per frame: how many of the zoned school are in flow vs calm.
frame_occ <- dt[, .(
    n_flow = sum(main_zone == "flow"),
    n_calm = sum(main_zone == "calm"),
    time   = time[1]
  ), by = .(trial_id, trial, timepoint, frame, treatment, tank, fish_density)]

# Assign each frame to a fixed time bin within its (trial,timepoint).
frame_occ[, bin := floor(time / BIN_S)]

# Per bin: pooled flow proportion = flow detections / (flow + calm detections).
bin_occ <- frame_occ[, .(
    n_flow = sum(n_flow),
    n_calm = sum(n_calm)
  ), by = .(trial_id, trial, timepoint, treatment, tank, fish_density, bin)]
bin_occ[, n_zoned := n_flow + n_calm]
bin_occ <- bin_occ[n_zoned > 0]
bin_occ[, prop_flow := n_flow / n_zoned]
data.table::setkey(bin_occ, trial_id, timepoint, bin)
ts_msg("Built ", format(nrow(bin_occ), big.mark = ","), " time-bins (", BIN_S, "s)")

# =============================================================================
# 3) STATE CLASSIFICATION (both alphabets)
# =============================================================================
# (1) GRADED: data-driven terciles of the pooled bin prop_flow distribution.
q <- stats::quantile(bin_occ$prop_flow, probs = c(1/3, 2/3), na.rm = TRUE)
cut_lo <- unname(q[1]); cut_hi <- unname(q[2])
ts_msg(sprintf("Graded tercile cut-points: Low<%.3f<=Med<%.3f<=High", cut_lo, cut_hi))
bin_occ[, state_graded := factor(
  fifelse(prop_flow < cut_lo, "Low",
    fifelse(prop_flow < cut_hi, "Med", "High")),
  levels = GRADED_LEVELS)]

# (2) BINARY: majority of the zoned school in flow.
bin_occ[, state_binary := factor(
  fifelse(prop_flow > 0.5, "Flow", "Calm"), levels = BINARY_LEVELS)]

# =============================================================================
# 4) SEQUENCE METRICS PER (trial x timepoint)
# =============================================================================
# Shannon entropy (natural log), returns 0 for degenerate distributions.
.shannon <- function(p) { p <- p[p > 0]; if (!length(p)) return(0); -sum(p * log(p)) }

# Compute transition matrix, dwell times, switch rate, entropy rate for ONE
# ordered state vector. Transitions counted only between temporally ADJACENT
# bins (bin index differs by exactly 1) so gaps never bridge.
.seq_metrics <- function(states, bins, levels_vec, bin_s) {
  states <- as.character(states)
  ord <- order(bins); states <- states[ord]; bins <- bins[ord]
  K <- length(levels_vec)
  n <- length(states)

  # --- transition counts (adjacent bins only) ---
  tc <- matrix(0, K, K, dimnames = list(levels_vec, levels_vec))
  if (n >= 2) {
    adj <- which(diff(bins) == 1L)
    for (i in adj) tc[states[i], states[i + 1L]] <- tc[states[i], states[i + 1L]] + 1
  }
  row_tot <- rowSums(tc)
  # Row-normalised transition probabilities. A row whose origin state is NEVER
  # visited has row_tot == 0, so P(to | from) is UNDEFINED (0/0), not zero.
  # Emitting 0 there would be a spurious measured value: it silently biases the
  # per-session persistence metrics downward for exactly those groups that most
  # often skip a state. Such rows are therefore NA, and the downstream models
  # drop them by listwise deletion (n reported per metric).
  tp <- tc / row_tot                              # row_tot == 0 -> NaN
  tp[!is.finite(tp)] <- NA_real_
  # n transitions observed FROM each state — the denominator behind each tp row,
  # carried through so sparse rows can be identified and reported.
  n_from <- row_tot

  # --- occupancy (time share of each state) ---
  occ <- table(factor(states, levels = levels_vec)) / n

  # --- switch rate: off-diagonal transitions per minute ---
  n_changes <- sum(tc) - sum(diag(tc))
  obs_min   <- (n * bin_s) / 60
  switch_rate <- if (obs_min > 0) n_changes / obs_min else NA_real_

  # --- entropy rate: occupancy-weighted mean row entropy, normalised to [0,1] ---
  # Rows for never-visited states are NA (see above) but also carry occupancy 0,
  # so they contribute nothing: treat their entropy as 0 rather than letting
  # 0 * NA propagate NA through the weighted sum.
  H_rows <- apply(tp, 1, function(p) .shannon(p[is.finite(p)]))
  H_rows[!is.finite(H_rows)] <- 0
  ent_rate <- if (K > 1) sum(as.numeric(occ) * H_rows) / log(K) else 0

  # --- mean dwell (bout) time per state, INTERIOR bouts only (uncensored) ---
  # Split the sequence at gaps, run-length encode each contiguous segment,
  # drop the first & last run of each segment (censored by boundaries/gaps).
  dwell <- setNames(rep(NA_real_, K), levels_vec)
  seg_id <- cumsum(c(1L, as.integer(diff(bins) != 1L)))
  bouts <- list()
  for (s in split(seq_len(n), seg_id)) {
    r <- rle(states[s])
    if (length(r$lengths) >= 3) {           # need >=1 interior run
      inner <- 2:(length(r$lengths) - 1)
      for (j in inner) bouts[[length(bouts) + 1L]] <-
        data.frame(state = r$values[j], len = r$lengths[j])
    }
  }
  if (length(bouts)) {
    bd <- do.call(rbind, bouts)
    md <- tapply(bd$len, bd$state, mean) * bin_s
    dwell[names(md)] <- as.numeric(md)
  }

  list(tc = tc, tp = tp, occ = occ, switch_rate = switch_rate,
       ent_rate = ent_rate, dwell = dwell, n_bins = n, obs_min = obs_min,
       n_changes = n_changes, n_from = n_from)
}

# Iterate over (trial,timepoint) for a given alphabet; return per-unit metrics
# table + accumulated per-treatment transition-count matrices.
.run_alphabet <- function(state_col, levels_vec) {
  keys <- unique(bin_occ[, .(trial_id, trial, timepoint, treatment, tank, fish_density)])
  rows <- vector("list", nrow(keys))
  pooled <- list()  # treatment -> KxK count matrix
  for (t in unique(keys$treatment))
    pooled[[t]] <- matrix(0, length(levels_vec), length(levels_vec),
                          dimnames = list(levels_vec, levels_vec))

  for (i in seq_len(nrow(keys))) {
    k <- keys[i]
    sub <- bin_occ[trial_id == k$trial_id & timepoint == k$timepoint]
    if (nrow(sub) < MIN_BINS) next
    m <- .seq_metrics(sub[[state_col]], sub$bin, levels_vec, BIN_S)
    pooled[[k$treatment]] <- pooled[[k$treatment]] + m$tc

    row <- data.table(
      trial_id = k$trial_id, trial = k$trial, timepoint = k$timepoint,
      treatment = k$treatment, tank = k$tank, fish_density = k$fish_density,
      n_bins = m$n_bins, obs_min = m$obs_min,
      switch_rate = m$switch_rate, entropy_rate = m$ent_rate,
      n_changes = m$n_changes)
    for (s in levels_vec) {
      row[[paste0("occ_", s)]]   <- as.numeric(m$occ[s])
      row[[paste0("dwell_", s)]] <- unname(m$dwell[s])
      # PERSISTENCE: P(stay in s | currently in s) — the per-session, per-state
      # self-transition probability. Computed for EVERY alphabet (previously
      # binary-only), so the "school stays committed to the flow zone" claim can
      # be tested inferentially rather than read off pooled descriptives. NA when
      # state s is never entered (n_from == 0); n_from_* records the denominator.
      row[[paste0("p_stay_", s)]]  <- unname(m$tp[s, s])
      row[[paste0("n_from_", s)]]  <- unname(m$n_from[s])
    }
    # Directional off-diagonal probabilities for the 2-state alphabet.
    if (length(levels_vec) == 2L) {
      row$p_flow_to_calm <- m$tp["Flow", "Calm"]
      row$p_calm_to_flow <- m$tp["Calm", "Flow"]
    }
    rows[[i]] <- row
  }
  metrics <- data.table::rbindlist(Filter(Negate(is.null), rows), fill = TRUE)

  # Pooled transition PROBABILITIES per treatment (row-normalised).
  #
  # IMPORTANT — these are COUNT-WEIGHTED and are a DESCRIPTIVE summary only.
  # They answer "at a random second spent in state s, how likely is the school
  # to still be in s next second?", so schools that dwell longest in a state
  # dominate the estimate. The inferential persistence statistic is the
  # per-session p_stay_* metric above, which weights every school equally and
  # is the quantity carried into the Level A / Level B models.
  # The two can differ in DIRECTION when occupancy itself differs by treatment
  # (it does here: exercise-choice schools occupy Flow far more), so the pooled
  # matrix must never be quoted as evidence of a treatment difference.
  pooled_tp <- lapply(names(pooled), function(t) {
    cm <- pooled[[t]]; rt <- rowSums(cm)
    tp <- cm / ifelse(rt == 0, 1, rt)
    d <- as.data.table(as.table(tp))
    setnames(d, c("from", "to", "prob")); d[, treatment := t][]
  })
  list(metrics = metrics, pooled = data.table::rbindlist(pooled_tp))
}

ts_msg("Computing sequence metrics — GRADED alphabet ...")
res_graded <- .run_alphabet("state_graded", GRADED_LEVELS)
ts_msg("Computing sequence metrics — BINARY alphabet ...")
res_binary <- .run_alphabet("state_binary", BINARY_LEVELS)

# =============================================================================
# 5) TREATMENT CONTRASTS — TWO SEPARATE LEVELS OF ANALYSIS
# =============================================================================
# Design: 16 physical schools (`trial`), 3 timepoints each = 48 sessions;
# treatment is between-school, balanced within each of 4 tanks, fully crossed
# (8 schools per treatment, 8 sessions per treatment x timepoint cell).
#
# The two levels below answer DIFFERENT questions and are kept separate in the
# code, in the exported tables, in the manuscript and in the figures:
#
#   LEVEL A — TRIAL LEVEL  (n = 16 schools; the experimental unit)
#       Each metric is first averaged over that school's three sessions, giving
#       ONE value per school, then
#           val ~ treatment + fish_density + factor(tank)
#       There are no repeated measures left after aggregation, so `trial` does
#       not appear as a random effect; `tank` (4 levels, below the random-effect
#       stability floor of ~5-8) is retained as a FIXED blocking factor instead
#       of a random intercept (design-based RE policy).
#       Q: does the treatment shift a school's OVERALL sequence dynamics?
#
#   LEVEL B — TRIAL x TIMEPOINT LEVEL  (n = 48 sessions)
#       All sessions retained, with timepoint as an explicit CATEGORICAL fixed
#       effect and its interaction with treatment, and repeated measures on the
#       school carried by a random intercept:
#           val ~ treatment * timepoint_f + fish_density + (1|trial)
#       Q: does the treatment effect DEVELOP across the three intervals, i.e.
#       is there a treatment x interval interaction?
#
# Level B is not a re-plot of Level A: it estimates two extra fixed-effect terms
# (interval, treatment x interval) that Level A cannot express, and Level A is
# not a subset of Level B because its response is the school-level mean.
#
# Both levels use Type III ANOVA. Effect size for every ANOVA term is partial
# eta-squared computed from that term's own F and degrees of freedom,
#     eta^2_p = (F * df1) / (F * df1 + df2),
# the estimator already used elsewhere in this project; partial omega-squared
#     omega^2_p = (df1 * (F - 1)) / (df1 * (F - 1) + df2 + 1)
# is reported alongside as the small-sample bias-corrected counterpart, clamped
# at 0. No multiple-comparison correction is applied (project-wide decision,
# 2026-08-07): all p-values are raw.
#
#   ROBUST (retained, unchanged): Wilcoxon on school-level means, 8 vs 8.

# ---- effect sizes from an F-test -------------------------------------------
.eta2_p <- function(F_val, df1, df2) {
  if (!is.finite(F_val) || !is.finite(df1) || !is.finite(df2)) return(NA_real_)
  (F_val * df1) / (F_val * df1 + df2)
}
.omega2_p <- function(F_val, df1, df2) {
  if (!is.finite(F_val) || !is.finite(df1) || !is.finite(df2)) return(NA_real_)
  max(0, (df1 * (F_val - 1)) / (df1 * (F_val - 1) + df2 + 1))
}
has_effectsize <- requireNamespace("effectsize", quietly = TRUE)
has_perf <- requireNamespace("performance", quietly = TRUE)
# R2m/R2c (Nakagawa & Schielzeth 2013) + ICC (adjusted) via 'performance'.
# R2c/ICC are NA for a plain lm fit (no random effect to partition).
.r2_icc <- function(fit) {
  na_out <- list(R2m = NA_real_, R2c = NA_real_, ICC = NA_real_)
  if (is.null(fit) || !has_perf) return(na_out)
  r2 <- tryCatch(performance::r2(fit), error = function(e) NULL)
  r2m <- if (!is.null(r2)) as.numeric(r2$R2_marginal %||% r2$R2[1]) else NA_real_
  r2c <- if (!is.null(r2)) as.numeric(r2$R2_conditional %||% NA_real_) else NA_real_
  icc <- if (inherits(fit, "merMod"))
    tryCatch(as.numeric(performance::icc(fit)$ICC_adjusted), error = function(e) NA_real_)
  else NA_real_
  list(R2m = r2m, R2c = r2c, ICC = icc)
}
# Noncentral-F 95% CI on eta2_p/omega2_p (Nakagawa & Cuthill 2007 §II.4).
# Returns c(lo, hi); c(NA, NA) when effectsize is unavailable or inputs invalid.
.eta2_p_ci <- function(F_val, df1, df2) {
  if (!has_effectsize || !is.finite(F_val) || !is.finite(df1) || !is.finite(df2))
    return(c(NA_real_, NA_real_))
  es <- tryCatch(effectsize::F_to_eta2(F_val, df1, df2, ci = 0.95, alternative = "two.sided"),
                 error = function(e) NULL)
  if (is.null(es) || !nrow(es)) return(c(NA_real_, NA_real_))
  c(max(es$CI_low[1], 0), min(es$CI_high[1], 1))
}
.omega2_p_ci <- function(F_val, df1, df2) {
  if (!has_effectsize || !is.finite(F_val) || !is.finite(df1) || !is.finite(df2))
    return(c(NA_real_, NA_real_))
  os <- tryCatch(effectsize::F_to_omega2(F_val, df1, df2, ci = 0.95, alternative = "two.sided"),
                 error = function(e) NULL)
  if (is.null(os) || !nrow(os)) return(c(NA_real_, NA_real_))
  c(max(os$CI_low[1], 0), min(os$CI_high[1], 1))
}
# Signed Hedges' g + noncentral-t 95% CI (Nakagawa & Cuthill 2007), from a
# fitted model's own treatment contrast t-ratio and denominator df (KR by
# default via emmeans' lmer.df option, verified 2026-08-08) -- not a Wald
# normal-approximation CI, and not |estimate| (sign is retained).
.hedges_g <- function(fit, group_var = "treatment") {
  na_out <- list(g = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_, method = "na")
  if (is.null(fit) || !requireNamespace("emmeans", quietly = TRUE)) return(na_out)
  tryCatch({
    em  <- emmeans::emmeans(fit, stats::as.formula(paste("~", group_var)))
    ctr <- as.data.frame(emmeans::contrast(em, method = "pairwise", infer = c(TRUE, TRUE)))
    if (!nrow(ctr)) return(na_out)
    df_ctr <- ctr$df[1]
    t_ctr  <- if ("t.ratio" %in% names(ctr)) ctr$t.ratio[1] else ctr$estimate[1] / ctr$SE[1]
    if (!is.finite(df_ctr) || !is.finite(t_ctr) || df_ctr <= 1 || !has_effectsize) return(na_out)
    J     <- 1 - 3 / (4 * df_ctr - 1)
    d_res <- effectsize::t_to_d(t_ctr, df_ctr, ci = 0.95)
    if (!nrow(d_res)) return(na_out)
    list(g = d_res$d[1] * J, ci_lo = d_res$CI_low[1] * J, ci_hi = d_res$CI_high[1] * J,
         method = "t_to_d_noncentral")
  }, error = function(e) na_out)
}

# ---- Type III ANOVA table from lmer OR lm, in one shape --------------------
# Returns data.frame(term, F, df1, df2, p, stat_type) or NULL.
#
# CHANGED 2026-08-08: denominator df is now Kenward-Roger by default, matching
# Methods 2.8 ("Kenward-Roger ... with Satterthwaite computed as a sensitivity
# check") and the STEP5 behavioural engine. Previously this called
# stats::anova(fit, type = 3) with no `ddf`, which takes lmerTest's default of
# Satterthwaite -- so the sequence and bout modules were the only place in the
# study not using the convention the manuscript states.
#
# `stat_type` records what was actually computed, so a KR failure can never be
# mistaken for a KR success: F-KR | F-SW | F-SW(KR-failed) | F-OLS.
.anova3 <- function(fit, ddf = c("Kenward-Roger", "Satterthwaite")) {
  ddf <- match.arg(ddf)
  if (inherits(fit, "merMod") || inherits(fit, "lmerModLmerTest")) {
    # lmerTest honours ddf=; a plain lme4 fit must be promoted first.
    ml <- if (inherits(fit, "lmerModLmerTest")) fit
          else tryCatch(lmerTest::as_lmerModLmerTest(fit), error = function(e) NULL)
    if (is.null(ml)) return(NULL)

    use_ddf <- if (ddf == "Kenward-Roger" && !has_pbkrtest) "Satterthwaite" else ddf
    stat_lab <- if (use_ddf == "Kenward-Roger") "F-KR" else "F-SW"
    av <- tryCatch(suppressWarnings(stats::anova(ml, type = 3, ddf = use_ddf)),
                   error = function(e) NULL)

    # KR can fail on near-singular fits; fall back to Satterthwaite but SAY SO
    # rather than silently reporting it as KR.
    if (is.null(av) && use_ddf == "Kenward-Roger") {
      av <- tryCatch(suppressWarnings(stats::anova(ml, type = 3, ddf = "Satterthwaite")),
                     error = function(e) NULL)
      if (!is.null(av)) {
        stat_lab <- "F-SW(KR-failed)"
        ts_msg("    WARNING: Kenward-Roger failed; fell back to Satterthwaite")
      }
    }
    if (is.null(av) || !nrow(av)) return(NULL)
    return(data.frame(term = rownames(av), F = av[["F value"]],
                      df1 = av[["NumDF"]], df2 = av[["DenDF"]],
                      p = av[["Pr(>F)"]], stat_type = stat_lab,
                      stringsAsFactors = FALSE))
  }
  # plain lm: use car::Anova type III when available (design is balanced, so
  # this coincides with the sequential decomposition for the treatment term),
  # else fall back to the sequential table.
  if (requireNamespace("car", quietly = TRUE)) {
    av <- tryCatch(suppressWarnings(car::Anova(fit, type = "III")),
                   error = function(e) NULL)
    if (!is.null(av)) {
      d <- as.data.frame(av)
      d$term <- rownames(d)
      d <- d[!d$term %in% c("(Intercept)", "Residuals"), , drop = FALSE]
      dfres <- stats::df.residual(fit)
      return(data.frame(term = d$term, F = d[["F value"]], df1 = d[["Df"]],
                        df2 = dfres, p = d[["Pr(>F)"]], stat_type = "F-OLS",
                        stringsAsFactors = FALSE))
    }
  }
  av <- tryCatch(stats::anova(fit), error = function(e) NULL)
  if (is.null(av)) return(NULL)
  d <- as.data.frame(av); d$term <- rownames(d)
  d <- d[d$term != "Residuals", , drop = FALSE]
  data.frame(term = d$term, F = d[["F value"]], df1 = d[["Df"]],
             df2 = stats::df.residual(fit), p = d[["Pr(>F)"]], stat_type = "F-OLS",
             stringsAsFactors = FALSE)
}

# ---- fit the first non-singular candidate model ----------------------------
# Returns list(fit, form, singular) or NULL. Mirrors the existing convention of
# preferring the most complex non-singular structure.
.fit_first_ok <- function(forms, d) {
  fallback <- NULL
  for (fm in forms) {
    is_mixed <- grepl("\\(1\\s*\\|", fm)
    fit <- tryCatch(suppressMessages(suppressWarnings(
      if (!is_mixed) stats::lm(stats::as.formula(fm), data = d)
      else if (has_lmerTest) lmerTest::lmer(stats::as.formula(fm), data = d)
      else lme4::lmer(stats::as.formula(fm), data = d, REML = TRUE))),
      error = function(e) NULL)
    if (is.null(fit)) next
    sing <- if (is_mixed)
      isTRUE(tryCatch(lme4::isSingular(fit), error = function(e) TRUE)) else FALSE
    # Overwrite each pass so that, if EVERY candidate is singular, the fallback
    # is the SIMPLEST fit tried rather than the most complex — matching the
    # convention used by the project's other engines.
    fallback <- list(fit = fit, form = fm, singular = sing)
    if (!sing) return(list(fit = fit, form = fm, singular = FALSE))
  }
  fallback
}

# ---- LEVEL A: TRIAL LEVEL (n = 16 schools) ---------------------------------
.level_trial <- function(dat, metric, alphabet) {
  if (!metric %in% names(dat)) return(NULL)
  d <- dat[is.finite(get(metric)),
           .(trial, treatment, tank, fish_density, val = get(metric))]
  if (!nrow(d)) return(NULL)
  # Aggregate each school's sessions to a single value: after this there is one
  # row per school, so no repeated-measures term is required.
  per_school <- d[, .(val = mean(val, na.rm = TRUE),
                      n_sessions = .N,
                      tank = tank[1], fish_density = fish_density[1]),
                  by = .(trial, treatment)]
  if (data.table::uniqueN(per_school$treatment) < 2) return(NULL)
  g_ctrl <- per_school[treatment == "control", val]
  g_exer <- per_school[treatment == "exercise choice", val]
  if (length(g_ctrl) < 2 || length(g_exer) < 2) return(NULL)

  w <- suppressWarnings(stats::wilcox.test(g_exer, g_ctrl))
  # Trial-level (one row per school): no random effect is identifiable, and the
  # previous (1|tank) had only 4 levels — below the 5-8 stability floor (Bolker
  # 2009; Harrison 2018). Fit OLS with tank as a FIXED covariate instead.
  pick <- .fit_first_ok(c("val ~ treatment + fish_density + factor(tank)",
                          "val ~ treatment + fish_density",
                          "val ~ treatment"), per_school)
  av <- if (!is.null(pick)) .anova3(pick$fit) else NULL
  tr <- if (!is.null(av)) av[grepl("^treatment", av$term), , drop = FALSE] else NULL
  # Genuine Satterthwaite sensitivity row (see .anova3 header).
  av_sw <- if (!is.null(pick)) .anova3(pick$fit, ddf = "Satterthwaite") else NULL
  tr_sw <- if (!is.null(av_sw)) av_sw[grepl("^treatment", av_sw$term), , drop = FALSE] else NULL
  hg    <- if (!is.null(pick)) .hedges_g(pick$fit) else list(g = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_)
  r2i   <- if (!is.null(pick)) .r2_icc(pick$fit) else list(R2m = NA_real_, R2c = NA_real_, ICC = NA_real_)

  data.table(
    level = "trial", alphabet = alphabet, metric = metric, term = "treatment",
    n_unit = nrow(per_school),
    n_ctrl = length(g_ctrl), n_exer = length(g_exer),
    mean_ctrl = mean(g_ctrl), mean_exer = mean(g_exer),
    F = if (!is.null(tr) && nrow(tr)) tr$F[1] else NA_real_,
    df1 = if (!is.null(tr) && nrow(tr)) tr$df1[1] else NA_real_,
    df2 = if (!is.null(tr) && nrow(tr)) tr$df2[1] else NA_real_,
    p = if (!is.null(tr) && nrow(tr)) tr$p[1] else NA_real_,
    stat_type = if (!is.null(tr) && nrow(tr)) tr$stat_type[1] else NA_character_,
    F_satt  = if (!is.null(tr_sw) && nrow(tr_sw)) tr_sw$F[1] else NA_real_,
    df2_satt = if (!is.null(tr_sw) && nrow(tr_sw)) tr_sw$df2[1] else NA_real_,
    p_satt  = if (!is.null(tr_sw) && nrow(tr_sw)) tr_sw$p[1] else NA_real_,
    hedges_g = hg$g, hedges_g_lo = hg$ci_lo, hedges_g_hi = hg$ci_hi,
    R2m = r2i$R2m, R2c = r2i$R2c, ICC = r2i$ICC,
    W = unname(w$statistic), p_wilcox = w$p.value,
    model = if (!is.null(pick)) pick$form else NA_character_,
    singular = if (!is.null(pick)) pick$singular else NA)
}

# ---- LEVEL B: TRIAL x TIMEPOINT LEVEL (n = 48 sessions) --------------------
.level_trial_tp <- function(dat, metric, alphabet) {
  if (!metric %in% names(dat)) return(NULL)
  d <- dat[is.finite(get(metric)),
           .(trial, treatment, tank, fish_density, timepoint, val = get(metric))]
  if (!nrow(d)) return(NULL)
  d[, timepoint_f := factor(timepoint)]
  if (data.table::uniqueN(d$timepoint_f) < 2) return(NULL)
  if (data.table::uniqueN(d$trial) >= nrow(d)) return(NULL)  # no repeats left

  # Interval-level (N=48, 3/trial): force (1|trial) unconditionally; tank (4
  # levels) is dropped from the candidate set since trial is nested within tank
  # and already absorbs tank-level variation (design-based RE policy).
  pick <- .fit_first_ok(c(
    "val ~ treatment * timepoint_f + fish_density + (1|trial)",
    "val ~ treatment * timepoint_f + (1|trial)"), d)
  if (is.null(pick)) return(NULL)
  av <- .anova3(pick$fit)
  if (is.null(av)) return(NULL)
  keep <- av[av$term %in% c("treatment", "timepoint_f", "treatment:timepoint_f"), ,
             drop = FALSE]
  if (!nrow(keep)) return(NULL)

  # Genuine Satterthwaite sensitivity, matched term-by-term.
  av_sw <- .anova3(pick$fit, ddf = "Satterthwaite")
  m_sw <- if (!is.null(av_sw)) match(keep$term, av_sw$term) else rep(NA_integer_, nrow(keep))

  # Hedges' g is only defined for a genuine 2-group contrast (the "treatment"
  # term); NA for timepoint_f (3 levels) and the interaction row.
  hg <- .hedges_g(pick$fit)
  hg_g <- ifelse(keep$term == "treatment", hg$g, NA_real_)
  hg_lo <- ifelse(keep$term == "treatment", hg$ci_lo, NA_real_)
  hg_hi <- ifelse(keep$term == "treatment", hg$ci_hi, NA_real_)
  r2i <- .r2_icc(pick$fit)

  data.table(
    level = "trial_x_timepoint", alphabet = alphabet, metric = metric,
    term = keep$term, n_unit = nrow(d),
    F = keep$F, df1 = keep$df1, df2 = keep$df2, p = keep$p,
    stat_type = keep$stat_type,
    F_satt   = if (!is.null(av_sw)) av_sw$F[m_sw]   else NA_real_,
    df2_satt = if (!is.null(av_sw)) av_sw$df2[m_sw] else NA_real_,
    p_satt   = if (!is.null(av_sw)) av_sw$p[m_sw]   else NA_real_,
    hedges_g = hg_g, hedges_g_lo = hg_lo, hedges_g_hi = hg_hi,
    R2m = r2i$R2m, R2c = r2i$R2c, ICC = r2i$ICC,
    model = pick$form, singular = pick$singular)
}

# ---- legacy pooled contrast (retained for the descriptive means) -----------
.contrast_one <- function(dat, metric, alphabet) {
  d <- dat[is.finite(get(metric)),
           .(trial, treatment, tank, fish_density, timepoint, val = get(metric))]
  if (!nrow(d)) return(NULL)

  # ROBUST: collapse to one value per physical school, then between-treatment test.
  per_school <- d[, .(val = mean(val, na.rm = TRUE)), by = .(trial, treatment)]
  g_ctrl <- per_school[treatment == "control", val]
  g_exer <- per_school[treatment == "exercise choice", val]
  if (length(g_ctrl) < 2 || length(g_exer) < 2) return(NULL)
  w <- suppressWarnings(stats::wilcox.test(g_exer, g_ctrl))

  out <- data.table(
    alphabet = alphabet, metric = metric,
    n_school_ctrl = length(g_ctrl), n_school_exer = length(g_exer),
    mean_ctrl = mean(g_ctrl), mean_exer = mean(g_exer),
    W = unname(w$statistic), p_wilcox = w$p.value,
    beta_lmm = NA_real_, F_lmm = NA_real_, df1_lmm = NA_real_,
    df2_lmm = NA_real_, p_lmm = NA_real_, lmm_form = NA_character_)

  # PRIMARY: nested LMM. Prefer the most complex NON-singular fit; only fall back
  # to a singular fit if every candidate is singular. Interval-level: force
  # (1|trial) unconditionally; tank (4 levels) is dropped from the candidate set
  # since trial is nested within tank and already absorbs its variation
  # (design-based RE policy).
  if (has_lme4 && data.table::uniqueN(d$trial) < nrow(d)) {
    forms <- c("val ~ treatment + fish_density + (1|trial)",
               "val ~ treatment + (1|trial)")
    chosen <- NULL; fallback <- NULL
    for (fm in forms) {
      fit <- tryCatch(suppressMessages(suppressWarnings(
        if (has_lmerTest) lmerTest::lmer(stats::as.formula(fm), data = d)
        else lme4::lmer(stats::as.formula(fm), data = d, REML = TRUE))),
        error = function(e) NULL)
      if (is.null(fit)) next
      fallback <- list(fit = fit, fm = fm)              # simplest valid fit seen
      sing <- isTRUE(tryCatch(lme4::isSingular(fit), error = function(e) TRUE))
      if (!sing) { chosen <- list(fit = fit, fm = fm); break }
    }
    pick <- if (!is.null(chosen)) chosen else fallback
    if (!is.null(pick)) {
      cf <- suppressWarnings(summary(pick$fit))$coefficients
      trow <- grep("^treatment", rownames(cf))
      if (length(trow)) {
        out$beta_lmm <- cf[trow[1], "Estimate"]   # retained for effect direction
        out$lmm_form <- pick$fm
      }
      # Type III ANOVA of the fitted LMM. Kenward-Roger F since 2026-08-08, via
      # the shared .anova3() helper, so this path cannot drift from the two
      # inferential levels above. For the two-level treatment factor it is the
      # exact equivalent of the coefficient t-test (F = t^2).
      if (has_lmerTest) {
        av <- .anova3(pick$fit)
        if (!is.null(av) && "treatment" %in% av$term) {
          ar <- av[av$term == "treatment", , drop = FALSE]
          out$F_lmm   <- ar$F[1]
          out$df1_lmm <- ar$df1[1]
          out$df2_lmm <- ar$df2[1]
          out$p_lmm   <- ar$p[1]
        }
      }
    }
  }
  out
}

# Metric sets analysed at BOTH levels. Persistence probabilities (p_stay_*) are
# now included for the graded alphabet too, so the "commitment to the flow zone"
# claim is testable rather than descriptive. obs_min is the observation-effort
# balance check (a confounder control for the dwell/occupancy metrics).
graded_metrics_to_test <- c("switch_rate", "entropy_rate",
                            "occ_High", "occ_Low", "dwell_High", "dwell_Low",
                            "p_stay_High", "p_stay_Low", "obs_min")
# NB p_flow_to_calm is the exact complement of p_stay_Flow in a two-state
# alphabet (p_flow_to_calm = 1 - p_stay_Flow), so it carries identical
# information and returns an identical F-test. Only the persistence form is
# tested, to avoid reporting the same result twice.
binary_metrics_to_test <- c("switch_rate", "entropy_rate", "occ_Flow",
                            "dwell_Flow", "dwell_Calm",
                            "p_stay_Flow", "p_stay_Calm",
                            "obs_min")

# ---- LEVEL A: trial level (n = 16 schools) ---------------------------------
ts_msg("Level A — trial-level models (n = 16 schools) ...")
seq_trial <- data.table::rbindlist(c(
  lapply(graded_metrics_to_test, function(m) .level_trial(res_graded$metrics, m, "graded")),
  lapply(binary_metrics_to_test, function(m) .level_trial(res_binary$metrics, m, "binary"))
), fill = TRUE)

# ---- LEVEL B: trial x timepoint level (n = 48 sessions) --------------------
ts_msg("Level B — trial x timepoint models (n = 48 sessions) ...")
seq_trial_tp <- data.table::rbindlist(c(
  lapply(graded_metrics_to_test, function(m) .level_trial_tp(res_graded$metrics, m, "graded")),
  lapply(binary_metrics_to_test, function(m) .level_trial_tp(res_binary$metrics, m, "binary"))
), fill = TRUE)

# Effect sizes for every ANOVA term at both levels, with noncentral-F 95% CIs
# (Nakagawa & Cuthill 2007) -- reported for every term regardless of
# significance, not gated on p < 0.05.
for (D in list(seq_trial, seq_trial_tp)) {
  D[, eta2_p   := mapply(.eta2_p,   F, df1, df2)]
  D[, omega2_p := mapply(.omega2_p, F, df1, df2)]
  D[, c("eta2_p_lo", "eta2_p_hi")     := data.table::as.data.table(t(mapply(.eta2_p_ci,   F, df1, df2)))]
  D[, c("omega2_p_lo", "omega2_p_hi") := data.table::as.data.table(t(mapply(.omega2_p_ci, F, df1, df2)))]
}

# ---- Legacy pooled contrast: retained ONLY for descriptive treatment means --
# (kept so the existing per-metric mean columns and the robust Wilcoxon remain
# available; the inferential statements now come from Level A / Level B above).
contrasts <- data.table::rbindlist(c(
  lapply(graded_metrics_to_test, function(m) .contrast_one(res_graded$metrics, m, "graded")),
  lapply(binary_metrics_to_test, function(m) .contrast_one(res_binary$metrics, m, "binary"))
), fill = TRUE)

# =============================================================================
# 6) OUTPUT DIR
# =============================================================================
out_dir <- file.path(PIPE_ROOT, "output", "SEQ_output",
                     paste0("SEQ_output_", format(Sys.time(), "%Y%m%d_%H%M%S")))
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

data.table::fwrite(res_graded$metrics, file.path(out_dir, "seq_metrics_per_trial_graded.csv"))
data.table::fwrite(res_binary$metrics, file.path(out_dir, "seq_metrics_per_trial_binary.csv"))
data.table::fwrite(res_graded$pooled,  file.path(out_dir, "transition_matrix_pooled_graded.csv"))
data.table::fwrite(res_binary$pooled,  file.path(out_dir, "transition_matrix_pooled_binary.csv"))
data.table::fwrite(contrasts,          file.path(out_dir, "treatment_contrasts.csv"))

# The two analysis levels are exported as SEPARATE result tables.
data.table::fwrite(seq_trial,    file.path(out_dir, "seq_anova_trial_level.csv"))
data.table::fwrite(seq_trial_tp, file.path(out_dir, "seq_anova_trial_x_timepoint.csv"))

# ---- POST-HOC: only where a Level-B treatment x interval interaction is sig --
# Interpreting a significant omnibus interaction requires the simple treatment
# effect within each interval; emmeans supplies those on the model scale.
posthoc_rows <- list()
sig_int <- seq_trial_tp[term == "treatment:timepoint_f" & is.finite(p) & p < 0.05]
if (nrow(sig_int) && requireNamespace("emmeans", quietly = TRUE)) {
  for (i in seq_len(nrow(sig_int))) {
    rr <- sig_int[i]
    src <- if (rr$alphabet == "graded") res_graded$metrics else res_binary$metrics
    d <- src[is.finite(get(rr$metric)),
             .(trial, treatment, tank, fish_density, timepoint,
               val = get(rr$metric))]
    d[, timepoint_f := factor(timepoint)]
    fit <- tryCatch(suppressMessages(suppressWarnings(
      lmerTest::lmer(stats::as.formula(rr$model), data = d))),
      error = function(e) NULL)
    if (is.null(fit)) next
    em <- tryCatch(suppressMessages(
      emmeans::emmeans(fit, ~ treatment | timepoint_f)), error = function(e) NULL)
    if (is.null(em)) next
    pr <- as.data.frame(pairs(em, adjust = "none"))   # no correction, per project
    posthoc_rows[[length(posthoc_rows) + 1L]] <- data.table(
      alphabet = rr$alphabet, metric = rr$metric,
      timepoint = as.character(pr$timepoint_f), contrast = as.character(pr$contrast),
      estimate = pr$estimate, SE = pr$SE, df = pr$df,
      t = pr$t.ratio, p = pr$p.value)
  }
}
seq_posthoc <- if (length(posthoc_rows))
  data.table::rbindlist(posthoc_rows, fill = TRUE) else
  data.table(alphabet = character(), metric = character(), timepoint = character(),
             contrast = character(), estimate = numeric(), SE = numeric(),
             df = numeric(), t = numeric(), p = numeric())
data.table::fwrite(seq_posthoc, file.path(out_dir, "seq_posthoc_trial_x_timepoint.csv"))

# ---- POLYNOMIAL DECOMPOSITION of the interval effect (EVERY metric) ---------
# The three tracked intervals are equally spaced (midpoints 15, 55, 95 min), so
# the 2-df interval effect splits ORTHOGONALLY into a linear (monotone trend)
# and a quadratic (mid-trial dip) component. This characterises the SHAPE of an
# effect; the omnibus interaction in seq_anova_trial_x_timepoint.csv remains the
# primary test and alone governs whether a metric qualifies for reporting.
#
# NOTE the contrast with the post-hoc block above: that one is deliberately
# gated on a significant omnibus, because simple effects are only interpretable
# once the interaction is established. This block is deliberately NOT gated --
# it runs for every Level-B metric, significant or not, and writes a row set for
# each. Reporting the decomposition exhaustively is what keeps it free of
# selection; emitting it only where it happened to be significant would make it
# a second, unadjusted significance filter. adjust = "none" is correct for the
# same reason: the two components are orthogonal and pre-specified, so there is
# no multiplicity to correct within the decomposition.
.poly_lab <- function(x) {
  lv <- unique(as.character(x))
  c("linear", "quadratic", "cubic", "quartic")[match(as.character(x), lv)]
}
poly_rows <- list()
if (requireNamespace("emmeans", quietly = TRUE)) {
  mods <- unique(seq_trial_tp[!is.na(model), .(alphabet, metric, model)])
  for (i in seq_len(nrow(mods))) {
    rr  <- mods[i]
    src <- if (rr$alphabet == "graded") res_graded$metrics else res_binary$metrics
    if (!rr$metric %in% names(src)) next
    d <- src[is.finite(get(rr$metric)),
             .(trial, treatment, tank, fish_density, timepoint,
               val = get(rr$metric))]
    if (!nrow(d)) next
    d[, timepoint_f := factor(timepoint)]
    if (data.table::uniqueN(d$timepoint_f) < 3) next
    fit <- tryCatch(suppressMessages(suppressWarnings(
      lmerTest::lmer(stats::as.formula(rr$model), data = d))),
      error = function(e) NULL)
    if (is.null(fit)) next

    # (a) treatment x interval, split into its linear and quadratic components
    ci <- tryCatch(suppressMessages(as.data.frame(emmeans::contrast(
      emmeans::emmeans(fit, ~ timepoint_f * treatment),
      interaction = c("poly", "pairwise"), adjust = "none"))),
      error = function(e) NULL)
    if (!is.null(ci) && nrow(ci)) {
      pc <- names(ci)[grepl("_poly$", names(ci))]
      wc <- names(ci)[grepl("_pairwise$", names(ci))]
      poly_rows[[length(poly_rows) + 1L]] <- data.table(
        alphabet = rr$alphabet, metric = rr$metric,
        block = "treatment_x_interval",
        component = if (length(pc)) .poly_lab(ci[[pc]]) else NA_character_,
        contrast  = if (length(wc)) as.character(ci[[wc]]) else NA_character_,
        treatment = NA_character_,
        estimate = ci$estimate, SE = ci$SE, df = ci$df,
        t = ci$t.ratio, p = ci$p.value)
    }

    # (b) interval shape WITHIN each treatment arm -- describes each plotted line
    cw <- tryCatch(suppressMessages(as.data.frame(emmeans::contrast(
      emmeans::emmeans(fit, ~ timepoint_f | treatment), "poly", adjust = "none"))),
      error = function(e) NULL)
    if (!is.null(cw) && nrow(cw)) {
      poly_rows[[length(poly_rows) + 1L]] <- data.table(
        alphabet = rr$alphabet, metric = rr$metric,
        block = "interval_within_treatment",
        component = as.character(cw$contrast),
        contrast  = NA_character_,
        treatment = as.character(cw$treatment),
        estimate = cw$estimate, SE = cw$SE, df = cw$df,
        t = cw$t.ratio, p = cw$p.value)
    }
  }
}
seq_poly <- if (length(poly_rows))
  data.table::rbindlist(poly_rows, fill = TRUE) else
  data.table(alphabet = character(), metric = character(), block = character(),
             component = character(), contrast = character(),
             treatment = character(), estimate = numeric(), SE = numeric(),
             df = numeric(), t = numeric(), p = numeric())
data.table::fwrite(seq_poly, file.path(out_dir, "seq_poly_trial_x_timepoint.csv"))

saveRDS(list(cut_lo = cut_lo, cut_hi = cut_hi, BIN_S = BIN_S),
        file.path(out_dir, "config.rds"))

# =============================================================================
# 7) FIGURES
# =============================================================================
theme_seq <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 13),
        plot.caption = ggtext::element_markdown(hjust = 0.5, size = 10,
                                                lineheight = 1.2),
        legend.position = "bottom")

pretty_treat <- function(x) stringr::str_to_title(x)

# --- STATEMENT builders -----------------------------------------------------
# Matches the manuscript convention of quoting the inferential statistic within
# the panel, as ggtext markdown so degrees of freedom render as subscripts.
# Every ANOVA statement now carries its partial eta-squared effect size in the
# project's X = Y form.
# Unified number-display convention (author decision, 2026-08-07) -- see
# 00_shared/number_formatting.R. .fmt_F: F-statistics and df, 4 sig figs,
# trailing zeros trimmed. .fmt_es/.fmt_p: everything else, 3dp, single
# trailing-zero trim to 2 when the third decimal digit is exactly zero.
.fmt_F <- function(x) {
  s <- format(signif(x, 4), scientific = FALSE, trim = TRUE)
  # Only trim trailing zeros when a decimal point is present -- otherwise an
  # integer denominator df like 30 is mangled to "3" (see project number-
  # formatting defect log).
  has_dot <- grepl(".", s, fixed = TRUE)
  s[has_dot] <- sub("\\.$", "", sub("0+$", "", s[has_dot]))
  ifelse(is.finite(x), s, "NA")
}
.fmt_df <- .fmt_F

# .fmt_Fstat / .fmt_es2 (2026-08-09): as .fmt_F / .fmt_es, but a value below
# 0.001 collapses to "< 0.001" instead of rendering either a long 4-sig-fig
# decimal ("0.0009258") or an uninformative rounded zero ("0.00"). Applied ONLY
# to the test statistic and to partial eta-squared -- denominator df keep
# .fmt_df (a df is never < 1) and p-values keep .fmt_p (own floor).
.fmt_Fstat <- function(x) {
  ifelse(!is.finite(x), "NA",
         ifelse(abs(x) < 0.001, "< 0.001", .fmt_F(x)))
}

.fmt_es <- function(x) {
  s <- sprintf("%.3f", x)
  s <- ifelse(grepl("0$", s), substr(s, 1, nchar(s) - 1), s)
  ifelse(is.finite(x), s, "NA")
}

.fmt_es2 <- function(x) {
  ifelse(!is.finite(x), "NA",
         ifelse(abs(x) < 0.001, "< 0.001", .fmt_es(x)))
}

.fmt_p <- function(p) ifelse(!is.finite(p), "NA",
                     ifelse(p < 0.001, "< 0.001", .fmt_es(p)))

# Level A (trial-level) treatment effect + the robust school-level Wilcoxon.
seq_statement <- function(alph, metr) {
  r <- seq_trial[alphabet == alph & metric == metr & term == "treatment"]
  if (!nrow(r)) return(NULL)
  sprintf(paste0("*treatment: F<sub>%s,%s</sub> = %s, p = %s, ",
                 "&eta;<sup>2</sup><sub>p</sub> = %s*<br>",
                 "*Wilcoxon: W = %g, p = %s*"),
          .fmt_df(r$df1[1]), .fmt_df(r$df2[1]), .fmt_Fstat(r$F[1]), .fmt_p(r$p[1]),
          .fmt_es2(r$eta2_p[1]), r$W[1], .fmt_p(r$p_wilcox[1]))
}

# Level B (trial x timepoint): treatment, interval, and their interaction.
seq_statement_tp <- function(alph, metr) {
  r <- seq_trial_tp[alphabet == alph & metric == metr]
  if (!nrow(r)) return(NULL)
  lab <- c(treatment = "treatment", timepoint_f = "interval",
           `treatment:timepoint_f` = "treatment &times; interval")
  parts <- vapply(c("treatment", "timepoint_f", "treatment:timepoint_f"),
    function(tm) {
      x <- r[term == tm]
      if (!nrow(x)) return(NA_character_)
      sprintf("*%s: F<sub>%s,%s</sub> = %s, p = %s, &eta;<sup>2</sup><sub>p</sub> = %s*",
              lab[[tm]], .fmt_df(x$df1[1]), .fmt_df(x$df2[1]), .fmt_Fstat(x$F[1]),
              .fmt_p(x$p[1]), .fmt_es2(x$eta2_p[1]))
    }, character(1))
  paste(parts[!is.na(parts)], collapse = "<br>")
}

# --- (A) Example state-ribbon timelines: 2 trials/treatment, timepoint 1 -----
.ribbon_df <- function(state_col, levels_vec) {
  ex_keys <- bin_occ[timepoint == min(timepoint),
                     .(nb = .N), by = .(trial, treatment)][nb >= MIN_BINS]
  ex_keys <- ex_keys[, head(.SD, 2), by = treatment]
  d <- bin_occ[timepoint == min(timepoint) &
               trial %in% ex_keys$trial,
               .(trial, treatment, bin, prop_flow, st = get(state_col))]
  d[, lbl := paste0(pretty_treat(treatment), " — school ", trial)]
  d
}
rb <- .ribbon_df("state_graded", GRADED_LEVELS)
p_ribbon <- ggplot(rb, aes(x = bin * BIN_S, y = 1, fill = st)) +
  geom_tile() +
  facet_wrap(~ lbl, ncol = 1, strip.position = "left") +
  scale_fill_manual(values = STATE_FILL, name = "Collective flow state", drop = FALSE) +
  scale_y_continuous(breaks = NULL) +
  labs(title = "A  Collective exercise-choice state over time (graded alphabet)",
       x = "Time (s)", y = NULL) +
  theme_seq + theme(strip.text.y.left = element_text(angle = 0, hjust = 1))
ggsave(file.path(fig_dir, "A_state_ribbons_graded.png"), p_ribbon,
       width = 250, height = 150, units = "mm", dpi = 300, bg = "white")

# --- (B) Transition-matrix heatmaps by treatment (both alphabets) ------------
.heat <- function(pooled, levels_vec, ttl) {
  pooled[, from := factor(from, levels = levels_vec)]
  pooled[, to   := factor(to,   levels = rev(levels_vec))]
  ggplot(pooled, aes(x = to, y = from, fill = prob)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", prob)), size = 3.4) +
    facet_wrap(~ pretty_treat(treatment)) +
    scale_fill_gradient(low = "#F7FBFF", high = "#08519C", limits = c(0, 1),
                        name = "P(to | from)") +
    coord_equal() +
    labs(title = ttl, x = "to state", y = "from state") +
    theme_seq
}
p_heat_g <- .heat(copy(res_graded$pooled), GRADED_LEVELS,
                  "B  Transition probabilities — graded alphabet")
p_heat_b <- .heat(copy(res_binary$pooled), BINARY_LEVELS,
                  "C  Transition probabilities — binary alphabet")
ggsave(file.path(fig_dir, "B_transition_heat_graded.png"), p_heat_g,
       width = 220, height = 120, units = "mm", dpi = 300, bg = "white")
ggsave(file.path(fig_dir, "C_transition_heat_binary.png"), p_heat_b,
       width = 200, height = 110, units = "mm", dpi = 300, bg = "white")

# --- (D) Per-trial metric comparisons exercise vs control --------------------
.metric_box <- function(metrics, metric, ylab, ttl, alph) {
  d <- metrics[is.finite(get(metric)),
               .(val = mean(get(metric), na.rm = TRUE)),
               by = .(trial, treatment)]
  ggplot(d, aes(x = pretty_treat(treatment), y = val, color = treatment)) +
    geom_boxplot(outlier.shape = NA, width = 0.5) +
    geom_jitter(width = 0.12, height = 0, size = 2, alpha = 0.8) +
    scale_color_manual(values = TREAT_COLORS, guide = "none") +
    labs(title = ttl, x = NULL, y = ylab, caption = seq_statement(alph, metric)) +
    theme_seq
}
p_sw  <- .metric_box(res_binary$metrics, "switch_rate",
                     "State switches / min", "D  Switch rate (binary)", "binary")
p_ent <- .metric_box(res_graded$metrics, "entropy_rate",
                     "Normalised entropy rate", "E  Sequence predictability (graded)", "graded")
p_flow<- .metric_box(res_binary$metrics, "occ_Flow",
                     "Time Flow-dominant", "F  Flow-dominance occupancy", "binary")
p_dwell <- .metric_box(res_binary$metrics, "dwell_Flow",
                       "Mean flow-bout (s)", "G  Flow commitment (bout length)", "binary")
p_metrics <- (p_sw | p_ent) / (p_flow | p_dwell)
ggsave(file.path(fig_dir, "D_metric_comparisons.png"), p_metrics,
       width = 240, height = 200, units = "mm", dpi = 300, bg = "white")

# --- (H) Dwell-time distributions by state x treatment (binary) --------------
dw <- rbind(
  res_binary$metrics[, .(trial, treatment, state = "Flow", dwell = dwell_Flow)],
  res_binary$metrics[, .(trial, treatment, state = "Calm", dwell = dwell_Calm)]
)[is.finite(dwell)]
dw <- dw[, .(dwell = mean(dwell, na.rm = TRUE)), by = .(trial, treatment, state)]
# Per-facet ANOVA statements (one per state; LMM + Wilcoxon on two lines).
dw_stmt <- data.frame(
  state = c("Flow", "Calm"),
  txt   = c(seq_statement("binary", "dwell_Flow"),
            seq_statement("binary", "dwell_Calm")),
  stringsAsFactors = FALSE)
p_dwdist <- ggplot(dw, aes(x = pretty_treat(treatment), y = dwell, color = treatment)) +
  geom_boxplot(outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.12, height = 0, size = 1.8, alpha = 0.7) +
  facet_wrap(~ state) +
  ggtext::geom_richtext(data = dw_stmt, aes(x = 1.5, y = Inf, label = txt),
            inherit.aes = FALSE, vjust = 1.05, size = 3.1,
            fill = NA, label.color = NA, lineheight = 1.05) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
  scale_color_manual(values = TREAT_COLORS, guide = "none") +
  labs(title = "H  Dwell time per collective state (binary alphabet)",
       x = NULL, y = "Mean bout length (s)") +
  theme_seq
ggsave(file.path(fig_dir, "H_dwell_by_state.png"), p_dwdist,
       width = 200, height = 120, units = "mm", dpi = 300, bg = "white")

# =============================================================================
# 7b) RESULTS-GRID FIGURES  (I = complete, J = significant-only)
# =============================================================================
# Two summary grids over ALL analysed sequence metrics, both alphabets and both
# analysis levels. Same aesthetics as the panels above: theme_seq, the shared
# TREAT_COLORS palette, ggtext markdown labels, 300 dpi PNG.
#
# Encoding: one row per metric, one column per model term, faceted by alphabet.
# Fill = partial eta-squared (the effect size), with the F/p/eta2 label printed
# in-cell. The trial level contributes a single "treatment" column; the trial x
# timepoint level contributes "treatment", "interval" and "treatment x interval".
# This keeps the two levels visually distinct while allowing direct comparison.

TERM_LAB <- c(treatment = "treatment",
              timepoint_f = "interval",
              `treatment:timepoint_f` = "treatment × interval")
LEVEL_LAB <- c(trial = "Level A — trial (n = 16)",
               trial_x_timepoint = "Level B — trial × interval (n = 48)")

grid_dat <- rbind(
  seq_trial[, .(level, alphabet, metric, term, F, df1, df2, p, eta2_p)],
  seq_trial_tp[, .(level, alphabet, metric, term, F, df1, df2, p, eta2_p)]
)
grid_dat[, term_lab := factor(TERM_LAB[term],
                              levels = unname(TERM_LAB))]
grid_dat[, level_lab := factor(LEVEL_LAB[level], levels = unname(LEVEL_LAB))]
grid_dat[, sig := is.finite(p) & p < 0.05]
# Cell label: F, p and the effect size in the project's X = Y form.
grid_dat[, cell_lab := sprintf(
  "F<sub>%s,%s</sub> = %s<br>p = %s<br>&eta;<sup>2</sup><sub>p</sub> = %s",
  .fmt_df(df1), .fmt_df(df2), .fmt_Fstat(F), .fmt_p(p), .fmt_es2(eta2_p))]
# Order metrics consistently: dynamics first, then occupancy/dwell/persistence.
metric_order <- c("switch_rate", "entropy_rate",
                  "occ_Flow", "occ_High", "occ_Low",
                  "dwell_Flow", "dwell_Calm", "dwell_High", "dwell_Low",
                  "p_stay_Flow", "p_stay_Calm", "p_stay_High", "p_stay_Low",
                  "obs_min")
grid_dat[, metric := factor(metric, levels = rev(intersect(metric_order,
                                                           unique(metric))))]

# free scales on BOTH axes so each alphabet shows only the metrics it actually
# has (binary has *_Flow/*_Calm, graded has *_High/*_Low/*_Med); `space = "free"`
# then keeps every cell the same physical size, which is what stops the
# multi-line in-cell labels from colliding.
.results_grid <- function(d, ttl, sub) {
  d <- droplevels(d)
  ggplot(d, aes(x = term_lab, y = metric, fill = eta2_p)) +
    geom_tile(aes(colour = sig), linewidth = 0.7, width = 0.97, height = 0.97) +
    ggtext::geom_richtext(aes(label = cell_lab), size = 2.6,
                          fill = NA, label.color = NA, lineheight = 1.15,
                          label.padding = grid::unit(rep(0, 4), "pt")) +
    facet_grid(alphabet ~ level_lab, scales = "free", space = "free") +
    scale_fill_gradient(low = "#F7FBFF", high = "#08519C",
                        limits = c(0, 1), name = expression(eta[p]^2)) +
    scale_colour_manual(values = c(`TRUE` = "#B8410E", `FALSE` = "grey85"),
                        guide = "none") +
    labs(title = ttl, subtitle = sub, x = NULL, y = NULL) +
    theme_seq +
    theme(axis.text.x = element_text(angle = 20, hjust = 1),
          panel.grid = element_blank(),
          panel.spacing = grid::unit(6, "pt"))
}

p_grid_all <- .results_grid(
  grid_dat,
  "I  Behavioural sequence results — complete grid",
  paste0("All analysed metrics, both state alphabets, both analysis levels. ",
         "Fill and label give partial η²; red outline marks p < 0.05 (raw)."))
ggsave(file.path(fig_dir, "I_results_grid_complete.png"), p_grid_all,
       width = 330, height = 300, units = "mm", dpi = 300, bg = "white",
       limitsize = FALSE)

# --- Significant-only grid, SAME structure and SAME criterion (raw p < 0.05) --
sig_dat <- grid_dat[sig == TRUE]
if (nrow(sig_dat)) {
  p_grid_sig <- .results_grid(
    sig_dat,
    "J  Behavioural sequence results — significant only",
    paste0("Terms reaching p < 0.05 (raw, uncorrected — the same criterion as\n",
           "the analyses). Combinations with no significant term are left empty."))
} else {
  p_grid_sig <- ggplot() +
    annotate("text", x = 0, y = 0,
             label = "No sequence term reached p < 0.05", size = 5) +
    theme_void()
}
ggsave(file.path(fig_dir, "J_results_grid_significant.png"), p_grid_sig,
       width = 300, height = 200, units = "mm", dpi = 300, bg = "white",
       limitsize = FALSE)
data.table::fwrite(grid_dat[, .(level, alphabet, metric, term, F, df1, df2,
                                p, eta2_p, sig)],
                   file.path(out_dir, "seq_results_grid.csv"))

# =============================================================================
# 8) CONSOLE SUMMARY
# =============================================================================
ts_msg("==================== SEQUENCE ANALYSIS DONE ====================")
ts_msg("Output dir: ", out_dir)
cat("\n--- Graded tercile cut-points ---\n")
cat(sprintf("  Low < %.3f <= Med < %.3f <= High   (BIN_S=%.1fs)\n", cut_lo, cut_hi, BIN_S))
cat("\n--- Units analysed (>= ", MIN_BINS, " bins) ---\n", sep = "")
cat("  graded:", nrow(res_graded$metrics), " binary:", nrow(res_binary$metrics), "\n")
cat("\n=== LEVEL A — TRIAL LEVEL (n = 16 schools) ===\n")
print(seq_trial[, .(alphabet, metric,
                    ctrl = round(mean_ctrl, 3), exer = round(mean_exer, 3),
                    F = round(F, 2), df1, df2 = round(df2, 1),
                    p = round(p, 4), eta2_p = round(eta2_p, 3),
                    p_wilcox = round(p_wilcox, 4))])

cat("\n=== LEVEL B — TRIAL x TIMEPOINT (n = 48 sessions) ===\n")
print(seq_trial_tp[, .(alphabet, metric, term,
                       F = round(F, 2), df1, df2 = round(df2, 1),
                       p = round(p, 4), eta2_p = round(eta2_p, 3))])

cat("\n--- Level-B post-hoc (only for significant treatment x interval) ---\n")
if (nrow(seq_posthoc)) print(seq_posthoc) else cat("  none required\n")

cat("\n--- Singular fits (RE variance at boundary) ---\n")
sing_a <- seq_trial[isTRUE(singular) | singular %in% TRUE]
sing_b <- unique(seq_trial_tp[singular %in% TRUE, .(alphabet, metric, model)])
if (nrow(sing_a)) print(sing_a[, .(level = "trial", alphabet, metric, model)])
if (nrow(sing_b)) print(sing_b[, .(level = "trial_x_timepoint", alphabet, metric, model)])
if (!nrow(sing_a) && !nrow(sing_b)) cat("  none\n")

cat("\n--- Sparse transition rows (state never entered -> p_stay = NA) ---\n")
for (nm in c("graded", "binary")) {
  M <- if (nm == "graded") res_graded$metrics else res_binary$metrics
  lv <- if (nm == "graded") GRADED_LEVELS else BINARY_LEVELS
  for (s in lv) {
    col <- paste0("p_stay_", s)
    if (col %in% names(M)) {
      nmiss <- sum(!is.finite(M[[col]]))
      if (nmiss) cat(sprintf("  %s %-12s : %d/%d sessions undefined\n",
                             nm, col, nmiss, nrow(M)))
    }
  }
}
cat("\nFigures written to: ", fig_dir, "\n")
