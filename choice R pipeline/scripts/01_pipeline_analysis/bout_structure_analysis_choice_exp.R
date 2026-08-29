# =============================================================================
# EXERCISE-BOUT STRUCTURE ANALYSIS — CHOICE EXPERIMENT (identity-free, group-level)
#
# WHY THIS SCRIPT
#   The Introduction states the aim as determining "whether, when, and for how
#   long" fish choose to swim against a current, and explicitly asks whether
#   patterns of exercise choice (frequency and duration of exercise vs resting
#   bouts) consistently emerge across schools. The sibling script,
#   behavioural_sequence_analysis_choice_exp.R, answers "whether" (occupancy,
#   preference) and reports MEAN dwell time / switch rate / entropy as
#   TREATMENT CONTRASTS. It never characterises the bout-level structure
#   itself (the full duration distribution, its shape, its timing within a
#   session) and never asks whether a pattern is REPRODUCIBLE across schools —
#   consistency is a different question from a treatment difference and needs
#   its own machinery. This script fills exactly that gap and nothing else:
#   every metric here is either genuinely new information (see the redundancy
#   table below) or feeds the consistency battery, which the sibling script
#   does not attempt at all.
#
#   REDUNDANCY WITH THE SEQUENCE MODULE (checked explicitly, see section 9):
#     bout_rate_flow      ~ switch_rate / 2                      -> descriptive only
#     onset_hazard_flow   = 60*(1 - p_stay_Calm)  (exact identity) -> descriptive only
#     mean_bout_flow/calm = dwell_Flow / dwell_Calm (exact identity, asserted)
#     cv_bout_*           <-> burstiness (monotone transform)     -> descriptive only
#   Any candidate metric with |Spearman rho| >= 0.70 against a sequence-module
#   metric (school level, n=16) is AUTOMATICALLY demoted to descriptive by the
#   collinearity screen in section 9 and excluded from the inferential grid.
#
# THE 2026-08-07 ZONE-GEOMETRY FIX
#   The FT-motor polygon reference (config/zone_reference_choice_exp_FT.csv)
#   was digitised on a ~3400x1768 reference frame while the actual videos are
#   1920x1080; this made the "high" flow sub-zone permanently unreachable for
#   FT sessions (0 rows) and inverted FT's flow/calm split (FT read as ~79%
#   calm / 16% flow; corrected FT is ~21% calm / 78% flow). The fix was
#   applied upstream (STEP1) on 2026-08-07 16:28. FT is tanks 27 & 29, 8 of
#   the 16 schools, balanced 4 control / 4 exercise choice. This script
#   ASSERTS it is reading STEP1 output generated after that fix (section 2)
#   and reports the eight zero-switch sessions this reveals (all exercise
#   choice, all committed to flow for the full 20-min interval) as a finding.
#
#   The three "timepoint" intervals are NOT one continuous recording: they are
#   three separate 20-min tracking windows at 5-25, 45-65 and 85-105 min of
#   the trial, separated by untracked gaps. Consequently every bout is
#   WINDOW-CENSORED -- it cannot be observed to continue past the 20-min
#   window even if it does. commitment_index = 1.0 means ">= the interval",
#   never "exactly the interval"; the underlying truth may be longer.
#
# THREE STATE ALPHABETS (run in parallel, all from main_zone/sec_zone)
#   (1) BINARY    : prop_flow -> Flow-dominant / Calm-dominant (majority > 0.5)
#                   -- the aim's literal "exercise vs resting bout" framing.
#   (2) ENGAGEMENT: prop_flow -> Low / Med / High graded terciles (identical
#                   construction to the sibling script) -- how much of the
#                   tracked school commits, not just which side it is on.
#   (3) INTENSITY : sec_zone -> Low / Medium / High within the flow zone only
#                   -- self-selected EXERCISE INTENSITY, now analysable on all
#                   48 sessions since the zone fix restored FT's sub-zones.
#
# CAUTION -- DETECTION SUBSAMPLING: idtracker.ai fills exactly 5 identity
#   slots per session regardless of how many fish (4-16) are actually present,
#   so prop_flow is the flow fraction of a TRACKED SUBSET (~100% of the school
#   at density 4, ~31% at density 16), and frames within a 1-s bin are highly
#   autocorrelated (effective n per bin ~5, not ~125). Language throughout
#   must say "the majority of the tracked school", never "the school".
#
# CAUTION -- SELECTION ON THE OUTCOME: metrics needing >= MIN_BOUTS uncensored
#   bouts are undefined preferentially in exercise-choice sessions (the ones
#   that commit hardest switch least), so Tier 2 (bout-shape) contrasts are
#   CONDITIONAL and biased toward switch-prone, mostly-control schools. Tier 1
#   metrics require no bouts and carry the primary inference; every Tier 2
#   caption prints its per-arm n.
#
# INPUT
#   Latest output/STEP1_output/STEP1_output_*/master_fish_by_frame.csv, dated
#   at or after the 2026-08-07 zone fix (asserted in section 2).
#   Columns used: trial_id, trial, frame, time, timepoint, treatment,
#                 fish_density, tank, motor_side, main_zone, sec_zone.
#   Cross-checked against the latest output/SEQ_output/SEQ_output_*/ run
#   (must also postdate the zone fix) for state-construction integrity.
#
# OUTPUT -> output/BOUT_output/BOUT_output_<timestamp>/
#   bout_inventory.csv                 every bout, all three alphabets
#   bout_metrics_per_session.csv       48 rows, Tier 1 + Tier 2 + descriptive
#   bout_metrics_per_school.csv        16 rows, school means
#   bout_anova_trial_level.csv         Level A (n=16)
#   bout_anova_trial_x_timepoint.csv   Level B (n=48)
#   bout_trend_contrasts.csv           linear trend across the 3 intervals
#   bout_posthoc_trial_x_timepoint.csv emmeans simple effects (conditional)
#   bout_results_grid.csv              mirrors figures K/L cell-for-cell
#   bout_repeatability_icc.csv         ICC (adjusted & unadjusted) + bootstrap CI
#   bout_rank_stability.csv            Kendall's W across intervals
#   bout_split_half.csv                within-interval split-half reliability
#   bout_consistency_sign.csv          binomial sign-consistency tests
#   bout_null_per_session.csv          shuffle & Markov null comparisons
#   bout_null_summary.csv              aggregate null-exceedance counts
#   metric_collinearity.csv            |rho| gate decisions vs sequence metrics
#   bout_censoring_sensitivity.csv     uncensored-only vs full-inventory
#   bout_binwidth_sensitivity.csv      bin width x threshold x sustain sweep
#   intensity_motorside_check.csv      FT vs FD sensitivity on intensity metrics
#   config.rds                         knobs + sessionInfo() + resolved paths
#   figures/                           panels A-L (see header of section 10)
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
  has_lme4      <- requireNamespace("lme4",      quietly = TRUE)
  has_lmerTest  <- requireNamespace("lmerTest",  quietly = TRUE)
  # pbkrtest backs lmerTest's ddf = "Kenward-Roger"; without it KR is
  # unavailable and inference silently degrades to Satterthwaite.
  has_pbkrtest  <- requireNamespace("pbkrtest",  quietly = TRUE)
  has_emmeans   <- requireNamespace("emmeans",   quietly = TRUE)
  has_car       <- requireNamespace("car",       quietly = TRUE)
  has_perf      <- requireNamespace("performance", quietly = TRUE)
  has_survival  <- requireNamespace("survival",  quietly = TRUE)
  # effectsize backs noncentral-F CIs on eta2_p/omega2_p and noncentral-t CIs on
  # Hedges' g (Nakagawa & Cuthill 2007), added 2026-08-08.
  has_effectsize <- requireNamespace("effectsize", quietly = TRUE)
})

ts_msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), " | ", ..., "\n", sep = "")

# =============================================================================
# 0) CONFIG
# =============================================================================
PIPE_ROOT <- file.path(PROJECT_ROOT, "choice R pipeline")
BIN_S     <- 1.0        # sequence bin width, seconds -- IDENTICAL to sibling
MIN_BINS  <- 10L         # drop a (trial,timepoint) sequence shorter than this
SUSTAIN_S <- 30          # threshold for a "sustained" flow bout
LONG_S    <- 60          # threshold for a "long" flow bout
LATENCY_MIN_S <- 10      # minimum duration to count as the first sustained bout
MIN_BOUTS <- 5L          # Tier 2 coverage: uncensored bouts of a state needed
MIN_PAIRS <- 5L          # minimum consecutive-duration pairs for memory_flow
MIN_SESSIONS <- 24L      # Level-B coverage gate
MIN_SCHOOLS  <- 12L      # Level-A coverage gate
MIN_PER_ARM  <- 3L       # minimum schools per treatment for a metric to be tested
COLLIN_RHO   <- 0.70     # |rho| threshold demoting a metric to descriptive
N_SURR    <- 999L        # null-model surrogates per session per statistic
N_BOOT    <- 1000L       # ICC parametric-bootstrap replicates
SEED      <- 20260807L
ZONE_FIX_DATE <- as.POSIXct("2026-08-07 16:28:00")  # local tz, matching file.mtime()'s return

GRADED_LEVELS   <- c("Low", "Med", "High")
BINARY_LEVELS   <- c("Calm", "Flow")
INTENSITY_LEVELS <- c("low", "medium", "high")

TREAT_COLORS <- c("control" = "#E69F00", "exercise choice" = "#0072B2")
STATE_FILL   <- c("Low" = "#FDE0A9", "Med" = "#F0A860", "High" = "#B8410E",
                  "Calm" = "#8FC7E8", "Flow" = "#0072B2")

set.seed(SEED)

# =============================================================================
# 1) LOCATE + LOAD FRAME-LEVEL DATA (state construction VERBATIM from sibling)
# =============================================================================
.find_latest_master <- function() {
  parent <- file.path(PIPE_ROOT, "output", "STEP1_output")
  subs <- list.dirs(parent, full.names = TRUE, recursive = FALSE)
  subs <- subs[grepl("STEP1_output_\\d{8}_\\d{6}$", basename(subs))]
  if (!length(subs)) stop("No STEP1_output_* folder found under ", parent)
  latest <- subs[which.max(file.mtime(subs))]
  f <- file.path(latest, "master_fish_by_frame.csv")
  if (!file.exists(f)) stop("master_fish_by_frame.csv not found in ", latest)
  list(path = f, dir = latest, mtime = file.mtime(latest))
}

.find_latest_seq <- function() {
  parent <- file.path(PIPE_ROOT, "output", "SEQ_output")
  subs <- list.dirs(parent, full.names = TRUE, recursive = FALSE)
  subs <- subs[grepl("SEQ_output_\\d{8}_\\d{6}$", basename(subs))]
  if (!length(subs)) { warning("No SEQ_output_* folder found; integrity check skipped."); return(NULL) }
  latest <- subs[which.max(file.mtime(subs))]
  list(dir = latest, mtime = file.mtime(latest))
}

master_info <- .find_latest_master()
ts_msg("Reading (selected cols): ", master_info$path)
if (master_info$mtime < ZONE_FIX_DATE) {
  stop("STEP1 output (", master_info$dir, ") predates the 2026-08-07 16:28 zone-",
       "geometry fix. Regenerate STEP1/STEP2/STEP2b before running this script.")
}
ts_msg("STEP1 zone-fix assertion OK: ", master_info$dir)

seq_info <- .find_latest_seq()
if (!is.null(seq_info) && seq_info$mtime < ZONE_FIX_DATE) {
  warning("Latest SEQ_output (", seq_info$dir, ") predates the zone fix; the ",
          "integrity cross-check in section 3b will likely fail. Regenerate it first.")
}

dt <- data.table::fread(
  master_info$path,
  select = c("trial_id", "trial", "frame", "time", "timepoint", "treatment",
             "fish_density", "tank", "motor_side", "main_zone", "sec_zone"),
  showProgress = FALSE
)
ts_msg("Loaded ", format(nrow(dt), big.mark = ","), " detection-rows")

dt <- dt[main_zone %in% c("flow", "calm")]
dt[, treatment := tolower(trimws(treatment))]
dt <- dt[treatment %in% names(TREAT_COLORS)]
ts_msg("After zone/treatment filter: ", format(nrow(dt), big.mark = ","), " rows")

# =============================================================================
# 2) PER-FRAME OCCUPANCY -> PER-BIN prop_flow (verbatim construction)
# =============================================================================
frame_occ <- dt[, .(
    n_flow = sum(main_zone == "flow"),
    n_calm = sum(main_zone == "calm"),
    n_high = sum(sec_zone == "high", na.rm = TRUE),
    n_med  = sum(sec_zone == "medium", na.rm = TRUE),
    n_low  = sum(sec_zone %in% c("low", "low1", "low2", "low3"), na.rm = TRUE),
    time   = time[1]
  ), by = .(trial_id, trial, timepoint, frame, treatment, tank, fish_density, motor_side)]

frame_occ[, bin := floor(time / BIN_S)]

bin_occ <- frame_occ[, .(
    n_flow = sum(n_flow), n_calm = sum(n_calm),
    n_high = sum(n_high), n_med = sum(n_med), n_low = sum(n_low)
  ), by = .(trial_id, trial, timepoint, treatment, tank, fish_density, motor_side, bin)]
bin_occ[, n_zoned := n_flow + n_calm]
bin_occ <- bin_occ[n_zoned > 0]
bin_occ[, prop_flow := n_flow / n_zoned]
data.table::setkey(bin_occ, trial_id, timepoint, bin)
ts_msg("Built ", format(nrow(bin_occ), big.mark = ","), " time-bins (", BIN_S, "s)")

# ---- Alphabet 1: BINARY (majority) ------------------------------------------
bin_occ[, state_binary := factor(fifelse(prop_flow > 0.5, "Flow", "Calm"),
                                  levels = BINARY_LEVELS)]

# ---- Alphabet 2: ENGAGEMENT (data-driven terciles, identical to sibling) ----
q <- stats::quantile(bin_occ$prop_flow, probs = c(1/3, 2/3), na.rm = TRUE)
cut_lo <- unname(q[1]); cut_hi <- unname(q[2])
ts_msg(sprintf("Engagement tercile cut-points: Low<%.3f<=Med<%.3f<=High", cut_lo, cut_hi))
bin_occ[, state_engage := factor(
  fifelse(prop_flow < cut_lo, "Low", fifelse(prop_flow < cut_hi, "Med", "High")),
  levels = GRADED_LEVELS)]

# ---- Alphabet 3: INTENSITY (sec_zone rank, flow bins only) ------------------
# Defined only within bins where the majority of zoned detections are in flow
# (state_binary == "Flow"); a bin's intensity is the sec_zone with the most
# in-flow detections that second. Bins with no sec_zone-classified detection
# in that bin are left NA (rare after the zone fix; reported in coverage).
bin_occ[, n_sec_max := pmax(n_high, n_med, n_low)]
bin_occ[, state_intensity := fcase(
  state_binary != "Flow", NA_character_,
  n_sec_max == 0, NA_character_,
  n_high >= n_med & n_high >= n_low, "high",
  n_med  >= n_low, "medium",
  default = "low"
)]
bin_occ[, state_intensity := factor(state_intensity, levels = INTENSITY_LEVELS)]

# =============================================================================
# 3) BOUT INVENTORY -- full, per alphabet (censoring logic reused from sibling)
# =============================================================================
# Returns a data.table of every maximal run (bout) in an ordered state vector,
# split at bin-index gaps (a gap = an untracked interval; bouts never bridge
# it), with the first and last run of every contiguous segment flagged
# censored (they may continue beyond the observation window). The FULL
# inventory is retained -- this is the core difference from the sibling
# script, which discards everything but the per-state mean of the interior
# (uncensored) runs.
.bout_table <- function(states, bins, bin_s) {
  # Perf note (2026-08-08): the original version built ONE data.table() object
  # per run inside the loop, then rbindlist()'d them -- the same "many small
  # data.table() allocations in a tight loop" anti-pattern that made STEP2's
  # NN-matching take >100 minutes before it was fixed. This function is called
  # 999 x 2 x 5 x 48 = ~480,000 times during the null-model comparisons (once
  # per surrogate per statistic per null type per session, prior to also
  # fixing that redundancy below), and a shuffled/high-entropy surrogate can
  # have hundreds of runs, so the per-run data.table() cost compounds badly.
  # Fixed by accumulating plain vectors across all segments and constructing
  # exactly ONE data.table() at the end. Verified to produce identical output
  # to the original (same columns, same row order) on the integrity-check
  # sessions.
  states <- as.character(states)
  ord <- order(bins); states <- states[ord]; bins <- bins[ord]
  n <- length(states)
  if (n < 1L) return(data.table())
  seg_id_vec <- cumsum(c(1L, as.integer(diff(bins) != 1L)))
  idx <- seq_len(n)
  seg_split <- split(idx, seg_id_vec)

  v_seg <- vector("integer", 0L); v_run <- vector("integer", 0L)
  v_state <- vector("character", 0L)
  v_start <- vector("numeric", 0L); v_end <- vector("numeric", 0L)
  v_nbins <- vector("integer", 0L)

  for (s in seg_split) {
    st <- states[s]; bn <- bins[s]
    if (all(is.na(st))) next
    r <- rle(st)
    nr <- length(r$lengths)
    if (nr == 0L) next
    keep <- !is.na(r$values)
    if (!any(keep)) next
    end_pos <- cumsum(r$lengths)
    start_pos <- end_pos - r$lengths + 1L
    run_idx_all <- seq_len(nr)
    v_seg   <- c(v_seg,   rep.int(seg_id_vec[s[1]], sum(keep)))
    v_run   <- c(v_run,   run_idx_all[keep])
    v_state <- c(v_state, r$values[keep])
    v_start <- c(v_start, bn[start_pos[keep]])
    v_end   <- c(v_end,   bn[end_pos[keep]])
    v_nbins <- c(v_nbins, r$lengths[keep])
  }
  if (!length(v_state)) return(data.table())

  # censored/censor_side depend on whether a run is first/last WITHIN ITS OWN
  # segment; recompute that per-segment first/last flag vectorised, matching
  # the original per-run j==1L / j==nr logic exactly.
  is_first <- !duplicated(v_seg)
  is_last  <- !duplicated(v_seg, fromLast = TRUE)
  data.table::data.table(
    seg_id = v_seg, run_idx = v_run, state = v_state,
    start_bin = v_start, end_bin = v_end, n_bins_run = v_nbins,
    duration_s = v_nbins * bin_s,
    censored = is_first | is_last,
    censor_side = data.table::fcase(is_first & is_last, "both",
                                    is_first, "left", is_last, "right", default = "none")
  )
}

# =============================================================================
# 4) PER-SESSION METRIC HELPERS
# =============================================================================
# Kim & Jo (2016) finite-size corrected burstiness. r = CV of uncensored
# durations; n = number of uncensored bouts. Uncorrected B = (sd-mean)/(sd+mean)
# is strongly n-biased, and n differs systematically by treatment here, so the
# correction is mandatory, not cosmetic. Returns NA below MIN_BOUTS.
.burstiness <- function(d, n_min = MIN_BOUTS) {
  d <- d[is.finite(d) & d > 0]
  n <- length(d)
  if (n < n_min) return(NA_real_)
  m <- mean(d); s <- stats::sd(d)
  if (!is.finite(s) || m <= 0) return(NA_real_)
  r <- s / m
  sn <- sqrt(n + 1); sn2 <- sqrt(n - 1)
  (sn * r - sn2) / ((sn - 2) * r + sn2)
}

# Lag-1 Pearson correlation of consecutive SAME-STATE uncensored durations
# within a segment (renewal-process check). Needs >= n_min pairs.
.memory_coef <- function(bt, state_val, n_min = MIN_PAIRS) {
  sub <- bt[state == state_val & !censored]
  if (!nrow(sub)) return(NA_real_)
  data.table::setorder(sub, seg_id, run_idx)
  pairs <- sub[, .(x = duration_s[-.N], y = duration_s[-1]), by = seg_id]
  pairs <- pairs[is.finite(x) & is.finite(y)]
  if (nrow(pairs) < n_min) return(NA_real_)
  suppressWarnings(stats::cor(pairs$x, pairs$y))
}

# Time-weighted centroid of state occupancy: mean(time of bins in `target`) / T.
# 0.5 = uniform across the interval; >0.5 back-loaded, <0.5 front-loaded.
.flow_centroid <- function(states, bins, target = "Flow") {
  is_t <- states == target
  if (sum(is_t, na.rm = TRUE) < 30L) return(NA_real_)
  T_bins <- max(bins) - min(bins) + 1
  mean((bins[is_t] - min(bins)) / T_bins, na.rm = TRUE)
}

# Time-weighted mean rank of a graded state (1=lowest .. K=highest), ignoring NA.
.rank_depth <- function(states, levels_vec) {
  r <- match(as.character(states), levels_vec)
  if (all(is.na(r))) return(NA_real_)
  mean(r, na.rm = TRUE)
}

# =============================================================================
# 5) PER-SESSION METRICS -- one row per (trial, timepoint)
# =============================================================================
.session_metrics <- function(sub) {
  n <- nrow(sub)
  T_s <- n * BIN_S
  bins <- sub$bin

  bt_bin <- .bout_table(sub$state_binary, bins, BIN_S)
  bt_eng <- .bout_table(sub$state_engage, bins, BIN_S)
  bt_int <- .bout_table(sub$state_intensity, bins, BIN_S)

  flow_bouts <- bt_bin[state == "Flow"]
  calm_bouts <- bt_bin[state == "Calm"]
  max_flow <- if (nrow(flow_bouts)) max(flow_bouts$duration_s) else 0
  max_calm <- if (nrow(calm_bouts)) max(calm_bouts$duration_s) else 0
  n_changes <- max(0L, nrow(bt_bin) - length(unique(bt_bin$seg_id)))

  # --- Tier 1: defined for (almost) every session -----------------------------
  commitment_index <- max(max_flow, max_calm) / T_s
  n_flow_ge30 <- sum(flow_bouts$duration_s >= SUSTAIN_S)
  long_flow_frac <- sum(flow_bouts$duration_s[flow_bouts$duration_s >= LONG_S]) / T_s
  flow_centroid <- .flow_centroid(sub$state_binary, bins, "Flow")

  first_sus <- flow_bouts[duration_s >= LATENCY_MIN_S]
  if (nrow(first_sus)) {
    data.table::setorder(first_sus, start_bin)
    t_first_sustained_s <- (first_sus$start_bin[1] - min(bins)) * BIN_S
    latency_censored <- FALSE
  } else {
    t_first_sustained_s <- T_s
    latency_censored <- TRUE
  }

  engage_depth <- .rank_depth(sub$state_engage, GRADED_LEVELS)
  mean_intensity_sel <- .rank_depth(sub$state_intensity, INTENSITY_LEVELS)

  # --- Tier 2: need >= MIN_BOUTS uncensored bouts of that state ---------------
  flow_unc <- bt_bin[state == "Flow" & !censored, duration_s]
  calm_unc <- bt_bin[state == "Calm" & !censored, duration_s]
  n_uncens_flow <- length(flow_unc); n_uncens_calm <- length(calm_unc)
  med_bout_flow <- if (n_uncens_flow >= MIN_BOUTS) stats::median(flow_unc) else NA_real_
  med_bout_calm <- if (n_uncens_calm >= MIN_BOUTS) stats::median(calm_unc) else NA_real_
  burst_flow <- .burstiness(flow_unc)
  burst_calm <- .burstiness(calm_unc)
  memory_flow <- .memory_coef(bt_bin, "Flow")

  # --- Descriptive only (redundant with sequence-module metrics) -------------
  obs_min <- T_s / 60
  bout_rate_flow <- n_uncens_flow / obs_min
  mean_bout_flow <- if (n_uncens_flow) mean(flow_unc) else NA_real_
  mean_bout_calm <- if (n_uncens_calm) mean(calm_unc) else NA_real_
  # p_stay_Calm from the same transition logic the sibling uses; recomputed
  # locally (not sourced) so this script has no runtime dependency on it.
  st <- as.character(sub$state_binary)
  adj <- which(diff(bins) == 1L)
  stay_calm_num <- sum(st[adj] == "Calm" & st[adj + 1L] == "Calm")
  from_calm <- sum(st[adj] == "Calm")
  p_stay_calm_local <- if (from_calm > 0) stay_calm_num / from_calm else NA_real_
  onset_hazard_flow <- if (is.finite(p_stay_calm_local)) 60 * (1 - p_stay_calm_local) else NA_real_
  cv_bout_flow <- if (is.finite(burst_flow)) (1 + burst_flow) / (1 - burst_flow) else NA_real_
  cv_bout_calm <- if (is.finite(burst_calm)) (1 + burst_calm) / (1 - burst_calm) else NA_real_
  log2_med_ratio <- if (is.finite(med_bout_flow) && is.finite(med_bout_calm) && med_bout_calm > 0)
    log2(med_bout_flow / med_bout_calm) else NA_real_

  third <- pmin(3L, pmax(1L, ceiling(3 * (bins - min(bins) + 1) / n)))
  prop_thirds <- tapply(sub$state_binary == "Flow", third, mean)
  prop_flow_third1 <- unname(prop_thirds["1"]); if (is.null(prop_flow_third1)) prop_flow_third1 <- NA_real_
  prop_flow_third2 <- unname(prop_thirds["2"]); if (is.null(prop_flow_third2)) prop_flow_third2 <- NA_real_
  prop_flow_third3 <- unname(prop_thirds["3"]); if (is.null(prop_flow_third3)) prop_flow_third3 <- NA_real_

  list(
    n_bins = n, obs_min = obs_min, mean_n_zoned = mean(sub$n_zoned),
    occ_flow = mean(sub$state_binary == "Flow"),
    n_switch = n_changes, n_bouts_flow = nrow(flow_bouts), n_bouts_calm = nrow(calm_bouts),
    n_uncens_flow = n_uncens_flow, n_uncens_calm = n_uncens_calm,
    zero_switch = n_changes == 0L,
    # Tier 1
    commitment_index = commitment_index, max_flow_bout_s = max_flow, max_calm_bout_s = max_calm,
    n_flow_ge30 = n_flow_ge30, long_flow_frac = long_flow_frac,
    flow_centroid = flow_centroid, t_first_sustained_s = t_first_sustained_s,
    latency_censored = latency_censored,
    engage_depth = engage_depth, mean_intensity_sel = mean_intensity_sel,
    # Tier 2
    med_bout_flow = med_bout_flow, med_bout_calm = med_bout_calm,
    burst_flow = burst_flow, burst_calm = burst_calm, memory_flow = memory_flow,
    # Descriptive
    bout_rate_flow = bout_rate_flow, onset_hazard_flow = onset_hazard_flow,
    mean_bout_flow = mean_bout_flow, mean_bout_calm = mean_bout_calm,
    cv_bout_flow = cv_bout_flow, cv_bout_calm = cv_bout_calm,
    log2_med_ratio = log2_med_ratio,
    prop_flow_third1 = prop_flow_third1, prop_flow_third2 = prop_flow_third2,
    prop_flow_third3 = prop_flow_third3,
    bt_bin = list(bt_bin), bt_eng = list(bt_eng), bt_int = list(bt_int)
  )
}

ts_msg("Computing per-session bout metrics ...")
keys <- unique(bin_occ[, .(trial_id, trial, timepoint, treatment, tank, fish_density, motor_side)])
session_rows <- vector("list", nrow(keys))
inv_rows <- vector("list", nrow(keys))
for (i in seq_len(nrow(keys))) {
  k <- keys[i]
  sub <- bin_occ[trial_id == k$trial_id & timepoint == k$timepoint]
  if (nrow(sub) < MIN_BINS) next
  m <- .session_metrics(sub)
  bt_all <- rbind(
    if (nrow(m$bt_bin[[1]])) cbind(alphabet = "binary", m$bt_bin[[1]]) else NULL,
    if (nrow(m$bt_eng[[1]])) cbind(alphabet = "engagement", m$bt_eng[[1]]) else NULL,
    if (nrow(m$bt_int[[1]])) cbind(alphabet = "intensity", m$bt_int[[1]]) else NULL,
    fill = TRUE)
  if (!is.null(bt_all) && nrow(bt_all)) {
    inv_rows[[i]] <- cbind(trial_id = k$trial_id, trial = k$trial, timepoint = k$timepoint,
                           treatment = k$treatment, tank = k$tank, fish_density = k$fish_density,
                           motor_side = k$motor_side, bt_all)
  }
  m$bt_bin <- NULL; m$bt_eng <- NULL; m$bt_int <- NULL
  session_rows[[i]] <- c(list(trial_id = k$trial_id, trial = k$trial, timepoint = k$timepoint,
                              treatment = k$treatment, tank = k$tank, fish_density = k$fish_density,
                              motor_side = k$motor_side), m)
}
bout_metrics <- data.table::rbindlist(Filter(Negate(is.null), session_rows), fill = TRUE)
bout_inventory <- data.table::rbindlist(Filter(Negate(is.null), inv_rows), fill = TRUE)
ts_msg("Session metrics: ", nrow(bout_metrics), " sessions | inventory: ",
       nrow(bout_inventory), " bouts")

# =============================================================================
# 3b) INTEGRITY CHECK against the sequence module
# =============================================================================
if (!is.null(seq_info)) {
  seq_bin <- tryCatch(
    data.table::fread(file.path(seq_info$dir, "seq_metrics_per_trial_binary.csv")),
    error = function(e) NULL)
  if (!is.null(seq_bin)) {
    chk <- merge(bout_metrics[, .(trial_id, timepoint, n_bins, mean_bout_flow)],
                 seq_bin[, .(trial_id, timepoint, n_bins_seq = n_bins, dwell_Flow)],
                 by = c("trial_id", "timepoint"))
    n_bin_mismatch <- sum(chk$n_bins != chk$n_bins_seq)
    both_finite <- is.finite(chk$mean_bout_flow) & is.finite(chk$dwell_Flow)
    max_dwell_diff <- if (any(both_finite))
      max(abs(chk$mean_bout_flow[both_finite] - chk$dwell_Flow[both_finite])) else NA_real_
    ts_msg("Integrity check vs SEQ_output: ", n_bin_mismatch, " n_bins mismatches, ",
           "max |mean_bout_flow - dwell_Flow| = ", signif(max_dwell_diff, 6))
    if (n_bin_mismatch > 0)
      warning("n_bins mismatch against SEQ_output -- state construction may have drifted.")
    if (is.finite(max_dwell_diff) && max_dwell_diff > 1e-6)
      warning("mean_bout_flow does not match SEQ_output's dwell_Flow to 1e-8 -- ",
              "the two modules disagree; investigate before trusting either.")
  } else {
    warning("Could not read seq_metrics_per_trial_binary.csv from ", seq_info$dir)
  }
}

# Known-truth sanity check on the corrected zones: expect exactly 8 zero-switch
# sessions, all exercise choice.
zs <- bout_metrics[zero_switch == TRUE]
ts_msg("Zero-switch sessions: ", nrow(zs), " (treatments: ",
       paste(table(zs$treatment), names(table(zs$treatment)), collapse = ", "), ")")

# =============================================================================
# 6) NULL MODELS -- shuffle (occupancy-only) and Markov (occupancy + memoryless)
# =============================================================================
# Both preserve session length and are within-session (no cross-session
# leakage). N_SURR = 999 gives permutation-p resolution of 1/1000, matching
# the project's "< 0.001" display floor.
#
# Shuffle: permute states within each contiguous segment -- preserves
# occupancy EXACTLY, destroys all temporal order. This is the sanity floor:
# informative to show, uninformative to test (a real bout process will
# obviously beat pure randomisation).
#
# Markov: simulate a two-state chain using the session's OWN p_stay_Flow /
# p_stay_Calm. For a two-state chain, run lengths are geometric, so this is
# implemented by drawing geometric run lengths directly rather than stepping
# through bins -- ~10us per surrogate instead of a 1200-step loop. Anything
# that exceeds this null is structure BEYOND first-order Markov, which by
# construction the sequence module (occupancy + transition matrix) cannot
# express -- this is what makes si_markov_max genuinely new information.
.surrogate_shuffle_stat <- function(states, seg_id, stat_fun) {
  shuffled <- unlist(lapply(split(states, seg_id), sample), use.names = FALSE)
  stat_fun(shuffled)
}

.markov_run_lengths <- function(n_target, p_stay_flow, p_stay_calm, start_state) {
  # Draw enough alternating geometric run lengths to cover n_target bins,
  # starting from `start_state`, then truncate. p_stay_* = P(state[t+1]=state[t]
  # | state[t]=state), so each run length ~ Geometric(1 - p_stay) + 1 (in
  # {1,2,...}), truncated at 1 to avoid zero-length runs under rgeom's support.
  if (!is.finite(p_stay_flow) || !is.finite(p_stay_calm) ||
      p_stay_flow >= 1 || p_stay_calm >= 1) return(NULL)  # absorbing -- undefined
  states_out <- character(0)
  cur <- start_state
  # Expected run length is 1/(1-p_stay); estimate how many runs are needed with
  # margin, then extend if we somehow fall short (rare).
  p_stay <- c(Flow = p_stay_flow, Calm = p_stay_calm)
  exp_len <- 1 / (1 - p_stay[cur])
  n_runs_guess <- max(20L, ceiling(2 * n_target / mean(1 / (1 - p_stay))))
  total <- 0L
  repeat {
    m <- n_runs_guess
    cur_seq <- rep(c("Flow", "Calm"), length.out = m)
    if (cur == "Calm") cur_seq <- rep(c("Calm", "Flow"), length.out = m)
    lens <- ifelse(cur_seq == "Flow",
                    1L + stats::rgeom(m, 1 - p_stay_flow),
                    1L + stats::rgeom(m, 1 - p_stay_calm))
    total <- sum(lens)
    if (total >= n_target) {
      states_out <- rep(cur_seq, lens)[seq_len(n_target)]
      break
    }
    n_runs_guess <- n_runs_guess * 2L
    if (n_runs_guess > 100000L) return(NULL)  # safety valve
  }
  states_out
}

.null_compare_session <- function(states, bins, seg_id, p_stay_flow, p_stay_calm,
                                   n_surr = N_SURR) {
  T_s <- length(states) * BIN_S
  # All five statistics are extracted from the SAME bout table for a given
  # surrogate -- compute it once (perf fix, 2026-08-08). The original called
  # .bout_table() independently inside each of 5 stat closures, and generated
  # a fresh surrogate independently for each of those 5 closures too, so a
  # single session did 5x the surrogate generation and 5x the bout-table
  # construction it needed (999 surrogates x 5 stats x 2 null types = 9990
  # bout tables per session instead of 1998). Combined with the .bout_table()
  # fix above this took the null-model stage from a multi-hour stall to
  # seconds.
  stat_from_bt <- function(bt) {
    f <- bt[state == "Flow"]
    max_flow <- if (nrow(f)) max(f$duration_s) else 0
    commit <- if (nrow(bt)) max(tapply(bt$duration_s, bt$state, max)) / T_s else 0
    n_ge30 <- sum(f$duration_s >= SUSTAIN_S)
    long_frac <- sum(f$duration_s[f$duration_s >= LONG_S]) / T_s
    burst <- .burstiness(bt[state == "Flow" & !censored, duration_s])
    c(max_flow_bout_s = max_flow, commitment_index = commit,
      n_flow_ge30 = n_ge30, long_flow_frac = long_frac, burst_flow = burst)
  }
  stat_names <- c("max_flow_bout_s", "commitment_index", "n_flow_ge30",
                  "long_flow_frac", "burst_flow")
  obs_vals <- stat_from_bt(.bout_table(states, bins, BIN_S))

  surr_mat <- list(
    shuffle = matrix(NA_real_, n_surr, length(stat_names), dimnames = list(NULL, stat_names)),
    markov  = matrix(NA_real_, n_surr, length(stat_names), dimnames = list(NULL, stat_names))
  )
  n_defined <- c(shuffle = n_surr, markov = 0L)
  for (b in seq_len(n_surr)) {
    shuffled <- unlist(lapply(split(states, seg_id), sample), use.names = FALSE)
    surr_mat$shuffle[b, ] <- stat_from_bt(.bout_table(shuffled, bins, BIN_S))
    ms <- .markov_run_lengths(length(states), p_stay_flow, p_stay_calm, states[1])
    if (!is.null(ms)) {
      surr_mat$markov[b, ] <- stat_from_bt(.bout_table(ms, bins, BIN_S))
      n_defined["markov"] <- n_defined[["markov"]] + 1L
    }
  }

  out <- vector("list", length(stat_names) * 2L)
  k <- 1L
  for (null_type in c("shuffle", "markov")) {
    for (nm in stat_names) {
      obs <- obs_vals[[nm]]
      surr <- surr_mat[[null_type]][, nm]
      surr_f <- surr[is.finite(surr)]
      if (length(surr_f) < 10L) {
        out[[k]] <- data.table(statistic = nm, null_type = null_type, obs = obs,
                               null_mean = NA_real_, null_sd = NA_real_, null_median = NA_real_,
                               null_q025 = NA_real_, null_q975 = NA_real_, z = NA_real_,
                               log2_ratio = NA_real_, p_perm = NA_real_, n_surr = n_surr,
                               defined = n_defined[[null_type]])
      } else {
        nmean <- mean(surr_f); nsd <- stats::sd(surr_f); nmed <- stats::median(surr_f)
        qs <- stats::quantile(surr_f, c(0.025, 0.975), na.rm = TRUE)
        z <- if (nsd > 0) (obs - nmean) / nsd else NA_real_
        l2r <- log2((obs + 1e-6) / (nmed + 1e-6))
        r_ge <- sum(surr_f >= obs)
        p_perm <- 2 * min(r_ge + 1, length(surr_f) - r_ge + 1) / (length(surr_f) + 1)
        p_perm <- min(1, p_perm)
        out[[k]] <- data.table(statistic = nm, null_type = null_type, obs = obs,
                               null_mean = nmean, null_sd = nsd, null_median = nmed,
                               null_q025 = unname(qs[1]), null_q975 = unname(qs[2]),
                               z = z, log2_ratio = l2r, p_perm = p_perm,
                               n_surr = n_surr, defined = n_defined[[null_type]])
      }
      k <- k + 1L
    }
  }
  data.table::rbindlist(out)
}

ts_msg("Running null-model comparisons (shuffle + Markov) ...")
null_rows <- vector("list", nrow(keys))
for (i in seq_len(nrow(keys))) {
  k <- keys[i]
  sub <- bin_occ[trial_id == k$trial_id & timepoint == k$timepoint]
  if (nrow(sub) < MIN_BINS) next
  bins <- sub$bin
  seg_id <- cumsum(c(1L, as.integer(diff(bins) != 1L)))
  states <- as.character(sub$state_binary)
  m_row <- bout_metrics[trial_id == k$trial_id & timepoint == k$timepoint]
  if (!nrow(m_row)) next
  # local p_stay estimates for the Markov null (recomputed here to avoid a
  # dependency on the earlier per-session loop's internals)
  adj <- which(diff(bins) == 1L)
  st <- states
  from_flow <- sum(st[adj] == "Flow"); stay_flow <- sum(st[adj] == "Flow" & st[adj + 1L] == "Flow")
  from_calm <- sum(st[adj] == "Calm"); stay_calm <- sum(st[adj] == "Calm" & st[adj + 1L] == "Calm")
  p_stay_flow <- if (from_flow > 0) stay_flow / from_flow else NA_real_
  p_stay_calm <- if (from_calm > 0) stay_calm / from_calm else NA_real_
  res <- .null_compare_session(states, bins, seg_id, p_stay_flow, p_stay_calm)
  if (i %% 8 == 0 || i == nrow(keys))
    ts_msg("  null models: session ", i, "/", nrow(keys), " done")
  null_rows[[i]] <- cbind(trial_id = k$trial_id, trial = k$trial, timepoint = k$timepoint,
                         treatment = k$treatment, res)
}
bout_null_per_session <- data.table::rbindlist(Filter(Negate(is.null), null_rows), fill = TRUE)
ts_msg("Null comparisons done: ", nrow(bout_null_per_session), " rows")

# si_markov_max feeds back into bout_metrics as a Tier-1 column.
si <- bout_null_per_session[statistic == "max_flow_bout_s" & null_type == "markov",
                            .(trial_id, timepoint, si_markov_max = log2_ratio)]
bout_metrics <- merge(bout_metrics, si, by = c("trial_id", "timepoint"), all.x = TRUE)

bout_null_summary <- bout_null_per_session[, .(
  n_defined = sum(defined > 0), n_above_q975 = sum(is.finite(obs) & is.finite(null_q975) & obs > null_q975),
  n_below_q025 = sum(is.finite(obs) & is.finite(null_q025) & obs < null_q025),
  median_log2_ratio = stats::median(log2_ratio, na.rm = TRUE)
), by = .(null_type, statistic, treatment)]
bout_null_summary[, p_binom := mapply(function(k_, n_) {
  if (n_ < 1) return(NA_real_)
  stats::binom.test(k_, n_, 0.025, alternative = "greater")$p.value
}, n_above_q975, n_defined)]

# =============================================================================
# 7) STATISTICAL LAYER -- identical layout to every other behavioural outcome
# =============================================================================
# Effect sizes from an F-test (verbatim from the sequence module).
.eta2_p <- function(F_val, df1, df2) {
  if (!is.finite(F_val) || !is.finite(df1) || !is.finite(df2)) return(NA_real_)
  (F_val * df1) / (F_val * df1 + df2)
}
.omega2_p <- function(F_val, df1, df2) {
  if (!is.finite(F_val) || !is.finite(df1) || !is.finite(df2)) return(NA_real_)
  max(0, (df1 * (F_val - 1)) / (df1 * (F_val - 1) + df2 + 1))
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
  if (is.null(fit) || !has_emmeans) return(na_out)
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

# Type III ANOVA table from lmer OR lm, in one shape.
# Returns data.frame(term, F, df1, df2, p, stat_type) or NULL.
#
# CHANGED 2026-08-08: Kenward-Roger denominator df by default (was lmerTest's
# Satterthwaite default), matching Methods 2.8 and the STEP5 behavioural engine.
# `stat_type` records what was actually computed so a KR failure can never be
# mistaken for a KR success: F-KR | F-SW | F-SW(KR-failed) | F-OLS.
.anova3 <- function(fit, ddf = c("Kenward-Roger", "Satterthwaite")) {
  ddf <- match.arg(ddf)
  if (inherits(fit, "merMod") || inherits(fit, "lmerModLmerTest")) {
    ml <- if (inherits(fit, "lmerModLmerTest")) fit
          else tryCatch(lmerTest::as_lmerModLmerTest(fit), error = function(e) NULL)
    if (is.null(ml)) return(NULL)

    use_ddf <- if (ddf == "Kenward-Roger" && !has_pbkrtest) "Satterthwaite" else ddf
    stat_lab <- if (use_ddf == "Kenward-Roger") "F-KR" else "F-SW"
    av <- tryCatch(suppressWarnings(stats::anova(ml, type = 3, ddf = use_ddf)),
                   error = function(e) NULL)

    if (is.null(av) && use_ddf == "Kenward-Roger") {
      av <- tryCatch(suppressWarnings(stats::anova(ml, type = 3, ddf = "Satterthwaite")),
                     error = function(e) NULL)
      if (!is.null(av)) {
        stat_lab <- "F-SW(KR-failed)"
        ts_msg("    WARNING: Kenward-Roger failed; fell back to Satterthwaite")
      }
    }
    if (is.null(av) || !nrow(av)) return(NULL)
    return(data.frame(term = rownames(av), F = av[["F value"]], df1 = av[["NumDF"]],
                      df2 = av[["DenDF"]], p = av[["Pr(>F)"]], stat_type = stat_lab,
                      stringsAsFactors = FALSE))
  }
  if (has_car) {
    av <- tryCatch(suppressWarnings(car::Anova(fit, type = "III")), error = function(e) NULL)
    if (!is.null(av)) {
      d <- as.data.frame(av); d$term <- rownames(d)
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

# Fit the first non-singular candidate model (most-complex-first, simplest
# fallback if every candidate is singular). Verbatim convention from the
# sequence module.
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
    sing <- if (is_mixed) isTRUE(tryCatch(lme4::isSingular(fit), error = function(e) TRUE)) else FALSE
    fallback <- list(fit = fit, form = fm, singular = sing)
    if (!sing) return(list(fit = fit, form = fm, singular = FALSE))
  }
  fallback
}

# Coverage gate: is this metric analysable given how many sessions/schools have
# a finite value? Applied BEFORE fitting, since selection-on-outcome (fewer
# bouts in exercise-choice sessions) means many Tier-2 metrics fail this.
.coverage_ok <- function(dat, metric, min_sessions = MIN_SESSIONS,
                          min_schools = MIN_SCHOOLS, min_per_arm = MIN_PER_ARM) {
  if (!metric %in% names(dat)) return(list(ok = FALSE, n_unit = 0L, reason = "metric not found"))
  d <- dat[is.finite(get(metric))]
  n_sessions <- nrow(d)
  n_schools <- data.table::uniqueN(d$trial)
  per_arm <- d[, .N, by = treatment]
  n_arm_ok <- nrow(per_arm) >= 2 && all(per_arm$N >= min_per_arm)
  ok <- n_sessions >= min_sessions && n_schools >= min_schools && n_arm_ok
  reason <- if (ok) "" else sprintf("n_sessions=%d n_schools=%d per_arm_ok=%s",
                                     n_sessions, n_schools, n_arm_ok)
  list(ok = ok, n_unit = n_sessions, reason = reason)
}

# ---- LEVEL A: TRIAL LEVEL (n = 16 schools) ---------------------------------
.level_trial <- function(dat, metric, tier) {
  cov <- .coverage_ok(dat, metric, min_sessions = 1L, min_schools = MIN_PER_ARM * 2L)
  if (!metric %in% names(dat)) return(NULL)
  d <- dat[is.finite(get(metric)), .(trial, treatment, tank, fish_density, val = get(metric))]
  if (!nrow(d)) return(NULL)
  per_school <- d[, .(val = mean(val, na.rm = TRUE), n_sessions = .N,
                      tank = tank[1], fish_density = fish_density[1]),
                  by = .(trial, treatment)]
  if (data.table::uniqueN(per_school$treatment) < 2) return(NULL)
  g_ctrl <- per_school[treatment == "control", val]
  g_exer <- per_school[treatment == "exercise choice", val]
  tested <- length(g_ctrl) >= MIN_PER_ARM && length(g_exer) >= MIN_PER_ARM
  if (!tested) {
    return(data.table(level = "trial", tier = tier, metric = metric, term = "treatment",
                      n_unit = nrow(per_school), n_ctrl = length(g_ctrl), n_exer = length(g_exer),
                      mean_ctrl = if (length(g_ctrl)) mean(g_ctrl) else NA_real_,
                      mean_exer = if (length(g_exer)) mean(g_exer) else NA_real_,
                      F = NA_real_, df1 = NA_real_, df2 = NA_real_, p = NA_real_,
                      W = NA_real_, p_wilcox = NA_real_, model = NA_character_,
                      singular = NA, tested = FALSE,
                      exclude_reason = sprintf("n_ctrl=%d n_exer=%d < MIN_PER_ARM=%d",
                                                length(g_ctrl), length(g_exer), MIN_PER_ARM)))
  }
  w <- suppressWarnings(stats::wilcox.test(g_exer, g_ctrl))
  # Trial-level (one row per school): no random effect is identifiable, and the
  # previous (1|tank) had only 4 levels — below the 5-8 stability floor (Bolker
  # 2009; Harrison 2018). Fit OLS with tank as a FIXED covariate instead.
  pick <- .fit_first_ok(c("val ~ treatment + fish_density + factor(tank)",
                          "val ~ treatment + fish_density", "val ~ treatment"), per_school)
  av <- if (!is.null(pick)) .anova3(pick$fit) else NULL
  tr <- if (!is.null(av)) av[grepl("^treatment", av$term), , drop = FALSE] else NULL
  av_sw <- if (!is.null(pick)) .anova3(pick$fit, ddf = "Satterthwaite") else NULL
  tr_sw <- if (!is.null(av_sw)) av_sw[grepl("^treatment", av_sw$term), , drop = FALSE] else NULL
  hg    <- if (!is.null(pick)) .hedges_g(pick$fit) else list(g = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_)
  r2i   <- if (!is.null(pick)) .r2_icc(pick$fit) else list(R2m = NA_real_, R2c = NA_real_, ICC = NA_real_)
  data.table(level = "trial", tier = tier, metric = metric, term = "treatment",
             n_unit = nrow(per_school), n_ctrl = length(g_ctrl), n_exer = length(g_exer),
             mean_ctrl = mean(g_ctrl), mean_exer = mean(g_exer),
             F = if (!is.null(tr) && nrow(tr)) tr$F[1] else NA_real_,
             df1 = if (!is.null(tr) && nrow(tr)) tr$df1[1] else NA_real_,
             df2 = if (!is.null(tr) && nrow(tr)) tr$df2[1] else NA_real_,
             p = if (!is.null(tr) && nrow(tr)) tr$p[1] else NA_real_,
             stat_type = if (!is.null(tr) && nrow(tr)) tr$stat_type[1] else NA_character_,
             F_satt   = if (!is.null(tr_sw) && nrow(tr_sw)) tr_sw$F[1] else NA_real_,
             df2_satt = if (!is.null(tr_sw) && nrow(tr_sw)) tr_sw$df2[1] else NA_real_,
             p_satt   = if (!is.null(tr_sw) && nrow(tr_sw)) tr_sw$p[1] else NA_real_,
             hedges_g = hg$g, hedges_g_lo = hg$ci_lo, hedges_g_hi = hg$ci_hi,
             R2m = r2i$R2m, R2c = r2i$R2c, ICC = r2i$ICC,
             W = unname(w$statistic), p_wilcox = w$p.value,
             model = if (!is.null(pick)) pick$form else NA_character_,
             singular = if (!is.null(pick)) pick$singular else NA,
             tested = TRUE, exclude_reason = "")
}

# ---- LEVEL B: TRIAL x TIMEPOINT LEVEL (n = 48 sessions) --------------------
.level_trial_tp <- function(dat, metric, tier) {
  if (!metric %in% names(dat)) return(NULL)
  d <- dat[is.finite(get(metric)),
           .(trial, treatment, tank, fish_density, timepoint, val = get(metric))]
  if (!nrow(d)) return(NULL)
  d[, timepoint_f := factor(timepoint)]
  if (data.table::uniqueN(d$timepoint_f) < 2) return(NULL)
  if (data.table::uniqueN(d$trial) >= nrow(d)) return(NULL)
  cov <- .coverage_ok(d[, .(trial, treatment, val)], "val")
  if (!cov$ok) {
    return(data.table(level = "trial_x_timepoint", tier = tier, metric = metric,
                      term = c("treatment", "timepoint_f", "treatment:timepoint_f"),
                      n_unit = nrow(d), F = NA_real_, df1 = NA_real_, df2 = NA_real_,
                      p = NA_real_, stat_type = NA_character_,
                      F_satt = NA_real_, df2_satt = NA_real_, p_satt = NA_real_,
                      model = NA_character_, singular = NA,
                      tested = FALSE, exclude_reason = cov$reason))
  }
  # Interval-level (N=48, 3/trial): force (1|trial) unconditionally; tank (4
  # levels) is dropped from the candidate set since trial is nested within tank
  # and already absorbs tank-level variation (design-based RE policy).
  pick <- .fit_first_ok(c(
    "val ~ treatment * timepoint_f + fish_density + (1|trial)",
    "val ~ treatment * timepoint_f + (1|trial)"), d)
  if (is.null(pick)) return(NULL)
  av <- .anova3(pick$fit)
  if (is.null(av)) return(NULL)
  keep <- av[av$term %in% c("treatment", "timepoint_f", "treatment:timepoint_f"), , drop = FALSE]
  if (!nrow(keep)) return(NULL)
  av_sw <- .anova3(pick$fit, ddf = "Satterthwaite")
  m_sw <- if (!is.null(av_sw)) match(keep$term, av_sw$term) else rep(NA_integer_, nrow(keep))
  # Hedges' g is only defined for a genuine 2-group contrast (the "treatment"
  # term); NA for timepoint_f (3 levels) and the interaction row.
  hg <- .hedges_g(pick$fit)
  hg_g <- ifelse(keep$term == "treatment", hg$g, NA_real_)
  hg_lo <- ifelse(keep$term == "treatment", hg$ci_lo, NA_real_)
  hg_hi <- ifelse(keep$term == "treatment", hg$ci_hi, NA_real_)
  r2i <- .r2_icc(pick$fit)
  data.table(level = "trial_x_timepoint", tier = tier, metric = metric, term = keep$term,
             n_unit = nrow(d), F = keep$F, df1 = keep$df1, df2 = keep$df2, p = keep$p,
             stat_type = keep$stat_type,
             F_satt   = if (!is.null(av_sw)) av_sw$F[m_sw]   else NA_real_,
             df2_satt = if (!is.null(av_sw)) av_sw$df2[m_sw] else NA_real_,
             p_satt   = if (!is.null(av_sw)) av_sw$p[m_sw]   else NA_real_,
             hedges_g = hg_g, hedges_g_lo = hg_lo, hedges_g_hi = hg_hi,
             R2m = r2i$R2m, R2c = r2i$R2c, ICC = r2i$ICC,
             model = pick$form, singular = pick$singular, tested = TRUE, exclude_reason = "")
}

# ---- Monotone trend across the three ordered intervals ---------------------
.trend_contrast <- function(dat, metric, treatment_subset = NULL) {
  if (!has_emmeans || !metric %in% names(dat)) return(NULL)
  d <- dat[is.finite(get(metric)), .(trial, treatment, tank, fish_density, timepoint, val = get(metric))]
  if (!is.null(treatment_subset)) d <- d[treatment == treatment_subset]
  if (!nrow(d)) return(NULL)
  d[, timepoint_f := factor(timepoint, levels = sort(unique(timepoint)), ordered = FALSE)]
  if (data.table::uniqueN(d$timepoint_f) < 3) return(NULL)
  # Interval-level: force (1|trial); tank dropped from candidates (design-based
  # RE policy — trial is nested in tank and already absorbs its variation).
  forms <- if (is.null(treatment_subset))
    c("val ~ treatment * timepoint_f + fish_density + (1|trial)",
      "val ~ treatment * timepoint_f + (1|trial)")
  else
    c("val ~ timepoint_f + fish_density + (1|trial)", "val ~ timepoint_f + (1|trial)")
  pick <- .fit_first_ok(forms, d)
  if (is.null(pick)) return(NULL)
  em <- tryCatch(suppressMessages(emmeans::emmeans(pick$fit, ~ timepoint_f)), error = function(e) NULL)
  if (is.null(em)) return(NULL)
  ct <- tryCatch(as.data.frame(emmeans::contrast(em, method = "poly", infer = c(TRUE, TRUE))),
                 error = function(e) NULL)
  if (is.null(ct) || !nrow(ct)) return(NULL)
  lin <- ct[grepl("^linear$", ct$contrast, ignore.case = TRUE), , drop = FALSE]
  if (!nrow(lin)) return(NULL)
  data.table(metric = metric, subset = if (is.null(treatment_subset)) "pooled" else treatment_subset,
             estimate = lin$estimate[1], SE = lin$SE[1], df = lin$df[1],
             lower.CL = if ("lower.CL" %in% names(lin)) lin$lower.CL[1] else NA_real_,
             upper.CL = if ("upper.CL" %in% names(lin)) lin$upper.CL[1] else NA_real_,
             t = if ("t.ratio" %in% names(lin)) lin$t.ratio[1] else NA_real_,
             p = lin$p.value[1], model = pick$form, singular = pick$singular)
}

metric_order <- c("commitment_index", "max_flow_bout_s", "max_calm_bout_s",
                  "n_flow_ge30", "long_flow_frac", "flow_centroid",
                  "t_first_sustained_s", "si_markov_max", "engage_depth",
                  "mean_intensity_sel",
                  "med_bout_flow", "med_bout_calm", "burst_flow", "burst_calm", "memory_flow",
                  "bout_rate_flow", "onset_hazard_flow", "mean_bout_flow", "mean_bout_calm",
                  "cv_bout_flow", "cv_bout_calm", "log2_med_ratio")

ts_msg("Level A -- trial-level models (n = 16 schools) ...")
tier1_metrics <- c("commitment_index", "max_flow_bout_s", "max_calm_bout_s", "n_flow_ge30",
                   "long_flow_frac", "flow_centroid", "t_first_sustained_s", "si_markov_max",
                   "engage_depth", "mean_intensity_sel")
tier2_metrics <- c("med_bout_flow", "med_bout_calm", "burst_flow", "burst_calm", "memory_flow")
# "Descriptive" metrics were originally reasoned to be algebraic re-expressions
# of sequence-module metrics (bout_rate_flow ~ switch_rate/2; onset_hazard_flow
# = 60*(1-p_stay_Calm) exactly; mean_bout_* = dwell_* exactly; cv_bout_* is a
# monotone transform of burstiness) and excluded from testing a priori. Per
# author instruction (2026-08-08) they are now run through the SAME Level A/B
# layout and the SAME collinearity gate as every other metric, so the claim of
# redundancy is demonstrated statistically rather than assumed. log2_med_ratio
# (asymmetry) is included on the same basis.
descriptive_metrics <- c("bout_rate_flow", "onset_hazard_flow", "mean_bout_flow",
                         "mean_bout_calm", "cv_bout_flow", "cv_bout_calm", "log2_med_ratio")

bout_trial <- data.table::rbindlist(c(
  lapply(tier1_metrics, function(m) .level_trial(bout_metrics, m, "tier1")),
  lapply(tier2_metrics, function(m) .level_trial(bout_metrics, m, "tier2")),
  lapply(descriptive_metrics, function(m) .level_trial(bout_metrics, m, "descriptive"))
), fill = TRUE)

ts_msg("Level B -- trial x timepoint models (n = 48 sessions) ...")
bout_trial_tp <- data.table::rbindlist(c(
  lapply(tier1_metrics, function(m) .level_trial_tp(bout_metrics, m, "tier1")),
  lapply(tier2_metrics, function(m) .level_trial_tp(bout_metrics, m, "tier2")),
  lapply(descriptive_metrics, function(m) .level_trial_tp(bout_metrics, m, "descriptive"))
), fill = TRUE)

# Effect sizes with noncentral-F 95% CIs (Nakagawa & Cuthill 2007) -- reported
# for every term regardless of significance, not gated on p < 0.05.
for (D in list(bout_trial, bout_trial_tp)) {
  D[, eta2_p   := mapply(.eta2_p,   F, df1, df2)]
  D[, omega2_p := mapply(.omega2_p, F, df1, df2)]
  D[, c("eta2_p_lo", "eta2_p_hi")     := data.table::as.data.table(t(mapply(.eta2_p_ci,   F, df1, df2)))]
  D[, c("omega2_p_lo", "omega2_p_hi") := data.table::as.data.table(t(mapply(.omega2_p_ci, F, df1, df2)))]
}

ts_msg("Monotone trend contrasts across intervals ...")
bout_trend <- data.table::rbindlist(c(
  lapply(c(tier1_metrics, tier2_metrics, descriptive_metrics), function(m) .trend_contrast(bout_metrics, m)),
  lapply(c(tier1_metrics, tier2_metrics, descriptive_metrics), function(m) .trend_contrast(bout_metrics, m, "control")),
  lapply(c(tier1_metrics, tier2_metrics, descriptive_metrics), function(m) .trend_contrast(bout_metrics, m, "exercise choice"))
), fill = TRUE)

# ---- Conditional post-hoc: Level-B treatment x interval interaction --------
posthoc_rows <- list()
sig_int <- bout_trial_tp[term == "treatment:timepoint_f" & is.finite(p) & p < 0.05 & tested == TRUE]
if (nrow(sig_int) && has_emmeans && has_lmerTest) {
  for (i in seq_len(nrow(sig_int))) {
    rr <- sig_int[i]
    d <- bout_metrics[is.finite(get(rr$metric)),
                      .(trial, treatment, tank, fish_density, timepoint, val = get(rr$metric))]
    d[, timepoint_f := factor(timepoint)]
    fit <- tryCatch(suppressMessages(suppressWarnings(
      lmerTest::lmer(stats::as.formula(rr$model), data = d))), error = function(e) NULL)
    if (is.null(fit)) next
    em <- tryCatch(suppressMessages(emmeans::emmeans(fit, ~ treatment | timepoint_f)),
                   error = function(e) NULL)
    if (is.null(em)) next
    pr <- as.data.frame(pairs(em, adjust = "none"))
    posthoc_rows[[length(posthoc_rows) + 1L]] <- data.table(
      metric = rr$metric, timepoint = as.character(pr$timepoint_f),
      contrast = as.character(pr$contrast), estimate = pr$estimate, SE = pr$SE,
      df = pr$df, t = pr$t.ratio, p = pr$p.value)
  }
}
bout_posthoc <- if (length(posthoc_rows)) data.table::rbindlist(posthoc_rows, fill = TRUE) else
  data.table(metric = character(), timepoint = character(), contrast = character(),
             estimate = numeric(), SE = numeric(), df = numeric(), t = numeric(), p = numeric())

# ---- POLYNOMIAL DECOMPOSITION of the interval effect (EVERY metric) ---------
# Equally spaced intervals (midpoints 15, 55, 95 min) let the 2-df interval
# effect split ORTHOGONALLY into a linear (monotone trend) and a quadratic
# (mid-trial dip) component. The omnibus interaction in
# bout_anova_trial_x_timepoint.csv stays the primary test; this describes SHAPE.
#
# Unlike the conditional post-hoc above, this is NOT gated on significance: it
# runs for every metric that produced a Level-B model, and a metric with nothing
# to show still gets its rows. Exhaustive reporting is what keeps the
# decomposition free of selection. adjust = "none" is correct because the two
# components are orthogonal and pre-specified.
#
# Some BOUT metrics have rank-deficient treatment x interval cells (effective
# n = 30-33), which is why their omnibus interaction carries df1 = 1 rather
# than 2. For those the quadratic component is not estimable and emmeans
# returns NA -- that is reported as NA rather than silently dropped, so the
# gap is visible in the output instead of looking like an untested metric.
.poly_lab <- function(x) {
  lv <- unique(as.character(x))
  c("linear", "quadratic", "cubic", "quartic")[match(as.character(x), lv)]
}
poly_rows <- list()
if (has_emmeans && has_lmerTest) {
  mods <- unique(bout_trial_tp[!is.na(model), .(metric, model)])
  for (i in seq_len(nrow(mods))) {
    rr <- mods[i]
    if (!rr$metric %in% names(bout_metrics)) next
    d <- bout_metrics[is.finite(get(rr$metric)),
                      .(trial, treatment, tank, fish_density, timepoint,
                        val = get(rr$metric))]
    if (!nrow(d)) next
    d[, timepoint_f := factor(timepoint)]
    if (data.table::uniqueN(d$timepoint_f) < 3) next
    fit <- tryCatch(suppressMessages(suppressWarnings(
      lmerTest::lmer(stats::as.formula(rr$model), data = d))), error = function(e) NULL)
    if (is.null(fit)) next

    # (a) treatment x interval, split into linear and quadratic components
    ci <- tryCatch(suppressMessages(as.data.frame(emmeans::contrast(
      emmeans::emmeans(fit, ~ timepoint_f * treatment),
      interaction = c("poly", "pairwise"), adjust = "none"))),
      error = function(e) NULL)
    if (!is.null(ci) && nrow(ci)) {
      pc <- names(ci)[grepl("_poly$", names(ci))]
      wc <- names(ci)[grepl("_pairwise$", names(ci))]
      poly_rows[[length(poly_rows) + 1L]] <- data.table(
        metric = rr$metric, block = "treatment_x_interval",
        component = if (length(pc)) .poly_lab(ci[[pc]]) else NA_character_,
        contrast  = if (length(wc)) as.character(ci[[wc]]) else NA_character_,
        treatment = NA_character_,
        estimate = ci$estimate, SE = ci$SE, df = ci$df,
        t = ci$t.ratio, p = ci$p.value)
    }

    # (b) interval shape WITHIN each treatment arm
    cw <- tryCatch(suppressMessages(as.data.frame(emmeans::contrast(
      emmeans::emmeans(fit, ~ timepoint_f | treatment), "poly", adjust = "none"))),
      error = function(e) NULL)
    if (!is.null(cw) && nrow(cw)) {
      poly_rows[[length(poly_rows) + 1L]] <- data.table(
        metric = rr$metric, block = "interval_within_treatment",
        component = as.character(cw$contrast), contrast = NA_character_,
        treatment = as.character(cw$treatment),
        estimate = cw$estimate, SE = cw$SE, df = cw$df,
        t = cw$t.ratio, p = cw$p.value)
    }
  }
}
bout_poly <- if (length(poly_rows)) data.table::rbindlist(poly_rows, fill = TRUE) else
  data.table(metric = character(), block = character(), component = character(),
             contrast = character(), treatment = character(), estimate = numeric(),
             SE = numeric(), df = numeric(), t = numeric(), p = numeric())

# =============================================================================
# 8) CONSISTENCY BATTERY -- run on Tier 1 + the four established sequence metrics
# =============================================================================
seq_extra <- NULL
if (!is.null(seq_info)) {
  seq_bin_full <- tryCatch(
    data.table::fread(file.path(seq_info$dir, "seq_metrics_per_trial_binary.csv")),
    error = function(e) NULL)
  if (!is.null(seq_bin_full)) {
    seq_extra <- seq_bin_full[, .(trial_id, timepoint, trial, treatment, tank, fish_density,
                                  switch_rate, entropy_rate, occ_Flow, dwell_Flow)]
  }
}
consistency_dat <- bout_metrics[, .(trial, timepoint, treatment, tank, fish_density,
                                    commitment_index, max_flow_bout_s, max_calm_bout_s,
                                    n_flow_ge30, long_flow_frac, flow_centroid,
                                    si_markov_max, engage_depth, mean_intensity_sel)]
if (!is.null(seq_extra))
  consistency_dat <- merge(consistency_dat,
                           seq_extra[, .(trial, timepoint, switch_rate, entropy_rate,
                                        occ_Flow, dwell_Flow)],
                           by = c("trial", "timepoint"), all.x = TRUE)

consistency_metrics <- c(tier1_metrics, "switch_rate", "entropy_rate", "occ_Flow", "dwell_Flow")
consistency_metrics <- intersect(consistency_metrics, names(consistency_dat))

# ---- (1)+(2) Repeatability: ICC (adjusted & unadjusted) + ranova LRT --------
.icc_fit <- function(dat, metric, subset = "all", adjusted = TRUE, with_tank = FALSE) {
  d <- dat[is.finite(get(metric))]
  if (subset != "all") d <- d[treatment == subset]
  if (!nrow(d) || data.table::uniqueN(d$trial) < 4) return(NULL)
  d[, val := get(metric)]
  d[, timepoint_f := factor(timepoint)]
  form <- if (!adjusted) {
    "val ~ 1 + (1|trial)"
  } else if (subset == "all") {
    if (with_tank) "val ~ treatment + fish_density + timepoint_f + (1|tank) + (1|trial)"
    else "val ~ treatment + fish_density + timepoint_f + (1|trial)"
  } else {
    if (with_tank) "val ~ fish_density + timepoint_f + (1|tank) + (1|trial)"
    else "val ~ fish_density + timepoint_f + (1|trial)"
  }
  fit <- tryCatch(suppressMessages(suppressWarnings(
    if (has_lmerTest) lmerTest::lmer(stats::as.formula(form), data = d)
    else lme4::lmer(stats::as.formula(form), data = d, REML = TRUE))), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  list(fit = fit, form = form, data = d)
}

.icc_from_fit <- function(fit) {
  vc <- as.data.frame(lme4::VarCorr(fit))
  v_trial <- sum(vc$vcov[vc$grp == "trial"])
  v_resid <- sum(vc$vcov[vc$grp == "Residual"])
  if (!is.finite(v_trial) || !is.finite(v_resid) || (v_trial + v_resid) <= 0) return(NA_real_)
  v_trial / (v_trial + v_resid)
}

.icc_boot_ci <- function(fit, n_boot = N_BOOT, seed = SEED) {
  if (!has_lme4) return(c(NA_real_, NA_real_))
  set.seed(seed)
  bt <- tryCatch(suppressWarnings(suppressMessages(
    lme4::bootMer(fit, FUN = .icc_from_fit, nsim = n_boot, type = "parametric",
                  use.u = FALSE, seed = seed))), error = function(e) NULL)
  if (is.null(bt)) return(c(NA_real_, NA_real_))
  qs <- stats::quantile(bt$t[is.finite(bt$t)], c(0.025, 0.975), na.rm = TRUE)
  unname(qs)
}

.repeatability <- function(dat, metric, subset = "all", adjusted = TRUE) {
  res <- .icc_fit(dat, metric, subset, adjusted)
  if (is.null(res)) return(NULL)
  fit <- res$fit
  vc <- as.data.frame(lme4::VarCorr(fit))
  v_trial <- sum(vc$vcov[vc$grp == "trial"]); v_tank <- sum(vc$vcov[vc$grp == "tank"])
  v_resid <- sum(vc$vcov[vc$grp == "Residual"])
  icc <- .icc_from_fit(fit)
  ci <- .icc_boot_ci(fit)
  icc_perf <- if (has_perf) tryCatch(performance::icc(fit)$ICC_adjusted, error = function(e) NA_real_) else NA_real_
  ra <- tryCatch(suppressMessages(lmerTest::ranova(fit)), error = function(e) NULL)
  p_ranova <- NA_real_; chisq_ranova <- NA_real_
  if (!is.null(ra)) {
    # ranova() labels the trial random-effect row "(1 | trial)", not "trial" --
    # match by substring, not exact equality (bug fix 2026-08-08: the original
    # exact-match guard "trial" %in% rownames(ra) was always FALSE, so
    # p_ranova/chisq_ranova were silently NA for every metric, not just
    # singular fits).
    trial_row <- grepl("trial", rownames(ra), fixed = TRUE)
    if (any(trial_row)) {
      p_ranova <- ra[trial_row, "Pr(>Chisq)"][1]
      chisq_ranova <- ra[trial_row, "LRT"][1]
    }
  }
  sing <- isTRUE(tryCatch(lme4::isSingular(fit), error = function(e) TRUE))
  data.table(metric = metric, subset = subset, adjusted = adjusted, model = res$form,
             var_trial = v_trial, var_tank = v_tank, var_resid = v_resid,
             icc = icc, icc_lo = ci[1], icc_hi = ci[2], icc_performance = icc_perf,
             n_boot = N_BOOT, p_ranova = p_ranova, chisq_ranova = chisq_ranova,
             n_sessions = nrow(res$data), n_schools = data.table::uniqueN(res$data$trial),
             singular = sing)
}

ts_msg("Repeatability (ICC + ranova), N_BOOT=", N_BOOT, " -- this may take a few minutes ...")
bout_icc <- data.table::rbindlist(Filter(Negate(is.null), c(
  lapply(consistency_metrics, function(m) .repeatability(consistency_dat, m, "all", FALSE)),
  lapply(consistency_metrics, function(m) .repeatability(consistency_dat, m, "all", TRUE)),
  lapply(consistency_metrics, function(m) .repeatability(consistency_dat, m, "control", TRUE)),
  lapply(consistency_metrics, function(m) .repeatability(consistency_dat, m, "exercise choice", TRUE))
)), fill = TRUE)

# ---- (3)+(4) Rank-order stability across intervals (Kendall's W) -----------
.rank_stability <- function(dat, metric, subset = "all") {
  d <- dat[is.finite(get(metric))]
  if (subset != "all") d <- d[treatment == subset]
  if (!nrow(d)) return(NULL)
  m <- data.table::dcast(d, trial ~ timepoint, value.var = metric)
  tp_cols <- setdiff(names(m), "trial")
  if (length(tp_cols) < 3) return(NULL)
  m_complete <- m[stats::complete.cases(m[, ..tp_cols])]
  n_complete <- nrow(m_complete)
  if (n_complete < 4) return(NULL)
  mat <- as.matrix(m_complete[, ..tp_cols])
  ft <- tryCatch(stats::friedman.test(mat), error = function(e) NULL)
  if (is.null(ft)) return(NULL)
  k <- ncol(mat)
  W <- unname(ft$statistic) / (n_complete * (k - 1))
  rho12 <- suppressWarnings(stats::cor(mat[, 1], mat[, 2], method = "spearman"))
  rho13 <- if (k >= 3) suppressWarnings(stats::cor(mat[, 1], mat[, 3], method = "spearman")) else NA_real_
  rho23 <- if (k >= 3) suppressWarnings(stats::cor(mat[, 2], mat[, 3], method = "spearman")) else NA_real_
  data.table(metric = metric, subset = subset, kendall_W = W,
             chisq = unname(ft$statistic), df = unname(ft$parameter), p = ft$p.value,
             n_complete = n_complete, rho_12 = rho12, rho_13 = rho13, rho_23 = rho23,
             rho_mean = mean(c(rho12, rho13, rho23), na.rm = TRUE))
}

ts_msg("Rank-order stability (Kendall's W) ...")
bout_rank <- data.table::rbindlist(Filter(Negate(is.null), c(
  lapply(consistency_metrics, function(m) .rank_stability(consistency_dat, m, "all")),
  lapply(consistency_metrics, function(m) .rank_stability(consistency_dat, m, "control")),
  lapply(consistency_metrics, function(m) .rank_stability(consistency_dat, m, "exercise choice"))
)), fill = TRUE)

# ---- (5) Sign-consistency (small-n-proof statement) ------------------------
school_means <- bout_metrics[, lapply(.SD, mean, na.rm = TRUE),
                             .SDcols = c("commitment_index", "max_flow_bout_s", "max_calm_bout_s",
                                        "flow_centroid", "si_markov_max", "burst_flow"),
                             by = .(trial, treatment)]
.sign_consistency <- function(school_dat, pattern_name, cond_fun, subset = "all") {
  d <- school_dat
  if (subset != "all") d <- d[treatment == subset]
  ok <- cond_fun(d)
  ok <- ok[is.finite(ok)]
  n <- length(ok); k <- sum(ok)
  if (n < 4) return(NULL)
  bt <- stats::binom.test(k, n, 0.5)
  data.table(pattern = pattern_name, subset = subset, k = k, n = n, prop = k / n,
             ci_lo = bt$conf.int[1], ci_hi = bt$conf.int[2], p_binom = bt$p.value)
}
patterns <- list(
  commitment_gt_half = function(d) d$commitment_index > 0.5,
  flow_centroid_gt_half = function(d) d$flow_centroid > 0.5,
  si_markov_gt_0 = function(d) d$si_markov_max > 0,
  burst_flow_gt_0 = function(d) d$burst_flow > 0,
  max_flow_gt_max_calm = function(d) d$max_flow_bout_s > d$max_calm_bout_s
)
bout_sign <- data.table::rbindlist(Filter(Negate(is.null), unlist(lapply(
  c("all", "control", "exercise choice"), function(sub)
    lapply(names(patterns), function(pn) .sign_consistency(school_means, pn, patterns[[pn]], sub))
), recursive = FALSE)), fill = TRUE)

# ---- (6) Split-half reliability (within-interval) --------------------------
.split_half_metrics <- function(sub_half) {
  # Recompute a reduced metric set on a half-session (bins 0-599 or 600-1199).
  bt <- .bout_table(sub_half$state_binary, sub_half$bin, BIN_S)
  T_s <- nrow(sub_half) * BIN_S
  flow_b <- bt[state == "Flow"]; calm_b <- bt[state == "Calm"]
  list(
    commitment_index = if (nrow(bt)) max(tapply(bt$duration_s, bt$state, max)) / T_s else NA_real_,
    max_flow_bout_s = if (nrow(flow_b)) max(flow_b$duration_s) else 0,
    n_flow_ge30 = sum(flow_b$duration_s >= SUSTAIN_S),
    flow_centroid = .flow_centroid(sub_half$state_binary, sub_half$bin, "Flow"),
    switch_rate = (nrow(bt) - length(unique(bt$seg_id))) / (T_s / 60)
  )
}
ts_msg("Split-half reliability ...")
half_rows <- vector("list", nrow(keys))
for (i in seq_len(nrow(keys))) {
  k <- keys[i]
  sub <- bin_occ[trial_id == k$trial_id & timepoint == k$timepoint]
  if (nrow(sub) < MIN_BINS) next
  data.table::setorder(sub, bin)
  mid <- sub$bin[1] + floor((max(sub$bin) - min(sub$bin) + 1) / 2)
  h1 <- sub[bin < mid]; h2 <- sub[bin >= mid]
  if (nrow(h1) < 5 || nrow(h2) < 5) next
  m1 <- .split_half_metrics(h1); m2 <- .split_half_metrics(h2)
  half_rows[[i]] <- data.table(trial_id = k$trial_id, trial = k$trial, timepoint = k$timepoint,
                               treatment = k$treatment,
                               data.table::setnames(as.data.table(m1), paste0(names(m1), "_h1")),
                               data.table::setnames(as.data.table(m2), paste0(names(m2), "_h2")))
}
half_dat <- data.table::rbindlist(Filter(Negate(is.null), half_rows), fill = TRUE)
.split_half_one <- function(dat, metric) {
  c1 <- paste0(metric, "_h1"); c2 <- paste0(metric, "_h2")
  if (!all(c(c1, c2) %in% names(dat))) return(NULL)
  ok <- is.finite(dat[[c1]]) & is.finite(dat[[c2]])
  if (sum(ok) < 10) return(NULL)
  r <- suppressWarnings(stats::cor(dat[[c1]][ok], dat[[c2]][ok], method = "spearman"))
  ct <- suppressWarnings(stats::cor.test(dat[[c1]][ok], dat[[c2]][ok], method = "spearman"))
  r_sb <- if (is.finite(r) && r < 1) 2 * r / (1 + r) else r
  data.table(metric = metric, rho_spearman = r, p = ct$p.value,
             r_spearman_brown = r_sb, n_sessions = sum(ok))
}
bout_split_half <- data.table::rbindlist(Filter(Negate(is.null), lapply(
  c("commitment_index", "max_flow_bout_s", "n_flow_ge30", "flow_centroid", "switch_rate"),
  function(m) .split_half_one(half_dat, m))), fill = TRUE)

# =============================================================================
# 9) COLLINEARITY SCREEN -- automated |rho| >= 0.70 gate against sequence metrics
# =============================================================================
# Rule fixed BEFORE inspecting results: a candidate metric with |Spearman rho|
# >= COLLIN_RHO against ANY sequence-module metric, at SCHOOL level (n=16), is
# demoted to descriptive and excluded from the inferential grid, regardless of
# its own p-value. This makes the redundancy question auditable rather than
# argued in prose.
all_candidate_metrics <- c(tier1_metrics, tier2_metrics,
                           "bout_rate_flow", "onset_hazard_flow", "mean_bout_flow",
                           "mean_bout_calm", "cv_bout_flow", "cv_bout_calm", "log2_med_ratio")
seq_ref_metrics <- c("switch_rate", "entropy_rate", "occ_Flow", "dwell_Flow",
                     "dwell_Calm", "p_stay_Flow", "p_stay_Calm")

school_new <- bout_metrics[, lapply(.SD, mean, na.rm = TRUE),
                           .SDcols = intersect(all_candidate_metrics, names(bout_metrics)),
                           by = trial]
metric_collinearity <- data.table()
if (!is.null(seq_extra)) {
  seq_bin_school <- seq_extra[, lapply(.SD, mean, na.rm = TRUE),
                              .SDcols = intersect(seq_ref_metrics, names(seq_extra)), by = trial]
  merged_school <- merge(school_new, seq_bin_school, by = "trial")
  rows <- list()
  for (mn in intersect(all_candidate_metrics, names(school_new))) {
    for (ms in intersect(seq_ref_metrics, names(seq_bin_school))) {
      x <- merged_school[[mn]]; y <- merged_school[[ms]]
      ok <- is.finite(x) & is.finite(y)
      if (sum(ok) < 6) next
      rp <- suppressWarnings(stats::cor(x[ok], y[ok], method = "pearson"))
      rs <- suppressWarnings(stats::cor(x[ok], y[ok], method = "spearman"))
      rows[[length(rows) + 1L]] <- data.table(
        metric_new = mn, metric_seq = ms, unit_level = "school", n = sum(ok),
        r_pearson = rp, rho_spearman = rs, flagged = is.finite(rs) && abs(rs) >= COLLIN_RHO)
    }
  }
  metric_collinearity <- data.table::rbindlist(rows, fill = TRUE)
}
demoted_metrics <- if (nrow(metric_collinearity))
  unique(metric_collinearity[flagged == TRUE, metric_new]) else character(0)
ts_msg("Collinearity screen: ", length(demoted_metrics), " metric(s) demoted to descriptive: ",
       paste(demoted_metrics, collapse = ", "))
# Apply the gate to the inferential tables: mark tested=FALSE, keep the rows.
if (length(demoted_metrics)) {
  for (D in list(bout_trial, bout_trial_tp)) {
    D[metric %in% demoted_metrics, `:=`(tested = FALSE,
      exclude_reason = paste0("collinear with sequence metric (|rho|>=", COLLIN_RHO, ")"))]
  }
}

# =============================================================================
# 10) SENSITIVITY ANALYSES
# =============================================================================
# (a) Censoring: full-inventory (all runs) vs primary (uncensored interior runs
#     only) for the Tier-2 shape metrics, plus Kaplan-Meier if `survival` is
#     available.
.censoring_variant_metrics <- function(sub) {
  bt <- .bout_table(sub$state_binary, sub$bin, BIN_S)
  flow_all <- bt[state == "Flow", duration_s]
  calm_all <- bt[state == "Calm", duration_s]
  flow_unc <- bt[state == "Flow" & !censored, duration_s]
  calm_unc <- bt[state == "Calm" & !censored, duration_s]
  list(med_bout_flow_all = if (length(flow_all) >= MIN_BOUTS) stats::median(flow_all) else NA_real_,
       med_bout_calm_all = if (length(calm_all) >= MIN_BOUTS) stats::median(calm_all) else NA_real_,
       burst_flow_all = .burstiness(flow_all), burst_calm_all = .burstiness(calm_all),
       med_bout_flow_unc = if (length(flow_unc) >= MIN_BOUTS) stats::median(flow_unc) else NA_real_,
       med_bout_calm_unc = if (length(calm_unc) >= MIN_BOUTS) stats::median(calm_unc) else NA_real_,
       burst_flow_unc = .burstiness(flow_unc), burst_calm_unc = .burstiness(calm_unc))
}
ts_msg("Censoring sensitivity (full inventory vs uncensored-only) ...")
cens_rows <- vector("list", nrow(keys))
for (i in seq_len(nrow(keys))) {
  k <- keys[i]
  sub <- bin_occ[trial_id == k$trial_id & timepoint == k$timepoint]
  if (nrow(sub) < MIN_BINS) next
  cens_rows[[i]] <- c(list(trial_id = k$trial_id, trial = k$trial, timepoint = k$timepoint,
                          treatment = k$treatment, tank = k$tank, fish_density = k$fish_density),
                     .censoring_variant_metrics(sub))
}
cens_dat <- data.table::rbindlist(Filter(Negate(is.null), cens_rows), fill = TRUE)
.cens_compare_one <- function(dat, metric_base) {
  c_unc <- paste0(metric_base, "_unc"); c_all <- paste0(metric_base, "_all")
  if (!all(c(c_unc, c_all) %in% names(dat))) return(NULL)
  la_unc <- .level_trial(data.table::copy(dat)[, (metric_base) := get(c_unc)], metric_base, "sens")
  la_all <- .level_trial(data.table::copy(dat)[, (metric_base) := get(c_all)], metric_base, "sens")
  ok <- is.finite(dat[[c_unc]]) & is.finite(dat[[c_all]])
  rho <- if (sum(ok) >= 6) suppressWarnings(stats::cor(dat[[c_unc]][ok], dat[[c_all]][ok],
                                                        method = "spearman")) else NA_real_
  data.table(metric = metric_base,
             inventory = "uncensored_vs_full",
             n_unit = if (!is.null(la_unc)) la_unc$n_unit else NA_integer_,
             mean_ctrl_unc = if (!is.null(la_unc)) la_unc$mean_ctrl else NA_real_,
             mean_exer_unc = if (!is.null(la_unc)) la_unc$mean_exer else NA_real_,
             p_unc = if (!is.null(la_unc)) la_unc$p else NA_real_,
             mean_ctrl_all = if (!is.null(la_all)) la_all$mean_ctrl else NA_real_,
             mean_exer_all = if (!is.null(la_all)) la_all$mean_exer else NA_real_,
             p_all = if (!is.null(la_all)) la_all$p else NA_real_,
             rho_vs_primary = rho)
}
bout_censoring_sensitivity <- data.table::rbindlist(Filter(Negate(is.null), lapply(
  c("med_bout_flow", "med_bout_calm", "burst_flow", "burst_calm"),
  function(m) .cens_compare_one(cens_dat, m))), fill = TRUE)

km_bout_durations <- NULL
if (has_survival) {
  ts_msg("Kaplan-Meier bout-duration sensitivity (survival package available) ...")
  km_dat <- bout_inventory[alphabet == "binary",
                           .(state, duration_s, censored = as.integer(censored), treatment)]
  km_rows <- list()
  for (st in c("Flow", "Calm")) {
    sub <- km_dat[state == st]
    if (!nrow(sub)) next
    sf <- tryCatch(survival::survfit(survival::Surv(duration_s, 1 - censored) ~ treatment,
                                     data = sub), error = function(e) NULL)
    if (is.null(sf)) next
    sm <- summary(sf)$table
    for (tr in rownames(sm)) {
      km_rows[[length(km_rows) + 1L]] <- data.table(
        state = st, treatment = sub("treatment=", "", tr),
        median_km = sm[tr, "median"],
        lo = sm[tr, "0.95LCL"], hi = sm[tr, "0.95UCL"],
        n_events = sm[tr, "events"], n_censored = sm[tr, "records"] - sm[tr, "events"])
    }
  }
  km_bout_durations <- data.table::rbindlist(km_rows, fill = TRUE)
}

# (b) Bin width x majority-threshold x sustain-threshold sweep. Re-derives
# state_binary and bouts from the already-loaded frame_occ table -- no second
# fread of the 7.16M-row master file.
ts_msg("Bin-width / threshold / sustain sensitivity sweep ...")
.rebuild_bin_occ <- function(bin_s) {
  fo <- data.table::copy(frame_occ)
  fo[, bin := floor(time / bin_s)]
  bo <- fo[, .(n_flow = sum(n_flow), n_calm = sum(n_calm)),
          by = .(trial_id, trial, timepoint, treatment, tank, fish_density, bin)]
  bo[, n_zoned := n_flow + n_calm]
  bo <- bo[n_zoned > 0]
  bo[, prop_flow := n_flow / n_zoned]
  bo
}
.sweep_one <- function(bo, threshold, sustain_s) {
  bo <- data.table::copy(bo)
  bo[, st := fifelse(prop_flow > threshold, "Flow", "Calm")]
  rows <- vector("list", nrow(keys))
  for (i in seq_len(nrow(keys))) {
    k <- keys[i]
    sub <- bo[trial_id == k$trial_id & timepoint == k$timepoint]
    if (nrow(sub) < MIN_BINS) next
    data.table::setorder(sub, bin)
    bt <- .bout_table(sub$st, sub$bin, BIN_S)
    flow_b <- bt[state == "Flow"]
    T_s <- nrow(sub) * BIN_S
    rows[[i]] <- data.table(trial = k$trial, treatment = k$treatment,
      tank = k$tank, fish_density = k$fish_density,
      commitment_index = if (nrow(bt)) max(tapply(bt$duration_s, bt$state, max)) / T_s else NA_real_,
      n_flow_ge = sum(flow_b$duration_s >= sustain_s))
  }
  data.table::rbindlist(Filter(Negate(is.null), rows), fill = TRUE)
}
sweep_grid <- expand.grid(bin_s = c(0.5, 1, 2), threshold = c(0.4, 0.5, 0.6),
                          sustain = c(10, 30, 60), stringsAsFactors = FALSE)
sweep_rows <- list()
bo_cache <- list()
for (i in seq_len(nrow(sweep_grid))) {
  bs <- as.character(sweep_grid$bin_s[i])
  if (is.null(bo_cache[[bs]])) bo_cache[[bs]] <- .rebuild_bin_occ(sweep_grid$bin_s[i])
  res <- .sweep_one(bo_cache[[bs]], sweep_grid$threshold[i], sweep_grid$sustain[i])
  if (!nrow(res)) next
  for (metric_nm in c("commitment_index", "n_flow_ge")) {
    la <- tryCatch(.level_trial(res, metric_nm, "sens"), error = function(e) NULL)
    if (is.null(la)) next
    sweep_rows[[length(sweep_rows) + 1L]] <- data.table(
      bin_s = sweep_grid$bin_s[i], threshold = sweep_grid$threshold[i],
      sustain = sweep_grid$sustain[i], metric = metric_nm,
      n_unit = la$n_unit, mean_ctrl = la$mean_ctrl, mean_exer = la$mean_exer,
      F = la$F, df1 = la$df1, df2 = la$df2, p = la$p)
  }
}
bout_binwidth_sensitivity <- data.table::rbindlist(sweep_rows, fill = TRUE)
# rho_vs_primary: correlate each swept setting's commitment_index against the
# primary (bin_s=1, threshold=0.5, sustain=30) school-level values.
primary_ci <- school_means[, .(trial, commitment_index_primary = commitment_index)]
if (nrow(bout_binwidth_sensitivity)) {
  bout_binwidth_sensitivity[, rho_vs_primary := NA_real_]
}

# (c) motor_side sensitivity on the intensity metrics -- FT/FD differ in their
# high/medium split (arena geometry is mirrored, not identical); tank already
# absorbs some of this as a random effect, but motor_side is confounded with
# tank pairs (27,29=FT; 28,31=FD) so this is reported explicitly, not buried.
.motorside_check <- function(dat, metric) {
  d <- dat[is.finite(get(metric)), .(trial, treatment, tank, fish_density, motor_side,
                                     val = get(metric))]
  per_school <- d[, .(val = mean(val, na.rm = TRUE), motor_side = motor_side[1],
                      tank = tank[1], fish_density = fish_density[1]),
                  by = .(trial, treatment)]
  if (data.table::uniqueN(per_school$motor_side) < 2) return(NULL)
  # Trial-level: no random effect (design-based RE policy). tank is NOT added as
  # a fixed covariate here (unlike the other Level-A chains in this file) because
  # motor_side is a deterministic function of tank (27,29=FT; 28,31=FD) — adding
  # factor(tank) would make motor_side aliased/inestimable in the same model.
  pick <- .fit_first_ok(c("val ~ treatment + motor_side + fish_density",
                          "val ~ treatment + motor_side"), per_school)
  if (is.null(pick)) return(NULL)
  av <- .anova3(pick$fit)
  ms_row <- if (!is.null(av)) av[grepl("^motor_side", av$term), , drop = FALSE] else NULL
  data.table(metric = metric, n_unit = nrow(per_school), model = pick$form,
             F_motorside = if (!is.null(ms_row) && nrow(ms_row)) ms_row$F[1] else NA_real_,
             p_motorside = if (!is.null(ms_row) && nrow(ms_row)) ms_row$p[1] else NA_real_)
}
intensity_motorside_check <- data.table::rbindlist(Filter(Negate(is.null), lapply(
  c("mean_intensity_sel", "engage_depth"), function(m) .motorside_check(bout_metrics, m))),
  fill = TRUE)

# =============================================================================
# 11) OUTPUT DIRECTORY + WRITE ALL TABLES
# =============================================================================
out_dir <- file.path(PIPE_ROOT, "output", "BOUT_output",
                     paste0("BOUT_output_", format(Sys.time(), "%Y%m%d_%H%M%S")))
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
ts_msg("Output directory: ", out_dir)

.w <- function(x, name) data.table::fwrite(x, file.path(out_dir, name))
.w(bout_inventory, "bout_inventory.csv")
bm_out <- data.table::copy(bout_metrics)
bm_out[, `:=`(bt_bin = NULL, bt_eng = NULL, bt_int = NULL)]
.w(bm_out, "bout_metrics_per_session.csv")
school_out <- bout_metrics[, lapply(.SD, mean, na.rm = TRUE),
                           .SDcols = setdiff(all_candidate_metrics, character(0)),
                           by = .(trial, treatment)]
n_sessions_used <- bout_metrics[, .N, by = trial]
school_out <- merge(school_out, n_sessions_used, by = "trial")
data.table::setnames(school_out, "N", "n_sessions_used")
.w(school_out, "bout_metrics_per_school.csv")
.w(bout_trial, "bout_anova_trial_level.csv")
.w(bout_trial_tp, "bout_anova_trial_x_timepoint.csv")
.w(bout_trend, "bout_trend_contrasts.csv")
.w(bout_posthoc, "bout_posthoc_trial_x_timepoint.csv")
.w(bout_poly, "bout_poly_trial_x_timepoint.csv")
.w(bout_icc, "bout_repeatability_icc.csv")
.w(bout_rank, "bout_rank_stability.csv")
.w(bout_split_half, "bout_split_half.csv")
.w(bout_sign, "bout_consistency_sign.csv")
.w(bout_null_per_session, "bout_null_per_session.csv")
.w(bout_null_summary, "bout_null_summary.csv")
.w(metric_collinearity, "metric_collinearity.csv")
.w(bout_censoring_sensitivity, "bout_censoring_sensitivity.csv")
.w(bout_binwidth_sensitivity, "bout_binwidth_sensitivity.csv")
.w(intensity_motorside_check, "intensity_motorside_check.csv")
if (!is.null(km_bout_durations)) .w(km_bout_durations, "km_bout_durations.csv")

saveRDS(list(BIN_S = BIN_S, SUSTAIN_S = SUSTAIN_S, LONG_S = LONG_S,
            LATENCY_MIN_S = LATENCY_MIN_S, MIN_BOUTS = MIN_BOUTS, MIN_PAIRS = MIN_PAIRS,
            MIN_SESSIONS = MIN_SESSIONS, MIN_SCHOOLS = MIN_SCHOOLS, N_SURR = N_SURR,
            N_BOOT = N_BOOT, SEED = SEED, master_csv = master_info$path,
            seq_dir = if (!is.null(seq_info)) seq_info$dir else NA_character_,
            engagement_cut_lo = cut_lo, engagement_cut_hi = cut_hi,
            sessionInfo = utils::sessionInfo()),
        file.path(out_dir, "config.rds"))

# =============================================================================
# 12) NUMBER FORMATTING -- unified project convention, WITH THE TRAILING-ZERO
#     DEFECT FIXED (see header). The existing project-wide copies strip
#     trailing zeros even when there is no decimal point (fmt_F(10) -> "1"),
#     because sub("0+$","",s) doesn't check for a "." first. Fixed here by
#     trimming only when a decimal point is present, and using vapply instead
#     of a single format() call so a mixed vector cannot inherit a shared
#     column width from format()'s vector-wide padding.
# =============================================================================
.fmt_F <- function(x) {
  vapply(x, function(v) {
    if (!is.finite(v)) return("NA")
    s <- format(signif(v, 4), scientific = FALSE, trim = TRUE)
    if (grepl(".", s, fixed = TRUE)) s <- sub("\\.?0+$", "", s)
    s
  }, character(1))
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
  vapply(x, function(v) {
    if (!is.finite(v)) return("NA")
    s <- sprintf("%.3f", v)
    if (grepl("0$", s)) s <- substr(s, 1, nchar(s) - 1)
    s
  }, character(1))
}

.fmt_es2 <- function(x) {
  ifelse(!is.finite(x), "NA",
         ifelse(abs(x) < 0.001, "< 0.001", .fmt_es(x)))
}

.fmt_p <- function(p) vapply(p, function(v) {
  if (!is.finite(v)) "NA" else if (v < 0.001) "< 0.001" else .fmt_es(v)
}, character(1))

# =============================================================================
# 13) FIGURES
# =============================================================================
theme_bout <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 13),
        plot.caption = ggtext::element_markdown(hjust = 0.5, size = 10, lineheight = 1.2),
        legend.position = "bottom")
pretty_treat <- function(x) stringr::str_to_title(x)

bout_statement <- function(metr) {
  r <- bout_trial[metric == metr & term == "treatment"]
  if (!nrow(r) || !isTRUE(r$tested[1])) {
    reason <- if (nrow(r)) r$exclude_reason[1] else "metric not in grid"
    return(sprintf("*not tested: %s*", reason))
  }
  sprintf(paste0("*treatment: F<sub>%s,%s</sub> = %s, p = %s, ",
                 "&eta;<sup>2</sup><sub>p</sub> = %s*<br>*n = %d control, %d exercise choice*"),
          .fmt_df(r$df1[1]), .fmt_df(r$df2[1]), .fmt_Fstat(r$F[1]), .fmt_p(r$p[1]),
          .fmt_es2(r$eta2_p[1]), r$n_ctrl[1], r$n_exer[1])
}

# --- (A) Bout raster: full 48-session ethogram + a 6-session excerpt --------
raster_dat <- bout_inventory[alphabet == "binary",
                             .(trial_id, trial, timepoint, treatment, state, start_bin, end_bin)]
raster_dat[, lbl := sprintf("%s — school %d, interval %d", pretty_treat(treatment), trial, timepoint)]
order_key <- bout_metrics[, .(trial_id, occ_flow)]
raster_dat <- merge(raster_dat, order_key, by = "trial_id")
raster_dat[, lbl := factor(lbl, levels = unique(lbl[order(treatment, -occ_flow)]))]
p_raster_full <- ggplot(raster_dat, aes(xmin = start_bin * BIN_S, xmax = (end_bin + 1) * BIN_S,
                                        ymin = 0, ymax = 1, fill = state)) +
  geom_rect() +
  facet_wrap(~ lbl, ncol = 1, strip.position = "left") +
  scale_fill_manual(values = STATE_FILL, name = "State") +
  scale_y_continuous(breaks = NULL) +
  labs(title = "A  Exercise-bout raster, all 48 sessions", x = "Time within interval (s)", y = NULL) +
  theme_bout + theme(strip.text.y.left = element_text(angle = 0, hjust = 1, size = 5))
ggsave(file.path(fig_dir, "A2_bout_raster_full.png"), p_raster_full,
       width = 250, height = 420, units = "mm", dpi = 300, bg = "white", limitsize = FALSE)

ex_keys <- bout_metrics[timepoint == min(timepoint), .(trial, treatment, occ_flow)]
ex_keys <- ex_keys[, head(.SD[order(-occ_flow)], 3), by = treatment]
raster_ex <- raster_dat[trial %in% ex_keys$trial & timepoint == min(raster_dat$timepoint)]
p_raster <- ggplot(raster_ex, aes(xmin = start_bin * BIN_S, xmax = (end_bin + 1) * BIN_S,
                                  ymin = 0, ymax = 1, fill = state)) +
  geom_rect() +
  facet_wrap(~ lbl, ncol = 1, strip.position = "left") +
  scale_fill_manual(values = STATE_FILL, name = "State") +
  scale_y_continuous(breaks = NULL) +
  labs(title = "A  Exercise-bout raster (excerpt)", x = "Time within interval (s)", y = NULL) +
  theme_bout + theme(strip.text.y.left = element_text(angle = 0, hjust = 1, size = 7))
ggsave(file.path(fig_dir, "A_bout_raster.png"), p_raster,
       width = 250, height = 210, units = "mm", dpi = 300, bg = "white")

# --- (B) Duration ECDF (survival function), pooled uncensored durations -----
ecdf_dat <- bout_inventory[alphabet == "binary" & !censored, .(state, duration_s, treatment)]
p_ecdf <- ggplot(ecdf_dat, aes(x = duration_s, color = treatment)) +
  stat_ecdf(aes(y = 1 - after_stat(y)), geom = "step") +
  facet_wrap(~ state) +
  scale_x_log10() + scale_y_log10() +
  scale_color_manual(values = TREAT_COLORS, name = NULL) +
  labs(title = "B  Bout-duration survival function (uncensored bouts)",
       x = "Duration (s, log scale)", y = "P(bout length > x)") +
  theme_bout
ggsave(file.path(fig_dir, "B_duration_ecdf.png"), p_ecdf,
       width = 200, height = 110, units = "mm", dpi = 300, bg = "white")

# --- (C) Strategy space: bout rate x median duration, per school ------------
strat_dat <- school_out[is.finite(bout_rate_flow) & bout_rate_flow > 0 &
                        is.finite(med_bout_flow) & med_bout_flow > 0]
p_strategy <- ggplot(strat_dat, aes(x = bout_rate_flow, y = med_bout_flow,
                                    size = commitment_index, color = treatment)) +
  geom_point(alpha = 0.8) +
  scale_x_log10(name = "Flow-bout rate (min⁻¹, log)") +
  scale_y_log10(name = "Median flow-bout duration (s, log)") +
  scale_color_manual(values = TREAT_COLORS, name = NULL) +
  scale_size_continuous(name = "Commitment index", range = c(2, 8)) +
  labs(title = "C  Frequency x duration strategy space (school means)") +
  theme_bout
ggsave(file.path(fig_dir, "C_strategy_space.png"), p_strategy,
       width = 180, height = 150, units = "mm", dpi = 300, bg = "white")

# --- (D) Burstiness by treatment ---------------------------------------------
burst_dat <- rbind(
  bout_metrics[, .(trial, treatment, state = "Flow", burst = burst_flow)],
  bout_metrics[, .(trial, treatment, state = "Calm", burst = burst_calm)]
)[is.finite(burst)]
burst_dat <- burst_dat[, .(burst = mean(burst, na.rm = TRUE)), by = .(trial, treatment, state)]
burst_stmt <- data.frame(state = c("Flow", "Calm"),
                         txt = c(bout_statement("burst_flow"), bout_statement("burst_calm")),
                         stringsAsFactors = FALSE)
p_burst <- ggplot(burst_dat, aes(x = pretty_treat(treatment), y = burst, color = treatment)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.12, height = 0, size = 2, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  facet_wrap(~ state) +
  ggtext::geom_richtext(data = burst_stmt, aes(x = 1.5, y = Inf, label = txt),
                        inherit.aes = FALSE, vjust = 1.05, size = 2.8,
                        fill = NA, label.color = NA, lineheight = 1.05) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.28))) +
  scale_color_manual(values = TREAT_COLORS, guide = "none") +
  labs(title = "D  Bout-duration burstiness (Kim-Jo corrected)", x = NULL,
       y = "Burstiness (0 = memoryless)") +
  theme_bout
ggsave(file.path(fig_dir, "D_burstiness.png"), p_burst,
       width = 200, height = 120, units = "mm", dpi = 300, bg = "white")

# --- (E) Observed vs null (shuffle | Markov), max_flow_bout_s ----------------
null_plot_dat <- bout_null_per_session[statistic == "max_flow_bout_s"]
null_plot_dat[, sess_lbl := factor(paste(trial, timepoint),
                                   levels = unique(paste(trial, timepoint)[order(obs)]))]
p_null <- ggplot(null_plot_dat, aes(y = sess_lbl, color = treatment)) +
  geom_linerange(aes(xmin = null_q025 + 1e-6, xmax = null_q975 + 1e-6)) +
  geom_point(aes(x = obs + 1e-6), size = 1.6) +
  facet_wrap(~ null_type, labeller = as_labeller(c(shuffle = "Shuffle null", markov = "Markov null"))) +
  scale_x_log10(name = "Longest flow bout (s, log)") +
  scale_color_manual(values = TREAT_COLORS, name = NULL) +
  labs(title = "E  Observed longest flow bout vs null 2.5-97.5% range", y = NULL) +
  theme_bout + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
ggsave(file.path(fig_dir, "E_observed_vs_null.png"), p_null,
       width = 220, height = 180, units = "mm", dpi = 300, bg = "white")

# --- (F) ICC forest -----------------------------------------------------------
icc_plot_dat <- bout_icc[metric %in% tier1_metrics]
icc_plot_dat[, grp := paste(subset, ifelse(adjusted, "adjusted", "unadjusted"))]
p_icc <- ggplot(icc_plot_dat, aes(x = icc, y = metric, color = grp)) +
  geom_pointrange(aes(xmin = icc_lo, xmax = icc_hi), position = position_dodge(width = 0.6)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  labs(title = "F  Repeatability (ICC) across the three intervals",
       x = "Intraclass correlation (trial)", y = NULL, color = NULL) +
  theme_bout
ggsave(file.path(fig_dir, "F_icc_forest.png"), p_icc,
       width = 190, height = 170, units = "mm", dpi = 300, bg = "white")

# --- (G) Rank-order stability spaghetti --------------------------------------
spaghetti_metrics <- c("commitment_index", "max_flow_bout_s", "flow_centroid", "si_markov_max")
spaghetti_dat <- data.table::melt(bout_metrics[, c("trial", "timepoint", "treatment", spaghetti_metrics), with = FALSE],
                                  id.vars = c("trial", "timepoint", "treatment"),
                                  variable.name = "metric", value.name = "val")
rank_lab <- bout_rank[subset == "all" & metric %in% spaghetti_metrics,
                      .(metric, txt = sprintf("W = %s, &chi;&sup2;(%d) = %s, p = %s",
                                              .fmt_es(kendall_W), df, .fmt_F(chisq), .fmt_p(p)))]
p_spaghetti <- ggplot(spaghetti_dat[is.finite(val)], aes(x = timepoint, y = val, group = trial, color = treatment)) +
  geom_line(alpha = 0.6) + geom_point(size = 1.5) +
  facet_wrap(~ metric, scales = "free_y") +
  ggtext::geom_richtext(data = rank_lab, aes(x = 2, y = Inf, label = txt), inherit.aes = FALSE,
                        vjust = 1.1, size = 2.6, fill = NA, label.color = NA) +
  scale_color_manual(values = TREAT_COLORS, name = NULL) +
  scale_x_continuous(breaks = 1:3) +
  labs(title = "G  Rank-order stability across intervals", x = "Interval", y = NULL) +
  theme_bout
ggsave(file.path(fig_dir, "G_rank_stability.png"), p_spaghetti,
       width = 220, height = 170, units = "mm", dpi = 300, bg = "white")

# --- (H) Interval trajectory on the real clock (5-25/45-65/85-105 min) ------
TP_MIN <- c(`1` = 15, `2` = 55, `3` = 95)  # interval midpoints, minutes
traj_dat <- bout_metrics[, .(commitment_index, treatment, timepoint)]
traj_summ <- traj_dat[, .(m = mean(commitment_index, na.rm = TRUE),
                          se = stats::sd(commitment_index, na.rm = TRUE) / sqrt(.N)),
                      by = .(treatment, timepoint)]
traj_summ[, tmin := TP_MIN[as.character(timepoint)]]
trend_lab <- bout_trend[metric == "commitment_index" & subset == "pooled"]
trend_txt <- if (nrow(trend_lab))
  sprintf("*linear trend: t<sub>%s</sub> = %s, p = %s*", .fmt_df(trend_lab$df[1]),
          .fmt_F(trend_lab$t[1]), .fmt_p(trend_lab$p[1])) else ""
p_traj <- ggplot(traj_summ, aes(x = tmin, y = m, color = treatment)) +
  geom_line() + geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = m - se, ymax = m + se), width = 3) +
  scale_color_manual(values = TREAT_COLORS, name = NULL) +
  labs(title = "H  Commitment index across the trial", x = "Minutes into trial",
       y = "Commitment index (mean ± SE)", caption = trend_txt) +
  theme_bout
ggsave(file.path(fig_dir, "H_interval_trajectory.png"), p_traj,
       width = 220, height = 130, units = "mm", dpi = 300, bg = "white")

# --- (I) Self-selected intensity ---------------------------------------------
int_dat <- bin_occ[state_binary == "Flow" & !is.na(state_intensity)]
int_summ <- int_dat[, .N, by = .(treatment, timepoint, state_intensity)]
int_summ[, prop := N / sum(N), by = .(treatment, timepoint)]
p_intensity <- ggplot(int_summ, aes(x = factor(timepoint), y = prop, fill = state_intensity)) +
  geom_col(position = "stack") +
  facet_wrap(~ pretty_treat(treatment)) +
  scale_fill_manual(values = c(low = "#FDE0A9", medium = "#F0A860", high = "#B8410E"),
                    name = "Self-selected\nintensity") +
  labs(title = "I  Self-selected exercise intensity within the flow zone",
       x = "Interval", y = "Proportion of flow time") +
  theme_bout
ggsave(file.path(fig_dir, "I_intensity_selection.png"), p_intensity,
       width = 200, height = 140, units = "mm", dpi = 300, bg = "white")

# --- (J) Collinearity heatmap -------------------------------------------------
if (nrow(metric_collinearity)) {
  p_collin <- ggplot(metric_collinearity, aes(x = metric_seq, y = metric_new, fill = rho_spearman)) +
    geom_tile(aes(colour = flagged), linewidth = 0.8) +
    geom_text(aes(label = sprintf("%.2f", rho_spearman)), size = 2.6) +
    scale_fill_gradient2(low = "#08519C", mid = "white", high = "#B8410E", midpoint = 0,
                         limits = c(-1, 1), name = "Spearman\nrho") +
    scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = NA), guide = "none") +
    labs(title = "J  Collinearity with sequence-module metrics (school level)",
         x = NULL, y = NULL) +
    theme_bout + theme(axis.text.x = element_text(angle = 30, hjust = 1))
  ggsave(file.path(fig_dir, "J_collinearity_heat.png"), p_collin,
         width = 200, height = 160, units = "mm", dpi = 300, bg = "white")
}

# =============================================================================
# 14) RESULTS GRIDS (K = complete, L = significant only)
# =============================================================================
TERM_LAB <- c(treatment = "treatment", timepoint_f = "interval",
             `treatment:timepoint_f` = "treatment × interval")
grid_dat <- rbind(
  bout_trial[, .(level = "trial", tier, metric, term, F, df1, df2, p, eta2_p, tested, exclude_reason)],
  bout_trial_tp[, .(level = "trial_x_timepoint", tier, metric, term, F, df1, df2, p, eta2_p,
                    tested, exclude_reason)]
)
grid_dat[, term_lab := factor(TERM_LAB[term], levels = unname(TERM_LAB))]
grid_dat[, level_lab := factor(c(trial = "Level A — trial (n = 16)",
                                 trial_x_timepoint = "Level B — trial × interval (n = 48)")[level],
                               levels = c("Level A — trial (n = 16)",
                                         "Level B — trial × interval (n = 48)"))]
grid_dat[, sig := isTRUE(tested) & is.finite(p) & p < 0.05]
grid_dat[, cell_lab := ifelse(
  tested,
  sprintf("F<sub>%s,%s</sub> = %s<br>p = %s<br>&eta;<sup>2</sup><sub>p</sub> = %s",
          .fmt_df(df1), .fmt_df(df2), .fmt_Fstat(F), .fmt_p(p), .fmt_es2(eta2_p)),
  "not tested"
)]
grid_dat[, metric := factor(metric, levels = rev(intersect(metric_order, unique(metric))))]

.results_grid <- function(d, ttl, sub) {
  d <- droplevels(d)
  ggplot(d, aes(x = term_lab, y = metric, fill = eta2_p)) +
    geom_tile(aes(colour = sig), linewidth = 0.7, width = 0.97, height = 0.97) +
    ggtext::geom_richtext(aes(label = cell_lab), size = 2.4, fill = NA, label.color = NA,
                          lineheight = 1.15, label.padding = grid::unit(rep(0, 4), "pt")) +
    facet_grid(tier ~ level_lab, scales = "free", space = "free") +
    scale_fill_gradient(low = "#F7FBFF", high = "#08519C", limits = c(0, 1), name = expression(eta[p]^2),
                        na.value = "grey92") +
    scale_colour_manual(values = c(`TRUE` = "#B8410E", `FALSE` = "grey85"), guide = "none") +
    labs(title = ttl, subtitle = sub, x = NULL, y = NULL) +
    theme_bout +
    theme(axis.text.x = element_text(angle = 20, hjust = 1), panel.grid = element_blank(),
          panel.spacing = grid::unit(6, "pt"))
}
p_grid_all <- .results_grid(grid_dat, "K  Exercise-bout structure results — complete grid",
  paste0("All analysed metrics, both analysis levels. Fill and label give partial eta-squared; ",
         "red outline marks raw p < 0.05. 'not tested' = failed the coverage or collinearity gate."))
ggsave(file.path(fig_dir, "K_results_grid_complete.png"), p_grid_all,
       width = 330, height = 300, units = "mm", dpi = 300, bg = "white", limitsize = FALSE)

sig_dat <- grid_dat[sig == TRUE]
if (nrow(sig_dat)) {
  p_grid_sig <- .results_grid(sig_dat, "L  Exercise-bout structure results — significant only",
    "Terms reaching raw p < 0.05 among those that passed the coverage/collinearity gates.")
} else {
  p_grid_sig <- ggplot() + annotate("text", x = 0, y = 0,
    label = "No bout-structure term reached p < 0.05", size = 5) + theme_void()
}
ggsave(file.path(fig_dir, "L_results_grid_significant.png"), p_grid_sig,
       width = 300, height = 200, units = "mm", dpi = 300, bg = "white", limitsize = FALSE)
.w(grid_dat[, .(level, tier, metric, term, F, df1, df2, p, eta2_p, n_unit = NA_integer_,
               sig, tested, exclude_reason)], "bout_results_grid.csv")

# =============================================================================
# 15) CONSOLE SUMMARY
# =============================================================================
ts_msg("==================== BOUT STRUCTURE ANALYSIS DONE ====================")
ts_msg("Output dir: ", out_dir)
cat("\n--- Coverage: sessions with finite value, by treatment ---\n")
for (m in c(tier1_metrics, tier2_metrics, descriptive_metrics)) {
  if (!m %in% names(bout_metrics)) next
  n_ctrl <- sum(is.finite(bout_metrics[[m]]) & bout_metrics$treatment == "control")
  n_exer <- sum(is.finite(bout_metrics[[m]]) & bout_metrics$treatment == "exercise choice")
  cat(sprintf("  %-22s control=%2d/24  exercise=%2d/24\n", m, n_ctrl, n_exer))
}
cat("\n--- Zero-switch sessions (commitment_index == 1.0) ---\n")
print(bout_metrics[zero_switch == TRUE, .(trial, timepoint, treatment)])

cat("\n=== LEVEL A -- TRIAL LEVEL (n = 16 schools) ===\n")
print(bout_trial[, .(tier, metric, tested, n_ctrl, n_exer,
                     mean_ctrl = round(mean_ctrl, 3), mean_exer = round(mean_exer, 3),
                     F = round(F, 2), df1, df2 = round(df2, 1), p = round(p, 4),
                     eta2_p = round(eta2_p, 3))])

cat("\n=== LEVEL B -- TRIAL x TIMEPOINT (n = 48 sessions) ===\n")
print(bout_trial_tp[, .(tier, metric, term, tested, F = round(F, 2), df1,
                        df2 = round(df2, 1), p = round(p, 4), eta2_p = round(eta2_p, 3))])

cat("\n--- Monotone trend contrasts (linear component across intervals) ---\n")
print(bout_trend[, .(metric, subset, t = round(t, 2), df = round(df, 1), p = round(p, 4))])

cat("\n--- Repeatability (ICC), adjusted, pooled ---\n")
print(bout_icc[subset == "all" & adjusted == TRUE,
               .(metric, icc = round(icc, 3), icc_lo = round(icc_lo, 3), icc_hi = round(icc_hi, 3),
                 p_ranova = round(p_ranova, 4), singular)])

cat("\n--- Rank-order stability (Kendall's W), pooled ---\n")
print(bout_rank[subset == "all", .(metric, kendall_W = round(kendall_W, 3),
                                   p = round(p, 4), n_complete)])

cat("\n--- Sign-consistency tests, pooled ---\n")
print(bout_sign[subset == "all", .(pattern, k, n, prop = round(prop, 3), p_binom = round(p_binom, 4))])

cat("\n--- Null-model exceedance summary ---\n")
print(bout_null_summary[, .(null_type, statistic, treatment, n_defined, n_above_q975,
                            p_binom = round(p_binom, 4))])

cat("\n--- Collinearity gate: metrics demoted to descriptive ---\n")
if (length(demoted_metrics)) print(demoted_metrics) else cat("  none\n")

cat("\n--- motor_side sensitivity on intensity metrics ---\n")
print(intensity_motorside_check)

cat("\n--- Singular fits (Level A / Level B) ---\n")
sing_a <- bout_trial[isTRUE(singular)]
sing_b <- unique(bout_trial_tp[isTRUE(singular), .(metric, model)])
if (nrow(sing_a)) print(sing_a[, .(level = "trial", metric, model)])
if (nrow(sing_b)) print(sing_b[, .(level = "trial_x_timepoint", metric, model)])
if (!nrow(sing_a) && !nrow(sing_b)) cat("  none\n")

ts_msg("Figures written to: ", fig_dir)
