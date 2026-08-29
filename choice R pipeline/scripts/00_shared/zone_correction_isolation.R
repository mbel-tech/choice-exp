# =============================================================================
# zone_correction_isolation.R
# =============================================================================
# Isolates the effect of the 2026-08 zone-assignment correction on the reported
# statistics (R2.6 memo 5.7).
#
# The problem the memo names: the analysis code was revised in the same interval
# as the correction, so comparing the old published statistics with the new ones
# confounds the two. The clean test is to put the RETAINED pre-correction
# frame-level dataset through the CURRENT analysis, and both datasets are on
# disk with the same 7,156,680 rows from the same 48 sessions.
#
# Derivation used here is STEP2's, transcribed:
#   per frame, over rows with a finite interpolated position,
#     n_total, n_in_flow, n_in_calm, n_in_high, n_in_medium, n_in_low,
#     n_in_calm_sec
#   dt_row = diff(time) within (trial_id, timepoint)
#   prop_time_in_<z> = weighted.mean(n_in_<z> / max(n_total,1), w = dt_row)
#   sub-zones divided by area units (high 10, medium 18, low 13, calm 42) and
#     rescaled to sum to 1                       -> prop_time_ac_*
#   logit_flow = log(p/(1-p)) on prop_time_in_flow squeezed to [1e-4, 1-1e-4]
#   lr_<z>     = log(prop_time_ac_<z> / prop_time_ac_calm), Haldane floor 1e-4
#
# The transcription is not taken on trust: applied to the POST-correction file
# it must reproduce the pipeline's own easy_scripts_dataset.csv to 1e-8, and the
# refitted models must reproduce the published F and p. Only then is the
# PRE-correction run interpretable.

suppressMessages({library(data.table); library(lme4); library(lmerTest)})
options(contrasts = c("contr.sum", "contr.poly"))
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

PIPE <- file.path(PROJECT_ROOT, "choice R pipeline")
RUN  <- file.path(PIPE, "output/STEP5_stats/STEP5_stats_20260818_144007")
PRE  <- file.path(PIPE, "output/STEP1_output/STEP1_output_20260426_230846/master_fish_by_frame.csv")
POST <- file.path(PIPE, "output/STEP1_output/STEP1_output_20260807_162823/master_fish_by_frame.csv")

AREA <- c(high = 10, medium = 18, low = 13, calm = 42)
EPS  <- 1e-4

derive <- function(path, label) {
  msg("reading ", label, " ...")
  d <- fread(path, select = c("trial_id","frame","time","timepoint",
                              "main_zone","sec_zone","x_interp","y_interp"),
             showProgress = FALSE)
  msg("  ", nrow(d), " rows")
  d <- d[is.finite(x_interp) & is.finite(y_interp)]
  fc <- d[, .(n_total = .N,
              n_in_flow     = sum(main_zone == "flow",   na.rm = TRUE),
              n_in_calm     = sum(main_zone == "calm",   na.rm = TRUE),
              n_in_high     = sum(sec_zone  == "high",   na.rm = TRUE),
              n_in_medium   = sum(sec_zone  == "medium", na.rm = TRUE),
              n_in_low      = sum(sec_zone  == "low",    na.rm = TRUE),
              n_in_calm_sec = sum(sec_zone  == "calm",   na.rm = TRUE),
              time = time[1]),
          by = .(trial_id, timepoint, frame)]
  setkey(fc, trial_id, timepoint, frame)
  fc[, dt_row := c(0, diff(time)), by = .(trial_id, timepoint)]
  fc[!is.finite(dt_row) | dt_row < 0, dt_row := 0]
  wm <- function(num, w) stats::weighted.mean(num / pmax(fc$n_total, 1L), w = w, na.rm = TRUE)
  out <- fc[, .(
      prop_time_in_flow     = stats::weighted.mean(n_in_flow     / pmax(n_total,1L), w = dt_row, na.rm = TRUE),
      prop_time_in_calm     = stats::weighted.mean(n_in_calm     / pmax(n_total,1L), w = dt_row, na.rm = TRUE),
      prop_time_in_high     = stats::weighted.mean(n_in_high     / pmax(n_total,1L), w = dt_row, na.rm = TRUE),
      prop_time_in_medium   = stats::weighted.mean(n_in_medium   / pmax(n_total,1L), w = dt_row, na.rm = TRUE),
      prop_time_in_low      = stats::weighted.mean(n_in_low      / pmax(n_total,1L), w = dt_row, na.rm = TRUE),
      prop_time_in_calm_sec = stats::weighted.mean(n_in_calm_sec / pmax(n_total,1L), w = dt_row, na.rm = TRUE)),
    by = .(trial_id, timepoint)]
  out[, `:=`(nh = prop_time_in_high/AREA[["high"]], nm = prop_time_in_medium/AREA[["medium"]],
             nl = prop_time_in_low/AREA[["low"]],   nc = prop_time_in_calm_sec/AREA[["calm"]])]
  out[, nt := nh + nm + nl + nc]
  out[, `:=`(prop_time_ac_high = nh/nt, prop_time_ac_medium = nm/nt,
             prop_time_ac_low  = nl/nt, prop_time_ac_calm   = nc/nt)]
  out[, `:=`(logit_flow = { p <- pmin(pmax(prop_time_in_flow, EPS), 1-EPS); log(p/(1-p)) },
             lr_high   = log(pmax(prop_time_ac_high,  EPS)/pmax(prop_time_ac_calm, EPS)),
             lr_medium = log(pmax(prop_time_ac_medium,EPS)/pmax(prop_time_ac_calm, EPS)),
             lr_low    = log(pmax(prop_time_ac_low,   EPS)/pmax(prop_time_ac_calm, EPS)))]
  out[, dataset := label]
  out[, c("nh","nm","nl","nc","nt") := NULL]
  out[]
}

post <- derive(POST, "post-correction (reported)")
pre  <- derive(PRE,  "pre-correction")

# --- validation 1: does the transcription reproduce the pipeline's own file? --
ez <- fread(file.path(PIPE, "easy_scripts/easy_scripts_dataset.csv"))
chk <- merge(post[, .(trial_id, timepoint, logit_flow, lr_high, lr_medium, lr_low)],
             ez[, .(trial_id, timepoint, ez_logit = logit_flow, ez_high = lr_high,
                    ez_med = lr_medium, ez_low = lr_low)],
             by = c("trial_id","timepoint"))
dev <- c(logit_flow = max(abs(chk$logit_flow - chk$ez_logit)),
         lr_high    = max(abs(chk$lr_high    - chk$ez_high)),
         lr_medium  = max(abs(chk$lr_medium  - chk$ez_med)),
         lr_low     = max(abs(chk$lr_low     - chk$ez_low)))
msg("max |deviation| from the pipeline's own dataset:")
print(signif(dev, 3))
TRANSCRIPTION_EXACT <- all(dev < 1e-8)
msg("transcription exact to 1e-8: ", TRANSCRIPTION_EXACT)

# --- refit the reported models on each dataset ------------------------------
meta <- unique(fread(file.path(RUN, "analysis_ready.csv"))[
  , .(trial_id, timepoint, treatment, tank, fish_density, phys_trial_id)])

fit_one <- function(dat, resp) {
  d <- merge(dat, meta, by = c("trial_id","timepoint"))
  d[, `:=`(treatment = factor(treatment, levels = c("control","exercise choice")),
           timepoint_f = factor(timepoint), trial_key = factor(phys_trial_id))]
  d[, y_t := as.numeric(get(resp))]
  m <- suppressMessages(suppressWarnings(
    lmerTest::lmer(y_t ~ treatment * timepoint_f + (1 | trial_key), data = d, REML = TRUE)))
  a <- suppressMessages(anova(m, type = 3, ddf = "Kenward-Roger"))
  em <- suppressMessages(emmeans::emmeans(m, ~ treatment))
  ct <- as.data.frame(emmeans::contrast(em, "pairwise"))
  list(F = a["treatment","F value"], df1 = a["treatment","NumDF"],
       df2 = a["treatment","DenDF"], p = a["treatment","Pr(>F)"],
       est = ct$estimate[1])
}

RESP <- c(`ALR(flow vs calm)` = "logit_flow", `ALR(high vs calm)` = "lr_high",
          `ALR(medium vs calm)` = "lr_medium", `ALR(low vs calm)` = "lr_low")
PUB  <- c(logit_flow = 17.503896, lr_high = 7.529668, lr_medium = 5.275043, lr_low = 2.178357)

cmp <- rbindlist(lapply(names(RESP), function(lab) {
  r <- RESP[[lab]]
  a <- fit_one(post, r); b <- fit_one(pre, r)
  data.table(outcome = lab, response = r,
             F_post = a$F, df2_post = a$df2, p_post = a$p, est_post = a$est,
             F_pre  = b$F, df2_pre  = b$df2, p_pre  = b$p, est_pre  = b$est,
             F_published = PUB[[r]],
             post_reproduces_published = abs(a$F - PUB[[r]]) <= 5e-3 * PUB[[r]],
             sig_post = a$p < 0.05, sig_pre = b$p < 0.05,
             conclusion_changed = (a$p < 0.05) != (b$p < 0.05) || sign(a$est) != sign(b$est))
}))
print(cmp)

# --- frame-level occupancy by configuration, the attributable quantity -------
occ <- rbindlist(lapply(list(post, pre), function(x) {
  d <- merge(x, unique(meta[, .(trial_id, tank)]), by = "trial_id")
  d[, motor := ifelse(tank %in% c(27, 29), "FT", "FD")]
  d[, .(prop_flow = mean(prop_time_in_flow), prop_calm = mean(prop_time_in_calm)),
    by = .(dataset, motor)]
}))
occ_w <- dcast(occ, dataset ~ motor, value.var = c("prop_flow","prop_calm"))
occ_w[, `:=`(disagreement_flow_pp = 100*abs(prop_flow_FT - prop_flow_FD),
             disagreement_calm_pp = 100*abs(prop_calm_FT - prop_calm_FD))]
print(occ_w)

fwrite(cmp,   file.path(RUN, "zone_correction_isolation.csv"))
fwrite(occ_w, file.path(RUN, "zone_correction_occupancy_by_config.csv"))
fwrite(rbind(post, pre), file.path(RUN, "zone_correction_derived_occupancy.csv"))
msg("written; transcription_exact=", TRANSCRIPTION_EXACT,
    " any_conclusion_changed=", any(cmp$conclusion_changed))
