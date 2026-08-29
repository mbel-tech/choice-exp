# =============================================================================
# collinearity_vs_ALR.R
# =============================================================================
# Screens every sequence-analysis and bout-structure metric for collinearity
# against the PRIMARY occupancy endpoint (the ALR log-ratio models, section
# 3.2.1), at the school (trial) level.
#
# WHY THIS EXISTS
# ----------------
# The bout engine's own collinearity screen (bout_structure_analysis_choice_
# exp.R section 9, exported as metric_collinearity.csv / Supplementary Table
# S13) tests candidate bout metrics against the ESTABLISHED SEQUENCE metrics
# (switch_rate, entropy_rate, occ_Flow, dwell_Flow). It never tests anything
# against the ALR primary endpoint. That is a gap: the sequence metrics it
# uses as the reference are THEMSELVES 74-88% redundant with ALR(flow) at the
# school level (see collinearity_vs_ALR.csv, this script's own output), so a
# metric can clear the existing screen and still be almost entirely a
# re-description of the paper's headline result. This script runs that
# second, higher-level screen explicitly.
#
# NO MODEL IS FITTED, NO STATISTIC IS RECOMPUTED. Every number here is either
# read from the engines' own exported CSVs or a Pearson correlation computed
# directly on those exported values.
#
# THE TRIAL-ID CROSSWALK -- READ THIS BEFORE TRUSTING ANY OUTPUT
# -----------------------------------------------------------------------------
# The sequence/bout engines (SEQ_output/, BOUT_output/) and the collective/
# zone-preference engine (easy_scripts_dataset.csv, phys_trial_id) label the
# 16 physical schools with DIFFERENT integers. Both are the SAME partition of
# the 48 sessions, just permuted: e.g. SEQ/BOUT school "1" = easy_scripts
# school "2", SEQ/BOUT "3" = easy_scripts "4", "13"<->"14", "15"<->"16".
# Joining naively on the numeric trial/school label silently pairs 8 of the
# 16 schools with the WRONG partner and produces uncorrelated-looking nonsense
# (occ_Flow vs ALR(flow) comes out at r ~ 0.08 instead of ~0.88, despite both
# being direct functions of the same prop_flow quantity). The only safe join
# key is the SESSION-level trial_id, which both sides share. This script
# builds and asserts that crosswalk before computing anything, and asserts
# the occ_Flow / ALR(flow) sanity check at the end -- if that assertion ever
# fails, the join is broken again; do not trust the output.
#
# INPUT
#   easy_scripts_dataset.csv                             (zone/ALR proportions)
#   SEQ_output_*/seq_metrics_per_trial_{binary,graded}.csv
#   BOUT_output_*/bout_metrics_per_school.csv
#
# OUTPUT
#   collinearity_vs_ALR.csv  -- one row per (metric x ALR outcome): r, r2, p, n
#   Console summary: per-metric max |r| across the four ALR outcomes + gate
#   decision (|r| >= 0.70 -> descriptive, matching the project's existing
#   collinearity convention, Table S13 / section 2.8.7).
#
# STANDALONE: run with
#   "C:/Program Files/R/R-4.6.0/bin/Rscript.exe" --vanilla <this file>
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
})

.PIPE <- file.path(PROJECT_ROOT, "choice R pipeline")
GATE  <- 0.70   # |r| threshold; matches section 2.8.7 / Table S13's own gate

.latest_run <- function(kind = c("SEQ", "BOUT")) {
  kind <- match.arg(kind)
  root <- file.path(.PIPE, "output", paste0(kind, "_output"))
  ds <- list.dirs(root, recursive = FALSE)
  ds <- ds[grepl(paste0(kind, "_output_\\d{8}_\\d{6}$"), ds)]
  if (!length(ds)) stop("No ", kind, "_output_* run found under ", root)
  sort(ds, decreasing = TRUE)[1]
}

SEQD  <- .latest_run("SEQ")
BOUTD <- .latest_run("BOUT")
cat("SEQ  run:", SEQD,  "\n")
cat("BOUT run:", BOUTD, "\n")

easy_path <- file.path(.PIPE, "easy_scripts", "easy_scripts_dataset.csv")
easy <- read.csv(easy_path, stringsAsFactors = FALSE)
seq_bin_trial <- read.csv(file.path(SEQD, "seq_metrics_per_trial_binary.csv"),
                          stringsAsFactors = FALSE)

# -----------------------------------------------------------------------------
# 1) CROSSWALK: SEQ/BOUT `trial` <-> easy_scripts `phys_trial_id`, via the
#    session-level `trial_id` both sides share. Asserted, not assumed.
# -----------------------------------------------------------------------------
cw <- easy %>%
  distinct(trial_id, phys_trial_id) %>%
  inner_join(seq_bin_trial %>% distinct(trial_id, trial), by = "trial_id") %>%
  distinct(trial, phys_trial_id)

.chk <- easy %>% distinct(trial_id, phys_trial_id, treatment) %>%
  inner_join(seq_bin_trial %>% distinct(trial_id, trial, treatment) %>%
               rename(trt_seq = treatment), by = "trial_id")
n_disagree <- sum(.chk$treatment != .chk$trt_seq)

stopifnot(
  "crosswalk must cover exactly 16 schools" = nrow(cw) == 16,
  "crosswalk `trial` must be unique"        = !anyDuplicated(cw$trial),
  "crosswalk `phys_trial_id` must be unique" = !anyDuplicated(cw$phys_trial_id),
  "0/48 sessions may disagree on treatment across the two labelings" =
    n_disagree == 0
)
cat(sprintf("Crosswalk OK: 16 schools, 1:1, 0/%d treatment disagreements.\n",
            nrow(.chk)))

# -----------------------------------------------------------------------------
# 2) ALR PRIMARY ENDPOINTS at trial level, from the SAME trial-aggregated
#    proportions the manuscript's own ALR models use (p_*_agg columns).
# -----------------------------------------------------------------------------
zone <- easy %>%
  group_by(phys_trial_id, treatment) %>%
  summarise(across(c(p_main_flow_agg, p_main_calm_agg,
                     p_an_sec_high_agg, p_an_sec_medium_agg,
                     p_an_sec_low_agg, p_an_sec_calm_agg),
                   ~ if (all(is.na(.x))) NA_real_ else first(na.omit(.x))),
            .groups = "drop") %>%
  mutate(
    `ALR flow`   = log(p_main_flow_agg   / p_main_calm_agg),
    `ALR high`   = log(p_an_sec_high_agg   / p_an_sec_calm_agg),
    `ALR medium` = log(p_an_sec_medium_agg / p_an_sec_calm_agg),
    `ALR low`    = log(p_an_sec_low_agg    / p_an_sec_calm_agg)
  )
ALR_COLS <- c("ALR flow", "ALR high", "ALR medium", "ALR low")

# -----------------------------------------------------------------------------
# 3) CANDIDATE METRICS at trial level.
#    Sequence: mean across the 3 intervals per trial (both alphabets).
#    Bout: already one row per school.
# -----------------------------------------------------------------------------
.agg_seq <- function(fname, suffix = "") {
  d <- read.csv(file.path(SEQD, fname), stringsAsFactors = FALSE)
  drop_cols <- c("n_bins", "n_changes", "obs_min", "trial", "trial_id",
                 "timepoint", "treatment", "tank", "fish_density")
  keep <- setdiff(names(d), drop_cols)
  keep <- keep[!grepl("^n_from", keep)]
  out <- d %>%
    group_by(trial) %>%
    summarise(across(all_of(keep),
                     ~ if (all(is.na(.x))) NA_real_ else mean(.x, na.rm = TRUE)),
              .groups = "drop")
  if (nzchar(suffix)) names(out)[-1] <- paste0(names(out)[-1], suffix)
  out
}
seq_binary <- .agg_seq("seq_metrics_per_trial_binary.csv")
seq_graded <- .agg_seq("seq_metrics_per_trial_graded.csv", "_gr")
bout <- read.csv(file.path(BOUTD, "bout_metrics_per_school.csv"),
                 stringsAsFactors = FALSE) %>%
  select(-any_of(c("treatment", "n_sessions_used")))

dat <- cw %>%
  left_join(seq_binary, by = "trial") %>%
  left_join(seq_graded, by = "trial") %>%
  left_join(bout,       by = "trial") %>%
  left_join(zone %>% select(phys_trial_id, all_of(ALR_COLS)), by = "phys_trial_id")

METRICS <- setdiff(names(dat), c("trial", "phys_trial_id", ALR_COLS))

# -----------------------------------------------------------------------------
# 4) SANITY ASSERTION -- the join is trustworthy only if this passes.
#    occ_Flow and ALR(flow) are both direct functions of the same prop_flow;
#    a low value here means the crosswalk (step 1) failed.
# -----------------------------------------------------------------------------
sanity_r <- suppressWarnings(cor(dat$occ_Flow, dat$`ALR flow`, use = "complete.obs"))
cat(sprintf("\nSANITY: occ_Flow vs ALR(flow) r = %.3f (expect ~0.88)\n", sanity_r))
stopifnot("occ_Flow/ALR(flow) sanity check failed -- crosswalk is likely broken" =
            is.finite(sanity_r) && sanity_r > 0.70)

# -----------------------------------------------------------------------------
# 5) PEARSON r FOR EVERY (metric x ALR outcome) PAIR
# -----------------------------------------------------------------------------
res <- do.call(rbind, lapply(METRICS, function(m) do.call(rbind, lapply(ALR_COLS, function(zz) {
  x <- dat[[m]]; y <- dat[[zz]]
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 4 || sd(x[ok]) == 0 || sd(y[ok]) == 0)
    return(data.frame(metric = m, alr_outcome = zz, n = sum(ok),
                       r = NA_real_, r2 = NA_real_, p = NA_real_))
  ct <- suppressWarnings(cor.test(x[ok], y[ok]))
  data.frame(metric = m, alr_outcome = zz, n = sum(ok),
             r = unname(ct$estimate), r2 = unname(ct$estimate)^2, p = ct$p.value)
}))))

summary_tbl <- res %>%
  filter(!is.na(r)) %>%
  group_by(metric) %>%
  slice_max(abs(r), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(metric, worst_alr_outcome = alr_outcome, max_abs_r = abs(r),
            r2_pct = 100 * r2, n,
            gate = ifelse(max_abs_r >= GATE, "descriptive", "independent")) %>%
  arrange(desc(max_abs_r))

OUT <- file.path(.PIPE, "easy_scripts", "standalone", "collinearity_vs_ALR.csv")
dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)
write.csv(res, OUT, row.names = FALSE)
write.csv(summary_tbl,
          file.path(dirname(OUT), "collinearity_vs_ALR_summary.csv"),
          row.names = FALSE)

cat(sprintf("\nWrote %s (%d rows) and collinearity_vs_ALR_summary.csv (%d metrics)\n",
            OUT, nrow(res), nrow(summary_tbl)))
cat(sprintf("\nGate: |r| >= %.2f vs any ALR outcome -> descriptive (n = %d schools)\n\n",
            GATE, nrow(dat)))
print(summary_tbl, n = nrow(summary_tbl))
