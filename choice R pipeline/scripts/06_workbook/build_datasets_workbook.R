# =============================================================================
# build_datasets_workbook.R
# -----------------------------------------------------------------------------
# Builds: D:/CHOICE R SCRIPTS/choice R pipeline/
#         Behavioural and neuroendocrine correlates datasets.xlsx
#
# Workbook structure (7 sheets) -- PREFERENCE EXPERIMENT ONLY:
#   00_README
#   Preference_exp_behavior          <- 48 rows (16 trials x 3 intervals), interval level
#   Preference_exp_behavior_legend
#   preference_exp_cortisol          <- 1 row per sampled fish
#   preference_exp_cortisol_legend
#   preference_exp_monoamines        <- 1 row per fish, 24 analyte x area cols
#   preferenceexp_monoamines_legend
#
# CHANGED 2026-08-29
#   - The observational (B1) study is removed. It was carried here as two
#     verbatim sheets copied out of "B1 recheck stuff/choice exp datasets.xlsx";
#     no script in this pipeline analyses it, and it is out of scope for the
#     deposit this workbook backs.
#   - Five reported outcomes were missing: flow bouts per minute, longest flow
#     bout, longest calm bout, mean flow bout duration and mean calm bout
#     duration (outcomes 6-10 of the fourteen in
#     scripts/00_shared/export_tukey_interval_module.R). They live in the BOUT
#     and SEQ engines, which this builder never read. Both are now joined in.
#   - zone_flux_per_session was being DROPPED while switches_per_session was
#     kept. That is backwards: the reported "flow-calm crossings" outcome
#     (F(1,14) = 45.00, p = 9.97e-06, flux_timepoint/anova.csv) is the flux
#     column. Both are now kept, correctly named and distinguished.
#   - Noradrenaline dropped: 24 analyte x area cells, not 28.
#   - Every column is per interval. No trial-aggregate (_agg) column survives.
#   - The within-trial interval is called `interval` everywhere in the workbook,
#     the deposit and this script. The source frames spell it `timepoint`,
#     `timepoint_f` and `timepoint_num`; all three are renamed on read, so the
#     old spelling survives only where a source header is literally being read.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(readxl)
  library(openxlsx)
  library(tibble)
})

# ---- Paths ------------------------------------------------------------------
.WB_ROOT     <- file.path(PROJECT_ROOT, "choice R pipeline")
.WB_OUT_PATH <- file.path(.WB_ROOT, "Behavioural and neuroendocrine correlates datasets.xlsx")
.CSV_BEH     <- file.path(.WB_ROOT, "easy_scripts/easy_scripts_dataset.csv")
.CSV_ENDO    <- file.path(.WB_ROOT, "easy_scripts/easy_scripts_endo_dataset.csv")
# BOUT and SEQ engines: outcomes 6-8 and 9-10 of the reported fourteen. Both are
# 48 rows keyed (tank, fish_density, interval) -- NOT on `trial`, whose values
# in these two files are the chronological trial_seq, not phys_trial_id.
.CSV_BOUT    <- file.path(.WB_ROOT, "output/BOUT_output/BOUT_output_20260816_134054/bout_metrics_per_session.csv")
.CSV_SEQ     <- file.path(.WB_ROOT, "output/SEQ_output/SEQ_output_20260816_132132/seq_metrics_per_trial_binary.csv")
.RUN_STEP5   <- file.path(.WB_ROOT, "output/STEP5_stats/STEP5_stats_20260818_144007")
.CSV_TRIALLOG<- file.path(.RUN_STEP5, "trial_log_16.csv")
# Upstream of the endocrine CSV. Carries all 80 sampled fish (the endo CSV keeps
# only the 78 with a cortisol value) plus the fish-level metadata -- weight,
# length, video ID, batch -- that the endo CSV drops.
.CSV_CORT_CLEAN <- file.path(PROJECT_ROOT, "choice exp for claude mono and cortisol/claude output REVISED/b3_cortisol_clean.csv")
.MANIFEST_R  <- file.path(.WB_ROOT, "scripts/04_reporting/manifest_indicators.R")
.MAP_R       <- file.path(.WB_ROOT, "scripts/06_workbook/indicator_column_map.R")
.LEGEND_R    <- file.path(.WB_ROOT, "scripts/06_workbook/legend_descriptions.R")

# ---- Source dependencies ----------------------------------------------------
.log <- function(...) cat("[workbook] ", ..., "\n", sep = "")

source(.MAP_R)        # INDICATOR_MAP_BEH, INDICATOR_MAP_ENDO, validate_indicator_map
source(.LEGEND_R)     # METADATA_DESCRIPTIONS, build_indicator_descriptions, describe_src_column
source(.MANIFEST_R)   # MANIFEST  (the 74-row tibble)

stopifnot(exists("MANIFEST"),
          exists("INDICATOR_MAP_BEH"), exists("INDICATOR_MAP_ENDO"),
          exists("METADATA_DESCRIPTIONS"))

# Split endo map into cortisol and monoamine sub-maps (for legend building)
INDICATOR_MAP_CORT <- INDICATOR_MAP_ENDO[INDICATOR_MAP_ENDO$wb_col == "cort_plasma", ,
                                          drop = FALSE]
INDICATOR_MAP_MONO <- INDICATOR_MAP_ENDO[INDICATOR_MAP_ENDO$wb_col != "cort_plasma", ,
                                          drop = FALSE]

# Validate map <-> manifest BEFORE building anything
validate_indicator_map(MANIFEST, INDICATOR_MAP_BEH, INDICATOR_MAP_ENDO, verbose = TRUE)

# Raw-only behaviour map: one column per indicator, no aggregated/broadcast versions.
# Rule: drop any wb_col ending in "_agg" (covers broadcast_mean, direct-agg Jacobs D,
# representative _agg, and missing _agg variants).
INDICATOR_MAP_BEH_RAW <- INDICATOR_MAP_BEH[!endsWith(INDICATOR_MAP_BEH$wb_col, "_agg"), ,
                                            drop = FALSE]
.log("Indicator map (raw): ", nrow(INDICATOR_MAP_BEH_RAW),
     " behaviour cols (dropped ",
     nrow(INDICATOR_MAP_BEH) - nrow(INDICATOR_MAP_BEH_RAW), " _agg cols)")

.IND_DESC <- build_indicator_descriptions(INDICATOR_MAP_BEH, INDICATOR_MAP_ENDO)

# ---- Build behaviour_wide ---------------------------------------------------
.log("Reading behaviour CSV: ", .CSV_BEH)
beh_raw <- readr::read_csv(.CSV_BEH, show_col_types = FALSE)
.log("  rows=", nrow(beh_raw), " cols=", ncol(beh_raw))

# Defensive trial_date format (idempotent: handles both ISO and DD.MM.YYYY)
if ("trial_date" %in% names(beh_raw)) {
  td <- tryCatch(as.Date(beh_raw$trial_date,
                         tryFormats = c("%Y-%m-%d", "%d.%m.%Y")),
                 error = function(e) NULL)
  if (!is.null(td) && !all(is.na(td)))
    beh_raw$trial_date <- format(td, "%d.%m.%Y")
}

# phys_trial: rename phys_trial_id -> phys_trial if needed
if (!"phys_trial" %in% names(beh_raw) && "phys_trial_id" %in% names(beh_raw))
  beh_raw$phys_trial <- as.character(beh_raw$phys_trial_id)

# ---- One name for the within-trial interval ---------------------------------
# easy_scripts_dataset.csv carries the same variable three times, as `timepoint`,
# `timepoint_f` and `timepoint_num`. Collapse them here so that everything below
# -- and everything in the workbook and the deposit -- says `interval`.
stopifnot(all(beh_raw$timepoint == as.integer(as.character(beh_raw$timepoint_f))),
          all(beh_raw$timepoint == beh_raw$timepoint_num))
beh_raw$interval <- as.integer(beh_raw$timepoint)
beh_raw[, c("timepoint", "timepoint_f", "timepoint_num")] <- NULL

# ---- Join the BOUT and SEQ engines ------------------------------------------
# The join key is (tank, fish_density, interval). It is NOT `trial`: in both
# engines that column holds the chronological trial_seq (1-16 in running order),
# which differs from phys_trial_id -- trials 1<->2, 13<->14 and 15<->16 are
# swapped. Joining on it would silently mislabel six trials.
# Both engines spell the interval `timepoint` in their own headers; it is
# renamed to `interval` on read, as above.
.join_engine <- function(beh, path, cols, tag) {
  if (!file.exists(path)) stop("[workbook] ", tag, " source not found: ", path, call. = FALSE)
  d <- readr::read_csv(path, show_col_types = FALSE)
  miss <- setdiff(cols, names(d))
  if (length(miss)) stop("[workbook] ", tag, " missing column(s): ",
                         paste(miss, collapse = ", "), call. = FALSE)
  d <- d[, c("tank", "fish_density", "timepoint", cols), drop = FALSE]
  names(d)[names(d) == "timepoint"] <- "interval"   # source header -> our name
  d$tank         <- as.numeric(d$tank)
  d$fish_density <- as.numeric(d$fish_density)
  d$interval     <- as.integer(d$interval)
  if (anyDuplicated(d[, c("tank", "fish_density", "interval")]))
    stop("[workbook] ", tag, ": join key is not unique", call. = FALSE)
  n0  <- nrow(beh)
  out <- dplyr::left_join(beh, d, by = c("tank", "fish_density", "interval"))
  if (nrow(out) != n0)
    stop("[workbook] ", tag, " join changed row count: ", n0, " -> ", nrow(out), call. = FALSE)
  unmatched <- sum(!Reduce(`|`, lapply(cols, function(k) !is.na(out[[k]]))))
  .log("  joined ", tag, ": ", length(cols), " cols, ",
       n0 - unmatched, "/", n0, " rows matched")
  out
}
beh_raw <- .join_engine(beh_raw, .CSV_BOUT,
                        c("bout_rate_flow", "max_flow_bout_s", "max_calm_bout_s"), "BOUT")
beh_raw <- .join_engine(beh_raw, .CSV_SEQ,
                        c("dwell_Flow", "dwell_Calm"), "SEQ")

# ---- trial_seq, from the design log -----------------------------------------
# Deposited alongside phys_trial because the two are different orderings of the
# same 16 trials and the difference is easy to misread.
.tl <- readr::read_csv(.CSV_TRIALLOG, show_col_types = FALSE)
beh_raw <- dplyr::left_join(
  beh_raw,
  dplyr::transmute(.tl, tank = as.numeric(tank), fish_density = as.numeric(n_fish),
                   trial_seq = as.integer(trial_seq)),
  by = c("tank", "fish_density"))
stopifnot(nrow(beh_raw) == 48L, !anyNA(beh_raw$trial_seq))

# switches_per_session is a count that arrives carrying float noise from a
# rate x duration round-trip (23.0012, 40.9997, 462.9969...). Restore the count.
if ("switches_per_session" %in% names(beh_raw)) {
  .sw <- suppressWarnings(as.numeric(beh_raw$switches_per_session))
  stopifnot(max(abs(.sw - round(.sw)), na.rm = TRUE) < 0.05)
  beh_raw$switches_per_session <- round(.sw)
}

.beh_meta_cols <- c("trial_seq", "interval", "trial_date", "treatment", "tank",
                    "fish_density", "motor_side")

# Build one indicator column per INDICATOR_MAP_BEH row
build_indicator_col <- function(beh, map_row) {
  es_col <- map_row$es_col
  derive <- map_row$derive
  if (derive == "missing" || is.na(es_col) || !nzchar(es_col))
    return(rep(NA_real_, nrow(beh)))
  if (!es_col %in% names(beh))
    return(rep(NA_real_, nrow(beh)))
  v <- beh[[es_col]]
  if (derive == "broadcast_mean") {
    if (!"phys_trial" %in% names(beh)) return(rep(NA_real_, nrow(beh)))
    pt  <- as.character(beh$phys_trial)
    agg <- tapply(suppressWarnings(as.numeric(v)), pt, mean, na.rm = TRUE)
    out <- unname(agg[pt])
    out[is.nan(out)] <- NA_real_
    return(out)
  }
  suppressWarnings(as.numeric(v))
}

ind_cols <- vector("list", nrow(INDICATOR_MAP_BEH_RAW))
names(ind_cols) <- INDICATOR_MAP_BEH_RAW$wb_col
for (i in seq_len(nrow(INDICATOR_MAP_BEH_RAW)))
  ind_cols[[i]] <- build_indicator_col(beh_raw, INDICATOR_MAP_BEH_RAW[i, ])
ind_df <- as.data.frame(ind_cols, check.names = FALSE)

# ---- Context columns -------------------------------------------------------
# Raw occupancy and its coordinate transforms. These are model INPUTS, not model
# outputs, so they have no manifest row and are taken straight from the source
# frame. Depositing them is what lets a reader recompute every ALR, CLR and
# Jacobs' D rather than having to trust ours.
#   prop_flow / prop_calm      main-zone occupancy, uncorrected. They do not sum
#                              to 1: the remainder is time not assignable to
#                              either zone (fish untracked or between polygons).
#   prop_*_ac                  sub-zone composition AFTER the area correction
#                              applied once in STEP2 section 6b; sums to 1.
#                              alr_high = log(prop_high_ac / prop_calm_ac) with a
#                              Haldane floor of 1e-4 -- verified exact.
#   clr_*                      centred log-ratio of that same 4-part composition.
.BEH_CONTEXT_MAP <- c(
  obs_seconds    = "obs_seconds",
  n_frames       = "n_frames",
  prop_flow      = "prop_flow",
  prop_calm      = "prop_calm",
  prop_high_ac   = "prop_high",
  prop_medium_ac = "prop_medium",
  prop_low_ac    = "prop_low",
  prop_calm_ac   = "prop_calm_sec",
  clr_high       = "clr_sec_high_tp",
  clr_medium     = "clr_sec_medium_tp",
  clr_low        = "clr_sec_low_tp",
  clr_calm       = "clr_sec_calm_tp"
)
.ctx_missing <- setdiff(unname(.BEH_CONTEXT_MAP), names(beh_raw))
if (length(.ctx_missing))
  stop("[workbook] context column(s) absent from source: ",
       paste(.ctx_missing, collapse = ", "), call. = FALSE)
ctx_df <- as.data.frame(
  lapply(.BEH_CONTEXT_MAP, function(cc) suppressWarnings(as.numeric(beh_raw[[cc]]))),
  check.names = FALSE)
names(ctx_df) <- names(.BEH_CONTEXT_MAP)

behaviour_wide <- cbind(
  data.frame(phys_trial = as.integer(beh_raw$phys_trial)),
  beh_raw[, intersect(.beh_meta_cols, names(beh_raw)), drop = FALSE],
  ctx_df,
  ind_df
)
# Numeric sort. The previous character sort produced 1, 10, 11, ... 16, 2, 3.
behaviour_wide <- behaviour_wide[order(behaviour_wide$phys_trial,
                                       behaviour_wide$interval), , drop = FALSE]
rownames(behaviour_wide) <- NULL

# ---- Post-process behaviour_wide: drop, rename, group, reorder --------------

# Columns dropped entirely.
#   prop_active            all 48 values NA -- bl_cm_trial was never computed, so
#                          the model was skipped and no result exists.
#   fish_fish_nnd,
#   centroid_distance,
#   turning_rate           derive = "missing"; all-NA placeholder columns.
#   subzone_alr_pairs_tp   derive = "missing"; all-NA.
#   main_cells_flow_tp     bit-identical to the prop_flow context column.
#   sub_cells_high_clr_tp  bit-identical to the clr_high context column.
# NOTE: zone_flux_per_session is NOT dropped any more -- it is the reported
# flow-calm crossings outcome. See the header note.
.BEH_DROP <- c("prop_active", "centroid_distance", "turning_rate", "fish_fish_nnd",
               "subzone_alr_pairs_tp", "main_cells_flow_tp", "sub_cells_high_clr_tp")
behaviour_wide <- behaviour_wide[, !names(behaviour_wide) %in% .BEH_DROP, drop = FALSE]

# Rename lookup (new name -> old name in indicator map / description tables).
# logit_flow and lr_* are renamed because they ARE the additive log-ratios the
# Results report; the old names invited the reader to treat alr_flow as a plain
# logit and lr_high as something other than ALR(high vs calm).
.BEH_RENAME_MAP <- c(alr_flow                = "logit_flow",
                     alr_high                = "lr_high",
                     alr_medium              = "lr_medium",
                     alr_low                 = "lr_low",
                     crossings_per_session   = "zone_flux_per_session",
                     bouts_per_min           = "bouts_per_min",
                     longest_flow_bout_s     = "longest_flow_bout_s",
                     longest_calm_bout_s     = "longest_calm_bout_s",
                     mean_flow_bout_s        = "mean_flow_bout_s",
                     mean_calm_bout_s        = "mean_calm_bout_s",
                     mean_school_area_cm2    = "mean_hull_area_cm2",
                     mean_school_speed_cm_s  = "mean_centroid_spd_cm")
for (.new in names(.BEH_RENAME_MAP)) {
  .old <- .BEH_RENAME_MAP[[.new]]
  if (.old != .new && .old %in% names(behaviour_wide))
    names(behaviour_wide)[names(behaviour_wide) == .old] <- .new
}

# Column groups -- defines both order and header colour
.BEH_META       <- c("phys_trial", "trial_seq", "interval", "trial_date", "treatment",
                     "tank", "fish_density", "motor_side", "obs_seconds", "n_frames")
# 2026-08-18: mean_polarisation removed. It is not analysed anywhere -- no such
# column in easy_scripts_dataset.csv, absent from Figure 6 and from the Results
# text -- so cross-validation Layer 1 could only ever report it as
# FAIL-ES-COLUMN-MISSING: a column in the distributed workbook with nothing
# behind it. STEP2b still computes it if it is ever wanted.
.BEH_FLOW       <- c("alr_flow", "alr_high", "alr_medium", "alr_low",
                     "prop_flow", "prop_calm",
                     "prop_high_ac", "prop_medium_ac", "prop_low_ac", "prop_calm_ac",
                     "clr_high", "clr_medium", "clr_low", "clr_calm")
.BEH_ENGAGE     <- c("crossings_per_session", "bouts_per_min",
                     "longest_flow_bout_s", "longest_calm_bout_s",
                     "mean_flow_bout_s", "mean_calm_bout_s",
                     "switches_per_session")
.BEH_COLLECTIVE <- c("mean_nnd_cm", "mean_iid_cm",
                     "mean_school_area_cm2", "mean_school_speed_cm_s")
.BEH_JACOBS     <- c("jacobs_d_flow", "jacobs_d_high", "jacobs_d_medium", "jacobs_d_low")

.beh_canonical_order <- c(.BEH_META, .BEH_FLOW, .BEH_ENGAGE, .BEH_COLLECTIVE, .BEH_JACOBS)
.beh_present <- intersect(.beh_canonical_order, names(behaviour_wide))
.beh_extra   <- setdiff(names(behaviour_wide), .beh_present)   # any stragglers
if (length(.beh_extra))
  stop("[workbook] behaviour_wide: ", length(.beh_extra),
       " column(s) outside the canonical order -- classify or drop them: ",
       paste(.beh_extra, collapse = ", "), call. = FALSE)
behaviour_wide <- behaviour_wide[, .beh_present, drop = FALSE]

# The eleven STEP5/BOUT/SEQ outcomes reported at interval level, by workbook name.
.BEH_REPORTED <- c("alr_flow", "alr_high", "alr_medium", "alr_low",
                   "crossings_per_session", "bouts_per_min",
                   "longest_flow_bout_s", "longest_calm_bout_s",
                   "mean_flow_bout_s", "mean_calm_bout_s",
                   "mean_nnd_cm", "mean_iid_cm",
                   "mean_school_area_cm2", "mean_school_speed_cm_s")

.BEH_FINAL_INDICATOR_COLS <- c(.BEH_FLOW, .BEH_ENGAGE, .BEH_COLLECTIVE, .BEH_JACOBS)

.log("behaviour_wide: ", nrow(behaviour_wide), " x ", ncol(behaviour_wide),
     "  (", length(intersect(.BEH_FINAL_INDICATOR_COLS, names(behaviour_wide))),
     " indicators: ",
     length(intersect(.BEH_FLOW, names(behaviour_wide))), " flow + ",
     length(intersect(.BEH_ENGAGE, names(behaviour_wide))), " engagement + ",
     length(intersect(.BEH_COLLECTIVE, names(behaviour_wide))), " collective + ",
     length(intersect(.BEH_JACOBS, names(behaviour_wide))), " jacobs)")

# ---- Build cortisol and monoamines sheets -----------------------------------
.log("Reading endo CSV: ", .CSV_ENDO)
endo_raw <- readr::read_csv(.CSV_ENDO, show_col_types = FALSE)
.log("  rows=", nrow(endo_raw), " cols=", ncol(endo_raw))

# The endocrine tables key the experimental unit as "<tank>_<fish count>" in the
# `trial` column, e.g. "27_16". Rename it trial_key and derive phys_trial from
# it, so a reader can join a fish to its behavioural trial.
#
# WHY THIS MATTERS. Until 2026-08-28 the endocrine source carried an integer
# `trial` running 1-19. That column is NOT a trial identifier: it is the N-index
# of the sampling video (trial_raw in b3_cortisol_clean.csv). Only 6 of its 29
# distinct (tank, trial, treatment) combinations were even consistent with the
# behavioural design, and it assigned both treatments to the same tank x trial.
# Any deposit that invited a reader to join on it would have produced silent
# nonsense. (tank, fish count) is the one key every dataset shares, and it maps
# 1:1 onto all 16 behavioural trials with treatment agreeing in every case.
.trial_key_to_phys <- behaviour_wide %>%
  dplyr::distinct(tank, fish_density, phys_trial, treatment, motor_side) %>%
  dplyr::mutate(trial_key = paste0(tank, "_", fish_density))
stopifnot(nrow(.trial_key_to_phys) == 16L,
          !anyDuplicated(.trial_key_to_phys$trial_key))

.attach_design <- function(d, tag) {
  if (!"trial" %in% names(d)) stop("[workbook] ", tag, ": no trial column", call. = FALSE)
  d$trial_key <- as.character(d$trial)
  bad <- setdiff(unique(d$trial_key), .trial_key_to_phys$trial_key)
  if (length(bad))
    stop("[workbook] ", tag, ": trial_key value(s) with no behavioural trial: ",
         paste(bad, collapse = ", "),
         ". Expected \"<tank>_<fish count>\"; if this is the old integer coding, ",
         "regenerate easy_scripts_endo_dataset.csv.", call. = FALSE)
  out <- dplyr::left_join(
    d, dplyr::select(.trial_key_to_phys, trial_key, phys_trial,
                     .design_treatment = treatment, motor_side),
    by = "trial_key")
  # treatment arrives from the assay sheet as control/treat; the workbook has
  # always shown "exercise choice". Recode here. Doing it as a find-and-replace
  # over the saved xlsx is what produced the corrupted header "exercise
  # choicement" in the distributed file -- the replacement hit the substring
  # "treat" inside "treatment".
  out$treatment <- dplyr::recode(as.character(out$treatment),
                                 treat = "exercise choice", .default = "control")
  mism <- sum(out$treatment != as.character(out$.design_treatment))
  if (mism > 0)
    stop("[workbook] ", tag, ": ", mism,
         " row(s) whose treatment disagrees with the behavioural design",
         call. = FALSE)
  out$.design_treatment <- NULL
  .log("  ", tag, ": ", nrow(out), " rows over ",
       dplyr::n_distinct(out$phys_trial), " trials; treatment agrees with design")
  out
}

.endo_meta_head <- c("sample_id", "phys_trial", "trial_key", "tank", "fish_density",
                     "treatment", "motor_side", "sex")

# --- cortisol -----------------------------------------------------------------
# All 80 fish that were sampled, not just the 78 that yielded a cortisol value.
# B3_19 and B3_70 have no cortisol; B3_19 nonetheless has monoamines, so leaving
# it out would put a fish in the monoamine sheet that appears nowhere else.
.cort_all <- readr::read_csv(.CSV_CORT_CLEAN, show_col_types = FALSE) %>%
  dplyr::transmute(sample_id = as.character(sample_id),
                   tank      = as.integer(tank),
                   trial     = as.character(trial),
                   treatment = as.character(condition),
                   sex       = ifelse(is.na(sex) | sex == "NA", NA_character_, as.character(sex)),
                   plate     = as.character(plate),
                   plate_date = as.Date(plate_date),
                   cort_plasma = suppressWarnings(as.numeric(plasma_cortisol)),
                   weight_g   = suppressWarnings(as.numeric(wt)),
                   length_cm  = suppressWarnings(as.numeric(lenght)),
                   video_id   = as.character(video_id),
                   batch      = suppressWarnings(as.integer(batch)))
stopifnot(nrow(.cort_all) == 80L, !anyDuplicated(.cort_all$sample_id))

# Cross-check against the analysis frame: every assayed value must agree.
.cort_endo <- endo_raw %>%
  filter(analyte == "cort") %>%
  dplyr::transmute(sample_id = as.character(sample_id),
                   cort_endo = suppressWarnings(as.numeric(value)))
.chk <- dplyr::inner_join(.cort_all, .cort_endo, by = "sample_id")
if (nrow(.chk) != 78L)
  stop("[workbook] cortisol: ", nrow(.chk),
       " samples shared with the analysis frame, expected 78", call. = FALSE)
if (max(abs(.chk$cort_plasma - .chk$cort_endo), na.rm = TRUE) > 1e-9)
  stop("[workbook] cortisol: values disagree between b3_cortisol_clean.csv and ",
       "easy_scripts_endo_dataset.csv", call. = FALSE)
.log("  cortisol: 80 sampled fish, ", sum(!is.na(.cort_all$cort_plasma)),
     " with a value; cross-checked against the analysis frame")

cortisol <- .cort_all %>%
  .attach_design("cortisol") %>%
  mutate(fish_density   = as.integer(sub("^.*_", "", trial_key)),
         plate          = as.character(plate),   # keep text: holds "Evg"
         missing_reason = ifelse(is.na(cort_plasma), "no_value_returned", NA_character_)) %>%
  arrange(sample_id)

# Format plate_date as DD.MM.YYYY
if ("plate_date" %in% names(cortisol)) {
  pd <- tryCatch(as.Date(cortisol$plate_date), error = function(e) NULL)
  if (!is.null(pd) && !all(is.na(pd)))
    cortisol$plate_date <- format(pd, "%d.%m.%Y")
}
cortisol <- cortisol[, c(.endo_meta_head, "weight_g", "length_cm", "video_id",
                         "batch", "plate", "plate_date", "cort_plasma",
                         "missing_reason"), drop = FALSE]
rownames(cortisol) <- NULL
.log("cortisol: ", nrow(cortisol), " x ", ncol(cortisol))

# --- monoamines ---------------------------------------------------------------
# No plate column: these are HPLC-ED, not plate assays.
monoamines <- endo_raw %>%
  filter(analyte != "cort", !is.na(area), nzchar(as.character(area))) %>%
  mutate(cell_name = paste0(analyte, "__", area)) %>%
  select(any_of(c("sample_id", "tank", "trial", "treatment", "sex")),
         cell_name, value) %>%
  group_by(across(any_of(c("sample_id", "tank", "trial", "treatment", "sex"))),
           cell_name) %>%
  summarise(value = dplyr::first(value), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = cell_name, values_from = value) %>%
  .attach_design("monoamines") %>%
  mutate(fish_density = as.integer(sub("^.*_", "", trial_key)))

# Enforce canonical column order from INDICATOR_MAP_MONO. NE is no longer in the
# map, so any ne__ cell still present in the source is dropped here.
ordered_mono_cols <- INDICATOR_MAP_MONO$wb_col
missing_mono      <- setdiff(ordered_mono_cols, names(monoamines))
for (mc in missing_mono) monoamines[[mc]] <- NA_real_
monoamines <- monoamines[, c(.endo_meta_head, ordered_mono_cols), drop = FALSE]
monoamines <- monoamines[order(as.integer(sub("^B3_", "", monoamines$sample_id))), ,
                         drop = FALSE]
rownames(monoamines) <- NULL
.log("monoamines: ", nrow(monoamines), " x ", ncol(monoamines),
     "  (", length(ordered_mono_cols), " analyte x area cells)")

# ---- Build legend tibbles ---------------------------------------------------

# Helper: map a behaviour column to its display-group label
.beh_col_group <- function(col) {
  if (col %in% .BEH_META)       "metadata"
  else if (col %in% .BEH_FLOW)  "flow"
  else if (col %in% .BEH_ENGAGE) "engagement"
  else if (col %in% .BEH_COLLECTIVE) "collective"
  else if (col %in% .BEH_JACOBS) "jacobs"
  else "other"
}

build_legend <- function(sheet_name, df, kind,
                          rename_map = NULL) {
  # rename_map: named character vector, new_col → old_col, for description lookups
  if (nrow(df) == 0)
    return(data.frame(sheet = character(0), column = character(0),
                      group = character(0), type = character(0),
                      indicator_name = character(0), unit = character(0),
                      description = character(0), notes = character(0),
                      stringsAsFactors = FALSE))
  rows <- vector("list", ncol(df))
  for (i in seq_len(ncol(df))) {
    cnm      <- names(df)[i]
    lookup   <- if (!is.null(rename_map) && cnm %in% names(rename_map))
                  rename_map[[cnm]] else cnm
    type     <- NA_character_; unit <- NA_character_; ind_nm <- NA_character_
    desc     <- NA_character_; notes <- ""; grp <- NA_character_

    if (kind == "behaviour") {
      grp <- .beh_col_group(cnm)
      if (lookup %in% names(METADATA_DESCRIPTIONS)) {
        e <- METADATA_DESCRIPTIONS[[lookup]]
        type <- e$type; unit <- e$unit; desc <- e$description
      } else if (lookup %in% names(.IND_DESC$beh)) {
        e <- .IND_DESC$beh[[lookup]]
        type <- e$type; unit <- e$unit; ind_nm <- e$indicator_name
        desc <- e$description
        notes <- if (!is.null(e$notes) && nzchar(e$notes)) e$notes else ""
        # Annotate renames in notes
        if (!is.null(rename_map) && cnm %in% names(rename_map))
          notes <- trimws(paste(notes,
            sprintf("(renamed from '%s')", rename_map[[cnm]])))
      } else {
        e <- describe_src_column(lookup)
        type <- e$type; desc <- e$description; notes <- e$notes
      }

    } else if (kind %in% c("cortisol", "monoamines")) {
      if (lookup %in% names(METADATA_DESCRIPTIONS)) {
        e <- METADATA_DESCRIPTIONS[[lookup]]
        type <- e$type; unit <- e$unit; desc <- e$description
      } else if (lookup %in% names(.IND_DESC$endo)) {
        e <- .IND_DESC$endo[[lookup]]
        type <- e$type; unit <- e$unit; ind_nm <- e$indicator_name
        desc <- e$description; notes <- e$notes
      } else {
        type <- "unknown"; desc <- "Needs review"
      }

    } else {
      stop("[workbook] build_legend: unknown kind '", kind, "'", call. = FALSE)
    }

    rows[[i]] <- data.frame(
      sheet          = sheet_name,
      column         = cnm,
      group          = if (is.na(grp)) NA_character_ else grp,
      type           = if (is.na(type)) "unknown" else type,
      indicator_name = if (is.na(ind_nm)) NA_character_ else ind_nm,
      unit           = if (is.null(unit) || is.na(unit)) NA_character_ else unit,
      description    = if (is.na(desc)) "Needs review" else desc,
      notes          = if (is.null(notes)) "" else notes,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

# ---- Validation -------------------------------------------------------------
validate_workbook <- function() {
  errs <- character(0)
  E <- function(...) errs <<- c(errs, sprintf(...))

  bw <- behaviour_wide

  # ---- shape ----------------------------------------------------------------
  if (names(bw)[1] != "phys_trial")
    E("behaviour_wide: col 1 is '%s', expected 'phys_trial'", names(bw)[1])
  if (nrow(bw) != 48L)
    E("behaviour_wide: %d rows, expected 48 (16 trials x 3 intervals)", nrow(bw))
  if (dplyr::n_distinct(bw$phys_trial) != 16L)
    E("behaviour_wide: %d distinct phys_trial, expected 16",
      dplyr::n_distinct(bw$phys_trial))
  if (!setequal(sort(unique(bw$phys_trial)), 1:16))
    E("behaviour_wide: phys_trial is not 1..16")
  dup_b <- names(bw)[duplicated(tolower(names(bw)))]
  if (length(dup_b) > 0)
    E("behaviour_wide: duplicate column names: %s", paste(dup_b, collapse = ", "))

  # ---- interval level only --------------------------------------------------
  # This is a trial x interval table. A trial-aggregate column here would be one
  # value broadcast across three rows, which reads as three observations.
  agg <- grep("_agg$", names(bw), value = TRUE)
  if (length(agg) > 0)
    E("behaviour_wide: %d trial-aggregate column(s) present: %s",
      length(agg), paste(agg, collapse = ", "))
  if (dplyr::n_distinct(bw$interval) != 3L)
    E("behaviour_wide: %d distinct intervals, expected 3",
      dplyr::n_distinct(bw$interval))
  if (!setequal(sort(unique(bw$interval)), 1:3))
    E("behaviour_wide: interval is not 1..3")
  n_per <- as.integer(table(bw$phys_trial))
  if (!all(n_per == 3L))
    E("behaviour_wide: not every trial has 3 intervals (min %d, max %d)",
      min(n_per), max(n_per))

  # ---- reported outcomes present and populated ------------------------------
  miss_rep <- setdiff(.BEH_REPORTED, names(bw))
  if (length(miss_rep) > 0)
    E("behaviour_wide: %d REPORTED outcome(s) missing: %s",
      length(miss_rep), paste(miss_rep, collapse = ", "))
  for (cc in intersect(.BEH_REPORTED, names(bw)))
    if (all(is.na(bw[[cc]]))) E("behaviour_wide: reported outcome '%s' is entirely NA", cc)
  miss_inds <- setdiff(.BEH_FINAL_INDICATOR_COLS, names(bw))
  if (length(miss_inds) > 0)
    E("behaviour_wide: %d expected indicator col(s) MISSING: %s",
      length(miss_inds), paste(miss_inds, collapse = ", "))
  allna <- names(bw)[vapply(bw, function(x) all(is.na(x)), logical(1))]
  if (length(allna) > 0)
    E("behaviour_wide: %d all-NA column(s) -- a column with nothing behind it: %s",
      length(allna), paste(allna, collapse = ", "))

  # ---- the values are the post-area-fix ones --------------------------------
  # The distributed workbook shipped pre-correction occupancy for months
  # (alr_flow -1.526 where the live run says +1.390 at trial 1 / interval 1).
  # Anchor on the source frame, not on a literal, so this stays true if the
  # pipeline is legitimately re-run.
  .src <- beh_raw[beh_raw$phys_trial == 1 & beh_raw$interval == 1, ]
  .wbr <- bw[bw$phys_trial == 1 & bw$interval == 1, ]
  if (nrow(.src) == 1L && nrow(.wbr) == 1L) {
    for (pr in list(c("alr_flow", "logit_flow"), c("alr_high", "lr_high"),
                    c("crossings_per_session", "zone_flux_per_session"),
                    c("bouts_per_min", "bout_rate_flow"),
                    c("longest_flow_bout_s", "max_flow_bout_s"),
                    c("mean_nnd_cm", "mean_nnd_cm"))) {
      if (pr[1] %in% names(.wbr) && pr[2] %in% names(.src)) {
        d <- abs(as.numeric(.wbr[[pr[1]]]) - as.numeric(.src[[pr[2]]]))
        if (!is.finite(d) || d > 1e-12)
          E("behaviour_wide: %s does not match source %s (delta %s)",
            pr[1], pr[2], format(d))
      }
    }
  } else E("behaviour_wide: could not locate trial 1 / interval 1 for the value check")

  # ---- occupancy composition ------------------------------------------------
  acs <- rowSums(bw[, c("prop_high_ac", "prop_medium_ac",
                        "prop_low_ac", "prop_calm_ac"), drop = FALSE])
  if (max(abs(acs - 1)) > 1e-9)
    E("behaviour_wide: area-corrected sub-zone proportions do not sum to 1 (max dev %s)",
      format(max(abs(acs - 1))))

  # ---- dates ----------------------------------------------------------------
  if ("trial_date" %in% names(bw)) {
    bad_fmt <- !grepl("^[0-9]{2}[.][0-9]{2}[.][0-9]{4}$", bw$trial_date) & !is.na(bw$trial_date)
    if (any(bad_fmt)) E("behaviour_wide.trial_date: %d rows not in DD.MM.YYYY", sum(bad_fmt))
    n_uniq <- length(unique(bw$trial_date[!is.na(bw$trial_date)]))
    if (n_uniq != 8L)
      E("behaviour_wide.trial_date: %d distinct dates, expected 8", n_uniq)
  } else E("behaviour_wide: missing trial_date column")

  # ---- endocrine ------------------------------------------------------------
  if (!"cort_plasma" %in% names(cortisol)) E("cortisol: missing cort_plasma column")
  if (nrow(cortisol) != 80L)
    E("cortisol: %d rows, expected 80 (16 trials x 5 fish)", nrow(cortisol))
  if (sum(!is.na(cortisol$cort_plasma)) != 78L)
    E("cortisol: %d assayed values, expected 78", sum(!is.na(cortisol$cort_plasma)))
  if (any(is.na(cortisol$cort_plasma) & is.na(cortisol$missing_reason)))
    E("cortisol: empty value(s) with no missing_reason")
  .cn <- as.integer(table(cortisol$phys_trial))
  if (!all(.cn == 5L))
    E("cortisol: sampling is not 5 fish per trial (min %d, max %d)", min(.cn), max(.cn))
  if (!"phys_trial" %in% names(cortisol))  E("cortisol: missing phys_trial link column")
  if ("treatment" %in% names(cortisol) &&
      !setequal(unique(cortisol$treatment), c("control", "exercise choice")))
    E("cortisol.treatment levels are %s, expected control / exercise choice",
      paste(unique(cortisol$treatment), collapse = " / "))

  n_mono_ind <- sum(INDICATOR_MAP_MONO$wb_col %in% names(monoamines))
  if (n_mono_ind != nrow(INDICATOR_MAP_MONO))
    E("monoamines: %d of %d analyte x area columns present",
      n_mono_ind, nrow(INDICATOR_MAP_MONO))
  if (nrow(INDICATOR_MAP_MONO) != 24L)
    E("monoamines: %d analyte x area cells, expected 24 (NE excluded)",
      nrow(INDICATOR_MAP_MONO))
  if (nrow(monoamines) == 0L) E("monoamines: 0 rows")
  if (!"phys_trial" %in% names(monoamines)) E("monoamines: missing phys_trial link column")

  # Every physiology fish must resolve to a real behavioural trial.
  for (nm in c("cortisol", "monoamines")) {
    d <- get(nm)
    if ("phys_trial" %in% names(d)) {
      orphan <- setdiff(unique(d$phys_trial), unique(bw$phys_trial))
      if (length(orphan) > 0)
        E("%s: phys_trial value(s) with no behavioural trial: %s",
          nm, paste(orphan, collapse = ", "))
    }
  }

  orph_fish <- setdiff(monoamines$sample_id, cortisol$sample_id)
  if (length(orph_fish) > 0)
    E("monoamines: sample_id(s) with no row in the cortisol sheet: %s",
      paste(orph_fish, collapse = ", "))

  # ---- no observational-study or NE residue ---------------------------------
  ne_cols <- grep("^ne__", c(names(monoamines), names(cortisol)), value = TRUE)
  if (length(ne_cols) > 0) E("NE column(s) still present: %s", paste(ne_cols, collapse = ", "))
  obs_hits <- grep("^(BS|HY|TEL|OT)_|FISH_ID|trial_2nd_exp",
                   c(names(bw), names(cortisol), names(monoamines)), value = TRUE)
  if (length(obs_hits) > 0)
    E("observational-study column(s) still present: %s", paste(obs_hits, collapse = ", "))

  # ---- no column name that breaks a CSV round-trip --------------------------
  badnm <- grep("[ ()/]", c(names(bw), names(cortisol), names(monoamines)), value = TRUE)
  if (length(badnm) > 0)
    E("column name(s) with space, parenthesis or slash: %s", paste(badnm, collapse = ", "))

  if (length(errs) > 0)
    stop("Workbook validation FAILED:\n  - ",
         paste(errs, collapse = "\n  - "), call. = FALSE)
  .log("Validation: OK  (48 x ", ncol(behaviour_wide), " behaviour, ",
       nrow(cortisol), " cortisol, ", nrow(monoamines), " monoamines)")
  invisible(TRUE)
}

# ---- Build workbook (openxlsx) ----------------------------------------------
.log("Building openxlsx workbook ...")
wb <- openxlsx::createWorkbook()

# Styles
hdr_style <- openxlsx::createStyle(
  textDecoration = "bold", fgFill = "#D9D9D9", border = "Bottom",
  borderColour = "#000000", halign = "left", valign = "center")
legend_hdr_style <- openxlsx::createStyle(
  textDecoration = "bold", fgFill = "#C6E0B4", border = "Bottom",
  borderColour = "#000000", halign = "left", valign = "center")

add_data_sheet <- function(wb, sheet, df, header_style = hdr_style,
                           col_widths = "auto") {
  openxlsx::addWorksheet(wb, sheetName = sheet, gridLines = TRUE)
  openxlsx::writeData(wb, sheet = sheet, x = df, colNames = TRUE,
                      withFilter = TRUE, rowNames = FALSE)
  openxlsx::addStyle(wb, sheet = sheet, style = header_style,
                     rows = 1, cols = seq_len(ncol(df)), gridExpand = TRUE)
  openxlsx::freezePane(wb, sheet = sheet, firstRow = TRUE, firstCol = TRUE)
  openxlsx::setColWidths(wb, sheet = sheet, cols = seq_len(ncol(df)),
                         widths = col_widths)
}

# README
readme_header <- data.frame(
  Item = c(
    "Workbook",
    "Scope",
    "Generated (datetime)",
    "Source: behaviour CSV",
    "Source: bout metrics (BOUT engine)",
    "Source: state sequence (SEQ engine)",
    "Source: endocrine CSV",
    "Source: trial design log",
    "Indicator manifest",
    "Unit of observation",
    "Behaviour indicators",
    "Reported outcomes",
    "Cortisol indicators",
    "Monoamine indicators",
    "Noradrenaline",
    "Trial numbering",
    "Cortisol-plate raw assays",
    "Sheet conventions"
  ),
  Value = c(
    basename(.WB_OUT_PATH),
    "Preference experiment only (control vs exercise choice). The observational study was removed on 2026-08-29.",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    .CSV_BEH,
    .CSV_BOUT,
    .CSV_SEQ,
    .CSV_ENDO,
    .CSV_TRIALLOG,
    .MANIFEST_R,
    "Trial x interval. 16 trials x 3 within-trial intervals = 48 rows. No trial-aggregate column is included.",
    sprintf("%d indicator columns (%d zone preference + %d engagement + %d collective + %d Jacobs D)",
            length(.BEH_FINAL_INDICATOR_COLS), length(.BEH_FLOW),
            length(.BEH_ENGAGE), length(.BEH_COLLECTIVE), length(.BEH_JACOBS)),
    sprintf("%d of the behaviour columns are outcomes reported in the manuscript: %s",
            length(.BEH_REPORTED), paste(.BEH_REPORTED, collapse = ", ")),
    sprintf("%d column in cortisol (cort_plasma, ng/mL plasma)", nrow(INDICATOR_MAP_CORT)),
    sprintf("%d columns in monoamines (%d analytes x 4 brain areas)",
            nrow(INDICATOR_MAP_MONO),
            length(unique(sub("__.*", "", INDICATOR_MAP_MONO$wb_col)))),
    "Quantified but not deposited: NE is not entered into any reported model, so the analysed grid is 24 analyte x region cells, not 28.",
    "phys_trial (1-16) is the manuscript numbering and matches the figures. trial_seq (1-16) is chronological running order. They differ: trials 1<->2, 13<->14 and 15<->16 are swapped. Both are given.",
    "Out of scope. See: D:/CHOICE R SCRIPTS/choice exp for claude mono and cortisol/cortisol/",
    "Every data sheet has a paired legend sheet documenting each column."
  ),
  stringsAsFactors = FALSE
)

openxlsx::addWorksheet(wb, "00_README", gridLines = TRUE)
openxlsx::writeData(wb, "00_README", readme_header,
                    startRow = 1, colNames = TRUE, withFilter = FALSE, rowNames = FALSE)
openxlsx::addStyle(wb, "00_README", style = hdr_style,
                   rows = 1, cols = 1:2, gridExpand = TRUE)
openxlsx::freezePane(wb, "00_README", firstRow = TRUE, firstCol = TRUE)
openxlsx::setColWidths(wb, "00_README", cols = 1:2, widths = c(32, 90))

# ---- Per-group header colour styles -----------------------------------------
# Metadata:   #D9D9D9  (grey)
# Flow:       #BDD7EE  (light blue)
# Collective: #E2EFDA  (light green)
# Jacobs:     #FCE4D6  (light peach/orange)
.mk_hdr <- function(fill)
  openxlsx::createStyle(textDecoration = "bold", fgFill = fill,
                        border = "Bottom", borderColour = "#000000",
                        halign = "left", valign = "center")
.hdr_meta  <- .mk_hdr("#D9D9D9")
.hdr_flow  <- .mk_hdr("#BDD7EE")
.hdr_coll  <- .mk_hdr("#E2EFDA")
.hdr_jac   <- .mk_hdr("#FCE4D6")
.hdr_eng   <- .mk_hdr("#FFF2CC")   # engagement: light yellow

# Sheet names, kept as the names the distributed workbook already used so the
# deposit and any downstream reference to them do not move.
.SH_BEH  <- "Preference_exp_behavior"
.SH_CORT <- "preference_exp_cortisol"
.SH_MONO <- "preference_exp_monoamines"

# behaviour
add_data_sheet(wb, .SH_BEH, behaviour_wide,
               header_style = .hdr_meta)   # baseline grey; overwritten below per group

# Apply per-group header colours (overrides the baseline row 1 style)
.bw_names <- names(behaviour_wide)
.apply_col_style <- function(wb, sheet, cols_vec, style) {
  idx <- which(.bw_names %in% cols_vec)
  if (length(idx)) openxlsx::addStyle(wb, sheet, style,
                                      rows = 1, cols = idx, gridExpand = TRUE)
}
.apply_col_style(wb, .SH_BEH, .BEH_META,       .hdr_meta)
.apply_col_style(wb, .SH_BEH, .BEH_FLOW,       .hdr_flow)
.apply_col_style(wb, .SH_BEH, .BEH_ENGAGE,     .hdr_eng)
.apply_col_style(wb, .SH_BEH, .BEH_COLLECTIVE, .hdr_coll)
.apply_col_style(wb, .SH_BEH, .BEH_JACOBS,     .hdr_jac)

add_data_sheet(wb, paste0(.SH_BEH, "_legend"),
               build_legend(.SH_BEH, behaviour_wide, kind = "behaviour",
                            rename_map = .BEH_RENAME_MAP),
               header_style = legend_hdr_style)

# cortisol
add_data_sheet(wb, .SH_CORT, cortisol)
add_data_sheet(wb, paste0(.SH_CORT, "_legend"),
               build_legend(.SH_CORT, cortisol, kind = "cortisol"),
               header_style = legend_hdr_style)

# monoamines
add_data_sheet(wb, .SH_MONO, monoamines)
add_data_sheet(wb, "preferenceexp_monoamines_legend",
               build_legend(.SH_MONO, monoamines, kind = "monoamines"),
               header_style = legend_hdr_style)

# Sheet-set check: 3 data sheets, each with its legend, plus the README.
.expected_sheets <- c("00_README",
                      .SH_BEH,  paste0(.SH_BEH,  "_legend"),
                      .SH_CORT, paste0(.SH_CORT, "_legend"),
                      .SH_MONO, "preferenceexp_monoamines_legend")
.got_sheets <- names(wb)
if (!setequal(.got_sheets, .expected_sheets))
  stop("Sheet set is wrong.
  expected: ", paste(.expected_sheets, collapse = ", "),
       "
  got     : ", paste(.got_sheets, collapse = ", "), call. = FALSE)

# Validate then save
validate_workbook()
openxlsx::saveWorkbook(wb, .WB_OUT_PATH, overwrite = TRUE)
.log("Wrote: ", .WB_OUT_PATH,
     "   (", round(file.info(.WB_OUT_PATH)$size / 1024, 1), " KB)")
.log("Done. Sheets: ", paste(names(wb), collapse = ", "))
