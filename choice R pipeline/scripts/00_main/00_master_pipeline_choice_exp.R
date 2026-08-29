# =============================================================================
# MASTER PIPELINE — CHOICE EXPERIMENT (B3 COHORT)
# =============================================================================
# Runs the full analysis pipeline for the exercise-choice experiment.
# Two treatments: "exercise choice" vs "control".
# Fish density (4 / 8 / 12 / 16) is housing-tank density (RE candidate only).
#
# DESIGN: one observation per TRIAL × TIMEPOINT (3 × 20-min segments). The
# choice arena cannot maintain persistent fish identities, so all indicators
# are identity-free school-level summaries (NN-matched velocities, group
# dynamics, frame-level zone occupancy). Tank is the repeated-measures RE.
#
# STEPS:
#   1  Parse raw idtracker.ai data; assign multi-level zones (main_zone + sec_zone);
#      attach treatment, fish_density, tank, motor_side from trial metadata;
#      parse timepoint from session folder names (segmentN pattern).
#      → STEP1_output/STEP1_output_<ts>/master_fish_by_frame.csv
#
#   2  Build trial × timepoint activity summary (identity-free).
#      Per trial: bodylength threshold (mean of per-fish median body_length_cm)
#      → speed via frame-to-frame nearest-neighbour matching → prop_active
#      (fraction moving > 1 BL/s). Frame-level zone counts aggregate to a
#      LONG-format occupancy table (zone_level ∈ main/sec; zone ∈ flow/calm/
#      high/medium/low) plus a wide companion. Identity-less main-zone switch
#      count from NN-matched flow↔calm crossings, normalised to a 20-min
#      session (switches_per_session). zone_flux_per_session retained as a
#      diagnostic.
#      → STEP2_output/STEP2_output_<ts>/trial_activity_summary.csv
#      → STEP2_output/STEP2_output_<ts>/trial_occupancy_long.csv
#
#   2b Group spatial dynamics per trial × timepoint: NND, polarisation, IID,
#      hull area, centroid speed. Joined into trial_activity_summary.
#      → STEP2b_output/STEP2b_output_<ts>/group_dynamics_summary.csv
#
#   3  Zone-switch diagnostics: Pearson + Spearman correlation between
#      zone_flux_per_session (diagnostic) and switches_per_session (statistical
#      response), faceted by treatment, with scatter plot, PowerPoint slide,
#      and a one-paragraph results-text.
#      → zone_switch_diagnostics/zone_switch_diagnostics_<ts>/
#
#   4  Scatter graphs and zone occupancy plots → editable PowerPoint files.
#      Two output sets (timepoint: faceted by segment; aggregated: plain):
#        Standard indicator plots: x = treatment, colour + shape = treatment
#          (pch 17 = control, pch 16 = exercise choice).
#          Timepoint plots faceted with "Segment N  (0–20 min)" strip labels.
#          y-axis labels use Unicode superscripts (cm s⁻¹).
#        Main zone plot: x = zone (flow, calm), colour = zone, shape = treatment;
#          treatments dodged within each zone; y = raw seconds.
#        Secondary zone plot: x = zone (high, medium, low, calm), same design;
#          y = area-corrected seconds.
#      → STEP4_graphs/STEP4_graphs_<ts>/
#          timepoint/
#            standard_indicator_plots_tp.pptx
#            zone_plots_tp.pptx  (main zones + secondary zones; faceted)
#            png_preview/
#          aggregated/
#            standard_indicator_plots_agg.pptx
#            zone_plots_agg.pptx
#            png_preview/
#
#   5  LMM statistics on trial_activity_summary + trial_occupancy_long.
#      Module A (trial × timepoint): main-zone & sub-zone occupancy
#      (treatment × zone × timepoint), switches_per_session, prop_active,
#      group dynamics (NND, polarisation, IID, hull area, centroid speed).
#      Workflow: normality (SW) → variance (Levene) → transform if needed →
#      AICc RE selection → Type III ANOVA → Tukey post-hoc + CLD →
#      pseudoreplication check (trial-level vs tank-level) → BH FDR.
#      → STEP5_stats/STEP5_stats_<ts>/
#
#   5b LMM statistics — Module B (trial-level, timepoints collapsed). Same
#      indicator list as Module A with treatment as the only fixed effect.
#      → STEP5_stats_trial/STEP5_stats_trial_<ts>/
#
# HOW TO USE:
#   1. Set PIPELINE_DATA_DIR (the folder containing your session raw data).
#   2. Make sure trial_summary_choice_exp.xlsx, zone_reference_choice_exp_FT.csv,
#      and zone_reference_choice_exp_FD.csv are in the same folder as this script.
#   3. In each zone reference CSV replace the placeholder x_px/y_px values with
#      polygon corner coordinates measured from a video frame of the correct
#      orientation in ImageJ (FT file from an FT-orientation frame; FD file from
#      an FD-orientation frame).  One row per corner vertex, any number of corners.
#   4. Source this file: click Source in RStudio.
# =============================================================================


# =============================================================================
# ==== CONFIG — THE ONLY SECTION YOU NEED TO EDIT =============================
# =============================================================================

# Full path to the folder that contains your session raw data (forward slashes).
# Raw idtracker.ai sessions live outside this repository; see config.R.
# The literal formerly hardcoded here pointed at D:/GOT R SCRIPTS/, a
# different project root - a stale path that never matched this tree.
if (!exists("PIPELINE_DATA_DIR")) source(file.path(PROJECT_ROOT, "config.R"))

# Which step to start from (1 = run everything from scratch).
# 2   = Step 1 done; re-run from activity summary onward.
# 2.5 = Steps 1-2 done; re-run group dynamics (Step 2b) + graphs + stats.
# 4   = Steps 1-2b done; re-run graphs only.
# 5   = Steps 1-4 done; re-run summary-level stats only.
START_FROM_STEP <- 2.5  # 2026-08-07 zone-geometry fix: re-run group dynamics + graphs + stats
                        # on the corrected STEP1/STEP2 outputs (Steps 1-2 already regenerated)

# Stop the pipeline immediately if any step throws an error (recommended).
STOP_ON_ERROR <- TRUE

# Frame height in pixels — informational; zone orientation is handled natively
# via the separate FT/FD polygon CSV files (no y-mirroring is applied).
# Change only if your video resolution differs from 1080p.
FRAME_HEIGHT_PX <- 1080L

# If TRUE, each side's zone-reference polygons are linearly rescaled (per axis,
# independently in x and y) to fit the actual observed x_interp/y_interp extent
# for that motor_side before point-in-polygon testing. Use this when a
# zone_reference CSV was digitised on a differently-sized frame/canvas than the
# actual tracking output (e.g. zone_reference_choice_exp_FT.csv measured on a
# ~3400x1768 canvas while FT tracking data only spans ~1810x1123 px).
RESCALE_ZONE_REF_TO_DATA <- TRUE

# Ordered treatment and density levels — controls factor ordering in plots/models.
TREATMENT_LEVELS  <- c("control", "exercise choice")
DENSITY_LEVELS    <- c(4L, 8L, 12L, 16L)

# Timepoint levels (1 = segment 1, 2 = segment 2, 3 = segment 3).
# Each segment is 20 minutes.  Timepoints are parsed from session folder names
# containing "segmentN" (e.g. "mydata_segment2").
TIMEPOINT_LEVELS  <- 1:3

# Duration of each 20-minute session segment.  Used by STEP2 to normalise
# switches_per_session.  Change only if your recording protocol differs.
SESSION_MINUTES_PER_TIMEPOINT <- 20

# Number of parallel workers for per-session processing in STEP1 and STEP2.
# NULL (default) = detectCores() - 1, using all but one logical core.
# Set to 1 to disable parallelisation entirely (useful for debugging).
# Hard ceiling of 8 applied automatically to avoid thrashing on shared machines.
N_PARALLEL_WORKERS <- NULL

# When TRUE: session folders already listed in STEP1_output/processed_sessions.txt
# are skipped.  Only new session_* folders are processed and merged with the
# existing master_fish_by_frame.csv from the previous run.
# When FALSE (default): all session_* folders in PIPELINE_DATA_DIR are processed.
SKIP_PROCESSED_SESSIONS <- FALSE

# Set to TRUE to treat fish_density as an ordered factor in models.
# Set to FALSE to treat it as a numeric covariate.
DENSITY_AS_FACTOR <- TRUE

# Zone areas in spatial units (1 unit = 10 × 25 cm section).
# Used to compute area-corrected seconds (time / zone_area) for secondary-zone
# models, and to rescale those to percentages summing to 100 for graphs.
# Values come from the flow-distribution map (Supplementary material, Figure 3).
# Change these values if you update the experimental arena layout.
ZONE_AREA_UNITS <- c(
  high   = 10,   # high-intensity flow zone
  medium = 18,   # medium-intensity flow zone
  low    = 13,   # low-intensity flow zone
  calm   = 42    # calm zone
)

# ---------------------------------------------------------------------------
# JUMP DETECTION — applied in Step 1 after trajectory read, before zone
# assignment.  A "jump" is a frame-to-frame displacement that exceeds a
# velocity threshold, indicating a tracking error (identity swap, lost track).
# The offending frame's coordinates are set to NA and then linearly
# interpolated from the nearest valid positions on either side.
#
# METHOD "sd_multiple"  (default — mirrors idtracker.ai validator):
#   threshold = mean(speed) + JUMP_SD_MULT × sd(speed)
#   Speeds are computed globally across all fish and all trials in the session.
#   The default multiplier of 10 matches the idtracker.ai validator GUI default.
#
# METHOD "percentile"  (mirrors idtracker.ai postprocessing):
#   threshold = JUMP_PCT_MULT × percentile(speed, JUMP_PCT)
#   Default: 2 × 99th percentile of all frame-to-frame speeds.
#
# Speed units: pixels per frame, divided by the actual frame gap so that
# tracks with missing frames are not penalised (gap-spanning displacements
# are normalised by the number of frames they span).
#
# Set JUMP_DETECTION_ENABLED <- FALSE to skip the filter entirely.
# ---------------------------------------------------------------------------
JUMP_DETECTION_ENABLED  <- TRUE
JUMP_THRESHOLD_METHOD   <- "sd_multiple"  # "sd_multiple" | "percentile"
JUMP_SD_MULT            <- 10.0           # sd_multiple: SDs above mean
JUMP_PCT                <- 99             # percentile: which percentile
JUMP_PCT_MULT           <- 2.0            # percentile: multiplier
JUMP_INTERP             <- TRUE           # TRUE = linear interpolation after removal
#                                         # FALSE = leave removed frames as NA

# ---------------------------------------------------------------------------
# IDENTITY SWITCH DETECTION — applied in Step 1 after jump detection, before
# zone assignment.  An identity switch occurs when the tracker assigns two fish
# the wrong labels — typically after a crossing event — so neither fish shows
# an individual velocity spike but the collective assignment is suboptimal.
#
# TIER 1 — Mutual nearest-neighbour (no extra packages):
#   For each consecutive frame pair, for each pair of fish (i, j):
#     ratio = (d_self_i + d_self_j) / (d_cross_i + d_cross_j)
#   where d_self  = distance from own previous position,
#         d_cross = distance from the other fish's previous position.
#   ratio > SWITCH_TIER1_RATIO means the swapped assignment is shorter
#   → flag as a switch candidate.  Default 1.5 requires the swap to be at
#   least 50% more parsimonious than the current assignment.
#
# TIER 2 — Optimal assignment via Hungarian algorithm (clue::solve_LSAP):
#   Builds the N×N squared-distance cost matrix across all fish and solves
#   the minimum-cost assignment.  Runs only on Tier 1 candidate frames to
#   confirm events; automatically disabled if clue is not installed.
#
# CONSOLIDATION:
#   Consecutive flagged frames for the same pair that are separated by at most
#   SWITCH_CONSOL_GAP frames are merged into one event.  Events shorter than
#   SWITCH_MIN_DURATION frames are discarded as noise.
#
# CORRECTION:
#   Events are processed chronologically per trial.  For each confirmed event
#   (fish_A, fish_B, frame_start), position columns are swapped between the
#   two fish for all frames >= frame_start.  Applying a second swap of the
#   same pair later automatically un-does the first one.
#
# PROXIMITY PRE-FILTER (optional):
#   Set SWITCH_PROXIMITY_PX > 0 to only check pairs whose inter-fish distance
#   in the previous frame was within that many pixels.  This can speed up
#   processing on long sessions since switches only occur when fish are close.
#   0 = disabled (check all pairs at all frames).
#
# Set SWITCH_DETECTION_ENABLED <- FALSE to skip entirely.
# ---------------------------------------------------------------------------
SWITCH_DETECTION_ENABLED <- TRUE
SWITCH_TIER1_RATIO       <- 1.5    # min improvement ratio to flag (Tier 1)
SWITCH_MIN_DURATION      <- 2L     # min consecutive flagged frames per event
SWITCH_CONSOL_GAP        <- 5L     # max frame gap to bridge into one event
SWITCH_USE_HUNGARIAN     <- TRUE   # confirm via Hungarian (requires clue pkg)
SWITCH_PROXIMITY_PX      <- 0      # proximity pre-filter in px (0 = off)


# =============================================================================
# ==== INTERNAL SETUP — do not edit below this line ===========================
# =============================================================================

# Resolve the directory that contains THIS script at source() time.
.pipeline_dir <- tryCatch({
  frames     <- sys.frames()
  ofile_envs <- Filter(function(f) exists("ofile", envir = f, inherits = FALSE), frames)
  if (length(ofile_envs) > 0) {
    normalizePath(
      dirname(get("ofile", envir = ofile_envs[[length(ofile_envs)]])),
      winslash = "/", mustWork = FALSE
    )
  } else {
    getwd()
  }
}, error = function(e) getwd())

# Validate PIPELINE_DATA_DIR
if (!nzchar(trimws(PIPELINE_DATA_DIR)) ||
    PIPELINE_DATA_DIR == "C:/path/to/your/choice_exp/session/folder") {
  stop(
    "PIPELINE_DATA_DIR has not been set.\n",
    "Set the CHOICE_EXP_DATA_DIR environment variable to the folder of checked ",
    "idtracker.ai sessions, or edit config.R. Only Step 1 reads it.",
    call. = FALSE
  )
}
PIPELINE_DATA_DIR <- normalizePath(PIPELINE_DATA_DIR, winslash = "/", mustWork = FALSE)
# Only Step 1 reads PIPELINE_DATA_DIR (raw idtracker.ai session folders), so
# only validate it when Step 1 will actually run. START_FROM_STEP > 1 means
# STEP1_output already exists and this path is never touched.
if (START_FROM_STEP <= 1 && !dir.exists(PIPELINE_DATA_DIR)) {
  stop("PIPELINE_DATA_DIR does not exist on disk:\n  ", PIPELINE_DATA_DIR, call. = FALSE)
}

if (!is.numeric(START_FROM_STEP) ||
    !START_FROM_STEP %in% c(1, 2, 2.5, 4, 4.5, 5, 5.5)) {
  stop("START_FROM_STEP must be one of: 1, 2, 2.5, 4, 4.5, 5, 5.5", call. = FALSE)
}

# Warn about missing prerequisite objects when skipping early steps
if (START_FROM_STEP >= 3 &&
    !exists("master_fish_by_frame", envir = .GlobalEnv, inherits = FALSE)) {
  warning(
    "START_FROM_STEP = ", START_FROM_STEP,
    " but 'master_fish_by_frame' is not in the R environment.\n",
    "Step 2 output will be searched in STEP2_output/ and STEP1_output/ automatically.",
    call. = FALSE
  )
}
if (START_FROM_STEP >= 4 &&
    !exists("fish_activity_summary", envir = .GlobalEnv, inherits = FALSE)) {
  warning(
    "START_FROM_STEP = ", START_FROM_STEP,
    " but 'fish_activity_summary' is not in the R environment.\n",
    "Steps 4-5 will search STEP2_output/ automatically.",
    call. = FALSE
  )
}


# =============================================================================
# ==== LOAD TRIAL METADATA ====================================================
# =============================================================================

suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(dplyr))

# Expect the Excel file in the same folder as this script; fall back to
# choice R pipeline/data/ (its actual location as of 2026-08-07) if absent.
.trial_excel_path <- file.path(.pipeline_dir, "trial_summary_choice_exp.xlsx")
if (!file.exists(.trial_excel_path)) {
  .trial_excel_path_fallback <- file.path(dirname(dirname(.pipeline_dir)), "data",
                                           "trial_summary_choice_exp.xlsx")
  if (file.exists(.trial_excel_path_fallback)) .trial_excel_path <- .trial_excel_path_fallback
}

if (!file.exists(.trial_excel_path)) {
  stop(
    "Trial metadata not found: ", .trial_excel_path,
    "\nMake sure trial_summary_choice_exp.xlsx is in the same folder as this script, ",
    "or in choice R pipeline/data/.",
    call. = FALSE
  )
}

cat("Loading trial metadata from:", .trial_excel_path, "\n")
TRIAL_META_RAW <- readxl::read_excel(.trial_excel_path)

# Rename 'condition' -> 'treatment' if present; tolerate existing 'treatment' column
if ("condition" %in% names(TRIAL_META_RAW) && !("treatment" %in% names(TRIAL_META_RAW))) {
  TRIAL_META_RAW <- dplyr::rename(TRIAL_META_RAW, treatment = condition)
  cat("  Renamed column 'condition' -> 'treatment'\n")
}

# Validate required columns
# Expected schema: video_ID, trial, tank, motor_side, fish_density, treatment, date
.required_meta_cols <- c("trial", "treatment", "fish_density", "tank", "motor_side")
.missing_meta       <- setdiff(.required_meta_cols, names(TRIAL_META_RAW))
if (length(.missing_meta) > 0) {
  stop(
    "trial_summary_choice_exp.xlsx is missing required column(s): ",
    paste(.missing_meta, collapse = ", "),
    "\nExpected columns: video_ID, trial, tank, motor_side, fish_density, treatment, date.",
    call. = FALSE
  )
}

# Normalise treatment strings; apply factor levels
TRIAL_META_RAW$treatment <- trimws(tolower(as.character(TRIAL_META_RAW$treatment)))

TRIAL_META <- TRIAL_META_RAW %>%
  dplyr::mutate(
    trial        = suppressWarnings(as.integer(trial)),
    fish_density = suppressWarnings(as.integer(fish_density)),
    tank         = as.character(tank),
    motor_side   = trimws(as.character(motor_side)),
    treatment    = factor(treatment, levels = TREATMENT_LEVELS)
  )

# Keep video_ID as character if present
if ("video_ID" %in% names(TRIAL_META))
  TRIAL_META$video_ID <- as.character(TRIAL_META$video_ID)

# Rename 'date' -> 'trial_date' (the column is called 'date' in the new schema)
if (!"trial_date" %in% names(TRIAL_META)) {
  date_candidates <- c("date", "Date", "trial_date", "DATE")
  date_col <- intersect(date_candidates, names(TRIAL_META))[1]
  if (!is.na(date_col)) {
    TRIAL_META <- dplyr::rename(TRIAL_META, trial_date = !!date_col)
    cat("  Renamed column '", date_col, "' -> 'trial_date'\n", sep = "")
  }
}

cat("  TRIAL_META loaded:", nrow(TRIAL_META), "rows x", ncol(TRIAL_META), "columns\n")
cat("  Treatments found :", paste(unique(as.character(TRIAL_META$treatment)), collapse = ", "), "\n")
cat("  Densities found  :", paste(sort(unique(TRIAL_META$fish_density)), collapse = ", "), "\n")
cat("  Tanks found      :", paste(sort(unique(TRIAL_META$tank)), collapse = ", "), "\n")
cat("  Motor sides found:", paste(unique(TRIAL_META$motor_side), collapse = ", "), "\n\n")

# Publish to GlobalEnv so step scripts can access it without re-loading
assign("TRIAL_META",   TRIAL_META,   envir = .GlobalEnv)
assign("TRIAL_META_RAW", TRIAL_META_RAW, envir = .GlobalEnv)


# =============================================================================
# ==== LOAD ZONE REFERENCES (FT and FD) =======================================
# =============================================================================
# Two separate polygon-vertex CSV files define zone boundaries:
#   zone_reference_choice_exp_FT.csv  — measured from an FT-orientation frame
#   zone_reference_choice_exp_FD.csv  — measured from an FD-orientation frame
#
# Accepts idtracker.ai ROI-export format (one row per vertex):
#   ROI_index, ROI_name, Vertex_index, X, Y
# or the legacy script format:
#   zone_label, vertex, x_px, y_px  (with optional pre-filled main_zone/sec_zone)
#
# Zone names are normalised to lowercase.  Discontinuous low-flow sub-zones
# exported as low1/low2/low3 (or any lowN) are treated as separate polygons
# but all assigned sec_zone = "low".  The broad "flow" polygon (if present)
# acts as a fallback for flow-area points not covered by any sub-zone polygon.

.load_zone_ref <- function(filename, label) {
  path <- file.path(.pipeline_dir, filename)
  if (!file.exists(path)) {
    stop("Zone reference CSV not found: ", path,
         "\nMake sure ", filename, " is in the same folder as this script.",
         call. = FALSE)
  }
  cat("Loading zone reference (", label, ") from: ", path, "\n", sep = "")
  zr <- read.csv(path, stringsAsFactors = FALSE)

  # --- Accept idtracker.ai column names OR legacy column names ---------------
  .rename_if <- function(df, candidates, target) {
    if (target %in% names(df)) return(df)
    hit <- intersect(candidates, names(df))[1L]
    if (!is.na(hit)) names(df)[names(df) == hit] <- target
    df
  }
  zr <- .rename_if(zr, c("ROI_name",     "roi_name"),     "zone_label")
  zr <- .rename_if(zr, c("Vertex_index", "vertex_index"), "vertex")
  zr <- .rename_if(zr, c("X"),                             "x_px")
  zr <- .rename_if(zr, c("Y"),                             "y_px")

  .req  <- c("zone_label", "vertex", "x_px", "y_px")
  .miss <- setdiff(.req, names(zr))
  if (length(.miss) > 0)
    stop(filename, " missing column(s): ", paste(.miss, collapse = ", "),
         ". Columns present: ", paste(names(zr), collapse = ", "), call. = FALSE)

  # --- Normalise zone names: lowercase; keep lowN labels separate -----------
  zr$zone_label <- tolower(trimws(as.character(zr$zone_label)))

  # Derive sec_zone: low/low1/low2/... all → "low"; flow polygon → NA
  zr$sec_zone <- dplyr::case_when(
    grepl("^low\\d*$", zr$zone_label) ~ "low",
    zr$zone_label == "calm"           ~ "calm",
    zr$zone_label == "high"           ~ "high",
    zr$zone_label == "medium"         ~ "medium",
    TRUE                              ~ NA_character_
  )

  # Derive main_zone: sub-zones (high/medium/lowN) are inside the flow area
  zr$main_zone <- dplyr::case_when(
    zr$zone_label == "calm"                          ~ "calm",
    zr$zone_label == "flow"                          ~ "flow",
    grepl("^(high|medium|low\\d*)$", zr$zone_label) ~ "flow",
    TRUE                                             ~ NA_character_
  )

  # Warn on unrecognised zone names
  .known  <- c("flow", "calm", "high", "medium",
               grep("^low", unique(zr$zone_label), value = TRUE))
  .unknown <- setdiff(unique(zr$zone_label), .known)
  if (length(.unknown) > 0)
    warning(filename, ": unrecognised zone label(s) — main_zone/sec_zone will be NA: ",
            paste(.unknown, collapse = ", "), call. = FALSE)

  # Warn on polygons with fewer than 3 vertices
  .vcount <- tapply(zr$vertex, zr$zone_label, length)
  .bad    <- names(.vcount)[.vcount < 3L]
  if (length(.bad) > 0)
    warning(filename, ": zone(s) with fewer than 3 vertices — check coordinates: ",
            paste(.bad, collapse = ", "), call. = FALSE)

  cat("  Zones (", label, "): ",
      paste(unique(zr$zone_label), collapse = ", "), "\n", sep = "")
  zr
}

ZONE_REF_FT <- .load_zone_ref("zone_reference_choice_exp_FT.csv", "FT")
ZONE_REF_FD <- .load_zone_ref("zone_reference_choice_exp_FD.csv", "FD")
cat("\n")

assign("ZONE_REF_FT", ZONE_REF_FT, envir = .GlobalEnv)
assign("ZONE_REF_FD", ZONE_REF_FD, envir = .GlobalEnv)


# =============================================================================
# ==== ZONE ASSIGNMENT FUNCTION ===============================================
# =============================================================================
# This function is called in Step 1 to attach main_zone + sec_zone to every
# row of master_fish_by_frame.  It is defined here (in the master pipeline) so
# that all step scripts automatically inherit it from GlobalEnv.

## Linearly rescale a zone-reference polygon set so its extent matches the
## observed data's extent. This corrects zone_reference CSVs that were
## digitised on a differently sized frame/canvas than the actual tracking
## output. A degenerate axis (zero-width data or polygon range) is left
## unscaled to avoid Inf/NaN.
##
## The data range is taken from trimmed quantiles (trim_prob on each tail),
## NOT raw min/max: single-frame tracking glitches (a handful of frames with
## wildly out-of-tank x/y, e.g. x > 5000 px or y < 0 px when the real arena
## is ~0-1800 px) can otherwise dominate the observed range and badly distort
## an already-correctly-calibrated reference. Polygon vertices themselves are
## still used at their full extent since they are precise digitised points,
## not noisy measurements.
##
## isotropic = TRUE applies a SINGLE scale factor (the geometric mean of the
## independently-fit x/y factors) to both axes instead of stretching x and y
## by different amounts. Tested against this dataset and found to INCREASE
## unassigned (NA) points overall (FD 3.5%->4.5%, FT 1.6%->4.5%) because the
## true required correction genuinely differs per axis; forcing one scale
## just fits both axes worse. Default is FALSE (independent per-axis fit).
## Kept as an option in case a future zone_reference file's mismatch really
## is a uniform canvas resize rather than a per-axis distortion.
.rescale_zone_ref_to_data <- function(zr, data_x, data_y, label = "",
                                       trim_prob = 0.001, isotropic = FALSE) {
  data_x <- data_x[is.finite(data_x)]
  data_y <- data_y[is.finite(data_y)]
  if (length(data_x) == 0L || length(data_y) == 0L) return(zr)

  poly_x_rng <- range(zr$x_px, na.rm = TRUE)
  poly_y_rng <- range(zr$y_px, na.rm = TRUE)
  data_x_rng <- stats::quantile(data_x, c(trim_prob, 1 - trim_prob),
                                 na.rm = TRUE, names = FALSE)
  data_y_rng <- stats::quantile(data_y, c(trim_prob, 1 - trim_prob),
                                 na.rm = TRUE, names = FALSE)

  poly_x_span <- diff(poly_x_rng)
  poly_y_span <- diff(poly_y_rng)
  scale_x <- if (poly_x_span > 0) diff(data_x_rng) / poly_x_span else 1
  scale_y <- if (poly_y_span > 0) diff(data_y_rng) / poly_y_span else 1

  if (isotropic) {
    scale_iso   <- sqrt(scale_x * scale_y)
    scale_x_use <- scale_iso
    scale_y_use <- scale_iso
  } else {
    scale_x_use <- scale_x
    scale_y_use <- scale_y
  }

  zr$x_px <- data_x_rng[1] + (zr$x_px - poly_x_rng[1]) * scale_x_use
  zr$y_px <- data_y_rng[1] + (zr$y_px - poly_y_rng[1]) * scale_y_use

  cat("  Rescaled zone reference (", label, "): x-scale=", round(scale_x_use, 4),
      " y-scale=", round(scale_y_use, 4),
      if (isotropic) paste0(" [isotropic; independent fit was x=", round(scale_x,4),
                             " y=", round(scale_y,4), "]") else "",
      " [data range trimmed at ", trim_prob * 100, "% / ",
      (1 - trim_prob) * 100, "%]\n", sep = "")
  zr
}

## Point-in-polygon test using the ray casting algorithm.
## Vectorised over points; iterates over polygon edges.
.pip_raycasting <- function(px, py, vx, vy) {
  n      <- length(vx)
  inside <- integer(length(px))
  j      <- n
  for (i in seq_len(n)) {
    yi <- vy[i]; yj <- vy[j]
    xi <- vx[i]; xj <- vx[j]
    cross <- ((yi > py) != (yj > py)) &
             (px < (xj - xi) * (py - yi) / (yj - yi) + xi)
    inside <- bitwXor(inside, as.integer(cross))
    j <- i
  }
  inside > 0L
}

assign_zones_choice_exp <- function(
    df,
    zone_ref_ft = get("ZONE_REF_FT", envir = .GlobalEnv),
    zone_ref_fd = get("ZONE_REF_FD", envir = .GlobalEnv),
    trial_meta  = get("TRIAL_META",  envir = .GlobalEnv)
) {
  # Returns df with new columns: zone_label, main_zone, sec_zone.
  # zone_ref_ft is applied to rows where motor_side == "FT";
  # zone_ref_fd is applied to rows where motor_side == "FD".
  # Both references are long-format polygon-vertex data.frames:
  #   zone_label, main_zone, sec_zone, vertex, x_px, y_px
  # Zone assignment uses point-in-polygon ray casting; first-match wins.

  requireNamespace("data.table", quietly = TRUE)
  dt <- data.table::as.data.table(df)

  # Attach motor_side per trial_id if not already present
  dt[, trial_id := suppressWarnings(as.integer(trial_id))]
  if (!"motor_side" %in% names(dt)) {
    .tm <- as.data.frame(trial_meta)
    # Accept both 'trial_id' and 'trial' as the key column
    .tid_col <- if ("trial_id" %in% names(.tm)) "trial_id" else "trial"
    meta_dt <- data.table::as.data.table(data.frame(
      trial_id   = as.integer(.tm[[.tid_col]]),
      motor_side = as.character(.tm[["motor_side"]])
    ))
    dt <- merge(dt, meta_dt, by = "trial_id", all.x = TRUE)
  }

  dt[, c("zone_label", "main_zone", "sec_zone") := list(
    NA_character_, NA_character_, NA_character_
  )]

  # Build a named list of polygon definitions from a zone reference data.frame.
  # Sub-zone polygons are tested first so they take priority over the broad
  # "flow" polygon that spatially contains them.
  .build_polys <- function(zr) {
    zr  <- as.data.frame(zr)
    .named_pri <- c(high = 1L, medium = 2L, calm = 4L, flow = 9L)
    zr$.pri <- ifelse(grepl("^low\\d*$", zr$zone_label), 3L,
                      .named_pri[zr$zone_label])
    zr$.pri[is.na(zr$.pri)] <- 5L
    zr  <- zr[order(zr$.pri, zr$zone_label, as.integer(zr$vertex)), ]
    zls <- unique(as.character(zr$zone_label))
    stats::setNames(
      lapply(zls, function(zl) {
        sub <- zr[zr$zone_label == zl, ]
        list(x_px      = as.numeric(sub$x_px),
             y_px      = as.numeric(sub$y_px),
             main_zone = as.character(sub$main_zone[1L]),
             sec_zone  = as.character(sub$sec_zone[1L]))
      }),
      zls
    )
  }

  .do_rescale <- exists("RESCALE_ZONE_REF_TO_DATA", envir = .GlobalEnv) &&
                  isTRUE(get("RESCALE_ZONE_REF_TO_DATA", envir = .GlobalEnv))

  for (side in c("FT", "FD")) {
    zr  <- as.data.frame(if (side == "FT") zone_ref_ft else zone_ref_fd)
    idx <- which(!is.na(dt$motor_side) & dt$motor_side == side)
    if (length(idx) == 0L) next

    px <- as.numeric(dt$x_interp[idx])
    py <- as.numeric(dt$y_interp[idx])

    if (.do_rescale) {
      zr <- .rescale_zone_ref_to_data(zr, px, py, label = side)
    }
    polys <- .build_polys(zr)

    assigned <- integer(length(idx))   # 0 = unassigned; positive = polygon index
    for (zi in seq_along(polys)) {
      inside <- .pip_raycasting(px, py, polys[[zi]]$x_px, polys[[zi]]$y_px)
      take   <- inside & assigned == 0L
      assigned[take] <- zi
    }

    assigned_na <- ifelse(assigned == 0L, NA_integer_, as.integer(assigned))
    zl_vec <- names(polys)[assigned_na]
    # Collapse discontinuous low sub-zone polygons (low1/low2/low3/...) into a
    # single "low" zone_label, mirroring the sec_zone grouping already done in
    # .load_zone_ref(). The separate polygons are still tested individually
    # above so their (possibly disjoint) shapes are preserved.
    zl_vec <- ifelse(grepl("^low\\d*$", zl_vec), "low", zl_vec)
    mz_vec <- vapply(assigned_na,
                     function(a) if (is.na(a)) NA_character_ else polys[[a]]$main_zone,
                     character(1))
    sz_vec <- vapply(assigned_na,
                     function(a) if (is.na(a)) NA_character_ else polys[[a]]$sec_zone,
                     character(1))

    data.table::set(dt, idx, "zone_label", zl_vec)
    data.table::set(dt, idx, "main_zone",  mz_vec)
    data.table::set(dt, idx, "sec_zone",   sz_vec)
  }

  as.data.frame(dt)
}

assign("assign_zones_choice_exp", assign_zones_choice_exp, envir = .GlobalEnv)


# =============================================================================
# ==== PUBLISH CONFIG TO GLOBALENV ============================================
# =============================================================================

assign("PIPELINE_DATA_DIR",          PIPELINE_DATA_DIR,          envir = .GlobalEnv)
assign("FRAME_HEIGHT_PX",            FRAME_HEIGHT_PX,            envir = .GlobalEnv)
assign("TREATMENT_LEVELS",           TREATMENT_LEVELS,           envir = .GlobalEnv)
assign("DENSITY_LEVELS",             DENSITY_LEVELS,             envir = .GlobalEnv)
assign("DENSITY_AS_FACTOR",          DENSITY_AS_FACTOR,          envir = .GlobalEnv)
assign("TIMEPOINT_LEVELS",           TIMEPOINT_LEVELS,           envir = .GlobalEnv)
assign("ZONE_AREA_UNITS",            ZONE_AREA_UNITS,            envir = .GlobalEnv)
assign("SKIP_PROCESSED_SESSIONS",   SKIP_PROCESSED_SESSIONS,    envir = .GlobalEnv)
assign(".pipeline_dir_choice",       .pipeline_dir,              envir = .GlobalEnv)
assign("JUMP_DETECTION_ENABLED",    JUMP_DETECTION_ENABLED,     envir = .GlobalEnv)
assign("JUMP_THRESHOLD_METHOD",     JUMP_THRESHOLD_METHOD,      envir = .GlobalEnv)
assign("JUMP_SD_MULT",              JUMP_SD_MULT,               envir = .GlobalEnv)
assign("JUMP_PCT",                  JUMP_PCT,                   envir = .GlobalEnv)
assign("JUMP_PCT_MULT",             JUMP_PCT_MULT,              envir = .GlobalEnv)
assign("JUMP_INTERP",               JUMP_INTERP,                envir = .GlobalEnv)
assign("SWITCH_DETECTION_ENABLED",  SWITCH_DETECTION_ENABLED,   envir = .GlobalEnv)
assign("SWITCH_TIER1_RATIO",        SWITCH_TIER1_RATIO,         envir = .GlobalEnv)
assign("SWITCH_MIN_DURATION",       SWITCH_MIN_DURATION,        envir = .GlobalEnv)
assign("SWITCH_CONSOL_GAP",         SWITCH_CONSOL_GAP,          envir = .GlobalEnv)
assign("SWITCH_USE_HUNGARIAN",      SWITCH_USE_HUNGARIAN,       envir = .GlobalEnv)
assign("SWITCH_PROXIMITY_PX",            SWITCH_PROXIMITY_PX,            envir = .GlobalEnv)
assign("SESSION_MINUTES_PER_TIMEPOINT",  SESSION_MINUTES_PER_TIMEPOINT,  envir = .GlobalEnv)
assign("N_PARALLEL_WORKERS",             N_PARALLEL_WORKERS,             envir = .GlobalEnv)


# =============================================================================
# ==== STEP RUNNER ============================================================
# =============================================================================

.pipeline_source <- function(step_number, script_name, start_from, stop_on_error) {
  if (step_number < start_from) {
    cat(strrep("-", 70), "\n")
    cat("SKIPPED (START_FROM_STEP=", start_from, "): STEP ", step_number,
        " — ", script_name, "\n", sep = "")
    cat(strrep("-", 70), "\n\n")
    return(invisible(NULL))
  }

  # Look for the script in the same folder as this master pipeline
  script_path <- file.path(.pipeline_dir, script_name)
  if (!file.exists(script_path)) {
    stop("Script not found: ", script_name,
         "\n  Searched: ", script_path, call. = FALSE)
  }

  cat(strrep("=", 70), "\n")
  cat("STEP ", step_number, " — ", script_name, "\n", sep = "")
  cat(strrep("=", 70), "\n\n")

  t_start <- proc.time()[["elapsed"]]

  result <- tryCatch(
    source(script_path, local = FALSE, echo = FALSE, verbose = FALSE),
    error = function(e) e
  )

  elapsed <- round(proc.time()[["elapsed"]] - t_start, 1)

  if (inherits(result, "error")) {
    cat("\n", strrep("!", 70), "\n", sep = "")
    cat("ERROR in STEP ", step_number, " — ", script_name, "\n", sep = "")
    cat("Message: ", conditionMessage(result), "\n", sep = "")
    cat("Elapsed: ", elapsed, " s\n", sep = "")
    cat(strrep("!", 70), "\n\n")
    if (isTRUE(stop_on_error)) {
      stop("Pipeline halted at step ", step_number, ".\n",
           "Fix the error and re-run with START_FROM_STEP = ", step_number, ".",
           call. = FALSE)
    } else {
      warning("Step ", step_number, " failed but STOP_ON_ERROR = FALSE; continuing.",
              call. = FALSE)
    }
  } else {
    cat("\n", strrep("-", 70), "\n", sep = "")
    cat("DONE: STEP ", step_number, " — ", script_name,
        "  [", elapsed, " s]\n\n", sep = "")
  }
  invisible(NULL)
}


# =============================================================================
# ==== PIPELINE SUMMARY =======================================================
# =============================================================================

.pipeline_start_time <- proc.time()[["elapsed"]]

cat(strrep("=", 70), "\n")
cat("MASTER PIPELINE — CHOICE EXPERIMENT\n")
cat("Data dir:        ", PIPELINE_DATA_DIR, "\n")
cat("Starting at step:", START_FROM_STEP, "\n")
cat("Stop on error:   ", STOP_ON_ERROR, "\n")
cat("Frame height px: ", FRAME_HEIGHT_PX, "\n")
cat("Density as factor:", DENSITY_AS_FACTOR, "\n")
cat("Timepoint levels :", paste(TIMEPOINT_LEVELS, collapse = ", "), "\n")
cat("Zone area units  : high=", ZONE_AREA_UNITS[["high"]],
    " medium=", ZONE_AREA_UNITS[["medium"]],
    " low=",    ZONE_AREA_UNITS[["low"]],
    " calm=",   ZONE_AREA_UNITS[["calm"]], "\n", sep = "")
cat(strrep("=", 70), "\n\n")


# =============================================================================
# ==== RUN THE STEPS ==========================================================
# =============================================================================

# Step 1: Parse raw idtracker data; assign zones (main + secondary); attach metadata.
#          Reads:    raw session CSVs (PIPELINE_DATA_DIR); TRIAL_META; ZONE_REF.
#          Produces: STEP1_output/STEP1_output_<ts>/master_fish_by_frame.csv
.pipeline_source(1, "idtracker_STEP1_choice_exp.R",              START_FROM_STEP, STOP_ON_ERROR)

# ---------------------------------------------------------------------------
# PARALLEL SETUP — initialise the future/furrr backend before STEP2 and STEP2b.
# Both scripts call furrr::future_map() over individual trials; the plan set
# here propagates to those calls automatically.
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  if (!requireNamespace("parallel", quietly = TRUE))
    install.packages("parallel")
  if (!requireNamespace("future",   quietly = TRUE))
    install.packages("future")
  if (!requireNamespace("furrr",    quietly = TRUE))
    install.packages("furrr")
  library(future)
  library(furrr)
})

.max_workers <- min(8L, parallel::detectCores(logical = TRUE) - 1L)
.n_workers   <- if (is.null(N_PARALLEL_WORKERS) || !is.numeric(N_PARALLEL_WORKERS)) {
  .max_workers
} else {
  as.integer(max(1L, min(as.integer(N_PARALLEL_WORKERS), .max_workers)))
}

if (.n_workers > 1L) {
  future::plan(future::multisession, workers = .n_workers)
  cat(sprintf(
    "Parallel backend: multisession  |  workers = %d  |  logical cores = %d\n\n",
    .n_workers, parallel::detectCores(logical = TRUE)
  ))
} else {
  future::plan(future::sequential)
  cat("Parallel backend: sequential (N_PARALLEL_WORKERS = 1 or single-core machine)\n\n")
}

assign("PIPELINE_WORKERS", .n_workers, envir = .GlobalEnv)

# Step 2: Trial-level activity summary (identity-free: prop_active, switches,
#          zone occupancy long + wide).
#          Reads:    master_fish_by_frame (GlobalEnv → STEP1_output/<latest>/).
#          Produces: STEP2_output/STEP2_output_<ts>/
#                      trial_activity_summary.csv   (one row per trial × timepoint)
#                      trial_occupancy_long.csv     (long: zone_level × zone)
.pipeline_source(2, "activity_analysis_STEP2_choice_exp.R",      START_FROM_STEP, STOP_ON_ERROR)

# Free RAM before STEP 2.5's heavy parallel frame computation.
# STEP 2.5 reloads master_fish_by_frame from STEP1_output if not in GlobalEnv.
# trial_activity_summary / trial_occupancy_long are already on disk.
# Only evict master_fish_by_frame (7M-row object). STEP2.5 reloads it from disk.
# fish_activity_summary / trial_* are tiny — kept in GlobalEnv for STEP2.5.
for (.obj in c("master_fish_by_frame")) {
  if (exists(.obj, envir = .GlobalEnv, inherits = FALSE))
    rm(list = .obj, envir = .GlobalEnv)
}
invisible(gc(verbose = FALSE, full = TRUE))
cat("[Master] RAM freed before STEP 2.5.\n")

# Step 3: DEPRECATED — frame-level LMM analysis removed.
#          Fish-level LMM analysis (the biological unit of interest) is in Step 5.
# .pipeline_source(3, "activity_analysis_STEP3_choice_exp.R",    START_FROM_STEP, STOP_ON_ERROR)

# Step 2b: Group dynamics — Tier 1 collective-behaviour metrics.
#           Uses swaRm (CRAN) for NND, polarisation, convex-hull area.
#           Group-level:      mean IID, NND, polarisation, hull area, centroid speed.
#           Individual-level: per-fish NND, centroid distance, turning rate.
#           Reads:    master_fish_by_frame + fish_activity_summary (GlobalEnv).
#           Produces: STEP2b_output/STEP2b_output_<ts>/
#                       group_dynamics_summary.csv  (1 row per trial × timepoint)
#                       fish_social_context.csv     (1 row per fish per timepoint)
#                       fish_activity_summary_with_social.csv
#           Updates:  fish_activity_summary in GlobalEnv (adds social columns).
#           Requires: install.packages("swaRm")   — falls back to base R if absent.
.pipeline_source(2.5, "group_dynamics_STEP2b_choice_exp.R",         START_FROM_STEP, STOP_ON_ERROR)

# Free master_fish_by_frame after STEP 2.5 — downstream steps load from disk.
for (.obj in c("master_fish_by_frame", "group_dynamics_summary")) {
  if (exists(.obj, envir = .GlobalEnv, inherits = FALSE))
    rm(list = .obj, envir = .GlobalEnv)
}
invisible(gc(verbose = FALSE, full = TRUE))
cat("[Master] RAM freed after STEP 2.5.\n")

# Step 4: Standard indicator scatter plots + zone occupancy combined plots.
#          Two sets: timepoint/ (faceted by segment) and aggregated/ (per-fish mean).
#          Standard: x = treatment, colour + shape = treatment (pch 17 / 16).
#          Zones:    x = zone, colour = zone, shape = treatment; dodged within zone.
#          y-axis labels use Unicode superscripts (e.g. cm s⁻¹).
#          Reads:    fish_activity_summary (GlobalEnv → STEP2_output/<latest>/).
#          Produces: STEP4_graphs/STEP4_graphs_<ts>/
#                      timepoint/
#                        standard_indicator_plots_tp.pptx
#                        zone_plots_tp.pptx  (faceted slides)
#                        png_preview/
#                      aggregated/
#                        standard_indicator_plots_agg.pptx
#                        zone_plots_agg.pptx
#                        png_preview/
.pipeline_source(4, "activity_analysis_GRAPHS_choice_exp.R",     START_FROM_STEP, STOP_ON_ERROR)

# Step 4b: Zone-switch diagnostics — correlation between zone_flux_per_session
#           (diagnostic) and switches_per_session (statistical response).
#           Reads:    trial_activity_summary (GlobalEnv → STEP2_output/<latest>/).
#           Produces: zone_switch_diagnostics/zone_switch_diagnostics_<ts>/
.pipeline_source(4.5, "zone_switch_diagnostics_choice_exp.R",    START_FROM_STEP, STOP_ON_ERROR)

# Step 5: LMM Module A — trial × timepoint statistics.
#          Indicator list: main-zone occupancy, sub-zone occupancy,
#          switches_per_session, prop_active, NND, polarisation, IID, hull
#          area, centroid speed. Normality (Shapiro-Wilk) → variance (Levene)
#          → transform if needed → AICc RE selection → Type III Wald χ² ANOVA
#          → Tukey post-hoc + CLD → pseudoreplication check → BH FDR.
#          Reads:    trial_activity_summary, trial_occupancy_long, group_dynamics_summary.
#          Produces: STEP5_stats/STEP5_stats_<ts>/
.pipeline_source(5, "activity_analysis_STATS_choice_exp.R",      START_FROM_STEP, STOP_ON_ERROR)

# Step 5b: LMM Module B — trial-level statistics (timepoints collapsed).
#           Treatment-only fixed effect; (1|tank) RE when tanks have multiple
#           trials, otherwise plain LM. Same indicator list as Module A.
#           Requires Module A helpers in GlobalEnv (run_lmm_analysis, fmt_p, etc.).
#           Produces: STEP5_stats_trial/STEP5_stats_trial_<ts>/
.pipeline_source(5.5, "activity_analysis_STATS_trial_choice_exp.R", START_FROM_STEP, STOP_ON_ERROR)


# =============================================================================
# ==== FINISHED ===============================================================
# =============================================================================

.total_elapsed <- round(proc.time()[["elapsed"]] - .pipeline_start_time, 1)

cat(strrep("=", 70), "\n")
cat("MASTER PIPELINE — CHOICE EXPERIMENT — COMPLETE\n")
cat("Total elapsed:   ", .total_elapsed, " s\n", sep = "")
cat("Output folders created inside: ", getwd(), "\n")
cat(strrep("=", 70), "\n")
