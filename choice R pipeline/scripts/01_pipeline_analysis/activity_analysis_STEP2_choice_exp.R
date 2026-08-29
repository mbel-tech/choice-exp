# =============================================================================
# STEP 2 — CHOICE EXPERIMENT (school-level, identity-free)
#
# WHY THIS REFACTOR
#   In the choice experiment, idtracker.ai cannot guarantee persistent fish
#   identities within a session (positions are correct, but identities get
#   reassigned). To keep all treatments comparable, the unit of analysis is
#   trial × timepoint and every indicator is identity-free:
#
#     - Zone occupancy aggregated frame-by-frame (no IDs needed).
#     - School-level activity (> 1 BL/s) recovered via frame-to-frame
#       nearest-neighbor matching of detections.
#     - Main-zone switch counts recovered via the same NN matching
#       (an "identity-less" switch each time a matched blob crosses
#       the flow↔calm boundary between consecutive frames).
#     - zone_flux_per_session retained as a diagnostic (correlates with
#       switches_per_session — see STEP3 zone_switch_diagnostics).
#
#   STEP1 trajectory preprocessing (gap interpolation, jump removal,
#   identity-switch correction) is retained as-is; this script only changes
#   how cleaned positions are aggregated.
#
# INPUT (from GlobalEnv or most recent STEP1_output/ run):
#   master_fish_by_frame  with columns:
#     trial_id, fish_id, frame, time, x_interp, y_interp,
#     main_zone, sec_zone, treatment, fish_density, tank, motor_side,
#     trial_date, timepoint, body_length_cm (per fish_id), length_unit (px/cm)
#
# OUTPUTS (GlobalEnv + disk):
#   trial_activity_summary       — one row per (trial_id, timepoint), wide
#   trial_occupancy_long         — long table for joint main+sub zone stats
#   master_fish_by_frame         — passed through (for STEP2b / heatmaps)
#   STEP2_output/STEP2_output_<timestamp>/
#     trial_activity_summary.csv
#     trial_occupancy_long.csv
#     diagnostics/
# =============================================================================


# =============================================================================
# ==== 0) PACKAGES + HELPERS ==================================================
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(readxl)
  library(data.table)
  library(stringr)
  library(ggplot2)
  library(purrr)
  library(furrr)
  library(future)
  library(jsonlite)
})

ts_msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), " | ", ..., "\n", sep = "")

.get_global <- function(name, default = NULL) {
  if (exists(name, envir = .GlobalEnv, inherits = FALSE))
    get(name, envir = .GlobalEnv, inherits = FALSE)
  else default
}

.step2_dir <- tryCatch({
  frames     <- sys.frames()
  ofile_envs <- Filter(function(f) exists("ofile", envir = f, inherits = FALSE), frames)
  if (length(ofile_envs) > 0)
    normalizePath(dirname(get("ofile", envir = ofile_envs[[length(ofile_envs)]])),
                  winslash = "/", mustWork = FALSE)
  else getwd()
}, error = function(e) getwd())

.find_latest_csv <- function(step_name, csv_filename) {
  parent <- file.path(getwd(), step_name)
  if (!dir.exists(parent)) return(NULL)
  subdirs <- list.dirs(parent, full.names = TRUE, recursive = FALSE)
  subdirs <- subdirs[grepl(paste0("^", step_name, "_\\d{8}_\\d{6}$"), basename(subdirs))]
  if (length(subdirs) == 0) return(NULL)
  latest <- subdirs[which.max(file.mtime(subdirs))]
  candidate <- file.path(latest, csv_filename)
  if (file.exists(candidate)) candidate else NULL
}

.make_output_dir <- function(step_name) {
  parent <- file.path(getwd(), step_name)
  if (!dir.exists(parent)) dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  run_dir <- file.path(parent, paste0(step_name, "_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  run_dir
}


# =============================================================================
# ==== 1) LOAD TRIAL_META + master_fish_by_frame ==============================
# =============================================================================

if (!exists("TRIAL_META", envir = .GlobalEnv, inherits = FALSE)) {
  ts_msg("TRIAL_META not in GlobalEnv; loading from Excel...")
  .xl <- file.path(.step2_dir, "trial_summary_choice_exp.xlsx")
  if (!file.exists(.xl)) stop("trial_summary_choice_exp.xlsx not found at: ", .xl)
  .tm <- readxl::read_excel(.xl)
  if ("condition" %in% names(.tm) && !("treatment" %in% names(.tm)))
    .tm <- dplyr::rename(.tm, treatment = condition)
  TRIAL_META <- .tm %>% dplyr::mutate(
    trial_id     = suppressWarnings(as.integer(trial_id)),
    fish_density = suppressWarnings(as.integer(fish_density)),
    tank         = as.character(tank),
    motor_side   = trimws(as.character(motor_side)),
    treatment    = trimws(tolower(as.character(treatment)))
  )
  assign("TRIAL_META", TRIAL_META, envir = .GlobalEnv)
}

if (!exists("master_fish_by_frame", envir = .GlobalEnv, inherits = FALSE)) {
  ts_msg("master_fish_by_frame not in GlobalEnv; searching disk...")
  .csv_path <- .find_latest_csv("STEP1_output", "master_fish_by_frame.csv")
  if (is.null(.csv_path)) {
    .data_dir <- .get_global("PIPELINE_DATA_DIR")
    if (!is.null(.data_dir)) {
      for (.sd in c(file.path(.data_dir, "master_outputs", "by_frame"),
                    file.path(.data_dir, "master_outputs"))) {
        .cand <- file.path(.sd, "master_fish_by_frame.csv")
        if (file.exists(.cand)) { .csv_path <- .cand; break }
      }
    }
  }
  if (is.null(.csv_path))
    stop("Cannot find master_fish_by_frame.csv. Run Step 1 first.", call. = FALSE)
  master_fish_by_frame <- readr::read_csv(.csv_path, show_col_types = FALSE)
  assign("master_fish_by_frame", master_fish_by_frame, envir = .GlobalEnv)
  ts_msg("Loaded from: ", .csv_path, " (", nrow(master_fish_by_frame), " rows)")
}

df <- get("master_fish_by_frame", envir = .GlobalEnv, inherits = FALSE)

if (!all(c("main_zone", "sec_zone") %in% names(df))) {
  if (exists("assign_zones_choice_exp", envir = .GlobalEnv, inherits = FALSE)) {
    ts_msg("main_zone/sec_zone missing; running assign_zones_choice_exp()...")
    df <- assign_zones_choice_exp(df)
  } else {
    warning("main_zone/sec_zone columns not found. Zone metrics will be NA.", call. = FALSE)
    df$main_zone <- NA_character_
    df$sec_zone  <- NA_character_
  }
}

.required <- c("trial_id", "frame", "time", "x_interp", "y_interp")
.miss <- setdiff(.required, names(df))
if (length(.miss) > 0)
  stop("master_fish_by_frame missing columns: ", paste(.miss, collapse = ", "), call. = FALSE)

if (!"timepoint" %in% names(df)) {
  warning("'timepoint' missing; setting NA.", call. = FALSE)
  df$timepoint <- NA_integer_
}
df$timepoint <- as.integer(df$timepoint)

ts_msg("Working data: ", nrow(df), " rows, ",
       dplyr::n_distinct(df$trial_id), " trials, ",
       dplyr::n_distinct(paste(df$trial_id, df$timepoint)),
       " (trial × timepoint) cells")


# =============================================================================
# ==== 2) OUTPUT DIRECTORY ====================================================
# =============================================================================

STEP2_OUT <- .make_output_dir("STEP2_output")
dir.create(file.path(STEP2_OUT, "diagnostics"), showWarnings = FALSE)
ts_msg("Output directory: ", STEP2_OUT)


# =============================================================================
# ==== 3) PER-TRIAL BODY LENGTH (cm) ==========================================
# =============================================================================
# Primary: read median_body_length (px) from idtracker.ai session.json files.
#   Session folders under PIPELINE_DATA_DIR are matched to trial_ids by the
#   same alphabetical sort that STEP1 uses when assigning trial_ids.
#   bl_cm_trial = median_body_length_px / length_unit_trial.
# Fallback: if session.json is absent or median_body_length unavailable,
#   compute per-fish median body_length_cm from master_fish_by_frame.
# =============================================================================

# --- length_unit per trial (px/cm) ------------------------------------------
.have_lu      <- "length_unit" %in% names(df)
.lu_per_trial <- if (.have_lu) {
  .lu_raw <- df %>%
    dplyr::filter(is.finite(length_unit), length_unit > 0) %>%
    dplyr::group_by(trial_id) %>%
    dplyr::summarise(length_unit_trial = stats::median(length_unit, na.rm = TRUE),
                     .groups = "drop")
  .lu_fallback <- mean(.lu_raw$length_unit_trial, na.rm = TRUE)
  .all_trials  <- data.frame(trial_id = sort(unique(df$trial_id)))
  .lu_out <- dplyr::left_join(.all_trials, .lu_raw, by = "trial_id") %>%
    dplyr::mutate(length_unit_trial = dplyr::if_else(
      is.finite(length_unit_trial), length_unit_trial, .lu_fallback))
  .n_imputed <- sum(!.lu_raw$trial_id %in% .all_trials$trial_id[
    .all_trials$trial_id %in% .lu_raw$trial_id])
  if (anyNA(.lu_out$length_unit_trial))
    warning("length_unit missing for some trials and no valid sessions to average from.",
            call. = FALSE)
  else {
    .n_imp <- nrow(.all_trials) - nrow(.lu_raw)
    if (.n_imp > 0)
      ts_msg("length_unit missing for ", .n_imp, " trial(s); imputed with cross-session mean (",
             round(.lu_fallback, 3), " px/cm)")
  }
  .lu_out
} else {
  data.frame(trial_id = sort(unique(df$trial_id)),
             length_unit_trial = NA_real_)
}

# --- Primary: session.json median_body_length --------------------------------
.pipeline_data_dir <- .get_global("PIPELINE_DATA_DIR")
.session_bl_px     <- NULL

if (!is.null(.pipeline_data_dir) && dir.exists(.pipeline_data_dir)) {
  .all_subdirs <- sort(list.dirs(.pipeline_data_dir, recursive = FALSE,
                                  full.names = TRUE))
  .sess_dirs   <- .all_subdirs[grepl("N\\d+_segment\\d+",
                                      basename(.all_subdirs), perl = TRUE)]
  if (length(.sess_dirs) > 0) {
    ts_msg("Found ", length(.sess_dirs),
           " session folder(s); reading session.json for median_body_length...")
    .rows <- lapply(seq_along(.sess_dirs), function(i) {
      .jf <- file.path(.sess_dirs[i], "session.json")
      if (!file.exists(.jf)) return(NULL)
      tryCatch({
        .j  <- jsonlite::read_json(.jf, simplifyVector = TRUE)
        bl  <- .j[["median_body_length"]]
        if (!is.null(bl) && is.numeric(bl) && is.finite(bl) && bl > 0)
          data.frame(trial_id = i, session_bl_px = as.numeric(bl))
        else NULL
      }, error = function(e) {
        warning("session.json parse error: ", basename(.sess_dirs[i]),
                " — ", conditionMessage(e), call. = FALSE)
        NULL
      })
    })
    .rows <- Filter(Negate(is.null), .rows)
    if (length(.rows) > 0) .session_bl_px <- dplyr::bind_rows(.rows)
  }
}

# --- Compute bl_cm_trial -----------------------------------------------------
.have_bl_cm <- "body_length_cm" %in% names(df)
.have_bl_px <- "body_length_px" %in% names(df)

if (!.have_bl_cm && .have_bl_px && .have_lu) {
  ts_msg("body_length_cm missing — recomputing from body_length_px / length_unit")
  df$body_length_cm <- df$body_length_px / df$length_unit
  .have_bl_cm <- TRUE
}

if (!is.null(.session_bl_px) && nrow(.session_bl_px) > 0) {
  # Primary path: session.json
  bl_per_trial <- .session_bl_px %>%
    dplyr::left_join(.lu_per_trial, by = "trial_id") %>%
    dplyr::mutate(
      bl_cm_trial   = dplyr::if_else(
        is.finite(length_unit_trial) & length_unit_trial > 0,
        session_bl_px / length_unit_trial, NA_real_),
      n_fish_for_bl = NA_integer_
    ) %>%
    dplyr::select(trial_id, bl_cm_trial, n_fish_for_bl)
  ts_msg("Body length (cm) from session.json: ",
         round(min(bl_per_trial$bl_cm_trial, na.rm = TRUE), 2), " \u2013 ",
         round(max(bl_per_trial$bl_cm_trial, na.rm = TRUE), 2),
         " (", sum(!is.na(bl_per_trial$bl_cm_trial)), " sessions matched)")
} else if (.have_bl_cm) {
  # Fallback: per-fish median from master_fish_by_frame
  ts_msg("session.json not available \u2014 computing body length from per-fish median.")
  bl_per_trial <- df %>%
    dplyr::filter(is.finite(body_length_cm), body_length_cm > 0) %>%
    dplyr::group_by(trial_id, fish_id) %>%
    dplyr::summarise(bl_fish_cm = stats::median(body_length_cm, na.rm = TRUE),
                     .groups = "drop") %>%
    dplyr::group_by(trial_id) %>%
    dplyr::summarise(bl_cm_trial   = mean(bl_fish_cm, na.rm = TRUE),
                     n_fish_for_bl = dplyr::n(),
                     .groups = "drop")
  ts_msg("Per-trial bodylength (cm, fallback): ",
         round(min(bl_per_trial$bl_cm_trial, na.rm = TRUE), 2), " \u2013 ",
         round(max(bl_per_trial$bl_cm_trial, na.rm = TRUE), 2),
         " (n trials = ", nrow(bl_per_trial), ")")
} else {
  warning("No body_length data available \u2014 bl_cm_trial = NA.", call. = FALSE)
  bl_per_trial <- data.frame(
    trial_id      = sort(unique(df$trial_id)),
    bl_cm_trial   = NA_real_,
    n_fish_for_bl = 0L
  )
}


# =============================================================================
# ==== 4) IDENTITY-FREE FRAME-TO-FRAME NN MATCHING ============================
# =============================================================================
# For each (trial_id, frame) collect all valid (x_interp, y_interp) detections.
# Match each detection to the nearest detection in the next frame (greedy NN,
# capped by max plausible displacement = 2 * bl_cm_trial * dt to reject
# crossings). Speed = displacement / dt, in cm/s after applying length_unit.
#
# Outputs:
#   - active_frame: per-frame activity counters (n_active, n_total) + main_zone
#                   transitions (n_main_switches in this transition).
#   - aggregated to (trial_id, timepoint) below in section 5.
# =============================================================================

ts_msg("Running identity-free NN matching for speed + main-zone switches...")

.nn_match_trial <- function(dt_trial, bl_cm, length_unit_trial) {
  # dt_trial: data.table for one trial, columns frame, time, x_interp, y_interp,
  #           main_zone, timepoint
  # bl_cm: scalar body length in cm (NA OK, threshold cap then disabled)
  # length_unit_trial: scalar pixels-per-cm calibration for this trial

  data.table::setkey(dt_trial, frame)
  frames <- sort(unique(dt_trial$frame))
  if (length(frames) < 2L) {
    return(data.table::data.table(
      frame_from = integer(), frame_to = integer(),
      timepoint = integer(),
      n_total_from = integer(), n_matched = integer(),
      n_active = integer(), n_main_switches = integer(),
      dt_s = numeric()
    ))
  }

  # Plausible displacement cap (cm). If bl_cm NA, use a generous fallback.
  max_disp_cm <- if (is.finite(bl_cm) && bl_cm > 0) 2 * bl_cm else Inf
  has_lu <- is.finite(length_unit_trial) && length_unit_trial > 0
  px_to_cm <- function(px) if (has_lu) px / length_unit_trial else NA_real_

  # Split into per-frame row-index lists ONCE (2026-08-07 perf fix): the
  # previous version re-ran dt_trial[frame == f] (a binary-search subset) on
  # every iteration of this loop -- ~30,000 times per trial -- which measured
  # at ~128s/trial (>100 min across 48 trials) and was the reason this step
  # never completed after the zone-geometry fix required a full pipeline
  # re-run. Precomputing the split is a pure lookup-strategy change: verified
  # to produce IDENTICAL per-frame subsets to the original, ~9x faster
  # (~14s/trial).
  idx_by_frame <- split(seq_len(nrow(dt_trial)), dt_trial$frame)
  idx_by_frame <- idx_by_frame[as.character(frames)]

  out_list <- vector("list", length(frames) - 1L)
  for (k in seq_len(length(frames) - 1L)) {
    f1 <- frames[k]; f2 <- frames[k + 1L]
    a <- dt_trial[idx_by_frame[[k]]]
    b <- dt_trial[idx_by_frame[[k + 1L]]]
    a <- a[is.finite(x_interp) & is.finite(y_interp)]
    b <- b[is.finite(x_interp) & is.finite(y_interp)]
    if (nrow(a) == 0L || nrow(b) == 0L) next

    dt_s <- as.numeric(b$time[1] - a$time[1])
    if (!is.finite(dt_s) || dt_s <= 0) next
    max_disp_px <- if (is.finite(max_disp_cm) && has_lu)
      max_disp_cm * length_unit_trial * dt_s + 1e-6 else Inf

    # Greedy NN: for each blob in `a`, pick nearest blob in `b` (without
    # replacement) within max_disp_px.
    pa <- as.matrix(a[, .(x_interp, y_interp)])
    pb <- as.matrix(b[, .(x_interp, y_interp)])
    used_b <- rep(FALSE, nrow(pb))
    speeds_cm_s <- rep(NA_real_, nrow(pa))
    zone_pairs_a <- a$main_zone
    zone_pairs_b <- character(nrow(pa))

    for (i in seq_len(nrow(pa))) {
      avail <- which(!used_b)
      if (length(avail) == 0L) break
      d2 <- (pb[avail, 1L] - pa[i, 1L])^2 + (pb[avail, 2L] - pa[i, 2L])^2
      jrel <- which.min(d2)
      d_px <- sqrt(d2[jrel])
      if (d_px > max_disp_px) next
      j <- avail[jrel]
      used_b[j] <- TRUE
      speeds_cm_s[i] <- px_to_cm(d_px) / dt_s
      zone_pairs_b[i] <- as.character(b$main_zone[j])
    }

    # Activity threshold: > 1 BL/s
    is_active <- if (is.finite(bl_cm) && bl_cm > 0)
      is.finite(speeds_cm_s) & speeds_cm_s > bl_cm else rep(NA, length(speeds_cm_s))

    # Main-zone switches: matched-pair zones differ between {flow,calm}
    valid_zone <- zone_pairs_a %in% c("flow", "calm") &
                  zone_pairs_b %in% c("flow", "calm")
    n_switches <- sum(valid_zone & zone_pairs_a != zone_pairs_b, na.rm = TRUE)

    out_list[[k]] <- data.table::data.table(
      frame_from = f1, frame_to = f2,
      timepoint  = a$timepoint[1],
      n_total_from = nrow(a),
      n_matched    = sum(!is.na(speeds_cm_s)),
      n_active     = sum(is_active, na.rm = TRUE),
      n_main_switches = n_switches,
      dt_s = dt_s
    )
  }
  data.table::rbindlist(out_list)
}

dt_all <- data.table::as.data.table(df[, c("trial_id", "fish_id", "frame", "time",
                                           "x_interp", "y_interp", "main_zone",
                                           "timepoint")])
data.table::setkey(dt_all, trial_id, frame)

bl_lookup <- setNames(bl_per_trial$bl_cm_trial,       bl_per_trial$trial_id)
lu_lookup <- setNames(.lu_per_trial$length_unit_trial, .lu_per_trial$trial_id)

# Pre-split by trial so each worker receives only its own slice.
# This avoids serialising the full dt_all to every parallel worker.
.dt_by_trial <- split(dt_all, by = "trial_id", keep.by = TRUE)

.n_workers_step2 <- .get_global("PIPELINE_WORKERS", 1L)
ts_msg("NN matching: ", length(.dt_by_trial),
       " trial(s) | ", .n_workers_step2, " worker(s)")

# Ensure a future plan is active; fall back to sequential if none was set by
# the master pipeline (e.g. when sourcing STEP2 standalone).
if (inherits(future::plan(), "SequentialFuture") && .n_workers_step2 > 1L) {
  future::plan(future::multisession, workers = .n_workers_step2)
}

nn_results <- furrr::future_map(
  .dt_by_trial,
  function(.dt_slice) {
    .tid_str <- as.character(.dt_slice$trial_id[1L])
    .tk      <- as.integer(.tid_str)
    .res <- .nn_match_trial(
      dt_trial          = .dt_slice,
      bl_cm             = bl_lookup[[.tid_str]],
      length_unit_trial = lu_lookup[[.tid_str]]
    )
    if (nrow(.res) > 0L) .res[, trial_id := .tk]
    .res
  },
  .options = furrr::furrr_options(
    seed     = TRUE,
    globals  = c("bl_lookup", "lu_lookup", ".nn_match_trial"),
    packages = "data.table"
  )
)

nn_dt <- data.table::rbindlist(nn_results, fill = TRUE)
ts_msg("NN matching complete: ", nrow(nn_dt), " inter-frame transitions across ",
       data.table::uniqueN(nn_dt$trial_id), " trials")


# =============================================================================
# ==== 5) PER-FRAME ZONE COUNTS ===============================================
# =============================================================================

ts_msg("Computing per-frame zone counts...")

dt_zones <- data.table::as.data.table(df[, c("trial_id", "frame", "time",
                                             "timepoint", "main_zone", "sec_zone",
                                             "x_interp", "y_interp")])
dt_zones[, valid_pos := is.finite(x_interp) & is.finite(y_interp)]

# Per-frame counts
frame_counts <- dt_zones[valid_pos == TRUE,
  .(n_total       = .N,
    n_in_flow     = sum(main_zone == "flow", na.rm = TRUE),
    n_in_calm     = sum(main_zone == "calm", na.rm = TRUE),
    n_in_high     = sum(sec_zone  == "high", na.rm = TRUE),
    n_in_medium   = sum(sec_zone  == "medium", na.rm = TRUE),
    n_in_low      = sum(sec_zone  == "low", na.rm = TRUE),
    n_in_calm_sec = sum(sec_zone  == "calm", na.rm = TRUE),
    time = time[1]),
  by = .(trial_id, timepoint, frame)]

# dt_row per frame within (trial × timepoint)
data.table::setkey(frame_counts, trial_id, timepoint, frame)
frame_counts[, dt_row := c(0, diff(time)), by = .(trial_id, timepoint)]
frame_counts[!is.finite(dt_row) | dt_row < 0, dt_row := 0]

# Frame-level zone-flux input (for diagnostic correlation later)
frame_counts[, prop_flow_frame := n_in_flow / pmax(n_total, 1L)]
frame_counts[, dprop_flow := c(0, abs(diff(prop_flow_frame))),
             by = .(trial_id, timepoint)]


# =============================================================================
# ==== 6) AGGREGATE TO (trial_id, timepoint) ==================================
# =============================================================================

ts_msg("Aggregating to trial × timepoint...")

# Wide occupancy + total time
trial_wide <- frame_counts[ , .(
    total_time_s  = sum(dt_row, na.rm = TRUE),
    prop_time_in_flow   = stats::weighted.mean(n_in_flow     / pmax(n_total, 1L),
                                                w = dt_row, na.rm = TRUE),
    prop_time_in_calm   = stats::weighted.mean(n_in_calm     / pmax(n_total, 1L),
                                                w = dt_row, na.rm = TRUE),
    prop_time_in_high   = stats::weighted.mean(n_in_high     / pmax(n_total, 1L),
                                                w = dt_row, na.rm = TRUE),
    prop_time_in_medium = stats::weighted.mean(n_in_medium   / pmax(n_total, 1L),
                                                w = dt_row, na.rm = TRUE),
    prop_time_in_low    = stats::weighted.mean(n_in_low      / pmax(n_total, 1L),
                                                w = dt_row, na.rm = TRUE),
    prop_time_in_calm_sec = stats::weighted.mean(n_in_calm_sec / pmax(n_total, 1L),
                                                w = dt_row, na.rm = TRUE),
    # Zone-flux diagnostic: mean |Δprop_flow| per frame
    mean_dprop_flow = mean(dprop_flow, na.rm = TRUE),
    n_frames        = .N
  ), by = .(trial_id, timepoint)]

# Switches and prop_active from NN matching
nn_summary <- if (nrow(nn_dt) > 0L) {
  nn_dt[ , .(
    n_main_switches = sum(n_main_switches, na.rm = TRUE),
    n_active_sum    = sum(n_active, na.rm = TRUE),
    n_matched_sum   = sum(n_matched, na.rm = TRUE),
    obs_seconds     = sum(dt_s, na.rm = TRUE)
  ), by = .(trial_id, timepoint)]
} else {
  data.table::data.table(trial_id = integer(), timepoint = integer(),
                         n_main_switches = integer(),
                         n_active_sum = integer(), n_matched_sum = integer(),
                         obs_seconds = numeric())
}

# Merge
trial_activity_summary <- merge(trial_wide, nn_summary,
                                 by = c("trial_id", "timepoint"), all.x = TRUE)

# Session normalisation: each timepoint is a 20-min segment of the trial
SESSION_MINUTES <- as.numeric(.get_global("SESSION_MINUTES_PER_TIMEPOINT", 20))
trial_activity_summary[ , `:=`(
  prop_active = data.table::fifelse(is.finite(n_matched_sum) & n_matched_sum > 0,
                                    n_active_sum / n_matched_sum, NA_real_),
  switches_per_session   = data.table::fifelse(is.finite(obs_seconds) & obs_seconds > 0,
                                               n_main_switches *
                                                 (SESSION_MINUTES * 60) / obs_seconds,
                                               NA_real_),
  zone_flux_per_session  = data.table::fifelse(is.finite(mean_dprop_flow) & n_frames > 0,
                                               mean_dprop_flow * n_frames *
                                                 (SESSION_MINUTES * 60) /
                                                 pmax(total_time_s, 1e-9),
                                               NA_real_)
)]

# Attach per-trial bodylength + meta
.trial_meta_cols <- intersect(c("trial_id", "trial_date", "treatment", "fish_density",
                                "tank", "motor_side"), names(df))
.trial_meta <- df %>%
  dplyr::select(dplyr::all_of(.trial_meta_cols)) %>%
  dplyr::distinct(trial_id, .keep_all = TRUE)

trial_activity_summary <- as.data.frame(trial_activity_summary) %>%
  dplyr::left_join(bl_per_trial, by = "trial_id") %>%
  dplyr::left_join(.trial_meta,  by = "trial_id") %>%
  dplyr::mutate(
    # prop_active is meaningful only when a body-length threshold was available.
    # When bl_cm_trial is NA, n_active = 0 (all is_active were NA → sum(NA,na.rm=T)=0),
    # giving a spurious prop_active = 0 instead of NA.
    prop_active    = dplyr::if_else(is.finite(bl_cm_trial), prop_active, NA_real_),
    timepoint_f = factor(timepoint, levels = .get_global("TIMEPOINT_LEVELS", 1:3)),
    fish_density_f = if (isTRUE(.get_global("DENSITY_AS_FACTOR", TRUE))) {
      factor(fish_density,
             levels = .get_global("DENSITY_LEVELS", c(4L, 8L, 12L, 16L)),
             ordered = TRUE)
    } else as.numeric(fish_density)
  )

ts_msg("trial_activity_summary: ", nrow(trial_activity_summary),
       " rows (one per trial × timepoint)")


# =============================================================================
# ==== 6b) AREA CORRECTION FOR SUB-ZONES ======================================
# =============================================================================
# Sub-zones (high / medium / low / calm) differ in physical area; raw
# occupancy proportions must be normalised before comparing across zones.
#
# Zone areas in "units" where 1 unit = 10 × 25 cm section (from arena maps):
#   high   = 10 units   medium = 18 units
#   low    = 13 units   calm   = 42 units
#
# Procedure (matches published methods):
#   1. Divide prop_time_in_<zone> by its area in units.
#   2. Rescale so the four area-corrected values sum to 1
#      (divide each by their sum → proportions reflect preference, not area).
#
# Main zones (flow / calm) are equal in area and are NOT area-corrected.
# =============================================================================

.AREA_HIGH     <- 10
.AREA_MEDIUM   <- 18
.AREA_LOW      <- 13
.AREA_CALM_SEC <- 42

trial_activity_summary <- trial_activity_summary %>%
  dplyr::mutate(
    .norm_high   = prop_time_in_high     / .AREA_HIGH,
    .norm_medium = prop_time_in_medium   / .AREA_MEDIUM,
    .norm_low    = prop_time_in_low      / .AREA_LOW,
    .norm_calm   = prop_time_in_calm_sec / .AREA_CALM_SEC,
    .norm_total  = .norm_high + .norm_medium + .norm_low + .norm_calm,
    prop_time_ac_high   = dplyr::if_else(.norm_total > 0, .norm_high   / .norm_total, NA_real_),
    prop_time_ac_medium = dplyr::if_else(.norm_total > 0, .norm_medium / .norm_total, NA_real_),
    prop_time_ac_low    = dplyr::if_else(.norm_total > 0, .norm_low    / .norm_total, NA_real_),
    prop_time_ac_calm   = dplyr::if_else(.norm_total > 0, .norm_calm   / .norm_total, NA_real_)
  ) %>%
  dplyr::select(-.norm_high, -.norm_medium, -.norm_low, -.norm_calm, -.norm_total)

ts_msg("Area correction applied to sub-zones (high=10, medium=18, low=13, calm=42 units)")


# =============================================================================
# ==== 7) LONG-FORMAT OCCUPANCY (for joint main + sub stats) ==================
# =============================================================================
# Main zones  → use raw prop_time_in_flow / prop_time_in_calm   (equal areas)
# Sub-zones   → use area-corrected prop_time_ac_*               (unequal areas)

trial_occupancy_long <- trial_activity_summary %>%
  dplyr::select(trial_id, timepoint, timepoint_f, treatment,
                fish_density, fish_density_f, tank, motor_side, trial_date,
                prop_time_in_flow, prop_time_in_calm,
                prop_time_ac_high, prop_time_ac_medium,
                prop_time_ac_low,  prop_time_ac_calm) %>%
  tidyr::pivot_longer(
    cols      = c(prop_time_in_flow, prop_time_in_calm,
                  prop_time_ac_high, prop_time_ac_medium,
                  prop_time_ac_low,  prop_time_ac_calm),
    names_to  = "zone_raw",
    values_to = "prop_time"
  ) %>%
  dplyr::mutate(
    zone_level = dplyr::case_when(
      zone_raw %in% c("prop_time_in_flow", "prop_time_in_calm") ~ "main",
      TRUE                                                       ~ "sec"
    ),
    zone = dplyr::case_when(
      zone_raw == "prop_time_in_flow"   ~ "flow",
      zone_raw == "prop_time_in_calm"   ~ "calm",
      zone_raw == "prop_time_ac_high"   ~ "high",
      zone_raw == "prop_time_ac_medium" ~ "medium",
      zone_raw == "prop_time_ac_low"    ~ "low",
      zone_raw == "prop_time_ac_calm"   ~ "calm",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::select(-zone_raw) %>%
  dplyr::mutate(
    zone_level = factor(zone_level, levels = c("main", "sec")),
    zone       = factor(zone, levels = c("flow", "calm", "high", "medium", "low"))
  )

ts_msg("trial_occupancy_long: ", nrow(trial_occupancy_long), " rows")


# =============================================================================
# ==== 8) SAVE OUTPUTS ========================================================
# =============================================================================

readr::write_csv(trial_activity_summary,
                 file.path(STEP2_OUT, "trial_activity_summary.csv"))
ts_msg("trial_activity_summary.csv written")

readr::write_csv(trial_occupancy_long,
                 file.path(STEP2_OUT, "trial_occupancy_long.csv"))
ts_msg("trial_occupancy_long.csv written")

readr::write_csv(bl_per_trial,
                 file.path(STEP2_OUT, "diagnostics", "bodylength_per_trial.csv"))

# Diagnostic: NN-matching coverage
if (nrow(nn_dt) > 0) {
  nn_diag <- nn_dt[ , .(
    n_transitions = .N,
    mean_match_rate = mean(n_matched / pmax(n_total_from, 1)),
    median_dt_s     = stats::median(dt_s, na.rm = TRUE)
  ), by = .(trial_id, timepoint)]
  readr::write_csv(nn_diag, file.path(STEP2_OUT, "diagnostics",
                                       "nn_matching_coverage.csv"))
}


# =============================================================================
# ==== 9) UPDATE GLOBALENV ====================================================
# =============================================================================

assign("trial_activity_summary", trial_activity_summary, envir = .GlobalEnv)
assign("trial_occupancy_long",   trial_occupancy_long,   envir = .GlobalEnv)
assign("bl_per_trial",           bl_per_trial,           envir = .GlobalEnv)
assign("STEP2_OUTPUT_DIR",       STEP2_OUT,              envir = .GlobalEnv)

# Backwards-compatibility alias for downstream scripts that still reference
# fish_activity_summary by name (graceful deprecation): point it at the new
# trial-level table. Anything that was filtering on fish_id will fail loudly,
# which is the desired signal that those callers need updating.
assign("fish_activity_summary", trial_activity_summary, envir = .GlobalEnv)

ts_msg("Step 2 complete — output: ", STEP2_OUT)
ts_msg("  trial_activity_summary in GlobalEnv (", nrow(trial_activity_summary), " rows)")
ts_msg("  trial_occupancy_long in GlobalEnv (", nrow(trial_occupancy_long), " rows)")
