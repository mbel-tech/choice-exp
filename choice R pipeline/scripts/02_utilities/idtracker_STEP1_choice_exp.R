# =============================================================================
# STEP 1 — CHOICE EXPERIMENT
# Read raw idtracker.ai trajectory CSVs; detect and correct jumps / identity
# switches; assign multi-level zones; attach treatment metadata.
# =============================================================================
#
# This script reads idtracker.ai output directly from PIPELINE_DATA_DIR:
#   1. Reads trajectories_csv/trajectories.csv from every session sub-folder,
#      interpolates NaN gaps, and assembles master_fish_by_frame.
#   2. Detects tracking jumps (section 2b) and identity switches (section 2c)
#      and corrects them in-place.
#   3. Parses trial number and timepoint from session folder names
#      (N{trial}_segment{timepoint} pattern).
#   4. Joins TRIAL_META to add treatment, fish_density, tank, motor_side,
#      trial_date to every row.
#   5. Assigns multi-level zones via assign_zones_choice_exp().
#   6. Writes master_fish_by_frame.csv to a timestamped STEP1_output/ folder.
#
# PRE-CONDITIONS (supplied by 00_master_pipeline_choice_exp.R):
#   - PIPELINE_DATA_DIR          path to session data folder
#   - TRIAL_META                 data.frame loaded from trial_summary Excel
#   - ZONE_REF_FT                data.frame loaded from zone_reference_choice_exp_FT.csv
#   - ZONE_REF_FD                data.frame loaded from zone_reference_choice_exp_FD.csv
#   - assign_zones_choice_exp()  polygon PIP zone-assignment function
#   - DENSITY_AS_FACTOR          logical
#   - TREATMENT_LEVELS           character vector
#   - DENSITY_LEVELS             integer vector
# =============================================================================


# =============================================================================
# ==== 0) PACKAGES + HELPERS ==================================================
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(readxl)
  library(readr)
  library(stringr)
  library(jsonlite)
  library(furrr)
  library(future)
})

ts_msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), " | ", ..., "\n", sep = "")

# Resolve the directory of THIS script at source() time (same trick as master)
.step1_dir <- tryCatch({
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

# Helper: safely get object from GlobalEnv
.get_global <- function(name, default = NULL) {
  if (exists(name, envir = .GlobalEnv, inherits = FALSE)) {
    get(name, envir = .GlobalEnv, inherits = FALSE)
  } else {
    default
  }
}

# Helper: read one idtracker.ai session folder → long-format data.frame
# Returns NULL (with a warning) if the trajectory CSV is absent or unreadable.
.read_one_session <- function(sess_folder, trial_id_seq, base_dir) {
  traj_csv  <- file.path(base_dir, sess_folder,
                          "trajectories", "trajectories_csv", "trajectories.csv")
  attr_json <- file.path(base_dir, sess_folder,
                          "trajectories", "trajectories_csv", "attributes.json")

  if (!file.exists(traj_csv)) {
    warning("trajectories.csv not found in session: ", sess_folder, call. = FALSE)
    return(NULL)
  }

  raw <- tryCatch(
    readr::read_csv(traj_csv, show_col_types = FALSE, progress = FALSE),
    error = function(e) {
      warning("Cannot read trajectories.csv for ", sess_folder, ": ",
              conditionMessage(e), call. = FALSE)
      NULL
    }
  )
  if (is.null(raw) || nrow(raw) == 0) return(NULL)

  # Detect fish count: columns are time, x1, y1, x2, y2, ...
  n_fish <- (ncol(raw) - 1L) %/% 2L
  if (n_fish < 1L) {
    warning("No fish columns found in: ", sess_folder, call. = FALSE)
    return(NULL)
  }

  # Read length_unit from attributes.json (pixels per cm calibration)
  length_unit <- NA_real_
  if (file.exists(attr_json)) {
    att <- tryCatch(jsonlite::read_json(attr_json), error = function(e) NULL)
    if (!is.null(att) && !is.null(att$length_unit)) {
      lu <- suppressWarnings(as.numeric(att$length_unit[[1]]))
      if (is.finite(lu) && lu > 0) length_unit <- lu
    }
  }

  frame_vec <- seq_len(nrow(raw))
  time_vec  <- if ("time" %in% names(raw)) raw[["time"]] else NA_real_

  long_list <- vector("list", n_fish)
  for (f in seq_len(n_fish)) {
    xcol <- paste0("x", f)
    ycol <- paste0("y", f)
    if (!xcol %in% names(raw) || !ycol %in% names(raw)) next

    x_raw <- as.numeric(raw[[xcol]])
    y_raw <- as.numeric(raw[[ycol]])

    # Linear interpolation over NaN/NA gaps (idtracker.ai uses NaN for missed frames)
    x_interp <- x_raw
    y_interp <- y_raw
    vld <- is.finite(x_raw) & is.finite(y_raw)
    if (sum(vld) >= 2L) {
      fi_vld   <- frame_vec[vld]
      x_interp <- approx(fi_vld, x_raw[vld], xout = frame_vec,
                          method = "linear", rule = 1L)$y
      y_interp <- approx(fi_vld, y_raw[vld], xout = frame_vec,
                          method = "linear", rule = 1L)$y
    }

    long_list[[f]] <- data.frame(
      trial_id    = trial_id_seq,
      fish_id     = as.character(f),
      frame       = frame_vec,
      time        = time_vec,
      x_interp    = x_interp,
      y_interp    = y_interp,
      length_unit = length_unit,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, Filter(Negate(is.null), long_list))
}

# Parallel plan — inherit from master pipeline; fall back to sequential
.n_workers_step1 <- .get_global("PIPELINE_WORKERS", 1L)
if (.n_workers_step1 > 1L && inherits(future::plan(), "SequentialFuture"))
  future::plan(future::multisession, workers = .n_workers_step1)


# =============================================================================
# ==== 1) LOAD TRIAL_META IF NOT ALREADY IN GLOBALENV =========================
# =============================================================================
# The master pipeline loads TRIAL_META before running steps.
# If Step 1 is run standalone (START_FROM_STEP = 1 after fresh R session),
# we load TRIAL_META here from the Excel file next to this script.

if (!exists("TRIAL_META", envir = .GlobalEnv, inherits = FALSE)) {
  ts_msg("TRIAL_META not found in GlobalEnv; loading from Excel...")

  .trial_excel <- file.path(.step1_dir, "trial_summary_choice_exp.xlsx")
  if (!file.exists(.trial_excel)) {
    stop("trial_summary_choice_exp.xlsx not found at: ", .trial_excel, call. = FALSE)
  }

  .tm_raw <- readxl::read_excel(.trial_excel)
  if ("condition" %in% names(.tm_raw) && !("treatment" %in% names(.tm_raw))) {
    .tm_raw <- dplyr::rename(.tm_raw, treatment = condition)
  }

  TRIAL_META <- .tm_raw %>%
    dplyr::mutate(
      trial        = suppressWarnings(as.integer(trial)),
      fish_density = suppressWarnings(as.integer(fish_density)),
      tank         = as.character(tank),
      motor_side   = trimws(as.character(motor_side)),
      treatment    = trimws(tolower(as.character(treatment)))
    )

  if ("video_ID" %in% names(TRIAL_META))
    TRIAL_META$video_ID <- as.character(TRIAL_META$video_ID)

  # Rename 'date' -> 'trial_date'
  if (!"trial_date" %in% names(TRIAL_META)) {
    .dc <- intersect(c("date", "Date", "DATE"), names(TRIAL_META))[1]
    if (!is.na(.dc)) TRIAL_META <- dplyr::rename(TRIAL_META, trial_date = !!.dc)
  }

  assign("TRIAL_META", TRIAL_META, envir = .GlobalEnv)
  ts_msg("TRIAL_META loaded: ", nrow(TRIAL_META), " rows")
}

TRIAL_META <- .get_global("TRIAL_META")


# =============================================================================
# ==== 1b) INCREMENTAL SESSION FILTER =========================================
# =============================================================================
# When SKIP_PROCESSED_SESSIONS = TRUE (set in 00_master_pipeline_choice_exp.R):
#   - Scans PIPELINE_DATA_DIR for session sub-folders.
#   - Reads a manifest (STEP1_output/processed_sessions.txt) listing folders
#     already processed in previous runs.
#   - If all sessions are already processed: loads master_fish_by_frame.csv
#     from disk (section 2 will then skip the trajectory read).
#   - If new sessions exist: reads trajectory CSVs for those sessions only,
#     assigns trial_ids based on global alphabetical sort (ensuring stable IDs
#     across incremental runs), then merges with the existing CSV.
#   - Updates the manifest after a successful run.
#
# When SKIP_PROCESSED_SESSIONS = FALSE (default): this block is skipped and
# section 2 reads all sessions from scratch.
# =============================================================================

.skip_processed <- isTRUE(.get_global("SKIP_PROCESSED_SESSIONS", FALSE))
.data_dir_incr  <- .get_global("PIPELINE_DATA_DIR")
.manifest_path  <- file.path(getwd(), "STEP1_output", "processed_sessions.txt")

# Helper: find the latest master_fish_by_frame.csv in STEP1_output/
.find_latest_mfbf_csv <- function() {
  parent <- file.path(getwd(), "STEP1_output")
  if (!dir.exists(parent)) return(NULL)
  subdirs <- list.dirs(parent, full.names = TRUE, recursive = FALSE)
  subdirs <- subdirs[grepl("^STEP1_output_\\d{8}_\\d{6}$", basename(subdirs))]
  if (length(subdirs) == 0) return(NULL)
  latest  <- subdirs[which.max(file.mtime(subdirs))]
  f <- file.path(latest, "master_fish_by_frame.csv")
  if (file.exists(f)) f else NULL
}

if (.skip_processed && !is.null(.data_dir_incr) && dir.exists(.data_dir_incr)) {

  # ---- Identify N{n}_segment{tp} folders in PIPELINE_DATA_DIR ---------------
  .all_sessions <- sort(list.dirs(.data_dir_incr, recursive = FALSE, full.names = FALSE))
  .all_sessions <- .all_sessions[
    grepl("N\\d+", .all_sessions) & grepl("seg(?:ment)?_?\\d+", .all_sessions, ignore.case = TRUE, perl = TRUE)
  ]

  if (length(.all_sessions) == 0) {
    warning(
      "SKIP_PROCESSED_SESSIONS = TRUE but no 'N{n}_segment{tp}' sub-folders found in:\n  ",
      .data_dir_incr, "\n",
      "Proceeding with normal full-pipeline run.",
      call. = FALSE
    )
  } else {

    # ---- Load manifest -------------------------------------------------------
    .done_sessions <- if (file.exists(.manifest_path))
      readLines(.manifest_path, warn = FALSE) else character(0)
    .done_sessions <- .done_sessions[nzchar(.done_sessions)]

    .new_sessions <- setdiff(.all_sessions, .done_sessions)

    ts_msg("Incremental mode: ", length(.all_sessions), " session folder(s) found | ",
           length(.done_sessions), " in manifest | ",
           length(.new_sessions), " new")

    if (length(.new_sessions) == 0) {
      # ---- All sessions already processed → load from disk ------------------
      ts_msg("All sessions already processed. Loading master_fish_by_frame from disk.")
      .existing_csv <- .find_latest_mfbf_csv()
      if (!is.null(.existing_csv)) {
        .mfbf_from_disk <- tryCatch(
          readr::read_csv(.existing_csv, show_col_types = FALSE),
          error = function(e) {
            warning("Could not read ", .existing_csv, ": ", conditionMessage(e), call. = FALSE)
            NULL
          }
        )
        if (!is.null(.mfbf_from_disk)) {
          assign("master_fish_by_frame", .mfbf_from_disk, envir = .GlobalEnv)
          ts_msg("Loaded ", nrow(.mfbf_from_disk), " rows from: ", .existing_csv)
        }
      } else {
        warning(
          "No master_fish_by_frame.csv found in STEP1_output/; ",
          "will run full trajectory read on all sessions.",
          call. = FALSE
        )
      }

    } else {
      # ---- New sessions: read trajectory CSVs using global trial_id order ---
      ts_msg("Reading trajectory CSVs for ", length(.new_sessions), " new session(s)...")

      # Global alphabetical sort → stable trial_ids regardless of when sessions
      # are added, so IDs match across incremental runs.
      .global_id_map <- data.frame(
        folder   = .all_sessions,
        trial_id = seq_along(.all_sessions),
        stringsAsFactors = FALSE
      )

      .new_data_list <- lapply(.new_sessions, function(.sf) {
        .tid <- .global_id_map$trial_id[.global_id_map$folder == .sf]
        .read_one_session(.sf, .tid, .data_dir_incr)
      })
      .new_mfbf <- do.call(rbind, Filter(Negate(is.null), .new_data_list))

      if (is.null(.new_mfbf) || nrow(.new_mfbf) == 0) {
        warning(
          "Could not read trajectory data from new session(s). ",
          "Falling back to full-pipeline run.",
          call. = FALSE
        )
      } else {
        ts_msg("New sessions produced ", nrow(.new_mfbf), " rows.")

        # Attempt to merge with existing CSV on disk
        .existing_csv <- .find_latest_mfbf_csv()
        if (!is.null(.existing_csv)) {
          .old_mfbf <- tryCatch(
            readr::read_csv(.existing_csv, show_col_types = FALSE),
            error = function(e) NULL
          )
          if (!is.null(.old_mfbf)) {
            .merged <- dplyr::bind_rows(.old_mfbf, .new_mfbf)
            assign("master_fish_by_frame", .merged, envir = .GlobalEnv)
            ts_msg("Merged: ", nrow(.old_mfbf), " existing + ",
                   nrow(.new_mfbf), " new = ", nrow(.merged), " total rows.")
          } else {
            assign("master_fish_by_frame", .new_mfbf, envir = .GlobalEnv)
          }
        } else {
          assign("master_fish_by_frame", .new_mfbf, envir = .GlobalEnv)
        }

        dir.create(dirname(.manifest_path), recursive = TRUE, showWarnings = FALSE)
        writeLines(c(.done_sessions, .new_sessions), .manifest_path)
        ts_msg("Manifest updated: ", .manifest_path)
      }
    }
  }
} else if (.skip_processed) {
  warning(
    "SKIP_PROCESSED_SESSIONS = TRUE but PIPELINE_DATA_DIR is not available or does not exist.",
    call. = FALSE
  )
}

# After section 5 completes we update the manifest for SKIP_PROCESSED_SESSIONS = FALSE too.
# That update is in section 5 below (see .manifest_update block).


# =============================================================================
# ==== 2) BUILD master_fish_by_frame FROM IDTRACKER.AI TRAJECTORY CSVs ========
# =============================================================================
# Reads trajectories_csv/trajectories.csv from every session sub-folder in
# PIPELINE_DATA_DIR, fills NaN gaps by linear interpolation, and assembles
# master_fish_by_frame in long format (one row per fish × frame).
#
# Session sub-folders must contain "N{n}" and "segment{n}" in their name.
# trial_id is assigned in alphabetical sort order (1 = first folder, etc.).
# length_unit (pixels per cm) is read from trajectories_csv/attributes.json.

if (!exists("master_fish_by_frame", envir = .GlobalEnv, inherits = FALSE)) {
  ts_msg("master_fish_by_frame not in GlobalEnv; reading trajectory CSV files...")

  .data_dir_s2 <- .get_global("PIPELINE_DATA_DIR")
  if (is.null(.data_dir_s2) || !dir.exists(.data_dir_s2))
    stop("PIPELINE_DATA_DIR is not set or does not exist.", call. = FALSE)

  .all_sess_s2 <- sort(list.dirs(.data_dir_s2, recursive = FALSE, full.names = FALSE))
  .all_sess_s2 <- .all_sess_s2[
    grepl("N\\d+",       .all_sess_s2, perl = TRUE) &
    grepl("seg(?:ment)?_?\\d+", .all_sess_s2, ignore.case = TRUE, perl = TRUE)
  ]

  if (length(.all_sess_s2) == 0)
    stop("No session folders matching 'N{n}...segment{tp}' found in:\n  ",
         .data_dir_s2, call. = FALSE)

  ts_msg("  Reading ", length(.all_sess_s2), " session folder(s)...")

  .mfbf_list <- vector("list", length(.all_sess_s2))
  for (.si in seq_along(.all_sess_s2)) {
    .mfbf_list[[.si]] <- .read_one_session(.all_sess_s2[.si], .si, .data_dir_s2)
    if (is.null(.mfbf_list[[.si]]))
      warning("No data read from session: ", .all_sess_s2[.si], call. = FALSE)
  }

  .mfbf_all <- do.call(rbind, Filter(Negate(is.null), .mfbf_list))
  if (is.null(.mfbf_all) || nrow(.mfbf_all) == 0)
    stop("No trajectory data could be read from any session folder.", call. = FALSE)

  assign("master_fish_by_frame", .mfbf_all, envir = .GlobalEnv)
  ts_msg("master_fish_by_frame assembled: ", nrow(.mfbf_all), " rows")

} else {
  ts_msg("master_fish_by_frame already in GlobalEnv; skipping trajectory read.")
}

# Confirm the object exists
if (!exists("master_fish_by_frame", envir = .GlobalEnv, inherits = FALSE))
  stop("master_fish_by_frame was not created. Check for errors above.", call. = FALSE)

mfbf <- get("master_fish_by_frame", envir = .GlobalEnv, inherits = FALSE)
ts_msg("master_fish_by_frame: ", nrow(mfbf), " rows x ", ncol(mfbf), " columns")


# =============================================================================
# ==== 2b) JUMP DETECTION AND LINEAR INTERPOLATION ============================
# =============================================================================
# Detects frame-to-frame position jumps caused by tracking errors (identity
# swaps, lost tracks) and replaces the offending frame coordinates with
# linearly interpolated positions.
#
# Algorithm mirrors the idtracker.ai validator (errors_explorer.py):
#   1. Compute Euclidean distance between consecutive valid positions for every
#      fish in every trial, divided by the frame-gap so that tracks with
#      missing frames are not penalised.  Units: pixels per frame.
#   2. Compute a global threshold from the pooled speed distribution.
#   3. Any frame whose incoming speed exceeds the threshold is marked as a jump.
#   4. Jump-frame coordinates are set to NA.
#   5. All NA positions are filled by linear interpolation between the nearest
#      valid frames before and after the gap (base R approx(), rule = 1:
#      edges that have no valid anchor on one side are left as NA).
#
# CONFIG (set in 00_master_pipeline_choice_exp.R):
#   JUMP_DETECTION_ENABLED   TRUE/FALSE
#   JUMP_THRESHOLD_METHOD    "sd_multiple" (mean + N×sd) | "percentile" (M×Pth)
#   JUMP_SD_MULT             N multiplier for sd_multiple (default 10, matching
#                            the idtracker.ai validator GUI default)
#   JUMP_PCT / JUMP_PCT_MULT percentile settings for percentile method
#   JUMP_INTERP              TRUE = fill gap; FALSE = leave as NA
#
# OUTPUTS (written to STEP1_output/<run>/):
#   jump_detection_summary.csv   one row per fish-track; columns:
#     trial_id, fish_id, n_frames, n_jumps, pct_jumps, threshold_px_per_frame
#   Columns added to master_fish_by_frame:
#     x_interp_raw, y_interp_raw  — original (pre-filter) coordinates
#     jump_flag                   — logical; TRUE on frames removed as jumps
# =============================================================================

.jump_enabled <- isTRUE(.get_global("JUMP_DETECTION_ENABLED", TRUE))

if (.jump_enabled) {
  ts_msg("Jump detection: running...")

  .jmethod   <- .get_global("JUMP_THRESHOLD_METHOD", "sd_multiple")
  .jsd_mult  <- .get_global("JUMP_SD_MULT",          10.0)
  .jpct      <- .get_global("JUMP_PCT",              99)
  .jpct_mult <- .get_global("JUMP_PCT_MULT",         2.0)
  .jinterp   <- isTRUE(.get_global("JUMP_INTERP",   TRUE))

  # Identify position and frame columns
  .x_col  <- if ("x_interp" %in% names(mfbf)) "x_interp" else
    intersect(c("x","X","pos_x"), names(mfbf))[1]
  .y_col  <- if ("y_interp" %in% names(mfbf)) "y_interp" else
    intersect(c("y","Y","pos_y"), names(mfbf))[1]
  .fr_col <- if ("frame"    %in% names(mfbf)) "frame"    else
    intersect(c("frame_number","Frame","t"), names(mfbf))[1]
  .id_col <- if ("fish_id"  %in% names(mfbf)) "fish_id"  else
    intersect(c("identity","animal_id","id"), names(mfbf))[1]

  if (is.na(.x_col) || is.na(.y_col) || is.na(.fr_col) || is.na(.id_col)) {
    warning(
      "Jump detection skipped: could not locate required columns.\n",
      "Expected: x_interp / y_interp / frame / fish_id in master_fish_by_frame.",
      call. = FALSE
    )
    .jump_enabled <- FALSE
  }
}

if (.jump_enabled) {

  # ---- Preserve original coordinates ----------------------------------------
  mfbf$x_interp_raw <- mfbf[[.x_col]]
  mfbf$y_interp_raw <- mfbf[[.y_col]]
  mfbf$jump_flag    <- FALSE

  # ---- Step 1: collect all per-frame speeds to compute the global threshold --
  .grp_key <- paste(mfbf$trial_id, mfbf[[.id_col]], sep = "__")
  .groups  <- split(seq_len(nrow(mfbf)), .grp_key)

  .all_speeds <- unlist(lapply(.groups, function(idx) {
    xi <- mfbf[[.x_col]][idx]
    yi <- mfbf[[.y_col]][idx]
    fi <- mfbf[[.fr_col]][idx]
    ord <- order(fi)
    xi <- xi[ord]; yi <- yi[ord]; fi <- fi[ord]
    vld <- is.finite(xi) & is.finite(yi)
    if (sum(vld) < 4L) return(numeric(0))
    vi <- which(vld)
    dx <- diff(xi[vi]); dy <- diff(yi[vi])
    fg <- pmax(diff(fi[vi]), 1L)
    spd <- sqrt(dx^2 + dy^2) / fg
    spd[is.finite(spd)]
  }))

  if (length(.all_speeds) < 20L) {
    warning("Jump detection: too few valid speed values; skipping.", call. = FALSE)
    .jump_enabled <- FALSE
  } else {
    .jthresh <- if (.jmethod == "percentile") {
      .jpct_mult * quantile(.all_speeds, .jpct / 100, na.rm = TRUE)
    } else {
      mean(.all_speeds, na.rm = TRUE) + .jsd_mult * sd(.all_speeds, na.rm = TRUE)
    }
    ts_msg("  method=", .jmethod,
           "  threshold=", round(.jthresh, 2), " px/frame",
           "  (", .jmethod,
           if (.jmethod == "sd_multiple")
             paste0(": mean=", round(mean(.all_speeds, na.rm=TRUE), 2),
                    " sd=", round(sd(.all_speeds, na.rm=TRUE), 2),
                    " mult=", .jsd_mult)
           else
             paste0(": ", .jpct, "th pct=",
                    round(quantile(.all_speeds, .jpct/100, na.rm=TRUE), 2),
                    " mult=", .jpct_mult),
           ")")
  }
}

if (.jump_enabled) {

  # ---- Step 2: detect and replace per track ----------------------------------
  # Pre-extract per-track slices so workers do not serialise the full mfbf.
  .jump_inputs <- lapply(names(.groups), function(.g) {
    .idx <- .groups[[.g]]
    list(
      g        = .g,
      idx      = .idx,
      xi       = mfbf[[.x_col]][.idx],
      yi       = mfbf[[.y_col]][.idx],
      fi       = mfbf[[.fr_col]][.idx],
      trial_id = mfbf$trial_id[.idx[1L]],
      fish_id  = mfbf[[.id_col]][.idx[1L]]
    )
  })
  names(.jump_inputs) <- names(.groups)

  .jump_par_results <- furrr::future_map(
    .jump_inputs,
    function(.inp) {
      if (length(.inp$idx) < 6L) return(NULL)

      .idx <- .inp$idx
      .xi  <- .inp$xi; .yi <- .inp$yi; .fi <- .inp$fi

      .ord <- order(.fi)
      .idx <- .idx[.ord]
      .xi  <- .xi[.ord]; .yi <- .yi[.ord]; .fi <- .fi[.ord]

      .vld <- is.finite(.xi) & is.finite(.yi)
      if (sum(.vld) < 4L) return(NULL)

      .vi  <- which(.vld)
      .dx  <- diff(.xi[.vi]); .dy <- diff(.yi[.vi])
      .fg  <- pmax(diff(.fi[.vi]), 1L)
      .spd <- sqrt(.dx^2 + .dy^2) / .fg

      .jmask <- logical(length(.xi))
      for (.k in seq_along(.spd)) {
        if (is.finite(.spd[.k]) && .spd[.k] > .jthresh)
          .jmask[.vi[.k + 1L]] <- TRUE
      }

      .n_jumps <- sum(.jmask)

      if (.n_jumps > 0L) {
        .xi[.jmask] <- NA_real_
        .yi[.jmask] <- NA_real_

        if (.jinterp) {
          .vf <- .fi[is.finite(.xi)]
          .vx <- .xi[is.finite(.xi)]
          .vy <- .yi[is.finite(.xi)]
          if (length(.vf) >= 2L) {
            .xi <- approx(.vf, .vx, xout = .fi, method = "linear", rule = 1L)$y
            .yi <- approx(.vf, .vy, xout = .fi, method = "linear", rule = 1L)$y
          }
        }
      }

      list(
        g        = .inp$g,
        idx      = .idx,
        xi       = .xi,
        yi       = .yi,
        jmask    = .jmask,
        n_frames = length(.idx),
        n_jumps  = .n_jumps,
        trial_id = .inp$trial_id,
        fish_id  = .inp$fish_id
      )
    },
    .options = furrr::furrr_options(
      seed     = TRUE,
      globals  = c(".jthresh", ".jinterp"),
      packages = character(0)
    )
  )
  names(.jump_par_results) <- names(.groups)

  # Apply corrections serially (each track's indices are disjoint in mfbf)
  .jump_summary_rows <- vector("list", length(.groups))
  names(.jump_summary_rows) <- names(.groups)

  for (.g in names(.groups)) {
    .r <- .jump_par_results[[.g]]
    if (is.null(.r)) next

    if (.r$n_jumps > 0L) {
      mfbf[[.x_col]][.r$idx]              <- .r$xi
      mfbf[[.y_col]][.r$idx]              <- .r$yi
      mfbf$jump_flag[.r$idx[.r$jmask]]   <- TRUE
    }

    .jump_summary_rows[[.g]] <- data.frame(
      group_key              = .r$g,
      trial_id               = .r$trial_id,
      fish_id                = .r$fish_id,
      n_frames               = .r$n_frames,
      n_jumps                = .r$n_jumps,
      pct_jumps              = round(100 * .r$n_jumps / .r$n_frames, 2),
      threshold_px_per_frame = round(.jthresh, 3),
      stringsAsFactors       = FALSE
    )
  }

  # ---- Build summary and report ----------------------------------------------
  .jump_df <- do.call(rbind, Filter(Negate(is.null), .jump_summary_rows))
  rownames(.jump_df) <- NULL

  .n_total_jumps  <- sum(.jump_df$n_jumps,  na.rm = TRUE)
  .n_total_frames <- sum(.jump_df$n_frames, na.rm = TRUE)
  .n_fish_affected <- sum(.jump_df$n_jumps > 0, na.rm = TRUE)

  ts_msg("  Jumps removed: ", .n_total_jumps,
         " frames across ", .n_fish_affected, " fish-tracks",
         " (", round(100 * .n_total_jumps / max(.n_total_frames, 1), 3), "% of all frames)")

  if (.n_total_jumps > 0.05 * .n_total_frames) {
    warning(
      round(100 * .n_total_jumps / .n_total_frames, 1),
      "% of frames removed as jumps — unusually high. ",
      "Consider increasing JUMP_SD_MULT or switching JUMP_THRESHOLD_METHOD to ",
      "'percentile' to use a less aggressive threshold.",
      call. = FALSE
    )
  }

  assign("jump_detection_summary", .jump_df, envir = .GlobalEnv)

} else {
  # Jump detection disabled or skipped: add placeholder columns for consistency
  if (!"x_interp_raw" %in% names(mfbf)) mfbf$x_interp_raw <- mfbf[[if ("x_interp" %in% names(mfbf)) "x_interp" else names(mfbf)[1]]]
  if (!"y_interp_raw" %in% names(mfbf)) mfbf$y_interp_raw <- mfbf[[if ("y_interp" %in% names(mfbf)) "y_interp" else names(mfbf)[2]]]
  if (!"jump_flag"    %in% names(mfbf)) mfbf$jump_flag    <- FALSE
  .jump_df <- data.frame()
  assign("jump_detection_summary", .jump_df, envir = .GlobalEnv)
}

ts_msg("Jump detection complete.")


# =============================================================================
# ==== 2c) IDENTITY SWITCH DETECTION AND CORRECTION ===========================
# =============================================================================
# Detects frames where the tracker swapped two fish identities (typically after
# a crossing event) and corrects the position columns by swapping them back.
#
# Runs on jump-corrected coordinates (x_interp / y_interp from section 2b).
# Produces:
#   switch_flag column in mfbf     — TRUE on frames whose coordinates were
#                                    swapped by correction
#   identity_switch_events object  — event table stored in GlobalEnv; written
#                                    to disk in section 5 (with timepoint info)
#
# See 00_master_pipeline_choice_exp.R (CONFIG block) for parameter descriptions.
# =============================================================================

.sw_enabled <- isTRUE(.get_global("SWITCH_DETECTION_ENABLED", TRUE))

if (.sw_enabled) {
  ts_msg("Switch detection: running...")

  .sw_ratio   <- .get_global("SWITCH_TIER1_RATIO",   1.5)
  .sw_min_dur <- as.integer(.get_global("SWITCH_MIN_DURATION", 2L))
  .sw_gap     <- as.integer(.get_global("SWITCH_CONSOL_GAP",   5L))
  .sw_hung    <- isTRUE(.get_global("SWITCH_USE_HUNGARIAN",  TRUE))
  .sw_prox    <- .get_global("SWITCH_PROXIMITY_PX",  0)

  # Reuse column names from section 2b when available
  if (!exists(".x_col",  inherits = FALSE))
    .x_col  <- if ("x_interp" %in% names(mfbf)) "x_interp" else names(mfbf)[1]
  if (!exists(".y_col",  inherits = FALSE))
    .y_col  <- if ("y_interp" %in% names(mfbf)) "y_interp" else names(mfbf)[2]
  if (!exists(".fr_col", inherits = FALSE))
    .fr_col <- if ("frame"    %in% names(mfbf)) "frame"    else names(mfbf)[3]
  if (!exists(".id_col", inherits = FALSE))
    .id_col <- if ("fish_id"  %in% names(mfbf)) "fish_id"  else names(mfbf)[4]

  .hung_avail <- .sw_hung && requireNamespace("clue", quietly = TRUE)
  if (.sw_hung && !.hung_avail)
    ts_msg("  clue not found — Tier 2 (Hungarian) disabled for switch detection")

  # Columns to swap during correction (position + associated derived columns)
  .sw_swap_cols <- intersect(
    c(.x_col, .y_col, "x_interp_raw", "y_interp_raw", "jump_flag"),
    names(mfbf)
  )

  # Save pre-correction positions to compute switch_flag accurately later
  .x_pre_sw <- mfbf[[.x_col]]
  .y_pre_sw <- mfbf[[.y_col]]

  # ---- DETECTION --------------------------------------------------------
  # Pre-extract per-trial slices; workers receive only ~30 K rows each.
  .tid_vec <- sort(unique(mfbf$trial_id))
  .sw_td_list <- lapply(.tid_vec, function(.tid) {
    .t_idx <- which(mfbf$trial_id == .tid)
    if (length(.t_idx) < 10L) return(NULL)
    df <- data.frame(
      frame    = mfbf[[.fr_col]][.t_idx],
      fish_id  = as.character(mfbf[[.id_col]][.t_idx]),
      x        = mfbf[[.x_col]][.t_idx],
      y        = mfbf[[.y_col]][.t_idx],
      trial_id = .tid,
      stringsAsFactors = FALSE
    )
    df[order(df$frame, df$fish_id), ]
  })
  names(.sw_td_list) <- as.character(.tid_vec)

  .sw_par_results <- furrr::future_map(
    .sw_td_list,
    function(.td) {
      if (is.null(.td)) return(list())
      .tid <- .td$trial_id[1L]

      .fish_ids <- sort(unique(.td$fish_id))
      .nf_ids   <- length(.fish_ids)
      if (.nf_ids < 2L) return(list())

      .frames <- sort(unique(.td$frame))
      .nf     <- length(.frames)
      if (.nf < 4L) return(list())

      .wx <- matrix(NA_real_, .nf, .nf_ids, dimnames = list(NULL, .fish_ids))
      .wy <- matrix(NA_real_, .nf, .nf_ids, dimnames = list(NULL, .fish_ids))
      for (.fi in seq_along(.fish_ids)) {
        .frows <- .td[.td$fish_id == .fish_ids[.fi], ]
        .fidx  <- match(.frows$frame, .frames)
        .wx[.fidx, .fi] <- .frows$x
        .wy[.fidx, .fi] <- .frows$y
      }

      .hung_conf_frame <- rep(FALSE, .nf - 1L)
      .pairs <- combn(.nf_ids, 2L, simplify = FALSE)
      .trial_events <- list()

      for (.p in .pairs) {
        .i <- .p[1]; .j <- .p[2]
        .fid_a <- .fish_ids[.i]; .fid_b <- .fish_ids[.j]

        .xi <- .wx[, .i]; .xj <- .wx[, .j]
        .yi <- .wy[, .i]; .yj <- .wy[, .j]

        .prox_dist <- sqrt((.xi[-.nf] - .xj[-.nf])^2 + (.yi[-.nf] - .yj[-.nf])^2)
        .prox_ok   <- if (.sw_prox > 0) .prox_dist <= .sw_prox else rep(TRUE, .nf - 1L)

        .d_si  <- sqrt((.xi[-1L] - .xi[-.nf])^2 + (.yi[-1L] - .yi[-.nf])^2)
        .d_sj  <- sqrt((.xj[-1L] - .xj[-.nf])^2 + (.yj[-1L] - .yj[-.nf])^2)
        .d_ci  <- sqrt((.xi[-1L] - .xj[-.nf])^2 + (.yi[-1L] - .yj[-.nf])^2)
        .d_cj  <- sqrt((.xj[-1L] - .xi[-.nf])^2 + (.yj[-1L] - .yi[-.nf])^2)

        .d_curr  <- .d_si + .d_sj
        .d_swap  <- .d_ci + .d_cj
        .ratio_v <- .d_curr / pmax(.d_swap, 1e-9)

        .any_na <- !is.finite(.xi[-.nf]) | !is.finite(.yi[-.nf]) |
                   !is.finite(.xi[-1L])  | !is.finite(.yi[-1L])  |
                   !is.finite(.xj[-.nf]) | !is.finite(.yj[-.nf]) |
                   !is.finite(.xj[-1L])  | !is.finite(.yj[-1L])
        .flag <- .ratio_v > .sw_ratio & !.any_na & .prox_ok

        if (!any(.flag, na.rm = TRUE)) next

        if (.hung_avail) {
          for (.k in which(.flag)) {
            .pc <- cbind(.wx[.k + 1L, ], .wy[.k + 1L, ])
            .pp <- cbind(.wx[.k,       ], .wy[.k,       ])
            if (any(!is.finite(.pc)) || any(!is.finite(.pp))) next
            .cost <- matrix(rowSums(outer(seq_len(.nf_ids), seq_len(.nf_ids),
              FUN = function(a, b) (.pc[a, 1] - .pp[b, 1])^2 +
                                   (.pc[a, 2] - .pp[b, 2])^2
            )), .nf_ids, .nf_ids)
            .asgn <- as.integer(clue::solve_LSAP(.cost))
            if (.asgn[.i] == .j && .asgn[.j] == .i)
              .hung_conf_frame[.k] <- TRUE
            else
              .flag[.k] <- FALSE
          }
        }

        if (!any(.flag, na.rm = TRUE)) next

        .rle  <- rle(.flag)
        .ends <- cumsum(.rle$lengths)
        .stts <- c(1L, .ends[-length(.ends)] + 1L)

        .evs <- list()
        .cur <- NULL

        for (.r in seq_along(.rle$values)) {
          if (.rle$values[.r]) {
            if (is.null(.cur)) {
              .cur <- list(s = .stts[.r], e = .ends[.r],
                           ratios = .ratio_v[.stts[.r]:.ends[.r]])
            } else {
              .cur$e      <- .ends[.r]
              .cur$ratios <- c(.cur$ratios, .ratio_v[.stts[.r]:.ends[.r]])
            }
          } else {
            if (!is.null(.cur)) {
              .gap <- .ends[.r] - .stts[.r] + 1L
              if (.gap > .sw_gap) {
                if (.cur$e - .cur$s + 1L >= .sw_min_dur)
                  .evs[[length(.evs) + 1L]] <- .cur
                .cur <- NULL
              }
            }
          }
        }
        if (!is.null(.cur) && .cur$e - .cur$s + 1L >= .sw_min_dur)
          .evs[[length(.evs) + 1L]] <- .cur

        for (.ev in .evs) {
          .sw_frame <- .frames[.ev$s + 1L]
          .tier     <- if (any(.hung_conf_frame[.ev$s:.ev$e])) 3L else 1L
          .trial_events[[length(.trial_events) + 1L]] <- data.frame(
            trial_id               = .tid,
            fish_id_A              = .fid_a,
            fish_id_B              = .fid_b,
            switch_frame           = .sw_frame,
            run_length_frames      = .ev$e - .ev$s + 1L,
            mean_improvement_ratio = round(mean(.ev$ratios, na.rm = TRUE), 3),
            tier                   = .tier,
            confirmed_hungarian    = (.tier == 3L),
            stringsAsFactors       = FALSE
          )
        }
      }   # end pairs loop

      .trial_events
    },
    .options = furrr::furrr_options(
      seed     = TRUE,
      globals  = c(".sw_ratio", ".sw_min_dur", ".sw_gap",
                   ".sw_prox", ".hung_avail"),
      packages = if (.hung_avail) "clue" else character(0)
    )
  )

  .raw_events <- unlist(.sw_par_results, recursive = FALSE)

  # ---- Build events table ------------------------------------------------
  .sw_ev <- if (length(.raw_events) > 0L)
    do.call(rbind, .raw_events)
  else
    data.frame(
      trial_id = integer(), fish_id_A = character(), fish_id_B = character(),
      switch_frame = integer(), run_length_frames = integer(),
      mean_improvement_ratio = numeric(), tier = integer(),
      confirmed_hungarian = logical(), stringsAsFactors = FALSE
    )

  # Add switch_frame_end and duration per (trial, pair) group
  if (nrow(.sw_ev) > 0L) {
    .sw_ev <- .sw_ev[order(.sw_ev$trial_id,
                           .sw_ev$fish_id_A,
                           .sw_ev$fish_id_B,
                           .sw_ev$switch_frame), ]
    rownames(.sw_ev) <- NULL

    .sw_ev$switch_frame_end <- NA_integer_
    .pair_key <- paste(.sw_ev$trial_id, .sw_ev$fish_id_A, .sw_ev$fish_id_B, sep = "__")

    for (.pk in unique(.pair_key)) {
      .pk_i  <- which(.pair_key == .pk)
      .pk_fr <- .sw_ev$switch_frame[.pk_i]
      .tid_p <- .sw_ev$trial_id[.pk_i[1]]
      .max_f <- max(mfbf[[.fr_col]][mfbf$trial_id == .tid_p], na.rm = TRUE)

      for (.e in seq_along(.pk_i)) {
        .sw_ev$switch_frame_end[.pk_i[.e]] <-
          if (.e < length(.pk_i)) .pk_fr[.e + 1L] - 1L else .max_f
      }
    }
    .sw_ev$duration_frames <- .sw_ev$switch_frame_end - .sw_ev$switch_frame + 1L
  }

  # ---- CORRECTION: swap coordinates per event, chronologically ----------
  if (nrow(.sw_ev) > 0L) {
    # Process in trial × time order so later events override earlier ones
    .sw_ev_ord <- .sw_ev[order(.sw_ev$trial_id, .sw_ev$switch_frame), ]

    for (.e in seq_len(nrow(.sw_ev_ord))) {
      .tid_c   <- .sw_ev_ord$trial_id[.e]
      .fid_ac  <- .sw_ev_ord$fish_id_A[.e]
      .fid_bc  <- .sw_ev_ord$fish_id_B[.e]
      .fr_st   <- .sw_ev_ord$switch_frame[.e]

      # Row indices for each fish from switch_frame onwards in this trial
      .ra <- which(mfbf$trial_id == .tid_c &
                   as.character(mfbf[[.id_col]]) == .fid_ac &
                   mfbf[[.fr_col]] >= .fr_st)
      .rb <- which(mfbf$trial_id == .tid_c &
                   as.character(mfbf[[.id_col]]) == .fid_bc &
                   mfbf[[.fr_col]] >= .fr_st)
      if (length(.ra) == 0L || length(.rb) == 0L) next

      # Align on matching frames (fish may have unequal frame coverage)
      .fa <- mfbf[[.fr_col]][.ra];  .fb <- mfbf[[.fr_col]][.rb]
      .common <- intersect(.fa, .fb)
      if (length(.common) == 0L) next

      .ia <- .ra[order(.fa)][.fa[order(.fa)] %in% .common]
      .ib <- .rb[order(.fb)][.fb[order(.fb)] %in% .common]

      # Atomic vectorised swap for each column
      for (.col in .sw_swap_cols) {
        .tmp              <- mfbf[[.col]][.ia]
        mfbf[[.col]][.ia] <- mfbf[[.col]][.ib]
        mfbf[[.col]][.ib] <- .tmp
      }
    }

    # switch_flag: frames whose position actually changed vs pre-correction
    mfbf$switch_flag <- (
      (abs(mfbf[[.x_col]] - .x_pre_sw) > 1e-6 |
       abs(mfbf[[.y_col]] - .y_pre_sw) > 1e-6) &
      is.finite(.x_pre_sw) & is.finite(mfbf[[.x_col]])
    )

    .n_sw_frames <- sum(mfbf$switch_flag, na.rm = TRUE)
    ts_msg("  Events: ", nrow(.sw_ev),
           " | Frames corrected: ", .n_sw_frames,
           " (", round(100 * .n_sw_frames / nrow(mfbf), 2), "%)")

    if (.n_sw_frames > 0.10 * nrow(mfbf))
      warning(round(100 * .n_sw_frames / nrow(mfbf), 1),
              "% of frames corrected — unusually high. ",
              "Consider raising SWITCH_TIER1_RATIO or enabling SWITCH_USE_HUNGARIAN.",
              call. = FALSE)
  } else {
    mfbf$switch_flag <- FALSE
    ts_msg("  No identity switches detected.")
  }

  assign("identity_switch_events", .sw_ev, envir = .GlobalEnv)

} else {
  if (!"switch_flag" %in% names(mfbf)) mfbf$switch_flag <- FALSE
  assign("identity_switch_events",
    data.frame(
      trial_id = integer(), fish_id_A = character(), fish_id_B = character(),
      switch_frame = integer(), run_length_frames = integer(),
      mean_improvement_ratio = numeric(), tier = integer(),
      confirmed_hungarian = logical(), switch_frame_end = integer(),
      duration_frames = integer(), stringsAsFactors = FALSE
    ),
    envir = .GlobalEnv
  )
}

ts_msg("Switch detection complete.")


# =============================================================================
# ==== 3) PARSE FOLDER NAMES → ATTACH trial, timepoint AND TRIAL METADATA =====
# =============================================================================
# Session folders in PIPELINE_DATA_DIR use the naming convention:
#   N{trial}_segment{timepoint}[_optional_suffix]
# where:
#   {trial}     — integer 1–16; matches the 'trial' column in TRIAL_META
#   {timepoint} — integer 1–3; the 20-min segment within that trial
#
# Folders are sorted alphabetically and trial_ids assigned sequentially (1, 2,
# 3, …) — the same order used by .read_one_session() in section 2.  trial and
# timepoint are parsed directly from each folder name.
#
# This block runs BEFORE zone assignment so that motor_side (needed by
# assign_zones_choice_exp) is already in mfbf when that function is called.
# TRIAL_META has one row per trial; all segments of a trial inherit its metadata.
# =============================================================================

ts_msg("Parsing session folder names (N{trial}_segment{tp})...")

.data_dir_tp <- .get_global("PIPELINE_DATA_DIR")

if (!is.null(.data_dir_tp) && dir.exists(.data_dir_tp)) {

  .all_fld   <- sort(list.dirs(.data_dir_tp, recursive = FALSE, full.names = FALSE))
  .valid_fld <- .all_fld[
    grepl("N\\d+",       .all_fld, perl = TRUE) &
    grepl("seg(?:ment)?_?\\d+", .all_fld, ignore.case = TRUE, perl = TRUE)
  ]

  if (length(.valid_fld) == 0) {
    warning(
      "No folders matching 'N{n}...segment{tp}' found in PIPELINE_DATA_DIR.\n",
      "Expected names like 'N3_segment1' or 'N12_segment2'.\n",
      "trial and timepoint will be NA.",
      call. = FALSE
    )
    mfbf$trial     <- NA_integer_
    mfbf$timepoint <- NA_integer_

  } else {

    # Parse trial number (digits immediately following "N") and segment number.
    # Handles: N2_segment_1, N6_Seg2, N1_segment_1_1 (underscore between word and digit).
    .parsed_trial <- suppressWarnings(as.integer(
      regmatches(.valid_fld, regexpr("(?<=N)\\d+", .valid_fld, perl = TRUE))
    ))
    .parsed_tp <- suppressWarnings(as.integer(
      gsub(".*[Ss]eg(?:ment)?_?(\\d+).*", "\\1", .valid_fld, perl = TRUE)
    ))

    # Sequential trial_ids correspond to the alphabetical folder order
    .folder_map <- data.frame(
      trial_id  = seq_along(.valid_fld),
      trial     = .parsed_trial,
      timepoint = .parsed_tp,
      folder    = .valid_fld,
      stringsAsFactors = FALSE
    )

    .bad_parse <- .folder_map[is.na(.folder_map$trial) | is.na(.folder_map$timepoint), ]
    if (nrow(.bad_parse) > 0)
      warning(
        "Could not parse trial/timepoint from: ",
        paste(.bad_parse$folder, collapse = ", "),
        "\nCheck that names contain 'N{number}' and 'segment{number}'.",
        call. = FALSE
      )

    ts_msg("  ", length(.valid_fld), " session folder(s): trial range ",
           min(.parsed_trial, na.rm = TRUE), "-", max(.parsed_trial, na.rm = TRUE),
           "  segments: ", paste(sort(unique(.parsed_tp)), collapse = ", "))

    # Attach trial and timepoint to mfbf via sequential trial_id
    mfbf$trial_id <- suppressWarnings(as.integer(mfbf$trial_id))
    mfbf <- dplyr::left_join(
      mfbf,
      .folder_map[, c("trial_id", "trial", "timepoint")],
      by = "trial_id"
    )

    .unmatched_fld <- sum(is.na(mfbf$trial))
    if (.unmatched_fld > 0)
      warning(
        .unmatched_fld, " rows could not be matched to a session folder.\n",
        "Verify that the number of unique trial_ids in master_fish_by_frame equals\n",
        "the number of N{n}_segment{tp} folders in PIPELINE_DATA_DIR.",
        call. = FALSE
      )

    ts_msg("  trial attached:     ",
           paste(sort(unique(mfbf$trial[!is.na(mfbf$trial)])), collapse = ", "))
    ts_msg("  timepoint attached: ",
           paste(sort(unique(mfbf$timepoint[!is.na(mfbf$timepoint)])), collapse = ", "))
  }

} else {
  warning("PIPELINE_DATA_DIR not available; trial and timepoint will be NA.", call. = FALSE)
  mfbf$trial     <- NA_integer_
  mfbf$timepoint <- NA_integer_
}

mfbf$trial     <- as.integer(mfbf$trial)
mfbf$timepoint <- as.integer(mfbf$timepoint)


# =============================================================================
# ==== 3b) JOIN TRIAL METADATA ================================================
# =============================================================================
# mfbf$trial now holds the parsed trial number matching TRIAL_META$trial.
# Join treatment, fish_density, tank, motor_side, trial_date on that key.
# motor_side is joined here so that assign_zones_choice_exp() (section 4) can
# use it directly from mfbf without needing to look it up in TRIAL_META.

ts_msg("Joining trial metadata (treatment, fish_density, tank, motor_side, trial_date)...")

.meta_want <- c("trial", "treatment", "fish_density", "tank", "motor_side", "trial_date")
.meta_want <- c("trial", setdiff(.meta_want[-1], names(mfbf)))
.meta_slim <- TRIAL_META[, intersect(.meta_want, names(TRIAL_META))]
.meta_slim <- dplyr::distinct(.meta_slim, trial, .keep_all = TRUE)
.meta_slim$trial <- suppressWarnings(as.integer(.meta_slim$trial))

mfbf <- dplyr::left_join(mfbf, .meta_slim, by = "trial")

TREATMENT_LEVELS_g  <- .get_global("TREATMENT_LEVELS",  c("control", "exercise choice"))
DENSITY_LEVELS_g    <- .get_global("DENSITY_LEVELS",    c(4L, 8L, 12L, 16L))
DENSITY_AS_FACTOR_g <- .get_global("DENSITY_AS_FACTOR", TRUE)

mfbf$treatment <- factor(
  trimws(tolower(as.character(mfbf$treatment))),
  levels = TREATMENT_LEVELS_g
)

if (isTRUE(DENSITY_AS_FACTOR_g)) {
  mfbf$fish_density_f <- factor(mfbf$fish_density, levels = DENSITY_LEVELS_g, ordered = TRUE)
} else {
  mfbf$fish_density_f <- as.numeric(mfbf$fish_density)
}

.unmatched <- sum(is.na(mfbf$treatment))
if (.unmatched > 0)
  warning(
    .unmatched, " rows have NA treatment after metadata join.\n",
    "Check that all parsed trial numbers appear in the 'trial' column of TRIAL_META.",
    call. = FALSE
  )

ts_msg("Metadata join complete. Rows with treatment: ",
       sum(!is.na(mfbf$treatment)), " / ", nrow(mfbf))


# =============================================================================
# ==== 4) ASSIGN MULTI-LEVEL ZONES ============================================
# =============================================================================
# Assigns main_zone (flow / calm) and sec_zone (high / medium / low / calm)
# using the polygon PIP function defined in the master pipeline.
# motor_side is now in mfbf (attached in 3b), so the function's internal
# TRIAL_META lookup is skipped by its guard: if (!"motor_side" %in% names(dt)).

if (!exists("assign_zones_choice_exp", envir = .GlobalEnv, inherits = FALSE)) {
  stop(
    "assign_zones_choice_exp() not found in GlobalEnv.\n",
    "Run via 00_master_pipeline_choice_exp.R which defines this function.",
    call. = FALSE
  )
}

ts_msg("Assigning multi-level zones (main_zone + sec_zone)...")

mfbf <- assign_zones_choice_exp(mfbf)

# Report zone coverage
.zone_na_frac <- mean(is.na(mfbf$main_zone))
if (.zone_na_frac > 0.05) {
  warning(
    round(.zone_na_frac * 100, 1), "% of rows have NA main_zone.\n",
    "Check that zone_reference_choice_exp_FT.csv and zone_reference_choice_exp_FD.csv\n",
    "cover the full frame area and that x_interp / y_interp values are within the\n",
    "expected pixel range for their respective arena orientations.",
    call. = FALSE
  )
}
ts_msg("Zone assignment done. NA main_zone: ",
       round(.zone_na_frac * 100, 1), "%")


# =============================================================================
# ==== 5) WRITE BACK TO GLOBALENV AND DISK ====================================
# =============================================================================

assign("master_fish_by_frame", mfbf, envir = .GlobalEnv)
ts_msg("master_fish_by_frame updated in GlobalEnv (",
       nrow(mfbf), " rows, ", ncol(mfbf), " cols)")

# ---- Write to timestamped run folder under STEP1_output/ --------------------
.step1_parent <- file.path(getwd(), "STEP1_output")
if (!dir.exists(.step1_parent))
  dir.create(.step1_parent, recursive = TRUE, showWarnings = FALSE)
.step1_run_dir <- file.path(.step1_parent,
                             paste0("STEP1_output_", format(Sys.time(), "%Y%m%d_%H%M%S")))
dir.create(.step1_run_dir, recursive = TRUE, showWarnings = FALSE)

tryCatch({
  readr::write_csv(mfbf, file.path(.step1_run_dir, "master_fish_by_frame.csv"))
  ts_msg("master_fish_by_frame.csv written to: ", .step1_run_dir)
}, error = function(e)
  warning("Could not write master_fish_by_frame.csv: ", conditionMessage(e), call. = FALSE))

# Write jump detection summary (empty data.frame if detection was disabled)
tryCatch({
  .jd_summary <- .get_global("jump_detection_summary", data.frame())
  if (nrow(.jd_summary) > 0) {
    readr::write_csv(.jd_summary,
                     file.path(.step1_run_dir, "jump_detection_summary.csv"))
    ts_msg("jump_detection_summary.csv written (",
           sum(.jd_summary$n_jumps > 0), " fish-tracks with at least one jump)")
  }
}, error = function(e)
  warning("Could not write jump_detection_summary.csv: ", conditionMessage(e), call. = FALSE))

# Write identity switch events (empty data.frame if detection was disabled)
tryCatch({
  .sw_ev_out <- .get_global("identity_switch_events", data.frame())
  if (nrow(.sw_ev_out) > 0) {
    # Join timepoint and n_frames_trial from mfbf (one timepoint per trial_id;
    # n_frames_trial = number of unique frame indices in that trial).
    .tp_lut    <- unique(mfbf[, c("trial_id", "timepoint")])
    .fr_col_sw <- if ("frame_idx" %in% names(mfbf)) "frame_idx" else
                  if ("frame"     %in% names(mfbf)) "frame"     else NA_character_
    if (!is.na(.fr_col_sw)) {
      .nf_lut <- tapply(mfbf[[.fr_col_sw]], mfbf$trial_id,
                        function(f) length(unique(f)))
      .nf_df  <- data.frame(trial_id       = names(.nf_lut),
                             n_frames_trial = as.integer(.nf_lut),
                             stringsAsFactors = FALSE)
    } else {
      .nf_df <- data.frame(trial_id = unique(mfbf$trial_id),
                           n_frames_trial = NA_integer_,
                           stringsAsFactors = FALSE)
    }

    .sw_ev_out <- merge(.sw_ev_out, .tp_lut,  by = "trial_id", all.x = TRUE)
    .sw_ev_out <- merge(.sw_ev_out, .nf_df,   by = "trial_id", all.x = TRUE)

    # canonical column order
    .sw_cols_ordered <- c(
      "trial_id", "timepoint",
      "fish_id_A", "fish_id_B",
      "switch_frame", "switch_frame_end", "duration_frames",
      "run_length_frames", "mean_improvement_ratio",
      "tier", "confirmed_hungarian",
      "n_frames_trial"
    )
    .sw_ev_out <- .sw_ev_out[, intersect(.sw_cols_ordered, names(.sw_ev_out)),
                               drop = FALSE]
    .sw_ev_out <- .sw_ev_out[order(.sw_ev_out$trial_id,
                                    .sw_ev_out$switch_frame), ]

    readr::write_csv(.sw_ev_out,
                     file.path(.step1_run_dir, "identity_switch_events.csv"))
    ts_msg("identity_switch_events.csv written (",
           nrow(.sw_ev_out), " switch events across ",
           length(unique(.sw_ev_out$trial_id)), " trials)")
  } else {
    ts_msg("No identity switch events detected; identity_switch_events.csv not written.")
  }
}, error = function(e)
  warning("Could not write identity_switch_events.csv: ", conditionMessage(e), call. = FALSE))

assign("STEP1_OUTPUT_DIR", .step1_run_dir, envir = .GlobalEnv)

# ---- Update processed-sessions manifest -------------------------------------
# Always update after a successful full run (whether SKIP_PROCESSED_SESSIONS is
# TRUE or FALSE) so that future incremental runs know which sessions are done.
tryCatch({
  .all_sess_now <- sort(
    list.dirs(.get_global("PIPELINE_DATA_DIR", ""), recursive = FALSE, full.names = FALSE)
  )
  .all_sess_now <- .all_sess_now[
    grepl("N\\d+", .all_sess_now, perl = TRUE) &
    grepl("seg(?:ment)?_?\\d+", .all_sess_now, ignore.case = TRUE, perl = TRUE)
  ]
  if (length(.all_sess_now) > 0) {
    dir.create(dirname(.manifest_path), recursive = TRUE, showWarnings = FALSE)
    # Merge with any entries already in the manifest (preserves incremental runs)
    .existing_manifest <- if (file.exists(.manifest_path))
      readLines(.manifest_path, warn = FALSE) else character(0)
    .updated_manifest  <- sort(unique(c(.existing_manifest, .all_sess_now)))
    .updated_manifest  <- .updated_manifest[nzchar(.updated_manifest)]
    writeLines(.updated_manifest, .manifest_path)
    ts_msg("Manifest updated (", length(.updated_manifest),
           " sessions): ", .manifest_path)
  }
}, error = function(e)
  warning("Could not update processed sessions manifest: ",
          conditionMessage(e), call. = FALSE))

ts_msg("Step 1 complete. Output: ", .step1_run_dir)
ts_msg("  Columns added: trial, timepoint, main_zone, sec_zone, zone_label, ",
       "treatment, fish_density, fish_density_f, tank, motor_side, trial_date")
ts_msg("  Jump filter columns: x_interp_raw, y_interp_raw, jump_flag ",
       "(x_interp / y_interp updated with jumps removed and linearly interpolated)")
ts_msg("  Switch filter columns: switch_flag ",
       "(positions corrected in-place; identity_switch_events.csv details each event)")
