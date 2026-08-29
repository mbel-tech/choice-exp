# =============================================================================
# STEP 2b — CHOICE EXPERIMENT  |  GROUP DYNAMICS (TIER 1)
# Computes collective-behaviour metrics from master_fish_by_frame.
# Uses swaRm (CRAN) for core calculations; base-R fallback where noted.
#
# TIER 1 GROUP-LEVEL (one row per trial × timepoint → group_dynamics_summary):
#   mean_iid_px          Mean inter-individual distance (all pairs), pixels
#   sd_iid_px            SD of IID across frames
#   mean_nnd_px          Mean nearest-neighbour distance (group mean), pixels
#   sd_nnd_px            SD of NND across frames
#   mean_polarisation    Mean polarisation order parameter P ∈ [0,1]
#   sd_polarisation      SD of P across frames
#   mean_hull_area_px2   Mean convex-hull area of group, px²
#   sd_hull_area_px2     SD of hull area across frames
#   mean_centroid_spd_px Mean centroid speed (px/s)
#   sd_centroid_spd_px   SD of centroid speed (px/s)
#   n_frames_used        Frames with ≥ 3 fish present (used in calculations)
#
# TIER 1 INDIVIDUAL-LEVEL (one row per fish per trial × timepoint → merged
#   into fish_activity_summary and written as fish_social_context.csv):
#   mean_nnd_fish_px     Mean NND for this fish across frames
#   sd_nnd_fish_px       SD of this fish's NND
#   mean_centdist_px     Mean distance of this fish to group centroid
#   sd_centdist_px       SD of centroid distance
#   mean_turning_rate    Mean |Δheading| per frame (rad/frame)
#   sd_turning_rate      SD of turning rate
#   prop_central         Proportion of frames this fish was closest to centroid
#
# INPUT  (from GlobalEnv or latest STEP2_output/ CSV):
#   fish_activity_summary  — one row per fish per timepoint
#   master_fish_by_frame   — one row per fish per frame
#   Both must have: trial_id, fish_id, timepoint, x_interp, y_interp,
#                   heading_rad (optional; computed from trajectory if present)
#
# OUTPUTS (disk: STEP2b_output/STEP2b_output_<ts>/):
#   group_dynamics_summary.csv
#   fish_social_context.csv
#   fish_activity_summary_with_social.csv   (enriched; replaces GlobalEnv copy)
# =============================================================================


# =============================================================================
# ==== 0) PACKAGES + HELPERS ==================================================
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(furrr)
  library(future)
})

ts_msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), " | ", ..., "\n", sep = "")

.get_global <- function(name, default = NULL) {
  if (exists(name, envir = .GlobalEnv, inherits = FALSE))
    get(name, envir = .GlobalEnv, inherits = FALSE)
  else default
}

.step2b_dir <- tryCatch({
  frames <- sys.frames()
  oe     <- Filter(function(f) exists("ofile", envir = f, inherits = FALSE), frames)
  if (length(oe) > 0)
    normalizePath(dirname(get("ofile", envir = oe[[length(oe)]])),
                  winslash = "/", mustWork = FALSE)
  else getwd()
}, error = function(e) getwd())

.find_latest_csv <- function(step_name, csv_name) {
  # Search getwd() first, then parent of getwd() — handles pipeline run from subdirectory
  for (.base in c(getwd(), dirname(getwd()))) {
    parent <- file.path(.base, step_name)
    if (!dir.exists(parent)) next
    subs <- list.dirs(parent, full.names = TRUE, recursive = FALSE)
    subs <- subs[grepl(paste0("^", step_name, "_\\d{8}_\\d{6}$"), basename(subs))]
    if (!length(subs)) next
    f <- file.path(subs[which.max(file.mtime(subs))], csv_name)
    if (file.exists(f)) return(f)
  }
  NULL
}


# =============================================================================
# ==== 1) CHECK swaRm AVAILABILITY ============================================
# =============================================================================
# swaRm (CRAN) provides: nnd(), pol(), chull_area() used in section 3.
# swaRmverse (2025, CRAN) wraps swaRm for workflow-level analysis; optional.
# Install if missing:  install.packages("swaRm")

.swaRm_ok <- requireNamespace("swaRm", quietly = TRUE)
if (!.swaRm_ok) {
  warning(
    "swaRm package not found.  Install with:  install.packages('swaRm')\n",
    "Falling back to base-R implementations for NND, polarisation, hull area.\n",
    "Results are numerically identical but slightly slower.",
    call. = FALSE
  )
}
ts_msg("swaRm available: ", .swaRm_ok)


# =============================================================================
# ==== 2) LOAD DATA ===========================================================
# =============================================================================

# ---- master_fish_by_frame ---------------------------------------------------
if (!exists("master_fish_by_frame", envir = .GlobalEnv, inherits = FALSE)) {
  .mfbf_path <- .find_latest_csv("STEP2_output", "master_fish_by_frame_step2.csv")
  if (is.null(.mfbf_path))
    .mfbf_path <- .find_latest_csv("STEP1_output", "master_fish_by_frame.csv")
  if (is.null(.mfbf_path))
    stop("master_fish_by_frame not found in GlobalEnv or on disk. Run Steps 1-2 first.",
         call. = FALSE)
  master_fish_by_frame <- readr::read_csv(.mfbf_path, show_col_types = FALSE)
  assign("master_fish_by_frame", master_fish_by_frame, envir = .GlobalEnv)
  ts_msg("Loaded master_fish_by_frame: ", .mfbf_path)
}
mfbf <- get("master_fish_by_frame", envir = .GlobalEnv, inherits = FALSE)

# ---- fish_activity_summary --------------------------------------------------
if (!exists("fish_activity_summary", envir = .GlobalEnv, inherits = FALSE)) {
  # STEP2 writes trial_activity_summary.csv; fish_activity_summary.csv is alias
  .fas_path <- .find_latest_csv("STEP2_output", "fish_activity_summary.csv")
  if (is.null(.fas_path))
    .fas_path <- .find_latest_csv("STEP2_output", "trial_activity_summary.csv")
  if (is.null(.fas_path))
    stop("fish_activity_summary not found. Run Step 2 first.", call. = FALSE)
  fish_activity_summary <- readr::read_csv(.fas_path, show_col_types = FALSE)
  assign("fish_activity_summary", fish_activity_summary, envir = .GlobalEnv)
  ts_msg("Loaded fish_activity_summary: ", .fas_path)
}
fas <- get("fish_activity_summary", envir = .GlobalEnv, inherits = FALSE)

# ---- Validate required columns ----------------------------------------------
.need_mfbf <- c("trial_id", "fish_id", "frame", "x_interp", "y_interp", "timepoint")
.miss <- setdiff(.need_mfbf, names(mfbf))
if (length(.miss))
  stop("master_fish_by_frame missing columns: ", paste(.miss, collapse = ", "), call. = FALSE)

# Compute heading_rad on the fly if absent: angle of the displacement vector
# between consecutive frames per (trial_id, fish_id), using atan2(dy, dx).
# The first frame of each track gets NA (no prior position to compare against).
if (!"heading_rad" %in% names(mfbf)) {
  ts_msg("heading_rad not in master_fish_by_frame; computing from consecutive positions...")
  mfbf <- mfbf %>%
    dplyr::arrange(trial_id, fish_id, frame) %>%
    dplyr::group_by(trial_id, fish_id) %>%
    dplyr::mutate(
      .dx = x_interp - dplyr::lag(x_interp),
      .dy = y_interp - dplyr::lag(y_interp),
      heading_rad = atan2(.dy, .dx)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-.dx, -.dy)
}
.has_heading <- "heading_rad" %in% names(mfbf)

ts_msg("Data: ", nrow(mfbf), " frame-rows | ",
       dplyr::n_distinct(mfbf$trial_id), " trials | ",
       dplyr::n_distinct(mfbf$timepoint[!is.na(mfbf$timepoint)]), " timepoints")


# =============================================================================
# ==== 3) BASE-R FALLBACK HELPERS =============================================
# =============================================================================
# Used when swaRm not installed; numerically identical to swaRm internals.

# NND: vector of per-fish nearest-neighbour distances
.nnd_base <- function(x, y) {
  n <- length(x)
  if (n < 2L) return(rep(NA_real_, n))
  dm <- as.matrix(dist(cbind(x, y)))
  diag(dm) <- NA_real_
  apply(dm, 1, min, na.rm = TRUE)
}

# Polarisation: |mean unit heading vector| ∈ [0,1]
.pol_base <- function(headings) {
  h <- headings[is.finite(headings)]
  if (length(h) < 2L) return(NA_real_)
  abs(mean(complex(modulus = 1, argument = h)))
}

# Convex hull area (px²) via shoelace formula on chull vertices
.hull_area_base <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L) return(NA_real_)
  pts <- cbind(x[ok], y[ok])
  h   <- chull(pts)
  if (length(h) < 3L) return(NA_real_)
  px <- pts[h, 1]; py <- pts[h, 2]
  n  <- length(px)
  abs(sum(px * c(py[-1], py[1]) - c(px[-1], px[1]) * py)) / 2
}

# Mean inter-individual distance (all pairwise, upper triangle)
.iid_mean <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 2L) return(NA_real_)
  mean(as.vector(dist(cbind(x[ok], y[ok]))))
}


# =============================================================================
# ==== 4) PER-FRAME GROUP METRICS =============================================
# =============================================================================
# Compute metrics for every (trial_id × timepoint × frame) combination.
# Strategy: split mfbf by trial×timepoint×frame → apply per-frame functions.

ts_msg("Computing per-frame group metrics...")

mfbf <- mfbf %>%
  dplyr::mutate(
    trial_id  = as.integer(trial_id),
    fish_id   = as.character(fish_id),
    frame     = as.integer(frame),
    timepoint = as.integer(timepoint),
    x_interp  = as.numeric(x_interp),
    y_interp  = as.numeric(y_interp)
  )

if (.has_heading)
  mfbf <- dplyr::mutate(mfbf, heading_rad = as.numeric(heading_rad))

# ---- Split: one element = one (trial_id, timepoint, frame) ------------------
.frame_list <- mfbf %>%
  dplyr::filter(!is.na(timepoint)) %>%
  dplyr::group_by(trial_id, timepoint, frame) %>%
  dplyr::group_split()

ts_msg("  ", length(.frame_list), " unique (trial × timepoint × frame) groups")

.n_workers_2b <- .get_global("PIPELINE_WORKERS", 1L)
ts_msg("  Workers: ", .n_workers_2b)

# Ensure a plan is active (fallback when sourcing STEP2b standalone)
if (inherits(future::plan(), "SequentialFuture") && .n_workers_2b > 1L)
  future::plan(future::multisession, workers = .n_workers_2b)

# Collect function names from GlobalEnv that per-frame workers need
.step2b_globals <- c(".nnd_base", ".pol_base", ".hull_area_base", ".iid_mean",
                 ".swaRm_ok", ".has_heading")

# ---- Per-frame worker function -----------------------------------------------
.per_frame_fn <- function(g) {
  x   <- g$x_interp
  y   <- g$y_interp
  tid <- g$trial_id[1L]
  tp  <- g$timepoint[1L]
  fr  <- g$frame[1L]
  fid <- g$fish_id
  n_ok <- sum(is.finite(x) & is.finite(y))

  if (n_ok < 2L) return(NULL)

  nnd_v <- if (.swaRm_ok) {
    tryCatch(swaRm::nnd(x, y), error = function(e) .nnd_base(x, y))
  } else .nnd_base(x, y)

  cx <- mean(x, na.rm = TRUE); cy <- mean(y, na.rm = TRUE)
  cdist_v <- sqrt((x - cx)^2 + (y - cy)^2)
  iid  <- .iid_mean(x, y)
  pol  <- if (.has_heading) {
    h <- g$heading_rad
    if (.swaRm_ok)
      tryCatch(swaRm::pol(x, y, h), error = function(e) .pol_base(h))
    else .pol_base(h)
  } else NA_real_
  hull <- if (.swaRm_ok) {
    tryCatch(swaRm::chull_area(x, y), error = function(e) .hull_area_base(x, y))
  } else .hull_area_base(x, y)

  list(
    group_row = data.frame(
      trial_id  = tid, timepoint = tp, frame = fr,
      centroid_x = cx, centroid_y = cy,
      iid       = iid,
      mean_nnd  = mean(nnd_v, na.rm = TRUE),
      pol       = pol,
      hull_area = hull,
      n_fish    = n_ok,
      stringsAsFactors = FALSE
    ),
    fish_rows = data.frame(
      trial_id    = tid, timepoint = tp, frame = fr,
      fish_id     = fid,
      nnd_fish    = nnd_v,
      cent_dist   = cdist_v,
      stringsAsFactors = FALSE
    )
  )
}

# ---- Apply per-frame computations in chunks to cap peak RAM ------------------
# Collecting all 1.4M results at once requires 3× peak RAM (results list +
# two lapply extractions + two rbind outputs). Chunking aggregates each slice
# immediately and discards it before loading the next.
.furrr_opts <- furrr::furrr_options(
  seed     = TRUE,
  globals  = c(.step2b_globals, ".per_frame_fn"),
  packages = if (.swaRm_ok) c("swaRm") else character(0)
)

.chunk_size        <- 100000L
.n_total_frames    <- length(.frame_list)
.n_chunks          <- ceiling(.n_total_frames / .chunk_size)
.group_frames_list <- vector("list", .n_chunks)
.fish_frames_list  <- vector("list", .n_chunks)
.frames_done       <- 0L

for (.ci in seq_len(.n_chunks)) {
  .idx <- seq((.ci - 1L) * .chunk_size + 1L,
              min(.ci * .chunk_size, .n_total_frames))
  .chunk_res <- furrr::future_map(.frame_list[.idx], .per_frame_fn,
                                   .options = .furrr_opts)
  .chunk_res <- Filter(Negate(is.null), .chunk_res)
  .group_frames_list[[.ci]] <- do.call(rbind, lapply(.chunk_res, `[[`, "group_row"))
  .fish_frames_list[[.ci]]  <- do.call(rbind, lapply(.chunk_res, `[[`, "fish_rows"))
  .frames_done <- .frames_done + length(.chunk_res)
  rm(.chunk_res); invisible(gc(verbose = FALSE))
}

ts_msg("  Frames computed: ", .frames_done)

.group_frames <- do.call(rbind, .group_frames_list)
.fish_frames  <- do.call(rbind, .fish_frames_list)
rm(.group_frames_list, .fish_frames_list); invisible(gc(verbose = FALSE))


# =============================================================================
# ==== 5) CENTROID SPEED (lag-based, px/s) ====================================
# =============================================================================
# Speed of group centroid between consecutive frames within same trial×timepoint.
# Result is per-second: displacement (px) divided by elapsed time (s) between
# frames. dt is derived from the `time` column in master_fish_by_frame (median
# per trial, robust to dropped frames); falls back to 25 fps when `time` is
# absent.

if ("time" %in% names(mfbf)) {
  .dt_per_trial <- mfbf %>%
    dplyr::filter(!is.na(timepoint), is.finite(time)) %>%
    dplyr::distinct(trial_id, frame, time) %>%
    dplyr::arrange(trial_id, frame) %>%
    dplyr::group_by(trial_id) %>%
    dplyr::mutate(dt = time - dplyr::lag(time)) %>%
    dplyr::summarise(
      dt_s_trial = stats::median(dt[is.finite(dt) & dt > 0], na.rm = TRUE),
      .groups    = "drop"
    ) %>%
    dplyr::mutate(trial_id = as.integer(trial_id))
  .dt_fallback <- stats::median(.dt_per_trial$dt_s_trial, na.rm = TRUE)
  if (!is.finite(.dt_fallback) || .dt_fallback <= 0) .dt_fallback <- 0.04
  ts_msg("Frame interval (median across trials): ",
         round(.dt_fallback, 4), " s (~",
         round(1 / .dt_fallback, 2), " fps)")
} else {
  ts_msg("WARNING: `time` not in master_fish_by_frame; assuming 25 fps (dt = 0.04 s).")
  .dt_per_trial <- data.frame(trial_id = integer(0), dt_s_trial = numeric(0))
  .dt_fallback  <- 0.04
}

.group_frames <- .group_frames %>%
  dplyr::mutate(trial_id = as.integer(trial_id)) %>%
  dplyr::left_join(.dt_per_trial, by = "trial_id") %>%
  dplyr::mutate(dt_s_trial = dplyr::if_else(
    is.finite(dt_s_trial) & dt_s_trial > 0, dt_s_trial, .dt_fallback)) %>%
  dplyr::arrange(trial_id, timepoint, frame) %>%
  dplyr::group_by(trial_id, timepoint) %>%
  dplyr::mutate(
    centroid_spd = sqrt((centroid_x - dplyr::lag(centroid_x))^2 +
                        (centroid_y - dplyr::lag(centroid_y))^2) / dt_s_trial
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(-dt_s_trial)


# =============================================================================
# ==== 6) TURNING RATE — REMOVED (identity-dependent) =========================
# =============================================================================
# Deleted in the school-level refactor: per-fish turning rate requires a stable
# identity across consecutive frames within a session. In the choice pipeline
# identity is unreliable across long horizons, so this metric is dropped for
# all treatments to keep them comparable.


# =============================================================================
# ==== 7) SUMMARISE → GROUP_DYNAMICS_SUMMARY ==================================
# =============================================================================

ts_msg("Summarising group-level metrics per trial × timepoint...")

group_dynamics_summary <- .group_frames %>%
  dplyr::group_by(trial_id, timepoint) %>%
  dplyr::summarise(
    n_frames_used      = sum(n_fish >= 3L, na.rm = TRUE),
    mean_iid_px        = mean(iid[n_fish >= 3L],          na.rm = TRUE),
    sd_iid_px          = sd(iid[n_fish >= 3L],            na.rm = TRUE),
    mean_nnd_px        = mean(mean_nnd[n_fish >= 3L],     na.rm = TRUE),
    sd_nnd_px          = sd(mean_nnd[n_fish >= 3L],       na.rm = TRUE),
    mean_polarisation  = mean(pol[n_fish >= 3L],          na.rm = TRUE),
    sd_polarisation    = sd(pol[n_fish >= 3L],            na.rm = TRUE),
    mean_hull_area_px2 = mean(hull_area[n_fish >= 3L],   na.rm = TRUE),
    sd_hull_area_px2   = sd(hull_area[n_fish >= 3L],     na.rm = TRUE),
    mean_centroid_spd_px = mean(centroid_spd, na.rm = TRUE),
    sd_centroid_spd_px   = sd(centroid_spd,  na.rm = TRUE),
    .groups = "drop"
  )

# Join metadata from fish_activity_summary (one row per trial×timepoint is enough)
.meta_cols <- c("trial_id", "timepoint", "trial", "treatment", "fish_density",
                "fish_density_f", "tank", "motor_side", "trial_date", "timepoint_f")
.meta_slim <- fas %>%
  dplyr::select(dplyr::any_of(.meta_cols)) %>%
  dplyr::distinct(trial_id, timepoint, .keep_all = TRUE) %>%
  dplyr::mutate(trial_id = as.integer(trial_id), timepoint = as.integer(timepoint))

group_dynamics_summary <- group_dynamics_summary %>%
  dplyr::mutate(trial_id = as.integer(trial_id), timepoint = as.integer(timepoint)) %>%
  dplyr::left_join(.meta_slim, by = c("trial_id", "timepoint"))

ts_msg("group_dynamics_summary: ", nrow(group_dynamics_summary), " rows")


# =============================================================================
# ==== 7b) PIXEL → CM CONVERSION ==============================================
# =============================================================================
# Get per-trial length_unit (px/cm) from master_fish_by_frame and compute cm
# columns alongside the existing pixel columns.
# =============================================================================

.mfbf_g <- tryCatch(get("master_fish_by_frame", envir = .GlobalEnv,
                          inherits = FALSE),
                     error = function(e) NULL)
.lu_g <- if (!is.null(.mfbf_g) && "length_unit" %in% names(.mfbf_g)) {
  .lu_raw_g <- .mfbf_g %>%
    dplyr::filter(is.finite(length_unit), length_unit > 0) %>%
    dplyr::group_by(trial_id) %>%
    dplyr::summarise(length_unit_trial = stats::median(length_unit, na.rm = TRUE),
                     .groups = "drop") %>%
    dplyr::mutate(trial_id = as.integer(trial_id))
  if (nrow(.lu_raw_g) == 0L) {
    ts_msg("WARNING: no session has a valid length_unit; cm columns will be NA.")
    NULL
  } else {
    .lu_fallback_g  <- mean(.lu_raw_g$length_unit_trial, na.rm = TRUE)
    .all_trials_g   <- data.frame(trial_id = sort(unique(as.integer(
                         group_dynamics_summary$trial_id))))
    .lu_out_g <- dplyr::left_join(.all_trials_g, .lu_raw_g, by = "trial_id") %>%
      dplyr::mutate(length_unit_trial = dplyr::if_else(
        is.finite(length_unit_trial), length_unit_trial, .lu_fallback_g))
    .n_imp_g <- sum(!.all_trials_g$trial_id %in% .lu_raw_g$trial_id)
    if (.n_imp_g > 0)
      ts_msg("length_unit missing for ", .n_imp_g, " trial(s) in STEP2b; ",
             "imputed with cross-session mean (", round(.lu_fallback_g, 3), " px/cm)")
    .lu_out_g
  }
} else {
  ts_msg("WARNING: length_unit not in master_fish_by_frame; cm columns will be NA.")
  NULL
}

if (!is.null(.lu_g)) {
  group_dynamics_summary <- group_dynamics_summary %>%
    dplyr::left_join(.lu_g, by = "trial_id") %>%
    dplyr::mutate(
      mean_nnd_cm          = mean_nnd_px         / length_unit_trial,
      sd_nnd_cm            = sd_nnd_px            / length_unit_trial,
      mean_iid_cm          = mean_iid_px          / length_unit_trial,
      sd_iid_cm            = sd_iid_px            / length_unit_trial,
      mean_hull_area_cm2   = mean_hull_area_px2   / length_unit_trial^2,
      sd_hull_area_cm2     = sd_hull_area_px2     / length_unit_trial^2,
      mean_centroid_spd_cm = mean_centroid_spd_px / length_unit_trial,
      sd_centroid_spd_cm   = sd_centroid_spd_px   / length_unit_trial
    ) %>%
    dplyr::select(-length_unit_trial)
}


# =============================================================================
# ==== 8) FISH-LEVEL SOCIAL CONTEXT — REMOVED =================================
# =============================================================================
# fish_social_context (mean_nnd_fish_px, mean_centdist_px, mean_turning_rate,
# prop_central) required persistent fish identities. Dropped in the school-
# level refactor. The school-level analogues (mean_nnd_px, mean_iid_px,
# mean_polarisation, mean_hull_area_px2) live in group_dynamics_summary above
# and are identity-free.


# =============================================================================
# ==== 9) MERGE WITH trial_activity_summary ===================================
# =============================================================================

ts_msg("Merging group dynamics into trial_activity_summary...")

# fas is now the trial-level table (one row per trial × timepoint) emitted by
# STEP2. Join group_dynamics_summary on (trial_id, timepoint).
.group_cols <- c("n_frames_used",
                 "mean_iid_px", "sd_iid_px",
                 "mean_nnd_px", "sd_nnd_px",
                 "mean_iid_cm", "sd_iid_cm",
                 "mean_nnd_cm", "sd_nnd_cm",
                 "mean_polarisation", "sd_polarisation",
                 "mean_hull_area_px2", "sd_hull_area_px2",
                 "mean_hull_area_cm2", "sd_hull_area_cm2",
                 "mean_centroid_spd_px", "sd_centroid_spd_px",
                 "mean_centroid_spd_cm", "sd_centroid_spd_cm")

fas_enriched <- fas %>%
  dplyr::select(-dplyr::any_of(.group_cols)) %>%
  dplyr::mutate(trial_id = as.integer(trial_id),
                timepoint = as.integer(timepoint)) %>%
  dplyr::left_join(
    group_dynamics_summary %>%
      dplyr::select(trial_id, timepoint, dplyr::any_of(.group_cols)) %>%
      dplyr::mutate(trial_id = as.integer(trial_id),
                    timepoint = as.integer(timepoint)),
    by = c("trial_id", "timepoint")
  )

assign("trial_activity_summary",  fas_enriched,           envir = .GlobalEnv)
assign("fish_activity_summary",   fas_enriched,           envir = .GlobalEnv)  # alias
assign("group_dynamics_summary",  group_dynamics_summary, envir = .GlobalEnv)


# =============================================================================
# ==== 10) WRITE TO DISK ======================================================
# =============================================================================

.step2b_parent  <- file.path(getwd(), "STEP2b_output")
dir.create(.step2b_parent, recursive = TRUE, showWarnings = FALSE)
.step2b_run_dir <- file.path(.step2b_parent,
                              paste0("STEP2b_output_",
                                     format(Sys.time(), "%Y%m%d_%H%M%S")))
dir.create(.step2b_run_dir, recursive = TRUE, showWarnings = FALSE)

readr::write_csv(group_dynamics_summary,
                 file.path(.step2b_run_dir, "group_dynamics_summary.csv"))
ts_msg("group_dynamics_summary.csv written")

readr::write_csv(fas_enriched,
                 file.path(.step2b_run_dir, "trial_activity_summary_with_group.csv"))
ts_msg("trial_activity_summary_with_group.csv written")

assign("STEP2b_OUTPUT_DIR", .step2b_run_dir, envir = .GlobalEnv)

ts_msg("Step 2b complete. Output: ", .step2b_run_dir)
ts_msg("  Group metrics:  mean_nnd_cm, mean_iid_cm, mean_polarisation,",
       " mean_hull_area_cm2, mean_centroid_spd_cm  (px columns also retained)")
ts_msg("  trial_activity_summary updated in GlobalEnv with group dynamics columns")
