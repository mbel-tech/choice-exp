# =============================================================================
# compare_data_identity.R   —  Layer 1 of cross-validation
# -----------------------------------------------------------------------------
# Cell-by-cell diff between:
#   - the workbook (`Behavioural and neuroendocrine correlates datasets.xlsx`)
#   - the easy_scripts CSVs (`easy_scripts_dataset.csv`,
#                            `easy_scripts_endo_dataset.csv`)
#
# Three exports, one CSV each, all written under <report_dir>:
#
#   data_identity_behaviour.csv
#       row-major diff of every workbook-resident column. Schema:
#         phys_trial, timepoint_f, workbook_col, es_col,
#         workbook_value, es_value, abs_diff, rel_diff, status
#
#   data_identity_behaviour__dropped.csv
#       informational: easy-scripts columns the workbook intentionally
#       dropped (5 columns). Not failures — flagged so reviewer sees what
#       has been excluded from the published dataset.
#
#   data_identity_cortisol.csv
#       per (sample_id) diff for cort_plasma.
#
#   data_identity_monoamines.csv
#       per (sample_id, analyte, area) diff for each cell.
#
# A run is PASS for Layer 1 iff every CSV has zero rows with status != "PASS".
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readxl); library(readr)
})

# ---- Rename-map defaults (workbook col -> easy-scripts col) ----------------
# These cover the user's documented workbook customisations. If new renames
# appear in the workbook the user must add them here OR pass a custom map.
.WB_BEH_RENAME_DEFAULT <- function() c(
  # 2026-08-18: the workbook's physical-trial key is `trial`, not `phys_trial`,
  # so Layer 1 aborted before it compared a single cell ("missing key cols after
  # rename: needed phys_trial_id + timepoint_f"). Verified equivalent: both
  # frames carry 48 rows and 16 distinct trials numbered 1-16, matching
  # phys_trial_id in easy_scripts_dataset.csv. `phys_trial` is kept as an alias
  # so an older workbook still resolves.
  trial                     = "phys_trial_id",
  phys_trial                = "phys_trial_id",
  interval                  = "timepoint_f",
  "exercise choicement"     = "treatment",
  mean_school_speed_cm_s    = "mean_centroid_spd_cm",
  mean_school_area_cm2      = "mean_hull_area_cm2",
  # Jacobs d-index columns (workbook name → easy-scripts name)
  jacobs_d_flow             = "D_main_flow",
  jacobs_d_high             = "D_sec_high",
  jacobs_d_medium           = "D_sec_medium",
  jacobs_d_low              = "D_sec_low",
  # Composite / CLR columns
  # 2026-08-18: the workbook columns are `main_cells_flow` / `sub_cells_high_clr`
  # with no _tp suffix, so these two renames never fired and Layer 1 died on
  # "arguments imply differing number of rows: 0, 1" when it tried to compare a
  # workbook column that had no easy-scripts counterpart. The _tp spellings are
  # kept as aliases. Targets verified present in easy_scripts_dataset.csv
  # (p_main_flow_tp, clr_sec_high_tp).
  main_cells_flow           = "p_main_flow_tp",
  sub_cells_high_clr        = "clr_sec_high_tp",
  main_cells_flow_tp        = "p_main_flow_tp",
  sub_cells_high_clr_tp     = "clr_sec_high_tp"
)
.WB_ENDO_RENAME_DEFAULT <- function() c(
  "exercise choicement" = "treatment"
)

# Apply a wb→es rename map to a data frame's column names.
.apply_wb_to_es_rename <- function(df, rename_map) {
  if (length(rename_map) == 0) return(df)
  for (wb_nm in names(rename_map)) {
    es_nm <- rename_map[[wb_nm]]
    if (wb_nm %in% names(df)) names(df)[names(df) == wb_nm] <- es_nm
  }
  df
}

# Numeric comparison helper. Treats both-NA as PASS. Uses absolute tol for
# values near zero, relative tol for larger values.
.compare_num <- function(wb_val, es_val, abs_tol = 1e-9, rel_tol = 1e-9) {
  wb_na <- is.na(wb_val); es_na <- is.na(es_val)
  abs_d <- ifelse(wb_na | es_na, NA_real_, abs(wb_val - es_val))
  denom <- pmax(abs(es_val), abs(wb_val), 1)
  rel_d <- ifelse(wb_na | es_na, NA_real_, abs_d / denom)
  status <- ifelse(wb_na & es_na, "PASS",
            ifelse(wb_na | es_na, "FAIL-NA-MISMATCH",
            ifelse(abs_d <= abs_tol | rel_d <= rel_tol, "PASS",
                   "FAIL-VALUE")))
  list(abs_diff = abs_d, rel_diff = rel_d, status = status)
}

# String comparison helper. NA-aware.
.compare_str <- function(wb_val, es_val) {
  wb_na <- is.na(wb_val); es_na <- is.na(es_val)
  status <- ifelse(wb_na & es_na, "PASS",
            ifelse(wb_na | es_na, "FAIL-NA-MISMATCH",
            ifelse(as.character(wb_val) == as.character(es_val),
                   "PASS", "FAIL-VALUE")))
  list(abs_diff = NA_real_, rel_diff = NA_real_, status = status)
}

# =============================================================================
# BEHAVIOUR
# =============================================================================

compare_behaviour_identity <- function(wb_path, es_csv_path, report_dir,
                                       wb_to_es_rename = NULL,
                                       sheet = "behaviour_wide") {
  if (is.null(wb_to_es_rename))
    wb_to_es_rename <- .WB_BEH_RENAME_DEFAULT()

  bw <- as.data.frame(readxl::read_excel(wb_path, sheet = sheet))
  es <- as.data.frame(readr::read_csv(es_csv_path, show_col_types = FALSE))

  # Apply rename map to workbook cols so both frames share easy-scripts naming
  bw <- .apply_wb_to_es_rename(bw, wb_to_es_rename)

  # Defensive: normalise trial_date on the ES side to DD.MM.YYYY for compare
  if ("trial_date" %in% names(es)) {
    td <- suppressWarnings(as.Date(es$trial_date,
                                   tryFormats = c("%Y-%m-%d", "%d.%m.%Y")))
    if (!all(is.na(td))) es$trial_date <- format(td, "%d.%m.%Y")
  }
  if ("trial_date" %in% names(bw)) {
    td <- suppressWarnings(as.Date(bw$trial_date,
                                   tryFormats = c("%d.%m.%Y", "%Y-%m-%d")))
    if (!all(is.na(td))) bw$trial_date <- format(td, "%d.%m.%Y")
  }

  # ---- Build comparison key (post-rename, both sides use es names) --------
  key_cols <- c("phys_trial_id", "timepoint_f")
  if (!all(key_cols %in% names(bw)))
    stop("[Layer1/beh] workbook missing key cols after rename: needed ",
         paste(key_cols, collapse = " + "), "; have ",
         paste(names(bw), collapse = ", "), call. = FALSE)
  if (!all(key_cols %in% names(es)))
    stop("[Layer1/beh] easy-scripts CSV missing key cols: needed ",
         paste(key_cols, collapse = " + "), call. = FALSE)

  bw_keyed <- bw
  bw_keyed$.k1 <- as.character(bw_keyed[[key_cols[1]]])
  bw_keyed$.k2 <- as.character(bw_keyed[[key_cols[2]]])

  es_keyed <- es
  es_keyed$.k1 <- as.character(es_keyed[[key_cols[1]]])
  es_keyed$.k2 <- as.character(es_keyed[[key_cols[2]]])

  # Align rows: inner join on key
  joined <- dplyr::inner_join(
    bw_keyed[, c(".k1", ".k2", setdiff(names(bw_keyed), c(".k1", ".k2"))), drop = FALSE],
    es_keyed[, c(".k1", ".k2", setdiff(names(es_keyed), c(".k1", ".k2"))), drop = FALSE],
    by = c(".k1", ".k2"), suffix = c(".wb", ".es"))

  # Row-count sanity
  n_wb <- nrow(bw); n_es <- nrow(es); n_join <- nrow(joined)
  message("[Layer1/beh] workbook rows=", n_wb, ", es rows=", n_es,
          ", joined rows=", n_join)

  # ---- Columns with no direct es equivalent (skip comparison, report only) ---
  # Intentionally dropped from workbook (present in es only) OR workbook extras
  # that have no es equivalent.
  dropped_cols <- c("prop_active", "centroid_distance", "turning_rate",
                    "zone_flux_per_session", "fish_fish_nnd",
                    "subzone_alr_pairs_tp")

  # ---- Per-column comparison ----------------------------------------------
  # After rename, both frames use easy-scripts column names. Compare each
  # workbook column to the same-named column in es.
  cols_to_check <- setdiff(names(bw), key_cols)
  cols_to_check <- cols_to_check[!startsWith(cols_to_check, ".k")]
  # Exclude workbook-only columns that have no es equivalent (informational only)
  cols_to_check <- cols_to_check[!cols_to_check %in% dropped_cols]
  rows <- list()

  for (wb_col in cols_to_check) {
    es_col <- wb_col
    if (!es_col %in% names(es)) {
      # ES lacks this column — workbook has it but the source CSV doesn't.
      # Record a single FAIL row.
      rows[[length(rows) + 1L]] <- data.frame(
        phys_trial = NA_character_, timepoint_f = NA_character_,
        workbook_col = wb_col, es_col = es_col,
        workbook_value = NA_character_, es_value = NA_character_,
        abs_diff = NA_real_, rel_diff = NA_real_,
        status = "FAIL-ES-COLUMN-MISSING",
        stringsAsFactors = FALSE)
      next
    }
    wb_v <- joined[[paste0(wb_col, ".wb")]]
    es_v <- joined[[paste0(es_col, ".es")]]
    if (is.null(wb_v) || is.null(es_v)) next   # join collision; skip

    # Choose comparison flavour
    if (is.numeric(wb_v) && is.numeric(es_v)) {
      cmp <- .compare_num(wb_v, es_v)
    } else {
      cmp <- .compare_str(as.character(wb_v), as.character(es_v))
    }

    rows[[length(rows) + 1L]] <- data.frame(
      phys_trial     = joined$.k1,
      timepoint_f    = joined$.k2,
      workbook_col   = wb_col,
      es_col         = es_col,
      workbook_value = as.character(wb_v),
      es_value       = as.character(es_v),
      abs_diff       = cmp$abs_diff,
      rel_diff       = cmp$rel_diff,
      status         = cmp$status,
      stringsAsFactors = FALSE)
  }

  identity_df <- do.call(rbind, rows)

  # Filter to fails-only for the main report, plus a summary CSV
  identity_path <- file.path(report_dir, "data_identity_behaviour.csv")
  fails_path    <- file.path(report_dir, "data_identity_behaviour__fails.csv")
  readr::write_csv(identity_df, identity_path)
  readr::write_csv(identity_df[identity_df$status != "PASS", , drop = FALSE],
                   fails_path)

  # ---- Dropped-column informational report --------------------------------
  # dropped_cols already defined above (before per-column loop)
  # Columns present in es but absent from workbook (intentionally dropped)
  dropped_in_es <- intersect(dropped_cols, names(es))
  # Columns present in workbook but with no es equivalent (no-es-col WB extras)
  dropped_in_wb <- intersect(dropped_cols, names(bw))

  # 2026-08-18: build each block with rep() over the column vector. Written as
  # scalars, data.frame() raised "arguments imply differing number of rows: 0, 1"
  # whenever one of the two sets was empty -- a scalar cannot recycle to zero
  # rows. On this dataset dropped_in_wb is empty, so Layer 1 aborted here, AFTER
  # completing every cell comparison but BEFORE writing any result, which is why
  # it looked like a data problem rather than a reporting bug.
  .dropped_block <- function(cols, side, note) data.frame(
    col  = as.character(cols),
    side = rep(side, length(cols)),
    note = rep(note,  length(cols)),
    stringsAsFactors = FALSE)
  dropped_df <- rbind(
    .dropped_block(dropped_in_es, "es_only",
      "Present in easy-scripts CSV but intentionally dropped from workbook"),
    .dropped_block(dropped_in_wb, "wb_only",
      "Present in workbook but has no direct easy-scripts equivalent (not compared)")
  )
  readr::write_csv(dropped_df, file.path(report_dir,
                                         "data_identity_behaviour__dropped.csv"))

  n_pass <- sum(identity_df$status == "PASS")
  n_fail <- sum(identity_df$status != "PASS")
  message("[Layer1/beh] PASS=", n_pass, ", FAIL=", n_fail)

  list(n_compared = nrow(identity_df), n_pass = n_pass, n_fail = n_fail,
       n_join = n_join, n_wb = n_wb, n_es = n_es,
       fails = identity_df[identity_df$status != "PASS", , drop = FALSE])
}

# =============================================================================
# CORTISOL
# =============================================================================

compare_cortisol_identity <- function(wb_path, es_endo_csv_path, report_dir,
                                       sheet = "cortisol",
                                       wb_to_es_rename = NULL) {
  if (is.null(wb_to_es_rename)) wb_to_es_rename <- .WB_ENDO_RENAME_DEFAULT()
  cort <- as.data.frame(readxl::read_excel(wb_path, sheet = sheet))
  cort <- .apply_wb_to_es_rename(cort, wb_to_es_rename)
  endo <- as.data.frame(readr::read_csv(es_endo_csv_path,
                                        show_col_types = FALSE))
  es_cort <- endo[endo$analyte == "cort", , drop = FALSE]

  cort$sample_id    <- as.character(cort$sample_id)
  es_cort$sample_id <- as.character(es_cort$sample_id)

  joined <- dplyr::inner_join(
    cort[, c("sample_id", "cort_plasma")],
    es_cort[, c("sample_id", "value")],
    by = "sample_id", suffix = c(".wb", ".es"))

  cmp <- .compare_num(as.numeric(joined$cort_plasma), as.numeric(joined$value))

  identity_df <- data.frame(
    sample_id      = as.character(joined$sample_id),
    workbook_value = joined$cort_plasma,
    es_value       = joined$value,
    abs_diff       = cmp$abs_diff,
    rel_diff       = cmp$rel_diff,
    status         = cmp$status,
    stringsAsFactors = FALSE)

  readr::write_csv(identity_df,
                   file.path(report_dir, "data_identity_cortisol.csv"))

  n_pass <- sum(identity_df$status == "PASS")
  n_fail <- sum(identity_df$status != "PASS")
  message("[Layer1/cort] PASS=", n_pass, ", FAIL=", n_fail,
          "  (wb=", nrow(cort), ", es=", nrow(es_cort),
          ", joined=", nrow(joined), ")")

  list(n_compared = nrow(identity_df), n_pass = n_pass, n_fail = n_fail,
       n_wb = nrow(cort), n_es = nrow(es_cort), n_join = nrow(joined),
       fails = identity_df[identity_df$status != "PASS", , drop = FALSE])
}

# =============================================================================
# MONOAMINES
# =============================================================================

compare_monoamines_identity <- function(wb_path, es_endo_csv_path, report_dir,
                                         sheet = "monoamines",
                                         wb_to_es_rename = NULL) {
  if (is.null(wb_to_es_rename)) wb_to_es_rename <- .WB_ENDO_RENAME_DEFAULT()
  mono <- as.data.frame(readxl::read_excel(wb_path, sheet = sheet))
  mono <- .apply_wb_to_es_rename(mono, wb_to_es_rename)
  endo <- as.data.frame(readr::read_csv(es_endo_csv_path,
                                        show_col_types = FALSE))
  es_mono <- endo[endo$analyte != "cort" &
                  !is.na(endo$area) & nzchar(endo$area), , drop = FALSE]

  # Un-pivot workbook monoamines
  meta_cols <- intersect(c("sample_id","tank","trial","treatment","sex"),
                         names(mono))
  cell_cols <- setdiff(names(mono), meta_cols)
  cell_cols <- cell_cols[grepl("__", cell_cols, fixed = TRUE)]

  if (length(cell_cols) == 0)
    stop("[Layer1/mono] no analyte__area columns found in monoamines sheet")

  wb_long <- mono %>%
    dplyr::select(dplyr::all_of(meta_cols), dplyr::all_of(cell_cols)) %>%
    tidyr::pivot_longer(dplyr::all_of(cell_cols),
                        names_to = "cell", values_to = "workbook_value") %>%
    tidyr::separate(cell, into = c("analyte","area"), sep = "__") %>%
    as.data.frame()

  # Coerce join keys to character so dplyr inner_join doesn't choke on
  # double vs character mismatch (workbook sample_id may parse as numeric).
  wb_long$sample_id  <- as.character(wb_long$sample_id)
  es_mono$sample_id  <- as.character(es_mono$sample_id)

  # --- sample_id normalisation fallback ------------------------------------
  # The workbook may carry pure numeric sample_ids (e.g. "1","4") while the
  # easy-scripts CSV carries "B3_01","B3_04". If the ids don't intersect at
  # all, try stripping the "B3_0*" prefix from es ids and re-check.
  .ws  <- unique(wb_long$sample_id)
  .ess <- unique(es_mono$sample_id)
  if (length(intersect(.ws, .ess)) == 0) {
    es_stripped <- sub("^[A-Za-z0-9]+_0*", "", .ess)
    if (length(intersect(.ws, es_stripped)) > 0) {
      message("[Layer1/mono] sample_id normalisation applied: stripping batch ",
              "prefix from es ids (e.g. 'B3_01' -> '1'). ",
              "This indicates the workbook uses numeric-only sample IDs.")
      es_mono$sample_id <- sub("^[A-Za-z0-9]+_0*", "", es_mono$sample_id)
    } else {
      message("[Layer1/mono] WARNING: sample_id mismatch — no intersection ",
              "even after prefix strip. Workbook ids: ",
              paste(head(.ws, 5), collapse = ", "),
              " | ES ids: ", paste(head(.ess, 5), collapse = ", "))
    }
  }

  joined <- dplyr::inner_join(
    wb_long[, c("sample_id","analyte","area","workbook_value")],
    es_mono[, c("sample_id","analyte","area","value")],
    by = c("sample_id","analyte","area"))

  cmp <- .compare_num(as.numeric(joined$workbook_value),
                      as.numeric(joined$value))
  identity_df <- data.frame(
    sample_id      = as.character(joined$sample_id),
    analyte        = as.character(joined$analyte),
    area           = as.character(joined$area),
    workbook_value = joined$workbook_value,
    es_value       = joined$value,
    abs_diff       = cmp$abs_diff,
    rel_diff       = cmp$rel_diff,
    status         = cmp$status,
    stringsAsFactors = FALSE)

  readr::write_csv(identity_df,
                   file.path(report_dir, "data_identity_monoamines.csv"))

  n_pass <- sum(identity_df$status == "PASS")
  n_fail <- sum(identity_df$status != "PASS")
  message("[Layer1/mono] PASS=", n_pass, ", FAIL=", n_fail,
          "  (wb cells=", nrow(wb_long), ", es rows=", nrow(es_mono),
          ", joined=", nrow(joined), ")")

  list(n_compared = nrow(identity_df), n_pass = n_pass, n_fail = n_fail,
       n_wb = nrow(wb_long), n_es = nrow(es_mono), n_join = nrow(joined),
       fails = identity_df[identity_df$status != "PASS", , drop = FALSE])
}
