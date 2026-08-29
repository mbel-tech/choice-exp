# =============================================================================
# reconstruct_from_workbook.R
# -----------------------------------------------------------------------------
# Take the consolidated workbook (`behaviour_wide`, `cortisol`, `monoamines`)
# and rebuild data frames whose schemas match the easy_scripts CSVs:
#
#   - reconstruct_behaviour_from_workbook(wb_path)
#       returns a data frame with the SAME columns as
#       easy_scripts/easy_scripts_dataset.csv (where derivable).
#       Columns that were intentionally dropped from the workbook
#       (`prop_active`, `centroid_distance`, `turning_rate`,
#       `zone_flux_per_session`, `fish_fish_nnd`) are emitted as NA so the
#       schema matches and easy-script `read_csv` succeeds.
#
#   - reconstruct_endo_from_workbook(wb_path)
#       returns a long-format data frame with the SAME columns as
#       easy_scripts/easy_scripts_endo_dataset.csv:
#       sample_id, tank, trial, treatment, sex, plate, plate_date,
#       area, analyte, value, log_value, sqrt_value
#       Rebuilt by un-pivoting the wide `cortisol` and `monoamines` sheets.
#
# This file defines functions only — no side effects on source().
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readxl); library(readr)
})

# ---- behaviour: workbook → easy-scripts CSV shape ---------------------------
# Map from workbook column name -> easy-scripts column name.
# Mirrors .WB_BEH_RENAME_DEFAULT in compare_data_identity.R; redefined here
# so this file is independently loadable.
.WB_TO_ES_RENAME <- c(
  # 2026-08-18: same drift as compare_data_identity.R -- the workbook's
  # physical-trial key is `trial`, not `phys_trial`. Without this, beh_recon
  # had no phys_trial_id column and Layer 2 died on
  # "replacement has 0 rows, data has 48" (as.character(NULL) is character(0)).
  trial                   = "phys_trial_id",
  phys_trial              = "phys_trial_id",
  interval                = "timepoint_f",
  "exercise choicement"   = "treatment",
  mean_school_speed_cm_s  = "mean_centroid_spd_cm",
  mean_school_area_cm2    = "mean_hull_area_cm2",
  # Jacobs d-index columns
  jacobs_d_flow           = "D_main_flow",
  jacobs_d_high           = "D_sec_high",
  jacobs_d_medium         = "D_sec_medium",
  jacobs_d_low            = "D_sec_low",
  # Composite / CLR columns
  # workbook spells these without the _tp suffix; both kept as aliases
  main_cells_flow         = "p_main_flow_tp",
  sub_cells_high_clr      = "clr_sec_high_tp",
  main_cells_flow_tp      = "p_main_flow_tp",
  sub_cells_high_clr_tp   = "clr_sec_high_tp"
)
.WB_TO_ES_RENAME_ENDO <- c(
  "exercise choicement" = "treatment"
)

reconstruct_behaviour_from_workbook <- function(wb_path,
                                                es_template_path = NULL,
                                                sheet = "behaviour_wide") {
  bw <- as.data.frame(readxl::read_excel(wb_path, sheet = sheet))

  # Apply rename map (workbook -> easy-scripts column convention)
  for (wb_nm in names(.WB_TO_ES_RENAME)) {
    es_nm <- .WB_TO_ES_RENAME[[wb_nm]]
    if (wb_nm %in% names(bw)) names(bw)[names(bw) == wb_nm] <- es_nm
  }

  # If template path given, align column set: any column in the easy-scripts
  # CSV that we don't have, emit as NA. Any column in our reconstruction
  # absent from the template is preserved (it shouldn't happen but stay safe).
  if (!is.null(es_template_path) && file.exists(es_template_path)) {
    es_cols <- names(readr::read_csv(es_template_path, n_max = 0,
                                     show_col_types = FALSE))
    for (cn in setdiff(es_cols, names(bw))) bw[[cn]] <- NA
    # Reorder to match template
    bw <- bw[, c(es_cols, setdiff(names(bw), es_cols)), drop = FALSE]
  }

  bw
}

# ---- endocrine: workbook → easy-scripts long CSV shape ----------------------
reconstruct_endo_from_workbook <- function(wb_path,
                                            sheet_cort = "cortisol",
                                            sheet_mono = "monoamines") {
  cort <- as.data.frame(readxl::read_excel(wb_path, sheet = sheet_cort))
  mono <- as.data.frame(readxl::read_excel(wb_path, sheet = sheet_mono))

  # Apply endo rename map (e.g. "exercise choicement" -> "treatment")
  .ap <- function(df) {
    for (wb_nm in names(.WB_TO_ES_RENAME_ENDO)) {
      es_nm <- .WB_TO_ES_RENAME_ENDO[[wb_nm]]
      if (wb_nm %in% names(df)) names(df)[names(df) == wb_nm] <- es_nm
    }
    df
  }
  cort <- .ap(cort); mono <- .ap(mono)

  # plate_date in workbook is DD.MM.YYYY (or NA); convert back to ISO for
  # parity with easy-scripts CSV (which carries YYYY-MM-DD)
  .as_iso <- function(x) {
    if (is.null(x)) return(NA_character_)
    out <- suppressWarnings(as.Date(x, tryFormats = c("%d.%m.%Y", "%Y-%m-%d")))
    format(out, "%Y-%m-%d")
  }
  if ("plate_date" %in% names(cort)) cort$plate_date <- .as_iso(cort$plate_date)

  # --- cortisol rows: area = NA, analyte = "cort" ---------------------------
  cort_long <- data.frame(
    sample_id  = as.character(cort$sample_id),
    tank       = as.character(cort$tank),
    trial      = as.character(cort$trial),
    treatment  = as.character(cort$treatment),
    sex        = as.character(cort$sex),
    plate      = as.character(cort$plate),
    plate_date = as.character(cort$plate_date),
    area       = NA_character_,
    analyte    = "cort",
    value      = as.numeric(cort$cort_plasma),
    stringsAsFactors = FALSE
  )

  # Note: workbook may carry numeric-only sample_ids (e.g. "1","4") while the
  # easy-scripts CSV carries "B3_01","B3_04". reconstruct_endo preserves the
  # workbook form; normalisation happens in compare_monoamines_identity().

  # --- monoamines: un-pivot 28 analyte__area cells --------------------------
  mono_meta <- c("sample_id", "tank", "trial", "treatment", "sex")
  mono_meta_present <- intersect(mono_meta, names(mono))
  cell_cols <- setdiff(names(mono), mono_meta_present)
  cell_cols <- cell_cols[grepl("__", cell_cols, fixed = TRUE)]

  if (length(cell_cols) == 0) {
    mono_long <- data.frame(
      sample_id = character(), tank = character(), trial = character(),
      treatment = character(), sex = character(), plate = character(),
      plate_date = character(), area = character(), analyte = character(),
      value = numeric(), stringsAsFactors = FALSE
    )
  } else {
    mono_long <- mono %>%
      dplyr::select(dplyr::all_of(mono_meta_present),
                    dplyr::all_of(cell_cols)) %>%
      tidyr::pivot_longer(cols = dplyr::all_of(cell_cols),
                          names_to  = "cell",
                          values_to = "value") %>%
      tidyr::separate(cell, into = c("analyte", "area"), sep = "__",
                      remove = TRUE) %>%
      dplyr::mutate(
        sample_id  = as.character(sample_id),
        tank       = as.character(tank),
        trial      = as.character(trial),
        treatment  = as.character(treatment),
        sex        = as.character(sex),
        plate      = NA_character_,           # monoamines: no plate
        plate_date = NA_character_,
        value      = suppressWarnings(as.numeric(value))
      ) %>%
      dplyr::select(sample_id, tank, trial, treatment, sex, plate, plate_date,
                    area, analyte, value) %>%
      as.data.frame()
  }

  combined <- dplyr::bind_rows(cort_long, mono_long)

  # log_value and sqrt_value: derived columns
  combined$log_value  <- ifelse(is.finite(combined$value) & combined$value > 0,
                                log(combined$value),  NA_real_)
  combined$sqrt_value <- ifelse(is.finite(combined$value) & combined$value >= 0,
                                sqrt(combined$value), NA_real_)

  # Match column order of the easy-scripts CSV
  combined <- combined[, c("sample_id", "tank", "trial", "treatment", "sex",
                           "plate", "plate_date", "area", "analyte",
                           "value", "log_value", "sqrt_value"), drop = FALSE]
  combined
}
