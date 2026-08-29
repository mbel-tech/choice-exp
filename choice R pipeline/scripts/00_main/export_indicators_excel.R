# Export all behavioural indicators per trial x timepoint to Excel
# Output: choice R pipeline/reports/indicators_by_trial_timepoint.xlsx

library(dplyr)
library(readr)
library(openxlsx)

# ── Paths ──────────────────────────────────────────────────────────────────────
STEP2B_DIR  <- file.path(PROJECT_ROOT, "choice R pipeline/output/STEP2b_output")
STEP5_CSV   <- file.path(PROJECT_ROOT, "choice R pipeline/STEP5_stats/STEP5_stats_20260511_175522/analysis_ready.csv")
OUT_FILE    <- file.path(PROJECT_ROOT, "choice R pipeline/reports/indicators_by_trial_timepoint.xlsx")

# ── Load data ──────────────────────────────────────────────────────────────────
# Latest STEP2b run (has all zone + group-dynamics indicators in one row per trial x tp)
step2b_latest <- sort(list.dirs(STEP2B_DIR, recursive = FALSE), decreasing = TRUE)[1]
message("Using STEP2b: ", basename(step2b_latest))
df_raw <- read_csv(file.path(step2b_latest, "trial_activity_summary_with_group.csv"),
                   show_col_types = FALSE)

# phys_trial_id from STEP5 (stable ID that collapses re-numbered idTracker trials)
df_step5 <- read_csv(STEP5_CSV, show_col_types = FALSE) |>
  select(trial_id, timepoint, phys_trial_id)

df <- left_join(df_raw, df_step5, by = c("trial_id", "timepoint"))

# ── Tidy & select columns ──────────────────────────────────────────────────────
df_out <- df |>
  transmute(
    # Identifiers
    phys_trial_id,
    trial_id,
    timepoint,
    timepoint_label = paste0("T", timepoint, " (", (timepoint - 1) * 20, "–", timepoint * 20, " min)"),
    treatment,
    tank,
    fish_density,
    motor_side,
    trial_date,

    # Zone occupancy — main zones
    pct_flow       = round(prop_time_in_flow  * 100, 2),
    pct_calm_main  = round(prop_time_in_calm  * 100, 2),

    # Zone occupancy — sub-zones (proportional to full session)
    pct_high       = round(prop_time_in_high   * 100, 2),
    pct_medium     = round(prop_time_in_medium * 100, 2),
    pct_low        = round(prop_time_in_low    * 100, 2),
    pct_calm_sec   = round(prop_time_in_calm_sec * 100, 2),

    # Activity
    prop_active,
    switches_per_session,
    n_main_switches,
    zone_flux_per_session,

    # Group dynamics
    mean_nnd_cm,
    mean_iid_cm,
    mean_polarisation,
    mean_hull_area_cm2,
    mean_centroid_spd_cm,

    # Sample size info
    n_frames,
    n_frames_used,
    obs_seconds
  ) |>
  arrange(treatment, phys_trial_id, timepoint)

# ── Build Excel workbook ───────────────────────────────────────────────────────
wb <- createWorkbook()
addWorksheet(wb, "All indicators")

# Column header style
hdr_style <- createStyle(
  fontName = "Arial", fontSize = 11, fontColour = "white",
  fgFill = "#2E4A7A", halign = "center", valign = "center",
  textDecoration = "bold", wrapText = TRUE,
  border = "Bottom", borderColour = "#1A2C50"
)

# Data styles
id_style <- createStyle(
  fontName = "Arial", fontSize = 10,
  fgFill = "#F2F2F2", halign = "left"
)
num_style <- createStyle(
  fontName = "Arial", fontSize = 10,
  numFmt = "0.00", halign = "right"
)
int_style <- createStyle(
  fontName = "Arial", fontSize = 10,
  numFmt = "0", halign = "right"
)
ec_style <- createStyle(
  fontName = "Arial", fontSize = 10,
  fgFill = "#E8F4FD", halign = "left"
)

# Write data
writeData(wb, "All indicators", df_out, startRow = 2, headerStyle = hdr_style)

# Row count (excludes header row we wrote at row 2)
n_rows  <- nrow(df_out)
n_cols  <- ncol(df_out)

# Freeze panes after columns A-C and row 2
freezePane(wb, "All indicators", firstActiveRow = 3, firstActiveCol = 4)

# Column widths
setColWidths(wb, "All indicators", cols = 1:n_cols,
             widths = c(12, 10, 10, 28, 17, 8, 12, 12, 16,
                        10, 14,
                        10, 10, 10, 12,
                        12, 20, 16, 22,
                        12, 12, 18, 20, 22,
                        10, 14, 12))

# Apply ID style to identifier columns (cols 1-9)
addStyle(wb, "All indicators", id_style,
         rows = 3:(n_rows + 2), cols = 1:9, gridExpand = TRUE)

# Apply numeric style to indicator columns
addStyle(wb, "All indicators", num_style,
         rows = 3:(n_rows + 2), cols = 10:n_cols, gridExpand = TRUE)

# Integer cols: timepoint, n_main_switches, switches, n_frames, n_frames_used, obs_sec
int_cols <- which(names(df_out) %in%
                    c("timepoint", "n_main_switches", "switches_per_session",
                      "zone_flux_per_session", "n_frames", "n_frames_used", "obs_seconds",
                      "fish_density"))
addStyle(wb, "All indicators", int_style,
         rows = 3:(n_rows + 2), cols = int_cols, gridExpand = TRUE)

# Tint exercise-choice rows
ec_rows <- which(df_out$treatment == "exercise choice") + 2
if (length(ec_rows) > 0) {
  addStyle(wb, "All indicators", ec_style,
           rows = ec_rows, cols = 1:9, gridExpand = TRUE)
}

# Row 1 = section label row
# Section header labels above the data columns
sections <- list(
  list(label = "IDENTIFIERS",           cols = 1:9,    fill = "#2E4A7A"),
  list(label = "MAIN ZONE OCCUPANCY %", cols = 10:11,  fill = "#1F6B42"),
  list(label = "SUB-ZONE OCCUPANCY %",  cols = 12:15,  fill = "#2D6A4F"),
  list(label = "ACTIVITY",              cols = 16:19,  fill = "#7B3F00"),
  list(label = "GROUP DYNAMICS",        cols = 20:24,  fill = "#4A235A"),
  list(label = "SAMPLE SIZE",           cols = 25:27,  fill = "#5D4037")
)

for (s in sections) {
  sec_style <- createStyle(
    fontName = "Arial", fontSize = 10, fontColour = "white",
    fgFill = s$fill, halign = "center", textDecoration = "bold",
    border = c("Top", "Bottom", "Left", "Right"), borderColour = "white"
  )
  # Merge cells across section columns in row 1
  if (length(s$cols) > 1) {
    mergeCells(wb, "All indicators",
               cols = s$cols, rows = 1)
  }
  writeData(wb, "All indicators", s$label,
            startCol = s$cols[1], startRow = 1)
  addStyle(wb, "All indicators", sec_style,
           rows = 1, cols = s$cols[1], stack = FALSE)
}

# Row height for section labels and header row
setRowHeights(wb, "All indicators", rows = 1, heights = 22)
setRowHeights(wb, "All indicators", rows = 2, heights = 40)

# Add a filter on the header row
addFilter(wb, "All indicators", row = 2, cols = 1:n_cols)

# ── Save ───────────────────────────────────────────────────────────────────────
saveWorkbook(wb, OUT_FILE, overwrite = TRUE)
message("Saved: ", OUT_FILE)
message("Rows: ", n_rows, "  |  Columns: ", n_cols)
