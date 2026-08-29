# Temporary driver: pre-load TRIAL_META (bypassing the buggy ofile-based path
# resolution that only works when this script is source()d, not Rscript'd
# directly) then source() STEP2 and STEP2b in sequence so each writes its
# output under choice R pipeline/output/ as usual.
setwd(file.path(PROJECT_ROOT, "choice R pipeline/output"))

ts_msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), " | ", ..., "\n", sep = "")

suppressPackageStartupMessages({ library(readxl); library(dplyr) })
# Exact replica of the master pipeline's TRIAL_META block (00_master_pipeline_
# choice_exp.R:296-373) -- the real schema is video_ID/trial/tank/motor_side/
# fish_density/treatment/date (no 'trial_id' column; that only appears later,
# assigned during STEP1's zone/session parsing).
TREATMENT_LEVELS <- c("control", "exercise choice")
.tm <- readxl::read_excel(file.path(PROJECT_ROOT, "choice R pipeline/data/trial_summary_choice_exp.xlsx"))
if ("condition" %in% names(.tm) && !("treatment" %in% names(.tm)))
  .tm <- dplyr::rename(.tm, treatment = condition)
.tm$treatment <- trimws(tolower(as.character(.tm$treatment)))
TRIAL_META <- .tm %>% dplyr::mutate(
  trial        = suppressWarnings(as.integer(trial)),
  fish_density = suppressWarnings(as.integer(fish_density)),
  tank         = as.character(tank),
  motor_side   = trimws(as.character(motor_side)),
  treatment    = factor(treatment, levels = TREATMENT_LEVELS)
)
if ("video_ID" %in% names(TRIAL_META)) TRIAL_META$video_ID <- as.character(TRIAL_META$video_ID)
if (!"trial_date" %in% names(TRIAL_META)) {
  .dc <- intersect(c("date", "Date", "trial_date", "DATE"), names(TRIAL_META))[1]
  if (!is.na(.dc)) TRIAL_META <- dplyr::rename(TRIAL_META, trial_date = !!.dc)
}
assign("TRIAL_META", TRIAL_META, envir = .GlobalEnv)
cat("TRIAL_META pre-loaded:", nrow(TRIAL_META), "rows\n")

source(file.path(PROJECT_ROOT, "choice R pipeline/scripts/01_pipeline_analysis/activity_analysis_STEP2_choice_exp.R"),
       local = FALSE, echo = FALSE)
cat("\n=== STEP2 DONE, running STEP2b ===\n\n")

# Free the 7M-row master frame before STEP2b re-derives what it needs, mirroring
# the master pipeline's own RAM-management convention.
if (exists("master_fish_by_frame", envir = .GlobalEnv, inherits = FALSE))
  rm(master_fish_by_frame, envir = .GlobalEnv)
invisible(gc(verbose = FALSE, full = TRUE))

source(file.path(PROJECT_ROOT, "choice R pipeline/scripts/02_utilities/group_dynamics_STEP2b_choice_exp.R"),
       local = FALSE, echo = FALSE)
cat("\n=== STEP2b DONE ===\n")
