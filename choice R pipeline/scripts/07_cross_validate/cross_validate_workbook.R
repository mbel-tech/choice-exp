# =============================================================================
# cross_validate_workbook.R   —  ENTRY POINT
# -----------------------------------------------------------------------------
# Cross-validates `Behavioural and neuroendocrine correlates datasets.xlsx`
# against the easy_scripts CSVs that the pipelines wrote.
#
# Two layers of evidence:
#   Layer 1  data identity check (compare_data_identity.R)
#            Every workbook cell must equal the matching easy-scripts cell.
#   Layer 2  refit equivalence check (refit_equivalence.R)
#            For a 50% random sample of indicators per dataset, refit the
#            model using the workbook-derived data and confirm identical
#            ANOVA statistics.
#
# Outputs land in:
#   D:/CHOICE R SCRIPTS/choice R pipeline/validation_reports/<timestamp>/
#
# Run standalone (after Stages 1–4 of `_run_full_session.R`):
#   source(file.path(PROJECT_ROOT, "choice R pipeline/scripts/07_cross_validate/cross_validate_workbook.R"))
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(readxl)
})

.CV_ROOT     <- file.path(PROJECT_ROOT, "choice R pipeline")
.CV_WB_PATH  <- file.path(.CV_ROOT,
                          "Behavioural and neuroendocrine correlates datasets.xlsx")
.CV_ES_BEH   <- file.path(.CV_ROOT, "easy_scripts/easy_scripts_dataset.csv")
.CV_ES_ENDO  <- file.path(.CV_ROOT, "easy_scripts/easy_scripts_endo_dataset.csv")
.CV_OUT_ROOT <- file.path(.CV_ROOT, "validation_reports")
.CV_SCRIPT_DIR <- file.path(.CV_ROOT, "scripts/07_cross_validate")

# ---- Pre-flight ------------------------------------------------------------
.cv_log <- function(...) cat("[cross-validate] ", ..., "\n", sep = "")

.cv_log("Workbook:        ", .CV_WB_PATH)
.cv_log("Easy-scripts beh: ", .CV_ES_BEH)
.cv_log("Easy-scripts endo: ", .CV_ES_ENDO)

stopifnot(file.exists(.CV_WB_PATH))
stopifnot(file.exists(.CV_ES_BEH))
stopifnot(file.exists(.CV_ES_ENDO))

# ---- Auto-resolve sheet names ---------------------------------------------
# The user may rename sheets between runs (e.g. behaviour_wide ->
# Preference_exp_behavior). Resolve by pattern, case-insensitive.
.resolve_sheet <- function(shts, patterns, role) {
  for (pat in patterns) {
    hit <- grep(pat, shts, ignore.case = TRUE, value = TRUE)
    # Exclude legend sheets
    hit <- hit[!grepl("legend", hit, ignore.case = TRUE)]
    # Exclude observation-study sheets (we only validate the choice exp)
    hit <- hit[!grepl("^obs[_]?study", hit, ignore.case = TRUE)]
    if (length(hit) >= 1) return(hit[1])
  }
  stop("[cross-validate] could not resolve sheet for role '", role,
       "' among: ", paste(shts, collapse = ", "), call. = FALSE)
}

.cv_all_sheets <- readxl::excel_sheets(.CV_WB_PATH)
.CV_SHEET_BEH  <- .resolve_sheet(.cv_all_sheets,
                                 c("^behaviour_wide$",
                                   "^preference[_ ]?exp[_ ]?behav",
                                   "behav"),
                                 "behaviour")
.CV_SHEET_CORT <- .resolve_sheet(.cv_all_sheets,
                                 c("^cortisol$",
                                   "^preference[_ ]?exp[_ ]?cort",
                                   "cort"),
                                 "cortisol")
.CV_SHEET_MONO <- .resolve_sheet(.cv_all_sheets,
                                 c("^monoamines$",
                                   "^preference[_ ]?exp[_ ]?mono",
                                   "mono"),
                                 "monoamines")
.cv_log("Resolved sheets:  beh=", .CV_SHEET_BEH,
        "  cort=", .CV_SHEET_CORT,
        "  mono=", .CV_SHEET_MONO)

# ---- Source helpers --------------------------------------------------------
source(file.path(.CV_SCRIPT_DIR, "reconstruct_from_workbook.R"))
source(file.path(.CV_SCRIPT_DIR, "compare_data_identity.R"))
source(file.path(.CV_SCRIPT_DIR, "refit_equivalence.R"))

# ---- Create timestamped output dir -----------------------------------------
.cv_ts        <- format(Sys.time(), "%Y%m%d_%H%M%S")
.cv_report    <- file.path(.CV_OUT_ROOT, .cv_ts)
dir.create(.cv_report, showWarnings = FALSE, recursive = TRUE)
.cv_log("Report dir:      ", .cv_report)

# ============================================================================
# LAYER 1 — Data identity
# ============================================================================
.cv_log("=== LAYER 1: data identity ===")
.cv_layer1_t0 <- Sys.time()

L1_beh <- compare_behaviour_identity(.CV_WB_PATH, .CV_ES_BEH,  .cv_report,
                                     sheet = .CV_SHEET_BEH)
L1_cor <- compare_cortisol_identity( .CV_WB_PATH, .CV_ES_ENDO, .cv_report,
                                     sheet = .CV_SHEET_CORT)
L1_mon <- compare_monoamines_identity(.CV_WB_PATH, .CV_ES_ENDO, .cv_report,
                                      sheet = .CV_SHEET_MONO)

.cv_layer1_secs <- round(as.numeric(difftime(Sys.time(), .cv_layer1_t0,
                                             units = "secs")), 1)

.layer1_pass <- (L1_beh$n_fail + L1_cor$n_fail + L1_mon$n_fail) == 0

# ============================================================================
# LAYER 2 — Refit equivalence
# ============================================================================
.cv_log("=== LAYER 2: refit equivalence (this takes a few minutes) ===")
.cv_layer2_t0 <- Sys.time()

L2 <- tryCatch(
  run_refit_equivalence(
    wb_path      = .CV_WB_PATH,
    es_beh_csv   = .CV_ES_BEH,
    es_endo_csv  = .CV_ES_ENDO,
    report_dir   = .cv_report,
    reconstruct_behaviour_fn = function(wb_path, es_template_path = NULL)
      reconstruct_behaviour_from_workbook(wb_path, es_template_path,
                                          sheet = .CV_SHEET_BEH),
    reconstruct_endo_fn      = function(wb_path)
      reconstruct_endo_from_workbook(wb_path,
                                     sheet_cort = .CV_SHEET_CORT,
                                     sheet_mono = .CV_SHEET_MONO),
    seed = 20260528),
  error = function(e) {
    .cv_log("[error] Layer 2 crashed: ", conditionMessage(e))
    list(n_compared = 0, n_pass = 0, n_fail = 0, n_skip = 0,
         sampling = NULL, fails = data.frame(), tmp_dir = NA_character_)
  })

.cv_layer2_secs <- round(as.numeric(difftime(Sys.time(), .cv_layer2_t0,
                                             units = "secs")), 1)

.layer2_pass <- (L2$n_fail == 0) && (L2$n_pass > 0)

# ============================================================================
# SUMMARY
# ============================================================================
.banner_overall <- if (.layer1_pass && .layer2_pass) "PASS" else "FAIL"

.summary_lines <- c(
  "============================================================================",
  paste0(" Cross-validation report — ", .cv_ts),
  paste0(" Workbook: ", basename(.CV_WB_PATH)),
  "============================================================================",
  "",
  paste0(" Overall verdict: ", .banner_overall),
  "",
  " Layer 1 — Data identity",
  paste0("   behaviour:   PASS=", L1_beh$n_pass, "  FAIL=", L1_beh$n_fail,
         "  (joined rows=", L1_beh$n_join, ")"),
  paste0("   cortisol:    PASS=", L1_cor$n_pass, "  FAIL=", L1_cor$n_fail,
         "  (joined rows=", L1_cor$n_join, ")"),
  paste0("   monoamines:  PASS=", L1_mon$n_pass, "  FAIL=", L1_mon$n_fail,
         "  (joined cells=", L1_mon$n_join, ")"),
  paste0("   elapsed: ", .cv_layer1_secs, " s"),
  "",
  " Layer 2 — Refit equivalence",
  paste0("   PASS=", L2$n_pass, "  FAIL=", L2$n_fail, "  SKIPPED=", L2$n_skip,
         "  (comparisons=", L2$n_compared, ")"),
  {
    if (!is.null(L2$sampling))
      c(paste0("   seed: ", L2$sampling$seed),
        paste0("   sampled behaviour (", length(L2$sampling$beh), "): ",
               paste(L2$sampling$beh, collapse = ", ")),
        paste0("   sampled cortisol  (1): ", L2$sampling$cort),
        paste0("   sampled monoamine (", length(L2$sampling$mono), "): ",
               paste(L2$sampling$mono, collapse = ", ")))
    else
      "   (Layer 2 produced no sampling)"
  },
  paste0("   elapsed: ", .cv_layer2_secs, " s"),
  "",
  " Drill-down CSVs:",
  paste0("   ", file.path(.cv_report, "data_identity_behaviour.csv")),
  paste0("   ", file.path(.cv_report, "data_identity_behaviour__fails.csv")),
  paste0("   ", file.path(.cv_report, "data_identity_behaviour__dropped.csv")),
  paste0("   ", file.path(.cv_report, "data_identity_cortisol.csv")),
  paste0("   ", file.path(.cv_report, "data_identity_monoamines.csv")),
  paste0("   ", file.path(.cv_report, "refit_equivalence.csv")),
  paste0("   ", file.path(.cv_report, "refit_sampling_log.csv")),
  "",
  "============================================================================"
)
writeLines(.summary_lines, file.path(.cv_report, "SUMMARY.txt"))
cat(paste(.summary_lines, collapse = "\n"), "\n", sep = "")

invisible(list(layer1 = list(behaviour = L1_beh, cortisol = L1_cor,
                              monoamines = L1_mon),
               layer2 = L2,
               report_dir = .cv_report,
               overall = .banner_overall))
