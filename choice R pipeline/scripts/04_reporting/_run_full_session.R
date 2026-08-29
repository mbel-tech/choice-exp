# Orchestrator: runs both pipelines + the reporter in one R session.
# Use with: Rscript _run_full_session.R
message("[FULL SESSION] === STAGE 1: Pipeline A (run_manu_graphs_only) ===")
source(file.path(PROJECT_ROOT, "choice R pipeline/scripts/00_main/run_manu_graphs_only.R"))

message("[FULL SESSION] === STAGE 2: Pipeline B (analysis_b3) ===")
# analysis_b3.R fits all models but may fail during its own Word-report
# rendering (officer temp-file issue / locked docx). Catch that error so
# we can still harvest the model objects that are already in globalenv.
tryCatch(
  source(file.path(PROJECT_ROOT, "choice exp for claude mono and cortisol/claude output/Scripts/analysis_b3_REVISED.R")),
  error = function(e) {
    message("[FULL SESSION][warn] Pipeline B threw an error (likely during report",
            " rendering, not model fitting):\n  ", conditionMessage(e))
    # Verify the objects we actually need are present
    missing <- character(0)
    if (!exists("fit_cort",         envir = globalenv())) missing <- c(missing, "fit_cort")
    if (!exists("mono_cell_models", envir = globalenv())) missing <- c(missing, "mono_cell_models")
    if (length(missing) > 0) {
      stop("[FULL SESSION] Critical model object(s) absent after Pipeline B error: ",
           paste(missing, collapse = ", "), call. = FALSE)
    }
    message("[FULL SESSION] fit_cort and mono_cell_models present — continuing.")
  }
)

message("[FULL SESSION] === STAGE 3: Reporter (build_ruth_statement) ===")
source(file.path(PROJECT_ROOT, "choice R pipeline/scripts/04_reporting/build_ruth_statement.R"))

message("[FULL SESSION] === STAGE 4: Datasets workbook (build_datasets_workbook) ===")
tryCatch(
  source(file.path(PROJECT_ROOT, "choice R pipeline/scripts/06_workbook/build_datasets_workbook.R")),
  error = function(e) {
    message("[FULL SESSION][warn] Datasets-workbook build failed:\n  ",
            conditionMessage(e))
  }
)

message("[FULL SESSION] === STAGE 5: Cross-validate workbook vs pipeline ===")
tryCatch(
  source(file.path(PROJECT_ROOT, "choice R pipeline/scripts/07_cross_validate/cross_validate_workbook.R")),
  error = function(e) {
    message("[FULL SESSION][warn] Cross-validation failed:\n  ",
            conditionMessage(e))
  }
)

message("[FULL SESSION] Done.")
