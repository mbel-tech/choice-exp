setwd(file.path(PROJECT_ROOT, "choice R pipeline"))
STEP5_OUT <- normalizePath(
  file.path(PROJECT_ROOT, "choice R pipeline/STEP5_stats/STEP5_stats_20260511_175522"),
  mustWork = FALSE
)
source(file.path(PROJECT_ROOT, "choice R pipeline/scripts/01_pipeline_analysis/activity_analysis_STATS_choice_exp.R"))
cat("\n=== DIAGNOSTICS ===\n")
cat("df_main_cells_agg: ", if(is.null(df_main_cells_agg)) "NULL" else paste(nrow(df_main_cells_agg), "rows,", paste(names(df_main_cells_agg), collapse=",")), "\n")
cat("df_sec_cells_agg: ", if(is.null(df_sec_cells_agg)) "NULL" else paste(nrow(df_sec_cells_agg), "rows,", paste(names(df_sec_cells_agg), collapse=",")), "\n")
