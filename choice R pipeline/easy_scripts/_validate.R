# easy_scripts/_validate.R
# -----------------------------------------------------------------------------
# Smoke-tests every mini-script and confirms it produces a PNG.
# Run as:
#   source(file.path(PROJECT_ROOT, "choice R pipeline/easy_scripts/_validate.R"))
#
# Outputs a green/red summary table — does NOT compare model outputs to the
# master pipeline yet (that's the next iteration). For now, the goal is
# "every mini-script runs without error and emits its expected PNG."
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr); library(stringr)
})

ROOT <- file.path(PROJECT_ROOT, "choice R pipeline/easy_scripts")
OUT_DIR <- file.path(ROOT, "outputs", "graphs")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Discover all mini-scripts
.scripts <- c(
  list.files(file.path(ROOT, "timepoint_aggregated"), pattern = "\\.R$",
             full.names = TRUE),
  list.files(file.path(ROOT, "by_timepoint"), pattern = "\\.R$",
             full.names = TRUE),
  list.files(file.path(ROOT, "endocrine"), pattern = "\\.R$",
             full.names = TRUE)
)
.scripts <- .scripts[!grepl("^_", basename(.scripts))]

message(sprintf("[validate] %d mini-scripts to test", length(.scripts)))

.results <- data.frame(
  folder = character(0),
  script = character(0),
  status = character(0),
  message = character(0),
  png_exists = logical(0),
  png_kb = numeric(0)
)

for (.s in .scripts) {
  .folder <- basename(dirname(.s))
  .name   <- basename(.s)
  .stem   <- sub("\\.R$", "", .name)
  # Expected PNG name follows mini-script convention
  .png_candidates <- c(
    file.path(OUT_DIR, paste0(.stem, ".png")),
    file.path(OUT_DIR, paste0(.stem, "_agg.png")),
    file.path(OUT_DIR, paste0(.stem, "_by_tp.png")),
    file.path(OUT_DIR, paste0(sub("_by_tp$", "", .stem), "_by_tp.png"))
  )

  .start <- Sys.time()
  .res <- tryCatch({
    .env <- new.env(parent = globalenv())
    sys.source(.s, envir = .env)
    list(status = "PASS", message = "")
  }, error = function(e) {
    list(status = "FAIL", message = conditionMessage(e))
  })
  .elapsed <- as.numeric(difftime(Sys.time(), .start, units = "secs"))

  .png_path <- .png_candidates[file.exists(.png_candidates)][1]
  .png_exists <- !is.na(.png_path)
  .png_kb <- if (.png_exists) round(file.info(.png_path)$size / 1024, 1) else NA_real_

  .results <- rbind(.results, data.frame(
    folder = .folder, script = .name,
    status = .res$status, message = .res$message,
    png_exists = .png_exists, png_kb = .png_kb,
    elapsed_s = round(.elapsed, 1)
  ))

  cat(sprintf("  [%-22s] %-40s %s (%.1fs, %s)\n",
              .folder, .name, .res$status, .elapsed,
              if (.png_exists) sprintf("%.1f KB PNG", .png_kb) else "no PNG"))
}

cat("\n=== Validation summary ===\n")
print(.results %>% dplyr::count(folder, status))
cat(sprintf("\nTotal: %d / %d passed\n",
            sum(.results$status == "PASS"), nrow(.results)))

# Write CSV summary
.csv <- file.path(ROOT, "outputs", "validation_summary.csv")
readr::write_csv(.results, .csv)
cat("Summary written to", .csv, "\n")

invisible(.results)
