# =============================================================================
# build_ruth_statement.R
# -----------------------------------------------------------------------------
# Builds D:/CHOICE R SCRIPTS/choice R pipeline/Ruth_statement.docx — a single
# manuscript-friendly table documenting the FINAL fitted model behind every
# presented inferential result across two pipelines.
#
# USAGE (must run in a single fresh R session):
#   1) source the main behavioural pipeline (or its run_manu_graphs_only wrapper)
#        source(file.path(PROJECT_ROOT, "choice R pipeline/scripts/00_main/run_manu_graphs_only.R"))
#   2) source the endocrine pipeline
#        source(file.path(PROJECT_ROOT, "choice exp for claude mono and cortisol/claude output/Scripts/analysis_b3_REVISED.R"))
#   3) source this file
#        source(file.path(PROJECT_ROOT, "choice R pipeline/scripts/04_reporting/build_ruth_statement.R"))
#
# The helper guard at line ~1758 of the main pipeline prevents sys.source()-based
# helper-only hydration of model objects, so single-session execution is the
# only viable path.
#
# OUTPUT:
#   D:/CHOICE R SCRIPTS/choice R pipeline/Ruth_statement.docx  (overwritten)
# =============================================================================

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(purrr)
  library(officer)
  library(flextable)
  library(lme4)
})

.RUTH_DIR <- file.path(PROJECT_ROOT, "choice R pipeline/scripts/04_reporting")
.RUTH_OUT <- file.path(PROJECT_ROOT, "choice R pipeline/Ruth_statement.docx")

source(file.path(.RUTH_DIR, "extractors.R"))
source(file.path(.RUTH_DIR, "manifest_indicators.R"))

# ---- 1. Harvest -------------------------------------------------------------
message("[Ruth_statement] Harvesting ", nrow(MANIFEST), " indicators...")
master <- harvest_manifest(MANIFEST, envir = globalenv())
message("[Ruth_statement] Harvested ", nrow(master), " rows.")

# ---- 2. Validation (soft) ---------------------------------------------------
.validation_log <- character(0)
.add_warn <- function(msg) {
  message("[Ruth_statement][warn] ", msg)
  .validation_log <<- c(.validation_log, msg)
}

# Check: any rows where the fit was NULL or class unsupported
.null_rows <- master[is.na(master$fixed_effects) |
                     grepl("NULL or does not exist", master$selection_basis), , drop = FALSE]
if (nrow(.null_rows) > 0) {
  for (i in seq_len(nrow(.null_rows))) {
    .add_warn(sprintf("Indicator '%s' (%s): model object missing or unsupported.",
                      .null_rows$indicator[i], .null_rows$block[i]))
  }
}

# Check: empty fixed-effects column (after the NULL filter)
.bad_fixed <- master[!is.na(master$fixed_effects) & nchar(master$fixed_effects) == 0, , drop = FALSE]
if (nrow(.bad_fixed) > 0) {
  for (i in seq_len(nrow(.bad_fixed))) {
    .add_warn(sprintf("Indicator '%s': empty Fixed effects column.",
                      .bad_fixed$indicator[i]))
  }
}

# Optional cross-check vs figure stems (soft):
.fig_dirs <- c(
  file.path(PROJECT_ROOT, "choice R pipeline/manu_graphs"),
  file.path(PROJECT_ROOT, "all_manu_graphs")
)
.fig_dirs <- .fig_dirs[dir.exists(.fig_dirs)]
if (length(.fig_dirs) > 0) {
  .pngs <- unlist(lapply(.fig_dirs, list.files,
                         pattern = "\\.png$", full.names = FALSE,
                         recursive = TRUE), use.names = FALSE)
  .pngs <- unique(tools::file_path_sans_ext(.pngs))
  if (length(.pngs) > 0)
    message(sprintf("[Ruth_statement] Discovered %d PNG stems in manu_graphs dirs (informational).",
                    length(.pngs)))
}

# ---- 3. Pretty-format the master table for Word -----------------------------
master_display <- master %>%
  transmute(
    `Analysis block`                  = block,
    `Indicator`                       = indicator,
    `Response variable`               = response_pretty,
    `Analysis unit`                   = unit,
    `Fixed effects`                   = ifelse(is.na(fixed_effects), "—", fixed_effects),
    `Random effects`                  = ifelse(is.na(random_effects), "—", random_effects),
    `Distribution / family`           = ifelse(is.na(family), "—", family),
    `Link function`                   = ifelse(is.na(link), "—", link),
    `Response scale / transformation` = ifelse(is.na(transform), "—", transform),
    `Final-model selection basis`     = ifelse(is.na(selection_basis), "—", selection_basis),
    `Source`                          = source
  ) %>%
  arrange(`Analysis block`, `Indicator`)

# ---- 4. Build the Word document ---------------------------------------------
.session_lines <- c(
  sprintf("R version: %s", paste(R.version$major, R.version$minor, sep = ".")),
  sprintf("Date built: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  sprintf("Total rows: %d  (primary: %d, descriptive: %d)",
          nrow(master),
          sum(MANIFEST$tier == "primary"),
          sum(MANIFEST$tier == "descriptive"))
)
.pkg_lines <- vapply(c("lme4","lmerTest","glmmTMB","emmeans","MuMIn","multcomp",
                       "officer","flextable","performance"),
                     function(p) {
                       v <- tryCatch(as.character(utils::packageVersion(p)),
                                     error = function(e) "not installed")
                       sprintf("  %s: %s", p, v)
                     }, character(1))

ft <- flextable(master_display) %>%
  set_table_properties(layout = "autofit", width = 1) %>%
  bold(part = "header") %>%
  bg(part = "header", bg = "#D9D9D9") %>%
  fontsize(part = "all", size = 8) %>%
  padding(part = "all", padding = 2) %>%
  border_outer(part = "all", border = officer::fp_border(color = "#888888", width = 0.5)) %>%
  border_inner(part = "all", border = officer::fp_border(color = "#CCCCCC", width = 0.25))

doc <- read_docx() %>%
  body_add_par("Ruth's statement — final-model documentation",
               style = "heading 1") %>%
  body_add_par("Final fitted model behind every presented inferential result in the choice experiment pipelines. One row per indicator. Models marked [Primary] are reported as the headline inference; [Descriptive] models are reported descriptively alongside.",
               style = "Normal") %>%
  body_add_par("Session information", style = "heading 2")
for (line in .session_lines) doc <- doc %>% body_add_par(line, style = "Normal")
doc <- doc %>% body_add_par("Package versions:", style = "Normal")
for (line in .pkg_lines)     doc <- doc %>% body_add_par(line, style = "Normal")

doc <- doc %>%
  body_add_par(" ", style = "Normal") %>%
  body_add_par("Master table", style = "heading 2") %>%
  body_add_flextable(ft) %>%
  body_add_par(" ", style = "Normal") %>%
  body_add_par("Glossary", style = "heading 2") %>%
  body_add_par("KR = Kenward-Roger; CLR = centred log-ratio; ALR = additive log-ratio; OLS = ordinary least squares; AICc = corrected Akaike Information Criterion; β-GLMM = beta generalised linear mixed model; RE = random effects.",
               style = "Normal") %>%
  body_add_par("Per-cell monoamine models may use heterogeneous Box-Cox λ values across cells of the same analyte; β estimates and EMMs are not directly comparable across cells without back-transformation.",
               style = "Normal")

if (length(.validation_log) > 0) {
  doc <- doc %>%
    body_add_par(" ", style = "Normal") %>%
    body_add_par("Validation notes", style = "heading 2")
  for (msg in .validation_log) doc <- doc %>% body_add_par(msg, style = "Normal")
} else {
  doc <- doc %>%
    body_add_par(" ", style = "Normal") %>%
    body_add_par("Validation: clean run — no warnings.", style = "Normal")
}

# Landscape page for the wide table
doc <- doc %>%
  body_end_section_continuous() %>%
  body_set_default_section(
    prop_section(
      page_size = page_size(width = 11.69, height = 8.27, orient = "landscape"),
      page_margins = page_mar(top = 0.4, bottom = 0.4, left = 0.4, right = 0.4)
    )
  )

print(doc, target = .RUTH_OUT)
message(sprintf("[Ruth_statement] Wrote: %s", .RUTH_OUT))
invisible(master)
