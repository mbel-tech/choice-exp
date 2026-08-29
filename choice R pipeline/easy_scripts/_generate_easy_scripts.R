# easy_scripts/_generate_easy_scripts.R
# -----------------------------------------------------------------------------
# Stamps out the ~50 routine mini-scripts from a curated spec list using
# _template_mini_script.R as the source-of-truth template.
#
# Hand-written mini-scripts (mean_nnd_cm.R, logit_flow_main_by_tp.R,
# plasma_cortisol.R) serve as reference; the generator covers the long tail
# (NND/Pol/IID/Hull/Speed/Switches/Active/Flux + zone alr + Jacobs + cells +
# 7 monoamines + 28 per-area cells).
#
# Idempotent — overwrites existing mini-scripts each run.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(glue); library(dplyr); library(readr); library(stringr)
})

ROOT <- file.path(PROJECT_ROOT, "choice R pipeline/easy_scripts")
TEMPLATE_FILE <- file.path(ROOT, "_template_mini_script.R")
stopifnot(file.exists(TEMPLATE_FILE))
.template_text <- paste(readLines(TEMPLATE_FILE), collapse = "\n")

# -----------------------------------------------------------------------------
# Spec table — one row per mini-script
# -----------------------------------------------------------------------------
# Columns:
#   folder, filename, res_obj, line_range, csv_name,
#   filter_expr, indicator_label, additional_factor_lines,
#   transform_diag_expr, re_cands, model_call, raw_grouping,
#   response, plot_builder_call, graph_filename, width_mm, height_mm

# ---- Behaviour aggregated (timepoint_aggregated/) ---------------------------
.spec_agg <- tibble::tribble(
  ~response,              ~indicator_label,                ~y_axis_label,
  "switches_per_session", "Zone switches per session",     "Switches/session",
  "zone_flux_per_session","Zone flux per session",         "Flux/session",
  "prop_active",          "Proportion of time active",     "Proportion active",
  "logit_flow",           "ALR(flow vs calm)",             "ALR(flow vs calm)",
  "lr_high",              "ALR(high vs calm)",             "ALR(high vs calm)",
  "lr_medium",            "ALR(medium vs calm)",           "ALR(medium vs calm)",
  "lr_low",               "ALR(low vs calm)",              "ALR(low vs calm)",
  "mean_nnd_cm",          "Mean nearest-neighbour distance (cm)", "Mean NND (cm)",
  "mean_polarisation",    "Mean polarisation",             "Polarisation",
  "mean_iid_cm",          "Mean inter-individual distance (cm)",  "Mean IID (cm)",
  "mean_hull_area_cm2",   "Mean convex hull area (cm²)",   "School area (cm²)",
  "mean_centroid_spd_cm", "Mean school centroid speed (cm/s)",    "School speed (cm/s)"
)

# ---- Behaviour by-timepoint (by_timepoint/) ---------------------------------
.spec_tp <- .spec_agg  # same indicators, plus timepoint_f interaction

# Filename + helper choice per response
.builder_for <- function(resp) {
  if (resp %in% c("logit_flow", "lr_high", "lr_medium", "lr_low")) "alr_scatter"
  else if (resp == "prop_active") "beta_collective"
  else "lmm_collective"
}

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
.csv_behavior <- "easy_scripts_dataset.csv"
.csv_endo     <- "easy_scripts_endo_dataset.csv"

.fmt_filter_agg <- function(resp) sprintf(
  'dat %%>%%
  dplyr::filter(is.finite(%s)) %%>%%
  dplyr::group_by(phys_trial_id, tank, trial_date, treatment, fish_density_f) %%>%%
  dplyr::summarise(%s = mean(%s, na.rm = TRUE), .groups = "drop")',
  resp, resp, resp
)

.fmt_filter_tp <- function(resp) sprintf(
  'dat %%>%%
  dplyr::filter(is.finite(%s)) %%>%%
  dplyr::select(trial_id, phys_trial_id, tank, trial_date, treatment,
                timepoint, timepoint_f, fish_density_f, %s)',
  resp, resp
)

.fmt_model_call_lmm_agg <- function(resp) sprintf(
  'run_lmm_analysis(
  label     = "%s_aggregated",
  data      = d,
  response  = "%s",
  fixed_str = "treatment",
  re_cands  = RE_WIDE_AGG,
  focal_terms = list("Treatment" = "^treatment$")
)', resp, resp
)

.fmt_model_call_lmm_tp <- function(resp) sprintf(
  'run_lmm_analysis(
  label          = "%s_timepoint",
  data           = d,
  response       = "%s",
  fixed_str      = "treatment * timepoint_f",
  cont_fixed_str = "treatment * timepoint",
  re_cands       = RE_WIDE,
  focal_terms = list(
    "Treatment"           = "^treatment$",
    "Timepoint"           = "^timepoint",
    "Treatment:Timepoint" = "treatment.*timepoint|timepoint.*treatment"
  )
)', resp, resp
)

.fmt_model_call_beta_agg <- function(resp) sprintf(
  'run_betaglmm_analysis(
  label     = "%s_aggregated",
  data      = d,
  response  = "%s",
  fixed_str = "treatment",
  re_cands  = RE_WIDE_AGG,
  focal_terms = list("Treatment" = "^treatment$")
)', resp, resp
)

.fmt_plot_call <- function(resp, y_label, builder) {
  if (builder == "alr_scatter") {
    sprintf(
      '.mg_make_alr_scatter(df_in = d, y_col = "%s", y_label = "%s",
                       res_obj = res, show_n = FALSE)',
      resp, y_label
    )
  } else {
    sprintf(
      '.mg_treatment_xaxis(.mg_mk_collective_panel(
  res, y_col = "%s", y_label = "%s",
  dat = d, show_n = FALSE
))',
      resp, y_label
    )
  }
}

.fmt_transform_diag <- function(resp) sprintf(
  'cat("Variance:        ", signif(var(d$%s, na.rm = TRUE), 3), "\\n")
sw_p <- tryCatch(shapiro.test(d$%s)$p.value, error = function(e) NA_real_)
cat("Shapiro-Wilk p:  ", signif(sw_p, 3), "\\n")
cat("N per treatment:\\n"); print(table(d$treatment))',
  resp, resp
)

# -----------------------------------------------------------------------------
# Stamp one mini-script
# -----------------------------------------------------------------------------
.stamp_one <- function(folder, filename, vars) {
  txt <- .template_text
  for (k in names(vars)) {
    txt <- gsub(paste0("{{", k, "}}"), vars[[k]], txt, fixed = TRUE)
  }
  out <- file.path(ROOT, folder, filename)
  dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
  writeLines(txt, out)
  invisible(out)
}

# -----------------------------------------------------------------------------
# Generate behaviour aggregated mini-scripts
# -----------------------------------------------------------------------------
.n_agg <- 0L
for (i in seq_len(nrow(.spec_agg))) {
  s <- .spec_agg[i, ]
  resp <- s$response
  builder <- .builder_for(resp)
  fname <- paste0(resp, ".R")
  vars <- list(
    FOLDER                  = "timepoint_aggregated",
    FILENAME                = fname,
    RES_OBJECT              = paste0("res_", sub("mean_|prop_|logit_|lr_", "", resp), "_agg"),
    LINE_RANGE              = "see §5/§5b",
    CSV_NAME                = .csv_behavior,
    FILTER_AND_AGGREGATE_EXPR = .fmt_filter_agg(resp),
    INDICATOR_LABEL         = s$indicator_label,
    ADDITIONAL_FACTOR_LINES = 'd$tank <- factor(d$tank)\nd$fish_density_f <- factor(d$fish_density_f, levels = as.character(DENSITY_LEVELS_g))',
    TRANSFORM_DIAGNOSTIC_EXPR = .fmt_transform_diag(resp),
    RE_CANDS                = "RE_WIDE_AGG",
    MODEL_CALL              = if (builder == "beta_collective") .fmt_model_call_beta_agg(resp)
                              else .fmt_model_call_lmm_agg(resp),
    RAW_GROUPING            = "treatment",
    RESPONSE                = resp,
    PLOT_BUILDER_CALL       = .fmt_plot_call(resp, s$y_axis_label, builder),
    GRAPH_FILENAME          = paste0(resp, "_agg"),
    WIDTH_MM                = "100",
    HEIGHT_MM               = "100"
  )
  .stamp_one("timepoint_aggregated", fname, vars)
  .n_agg <- .n_agg + 1L
}
message(sprintf("[generator] timepoint_aggregated/ : %d scripts written", .n_agg))

# -----------------------------------------------------------------------------
# Generate behaviour by-timepoint mini-scripts
# -----------------------------------------------------------------------------
.n_tp <- 0L
for (i in seq_len(nrow(.spec_tp))) {
  s <- .spec_tp[i, ]
  resp <- s$response
  builder <- .builder_for(resp)
  fname <- paste0(resp, "_by_tp.R")
  vars <- list(
    FOLDER                  = "by_timepoint",
    FILENAME                = fname,
    RES_OBJECT              = paste0("res_", sub("mean_|prop_|logit_|lr_", "", resp), "_tp"),
    LINE_RANGE              = "see §5/§5b",
    CSV_NAME                = .csv_behavior,
    FILTER_AND_AGGREGATE_EXPR = .fmt_filter_tp(resp),
    INDICATOR_LABEL         = paste(s$indicator_label, "(by interval)"),
    ADDITIONAL_FACTOR_LINES = 'd$tank <- factor(d$tank)\nd$timepoint_f <- factor(d$timepoint_f, levels = TIMEPOINT_LEVELS_g)\nd$fish_density_f <- factor(d$fish_density_f, levels = as.character(DENSITY_LEVELS_g))',
    TRANSFORM_DIAGNOSTIC_EXPR = .fmt_transform_diag(resp),
    RE_CANDS                = "RE_WIDE",
    MODEL_CALL              = if (builder == "beta_collective")
                                gsub("RE_WIDE_AGG", "RE_WIDE",
                                     .fmt_model_call_beta_agg(resp), fixed = TRUE)
                              else .fmt_model_call_lmm_tp(resp),
    RAW_GROUPING            = "treatment, timepoint_f",
    RESPONSE                = resp,
    PLOT_BUILDER_CALL       = .fmt_plot_call(resp, s$y_axis_label, builder),
    GRAPH_FILENAME          = paste0(resp, "_by_tp"),
    WIDTH_MM                = "140",
    HEIGHT_MM               = "100"
  )
  .stamp_one("by_timepoint", fname, vars)
  .n_tp <- .n_tp + 1L
}
message(sprintf("[generator] by_timepoint/         : %d scripts written", .n_tp))

# -----------------------------------------------------------------------------
# Generate per-analyte monoamine mini-scripts (treatment × area)
# -----------------------------------------------------------------------------
.analytes <- c(
  ht_5            = "5-HT",
  hiaa_5          = "5-HIAA",
  hiaa_5_ratio    = "5-HIAA/5-HT",
  da              = "DA",
  dopac           = "DOPAC",
  dopac_da_ratio  = "DOPAC/DA",
  ne              = "NE"
)
.endo_re_text <- 'list("(1|sample_id)"=lmer(log_v~treatment*area+(1|sample_id),data=d,REML=TRUE))'
# Endocrine mini-script template (different from behaviour — uses lmer directly)
.endo_template <- '# easy_scripts/endocrine/{{FILENAME}}
# Reproduces §4/§5 per-analyte monoamine analysis from analysis_b3.R.
# Indicator: {{INDICATOR_LABEL}}
# Response: log(value); Fixed: treatment * area.
# Self-contained — does not source the main STATS bridge.
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(lme4); library(lmerTest)
  library(MuMIn); library(emmeans); library(multcomp); library(ggplot2)
})
ROOT <- file.path(PROJECT_ROOT, "choice R pipeline/easy_scripts")

# Step 1: Load
dat <- readr::read_csv(file.path(ROOT, "easy_scripts_endo_dataset.csv"),
                       show_col_types = FALSE)
message("Step 1 - Load: ", nrow(dat), " rows x ", ncol(dat), " cols")

# Step 2: Filter
d <- dat %>%
  dplyr::filter(analyte == "{{ANALYTE_KEY}}", is.finite(value), value > 0) %>%
  dplyr::mutate(log_v = log(value))
message("Step 2 - Filter: ", nrow(d), " {{INDICATOR_LABEL}} measurements")

# Step 3: Factors
d$treatment <- factor(d$treatment, levels = c("control","treat"))
d$area      <- factor(d$area)
d$tank      <- factor(d$tank)
d$sample_id <- factor(d$sample_id)
options(contrasts = c("contr.sum","contr.poly"))
message("Step 3 - Factors set; areas: ", paste(levels(d$area), collapse=", "))

# Step 4: Diagnostics
cat("Variance raw:  ", signif(var(d$value), 3), "\\n")
cat("Variance log:  ", signif(var(d$log_v), 3), "\\n")
cat("N per cell:\\n"); print(table(d$treatment, d$area))
message("Step 4 - Diagnostics emitted above")

# Step 5: RE candidates
.re_cands <- c(
  "(1|sample_id)" = "(1|sample_id)",
  "(1|tank)"      = "(1|tank)",
  "(1|sample_id)+(1|tank)" = "(1|sample_id) + (1|tank)"
)
message("Step 5 - RE candidates: ", paste(names(.re_cands), collapse="; "))

# Step 6: AICc-select
.fits <- list()
for (k in names(.re_cands)) {
  f <- as.formula(paste("log_v ~ treatment * area +", .re_cands[[k]]))
  fit <- tryCatch(suppressMessages(lmer(f, data = d, REML = TRUE)),
                  error = function(e) NULL)
  if (!is.null(fit) && !isSingular(fit)) {
    .fits[[k]] <- list(fit = fit, aicc = MuMIn::AICc(fit))
  }
}
if (length(.fits) == 0) {
  fit <- lm(log_v ~ treatment * area, data = d); re_label <- "lm (singular)"
} else {
  .aicc <- sapply(.fits, function(x) x$aicc)
  fit   <- .fits[[which.min(.aicc)]]$fit; re_label <- names(.fits)[which.min(.aicc)]
  cat("\\n-- AICc --\\n"); print(.aicc)
}
message("Step 6 - Selected RE: ", re_label)

# Step 7: ANOVA
anv <- tryCatch(anova(fit, ddf="Kenward-Roger"), error = function(e) anova(fit))
cat("\\n-- Type III ANOVA --\\n"); print(anv)
message("Step 7 - ANOVA printed above")

# Step 8: Tukey + CLD per area
emm <- emmeans::emmeans(fit, ~ treatment | area)
cat("\\n-- emmeans (log scale) --\\n"); print(emm)
contr <- pairs(emm, adjust = "tukey")
cat("\\n-- Tukey contrasts --\\n"); print(contr)
cld <- multcomp::cld(emm, Letters = letters, adjust = "tukey")
cat("\\n-- CLD --\\n"); print(cld)
message("Step 8 - Tukey + CLD printed above")

# Step 8b: Significant Tukey contrasts only (p < 0.05)
cat("\\n=== SIGNIFICANT CONTRASTS (p < 0.05) ===\\n")
.ct8b <- tryCatch(as.data.frame(summary(contr)), error = function(e) NULL)
if (!is.null(.ct8b) && "p.value" %in% names(.ct8b)) {
  .sig8b <- .ct8b[.ct8b$p.value < 0.05, , drop = FALSE]
  if (nrow(.sig8b) > 0) print(.sig8b) else cat("  (no significant contrasts at p < 0.05)\\n")
} else {
  cat("  (could not extract p.value from contrasts)\\n")
}
message("Step 8b - Significant Tukey contrasts (p < 0.05) printed above")

# Step 9: Raw summary
raw <- d %>% dplyr::group_by(treatment, area) %>%
  dplyr::summarise(n=dplyr::n(), mean=mean(value), sd=sd(value), .groups="drop")
cat("\\n-- raw mean+/-SD --\\n"); print(raw)
message("Step 9 - Raw summary printed above")

# Step 10: Save graph
pal <- c(control = "#2166AC", treat = "#D6604D")
p <- ggplot(d, aes(area, value, color = treatment)) +
  geom_point(size = 2.5, alpha = 0.7,
             position = position_jitterdodge(jitter.width = 0.20,
                                              dodge.width = 0.6)) +
  stat_summary(fun = mean, geom = "point", size = 3, shape = 18,
               position = position_dodge(width = 0.6), color = "black") +
  scale_color_manual(values = pal) +
  labs(x = "Brain area", y = "{{INDICATOR_LABEL}}", color = "Treatment") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
out_dir <- file.path(ROOT, "outputs", "graphs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
ggplot2::ggsave(file.path(out_dir, "monoamine_{{ANALYTE_SAFE}}.png"),
                p, width = 140, height = 110, units = "mm", dpi = 300)
message("Step 10 - Saved graph")
'

.n_endo <- 0L
for (akey in names(.analytes)) {
  label <- .analytes[[akey]]
  safe <- gsub("[/-]", "_", tolower(label))
  fname <- paste0("monoamine_", safe, ".R")
  txt <- .endo_template
  txt <- gsub("{{FILENAME}}",        fname,  txt, fixed = TRUE)
  txt <- gsub("{{INDICATOR_LABEL}}", label,  txt, fixed = TRUE)
  txt <- gsub("{{ANALYTE_KEY}}",     akey,   txt, fixed = TRUE)
  txt <- gsub("{{ANALYTE_SAFE}}",    safe,   txt, fixed = TRUE)
  out <- file.path(ROOT, "endocrine", fname)
  writeLines(txt, out)
  .n_endo <- .n_endo + 1L
}
message(sprintf("[generator] endocrine/            : %d scripts written", .n_endo))

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
.total <- .n_agg + .n_tp + .n_endo
message(sprintf("[generator] Total                 : %d scripts written", .total))
invisible(NULL)
