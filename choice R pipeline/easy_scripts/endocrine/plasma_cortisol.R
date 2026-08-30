# =============================================================================
# NON-AUTHORITATIVE / SUPERSEDED (stats-revision, 2026-07-06)
# This is a simplified mirror, not a source of manuscript numbers. The
# authoritative engines are:
#   Behaviour : choice R pipeline/scripts/01_pipeline_analysis/
#               activity_analysis_STATS_choice_exp.R
#   Endocrine : choice exp for claude mono and cortisol/claude output/Scripts/
#               analysis_b3.R
# Those two scripts carry the current M1-M5 statistical-revision fixes
# (BH-FDR by brain region, trial-level RE nesting, order/day sensitivity
# test, behaviour<->neurochemistry association). This file has NOT been
# updated to match and may report stale numbers. Do not cite its output in
# the manuscript. See D:\CHOICE R SCRIPTS\STATS_REVISION_INSTRUCTIONS.md
# and stat_map.md for the current source-of-truth mapping.
# =============================================================================

# easy_scripts/endocrine/plasma_cortisol.R
# -----------------------------------------------------------------------------
# Reproduces the §4 plasma cortisol primary model from analysis_b3.R.
# Response: log(plasma_cortisol). Fixed: treatment (=condition).
# RE candidates: 6 combinations of (1|tank), (1|trial), (1|plate), AICc-selected.
# Self-contained — does not source the main STATS file (cortisol has its own
# RE design separate from the behaviour models).
# Line-by-line runnable in RStudio.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr)
  library(lme4); library(lmerTest); library(MuMIn)
  library(emmeans); library(multcomp); library(ggplot2)
})

ROOT <- file.path(PROJECT_ROOT, "choice R pipeline/easy_scripts")
source(file.path(ROOT, "_data_access.R"))
source(file.path(ROOT, "_re_candidates.R"))   # gives RE_CORT_CAND

# ---- Step 1: Load -----------------------------------------------------------
dat <- load_endocrine_dataset()
message("Step 1 — Load: ", nrow(dat), " rows x ", ncol(dat), " cols")

# ---- Step 2: Filter to cortisol rows ----------------------------------------
d <- dat %>%
  dplyr::filter(analyte == "cort",
                is.finite(value),
                value > 0) %>%
  dplyr::mutate(log_cort = log(value))
message("Step 2 — Filter+log: ", nrow(d), " cortisol samples")

# ---- Step 3: Factor levels --------------------------------------------------
d$treatment <- factor(d$treatment, levels = c("control", "treat"))
d$tank      <- factor(d$tank)
d$trial     <- factor(d$trial)
d$plate     <- factor(d$plate)
message("Step 3 — Factors set; treatment levels: ",
        paste(levels(d$treatment), collapse = ", "))

# ---- Step 4: Transformation diagnostics -------------------------------------
cat("Raw plasma_cortisol (ng/mL):\n")
cat("  mean:           ", signif(mean(d$value), 3), "\n")
cat("  sd:             ", signif(sd(d$value), 3), "\n")
cat("  Shapiro raw p:  ", signif(shapiro.test(d$value)$p.value, 3), "\n")
cat("  Shapiro log p:  ", signif(shapiro.test(d$log_cort)$p.value, 3), "\n")
cat("  → transform applied: log()\n")
message("Step 4 — Diagnostics emitted above")

# ---- Step 5: Candidate RE structures ----------------------------------------
message("Step 5 — RE candidates (RE_CORT_CAND): ",
        paste(names(RE_CORT_CAND), collapse = "; "))

# ---- Step 6: Fit + AICc-select RE -------------------------------------------
# Mirrors select_re_aicc() from analysis_b3.R: fit each candidate, drop
# singular fits, pick lowest AICc.
.fit_one <- function(re_text) {
  formula_str <- paste("log_cort ~ condition +", re_text)
  # Use original column name 'condition' that lmer expects:
  d2 <- d %>% dplyr::rename(condition = treatment)
  f <- as.formula(formula_str)
  fit <- tryCatch(suppressMessages(lmer(f, data = d2, REML = TRUE)),
                  error = function(e) NULL)
  if (is.null(fit) || isSingular(fit)) return(NULL)
  list(re = re_text, fit = fit, aicc = MuMIn::AICc(fit))
}
.fits <- purrr::compact(lapply(RE_CORT_CAND, .fit_one))
if (length(.fits) == 0) {
  message("  All RE candidates singular — falling back to lm")
  fit_cort  <- lm(log_cort ~ condition, data = d %>% dplyr::rename(condition = treatment))
  re_label  <- "lm (no RE)"
} else {
  .aicc_tbl <- data.frame(
    re   = sapply(.fits, `[[`, "re"),
    AICc = sapply(.fits, `[[`, "aicc")
  ) %>% dplyr::arrange(AICc)
  cat("\n--- AICc table ---\n"); print(.aicc_tbl)
  .winner <- .fits[[which.min(sapply(.fits, `[[`, "aicc"))]]
  fit_cort <- .winner$fit
  re_label  <- .winner$re
}
message("Step 6 — Selected RE: ", re_label)

# ---- Step 7: Type-III ANOVA (Kenward-Roger) ---------------------------------
if (inherits(fit_cort, "lmerMod")) {
  cort_anova <- tryCatch(anova(fit_cort, ddf = "Kenward-Roger"),
                         error = function(e) anova(fit_cort))
} else {
  cort_anova <- anova(fit_cort)
}
cat("\n--- Type III ANOVA (Kenward-Roger) ---\n"); print(cort_anova)
message("Step 7 — Type III ANOVA printed above")

# ---- Step 8: Tukey + CLD ----------------------------------------------------
emm <- emmeans::emmeans(fit_cort, ~ condition)
cat("\n--- emmeans (log scale) ---\n"); print(emm)
contr <- pairs(emm, adjust = "tukey")
cat("\n--- Tukey contrast ---\n"); print(contr)
cld <- multcomp::cld(emm, Letters = letters, adjust = "tukey")
cat("\n--- CLD ---\n"); print(cld)
message("Step 8 — Tukey + CLD printed above")

# ---- Step 8b: Significant Tukey contrasts only (p < 0.05) ------------------
cat("\n=== SIGNIFICANT CONTRASTS (p < 0.05) ===\n")
.ct8b <- tryCatch(as.data.frame(summary(contr)), error = function(e) NULL)
if (!is.null(.ct8b) && "p.value" %in% names(.ct8b)) {
  .sig8b <- .ct8b[.ct8b$p.value < 0.05, , drop = FALSE]
  if (nrow(.sig8b) > 0) print(.sig8b) else cat("  (no significant contrasts at p < 0.05)\n")
} else {
  cat("  (could not extract p.value from contrasts)\n")
}
message("Step 8b — Significant Tukey contrasts (p < 0.05) printed above")

# ---- Step 9: Raw summary ----------------------------------------------------
raw <- d %>%
  dplyr::group_by(treatment) %>%
  dplyr::summarise(n = dplyr::n(),
                   mean_ng_ml = mean(value, na.rm = TRUE),
                   sd_ng_ml   = sd(value,   na.rm = TRUE),
                   sem_ng_ml  = sd_ng_ml / sqrt(n),
                   .groups = "drop")
cat("\n--- raw mean +/- SD (treatment, ng/mL) ---\n"); print(raw)
message("Step 9 — Raw summary printed above")

# ---- Step 10: Save graph (legacy jitter style) ------------------------------
# Mirrors p_cort_leg recipe from analysis_b3.R but simplified (no significance bar).
pal_treat <- c(control = "#2166AC", treat = "#D6604D")
p_cort <- ggplot(d, aes(treatment, value, color = treatment)) +
  geom_jitter(width = 0.30, size = 3.0, alpha = 0.85) +
  geom_errorbar(data = raw,
                aes(x = treatment, y = mean_ng_ml,
                    ymin = mean_ng_ml - sem_ng_ml,
                    ymax = mean_ng_ml + sem_ng_ml),
                width = 0.15, linewidth = 1.0, colour = "black",
                inherit.aes = FALSE) +
  geom_point(data = raw, aes(x = treatment, y = mean_ng_ml),
             size = 2.5, colour = "black", inherit.aes = FALSE) +
  scale_color_manual(values = pal_treat, guide = "none") +
  labs(x = NULL, y = "[Cortisol] (ng/mL)") +
  theme_bw(base_size = 12) +
  theme(aspect.ratio = 1.2)

out_dir <- file.path(ROOT, "outputs", "graphs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_png <- file.path(out_dir, "plasma_cortisol.png")
ggplot2::ggsave(out_png, p_cort, width = 110, height = 130, units = "mm", dpi = 300)
message("Step 10 — Saved graph to ", out_png)
