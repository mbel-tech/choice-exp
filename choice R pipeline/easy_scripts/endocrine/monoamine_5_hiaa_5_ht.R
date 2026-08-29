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

# easy_scripts/endocrine/monoamine_5_hiaa_5_ht.R
# Reproduces §4/§5 per-analyte monoamine analysis from analysis_b3.R.
# Indicator: 5-HIAA/5-HT
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
  dplyr::filter(analyte == "hiaa_5_ratio", is.finite(value), value > 0) %>%
  dplyr::mutate(log_v = log(value))
message("Step 2 - Filter: ", nrow(d), " 5-HIAA/5-HT measurements")

# Step 3: Factors
d$treatment <- factor(d$treatment, levels = c("control","treat"))
d$area      <- factor(d$area)
d$tank      <- factor(d$tank)
d$sample_id <- factor(d$sample_id)
options(contrasts = c("contr.sum","contr.poly"))
message("Step 3 - Factors set; areas: ", paste(levels(d$area), collapse=", "))

# Step 4: Diagnostics
cat("Variance raw:  ", signif(var(d$value), 3), "\n")
cat("Variance log:  ", signif(var(d$log_v), 3), "\n")
cat("N per cell:\n"); print(table(d$treatment, d$area))
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
  cat("\n-- AICc --\n"); print(.aicc)
}
message("Step 6 - Selected RE: ", re_label)

# Step 7: ANOVA
anv <- tryCatch(anova(fit, ddf="Kenward-Roger"), error = function(e) anova(fit))
cat("\n-- Type III ANOVA --\n"); print(anv)
message("Step 7 - ANOVA printed above")

# Step 8: Tukey + CLD per area
emm <- emmeans::emmeans(fit, ~ treatment | area)
cat("\n-- emmeans (log scale) --\n"); print(emm)
contr <- pairs(emm, adjust = "tukey")
cat("\n-- Tukey contrasts --\n"); print(contr)
cld <- multcomp::cld(emm, Letters = letters, adjust = "tukey")
cat("\n-- CLD --\n"); print(cld)
message("Step 8 - Tukey + CLD printed above")

# Step 8b: Significant Tukey contrasts only (p < 0.05)
cat("\n=== SIGNIFICANT CONTRASTS (p < 0.05) ===\n")
.ct8b <- tryCatch(as.data.frame(summary(contr)), error = function(e) NULL)
if (!is.null(.ct8b) && "p.value" %in% names(.ct8b)) {
  .sig8b <- .ct8b[.ct8b$p.value < 0.05, , drop = FALSE]
  if (nrow(.sig8b) > 0) print(.sig8b) else cat("  (no significant contrasts at p < 0.05)\n")
} else {
  cat("  (could not extract p.value from contrasts)\n")
}
message("Step 8b - Significant Tukey contrasts (p < 0.05) printed above")

# Step 9: Raw summary
raw <- d %>% dplyr::group_by(treatment, area) %>%
  dplyr::summarise(n=dplyr::n(), mean=mean(value), sd=sd(value), .groups="drop")
cat("\n-- raw mean+/-SD --\n"); print(raw)
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
  labs(x = "Brain area", y = "5-HIAA/5-HT", color = "Treatment") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
out_dir <- file.path(ROOT, "outputs", "graphs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
ggplot2::ggsave(file.path(out_dir, "monoamine_5_hiaa_5_ht.png"),
                p, width = 140, height = 110, units = "mm", dpi = 300)
message("Step 10 - Saved graph")

