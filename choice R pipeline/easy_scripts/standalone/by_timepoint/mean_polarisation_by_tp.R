# =============================================================================
# Standalone easy-script: mean polarisation by interval
# -----------------------------------------------------------------------------
# Self-contained reproduction of the Pipeline-A model for res_pol_tp.
# Source no project file. Uses only CRAN packages.
#
# Why this exists — see the header of mean_iid_cm_by_tp.R for the full story.
# Short version: the figure, the CSV, and the post-hocs differ because
# three independent transformations sit between raw data and the letters
# you see on Figure_10 (panel B):
#   1) response transform = sqrt
#   2) timepoint is treated as CONTINUOUS in the chosen model
#   3) display CLD letters are remapped by precedence — for polarisation
#      this remap is NON-IDENTITY (raw `c` -> display `a`, etc.)
#
# Pipeline source: scripts/01_pipeline_analysis/activity_analysis_STATS_choice_exp.R
#   model fit: line 2942  (res_pol_tp <- run_lmm_analysis(...))
#   AICc table: STEP5_stats/STEP5_stats_20260511_175522/polarisation_timepoint/aicc_selection.csv
#   normality: STEP5_stats/STEP5_stats_20260511_175522/polarisation_timepoint/normality_check.csv
#
# Hard-coded pipeline choices (from the May-11 outputs):
#   transform     = sqrt
#   fixed formula = sqrt(y) ~ treatment * timepoint        # continuous t
#   random effect = (1 | trial_date)
#   selected AICc = -189.36
# Note: six RE candidates TIE at AICc = -189.36 (within 0.005). The wrapper
# picks (1|trial_date) by its position in the candidate list — the choice
# is essentially arbitrary at the model-comparison level; any of the six
# would give equivalent inference.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(lme4); library(lmerTest); library(pbkrtest)
  library(emmeans); library(multcomp); library(ggplot2); library(car)
})

# ---- 1. Resolve paths -------------------------------------------------------
.find_csv <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grep("^--file=", args)])
  here <- if (length(f)) dirname(normalizePath(f, mustWork = FALSE)) else
    tryCatch(dirname(rstudioapi::getSourceEditorContext()$path),
             error = function(e) ".")
  cand <- normalizePath(file.path(here, "..", "..", "easy_scripts_dataset.csv"),
                        mustWork = FALSE)
  if (file.exists(cand)) return(cand)
  file.path(PROJECT_ROOT, "choice R pipeline/easy_scripts/easy_scripts_dataset.csv")
}
.CSV     <- .find_csv()
.OUT_DIR <- file.path(PROJECT_ROOT, "choice R pipeline/easy_scripts/standalone/outputs")
dir.create(.OUT_DIR, showWarnings = FALSE, recursive = TRUE)
cat("Loading:", .CSV, "\n")

# ---- 2. Load + filter -------------------------------------------------------
dat <- read.csv(.CSV, stringsAsFactors = FALSE)
d <- dat %>%
  dplyr::filter(is.finite(mean_polarisation)) %>%
  dplyr::select(phys_trial_id, tank, trial_date, treatment,
                timepoint, mean_polarisation)
cat("Rows:", nrow(d), "\n")

# ---- 3. Factors -------------------------------------------------------------
# timepoint stays NUMERIC (1, 2, 3) — pipeline-selected model uses
# continuous time.
d$treatment  <- factor(d$treatment, levels = c("control", "exercise choice"))
d$trial_date <- factor(d$trial_date)
stopifnot(all(d$timepoint %in% 1:3))
cat("Treatments:", paste(levels(d$treatment), collapse = " / "),
    "  Trial dates:", nlevels(d$trial_date), "\n")
cat("N per cell:\n"); print(table(d$treatment, d$timepoint))

# ---- 4. Why sqrt was selected ----------------------------------------------
# Pipeline's saved values (for reference): sw_p_raw=0.0201, sw_p_trans=0.0301
# (sqrt helped only marginally; pipeline still picked it because Levene
# stayed clean at p>0.7 both raw and transformed).
.sw  <- function(x)    suppressWarnings(shapiro.test(x)$p.value)
.lev <- function(x, g) suppressWarnings(car::leveneTest(x ~ g)$`Pr(>F)`[1])
cat("\n[normality]  Shapiro-Wilk p:  raw =", signif(.sw(d$mean_polarisation), 4),
    "  sqrt =", signif(.sw(sqrt(d$mean_polarisation)), 4), "\n")
cat("[variance]   Levene p (by trt): raw =",
    signif(.lev(d$mean_polarisation,       d$treatment), 4),
    "  sqrt =", signif(.lev(sqrt(d$mean_polarisation), d$treatment), 4), "\n")
cat("# Pipeline picked sqrt; the SW gain is small but variance is preserved.\n")

# ---- 5. Why this random effect wins ----------------------------------------
# Pipeline's AICc table (12 candidates) has SIX candidates tied at -189.36.
# The wrapper picks (1|trial_date) by candidate-list order. Show the tie.
.fit_re <- function(re) lmer(
  as.formula(paste0("sqrt(mean_polarisation) ~ treatment * timepoint + ", re)),
  data = d, REML = TRUE)
re_cands <- list(
  "(1|trial_date)"     = "(1|trial_date)",
  "(1|phys_trial_id)"  = "(1|phys_trial_id)",
  "(1|tank)"           = "(1|tank)"
)
aicc <- sapply(re_cands, function(f) AIC(.fit_re(f)))
aicc_tbl <- data.frame(re = names(re_cands), AIC = round(aicc, 2),
                       dAIC = round(aicc - min(aicc), 3))
cat("\n[RE selection — top candidates (six-way tie in full pipeline table)]\n")
print(aicc_tbl, row.names = FALSE)
cat("# Six candidates tie at AICc = -189.36 -- choice of RE is\n",
    "# essentially arbitrary.  Pipeline picks (1|trial_date) by order.\n",
    sep = "")

# ---- 6. Fit the chosen model ------------------------------------------------
fit <- lmer(sqrt(mean_polarisation) ~ treatment * timepoint + (1 | trial_date),
            data = d, REML = TRUE)

# ---- 7. Type-III Kenward-Roger ANOVA ----------------------------------------
av <- anova(fit, type = 3, ddf = "Kenward-Roger")
cat("\n[Type-III ANOVA — Kenward-Roger]\n"); print(av)
cat("# Pipeline anova.csv: treatment F(1,37) = 4.343, p = 0.044;\n",
    "# timepoint p = 0.160; treatment:timepoint p = 0.337 (no interaction).\n",
    sep = "")

# ---- 8. emmeans + Tukey CLD (raw letters; match the saved CSV) -------------
emm <- emmeans(fit, ~ treatment * timepoint,
               at = list(timepoint = 1:3), type = "response")
cld_raw <- multcomp::cld(emm, adjust = "tukey", Letters = letters,
                          decreasing = FALSE) %>% as.data.frame()
cld_raw$.group <- trimws(cld_raw$.group)
cat("\n[Tukey CLD — RAW letters; matches",
    "STEP5_stats/<ts>/polarisation_timepoint/cld_treatmentxtimepoint.csv]\n")
print(cld_raw[, c("treatment", "timepoint", "response", "SE", ".group")],
      row.names = FALSE)

# ---- 8b. Pairwise contrasts — significant pairs only (p < 0.05) ------------
# emmeans auto-switches "tukey" -> "sidak" for a crossed treatment*timepoint
# grid (Tukey HSD is only valid for a single balanced set of means).
# Estimates are on the SQRT scale (the model's response scale), not raw units.
con <- as.data.frame(summary(emmeans::contrast(emm, method = "pairwise",
                                               adjust = "tukey")))
con_sig <- con[!is.na(con$p.value) & con$p.value < 0.05, , drop = FALSE]
cat("\n[Pairwise contrasts (Sidak-adjusted) — significant pairs only (p < 0.05)]\n")
cat("# Estimates on sqrt scale; emmeans auto-switched Tukey -> Sidak.\n")
if (nrow(con_sig) == 0) {
  cat("# No significant pairs.\n")
} else {
  show_cols <- intersect(c("contrast","estimate","SE","df","t.ratio","z.ratio","p.value"),
                         names(con_sig))
  print(con_sig[, show_cols, drop = FALSE], row.names = FALSE, digits = 4)
}

# ---- 9. CLD display remap (what Figure_10 panel B shows) -------------------
# Inline port of .mg_remap_cld() — same algorithm as the IID script.
# For polarisation the remap is NON-IDENTITY because the raw 'c' letter
# is on control/interval-1, but the display walks control first.
remap_cld <- function(df) {
  df <- df[order(match(df$treatment, c("control", "exercise choice")),
                 df$timepoint), , drop = FALSE]
  seen <- character(0)
  for (lbl in df$.group)
    for (l in strsplit(lbl, "")[[1]])
      if (!l %in% seen) seen <- c(seen, l)
  m <- setNames(letters[seq_along(seen)], seen)
  df$.display <- vapply(df$.group, function(g) {
    new <- m[strsplit(g, "")[[1]]]
    paste(sort(new[!is.na(new)]), collapse = "")
  }, character(1))
  df
}
cld_show <- remap_cld(cld_raw)
cat("\n[CLD — RAW vs FIGURE-REMAPPED letters]\n")
print(cld_show[, c("treatment", "timepoint", ".group", ".display")],
      row.names = FALSE)
cat("# Expected display mapping for polarisation:\n",
    "#   ctrl/1: c -> a   ctrl/2: bc -> ab   ctrl/3: abc -> abc\n",
    "#   ec/1: abc -> abc ec/2: a  -> c     ec/3: ab  -> bc\n",
    sep = "")
cat("# Figure_10 panel B further applies collapse_tp=2L (single label\n",
    "# above interval 2) and force_separate at interval 3 for cosmetics --\n",
    "# that's a display step, not a stats step.\n", sep = "")

# ---- 10. Observed cell means + SE (reality check) --------------------------
obs <- d %>%
  dplyr::group_by(treatment, timepoint) %>%
  dplyr::summarise(n = dplyr::n(),
                   mean = mean(mean_polarisation),
                   sem  = sd(mean_polarisation) / sqrt(dplyr::n()),
                   .groups = "drop")
cat("\n[Observed cell means + SE]\n"); print(as.data.frame(obs))

# ---- 11. Minimal plot ------------------------------------------------------
plot_df <- merge(obs, cld_show, by = c("treatment", "timepoint"))
plot_df$top <- plot_df$mean + plot_df$sem
yr <- diff(range(c(plot_df$mean - plot_df$sem, plot_df$mean + plot_df$sem)))
p <- ggplot(plot_df, aes(timepoint, mean, colour = treatment,
                          shape = treatment, group = treatment)) +
  geom_line() + geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem),
                width = 0.10, colour = "black") +
  geom_text(aes(y = top + 0.05 * yr, label = .group),
            fontface = "bold", show.legend = FALSE) +
  geom_text(aes(y = top + 0.13 * yr, label = paste0("(", .display, ")")),
            colour = "grey35", size = 3.2, show.legend = FALSE) +
  scale_x_continuous(breaks = 1:3, labels = c("5-25", "45-65", "85-105")) +
  labs(x = "Interval (min)", y = "Polarisation",
       title = "Polarisation by interval - observed means + Tukey CLD",
       subtitle = "Top label = RAW letters (CSV);  (paren) = figure remap") +
  theme_minimal(base_size = 11)
out_png <- file.path(.OUT_DIR, "mean_polarisation_by_tp.png")
ggsave(out_png, p, width = 150, height = 110, units = "mm", dpi = 200)
cat("\nSaved plot:", out_png, "\n")
