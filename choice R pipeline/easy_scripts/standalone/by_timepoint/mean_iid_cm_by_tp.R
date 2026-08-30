# =============================================================================
# Standalone easy-script: mean inter-individual distance (IID, cm) by interval
# -----------------------------------------------------------------------------
# Self-contained reproduction of the Pipeline-A model for res_iid_tp.
# Source no project file. Uses only CRAN packages.
#
# Why this exists:
#   The auto-generated counterpart (../by_timepoint/mean_iid_cm_by_tp.R)
#   sources _helpers.R, which pulls in the full pipeline machinery
#   (run_lmm_analysis + figure helpers).  That's faithful, but opaque.
#   This standalone version makes EVERY step visible so the reader can see
#   why the figure, the CSV and the post-hocs don't trivially line up.
#
# Three layers sit between the raw data and Figure_10's letters:
#   1) response transform = sqrt
#   2) timepoint is treated as CONTINUOUS in the chosen model (df_num=1
#      for the interaction; emmeans evaluated at t in {1,2,3} on a line)
#   3) display CLD letters are remapped by precedence (control before
#      exercise choice; interval 1 -> 2 -> 3)
#
# Pipeline source: scripts/01_pipeline_analysis/activity_analysis_STATS_choice_exp.R
#   model fit: line 2967  (res_iid_tp <- run_lmm_analysis(...))
#   AICc table: STEP5_stats/STEP5_stats_20260511_175522/iid_timepoint/aicc_selection.csv
#   normality: STEP5_stats/STEP5_stats_20260511_175522/iid_timepoint/normality_check.csv
#
# Hard-coded pipeline choices (from the May-11 outputs above):
#   transform     = sqrt          (best of {identity, log1p, sqrt})
#   fixed formula = sqrt(y) ~ treatment * timepoint        # continuous t
#   random effect = (1 | phys_trial_id)
#   selected AICc = 176.31         (next-best 179.06, +2.75 over winner)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(lme4); library(lmerTest); library(pbkrtest)
  library(emmeans); library(multcomp); library(ggplot2); library(car)
})

# ---- 1. Resolve paths -------------------------------------------------------
# Works whether sourced in RStudio or run via Rscript.  Fallback to the
# known absolute path so a tweaked working dir doesn't break the script.
.find_csv <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grep("^--file=", args)])
  here <- if (length(f)) dirname(normalizePath(f, mustWork = FALSE)) else
    tryCatch(dirname(rstudioapi::getSourceEditorContext()$path),
             error = function(e) ".")
  cand <- normalizePath(file.path(here, "..", "..", "easy_scripts_dataset.csv"),
                        mustWork = FALSE)
  if (file.exists(cand)) return(list(path = cand, kind = "internal"))
  # The Zenodo deposit, if it has been unpacked. Preferred when present, so
  # that this script exercises the same data a reader downloads.
  dep <- file.path(Sys.getenv("CHOICE_EXP_DATA_ROOT",
                              unset = file.path(PROJECT_ROOT, "zenodo_dataset")),
                   "pref_trials.csv")
  if (file.exists(dep)) return(list(path = dep, kind = "deposit"))
  list(path = file.path(PROJECT_ROOT,
         "choice R pipeline/easy_scripts/easy_scripts_dataset.csv"),
       kind = "internal")
}
.SRC     <- .find_csv()
.CSV     <- .SRC$path
.OUT_DIR <- file.path(PROJECT_ROOT, "choice R pipeline/easy_scripts/standalone/outputs")
dir.create(.OUT_DIR, showWarnings = FALSE, recursive = TRUE)
cat("Loading:", .CSV, "\n")

# ---- 2. Load + filter -------------------------------------------------------
dat <- read.csv(.CSV, stringsAsFactors = FALSE)
if (identical(.SRC$kind, "deposit")) {
  names(dat)[names(dat) == "phys_trial"]           <- "phys_trial_id"
  names(dat)[names(dat) == "interval"]             <- "timepoint"
  names(dat)[names(dat) == "mean_school_area_cm2"] <- "mean_hull_area_cm2"
  dat$trial_date <- format(as.Date(dat$trial_date, "%Y-%m-%d"), "%d.%m.%Y")
}
d <- dat %>%
  dplyr::filter(is.finite(mean_iid_cm)) %>%
  dplyr::select(phys_trial_id, tank, trial_date, treatment,
                timepoint, mean_iid_cm)
cat("Rows:", nrow(d), "\n")

# ---- 3. Factors -------------------------------------------------------------
# treatment as ordered factor so control is the reference;
# timepoint stays NUMERIC (1, 2, 3) — the pipeline-selected model uses
# continuous time, NOT the 3-level factor `timepoint_f`.
d$treatment <- factor(d$treatment, levels = c("control", "exercise choice"))
stopifnot(all(d$timepoint %in% 1:3))
cat("Treatments:", paste(levels(d$treatment), collapse = " / "),
    "  N per cell:\n"); print(table(d$treatment, d$timepoint))

# ---- 4. Why sqrt was selected ----------------------------------------------
# Reproduce the SW/Levene comparison stored in normality_check.csv.
# Pipeline's saved values (for reference): sw_p_raw=0.0197, sw_p_trans=0.1317.
.sw <- function(x) suppressWarnings(shapiro.test(x)$p.value)
.lev <- function(x, g) suppressWarnings(
  car::leveneTest(x ~ g)$`Pr(>F)`[1])
cat("\n[normality]  Shapiro-Wilk p:  raw =", signif(.sw(d$mean_iid_cm), 4),
    "  sqrt =", signif(.sw(sqrt(d$mean_iid_cm)), 4), "\n")
cat("[variance]   Levene p (by trt): raw =",
    signif(.lev(d$mean_iid_cm,       d$treatment), 4),
    "  sqrt =", signif(.lev(sqrt(d$mean_iid_cm), d$treatment), 4), "\n")
cat("# Pipeline picked sqrt because SW improves from 0.020 -> 0.132.\n")

# ---- 5. Why this random effect wins ----------------------------------------
# Fit the top-3 RE candidates from the pipeline's AICc table and show the
# winner beats the others.  Full table (all 12 candidates) is at
#   STEP5_stats/STEP5_stats_20260511_175522/iid_timepoint/aicc_selection.csv
.fit_re <- function(re) lmer(
  as.formula(paste0("sqrt(mean_iid_cm) ~ treatment * timepoint + ", re)),
  data = d, REML = TRUE)
re_cands <- list(
  "(1|phys_trial_id)"                  = "(1|phys_trial_id)",
  "(1|phys_trial_id) + (1|tank)"       = "(1|phys_trial_id) + (1|tank)",
  "(1|phys_trial_id) + (1|trial_date)" = "(1|phys_trial_id) + (1|trial_date)"
)
aicc <- sapply(re_cands, function(f) AIC(.fit_re(f)))
aicc_tbl <- data.frame(re = names(re_cands), AIC = round(aicc, 2),
                       dAIC = round(aicc - min(aicc), 2))
cat("\n[RE selection — top 3 candidates from pipeline AICc table]\n")
print(aicc_tbl, row.names = FALSE)
cat("# Pipeline picked (1|phys_trial_id), AICc = 176.31",
    " (next-best +2.75).\n")

# ---- 6. Fit the chosen model ------------------------------------------------
fit <- lmer(sqrt(mean_iid_cm) ~ treatment * timepoint + (1 | phys_trial_id),
            data = d, REML = TRUE)

# ---- 7. Type-III Kenward-Roger ANOVA ----------------------------------------
# treatment:timepoint has df_num = 1 because timepoint is continuous, so
# the interaction is the SLOPE difference between treatments (one scalar),
# not the {tp2 vs tp1, tp3 vs tp1}-pair test you'd get from timepoint_f.
av <- anova(fit, type = 3, ddf = "Kenward-Roger")
cat("\n[Type-III ANOVA — Kenward-Roger]\n"); print(av)
cat("# Pipeline anova.csv: treatment p=0.740, timepoint p=0.069,\n",
    "# treatment:timepoint F(1,30) = 6.096, p = 0.019.\n", sep = "")

# ---- 8. emmeans + Tukey CLD (raw letters; match the saved CSV) -------------
# Continuous-time model -> evaluate predicted means at integer t = 1, 2, 3.
# type='response' back-transforms from sqrt scale.
emm <- emmeans(fit, ~ treatment * timepoint,
               at = list(timepoint = 1:3), type = "response")
cld_raw <- multcomp::cld(emm, adjust = "tukey", Letters = letters,
                          decreasing = FALSE) %>% as.data.frame()
cld_raw$.group <- trimws(cld_raw$.group)
cat("\n[Tukey CLD — RAW letters; this matches",
    "STEP5_stats/<ts>/iid_timepoint/cld_treatmentxtimepoint.csv]\n")
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

# ---- 9. CLD display remap (what Figure_10 shows) ---------------------------
# Inline port of .mg_remap_cld() from
#   scripts/01_pipeline_analysis/activity_analysis_STATS_choice_exp.R (~L6056).
# Walk cells in precedence order (control then ec; interval 1 -> 2 -> 3),
# build a {raw letter -> new letter} map in first-seen order, re-sort each
# label alphabetically.
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
cat("# For IID this remap is identity (raw letters already ordered).\n")

# ---- 10. Observed cell means + SE (reality check) --------------------------
# Compare these to Step 8's predicted emmeans: differences = the linear-fit
# 'smoothing' effect.  This is the single biggest reason figure dots can
# look 'inconsistent' with the CLD letters.
obs <- d %>%
  dplyr::group_by(treatment, timepoint) %>%
  dplyr::summarise(n = dplyr::n(),
                   mean = mean(mean_iid_cm),
                   sem  = sd(mean_iid_cm) / sqrt(dplyr::n()),
                   .groups = "drop")
cat("\n[Observed cell means + SE]\n"); print(as.data.frame(obs))

# ---- 11. Minimal plot (observed means, raw + remapped CLD overlaid) --------
# Plain ggplot.  Letters are placed above each error-bar tip.
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
  labs(x = "Interval (min)", y = "Mean IID (cm)",
       title = "IID by interval - observed means + Tukey CLD",
       subtitle = "Top label = RAW letters (CSV);  (paren) = figure remap") +
  theme_minimal(base_size = 11)
out_png <- file.path(.OUT_DIR, "mean_iid_cm_by_tp.png")
ggsave(out_png, p, width = 150, height = 110, units = "mm", dpi = 200)
cat("\nSaved plot:", out_png, "\n")
