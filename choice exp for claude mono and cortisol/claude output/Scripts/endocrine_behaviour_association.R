# =============================================================================
# M5 (R1-4): Trial-level association between behavioural flow-preference and
# trial-aggregated brain monoamines — NOT a directional/causal analysis.
#
# WHY THIS SCRIPT
#   Neurochemistry (analysis_b3.R) is sampled terminally, once per fish, at the
#   end of a trial. Behaviour (activity_analysis_STATS_choice_exp.R) is a
#   whole-session, group-level measure. The manuscript must not imply that
#   monoamine levels "underlie" or "drive" the behavioural preference without
#   evidence. This script checks whether a trial-level association even exists
#   AT ALL, and reports it (r, p, q) with the sampling-design limitation
#   stated explicitly — it cannot and does not establish direction or causality.
#
# LINKAGE
#   Behaviour and neurochemistry share the same 16 physical trials (same
#   tank / testing date / treatment triple; confirmed by direct inspection of
#   trial_summary_choice_exp.xlsx vs the B3 cortisol/monoamine sheets — the
#   raw `trial` INTEGER differs in numbering between the two datasets, so the
#   robust join key is (tank, date, treatment), not the raw trial id.
#
# INPUTS
#   Behaviour : latest choice R pipeline/output/STEP5_stats/*/analysis_ready.csv
#   Endocrine : choice exp for claude mono and cortisol/claude output/Data/monoamine_long.csv
#               (written by analysis_b3.R; run that script first)
#
# OUTPUT
#   choice exp for claude mono and cortisol/claude output/Data/
#     endocrine_behaviour_association.csv   (r, p, q per analyte x area)
#     endocrine_behaviour_trial_table.csv   (the joined trial-level table used)
#   choice exp for claude mono and cortisol/claude output/Figures/
#     assoc_DM_DA_vs_flow.png, assoc_DM_5HT_vs_flow.png  (headline scatter pair)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(ggplot2); library(purrr)
})
set.seed(20260706)

ROOT_ENDO <- file.path(PROJECT_ROOT, "choice exp for claude mono and cortisol")
OUT_ENDO  <- file.path(ROOT_ENDO, "claude output")
DATA_DIR  <- file.path(OUT_ENDO, "Data")
FIG_DIR   <- file.path(OUT_ENDO, "Figures")

fmt_p <- function(p) {
  vapply(p, function(pi) {
    if (is.na(pi)) return("NA")
    base <- if (pi < 0.0001) "< 0.0001" else format(signif(pi, 3), scientific = FALSE)
    star <- if (pi < 0.001) " ***" else if (pi < 0.01) " **" else if (pi < 0.05) " *" else ""
    paste0(base, star)
  }, character(1))
}

# ---- 1. Behaviour: latest STEP5_stats analysis_ready.csv, aggregated to trial level ----
.find_latest_analysis_ready <- function() {
  parent <- file.path(PROJECT_ROOT, "choice R pipeline/output/STEP5_stats")
  subs <- list.dirs(parent, full.names = TRUE, recursive = FALSE)
  subs <- subs[file.exists(file.path(subs, "analysis_ready.csv"))]
  if (!length(subs)) stop("No analysis_ready.csv found under ", parent)
  subs[which.max(file.mtime(subs))]
}
beh_dir  <- .find_latest_analysis_ready()
beh_path <- file.path(beh_dir, "analysis_ready.csv")
message("Behaviour source: ", beh_path)
beh <- read_csv(beh_path, show_col_types = FALSE)

beh_trial <- beh %>%
  mutate(date_key = as.Date(substr(trial_date, 1, 10)),
         tank_key = as.character(tank),
         # Harmonise to the endocrine sheet's vocabulary ("control"/"treat").
         treat_key = dplyr::recode(treatment, "exercise choice" = "treat", "control" = "control")) %>%
  group_by(phys_trial_id, tank_key, date_key, treat_key) %>%
  summarise(
    prop_flow_trial   = mean(prop_time_in_flow, na.rm = TRUE),
    n_timepoints      = n(),
    treatment_behaviour = dplyr::first(treatment),
    .groups = "drop"
  )
message(sprintf("Behaviour: %d physical trials aggregated.", nrow(beh_trial)))

# ---- 2. Endocrine: monoamine_long.csv, aggregated to trial level ----------
mono_path <- file.path(DATA_DIR, "monoamine_long.csv")
if (!file.exists(mono_path))
  stop("monoamine_long.csv not found — run analysis_b3.R first: ", mono_path)
mono_long <- read_csv(mono_path, show_col_types = FALSE)

mono_trial <- mono_long %>%
  filter(is.finite(value), value > 0) %>%
  mutate(date_key = as.Date(date), tank_key = as.character(tank),
         treat_key = as.character(treatment)) %>%
  group_by(tank_key, date_key, treat_key, area, analyte) %>%
  summarise(mono_mean = mean(value, na.rm = TRUE), n_fish = n(), .groups = "drop")

# ---- 3. Join on (tank, date, treatment) — each tank tests BOTH treatments ----
# on the same day, so (tank, date) alone is not a unique key; treatment must
# be part of the join (verified: every tank x date combination hosts exactly
# one control and one exercise-choice trial).
joined <- mono_trial %>%
  inner_join(beh_trial, by = c("tank_key", "date_key", "treat_key"),
            relationship = "many-to-one")

n_trials_linked <- dplyr::n_distinct(joined$phys_trial_id)
message(sprintf("Linked %d of %d behavioural trials to endocrine samples via (tank, date).",
                n_trials_linked, nrow(beh_trial)))

write_csv(joined %>% select(phys_trial_id, tank = tank_key, date = date_key,
                            treatment = treat_key, treatment_behaviour,
                            area, analyte, mono_mean, n_fish, prop_flow_trial, n_timepoints),
          file.path(DATA_DIR, "endocrine_behaviour_trial_table.csv"))

# ---- 4. Feasibility gate ---------------------------------------------------
MIN_TRIALS_FOR_ASSOC <- 5
feasible <- n_trials_linked >= MIN_TRIALS_FOR_ASSOC

if (!feasible) {
  msg <- sprintf(paste0(
    "M5 ASSOCIATION NOT FEASIBLE: only %d physical trials link behaviour to ",
    "endocrine samples via (tank, date) — below the %d-trial floor for a ",
    "meaningful correlation. Reporting explicit non-feasibility; directional ",
    "wording linking neurochemistry to behavioural preference must be stripped ",
    "from Results/Discussion/Abstract regardless of this script's other output."),
    n_trials_linked, MIN_TRIALS_FOR_ASSOC)
  message(msg)
  writeLines(msg, file.path(DATA_DIR, "endocrine_behaviour_association_LIMITATION.txt"))
  quit(save = "no", status = 0)
}

# ---- 5. Per-analyte x area trial-level correlation ------------------------
assoc_one <- function(df) {
  d <- df %>% filter(is.finite(mono_mean), is.finite(prop_flow_trial))
  if (nrow(d) < MIN_TRIALS_FOR_ASSOC) return(NULL)
  ct_pearson  <- suppressWarnings(cor.test(d$mono_mean, d$prop_flow_trial, method = "pearson"))
  ct_spearman <- suppressWarnings(cor.test(d$mono_mean, d$prop_flow_trial, method = "spearman"))
  data.frame(
    n         = nrow(d),
    r_pearson = unname(ct_pearson$estimate),
    p_pearson = ct_pearson$p.value,
    rho_spearman = unname(ct_spearman$estimate),
    p_spearman   = ct_spearman$p.value
  )
}

assoc_results <- joined %>%
  group_by(area, analyte) %>%
  group_modify(~ assoc_one(.x)) %>%
  ungroup() %>%
  filter(!is.na(n))

# BH-FDR per brain region — mirrors the M1 family definition used in analysis_b3.R.
assoc_results <- assoc_results %>%
  group_by(area) %>%
  mutate(q_pearson = p.adjust(p_pearson, method = "BH")) %>%
  ungroup() %>%
  mutate(sig_raw = fmt_p(p_pearson), sig_BH = fmt_p(q_pearson)) %>%
  arrange(area, analyte)

write_csv(assoc_results, file.path(DATA_DIR, "endocrine_behaviour_association.csv"))
message(sprintf("Association results: %d analyte x area tests (n trials range %d-%d).",
                nrow(assoc_results), min(assoc_results$n), max(assoc_results$n)))

# ---- 6. Headline pair: DM DA and DM 5-HT (the manuscript's primary claims) ----
make_assoc_plot <- function(area_sel, analyte_sel, fname) {
  d <- joined %>% filter(area == area_sel, analyte == analyte_sel,
                        is.finite(mono_mean), is.finite(prop_flow_trial))
  if (nrow(d) < MIN_TRIALS_FOR_ASSOC) return(invisible(NULL))
  row <- assoc_results %>% filter(area == area_sel, analyte == analyte_sel)
  cap <- if (nrow(row))
    sprintf("Pearson r = %.2f, p = %.3f, q = %.3f (n = %d trials)",
            row$r_pearson[1], row$p_pearson[1], row$q_pearson[1], row$n[1])
  else ""
  p <- ggplot(d, aes(x = prop_flow_trial, y = mono_mean)) +
    geom_point(aes(colour = treat_key), size = 2.5) +
    geom_smooth(method = "lm", se = TRUE, colour = "black", linewidth = 0.6) +
    labs(title = sprintf("%s %s vs. trial flow-preference", area_sel, analyte_sel),
         x = "Trial-mean proportion of time in flow zone (behaviour)",
         y = sprintf("Trial-mean [%s] in %s (endocrine)", analyte_sel, area_sel),
         colour = "Treatment", caption = cap) +
    theme_bw()
  ggsave(file.path(FIG_DIR, fname), p, width = 6, height = 4.5, dpi = 200)
}
make_assoc_plot("DM", "DA",   "assoc_DM_DA_vs_flow.png")
make_assoc_plot("DM", "5-HT", "assoc_DM_5HT_vs_flow.png")

# ---- 7. Explicit limitation statement (always written) -------------------
limitation_txt <- paste0(
  "M5 (R1-4) trial-level association — limitation statement.\n\n",
  sprintf("n = %d physical trials had both behavioural and endocrine data, linked via ", n_trials_linked),
  "(tank, testing date) rather than a shared trial-identifier scheme (the raw `trial` ",
  "integer differs in numbering between the behaviour and endocrine worksheets).\n\n",
  "This is a TRIAL-LEVEL, CROSS-SECTIONAL association only. Monoamines were sampled ",
  "terminally from individual fish at a single timepoint (end of session); behaviour ",
  "is a whole-session, group-level aggregate. The design cannot distinguish direction ",
  "(behaviour -> neurochemistry, neurochemistry -> behaviour, or a shared upstream cause), ",
  "and the small number of trials (n = ", n_trials_linked, ") limits power to detect all ",
  "but large effects. Any reported correlation is evidence of covariation only; ",
  "directional language ('underlies', 'drives', 'causes') is not supported by this ",
  "analysis and should not be used in Results/Discussion/Abstract for these effects."
)
writeLines(limitation_txt, file.path(DATA_DIR, "endocrine_behaviour_association_LIMITATION.txt"))
cat(limitation_txt, "\n")
