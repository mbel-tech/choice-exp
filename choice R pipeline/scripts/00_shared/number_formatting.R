# =============================================================================
# number_formatting.R -- the project's single number-display convention
# =============================================================================
# Author decision, 2026-08-07: ALL text and figures use ONE formatting rule,
# replacing the earlier per-section conventions (behaviour text at 3 decimals,
# endocrine text at 2, figures at signif(4) or fixed decimals depending on
# which engine produced them). Two functions:
#
#   fmt_F(x)   -- F-statistics and Wald chi-square statistics: FOUR SIGNIFICANT
#                 FIGURES, trailing zeros trimmed. E.g. 10.348 -> "10.35",
#                 9.724 -> "9.724", 0.4585 -> "0.4585", 3.930 -> "3.93".
#                 Because it is significant figures (not fixed decimal places),
#                 the number of digits shown after the decimal point varies
#                 with the magnitude of the statistic -- this is intentional.
#
#   fmt3(x)    -- every OTHER statistical value introduced by this analysis
#                 pass: p-values (once below the "< 0.001" threshold), partial
#                 eta-squared, partial omega-squared, Cohen's d, R^2m/R^2c,
#                 ICC. THREE decimal places, dropped to TWO only when the
#                 third decimal digit is exactly zero. E.g. 0.019 -> "0.019",
#                 0.150 -> "0.15", 0.500 -> "0.50". This is a single trailing-
#                 zero trim, not full trailing-zero stripping (0.500 -> "0.50",
#                 never "0.5").
#
# OUT OF SCOPE: descriptive means, SDs/SEMs and percentages that predate this
# analysis pass (concentrations, masses, cm/s values, etc.) keep the precision
# already established in the manuscript for that quantity -- reformatting them
# was not requested and would be an unrelated change.
#
# This file is a REFERENCE copy for anyone auditing the convention. Because
# every analysis engine in this project is written to run standalone
# (see each script's own header), the two functions below are also DUPLICATED
# verbatim, not sourced, into every engine that renders a figure caption:
#   - choice R pipeline/scripts/01_pipeline_analysis/
#       activity_analysis_STATS_choice_exp.R (+ _REVISED.R)
#       behavioural_sequence_analysis_choice_exp.R
#   - choice exp for claude mono and cortisol/claude output/Scripts/
#       analysis_b3_REVISED.R
#   - choice R pipeline/scripts/00_shared/
#       rebuild_figure_caption_table.R
#       effect_size_and_sample_size.R
# If this rule ever changes again, all seven copies must be updated together.
# =============================================================================

fmt_F <- function(x) {
  if (!is.finite(x)) return("NA")
  s <- format(signif(x, 4), scientific = FALSE, trim = TRUE)
  # Only trim trailing zeros when a decimal point is present -- otherwise an
  # integer denominator df like 30 is mangled to "3" (see project number-
  # formatting defect log).
  if (grepl(".", s, fixed = TRUE)) {
    s <- sub("0+$", "", s)
    s <- sub("\\.$", "", s)
  }
  s
}

fmt3 <- function(x) {
  if (!is.finite(x)) return("NA")
  s <- sprintf("%.3f", x)
  if (grepl("0$", s)) s <- substr(s, 1, nchar(s) - 1)
  s
}
