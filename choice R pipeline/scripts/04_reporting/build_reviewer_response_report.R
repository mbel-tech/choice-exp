# =============================================================================
# build_reviewer_response_report.R
# -----------------------------------------------------------------------------
# Consolidates the current headline results (behaviour, cortisol, monoamines)
# together with the M1-M5 statistical-revision evidence into ONE reviewer/
# editor response report, in two formats:
#   - Reviewer_Response_Report.md   (plain structured markdown, easy for an
#                                    LLM/Claude to read back in a later session)
#   - Reviewer_Response_Report.docx (formatted Word doc for humans)
#
# Reads only pre-computed CSV/txt outputs already sitting on disk (produced by
# activity_analysis_STATS_choice_exp.R and analysis_b3.R on 2026-07-06, the
# M1-M5 statistical-revision pass) -- it does not refit any model itself.
# Re-run those two engines first if you need fresher numbers.
#
# Covers every reviewer/editor code from the annotated manuscript EXCEPT
# R1-1 (rheotaxis-confound framing) and R2-1 (split into two manuscripts),
# excluded on explicit user instruction -- both are prose/strategy decisions,
# not something a script can resolve.
#
# USAGE:
#   Rscript file.path(PROJECT_ROOT, "choice R pipeline/scripts/04_reporting/build_reviewer_response_report.R")
#
# OUTPUT:
#   D:/CHOICE R SCRIPTS/Reviewer_Response_Report.md    (overwritten)
#   D:/CHOICE R SCRIPTS/Reviewer_Response_Report.docx  (overwritten)
# =============================================================================

suppressPackageStartupMessages({
  library(officer)
  library(flextable)
})

# ---- 0. Paths ---------------------------------------------------------------

BEHAV_DIR <- file.path(PROJECT_ROOT, "choice R pipeline/STEP5_stats/STEP5_stats_20260511_175522")
ENDO_DIR  <- file.path(PROJECT_ROOT, "choice exp for claude mono and cortisol/claude output/Data")
OUT_MD    <- file.path(PROJECT_ROOT, "Reviewer_Response_Report.md")
OUT_DOCX  <- file.path(PROJECT_ROOT, "Reviewer_Response_Report.docx")

stopifnot(
  "Behaviour output dir not found -- re-run run_manu_graphs_only.R first" = dir.exists(BEHAV_DIR),
  "Endocrine output dir not found -- re-run analysis_b3.R first"          = dir.exists(ENDO_DIR)
)

# ---- 1. Helpers --------------------------------------------------------------

read_csv_safe <- function(path, colClasses = NA) {
  if (!file.exists(path)) {
    warning("[build_reviewer_response_report] Missing file: ", path)
    return(NULL)
  }
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, colClasses = colClasses)
}

round_df <- function(df, digits = 3) {
  if (is.null(df)) return(df)
  df[] <- lapply(df, function(col) if (is.numeric(col)) signif(col, digits) else col)
  df
}

na_to_dash <- function(df) {
  if (is.null(df)) return(df)
  df[] <- lapply(df, function(col) {
    col <- as.character(col)
    col[is.na(col) | col == "NA"] <- "\u2014"
    col
  })
  df
}

# Effect sizes derived directly from an F-test's own F/df1/df2 -- no refit
# needed, valid only for F-tests (not the Wald-chisq beta-GLMM rows, where
# R2m is the reported variance-explained metric instead).
partial_eta_sq <- function(F, df1, df2) (F * df1) / (F * df1 + df2)

partial_omega_sq <- function(F, df1, df2) {
  num <- df1 * (F - 1)
  den <- F * df1 + df2 + 1
  pmax(num / den, 0)  # clamp negative estimates (near-null effects) to 0
}

# ---- Sensitivity (minimum-detectable-effect) helpers -------------------------
# Retrospective sample-size justification for an already-collected dataset.
# NOT post-hoc/observed power (which is a monotone function of p and carries no
# information); instead: hold N, alpha and target power fixed and solve for the
# smallest effect the design could have detected.
#
# N is taken as df1 + df2 + 1 from the model's OWN Kenward-Roger denominator df,
# so the experimental unit (trial for behaviour, fish for endocrine) and the
# random-effect structure are inherited automatically -- no separate ICC or
# design-effect assumption is introduced. Verified against power.t.test().

lambda_for_power <- function(df1, df2, alpha, power) {
  if (is.na(df1) || is.na(df2)) return(NA_real_)
  Fcrit <- qf(1 - alpha, df1, df2)
  tryCatch(
    uniroot(function(lam) pf(Fcrit, df1, df2, ncp = lam, lower.tail = FALSE) - power,
            c(1e-6, 1000), tol = 1e-9)$root,
    error = function(e) NA_real_)
}

# Minimum detectable partial eta-sq. G*Power "as in Cohen" convention:
# lambda = f2 * N and eta2p = f2/(1+f2)  =>  eta2p = lambda/(N + lambda).
mdes_eta2 <- function(df1, df2, alpha = 0.05, power = 0.80) {
  lam <- lambda_for_power(df1, df2, alpha, power)
  if (is.na(lam)) return(NA_real_)
  lam / (df1 + df2 + 1 + lam)
}

# Same, expressed as a two-group standardised mean difference (d = 2*sqrt(f2)).
mdes_d <- function(df1, df2, alpha = 0.05, power = 0.80) {
  lam <- lambda_for_power(df1, df2, alpha, power)
  if (is.na(lam)) return(NA_real_)
  2 * sqrt(lam / (df1 + df2 + 1))
}

# Minimum detectable |r|, solved on the exact noncentral t (the Fisher-z
# approximation runs ~5 percentage points conservative at n = 15).
mdes_r <- function(n, alpha = 0.05, power = 0.80) {
  tryCatch(
    uniroot(function(r) {
      ncp <- r * sqrt(n - 2) / sqrt(1 - r^2); tc <- qt(1 - alpha/2, n - 2)
      pt(tc, n - 2, ncp = ncp, lower.tail = FALSE) + pt(-tc, n - 2, ncp = ncp) - power
    }, c(0.01, 0.999))$root,
    error = function(e) NA_real_)
}

# Power actually available to detect an observed correlation (this IS a
# legitimate use: it quantifies how uninformative a given null result is).
power_at_r <- function(r, n, alpha = 0.05) {
  ncp <- abs(r) * sqrt(n - 2) / sqrt(1 - r^2); tc <- qt(1 - alpha/2, n - 2)
  pt(tc, n - 2, ncp = ncp, lower.tail = FALSE) + pt(-tc, n - 2, ncp = ncp)
}

# Raw CSV column name -> human-readable label. Applied to every table in both
# deliverables so headers never show a raw R/CSV identifier; full definitions
# of each term live in the Glossary section at the top of the report.
COLUMN_LABELS <- c(
  "Analysis"        = "Analysis / model",
  "p_raw"           = "Raw p",
  "p_BH"            = "BH-adjusted q",
  "sig_BH"          = "Significant (q<0.05)",
  "sig_raw"         = "Significant (raw p)",
  "stat_type"       = "Test type",
  "df1"             = "df (numerator)",
  "df2"             = "df (denominator)",
  "model"           = "Model formula",
  "term"            = "Term tested",
  "p_value"         = "p-value",
  "Trial_level_p"   = "p (trial-level refit)",
  "Tank_level_p"    = "p (tank-level refit)",
  "Trial_beta"      = "β (trial-level)",
  "Tank_beta"       = "β (tank-level)",
  "Beta_ratio"      = "β ratio (trial/tank)",
  "Trial_sig"       = "Significant (trial level)",
  "Tank_sig"        = "Significant (tank level)",
  "Direction_agree" = "Trial & tank agree on direction",
  "outcome"         = "Outcome variable",
  "F"               = "F statistic",
  "Df"              = "df (numerator)",
  "Df.res"          = "df (residual)",
  "Pr(>F)"          = "p-value",
  "Sum Sq"          = "Sum of squares",
  "F value"         = "F statistic (Wald)",
  "spec"            = "Random-effect specification",
  "beta"            = "β (effect estimate)",
  "SE"              = "Standard error",
  "p"               = "p-value",
  "dSE_pct"         = "ΔSE vs winner (%)",
  "dBeta"           = "Δβ vs winner",
  "flag"            = ">30% shift flag",
  "area"            = "Brain region",
  "analyte"         = "Analyte",
  "n"               = "Sample size (n)",
  "re_spec"         = "Random-effect structure",
  "r_pearson"       = "Pearson r",
  "p_pearson"       = "p (Pearson)",
  "rho_spearman"    = "Spearman ρ",
  "p_spearman"      = "p (Spearman)",
  "q_pearson"       = "BH-adjusted q (Pearson)",
  "step"            = "Decision #",
  "decision"        = "Decision",
  "basis"           = "Basis",
  "check"           = "Check",
  "value"           = "Value",
  "stat_val"        = "Test statistic (F or χ²)",
  "F_stat"          = "F statistic",
  "eta2_partial"    = "Partial η²",
  "omega2_partial"  = "Partial ω²",
  "R2m"             = "R²m (marginal)",
  "R2_marginal"     = "R²m (marginal)",
  "R2_conditional"  = "R²c (conditional)",
  "g_ci"            = "Hedges' g [95% CI]",
  "Stratum"         = "Stratum",
  "N_eff"           = "Effective N",
  "eta2_obs"        = "Observed partial η²",
  "mdes_eta2_05"    = "Min. detectable η² (α=.05)",
  "mdes_d_05"       = "Min. detectable d (α=.05)",
  "mdes_eta2_bh"    = "Min. detectable η² (strictest BH α)",
  "verdict"         = "Observed vs detection floor",
  "mdes_r_05"       = "Min. detectable |r| (α=.05)",
  "mdes_r_bh"       = "Min. detectable |r| (strictest BH α)",
  "power_at_obs"    = "Power at the observed r"
)

pretty_cols <- function(df) {
  if (is.null(df)) return(df)
  nm <- names(df)
  matched <- nm %in% names(COLUMN_LABELS)
  nm[matched] <- COLUMN_LABELS[nm[matched]]
  names(df) <- nm
  df
}

# One markdown pipe-table per data.frame; used for the .md deliverable.
# Literal "|" must be escaped or it is parsed as a cell delimiter and silently
# corrupts the table -- this bites the random-effect columns, which are full of
# lme4 "(1|trial)" notation, and any header using |r| for absolute correlation.
md_escape <- function(x) gsub("|", "\\|", as.character(x), fixed = TRUE)

md_table <- function(df, digits = 3) {
  if (is.null(df) || nrow(df) == 0) return("*(no data)*\n")
  df <- pretty_cols(na_to_dash(round_df(df, digits)))
  df[] <- lapply(df, md_escape)
  hdr <- paste0("| ", paste(md_escape(names(df)), collapse = " | "), " |")
  sep <- paste0("|", paste(rep(":---", ncol(df)), collapse = "|"), "|")
  rows <- apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
  paste(c(hdr, sep, rows), collapse = "\n")
}

# Matching flextable for the .docx deliverable.
ft_table <- function(df, digits = 3, fontsize = 8) {
  df <- pretty_cols(na_to_dash(round_df(df, digits)))
  flextable(df) |>
    set_table_properties(layout = "autofit", width = 1) |>
    bold(part = "header") |>
    bg(part = "header", bg = "#D9D9D9") |>
    fontsize(part = "all", size = fontsize) |>
    padding(part = "all", padding = 2) |>
    border_outer(part = "all", border = fp_border(color = "#888888", width = 0.5)) |>
    border_inner(part = "all", border = fp_border(color = "#CCCCCC", width = 0.25))
}

# ---- 2. Load source data -----------------------------------------------------

behav <- list(
  primary_nonzone      = read_csv_safe(file.path(BEHAV_DIR, "bh_primary_nonzone.csv")),
  jacobs_beta_primary  = read_csv_safe(file.path(BEHAV_DIR, "bh_correction_summary.csv")),
  supplementary        = read_csv_safe(file.path(BEHAV_DIR, "bh_supplementary.csv")),
  jacobs_sensitivity   = read_csv_safe(file.path(BEHAV_DIR, "bh_jacobs_sensitivity.csv")),
  exploratory          = read_csv_safe(file.path(BEHAV_DIR, "bh_exploratory.csv")),
  order_effect         = read_csv_safe(file.path(BEHAV_DIR, "order_effect_summary.csv")),
  pseudoreplication    = read_csv_safe(file.path(BEHAV_DIR, "pseudoreplication_check.csv")),
  design_balance       = read_csv_safe(file.path(BEHAV_DIR, "design_balance.csv")),
  hedges_g             = read_csv_safe(file.path(BEHAV_DIR, "hedges_g_summary.csv"))
)

endo <- list(
  cell_cross_bh        = read_csv_safe(file.path(ENDO_DIR, "cell_cross_summary_bh.csv")),
  cortisol_anova       = read_csv_safe(file.path(ENDO_DIR, "anova_kr_cortisol.csv")),
  pseudorep_cortisol   = read_csv_safe(file.path(ENDO_DIR, "pseudorep_check.csv")),
  pseudorep_mono       = read_csv_safe(file.path(ENDO_DIR, "pseudorep_check_monoamines.csv")),
  decision_log         = read_csv_safe(file.path(ENDO_DIR, "decision_log.csv"),
                                        colClasses = c(step = "character")),
  order_effect         = read_csv_safe(file.path(ENDO_DIR, "order_effect_summary.csv")),
  assoc                = read_csv_safe(file.path(ENDO_DIR, "endocrine_behaviour_association.csv")),
  cell_r2              = read_csv_safe(file.path(ENDO_DIR, "cell_r2.csv"))
)

assoc_limitation_path <- file.path(ENDO_DIR, "endocrine_behaviour_association_LIMITATION.txt")
assoc_limitation <- if (file.exists(assoc_limitation_path)) {
  paste(readLines(assoc_limitation_path, warn = FALSE), collapse = "\n")
} else {
  "(limitation file not found)"
}

# Derived views ----------------------------------------------------------------

behav_headline <- rbind(behav$primary_nonzone, behav$jacobs_beta_primary)

# Partial eta-sq / omega-sq are only defined for an F-test; the Jacobs' D
# beta-GLMM rows use a Wald chi-square instead, so those stay NA and are
# read via R2m (already computed for every row, F-KR or Wald-chisq alike).
.is_F <- behav_headline$stat_type == "F-KR"
behav_headline$eta2_partial   <- ifelse(.is_F, partial_eta_sq(behav_headline$stat_val, behav_headline$df1, behav_headline$df2), NA_real_)
behav_headline$omega2_partial <- ifelse(.is_F, partial_omega_sq(behav_headline$stat_val, behav_headline$df1, behav_headline$df2), NA_real_)

# Hedges' g (bias-corrected standardised mean difference) -- only computed
# upstream for the Gaussian (F-KR) primary-family models; join by exact
# Analysis-label match, left as NA for the Jacobs' beta-GLMM rows.
.hg_idx <- match(behav_headline$Analysis, behav$hedges_g$Analysis)
behav_headline$g_ci <- ifelse(
  !is.na(.hg_idx) & !is.na(behav$hedges_g$g[.hg_idx]),
  sprintf("%.2f [%.2f, %.2f]", behav$hedges_g$g[.hg_idx],
          behav$hedges_g$ci_lo[.hg_idx], behav$hedges_g$ci_hi[.hg_idx]),
  NA_character_
)

behav_headline <- behav_headline[, c("Analysis", "p_raw", "p_BH", "sig_BH", "stat_type",
                                      "stat_val", "df1", "df2", "R2m",
                                      "eta2_partial", "omega2_partial", "g_ci")]

behav_order_interaction <- behav$order_effect[
  grepl("^treatment:order_idx$", behav$order_effect$term), c("model", "term", "p_value")
]
names(behav_order_interaction) <- c("Model", "Term", "p_value")
behav_order_interaction$Flag <- ifelse(behav_order_interaction$p_value < 0.05,
                                        "order-sensitive", "no confound detected")

endo_full <- endo$cell_cross_bh[, c("area", "analyte", "n", "re_spec", "beta", "SE",
                                    "F_stat", "df1", "df2", "p_raw", "p_BH")]
endo_full$sig_BH        <- endo$cell_cross_bh$p_BH < 0.05
endo_full$eta2_partial  <- partial_eta_sq(endo$cell_cross_bh$F_stat, endo$cell_cross_bh$df1, endo$cell_cross_bh$df2)
endo_full$omega2_partial <- partial_omega_sq(endo$cell_cross_bh$F_stat, endo$cell_cross_bh$df1, endo$cell_cross_bh$df2)

.r2_idx <- match(paste(endo$cell_cross_bh$area, endo$cell_cross_bh$analyte),
                 paste(endo$cell_r2$area, endo$cell_r2$analyte))
endo_full$R2_marginal    <- endo$cell_r2$R2_marginal[.r2_idx]
endo_full$R2_conditional <- endo$cell_r2$R2_conditional[.r2_idx]

endo_full <- endo_full[, c("area", "analyte", "n", "re_spec", "beta", "SE",
                           "F_stat", "df1", "df2", "eta2_partial", "omega2_partial",
                           "R2_marginal", "R2_conditional", "p_raw", "p_BH", "sig_BH")]
endo_dm <- endo_full[endo_full$area == "DM", ]

endo_order_dm <- endo$order_effect[grepl("\\[DM\\]|plasma_cortisol", endo$order_effect$outcome), ]

# Cortisol has no standalone R2 export (that's monoamine-only); partial
# eta-sq/omega-sq -- derivable straight from its own F/Df/Df.res -- serve as
# its variance-explained metric instead.
cortisol_headline <- endo$cortisol_anova
cortisol_headline$eta2_partial   <- partial_eta_sq(cortisol_headline$F, cortisol_headline$Df, cortisol_headline$Df.res)
cortisol_headline$omega2_partial <- partial_omega_sq(cortisol_headline$F, cortisol_headline$Df, cortisol_headline$Df.res)

# ---- Sensitivity / minimum-detectable-effect table --------------------------
# One row per F-test-based analysis. alpha_bh = the STRICTEST Benjamini-Hochberg
# threshold in that family (rank-1, i.e. 0.05/m), which brackets the worst case.

.sens <- rbind(
  data.frame(Stratum = "Behaviour (trial-level)",
             Analysis = behav_headline$Analysis,
             df1 = behav_headline$df1, df2 = behav_headline$df2,
             eta2_obs = behav_headline$eta2_partial,
             m_family = 4, stringsAsFactors = FALSE)[behav_headline$stat_type == "F-KR", ],
  data.frame(Stratum = "Cortisol (fish-level)",
             Analysis = "Plasma cortisol x Treatment",
             df1 = cortisol_headline$Df, df2 = cortisol_headline$Df.res,
             eta2_obs = cortisol_headline$eta2_partial,
             m_family = 1, stringsAsFactors = FALSE),
  data.frame(Stratum = "Monoamines, Dm (fish-level)",
             Analysis = paste0("Dm ", endo_dm$analyte, " x Treatment"),
             df1 = endo_dm$df1, df2 = endo_dm$df2,
             eta2_obs = endo_dm$eta2_partial,
             m_family = 7, stringsAsFactors = FALSE)
)

.sens$N_eff        <- .sens$df1 + .sens$df2 + 1
.sens$mdes_eta2_05 <- mapply(mdes_eta2, .sens$df1, .sens$df2, MoreArgs = list(alpha = 0.05, power = 0.80))
.sens$mdes_d_05    <- mapply(mdes_d,    .sens$df1, .sens$df2, MoreArgs = list(alpha = 0.05, power = 0.80))
.sens$mdes_eta2_bh <- mapply(function(a, b, m) mdes_eta2(a, b, alpha = 0.05 / m, power = 0.80),
                             .sens$df1, .sens$df2, .sens$m_family)
.sens$verdict <- ifelse(
  is.na(.sens$mdes_eta2_05), "—",
  ifelse(.sens$eta2_obs >= .sens$mdes_eta2_05,
         "above detection floor",
         "BELOW floor — null uninformative"))

sensitivity_tbl <- .sens[, c("Stratum", "Analysis", "df1", "df2", "N_eff",
                             "eta2_obs", "mdes_eta2_05", "mdes_d_05",
                             "mdes_eta2_bh", "verdict")]

# M5 association: minimum detectable |r| and the power actually available at
# each observed r -- the quantity that says how much the null result is worth.
.m5 <- endo$assoc[endo$assoc$area == "DM", c("area", "analyte", "n", "r_pearson", "p_pearson")]
.m5$mdes_r_05    <- mdes_r(15, alpha = 0.05,     power = 0.80)
.m5$mdes_r_bh    <- mdes_r(15, alpha = 0.05 / 7, power = 0.80)
.m5$power_at_obs <- power_at_r(.m5$r_pearson, n = 15, alpha = 0.05)
m5_sensitivity <- .m5

mono_flagged <- endo$pseudorep_mono[endo$pseudorep_mono$flag == TRUE, ]

# ---- 3. Reviewer/editor coverage checklist -----------------------------------
# Source: STATS_REVISION_INSTRUCTIONS.md + "Revision plan - to-do and
# roadmap.docx" (both in the project root). R1-1 and R2-1 intentionally
# omitted per explicit instruction.

coverage <- data.frame(
  Code = c("ED", "R1-2", "R1-3", "R1-4", "R1-5", "R1-6", "R1-7", "R1-8",
           "R2-2", "R2-3", "R2-4", "R2-5", "R2-6", "R2-7"),
  Source = c("Editor", rep("Reviewer 1", 7), rep("Reviewer 2", 6)),
  Requirement = c(
    "Remove/qualify every claim that swimming is rewarding; foreground the choice-vs-exercise-vs-rheotaxis limitation early.",
    "Soften 'intrinsically rewarding' language throughout (title, abstract, intro, discussion, conclusions).",
    "Frame the single 2-h exposure as an acute challenge, not evidence of stable exercise-reinforcement.",
    "Do not imply a directional/causal link between neurochemistry and behavioural preference.",
    "Acknowledge handling + terminal-sampling confound in cortisol/monoamine profiles.",
    "State the experimental unit (trial/group vs individual fish) explicitly for every analysis.",
    "Justify the multi-step analytical strategy (AICc, Box-Cox, ALR/CLR, Jacobs, beta-GLMM) and show robustness.",
    "Control for multiple comparisons across the analyte x brain-region tests.",
    "Add a rearing-history / trial-chronology figure; reconsider Figure 5.",
    "Justify assessing 'prior experience' without individually tagging fish.",
    "Show measured flow velocities by arena area and mean body length in Figure 3; show within-trial variability.",
    "Explain high between-trial variability (Fig 9) and test the treatment-order confound.",
    "Justify the macro-to-targeted brain-region approach; state neurotransmitter normalisation units.",
    "State how exercise-choice fish (which may switch zones) were captured and handled at sampling."
  ),
  Status = c(
    "Open", "Open", "Open", "Resolved", "Open", "Resolved", "Resolved", "Resolved",
    "Open", "Open", "Open", "Resolved", "Open", "Open"
  ),
  `Addressed in` = c(
    "Manuscript text only (Abstract/Discussion rewrite) -- depends on M1/M5 outcome below",
    "Manuscript text only (global find/replace pass)",
    "Manuscript text only (Discussion caveat)",
    "Part A4 below -- endocrine_behaviour_association.R (M5)",
    "Manuscript text only (Methods 2.4.3 + Discussion sentence)",
    "Part A1/A2/A3 + Part B below -- RE structure + pseudoreplication checks (M2)",
    "Part B below -- decision_log.csv + sensitivity tables (M3)",
    "Part A1/A3 below -- BH-FDR q-values (M1)",
    "New figure -- no script produces this; needs rearing/trial metadata compiled by hand",
    "Manuscript text only (Methods/Discussion) or drop the objective",
    "New figure annotation -- needs flow-velocity + body-length data overlaid on Fig 3",
    "Part A1/A3 below -- order_idx x treatment test (M4)",
    "Manuscript text only (Methods rationale + units clarification)",
    "Manuscript text only (Methods 2.4 sampling sentence)"
  ),
  check.names = FALSE
)

excluded_note <- paste(
  "Excluded on explicit instruction (prose/strategic, not statistical):",
  "R1-1 (innate-rheotaxis confound framing) and R2-1 (whether to split into two manuscripts)."
)

# ---- 3b. Glossary -------------------------------------------------------------
# Every abbreviation/statistical term used anywhere in the tables or prose
# below, defined once here so neither deliverable requires outside knowledge.

glossary <- data.frame(
  Term = c(
    "q  /  BH-adjusted q  /  p_BH",
    "Raw p  /  p_raw  /  Pr(>F)",
    "sig_BH  /  sig_raw",
    "BH-FDR  /  family",
    "β (beta)",
    "SE",
    "F statistic  /  F value",
    "Wald-chisq",
    "Partial η² (eta-squared)",
    "Partial ω² (omega-squared)",
    "R²m / R²c (marginal / conditional R²)",
    "Hedges' g",
    "Sensitivity analysis / MDES (minimum detectable effect size)",
    "Effective N",
    "Detection floor",
    "Post-hoc (observed) power — and why it is NOT reported",
    "df1 / df2  /  Df / Df.res",
    "F-KR  /  KR (Kenward-Roger)",
    "Satterthwaite",
    "AICc",
    "Random effect (RE)  /  (1|x) notation",
    "re_spec",
    "BC(λ = ...)",
    "n",
    "Brain region: Dm / POA / VV / VD",
    "Analyte",
    "order_idx",
    "treatment:order_idx",
    "Trial-level vs tank-level refit",
    "Direction_agree",
    "Jacobs' D",
    "CLR  /  ALR",
    "M1–M5",
    "ED  /  R1-x  /  R2-x",
    "r_pearson  /  rho_spearman",
    "Decision log step (a priori vs data-driven)"
  ),
  Meaning = c(
    "Benjamini-Hochberg false-discovery-rate-adjusted p-value. This is the significance measure to read wherever both raw p and q are shown — effects with q ≥ 0.05 are not significant, even if the raw p was < 0.05.",
    "The unadjusted p-value straight from the model, before any multiple-comparison correction.",
    "Significance flag at the BH-adjusted / raw threshold; an asterisk (*, **, ***) marks p or q < 0.05 / 0.01 / 0.001.",
    "The pre-declared group of tests corrected together (e.g. the 7 monoamines within one brain region). BH-FDR = the Benjamini-Hochberg correction procedure itself.",
    "The fitted model's effect-size estimate for the term in question (e.g. the treatment effect), on the model's response scale (which may be log- or Box-Cox-transformed — see BC(λ) below).",
    "Standard error of the β estimate.",
    "The test statistic from an ANOVA-type F-test on the fitted model.",
    "Wald chi-squared test statistic, used instead of an F-test where denominator degrees of freedom aren't defined (chiefly the beta-GLMMs for Jacobs' D).",
    "Proportion of variance in the response attributable to this term, holding the model's other terms constant. Computed here directly from the term's own F-statistic and degrees of freedom (η²p = F·df1 / (F·df1 + df2)) -- no refit needed. Cohen's rough benchmarks: ~0.01 small, ~0.06 medium, ~0.14 large. Only defined for an F-test; NA on the Wald-chisq (beta-GLMM) rows -- read R²m there instead.",
    "A small-sample-corrected version of partial η² that runs slightly lower and is less biased upward for small df2. Negative raw estimates (which can occur for near-null effects) are clamped to 0. Same F-test-only scope as partial η².",
    "Proportion of variance explained by the fixed effects alone (marginal, R²m) vs fixed + random effects together (conditional, R²c), via `performance::r2()`. Cortisol has no standalone R² export upstream -- partial η²/ω² (computed here from its own F/df) serve as its variance-explained metric instead.",
    "Bias-corrected standardised mean difference between treatment groups (a small-sample-corrected Cohen's d), computed from estimated marginal means and pooled SE, reported with a 95% CI. Only computed upstream for the Gaussian (F-KR) behavioural models, not the Jacobs' beta-GLMM rows (NA there -- read R²m instead). The underlying engine also reports a plain Cohen's d / odds ratio / incidence-rate ratio inline in its Word sub-report when the omnibus treatment test is significant, but only as narrative text, not an exportable table, so it is not reproduced here.",
    "The retrospective sample-size justification appropriate for an already-collected dataset: hold N, α and target power (0.80) fixed, and solve for the SMALLEST effect the design could have detected. Answers 'what could this study have found?' rather than the circular 'how much power did we have for the effect we happened to observe?'. Computed here on the noncentral F, using each model's own df1/df2.",
    "df1 + df2 + 1, taken from the model's own Kenward-Roger denominator df. Because KR df already reflect the random-effect structure, this inherits the correct experimental unit automatically -- trial (not fish) for the behavioural models, fish for the endocrine models -- without introducing any separate ICC or design-effect assumption. This is the same experimental-unit question Reviewer 1 raised in R1-6.",
    "The MDES for a given analysis. An observed effect ABOVE the floor was detectable by this design. An observed effect BELOW the floor means a non-significant result is uninformative -- it is evidence of absence of power, not evidence of absence of effect.",
    "Post-hoc power (recomputing power from the effect actually observed) is deliberately NOT reported anywhere here. It is a deterministic monotone function of the p-value, so it adds no information beyond p, and is widely regarded as a statistical fallacy. The sensitivity/MDES analysis above is the defensible retrospective alternative. The one legitimate exception, used in Part A5b, is quantifying the power available at an observed effect specifically in order to show how little a null result constrains -- which is an argument about the design, not a re-test of the finding.",
    "Numerator / denominator degrees of freedom for the F-test.",
    "Kenward-Roger degrees-of-freedom approximation for mixed-model F-tests — this project's primary/default inference method (less anti-conservative than a Wald test).",
    "Alternative denominator-df approximation, reported as a sensitivity check alongside KR.",
    "Corrected (small-sample) Akaike Information Criterion, used to select which random-effect structure best fits each model.",
    "A random intercept for the named grouping variable. `(1|trial)` = fish sampled from the same trial share a random intercept, so they are not treated as independent replicates.",
    "Which random-effect structure won AICc selection for that specific model/cell.",
    "Box-Cox power transform applied to the response variable before fitting, with the fitted power parameter λ.",
    "Sample size. What it counts (fish, trials, or groups) depends on the table — stated in the prose immediately above each table.",
    "Dm = dorsomedial pallium; POA = preoptic area; VV = ventral telencephalon (Vv); VD = dorsal part of the ventral telencephalon (Vd).",
    "The monoamine or metabolite measured: 5-HT, 5-HIAA, DA, DOPAC, NE, and the ratios 5-HIAA/5-HT and DOPAC/DA.",
    "Chronological rank (1–8) of the testing day/trial — the covariate used to test the treatment-order confound (R2-5 / M4).",
    "The interaction term testing whether the *treatment effect itself* changes across the order sequence — the actual confound test (not just an order main effect).",
    "The same model refit with the replicate defined at the trial level vs the tank level, to check that the experimental-unit assumption isn't driving the result (R1-6 / M2).",
    "Whether the trial-level and tank-level refits agree on the sign (direction) of the effect.",
    "Jacobs' electivity/preference index: a bounded (−1 to 1) measure of zone preference relative to zone availability.",
    "Centred / additive log-ratio transforms, used for compositional (proportion-of-zones) occupancy data.",
    "The five statistical-revision fixes from the rejection: M1 = multiple-comparison correction, M2 = pseudoreplication/experimental-unit fix, M3 = decisions log & robustness, M4 = order/day confound test, M5 = neurochemistry↔behaviour association.",
    "Reviewer/editor comment codes from the annotated manuscript (ED = editor; R1 = Reviewer 1; R2 = Reviewer 2), cross-referenced to the M1-M5 codes where a statistical fix applies.",
    "Pearson linear correlation coefficient / Spearman rank correlation coefficient.",
    "Each row of the Statistical Decisions Log (Part B) is marked a priori (decided before seeing the data) or data-driven (chosen after inspecting it), with its rationale."
  ),
  check.names = FALSE
)

# ---- 4. Compose the report as a list of sections -----------------------------
# Each section: list(level, title, paragraphs = chr vector, table = df or NULL,
# caption = chr or NULL)

sections <- list()
add_section <- function(level, title, paragraphs = NULL, table = NULL, caption = NULL, digits = 3) {
  sections[[length(sections) + 1]] <<- list(
    level = level, title = title, paragraphs = paragraphs,
    table = table, caption = caption, digits = digits
  )
}

add_section(1, "Reviewer & Editor Response \u2014 Consolidated Report", paragraphs = c(
  "**Manuscript:** Behavioural and neuroendocrine correlates of swimming exercise choice in juvenile Atlantic salmon (Bellio et al.), HB-D-26-00125.",
  "**Decision:** Rejected by Hormones and Behavior, 28 June 2026 (editor Cheryl McCormick; two reviewers).",
  sprintf("**Report generated:** %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  excluded_note
))

add_section(2, "Glossary of terms", paragraphs = c(
  "Every abbreviation and statistical term used anywhere in this report, defined once here. Column headers throughout have already been converted to plain-English labels; consult this table if a label is still unclear or you want the precise statistical definition behind it."
), table = glossary, caption = "Table G1. Glossary.", digits = 3)

add_section(2, "How this report is organised", paragraphs = c(
  "**Part A** \u2014 current headline results (behaviour, cortisol, monoamines, cross-domain association), each annotated with the relevant statistical-revision fix (multiple-comparison correction, pseudoreplication structure, order/day sensitivity).",
  "**Part B** \u2014 rigor & robustness documentation backing R1-7/M3 (a-priori-vs-data-driven decisions log, design-balance audit).",
  "**Part C** \u2014 the full reviewer/editor coverage checklist: every code from the annotated manuscript, whether it was resolved by a script (with a pointer into Part A/B) or remains an open prose/figure/design task no script can produce.",
  "**Part D** \u2014 provenance: source files, generation timestamps, companion documents.",
  "All numbers below are read directly from the CSV/txt outputs of `activity_analysis_STATS_choice_exp.R` (behaviour) and `analysis_b3.R` (cortisol + monoamines), last regenerated 2026-07-06 during the M1-M5 statistical-revision pass. Nothing here is hand-typed."
))

## --- Part A1: Behaviour -------------------------------------------------
add_section(2, "Part A \u2014 Current headline results")
add_section(3, "A1. Behavioural preference (Figures 6\u201310, Table 1 behaviour rows)", paragraphs = c(
  "Primary treatment effects (zone/flow preference, Jacobs' D by sub-zone, transitions, collective metrics), with Benjamini-Hochberg FDR applied within family \u2014 this is **R1-8 / M1**.",
  "Read q (`p_BH`), not raw p, for significance: effects with q \u2265 0.05 are not significant.",
  "**Effect sizes** (magnitude, independent of sample size): partial \u03b7\u00b2 and partial \u03c9\u00b2 are computed directly from each row's own F/df1/df2 (Cohen's benchmarks: ~0.01 small, ~0.06 medium, ~0.14 large) \u2014 defined only for the F-KR rows, not the Jacobs' D Wald-chisq rows, where R2m (marginal R\u00b2, already computed upstream for every row) is the variance-explained metric to read instead. Hedges' g [95% CI] (bias-corrected standardised mean difference) is shown where the upstream engine computed it \u2014 also F-KR-only."
), table = behav_headline, caption = "Table A1a. Primary behavioural family: raw p, BH-adjusted q, and effect sizes.")

add_section(3, "A1b. Treatment-order / day confound (R2-5 / M4)", paragraphs = c(
  "`order_idx` (chronological testing-day rank) x treatment interaction, tested per behavioural outcome. Non-significant interaction = no evidence the order/day sequence confounds the treatment effect."
), table = behav_order_interaction, caption = "Table A1b. treatment:order_idx interaction p-values, behaviour engine.")

add_section(3, "A1c. Experimental-unit / pseudoreplication check (R1-6 / M2)", paragraphs = c(
  "Same model refit at trial level vs tank level; `Direction_agree = FALSE` flags where the two levels disagree on effect direction \u2014 read alongside the primary model's own random-effect structure, not as a replacement for it."
), table = behav$pseudoreplication, caption = "Table A1c. Trial-level vs tank-level refits, behaviour engine.")

## --- Part A2: Cortisol ---------------------------------------------------
add_section(3, "A2. Cortisol (Table 1 observational, Figure 13)", paragraphs = c(
  "Kenward-Roger F-test, treatment (`condition`) effect on log-transformed plasma cortisol.",
  "Partial η²/ω² computed directly from this F/Df/Df.res (no separate R² export exists for cortisol upstream — see glossary)."
), table = cortisol_headline, caption = "Table A2a. Cortisol treatment effect and effect size.")

add_section(3, "A2b. Cortisol order/day confound (R2-5 / M4)", paragraphs = NULL,
  table = endo$order_effect[endo$order_effect$outcome == "plasma_cortisol", ],
  caption = "Table A2b. order_idx x condition test, cortisol.")

add_section(3, "A2c. Cortisol pseudoreplication / RE sensitivity (R1-6 / M2 & M3)", paragraphs = c(
  "Winning random-effect structure vs alternatives; `flag = TRUE` would mark a >30% shift in beta/SE between specifications (none here)."
), table = endo$pseudorep_cortisol, caption = "Table A2c. Cortisol RE sensitivity.")

## --- Part A3: Monoamines -------------------------------------------------
add_section(3, "A3. Monoamines (Table 2, Figures 11\u201312)", paragraphs = c(
  "All 28 analyte x brain-region cells, treatment effect, with BH-FDR applied **per region** (4 families of 7 analytes) \u2014 **R1-8 / M1**. `re_spec` shows the winning random-effect structure per cell (`(1|trial)` where selected) \u2014 **R1-6 / M2**.",
  "**Effect sizes:** partial \u03b7\u00b2/\u03c9\u00b2 computed from each cell's own F_stat/df1/df2 (Cohen's benchmarks: ~0.01 small, ~0.06 medium, ~0.14 large); R\u00b2m/R\u00b2c (marginal/conditional R\u00b2, from `cell_r2.csv`) show the proportion of variance explained by treatment alone vs treatment + random effects together.",
  sprintf("Headline dorsomedial pallium (Dm) analytes: %s",
          paste(sprintf("%s (\u03b2=%.3g, partial \u03b7\u00b2=%.3g, p_BH=%.3g, %s)",
                        endo_dm$analyte, endo_dm$beta, endo_dm$eta2_partial, endo_dm$p_BH,
                        ifelse(endo_dm$sig_BH, "significant (q<0.05)", "not significant")),
                collapse = "; "))
), table = endo_full, caption = "Table A3a. Full 28-cell monoamine treatment-effect table, all regions, with effect sizes.")

add_section(3, "A3b. Monoamine order/day confound, Dm region (R2-5 / M4)", paragraphs = c(
  "Order/day x treatment interaction for the three headline Dm analytes plus cortisol -- none reach significance, i.e. no evidence the order/day sequence explains the Dm effects.",
  "Full 4-region, 7-analyte table (58 rows) available in `order_effect_summary.csv`; one incidental order effect elsewhere in the table is worth noting for completeness: NE in area VD shows a significant main order effect (p = 0.0033) -- not one of the manuscript's headline claims, flagged here for transparency."
), table = endo_order_dm, caption = "Table A3b. order_idx test, Dm analytes + cortisol.")

add_section(3, "A3c. Monoamine pseudoreplication sensitivity (R1-6 / M2 & M3)", paragraphs = c(
  sprintf("%d of %d analyte x region cells had a >30%% beta/SE shift across random-effect specifications (flagged below); the remainder were stable.",
          nrow(mono_flagged), length(unique(paste(endo$pseudorep_mono$analyte, endo$pseudorep_mono$area)))),
  "Full sensitivity table (66 rows, all specifications x all cells) available in `pseudorep_check_monoamines.csv`."
), table = mono_flagged, caption = "Table A3c. Flagged (>30% shift) random-effect sensitivity rows.")

## --- Part A4: Cross-domain association -----------------------------------
add_section(3, "A4. Neurochemistry \u2194 behaviour association (R1-4 / M5)", paragraphs = c(
  "Trial-level Pearson/Spearman correlation between behavioural flow-preference and trial-aggregated monoamines, computed only where both exist for the same physical trial (n = 15/16 trials linked via tank + testing date + treatment).",
  "**None of the three manuscript-headline Dm analytes (DA, 5-HT, DOPAC) reaches significance at the trial level.** Crucially, this is **not** evidence that no association exists: the sensitivity analysis in Part A5b shows this design (n = 15 trials) had only 5–21% power at the observed correlation magnitudes and could not have detected anything below |r| ≈ 0.64. The correct statement is that the study is **uninformative** about a neurochemistry↔behaviour association of small-to-moderate size — not that it rules one out. Directional language remains unsupported, but so does any claim of a demonstrated null.",
  paste0("**Limitation statement (from `endocrine_behaviour_association_LIMITATION.txt`):**\n\n> ", gsub("\n", "\n> ", assoc_limitation))
), table = endo$assoc, caption = "Table A4. Trial-level association, behaviour vs monoamines, by region x analyte.")

## --- Part A5: sensitivity / sample size ------------------------------------
add_section(3, "A5. Sample-size justification — sensitivity analysis", paragraphs = c(
  "**What this is.** A retrospective sample-size justification for an already-collected dataset. N, α and target power (0.80) are held fixed and the **smallest detectable effect** is solved for. This answers \"what could this design have found?\" — which is answerable — rather than \"how much power did we have for the effect we happened to observe?\", which is circular.",
  "**Post-hoc (observed) power is deliberately not reported.** It is a monotone function of the p-value and therefore adds nothing beyond p. Reporting it would invite exactly the R1-7 criticism about analytical rigour.",
  "**Experimental unit (ties to R1-6).** Each row uses its own model's Kenward-Roger denominator df, so the effective N inherits the correct replicate automatically — **trial** for the behavioural models (N = 16 trials, 8/treatment, giving df₂ = 11 and an effective N of 13), **fish** for the endocrine models. No separate ICC or design-effect assumption is introduced. Computed on the noncentral F and cross-checked against `power.t.test()`.",
  "**PLACEHOLDER — requires author confirmation:** whether an *a priori* power analysis was conducted before data collection could not be determined from the analysis outputs. If none was, Methods should say so plainly and cite the practical constraint that set N (facility capacity, fish availability, ethical limits), with this sensitivity analysis as the retrospective justification. **Do not assert an a priori analysis that was not performed.**"
), table = sensitivity_tbl, caption = "Table A5a. Minimum detectable effect size by analysis stratum (power = 0.80).")

add_section(3, "A5b. What the M5 null result can and cannot support (R1-4)", paragraphs = c(
  "Minimum detectable |r| at n = 15 trials, and the power actually available at each observed correlation. Solved on the exact noncentral t (the Fisher-z approximation runs ~5 percentage points conservative at this n).",
  "**Reading:** the three headline Dm analytes had **5–21% power** at their observed effect magnitudes. A non-significant result under those conditions is uninformative — it constrains nothing about a true association of small-to-moderate size. This strengthens rather than weakens the R1-4 response: it removes directional language *and* declines to over-claim a null."
), table = m5_sensitivity, caption = "Table A5b. M5 association — detectability at n = 15 trials.")

add_section(3, "A5c. Consequences for interpretation", paragraphs = c(
  "**1. Two non-significant behavioural results sit below the detection floor.** Switches/session (observed η² = 0.285) and school speed (0.276) fall under the α = .05 floor of η² = 0.42, so \"not significant\" is uninformative for these, not null.",
  "**2. The significant behavioural effects clear the floor but sit close to it.** A design that can only detect η² ≥ 0.42 (d ≈ 1.7) will, when it does detect something, tend to **overestimate** its magnitude (the winner's-curse / effect-inflation problem). NND (0.485) and polarisation (0.705) are genuine detections, but their point estimates should be presented as likely upper bounds rather than precise magnitudes. Stating this proactively speaks directly to the editor's core objection that the data overstate the case.",
  "**3. The endocrine analyses are better powered than the behavioural ones**, because their replicate is the fish rather than the trial — cortisol detects down to η² ≈ 0.10 (d ≈ 0.66) and the per-cell monoamine models to η² ≈ 0.19–0.23. The BH-FDR correction raises each of these floors further, as the strictest-α column shows.",
  "**4. Reporting recommendation.** Report partial η² with a **90%** CI (not 95%): F-tests are one-sided, and the 90% interval excludes zero exactly when the F-test is significant. Keep reporting ω² alongside η², since η² is biased upward and a future researcher using it for an a-priori power calculation would systematically **underpower** their replication."
))

## --- Part B: rigor ---------------------------------------------------------
add_section(2, "Part B \u2014 Rigor & robustness documentation (R1-7 / M3)", paragraphs = c(
  "**Endocrine engine:** full a-priori-vs-data-driven Statistical Decisions Log (15 rows) below.",
  "**Behaviour engine:** has no standalone decision-log file; the equivalent robustness evidence is `bh_correction_summary.csv` + `aicc_selection.csv` (per model, inside each `*_aggregated*` output folder) -- pre-existing AICc random-effect selection, carried over unchanged, cross-referenced rather than duplicated here.",
  "**Known pre-existing inconsistency (flagged, not fixed in this pass):** `design_balance.csv` reports `tank_density_1to1 = FALSE` and `cor_tank_density = 0` -- fish_density is not 1:1 with tank (each tank cycles through all 4 densities). This contradicts an existing in-code comment used to justify dropping fish_density as an aggregated-level random-effect candidate. Out of M1-M5 scope; noted here for the manuscript team's judgement."
), table = endo$decision_log, caption = "Table B1. Endocrine engine -- statistical decisions log.")

add_section(3, "B2. Design-balance audit (behaviour engine)", table = behav$design_balance,
            caption = "Table B2. tank/density/date collinearity checks.")

## --- Part C: coverage checklist --------------------------------------------
add_section(2, "Part C \u2014 Full reviewer/editor coverage checklist", paragraphs = c(
  "Every code from the annotated manuscript except R1-1 and R2-1 (see scope note above). 'Resolved' = a script now produces the evidence (see the Part A/B pointer); 'Open' = a prose, figure, or methods-description change with no analysis behind it -- no script can close these, they need a manuscript edit."
), table = coverage, caption = "Table C1. Reviewer/editor coverage checklist.", digits = 6)

## --- Part D: provenance -----------------------------------------------------
add_section(2, "Part D \u2014 Provenance", paragraphs = c(
  sprintf("Behaviour source directory: `%s` (files dated 2026-07-06).", BEHAV_DIR),
  sprintf("Endocrine source directory: `%s` (files dated 2026-07-06).", ENDO_DIR),
  "Companion documents (project root): `stat_map.md` (statistic -> script traceability), `CHANGELOG_stats.md` (old value -> new value per M1-M5 change), `reviewer_coverage_table.md` (the original 5-code M1-M5 coverage table this report supersedes/extends), `STATS_REVISION_INSTRUCTIONS.md` (the M1-M5 task spec), `Revision plan - to-do and roadmap.docx` (the ED/R1-x/R2-x to-do list and phased roadmap this report's Part C is drawn from).",
  sprintf("This report generated by `choice R pipeline/scripts/04_reporting/build_reviewer_response_report.R`, R %s.",
          paste(R.version$major, R.version$minor, sep = "."))
))

# ---- 5. Render: Markdown -----------------------------------------------------

md_lines <- character(0)
for (s in sections) {
  md_lines <- c(md_lines, paste0(strrep("#", s$level), " ", s$title), "")
  if (!is.null(s$paragraphs)) {
    for (p in s$paragraphs) md_lines <- c(md_lines, p, "")
  }
  if (!is.null(s$table)) {
    if (!is.null(s$caption)) md_lines <- c(md_lines, paste0("*", s$caption, "*"), "")
    md_lines <- c(md_lines, md_table(s$table, digits = s$digits), "")
  }
}
writeLines(md_lines, OUT_MD, useBytes = TRUE)
message("[build_reviewer_response_report] Wrote: ", OUT_MD)

# ---- 6. Render: Word ---------------------------------------------------------

doc <- read_docx()
style_for <- function(level) if (level == 1) "heading 1" else if (level == 2) "heading 2" else "heading 3"

for (s in sections) {
  doc <- body_add_par(doc, s$title, style = style_for(s$level))
  if (!is.null(s$paragraphs)) {
    for (p in s$paragraphs) {
      p_plain <- gsub("\\*\\*|`", "", p)
      doc <- body_add_par(doc, p_plain, style = "Normal")
    }
  }
  if (!is.null(s$table)) {
    if (!is.null(s$caption)) doc <- body_add_par(doc, s$caption, style = "Normal")
    doc <- body_add_flextable(doc, ft_table(s$table, digits = s$digits))
    doc <- body_add_par(doc, " ", style = "Normal")
  }
}

doc <- doc |>
  body_end_section_continuous() |>
  body_set_default_section(
    prop_section(
      page_size = page_size(width = 11.69, height = 8.27, orient = "landscape"),
      page_margins = page_mar(top = 0.4, bottom = 0.4, left = 0.4, right = 0.4)
    )
  )

print(doc, target = OUT_DOCX)
message("[build_reviewer_response_report] Wrote: ", OUT_DOCX)

invisible(list(md = OUT_MD, docx = OUT_DOCX))
