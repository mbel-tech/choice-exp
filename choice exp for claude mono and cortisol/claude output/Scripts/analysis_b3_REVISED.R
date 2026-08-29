# =============================================================================
# B3 plasma cortisol + monoamines — assay matching, CV%, and sex analysis
# Plan reference: C:\Users\marti\.claude\plans\you-are-claude-opus-cheeky-badger.md
# Stats-revision reference: D:\CHOICE R SCRIPTS\STATS_REVISION_INSTRUCTIONS.md (M1-M5)
# Scope: B3 / batch 3 only. Mucus, MU*, batch 1-2 plates excluded.
# Phase A: DHARMa diagnostics, performance::check_model, BayesFactor, omega_squared
# Phase B: tank/trial nesting in cortisol model; (1|tank) added to monoamine models
# Phase C (stats revision): (1|trial) added to monoamine RE candidates (M2);
#   trial-order/day fixed-effect test added (M4); reproducibility seed + sessionInfo (M3).
# =============================================================================

set.seed(20260706)

# ---- 0. Packages ------------------------------------------------------------
need <- c("readxl","writexl","openxlsx",
          "multcomp","MuMIn","MASS",  # load BEFORE dplyr so dplyr::select masks them
          "dplyr","tidyr","stringr","purrr",
          "forcats","lubridate","janitor","broom","broom.mixed","car",
          "emmeans","lme4","lmerTest","performance","rstatix",
          "ggplot2","patchwork","ggpubr","ggtext","officer","flextable","rvg",
          "DHARMa","effectsize","BayesFactor","cols4all","cowplot",
          "tibble")
inst <- rownames(installed.packages())
miss <- setdiff(need, inst)
if (length(miss)) {
  user_lib <- Sys.getenv("R_LIBS_USER", file.path(Sys.getenv("APPDATA"), "R", "win-library", getRversion()))
  dir.create(user_lib, showWarnings = FALSE, recursive = TRUE)
  install.packages(miss, repos = "https://cloud.r-project.org", lib = user_lib)
}
invisible(lapply(need, function(p) suppressPackageStartupMessages(library(p, character.only = TRUE))))

# ---- 1. Paths ---------------------------------------------------------------
ROOT     <- file.path(PROJECT_ROOT, "choice exp for claude mono and cortisol")
# Revised run writes to its OWN output tree so the 2026-07-06 originals under
# "claude output" are preserved for comparison. Inputs are still read from ROOT.
OUT      <- file.path(ROOT, "claude output REVISED")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

DATA_XLSX <- file.path(ROOT, "choice exp datasets.xlsx")
MONO_XLSX <- file.path(ROOT, "MD monoamine concentrations choice exp MAY 2025.xlsx")
NICE_XLSX <- file.path(ROOT, "cortisol", "cortisol choice exp nice plates.xlsx")

PLATE_FILES <- list(
  `5`   = list(file = NICE_XLSX,                                              sheet = "my assays plate 5",
               plate_date = as.Date("2024-07-25")),
  `6`   = list(file = file.path(ROOT, "cortisol", "plate 6 my assays.xlsx"),  sheet = "Sheet1",
               plate_date = as.Date("2024-08-02")),
  `7`   = list(file = file.path(ROOT, "cortisol", "my assays plate 7.xlsx"),  sheet = "Sheet1",
               plate_date = as.Date("2024-08-07")),
  `8`   = list(file = file.path(ROOT, "cortisol", "plate 8 my assays.xlsx"),  sheet = "plate 8 my assays",
               plate_date = as.Date("2024-08-08")),
  `Evg` = list(file = NICE_XLSX,                                              sheet = "my assay Evgenia",
               plate_date = as.Date("2024-08-21"))
)

NOMINAL_STD <- c(Standard1 = 3200, Standard2 = 1600, Standard3 = 800,
                 Standard4 = 400,  Standard5 = 200,  Standard6 = 100, Standard7 = 50)

# ---- 2. Helpers -------------------------------------------------------------
pad_b3 <- function(x) {
  x <- as.character(x)
  m <- stringr::str_match(x, "^B3[_ ]?([0-9]+)$")
  ifelse(!is.na(m[,2]), sprintf("B3_%02d", as.integer(m[,2])), x)
}

to_num <- function(x) suppressWarnings(as.numeric(x))

tidy_model <- function(fit, ...) {
  if (inherits(fit, "lmerMod") || inherits(fit, "lmerModLmerTest"))
    broom.mixed::tidy(fit, conf.int = TRUE, ...)
  else
    broom::tidy(fit, conf.int = TRUE, ...)
}

fmt_b <- function(b) sprintf("%.2f", b)

# ---- Unified number-display convention (author decision, 2026-08-07) -------
# fmt_F: F-statistics -> 4 significant figures, trailing zeros trimmed
#        (10.348 -> "10.35"; 9.724 stays "9.724"). Duplicated verbatim from
#        00_shared/number_formatting.R because this engine is run standalone.
# fmt3 : p (once below the "< 0.001" floor), partial eta-squared, partial
#        omega-squared, R^2 -> 3 decimals, dropped to 2 ONLY when the 3rd
#        decimal digit is exactly zero (0.150 -> "0.15"; 0.019 stays "0.019").
# Both are fully vectorised (unlike a naive scalar version with `if`), since
# this engine applies them inside dplyr::mutate() over whole columns.
fmt_F <- function(x) {
  s <- format(signif(x, 4), scientific = FALSE, trim = TRUE)
  # Only trim trailing zeros when a decimal point is present -- otherwise an
  # integer denominator df like 30 is mangled to "3" (see project number-
  # formatting defect log).
  has_dot <- grepl(".", s, fixed = TRUE)
  s[has_dot] <- sub("\\.$", "", sub("0+$", "", s[has_dot]))
  ifelse(is.finite(x), s, "NA")
}
fmt3 <- function(x) {
  s <- sprintf("%.3f", x)
  s <- ifelse(grepl("0$", s), substr(s, 1, nchar(s) - 1), s)
  ifelse(is.finite(x), s, "NA")
}

# fmt_Fstat / fmt_es3 (2026-08-09): as fmt_F / fmt3, but a value below 0.001
# collapses to "< 0.001" instead of rendering either a long 4-sig-fig decimal
# ("0.0009258") or an uninformative rounded zero ("0.00"). Applied ONLY to the
# F statistic and to partial eta-squared -- denominator df keep fmt_F (a df is
# never < 1) and p-values keep fmt_p (which owns its own floor).
# 2026-08-28: two decimals, not four significant figures. Table 2 in the manuscript
# reports F to 2 dp, and a figure caption reading 5.927 next to a table reading 5.93 is
# the kind of mismatch a reviewer notices. fmt_F is still used for everything else.
fmt_Fstat <- function(x) {
  ifelse(!is.finite(x), "NA",
         ifelse(abs(x) < 0.001, "< 0.001", formatC(x, format = "f", digits = 2)))
}
fmt_es3 <- function(x) {
  ifelse(!is.finite(x), "NA",
         ifelse(abs(x) < 0.001, "< 0.001", fmt3(x)))
}

# Template-style p-value formatter: 3 sig figs when p > 0.0001 else "< 0.0001"; appends sig stars
fmt_p <- function(p) {
  if (length(p) == 0) return("")
  vapply(p, function(pi) {
    if (is.na(pi)) return("NA")
    base <- if (pi < 0.0001) "< 0.0001" else format(signif(pi, 3), scientific = FALSE)
    star <- if (is.na(pi)) "" else if (pi < 0.001) " ***" else if (pi < 0.01) " **" else if (pi < 0.05) " *" else ""
    paste0(base, star)
  }, character(1))
}

fmt_estse <- function(est, se, digits = 3) {
  if (length(est) == 0 || is.na(est)) return("NA")
  sprintf(sprintf("%%.%df ± %%.%df", digits, digits), est, se)
}

# M4 (R2-5): fixed vocabulary of the 8 preference-experiment testing dates, shared
# by the cortisol and monoamine datasets, so order_idx (1..8) encodes identically
# in both regardless of which subset of dates each table happens to contain.
ALL_TRIAL_DATES <- as.Date(sprintf("2024-06-%02d", 15:22))
make_order_idx  <- function(date_col) as.numeric(factor(as.Date(date_col), levels = ALL_TRIAL_DATES))

# ---- 2b. Statistical helpers (template-style) ------------------------------
# Sum-to-zero contrasts on listed factor columns of a data frame.
set_sum_contrasts <- function(df, vars) {
  for (v in vars) if (v %in% names(df) && is.factor(df[[v]])) {
    contrasts(df[[v]]) <- contr.sum(nlevels(df[[v]]))
  }
  df
}

# AICc-based RE selection. Returns list(table=AICc data.frame, winner=refit-REML model, label=winner-label)
#
# CHANGED 2026-08-08 (Reviewer 1 comments 6 & 7; Barr et al. 2013; Hurlbert 1984):
# `force_re` names a random-effect term that the DESIGN imposes -- here (1|trial),
# because treatment was applied to the arena trial while the response is measured
# on the individual fish. It is appended to every candidate, so AICc can now only
# choose among OPTIONAL nuisance terms (plate, date, tank) and can no longer
# delete the design term. Previously (1|trial) competed head-to-head against
# those nuisance terms and lost: it was dropped from the cortisol model, from
# every region-aggregated monoamine model, and from 15 of 24 per-cell models
# (5 of which ended up with no clustering term at all).
#
# Singular candidates are also no longer excluded outright. That rule
# systematically penalised (1|trial) -- the very term the design requires --
# because a small between-trial variance frequently estimates to the boundary.
# Since the forced term is now present in every candidate, singularity carries
# information only about the optional terms, so non-singular fits are still
# preferred, but an all-singular set no longer discards the design structure.
select_re_aicc <- function(fixed_formula, re_list, data, family = NULL,
                           force_re = NULL) {
  # Build the candidate set with the design-imposed term folded in (de-duplicated).
  if (!is.null(force_re) && nzchar(force_re)) {
    lbls <- names(re_list)
    new_list <- list(); new_names <- character(0)
    # the design term on its own is always a candidate
    new_list[[force_re]] <- force_re; new_names <- c(new_names, force_re)
    for (lbl in lbls) {
      trm <- re_list[[lbl]]
      if (identical(trimws(trm), trimws(force_re))) next        # already covered
      if (grepl(force_re, trm, fixed = TRUE)) {                  # already contains it
        new_list[[lbl]] <- trm; new_names <- c(new_names, lbl); next
      }
      combo_lbl <- paste0(lbl, "+", force_re)
      new_list[[combo_lbl]] <- paste(trm, "+", force_re)
      new_names <- c(new_names, combo_lbl)
    }
    # De-duplicate models that differ only in the ORDER of their random terms
    # (e.g. "(1|trial) + (1|plate)" vs "(1|plate) + (1|trial)"), which would
    # otherwise be fitted and reported twice.
    .canon <- function(s) paste(sort(trimws(strsplit(s, "+", fixed = TRUE)[[1]])),
                                collapse = " + ")
    keep <- !duplicated(vapply(new_list, .canon, character(1)))
    re_list <- new_list[keep]
  }

  # REML (not ML): this loop compares RANDOM-effect structures with the fixed
  # formula held constant across every candidate -- exactly the case where a
  # REML-based AICc comparison is valid (Pinheiro & Bates 2000; Zuur et al.
  # 2009). ML is required instead when candidates differ in their FIXED
  # effects (see the corresponding STEP5 fixed-structure fix, df_method_memo
  # cross-reference); this function never does that, so REML applies uniformly.
  fits <- list(); aics <- c(); singular <- c()
  for (lbl in names(re_list)) {
    re_term <- re_list[[lbl]]
    full_form <- update(fixed_formula, paste(". ~ . +", re_term))
    fit <- tryCatch(
      suppressMessages(suppressWarnings(
        lmerTest::lmer(full_form, data = data, REML = TRUE)
      )),
      error = function(e) NULL
    )
    if (is.null(fit)) { aics[lbl] <- NA_real_; singular[lbl] <- NA; next }
    aics[lbl]     <- tryCatch(MuMIn::AICc(fit), error = function(e) AIC(fit))
    singular[lbl] <- lme4::isSingular(fit, tol = 1e-4)
    fits[[lbl]]   <- fit
  }
  ok <- !is.na(aics) & !singular
  # Prefer non-singular; but never discard the design term just because every
  # candidate containing it is singular.
  if (!any(ok) && !is.null(force_re)) ok <- !is.na(aics)
  if (!any(ok)) {
    return(list(table = data.frame(re = names(re_list), AICc = aics, singular = singular,
                                   selected = FALSE), winner = NULL, label = NA_character_))
  }
  win_lbl <- names(aics)[ok][which.min(aics[ok])]
  # Candidates are already fit REML (see above); reuse the winner directly
  # rather than refitting, which would risk a spurious divergence between the
  # AICc-selected fit and the fit actually used for inference.
  win_reml <- fits[[win_lbl]]
  tab <- data.frame(re = names(re_list), AICc = aics, singular = singular,
                    selected = names(re_list) == win_lbl,
                    dAICc = aics - min(aics, na.rm = TRUE),
                    row.names = NULL)
  # `re_formula` is the winning RE term string FROM THIS FUNCTION'S OWN
  # (possibly force_re-expanded) candidate list, e.g. re_list[["(1|plate)+(1|trial)"]].
  # Callers must use this, not re-index their original pre-expansion candidate
  # list by win_lbl -- when force_re generates a new combo label (e.g.
  # "(1|plate)+(1|trial)" from a lone "(1|plate)" candidate), that label need
  # not exist as a literal key in the caller's own list (which may instead
  # have the same combination spelled "(1|trial)+(1|plate)"), and indexing
  # with it silently returns NULL rather than erroring.
  list(table = tab, winner = win_reml, label = win_lbl, re_formula = re_list[[win_lbl]])
}

# Box-Cox + Shapiro on raw / log / boxcox-transformed
transform_check <- function(y, label = "") {
  y <- y[is.finite(y) & y > 0]
  if (length(y) < 3) return(data.frame(outcome = label, n = length(y),
                                       sw_raw = NA, sw_log = NA, sw_bc = NA, bc_lambda = NA))
  sw_raw <- tryCatch(shapiro.test(y)$p.value, error = function(e) NA_real_)
  sw_log <- tryCatch(shapiro.test(log(y))$p.value, error = function(e) NA_real_)
  # Grid search over lambda: optimise Shapiro-Wilk p directly (avoids MASS::boxcox scoping issues)
  bc_lambdas <- seq(-2, 2, by = 0.1)
  bc_sw <- vapply(bc_lambdas, function(lam) {
    bc_y <- if (abs(lam) < 1e-6) log(y) else (y^lam - 1) / lam
    if (!all(is.finite(bc_y))) return(NA_real_)
    tryCatch(shapiro.test(bc_y)$p.value, error = function(e) NA_real_)
  }, numeric(1))
  best_i <- which.max(bc_sw)
  if (length(best_i) == 1 && !is.na(bc_sw[best_i])) {
    lam  <- bc_lambdas[best_i]; sw_bc <- bc_sw[best_i]
  } else { lam <- NA_real_; sw_bc <- NA_real_ }
  data.frame(outcome = label, n = length(y),
             sw_raw = sw_raw, sw_log = sw_log, sw_bc = sw_bc, bc_lambda = lam)
}

# Levene: takes a formula y ~ group; returns data.frame row
levene_check <- function(formula, data, label = "") {
  res <- tryCatch(car::leveneTest(formula, data = data), error = function(e) NULL)
  if (is.null(res)) return(data.frame(outcome = label, levene_F = NA, df1 = NA, df2 = NA, p = NA))
  data.frame(outcome = label, levene_F = res$`F value`[1],
             df1 = res$Df[1], df2 = res$Df[2], p = res$`Pr(>F)`[1])
}

# Adaptive transform selection based on transform_check() output.
# Prefers Box-Cox when available and Shapiro p > log() by factor ≥ 1.5;
# collapses to log() when lambda near 0. Falls back to log() if BC unavailable.
select_transform <- function(norm_row) {
  sw_log <- norm_row$sw_log
  sw_bc  <- norm_row$sw_bc
  lam    <- norm_row$bc_lambda
  if (!is.na(sw_bc) && !is.na(sw_log) && is.finite(sw_bc) && sw_bc > sw_log * 1.5 &&
      !is.na(lam) && is.finite(lam)) {
    if (abs(lam) < 0.1) {
      return(list(type = "log", lambda = NA_real_,
                  label = "log(value)",
                  rationale = sprintf(
                    "Box-Cox λ=%.3f (~0); log() used as equivalent. Shapiro-Wilk: raw=%.4f, Shapiro-Wilk (log)=%.4f, Shapiro-Wilk (Box-Cox)=%.4f.",
                    lam, norm_row$sw_raw, sw_log, sw_bc)))
    }
    return(list(type = "boxcox", lambda = lam,
                label = sprintf("BC(λ=%.3f)", lam),
                rationale = sprintf(
                  "Box-Cox (λ=%.3f) selected: Shapiro-Wilk (Box-Cox)=%.4f > Shapiro-Wilk (log)=%.4f (ratio=%.2f).",
                  lam, sw_bc, sw_log, sw_bc / sw_log)))
  }
  rat <- if (!is.na(sw_log) && is.finite(sw_log))
    sprintf("log() selected: Shapiro-Wilk (log)=%.4f, Shapiro-Wilk (Box-Cox)=%s (ratio<1.5).",
            sw_log,
            if (!is.na(sw_bc)) sprintf("%.4f", sw_bc) else "not available")
  else "log() selected (default for positive concentration data)."
  list(type = "log", lambda = NA_real_, label = "log(value)", rationale = rat)
}

# Apply selected transform; output assigned to d$log_v so all formulas stay unchanged.
apply_transform <- function(y, trans) {
  if (trans$type == "log") return(log(y))
  lam <- trans$lambda
  if (abs(lam) < 1e-6) log(y) else (y^lam - 1) / lam
}

# DHARMa residual panel saved to PNG path (returns path or NA)
dharma_panel <- function(model, path, title = "") {
  ok <- tryCatch({
    sim <- DHARMa::simulateResiduals(model, plot = FALSE)
    png(path, width = 900, height = 500, res = 120)
    plot(sim, title = title)
    dev.off()
    TRUE
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    message("DHARMa failed [", title, "]: ", conditionMessage(e))
    FALSE
  })
  if (ok) path else NA_character_
}

# R2m / R2c via performance::r2; falls back gracefully
r2_row <- function(model) {
  r2 <- tryCatch(performance::r2(model), error = function(e) NULL)
  icc_val <- if (inherits(model, "merMod"))
    tryCatch(as.numeric(performance::icc(model)$ICC_adjusted), error = function(e) NA_real_)
  else NA_real_
  if (is.null(r2))
    return(data.frame(R2_marginal = NA_real_, R2_conditional = NA_real_, ICC_adjusted = icc_val))
  if (!is.null(r2$R2_marginal))
    data.frame(R2_marginal = as.numeric(r2$R2_marginal),
               R2_conditional = as.numeric(r2$R2_conditional),
               ICC_adjusted = icc_val)
  else
    data.frame(R2_marginal = as.numeric(r2$R2[1]), R2_conditional = NA_real_,
               ICC_adjusted = icc_val)
}

# ---------------------------------------------------------------------------
# WARNING (fixed 2026-08-08): car:::Anova.merMod has NO `ddf` argument -- its
# formals are (mod, type, test.statistic, vcov., singular.ok, ...), so a `ddf=`
# passed here is silently swallowed by `...`, and `test = "F"` ALWAYS routes to
# Kenward-Roger via pbkrtest::vcovAdj. The previous anova_satt() therefore
# returned a byte-identical copy of anova_kr(), and the "Satterthwaite
# sensitivity row" advertised in decision_log 2.5 never actually existed.
# Verified: car KR p = 0.0029585743 == car "Satterthwaite" p = 0.0029585743,
# whereas a genuine lmerTest Satterthwaite gives p = 0.0028638571
# (F 9.5207, df2 73.4604 vs KR F 9.4561, df2 73.1275).
# Satterthwaite must therefore go through lmerTest, which does honour `ddf`.
# ---------------------------------------------------------------------------

# Type III ANOVA with Kenward-Roger F (preferred for unbalanced LMMs).
# car::Anova(test = "F") on a merMod IS the Kenward-Roger test; no ddf arg.
anova_kr <- function(model, type = 3) {
  if (inherits(model, "lmerMod") || inherits(model, "lmerModLmerTest")) {
    res <- car::Anova(model, type = type, test = "F")
  } else {
    res <- car::Anova(model, type = type)
  }
  res
}

# Type III ANOVA with Satterthwaite denominator df, for genuine sensitivity.
# Returns a car-shaped table (rownames = terms; columns F, Df, Df.res, Pr(>F))
# so downstream rownames_to_column("term") and column lookups keep working.
anova_satt <- function(model, type = 3) {
  if (!(inherits(model, "lmerMod") || inherits(model, "lmerModLmerTest")))
    return(car::Anova(model, type = type))

  # Re-wrapping a model that is ALREADY class lmerModLmerTest via
  # as_lmerModLmerTest() again corrupts state that anova()'s internal
  # Satterthwaite refit relies on (silently returns a 0-row table, no error).
  # Only promote plain lme4 fits; pass lmerModLmerTest fits through untouched.
  ml <- if (inherits(model, "lmerModLmerTest")) model
        else tryCatch(lmerTest::as_lmerModLmerTest(model), error = function(e) NULL)
  if (is.null(ml)) return(NULL)
  ty <- if (is.numeric(type)) c("I", "II", "III")[type] else type
  a <- tryCatch(stats::anova(ml, type = ty, ddf = "Satterthwaite"),
                error = function(e) NULL)
  if (is.null(a) || !nrow(a)) return(NULL)

  out <- data.frame(
    F         = a[["F value"]],
    Df        = a[["NumDF"]],
    Df.res    = a[["DenDF"]],
    `Pr(>F)`  = a[["Pr(>F)"]],
    check.names = FALSE
  )
  rownames(out) <- rownames(a)
  out
}

# emmeans + Tukey + multcomp::cld; returns named list of tables.
# When `specs` is an interaction (e.g. ~ condition * timepoint), also compute
# the SIMPLE-EFFECTS view: pairwise within each level of each stratifier.
# Stored as $simple = list(<focal>_within_<stratifier> = data.frame, ...).
# Simple effects are typically what aligns with patterns shown in graphs.
emm_tukey_cld <- function(model, specs, by = NULL) {
  emm <- tryCatch(emmeans::emmeans(model, specs = specs, by = by), error = function(e) NULL)
  if (is.null(emm)) return(NULL)
  emm_tab <- as.data.frame(emm)
  pair_tab <- tryCatch(as.data.frame(pairs(emm, adjust = "tukey", infer = c(TRUE, TRUE))),
                       error = function(e) NULL)
  cld_tab <- tryCatch({
    cld <- multcomp::cld(emm, Letters = letters, adjust = "tukey")
    as.data.frame(cld)
  }, error = function(e) NULL)

  # Simple-effects (per-stratum Tukey). Only when no `by` was passed AND
  # specs is an interaction (i.e. RHS contains '*' or ':').
  simple <- list()
  if (is.null(by)) {
    rhs <- if (inherits(specs, "formula"))
             paste(deparse(specs[[length(specs)]]), collapse = " ")
           else as.character(specs)[length(as.character(specs))]
    vars <- trimws(unlist(strsplit(rhs, "\\s*[\\*:]\\s*")))
    vars <- vars[nzchar(vars) & vars != "1"]
    if (length(vars) >= 2) {
      for (stratifier in vars) {
        focal <- setdiff(vars, stratifier)
        focal_spec <- paste(focal, collapse = " * ")
        em_s <- tryCatch(
          emmeans::emmeans(model,
                           as.formula(paste("~", focal_spec, "|", stratifier))),
          error = function(e) NULL)
        if (is.null(em_s)) next
        ctr_s <- tryCatch(as.data.frame(pairs(em_s, adjust = "tukey", infer = c(TRUE, TRUE))),
                          error = function(e) NULL)
        if (!is.null(ctr_s)) {
          slot <- paste0(paste(focal, collapse = "_"), "_within_", stratifier)
          simple[[slot]] <- ctr_s
        }
      }
    }
  }

  list(emm = emm_tab, pairs = pair_tab, cld = cld_tab,
       simple = if (length(simple) > 0) simple else NULL)
}

# (bh_family() helper removed 2026-08-07: unused, and no correction is applied
# in this script -- see 00_shared/correction_methods.R if it's needed again.)

# Pseudoreplication probe: refit dropping each RE; compare term β & SE.
pseudorep_compare <- function(fixed_formula, re_winner_term, re_list, data, term_name) {
  rows <- list()
  # Winner row
  win_form <- update(fixed_formula, paste(". ~ . +", re_winner_term))
  m_win <- suppressMessages(suppressWarnings(lmerTest::lmer(win_form, data = data, REML = TRUE)))
  fe <- summary(m_win)$coefficients
  if (term_name %in% rownames(fe)) {
    rows[["winner"]] <- data.frame(
      spec = paste0("WINNER: ", re_winner_term),
      beta = fe[term_name, "Estimate"],
      SE   = fe[term_name, "Std. Error"],
      p    = fe[term_name, ncol(fe)]
    )
  }
  for (lbl in names(re_list)) {
    if (identical(re_list[[lbl]], re_winner_term)) next
    f <- update(fixed_formula, paste(". ~ . +", re_list[[lbl]]))
    m <- tryCatch(
      suppressMessages(suppressWarnings(lmerTest::lmer(f, data = data, REML = TRUE))),
      error = function(e) NULL)
    if (is.null(m)) next
    fe2 <- summary(m)$coefficients
    if (!(term_name %in% rownames(fe2))) next
    rows[[lbl]] <- data.frame(
      spec = lbl, beta = fe2[term_name, "Estimate"],
      SE   = fe2[term_name, "Std. Error"], p = fe2[term_name, ncol(fe2)]
    )
  }
  out <- do.call(rbind, rows)
  if (is.null(out) || nrow(out) == 0) return(NULL)
  out$dSE_pct  <- (out$SE - out$SE[1]) / out$SE[1] * 100
  out$dBeta    <- out$beta - out$beta[1]
  out$flag     <- abs(out$dSE_pct) > 30 | sign(out$beta) != sign(out$beta[1])
  rownames(out) <- NULL
  out
}

# ---- 3. Read worksheets -----------------------------------------------------
b3_cort_raw <- read_excel(DATA_XLSX, sheet = "plasma cortisol B3", .name_repair = "minimal")
b3_samp_raw <- read_excel(DATA_XLSX, sheet = "batch 3 sampling",   .name_repair = "minimal")
mono_raw    <- read_excel(MONO_XLSX, sheet = "Sheet3",              .name_repair = "minimal")

# Save column dictionary
col_dict <- bind_rows(
  tibble(source = "plasma cortisol B3", original = names(b3_cort_raw)),
  tibble(source = "batch 3 sampling",   original = names(b3_samp_raw)),
  tibble(source = "monoamine Sheet3",   original = names(mono_raw))
)
write.csv(col_dict, file.path(OUT, "column_dictionary.csv"), row.names = FALSE)

# Clean B3 cortisol
b3_cort <- b3_cort_raw %>%
  janitor::clean_names() %>%
  select(any_of(c("sample_id","md","batch","wt","lenght","sex","video_id","trial",
                  "tank","side","density","condition","date","plasma_cortisol",
                  "plate","plate_date","microdissecred"))) %>%
  mutate(sample_id = pad_b3(sample_id),
         plate     = as.character(plate),
         plate_date = suppressWarnings(as.Date(plate_date)),
         sex       = as.character(sex),
         condition = as.character(condition),
         plasma_cortisol = to_num(plasma_cortisol),
         wt        = to_num(wt),
         lenght    = to_num(lenght),
         density   = to_num(density))

# Clean batch 3 sampling
b3_samp <- b3_samp_raw %>%
  janitor::clean_names() %>%
  filter(!is.na(sample_id)) %>%
  mutate(sample_n  = suppressWarnings(as.integer(sample_id)),
         sample_id = sprintf("B3_%02d", sample_n))

# Clean monoamines Sheet3
# Detect duplicate 'area' columns; keep the first
mono_names <- names(mono_raw)
keep_cols  <- c("area","sample","5-HT","5-HIAA","5-HIAA/5-HT","DA","DOPAC","DOPAC/DA",
                "NE","mg","NOTES","treatment","date")
mono <- mono_raw[, intersect(keep_cols, mono_names), drop = FALSE]
names(mono) <- janitor::make_clean_names(names(mono))
# robust column references
mono <- mono %>%
  rename_with(~ sub("^x5_ht$",          "ht_5",         .x)) %>%
  rename_with(~ sub("^x5_hiaa$",        "hiaa_5",       .x)) %>%
  rename_with(~ sub("^x5_hiaa_5_ht$",   "hiaa_5_ratio", .x)) %>%
  rename_with(~ sub("^dopac_da$",       "dopac_da_ratio", .x))
mono <- mono %>%
  filter(!is.na(sample) & !is.na(area)) %>%
  mutate(sample_n  = suppressWarnings(as.integer(sample)),
         sample_id = sprintf("B3_%02d", sample_n),
         area      = toupper(trimws(as.character(area))))

# Warn if area values do not match expected set
expected_areas <- c("DM","POA","VV","VD")
observed_areas <- sort(unique(mono$area))
if (!all(observed_areas %in% expected_areas)) {
  warning(sprintf("Monoamine areas observed = {%s}; expected = {%s}",
                  paste(observed_areas, collapse=", "),
                  paste(expected_areas, collapse=", ")))
}

# ---- 4. Parse MyAssays sheets into long well-level data ---------------------
parse_myassays <- function(plate_id, file, sheet, plate_date) {
  cat(sprintf("Parsing plate %s -> %s :: %s\n", plate_id, basename(file), sheet))
  ws <- suppressMessages(read_excel(file, sheet = sheet, col_names = FALSE,
                                    .name_repair = "minimal"))
  ws <- as.data.frame(ws, stringsAsFactors = FALSE)

  # Locate the two header rows (Calibrator / Sample) in column B (index 2)
  colB <- as.character(ws[[2]])
  cal_row <- which(colB %in% c("Calibrator"))[1]
  smp_row <- which(colB %in% c("Sample"))[1]
  if (is.na(cal_row) || is.na(smp_row)) {
    stop(sprintf("Header rows not found in plate %s", plate_id))
  }

  parse_block <- function(start_row, kind) {
    end_row <- nrow(ws)
    rows <- (start_row + 1):end_row
    out  <- list()
    cur_sample <- NA_character_
    cur_dilution <- NA_real_
    for (i in rows) {
      r <- ws[i, , drop = TRUE]
      if (all(is.na(r) | r == "")) next
      if (kind == "standard" && !is.na(r[[2]]) && r[[2]] == "Sample") break
      smp <- r[[2]]; dil <- r[[3]]; well <- r[[4]]; raw <- r[[5]]
      pct <- r[[6]]; conc <- r[[7]]; conc_avg <- r[[8]]; my_cv <- r[[9]]
      if (!is.na(smp) && nzchar(trimws(as.character(smp)))) {
        cur_sample   <- trimws(as.character(smp))
        cur_dilution <- to_num(dil)
      }
      if (is.na(well) || !nzchar(trimws(as.character(well)))) next
      out[[length(out)+1]] <- tibble(
        plate         = plate_id,
        plate_date    = plate_date,
        sample_raw    = cur_sample,
        dilution      = cur_dilution,
        well          = trimws(as.character(well)),
        raw_OD        = to_num(raw),
        pct_B_B0      = to_num(pct),
        conc_pgmL_str = as.character(conc),
        conc_avg_str  = as.character(conc_avg),
        myassays_pct_cv_str = as.character(my_cv),
        kind          = kind
      )
    }
    bind_rows(out)
  }

  std_block <- parse_block(cal_row, kind = "standard")
  smp_block <- parse_block(smp_row, kind = "sample")
  out <- bind_rows(std_block, smp_block)

  out <- out %>%
    mutate(
      out_of_range_high = grepl("> *Curve", conc_pgmL_str, ignore.case = TRUE),
      out_of_range_low  = grepl("< *Curve", conc_pgmL_str, ignore.case = TRUE),
      conc_pgmL         = to_num(conc_pgmL_str),
      conc_avg_pgmL     = to_num(conc_avg_str),
      myassays_pct_cv   = to_num(myassays_pct_cv_str),
      is_standard       = kind == "standard" | grepl("^Standard", sample_raw),
      is_b0             = sample_raw == "B0",
      is_nsb            = sample_raw == "NSB",
      sample_id         = ifelse(is_standard | is_b0 | is_nsb, sample_raw,
                                 pad_b3(sample_raw))
    ) %>%
    group_by(plate, sample_raw) %>%
    mutate(replicate_index = row_number()) %>%
    ungroup()
  out
}

assay_long <- purrr::imap_dfr(PLATE_FILES, function(spec, plate_id)
  parse_myassays(plate_id, spec$file, spec$sheet, spec$plate_date))

write.csv(assay_long, file.path(OUT, "assay_long.csv"), row.names = FALSE)

# ---- 5. Build sample-level summary (assay_wide) -----------------------------
assay_wide <- assay_long %>%
  filter(!is_standard & !is_b0 & !is_nsb) %>%
  group_by(plate, plate_date, sample_id) %>%
  summarise(
    n_replicates           = sum(!is.na(conc_pgmL)),
    mean_conc              = mean(conc_pgmL, na.rm = TRUE),
    sd_conc                = sd(conc_pgmL, na.rm = TRUE),
    cv_intra_pct_recomputed = ifelse(mean_conc > 0 & n_replicates >= 2,
                                     (sd_conc / mean_conc) * 100, NA_real_),
    myassays_pct_cv        = suppressWarnings(mean(myassays_pct_cv, na.rm = TRUE)),
    n_out_of_range         = sum(out_of_range_high | out_of_range_low),
    .groups = "drop"
  ) %>%
  mutate(
    cv_source        = ifelse(is.finite(myassays_pct_cv), "myassays", "recomputed"),
    cv_intra_pct     = ifelse(is.finite(myassays_pct_cv), myassays_pct_cv, cv_intra_pct_recomputed)
  )

# ---- 6. Matching audit -------------------------------------------------------
matching_audit <- b3_cort %>%
  select(sample_id, worksheet_plate = plate, worksheet_concentration = plasma_cortisol) %>%
  full_join(
    assay_wide %>%
      filter(grepl("^B3_", sample_id)) %>%
      select(sample_id, assay_plate = plate, derived_concentration_pgmL = mean_conc),
    by = "sample_id"
  ) %>%
  mutate(
    derived_concentration = derived_concentration_pgmL / 1000,
    concentration_diff    = worksheet_concentration - derived_concentration,
    match_status = case_when(
      is.na(worksheet_plate)               ~ "missing_in_worksheet",
      is.na(assay_plate)                   ~ "missing_in_assay",
      worksheet_plate == assay_plate       ~ "match",
      worksheet_plate != assay_plate       ~ "plate_mismatch",
      TRUE                                 ~ "unknown"
    )
  )

# ---- 7. Intra-assay CV tables -----------------------------------------------
intra_assay_cv <- assay_wide %>%
  arrange(plate, sample_id)

intra_assay_summary <- intra_assay_cv %>%
  filter(is.finite(cv_intra_pct)) %>%
  group_by(plate) %>%
  summarise(
    n_samples       = n(),
    median_cv       = median(cv_intra_pct),
    weighted_mean_cv = sum(cv_intra_pct * n_replicates) / sum(n_replicates),
    iqr_cv          = IQR(cv_intra_pct),
    pct_gt_15       = mean(cv_intra_pct > 15) * 100,
    pct_gt_20       = mean(cv_intra_pct > 20) * 100,
    .groups = "drop"
  )

# ---- 8. Inter-assay CV ------------------------------------------------------
# 8a. B0 / NSB controls across plates
ctrl_long <- assay_long %>% filter(is_b0 | is_nsb) %>%
  mutate(control = ifelse(is_b0, "B0", "NSB"))
inter_assay_cv_controls <- ctrl_long %>%
  group_by(control) %>%
  summarise(
    n_plates  = n_distinct(plate),
    n_wells   = sum(!is.na(raw_OD)),
    mean_OD   = mean(raw_OD, na.rm = TRUE),
    sd_OD     = sd(raw_OD, na.rm = TRUE),
    cv_inter_pct_OD = ifelse(mean_OD != 0, sd_OD / mean_OD * 100, NA_real_),
    mean_pgmL = mean(conc_pgmL, na.rm = TRUE),
    sd_pgmL   = sd(conc_pgmL, na.rm = TRUE),
    cv_inter_pct_pgmL = ifelse(is.finite(mean_pgmL) & mean_pgmL > 0,
                               sd_pgmL / mean_pgmL * 100, NA_real_),
    .groups = "drop"
  )

# 8b. Standards across plates
std_long <- assay_long %>% filter(is_standard) %>%
  mutate(nominal = NOMINAL_STD[sample_raw])
inter_assay_cv_standards <- std_long %>%
  filter(!is.na(nominal) & !is.na(conc_pgmL)) %>%
  group_by(sample_raw, nominal) %>%
  summarise(
    n_plates = n_distinct(plate),
    n_wells  = n(),
    mean_backfit = mean(conc_pgmL),
    sd_backfit   = sd(conc_pgmL),
    cv_inter_pct = ifelse(mean_backfit > 0, sd_backfit / mean_backfit * 100, NA_real_),
    mean_recovery_pct = mean(conc_pgmL / nominal * 100),
    .groups = "drop"
  )

# 8c. B3 unknowns repeated across plates
b3_repeat <- assay_wide %>% filter(grepl("^B3_", sample_id)) %>%
  group_by(sample_id) %>% filter(n_distinct(plate) >= 2) %>% ungroup()
inter_assay_cv_repeats <- b3_repeat %>%
  group_by(sample_id) %>%
  summarise(
    n_plates = n_distinct(plate),
    plates   = paste(sort(unique(plate)), collapse = ","),
    mean_conc_pgmL = mean(mean_conc, na.rm = TRUE),
    sd_conc_pgmL   = sd(mean_conc, na.rm = TRUE),
    cv_inter_pct   = ifelse(mean_conc_pgmL > 0,
                            sd_conc_pgmL / mean_conc_pgmL * 100, NA_real_),
    .groups = "drop"
  )

# 8d. Per-plate standard curve quality
plate_curve_quality <- purrr::imap_dfr(PLATE_FILES, function(spec, plate_id) {
  ws <- suppressMessages(read_excel(spec$file, sheet = spec$sheet, col_names = FALSE,
                                    .name_repair = "minimal"))
  ws <- as.data.frame(ws, stringsAsFactors = FALSE)
  pick <- function(label) {
    idx <- which(as.character(ws[[10]]) == label)
    if (length(idx) == 0) return(NA_real_)
    to_num(ws[[11]][idx[1]])
  }
  tibble(plate = plate_id,
         a   = pick("a"), b = pick("b"), c = pick("c"), d = pick("d"),
         MSE = pick("MSE"), R2 = pick("R²"), SS = pick("SS"), SYX = pick("SYX"))
})
if (all(is.na(plate_curve_quality$R2))) {
  plate_curve_quality$R2 <- purrr::imap_dbl(PLATE_FILES, function(spec, plate_id) {
    ws <- suppressMessages(read_excel(spec$file, sheet = spec$sheet, col_names = FALSE,
                                      .name_repair = "minimal"))
    ws <- as.data.frame(ws, stringsAsFactors = FALSE)
    idx <- grep("^R", as.character(ws[[10]]))
    if (length(idx) == 0) NA_real_ else to_num(ws[[11]][idx[1]])
  })
}

# 8e. Replicate well summary
replicate_well_summary <- assay_long %>%
  group_by(plate) %>%
  summarise(
    n_wells_total      = n(),
    n_wells_standards  = sum(is_standard),
    n_wells_b0         = sum(is_b0),
    n_wells_nsb        = sum(is_nsb),
    n_wells_unknown    = sum(!is_standard & !is_b0 & !is_nsb),
    n_out_of_range     = sum(out_of_range_high | out_of_range_low, na.rm = TRUE),
    .groups = "drop"
  )

# Write all CV tables to one workbook
write_xlsx(list(
  intra_assay_cv         = intra_assay_cv,
  intra_assay_summary    = intra_assay_summary,
  inter_assay_cv_controls= inter_assay_cv_controls,
  inter_assay_cv_standards = inter_assay_cv_standards,
  inter_assay_cv_repeats = inter_assay_cv_repeats,
  plate_curve_quality    = plate_curve_quality,
  replicate_well_summary = replicate_well_summary,
  matching_audit         = matching_audit
), file.path(OUT, "assay_cv_tables.xlsx"))

# ---- 9. Build clean cortisol dataset and exports ----------------------------
b3_cort_clean <- b3_cort %>%
  left_join(assay_wide %>% select(sample_id, plate, mean_conc, cv_intra_pct, cv_source),
            by = c("sample_id","plate"))

# ---- Trial-ID repair (2026-08-08, df_method_memo section 6 check) ----------
# The raw `trial` column does NOT identify the arena trial: e.g. raw trial 1
# holds 6 control + 5 treat fish (11 total), but a real trial is 5 fish of ONE
# condition. This was caught by the memo's own diagnostic ("confirm anova()
# reports ~14 df for a between-trial contrast; if it reports 70+, the trial RE
# is missing") -- the cortisol model returned df2 = 50.19, neither figure,
# which is what happens when the grouping factor mixes both arms.
#
# Verified: `tank x density` reproduces the documented 4x4 grid exactly -- 16
# cells, 5 fish each, a SINGLE condition each (8 control / 8 treat trials),
# matching Table S1. The raw column is kept as `trial_raw` for audit; `trial`
# itself is overwritten with the corrected grouping so every existing
# downstream join/model that already reads `trial` (cortisol RE, the
# monoamine sample_id join, sex-balance checks) inherits the fix without
# needing to be individually rewritten.
b3_cort_clean <- b3_cort_clean %>%
  mutate(trial_raw = trial,
         trial     = factor(paste(tank, density, sep = "_")))

# Regression guard: every trial must contain exactly one condition. Fail loud
# rather than silently reporting inference on a mixed-condition cluster again.
.trial_integrity <- b3_cort_clean %>%
  filter(!is.na(trial), !is.na(condition)) %>%
  dplyr::distinct(trial, condition) %>%
  dplyr::count(trial) %>%
  dplyr::filter(n > 1)
if (nrow(.trial_integrity) > 0)
  stop("Trial-ID integrity check failed: ", nrow(.trial_integrity),
       " trial(s) span more than one condition: ",
       paste(.trial_integrity$trial, collapse = ", "))
message(sprintf("Trial-ID integrity OK: %d trials (tank x density), 1 condition each.",
                length(unique(b3_cort_clean$trial))))

write.csv(b3_cort_clean, file.path(OUT, "b3_cortisol_clean.csv"), row.names = FALSE)

# ---- 10. Section §3 + §4 cortisol -----------------------------------------
# §3 = sex-balance audit (full cortisol cohort + monoamine subsample)
# §4 = treatment main effect on cortisol (Part 1, template-style block)
# ---------------------------------------------------------------------------

# Output directories per plan
DATA_DIR  <- file.path(OUT, "Data");        dir.create(DATA_DIR,  showWarnings = FALSE, recursive = TRUE)
FIG_DIR   <- file.path(OUT, "Figures");     dir.create(FIG_DIR,   showWarnings = FALSE, recursive = TRUE)
FINAL_DIR <- file.path(OUT, "B3_Final_Report"); dir.create(FINAL_DIR, showWarnings = FALSE, recursive = TRUE)

# Cleaned cortisol working data — ALL fish (sex used for §3 + §9 only)
# M4 (R2-5): order_idx = chronological rank of testing date (1..8), the trial-order/
# day covariate. Derived from `date` (the behavioural testing day, shared with the
# monoamine sheet and the behaviour pipeline's trial_date).
cort_full <- b3_cort_clean %>%
  filter(!is.na(plasma_cortisol)) %>%
  mutate(
    log_cort  = log(plasma_cortisol),
    order_idx = make_order_idx(date),
    trial     = factor(trial),
    tank      = factor(tank),
    plate     = factor(plate),
    condition = factor(condition, levels = c("control","treat")),
    sex       = factor(sex, levels = c("F","M"))
  )

# Monoamine subsample identifier (sex-balance check § b uses this)
mono_subsample_ids <- unique(mono$sample_id)

cort_full <- set_sum_contrasts(cort_full, c("condition","sex"))

# `cort_sex` retained for §14b legacy plots and §9 secondary analyses
cort_sex <- cort_full %>% filter(sex %in% c("M","F"))

# ---- §3 SEX-BALANCE CHECK -------------------------------------------------
# Test whether sex composition differs across treatment groups.
# Cohort (a): full cortisol n=80 ; Cohort (b): monoamine subsample n≈37.
sex_balance_run <- function(df, cohort_label) {
  df <- df %>% filter(!is.na(sex), sex %in% c("F","M"), !is.na(condition)) %>%
    droplevels()
  if (nrow(df) < 5 || length(unique(df$sex)) < 2 || length(unique(df$condition)) < 2)
    return(data.frame(cohort = cohort_label, test = "skipped",
                      statistic = NA, df = NA, p = NA, interpretation = "insufficient data"))
  out <- list()
  # 1. Overall χ² / Fisher
  tab <- table(df$sex, df$condition)
  use_fisher <- any(chisq.test(tab)$expected < 5)
  if (use_fisher) {
    fr <- fisher.test(tab)
    out[[1]] <- data.frame(cohort = cohort_label, test = "Fisher exact (overall)",
                           statistic = NA, df = NA, p = fr$p.value,
                           interpretation = ifelse(fr$p.value < 0.05,
                             "imbalance detected", "balanced"))
  } else {
    cr <- chisq.test(tab)
    out[[1]] <- data.frame(cohort = cohort_label, test = "Chi-sq (overall)",
                           statistic = unname(cr$statistic), df = unname(cr$parameter),
                           p = cr$p.value,
                           interpretation = ifelse(cr$p.value < 0.05,
                             "imbalance detected", "balanced"))
  }
  # 2. Cochran-Mantel-Haenszel stratified by trial
  if (length(unique(df$trial)) >= 2) {
    cmh <- tryCatch(mantelhaen.test(df$sex, df$condition, df$trial), error = function(e) NULL)
    if (!is.null(cmh)) {
      out[[2]] <- data.frame(cohort = cohort_label, test = "CMH (stratified by trial)",
                             statistic = unname(cmh$statistic),
                             df = unname(cmh$parameter)[1],
                             p = cmh$p.value,
                             interpretation = ifelse(cmh$p.value < 0.05,
                               "imbalance after trial stratification",
                               "balanced after trial stratification"))
    }
  }
  # 3. GLMM: P(M) ~ condition + (1|trial)
  m_sb <- tryCatch(suppressMessages(suppressWarnings(
    lme4::glmer(I(sex == "M") ~ condition + (1|trial), data = df, family = binomial))),
    error = function(e) NULL)
  if (is.null(m_sb) || lme4::isSingular(m_sb, tol = 1e-4)) {
    m_sb_g <- tryCatch(glm(I(sex == "M") ~ condition + factor(trial), data = df, family = binomial),
                       error = function(e) NULL)
    if (!is.null(m_sb_g)) {
      a <- car::Anova(m_sb_g, type = 2)
      idx <- grep("^condition$", rownames(a))[1]
      if (!is.na(idx)) {
        out[[3]] <- data.frame(cohort = cohort_label, test = "GLM (trial fixed; lmm singular)",
                               statistic = a$`LR Chisq`[idx], df = a$Df[idx],
                               p = a$`Pr(>Chisq)`[idx],
                               interpretation = ifelse(a$`Pr(>Chisq)`[idx] < 0.05,
                                 "treatment shifts sex prob.",
                                 "no treatment shift in sex prob."))
      }
    }
  } else {
    a <- car::Anova(m_sb, type = 2)
    idx <- grep("^condition$", rownames(a))[1]
    if (!is.na(idx)) {
      out[[3]] <- data.frame(cohort = cohort_label, test = "GLMM I(sex=M)~condition+(1|trial)",
                             statistic = a$Chisq[idx], df = a$Df[idx],
                             p = a$`Pr(>Chisq)`[idx],
                             interpretation = ifelse(a$`Pr(>Chisq)`[idx] < 0.05,
                               "treatment shifts sex prob.",
                               "no treatment shift in sex prob."))
    }
  }
  do.call(rbind, out)
}

sex_balance_full   <- sex_balance_run(cort_full, "full cortisol cohort (n=80)")
sex_balance_subset <- sex_balance_run(
  cort_full %>% filter(sample_id %in% mono_subsample_ids),
  sprintf("monoamine subsample (n=%d)", sum(unique(cort_full$sample_id) %in% mono_subsample_ids)))
sex_balance <- rbind(sex_balance_full, sex_balance_subset)
write.csv(sex_balance, file.path(DATA_DIR, "sex_balance_test.csv"), row.names = FALSE)
message(sprintf("§3 sex-balance: %d tests across %d cohorts; min p = %s.",
                nrow(sex_balance), 2, fmt_p(min(sex_balance$p, na.rm = TRUE))))

# ---- §4 CORTISOL: Part 1 (Treatment Main Effect) -------------------------
cort_dat <- cort_full %>% filter(is.finite(log_cort))

# 4.x.1 Normality / Levene
cort_norm  <- transform_check(cort_dat$plasma_cortisol, label = "plasma_cortisol")
cort_lev   <- levene_check(plasma_cortisol ~ condition, data = cort_dat,
                            label = "plasma_cortisol ~ condition")
norm_lev_cort <- merge(cort_norm, cort_lev, by = "outcome", all = TRUE)

# 4.x.2 RE selection via AICc
fixed_cort_pri <- log_cort ~ condition
re_cort_cand <- list(
  "(1|tank)"            = "(1|tank)",
  "(1|trial)"           = "(1|trial)",
  "(1|plate)"           = "(1|plate)",
  "(1|tank)+(1|trial)"  = "(1|tank) + (1|trial)",
  "(1|tank)+(1|plate)"  = "(1|tank) + (1|plate)",
  "(1|trial)+(1|plate)" = "(1|trial) + (1|plate)"
)
# (1|trial) is forced: treatment is applied to the arena trial, cortisol is
# measured on the individual fish, so trial is the unit of replication. AICc
# now only decides whether the optional assay/housing nuisance terms help.
sel_cort <- select_re_aicc(fixed_cort_pri, re_cort_cand, data = cort_dat,
                           force_re = "(1|trial)")
cort_re_aicc_table <- sel_cort$table %>% mutate(outcome = "plasma_cortisol", .before = 1)

if (!is.null(sel_cort$winner)) {
  fit_cort <- sel_cort$winner
  cort_re_label   <- sel_cort$label
  # select_re_aicc's OWN winning formula string -- NOT re-derived via
  # re_cort_cand[[cort_re_label]] downstream. force_re can mint a combo label
  # (e.g. "(1|plate)+(1|trial)" from a lone "(1|plate)" candidate) that need
  # not exist as a literal key in re_cort_cand (which spells the same
  # combination "(1|trial)+(1|plate)"), and indexing with it there silently
  # returned NULL -- this crashed pseudorep_compare() 2026-08-08 the first
  # time the corrected trial grouping changed which candidate actually won.
  cort_re_formula <- sel_cort$re_formula
  cort_model_type <- paste("lmer:", cort_re_label)
  message(sprintf("§4 cortisol primary model: %s", cort_model_type))
} else {
  message("§4 all RE candidates singular — falling back to lm")
  fit_cort <- lm(log_cort ~ condition, data = cort_dat)
  cort_re_label   <- "lm (no RE)"
  cort_re_formula <- NULL
  cort_model_type <- "lm (no random effects)"
}

# 4.x.3 R²
cort_r2 <- r2_row(fit_cort) %>% mutate(outcome = "plasma_cortisol", .before = 1)

# 4.x.4 KR ANOVA (primary)
cort_anova_kr <- tryCatch({
  a <- as.data.frame(anova_kr(fit_cort, type = 2))  # type 2 = type 3 since 1 predictor
  a %>% tibble::rownames_to_column("term") %>% mutate(outcome = "plasma_cortisol", .before = 1)
}, error = function(e) NULL)

# 4.x.5 Satterthwaite ANOVA (sensitivity)
cort_anova_satt <- tryCatch({
  a <- as.data.frame(anova_satt(fit_cort, type = 2))
  a %>% tibble::rownames_to_column("term") %>% mutate(outcome = "plasma_cortisol", .before = 1)
}, error = function(e) NULL)

# 4.x.6 EMM + Tukey + CLD
cort_emm_list <- emm_tukey_cld(fit_cort, ~ condition)
if (!is.null(cort_emm_list)) {
  cort_emm_df <- cort_emm_list$emm    %>% mutate(outcome = "plasma_cortisol", .before = 1)
  cort_pairs_df <- cort_emm_list$pairs %>% mutate(outcome = "plasma_cortisol", .before = 1)
  cort_cld_df <- cort_emm_list$cld    %>% mutate(outcome = "plasma_cortisol", .before = 1)
} else {
  cort_emm_df <- cort_pairs_df <- cort_cld_df <- NULL
}

# 4.x.7 DHARMa diagnostic
diag_cort_dharma_path <- file.path(FIG_DIR, "diag_cortisol_dharma.png")
dharma_panel(fit_cort, diag_cort_dharma_path,
             title = paste("Cortisol DHARMa —", cort_model_type))

# Pull treatment β + SE for cross-summary §6
cort_summary <- summary(fit_cort)
cort_fe <- cort_summary$coefficients
cort_treat_idx <- grep("^condition", rownames(cort_fe))[1]
cort_treat_beta <- if (!is.na(cort_treat_idx)) cort_fe[cort_treat_idx, "Estimate"] else NA_real_
cort_treat_se   <- if (!is.na(cort_treat_idx)) cort_fe[cort_treat_idx, "Std. Error"] else NA_real_

# Save §4-cortisol exports
write.csv(cort_re_aicc_table, file.path(DATA_DIR, "re_selection_aicc_cortisol.csv"), row.names = FALSE)
write.csv(cort_r2,            file.path(DATA_DIR, "r2_icc_cortisol.csv"),             row.names = FALSE)
write.csv(cort_anova_kr,      file.path(DATA_DIR, "anova_kr_cortisol.csv"),           row.names = FALSE)
write.csv(cort_anova_satt,    file.path(DATA_DIR, "anova_satt_cortisol.csv"),         row.names = FALSE)
if (!is.null(cort_emm_df))   write.csv(cort_emm_df,   file.path(DATA_DIR, "emm_cortisol.csv"),   row.names = FALSE)
if (!is.null(cort_pairs_df)) write.csv(cort_pairs_df, file.path(DATA_DIR, "pairs_cortisol.csv"), row.names = FALSE)
if (!is.null(cort_cld_df))   write.csv(cort_cld_df,   file.path(DATA_DIR, "cld_cortisol.csv"),   row.names = FALSE)

# Legacy outputs preserved (used elsewhere)
write.csv(cort_dat %>% group_by(condition) %>%
            summarise(n = n(),
                      mean = mean(plasma_cortisol, na.rm = TRUE),
                      sd   = sd(plasma_cortisol,   na.rm = TRUE),
                      median = median(plasma_cortisol, na.rm = TRUE),
                      q1 = quantile(plasma_cortisol, 0.25, na.rm = TRUE),
                      q3 = quantile(plasma_cortisol, 0.75, na.rm = TRUE),
                      .groups = "drop"),
          file.path(DATA_DIR, "descriptive_cortisol.csv"), row.names = FALSE)

# Legacy violin plot kept (used by Word report appendix)
p_cort <- ggplot(cort_dat, aes(condition, plasma_cortisol, fill = condition)) +
  geom_violin(alpha = 0.4) +
  geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.7) +
  labs(title = "B3 plasma cortisol by treatment",
       y = "Plasma cortisol (ng/mL)", x = NULL) +
  theme_bw() + theme(legend.position = "none")
ggsave(file.path(FIG_DIR, "fig_cortisol_treatment.png"), p_cort, width = 7, height = 5, dpi = 200)
ggsave(file.path(FIG_DIR, "fig_cortisol_treatment.pdf"), p_cort, width = 7, height = 5)

# ---- 11. Section §4 + §5 monoamines ----------------------------------------
# NE (noradrenaline) is excluded from all analyses by author decision (2026-08-04):
# 6 analytes x 4 brain areas = 24 analyte x region cells, not 28.
analyte_cols  <- c("ht_5","hiaa_5","da","dopac","hiaa_5_ratio","dopac_da_ratio")
analyte_label <- c(ht_5 = "5-HT", hiaa_5 = "5-HIAA", da = "DA", dopac = "DOPAC",
                   hiaa_5_ratio = "5-HIAA/5-HT", dopac_da_ratio = "DOPAC/DA")

# No multiple-comparison correction is applied (author decision, 2026-08-07):
# every treatment-effect p-value below is the raw, unadjusted p from its own
# model. The BH/Holm family-based correction previously computed here (and the
# 7-procedure sensitivity comparison) has been moved, unmodified, to
# 00_shared/correction_methods.R -- it is no longer sourced by this script.

for (a in analyte_cols) {
  if (a %in% names(mono)) mono[[a]] <- to_num(mono[[a]])
}

# Join condition/tank/trial/sex/cortisol from b3_cort_clean (sex used only for §9)
mono_join <- mono %>%
  left_join(b3_cort_clean %>% select(sample_id, sex, condition, tank, trial, plasma_cortisol),
            by = "sample_id") %>%
  mutate(
    sex       = factor(sex, levels = c("F","M")),
    condition = factor(condition, levels = c("control","treat")),
    tank      = factor(tank),
    trial     = factor(trial),
    area      = factor(area, levels = expected_areas),
    # M4 (R2-5): order_idx from the monoamine sheet's own `date` (same testing-day
    # domain as cortisol/behaviour — see ALL_TRIAL_DATES).
    order_idx = make_order_idx(date)
  )

# Treatment column harmonisation (mono sheet has 'treatment' col; cortisol has 'condition')
if (!"treatment" %in% names(mono_join))  mono_join$treatment <- mono_join$condition
mono_join$treatment <- factor(mono_join$treatment, levels = c("control","treat"))
mono_join <- set_sum_contrasts(mono_join, c("treatment","area"))

mono_long <- mono_join %>%
  pivot_longer(any_of(analyte_cols), names_to = "analyte_key", values_to = "value") %>%
  mutate(analyte = analyte_label[analyte_key])
write.csv(mono_long, file.path(DATA_DIR, "monoamine_long.csv"), row.names = FALSE)

# §8 descriptive (treatment × area)
mono_descr <- mono_long %>%
  filter(is.finite(value)) %>%
  group_by(analyte, area, treatment) %>%
  summarise(n = n(),
            mean = mean(value), sd = sd(value), sem = sd/sqrt(pmax(n,1)),
            median = median(value),
            q1 = quantile(value, 0.25), q3 = quantile(value, 0.75),
            .groups = "drop")
write.csv(mono_descr, file.path(DATA_DIR, "descriptive_monoamine.csv"), row.names = FALSE)

# Per-analyte fit pipeline (§4 + §5 template-style)
fit_one_mono <- function(akey) {
  d <- mono_long %>% filter(analyte_key == akey, is.finite(value), value > 0) %>%
    droplevels()
  if (nrow(d) < 10 || length(unique(d$treatment)) < 2 || length(unique(d$area)) < 2)
    return(NULL)
  d$log_v <- log(d$value)
  d <- set_sum_contrasts(d, c("treatment","area"))

  # ----- 4.x.1 Normality + Levene
  norm_row <- transform_check(d$value, label = analyte_label[akey])
  lev_row  <- levene_check(value ~ treatment * area, data = d,
                            label = analyte_label[akey])
  norm_lev <- merge(norm_row, lev_row, by = "outcome", all = TRUE)

  # ----- 4.x.2 RE-AICc selection
  # M2 (pseudoreplication, R1-6): `trial` is now FORCED, not merely offered as a
  # candidate. Multiple fish are drawn from the same trial (video), so (1|trial)
  # is the correct higher-level nesting above (1|sample_id). Before 2026-08-08 it
  # was a competing candidate and AICc never selected it for any analyte.
  re_cand <- list(
    "(1|sample_id)"            = "(1|sample_id)",
    "(1|sample_id)+(1|date)"   = "(1|sample_id) + (1|date)",
    "(1|sample_id)+(1|tank)"   = "(1|sample_id) + (1|tank)"
  )
  fixed_form <- log_v ~ treatment * area
  sel <- select_re_aicc(fixed_form, re_cand, data = d, force_re = "(1|trial)")
  if (is.null(sel$winner)) {
    # fallback keeps BOTH design terms: fish within trial.
    fit <- tryCatch(lmerTest::lmer(log_v ~ treatment * area + (1|sample_id) + (1|trial),
                                   data = d, REML = TRUE),
                    error = function(e) NULL)
    re_label   <- "(1|sample_id) + (1|trial) [forced]"
    re_formula <- "(1|sample_id) + (1|trial)"
  } else {
    fit        <- sel$winner
    re_label   <- sel$label
    # Use select_re_aicc's OWN winning formula string, not a re-derivation via
    # re_label indexed into a caller-side list -- force_re can mint a combo
    # label (e.g. "(1|tank)+(1|trial)") that need not exist as a literal key
    # in any list outside this function.
    re_formula <- sel$re_formula
  }
  if (is.null(fit)) return(NULL)
  re_aicc_table <- sel$table %>% mutate(outcome = analyte_label[akey], .before = 1)

  # ----- 4.x.3 R²
  r2 <- r2_row(fit) %>% mutate(outcome = analyte_label[akey], .before = 1)

  # ----- 4.x.4 KR ANOVA
  kr <- tryCatch(as.data.frame(anova_kr(fit, type = 3)) %>%
                   tibble::rownames_to_column("term") %>%
                   mutate(outcome = analyte_label[akey], .before = 1),
                 error = function(e) NULL)

  # ----- 4.x.5 Satterthwaite ANOVA
  satt <- tryCatch(as.data.frame(anova_satt(fit, type = 3)) %>%
                     tibble::rownames_to_column("term") %>%
                     mutate(outcome = analyte_label[akey], .before = 1),
                   error = function(e) NULL)

  # ----- 4.x.6 EMM aggregated (~treatment) and stratified (~treatment|area)
  emm_agg <- emm_tukey_cld(fit, ~ treatment)
  emm_str <- emm_tukey_cld(fit, ~ treatment, by = "area")
  add_label <- function(x) if (!is.null(x))
    purrr::map(x, function(d) {
      if (is.data.frame(d)) d %>% mutate(outcome = analyte_label[akey], .before = 1)
      else d   # pass-through for non-data-frame slots (e.g. $simple list)
    })
  else NULL
  emm_agg <- add_label(emm_agg)
  emm_str <- add_label(emm_str)

  # ----- 4.x.7 DHARMa
  dharma_path <- file.path(FIG_DIR, sprintf("diag_monoamine_dharma_%s.png", akey))
  dharma_panel(fit, dharma_path,
               title = sprintf("DHARMa: %s [%s]", analyte_label[akey], re_label))

  # treatment β + SE for cross-summary §6 (use coefficient on contrast-coded `treatment1`)
  sm <- summary(fit)$coefficients
  treat_idx <- grep("^treatment", rownames(sm))[1]
  beta <- if (!is.na(treat_idx)) sm[treat_idx, "Estimate"] else NA_real_
  se   <- if (!is.na(treat_idx)) sm[treat_idx, "Std. Error"] else NA_real_

  # KR p for treatment main effect (aggregated)
  kr_p <- NA_real_; kr_F <- NA_real_; kr_df1 <- NA_real_; kr_df2 <- NA_real_
  if (!is.null(kr)) {
    idx <- grep("^treatment$", kr$term)[1]
    if (!is.na(idx)) {
      kr_F   <- kr[["F"]][idx]
      if (is.null(kr_F) || !is.finite(kr_F)) kr_F <- kr[["F value"]][idx]
      kr_p   <- kr[["Pr(>F)"]][idx]
      kr_df1 <- if ("Df" %in% names(kr)) kr[["Df"]][idx] else kr[["NumDF"]][idx]
      kr_df2 <- if ("Df.res" %in% names(kr)) kr[["Df.res"]][idx]
                 else if ("DenDF"  %in% names(kr)) kr[["DenDF"]][idx]
                 else if ("Den Df" %in% names(kr)) kr[["Den Df"]][idx]
                 else {
                   res_row <- kr[grepl("^Residual", kr$term, ignore.case = TRUE), ]
                   if (nrow(res_row) > 0) res_row[["Df"]][1] else NA_real_
                 }
    }
  }

  list(
    analyte_key   = akey,
    label         = analyte_label[akey],
    re_label      = re_label,
    re_formula    = re_formula,
    fit           = fit,
    n             = nrow(d),
    norm_lev      = norm_lev,
    re_aicc_table = re_aicc_table,
    r2            = r2,
    anova_kr      = kr,
    anova_satt    = satt,
    emm_agg       = emm_agg,
    emm_str       = emm_str,
    dharma_path   = dharma_path,
    treat_beta    = beta,
    treat_se      = se,
    treat_F       = kr_F,
    treat_df1     = kr_df1,
    treat_df2     = kr_df2,
    treat_p       = kr_p,
    tidy          = broom.mixed::tidy(fit, conf.int = TRUE)
  )
}

mono_models <- purrr::compact(purrr::map(analyte_cols, fit_one_mono))
names(mono_models) <- vapply(mono_models, function(x) x$analyte_key, character(1))

# Bind per-outcome tables
bind_named <- function(field) purrr::map_dfr(mono_models, function(x) x[[field]])

mono_re_aicc      <- bind_named("re_aicc_table")
mono_r2           <- bind_named("r2")
mono_anova_kr     <- bind_named("anova_kr")
mono_anova_satt   <- bind_named("anova_satt")
mono_norm_lev     <- bind_named("norm_lev")

mono_emm_agg <- purrr::map_dfr(mono_models, function(x) if (!is.null(x$emm_agg)) x$emm_agg$emm)
mono_pairs_agg <- purrr::map_dfr(mono_models, function(x) if (!is.null(x$emm_agg)) x$emm_agg$pairs)
mono_cld_agg <- purrr::map_dfr(mono_models, function(x) if (!is.null(x$emm_agg)) x$emm_agg$cld)
mono_emm_str <- purrr::map_dfr(mono_models, function(x) if (!is.null(x$emm_str)) x$emm_str$emm)
mono_pairs_str <- purrr::map_dfr(mono_models, function(x) if (!is.null(x$emm_str)) x$emm_str$pairs)
mono_cld_str <- purrr::map_dfr(mono_models, function(x) if (!is.null(x$emm_str)) x$emm_str$cld)

mono_re_summary <- purrr::map_dfr(mono_models, function(x)
  tibble(analyte_key = x$analyte_key, analyte = x$label, re_spec = x$re_label,
         singular = lme4::isSingular(x$fit, tol = 1e-4)))

mono_models_tidy <- purrr::map_dfr(mono_models, function(x)
  x$tidy %>% mutate(analyte_key = x$analyte_key, analyte = x$label, re_spec = x$re_label))

write.csv(mono_re_aicc,       file.path(DATA_DIR, "re_selection_aicc_monoamines.csv"), row.names = FALSE)
write.csv(mono_r2,            file.path(DATA_DIR, "r2_monoamines.csv"),                 row.names = FALSE)
write.csv(mono_anova_kr,      file.path(DATA_DIR, "anova_kr_monoamines.csv"),           row.names = FALSE)
write.csv(mono_anova_satt,    file.path(DATA_DIR, "anova_satt_monoamines.csv"),         row.names = FALSE)
write.csv(mono_norm_lev,      file.path(DATA_DIR, "normality_levene_monoamines.csv"),   row.names = FALSE)
write.csv(mono_emm_agg,       file.path(DATA_DIR, "emm_aggregated_monoamines.csv"),     row.names = FALSE)
write.csv(mono_pairs_agg,     file.path(DATA_DIR, "pairs_aggregated_monoamines.csv"),   row.names = FALSE)
write.csv(mono_cld_agg,       file.path(DATA_DIR, "cld_aggregated_monoamines.csv"),     row.names = FALSE)
write.csv(mono_emm_str,       file.path(DATA_DIR, "emm_stratified_monoamines.csv"),     row.names = FALSE)
write.csv(mono_pairs_str,     file.path(DATA_DIR, "pairs_stratified_monoamines.csv"),   row.names = FALSE)
write.csv(mono_cld_str,       file.path(DATA_DIR, "cld_stratified_monoamines.csv"),     row.names = FALSE)
write.csv(mono_re_summary,    file.path(DATA_DIR, "monoamine_model_re_summary.csv"),    row.names = FALSE)
write.csv(mono_models_tidy,   file.path(DATA_DIR, "monoamine_models_tidy.csv"),         row.names = FALSE)

# ---- §4 PRIMARY: Per-analyte × per-area models (treatment as sole fixed effect) ----
# Inferential results: log_v ~ treatment + RE for each of 7 analytes × 4 areas = 28 models.
# RE candidates: (1|date), (1|tank), (1|trial). sample_id cannot be RE here — one obs per
# fish per area-subset. If all RE candidates singular, falls back to lm; flagged in output.
# M2 (pseudoreplication, R1-6): `trial` added as a candidate — several fish are sampled
# from the same trial/video, so treating them as independent overstates precision unless
# trial is either a random intercept (here) or absorbed by the winning RE.

fit_one_mono_area <- function(akey, ar) {
  d <- mono_long %>%
    filter(analyte_key == akey, area == ar,
           is.finite(value), value > 0, !is.na(treatment)) %>%
    droplevels()
  n_total <- nrow(d)
  n_per_treat <- if (length(unique(d$treatment)) >= 2) min(table(d$treatment)) else 0L
  if (n_total < 6 || n_per_treat < 2) return(NULL)

  # Normality diagnostics on raw values → adaptive transform selection
  cell_lbl <- sprintf("%s [%s]", analyte_label[akey], ar)
  norm_row  <- transform_check(d$value, label = cell_lbl)
  trans     <- select_transform(norm_row)

  # Apply selected transform; stored as log_v so all downstream formulas unchanged
  d$log_v <- apply_transform(d$value, trans)
  d <- set_sum_contrasts(d, "treatment")

  # Levene (raw and post-transform)
  lev_raw_row <- tryCatch(
    levene_check(value ~ treatment, data = d, label = cell_lbl),
    error = function(e) data.frame(outcome = cell_lbl,
                                   levene_F = NA_real_, df1 = NA, df2 = NA, p = NA_real_))
  lev_log_row <- tryCatch(
    levene_check(log_v ~ treatment, data = d, label = cell_lbl),
    error = function(e) data.frame(outcome = cell_lbl,
                                   levene_F = NA_real_, df1 = NA, df2 = NA, p = NA_real_))

  # RE selection: (1|trial) is FORCED (M2/R1-6) because the fish sampled together
  # from the same video are not independent replicates of treatment; (1|date) and
  # (1|tank) remain optional nuisance terms decided by AICc. sample_id is 1 obs
  # per fish in this subset, so it carries no information here.
  # Before 2026-08-08 trial merely competed against date/tank and was selected in
  # only 9 of 24 cells, with 5 cells falling through to an lm with no clustering.
  fixed_form <- log_v ~ treatment
  re_cand <- list("(1|date)" = "(1|date)", "(1|tank)" = "(1|tank)")
  sel <- select_re_aicc(fixed_form, re_cand, data = d, force_re = "(1|trial)")

  if (!is.null(sel$winner)) {
    fit           <- sel$winner
    re_label      <- sel$label
    # select_re_aicc's OWN winning formula string -- see the analogous fix and
    # comment on the region-aggregated model above; re_label (a display label)
    # is not guaranteed to be indexable in any caller-side list.
    re_formula    <- sel$re_formula
    model_formula <- paste0(trans$label, " ~ treatment + ", re_label)
  } else {
    # Design term retained even here: fall back to (1|trial) alone, NOT to lm.
    fit <- tryCatch(lmerTest::lmer(log_v ~ treatment + (1|trial), data = d, REML = TRUE),
                    error = function(e) NULL)
    if (is.null(fit)) {
      fit           <- tryCatch(lm(log_v ~ treatment, data = d), error = function(e) NULL)
      re_label      <- "lm (no RE — (1|trial) could not be fitted)"
      re_formula    <- NULL
      model_formula <- paste0(trans$label, " ~ treatment [lm, no RE]")
      message(sprintf("[%s x %s] (1|trial) unfittable -> lm", analyte_label[akey], ar))
    } else {
      re_label      <- "(1|trial) [forced]"
      re_formula    <- "(1|trial)"
      model_formula <- paste0(trans$label, " ~ treatment + (1|trial)")
      message(sprintf("[%s x %s] AICc set empty -> forced (1|trial)",
                      analyte_label[akey], ar))
    }
  }
  if (is.null(fit)) return(NULL)

  re_aicc_table <- sel$table %>%
    mutate(analyte = analyte_label[akey], area = ar, .before = 1)

  r2 <- r2_row(fit) %>% mutate(analyte = analyte_label[akey], area = ar, .before = 1)

  kr <- tryCatch(
    as.data.frame(anova_kr(fit, type = 2)) %>%
      tibble::rownames_to_column("term") %>%
      mutate(analyte = analyte_label[akey], area = ar, .before = 1),
    error = function(e) NULL)

  satt <- tryCatch(
    as.data.frame(anova_satt(fit, type = 2)) %>%
      tibble::rownames_to_column("term") %>%
      mutate(analyte = analyte_label[akey], area = ar, .before = 1),
    error = function(e) NULL)

  emm_list <- emm_tukey_cld(fit, ~ treatment)
  if (!is.null(emm_list))
    emm_list <- purrr::map(emm_list, function(t)
      if (is.data.frame(t)) t %>% mutate(analyte = analyte_label[akey], area = ar, .before = 1)
      else t)  # pass-through for non-data-frame slots (e.g. $simple list)

  # DHARMa — skip for lm (DHARMa works but is less informative)
  dharma_path_cell <- file.path(FIG_DIR,
    sprintf("diag_cell_%s_%s_dharma.png",
            gsub("[^A-Za-z0-9]", "_", analyte_label[akey]), ar))
  if (inherits(fit, c("lmerMod","lmerModLmerTest")))
    dharma_panel(fit, dharma_path_cell,
                 title = sprintf("DHARMa: %s × %s [%s]", analyte_label[akey], ar, re_label))

  # Extract treatment effect quantities for cross-summary
  kr_p <- NA_real_; kr_F <- NA_real_; kr_df1 <- NA_real_; kr_df2 <- NA_real_
  if (!is.null(kr)) {
    idx <- grep("^treatment$", kr$term, ignore.case = TRUE)[1]
    if (!is.na(idx)) {
      fv     <- kr[["F"]][idx]
      if (is.null(fv) || !is.finite(fv)) fv <- kr[["F value"]][idx]
      kr_F   <- fv
      kr_p   <- kr[["Pr(>F)"]][idx]
      kr_df1 <- if ("Df"    %in% names(kr)) kr[["Df"]][idx]    else kr[["NumDF"]][idx]
      kr_df2 <- if ("Df.res" %in% names(kr)) kr[["Df.res"]][idx]
               else if ("DenDF"  %in% names(kr)) kr[["DenDF"]][idx]
               else if ("Den Df" %in% names(kr)) kr[["Den Df"]][idx]
               else {
                 res_row <- kr[grepl("^Residual", kr$term, ignore.case = TRUE), ]
                 if (nrow(res_row) > 0) res_row[["Df"]][1] else NA_real_
               }
    }
  }
  sm        <- summary(fit)$coefficients
  treat_idx <- grep("^treatment", rownames(sm))[1]
  beta      <- if (!is.na(treat_idx)) sm[treat_idx, "Estimate"]   else NA_real_
  se        <- if (!is.na(treat_idx)) sm[treat_idx, "Std. Error"] else NA_real_

  # 95% CI on treatment contrast (log scale), t-based on the model's own KR
  # denominator df rather than a z/1.96 Wald interval -- some cells have
  # df as low as ~10-15 (few trials contribute), where the normal approximation
  # is visibly too narrow relative to the t distribution.
  .ci_df <- if (is.finite(kr_df2)) kr_df2 else if (inherits(fit, "lm")) fit$df.residual else Inf
  .tcrit <- suppressWarnings(qt(0.975, .ci_df))
  if (!is.finite(.tcrit)) .tcrit <- 1.96
  ci_lo <- beta - .tcrit * se
  ci_hi <- beta + .tcrit * se

  # Effect size: signed Hedges' g from the model's own treatment contrast
  # (noncentral-t CI via effectsize::t_to_d), NOT effectsize::cohens_d() on the
  # raw two-group data -- the latter treats every fish as an independent
  # replicate and ignores the (1|trial) clustering that Part 2c of the
  # random-effects policy forces into this exact model, understating the
  # contrast's true uncertainty (Nakagawa & Cuthill 2007 §III.5).
  eff_d <- tryCatch({
    ctr1 <- emm_list$pairs
    if (is.null(ctr1) || !nrow(ctr1)) NULL else {
      df_ctr <- ctr1$df[1]
      t_ctr  <- if ("t.ratio" %in% names(ctr1)) ctr1$t.ratio[1] else ctr1$estimate[1] / ctr1$SE[1]
      if (!is.finite(df_ctr) || !is.finite(t_ctr) || df_ctr <= 1) NULL else {
        J     <- 1 - 3 / (4 * df_ctr - 1)
        d_res <- effectsize::t_to_d(t_ctr, df_ctr, ci = 0.95)
        data.frame(Cohens_d = d_res$d[1] * J, CI_low = d_res$CI_low[1] * J,
                   CI_high = d_res$CI_high[1] * J, method = "hedges_g_t_to_d")
      }
    }
  }, error = function(e) NULL)

  # Descriptive by treatment (raw scale)
  descr <- d %>%
    group_by(treatment) %>%
    summarise(n      = n(),
              mean   = mean(value,  na.rm = TRUE),
              sd     = sd(value,    na.rm = TRUE),
              sem    = sd / sqrt(n),
              median = median(value),
              q1     = quantile(value, 0.25),
              q3     = quantile(value, 0.75),
              .groups = "drop") %>%
    mutate(analyte = analyte_label[akey], area = ar, .before = 1)

  list(
    analyte_key       = akey,
    analyte           = analyte_label[akey],
    area              = ar,
    cell_id           = paste(akey, ar, sep = "__"),
    re_label          = re_label,
    re_formula        = re_formula,
    model_formula     = model_formula,
    trans_type        = trans$type,
    trans_label       = trans$label,
    trans_lambda      = trans$lambda,
    trans_rationale   = trans$rationale,
    fit               = fit,
    n                 = n_total,
    norm_lev      = {
      base <- merge(norm_row, lev_raw_row, by = "outcome", all = TRUE)
      names(base)[names(base) %in% c("levene_F","df1","df2","p")] <-
        c("levene_F_raw","df1_raw","df2_raw","levene_p_raw")
      log_sub <- lev_log_row[, c("outcome","levene_F","df1","df2","p"), drop = FALSE]
      names(log_sub) <- c("outcome","levene_F_log","df1_log","df2_log","levene_p_log")
      merge(base, log_sub, by = "outcome", all = TRUE)
    },
    re_aicc_table = re_aicc_table,
    r2            = r2,
    anova_kr      = kr,
    anova_satt    = satt,
    emm           = emm_list,
    eff_d         = eff_d,
    dharma_path   = dharma_path_cell,
    descr         = descr,
    treat_beta    = beta,
    treat_se      = se,
    treat_ci_lo   = ci_lo,
    treat_ci_hi   = ci_hi,
    treat_F       = kr_F,
    treat_df1     = kr_df1,
    treat_df2     = kr_df2,
    treat_p       = kr_p
  )
}

# Run 7 analytes × 4 areas = up to 28 per-cell models
mono_cell_combos <- expand.grid(akey = analyte_cols, area = expected_areas,
                                stringsAsFactors = FALSE)
mono_cell_models <- purrr::compact(purrr::map2(
  mono_cell_combos$akey, mono_cell_combos$area,
  function(ak, ar) tryCatch(
    fit_one_mono_area(ak, ar),
    error = function(e) {
      message(sprintf("fit_one_mono_area failed [%s × %s]: %s",
                      analyte_label[ak], ar, conditionMessage(e)))
      NULL
    })
))
names(mono_cell_models) <- vapply(mono_cell_models, function(x) x$cell_id, character(1))
message(sprintf("Fitted %d monoamine × area cell models (primary).", length(mono_cell_models)))

# Export per-cell results
cell_bind <- function(field) purrr::map_dfr(mono_cell_models, function(x) x[[field]])

write.csv(cell_bind("norm_lev"),      file.path(DATA_DIR, "cell_normality_levene.csv"),  row.names = FALSE)
write.csv(cell_bind("re_aicc_table"), file.path(DATA_DIR, "cell_re_selection_aicc.csv"), row.names = FALSE)
write.csv(cell_bind("r2"),            file.path(DATA_DIR, "cell_r2.csv"),                row.names = FALSE)
write.csv(cell_bind("anova_kr"),      file.path(DATA_DIR, "cell_anova_kr.csv"),          row.names = FALSE)
write.csv(cell_bind("anova_satt"),    file.path(DATA_DIR, "cell_anova_satt.csv"),        row.names = FALSE)
write.csv(cell_bind("descr"),         file.path(DATA_DIR, "cell_descriptive.csv"),       row.names = FALSE)

cell_emm   <- purrr::map_dfr(mono_cell_models, function(x) if (!is.null(x$emm)) x$emm$emm)
cell_pairs <- purrr::map_dfr(mono_cell_models, function(x) if (!is.null(x$emm)) x$emm$pairs)
cell_cld   <- purrr::map_dfr(mono_cell_models, function(x) if (!is.null(x$emm)) x$emm$cld)
if (nrow(cell_emm)   > 0) write.csv(cell_emm,   file.path(DATA_DIR, "cell_emm.csv"),   row.names = FALSE)
if (nrow(cell_pairs) > 0) write.csv(cell_pairs, file.path(DATA_DIR, "cell_pairs.csv"), row.names = FALSE)
if (nrow(cell_cld)   > 0) write.csv(cell_cld,   file.path(DATA_DIR, "cell_cld.csv"),   row.names = FALSE)

# BH-FDR within each declared neurotransmitter-system family, FLAT across all four
# brain areas (M1 / R1-8). Replaces the earlier per-area scheme (4 families of 7):
# the family definition, not the correction method, determines which effects survive
# here, so it is pre-declared above and reported in Methods.
cell_cross_summary <- local({
  df <- purrr::map_dfr(mono_cell_models, function(x) data.frame(
    area          = x$area,
    analyte       = x$analyte,
    n             = x$n,
    re_spec       = x$re_label,
    model         = x$model_formula,
    beta          = x$treat_beta,
    SE            = x$treat_se,
    CI_lo_95      = x$treat_ci_lo,
    CI_hi_95      = x$treat_ci_hi,
    F_stat        = x$treat_F,
    df1           = x$treat_df1,
    df2           = x$treat_df2,
    p_raw         = x$treat_p,
    R2m           = if (!is.null(x$r2$R2_marginal))    x$r2$R2_marginal[1]    else NA_real_,
    R2c           = if (!is.null(x$r2$R2_conditional)) x$r2$R2_conditional[1] else NA_real_,
    ICC_adjusted  = if (!is.null(x$r2$ICC_adjusted))   x$r2$ICC_adjusted[1]   else NA_real_
  ))
  # No correction: significance is read directly from the raw p (alpha = 0.05).
  df$sig_raw <- vapply(df$p_raw, fmt_p, character(1))
  df[order(df$area, df$analyte), ]
})
message(sprintf("Cell-level treatment effects (uncorrected, raw p): %d of %d cells significant at p<0.05",
                sum(cell_cross_summary$p_raw < 0.05), nrow(cell_cross_summary)))
write.csv(cell_cross_summary, file.path(DATA_DIR, "cell_cross_summary_bh.csv"), row.names = FALSE)
message(sprintf("Per-cell cross-summary written: %d rows.", nrow(cell_cross_summary)))

# Per-analyte × area caption lookup from per-cell models (used by legacy plots).
# No correction applied: captions report the raw p only, and figure asterisks
# (`is_sig_q`, name retained for downstream compatibility) are driven purely by
# raw p < 0.05.
cell_anova_caps <- cell_cross_summary %>%
  mutate(
    analyte_key = names(analyte_label)[match(analyte, analyte_label)],
    # Partial eta-squared from the term's own F and df: eta^2_p = F*df1/(F*df1+df2).
    eta2_p = ifelse(is.finite(F_stat) & is.finite(df1) & is.finite(df2),
                    (F_stat * df1) / (F_stat * df1 + df2), NA_real_),
    # Noncentral-F 95% CI on eta2_p (Nakagawa & Cuthill 2007 SII.4), added
    # 2026-08-08 alongside the same fix in the STEP5 behavioural engine.
    eta2_p_ci = purrr::pmap(list(F_stat, df1, df2), function(f, d1, d2) {
      if (!is.finite(f) || !is.finite(d1) || !is.finite(d2))
        return(c(NA_real_, NA_real_))
      es <- tryCatch(effectsize::F_to_eta2(f, d1, d2, ci = 0.95, alternative = "two.sided"),
                     error = function(e) NULL)
      if (is.null(es) || !nrow(es)) return(c(NA_real_, NA_real_))
      c(max(es$CI_low[1], 0), min(es$CI_high[1], 1))
    }),
    eta2_p_lo = vapply(eta2_p_ci, `[`, numeric(1), 1),
    eta2_p_hi = vapply(eta2_p_ci, `[`, numeric(1), 2),
    # Single-line caption: F, p and the effect size point estimate together,
    # comma-separated, matching the project-wide convention (2026-08-09).
    # The CI bracket is dropped here (kept in the underlying CSV/Results
    # text) specifically because it was what pushed this past one line at
    # the Figure 11/12 caption size (34, 3740 px canvas); the point estimate
    # alone fits comfortably.
    # F at 4 significant figures (fmt_F) and df2 shown exactly as the
    # manuscript text does: integer when whole, one decimal otherwise
    # (26.8, 27.2, 27.9) -- the project's single number-display convention.
    df2_lab = fmt_F(df2),
    caption = dplyr::case_when(
      !is.na(F_stat) & !is.na(df2) ~ sprintf(
        "F<sub>%d,%s</sub> = %s, p = %s, &eta;<sup>2</sup><sub>p</sub> = %s",
        round(df1), suppressWarnings(round(as.numeric(df2_lab))), fmt_Fstat(F_stat),
        ifelse(p_raw < 0.001, "< 0.001", fmt3(p_raw)), fmt_es3(eta2_p)),
      !is.na(F_stat) ~ sprintf(
        "F<sub>%d</sub> = %s, p = %s",
        round(df1), fmt_Fstat(F_stat),
        ifelse(p_raw < 0.001, "< 0.001", fmt3(p_raw))),
      TRUE ~ ""
    ),
    is_sig_q = !is.na(p_raw) & p_raw < 0.05
  ) %>%
  dplyr::select(analyte_key, area, analyte, p_val = p_raw, eta2_p, eta2_p_lo, eta2_p_hi,
                is_sig_q, caption)

# ---- §6 CROSS-ANALYSIS SUMMARY (region-aggregated per monoamine, raw p) ---
# No correction applied (see decision log). File name kept as
# cross_summary_bh.csv for downstream compatibility.
cross_summary_rows <- purrr::map_dfr(mono_models, function(x) {
  data.frame(
    family   = "monoamine_aggregated",
    outcome  = x$label, n = x$n,
    re_spec  = x$re_label,
    beta     = x$treat_beta, SE = x$treat_se,
    F        = x$treat_F, df1 = x$treat_df1, df2 = x$treat_df2,
    p_raw    = x$treat_p
  )
})
cross_summary_rows$sig_raw <- vapply(cross_summary_rows$p_raw, fmt_p, character(1))

cort_kr_p <- if (!is.null(cort_anova_kr)) {
  idx <- grep("^condition$", cort_anova_kr$term)[1]
  if (!is.na(idx)) cort_anova_kr[["Pr(>F)"]][idx] else NA_real_
} else NA_real_
cort_kr_F <- if (!is.null(cort_anova_kr)) {
  idx <- grep("^condition$", cort_anova_kr$term)[1]
  if (!is.na(idx)) {
    fv <- cort_anova_kr[["F"]][idx]; if (is.null(fv) || !is.finite(fv)) cort_anova_kr[["F value"]][idx] else fv
  } else NA_real_
} else NA_real_
# df1/df2 (previously hardcoded NA here, which forced the cross-summary forest
# plot to fall back to a 1.96 z-approximation for the cortisol row specifically
# even though the KR fit has a perfectly good denominator df). Extracted the
# same way as the per-cell monoamine models above.
cort_kr_df1 <- NA_real_; cort_kr_df2 <- NA_real_
if (!is.null(cort_anova_kr)) {
  idx <- grep("^condition$", cort_anova_kr$term)[1]
  if (!is.na(idx)) {
    cort_kr_df1 <- if ("Df" %in% names(cort_anova_kr)) cort_anova_kr[["Df"]][idx]
                   else cort_anova_kr[["NumDF"]][idx]
    cort_kr_df2 <- if ("Df.res" %in% names(cort_anova_kr)) cort_anova_kr[["Df.res"]][idx]
                   else if ("DenDF"  %in% names(cort_anova_kr)) cort_anova_kr[["DenDF"]][idx]
                   else if ("Den Df" %in% names(cort_anova_kr)) cort_anova_kr[["Den Df"]][idx]
                   else NA_real_
  }
}
cortisol_row <- data.frame(
  family = "cortisol", outcome = "plasma_Cortisol",
  n = nrow(cort_dat), re_spec = cort_re_label,
  beta = cort_treat_beta, SE = cort_treat_se,
  F = cort_kr_F, df1 = cort_kr_df1, df2 = cort_kr_df2,
  p_raw = cort_kr_p, sig_raw = fmt_p(cort_kr_p)
)
cross_summary <- rbind(cortisol_row, cross_summary_rows)
write.csv(cross_summary, file.path(DATA_DIR, "cross_summary_bh.csv"), row.names = FALSE)

# Stratified treatment-within-area simple effects (7 x 4 = 28), raw p only.
stratified_summary <- if (!is.null(mono_pairs_str) && nrow(mono_pairs_str) > 0) {
  mono_pairs_str %>%
    select(any_of(c("outcome","contrast","area","estimate","SE","df","t.ratio","p.value"))) %>%
    mutate(sig_raw = vapply(p.value, fmt_p, character(1)))
} else NULL
if (!is.null(stratified_summary))
  write.csv(stratified_summary, file.path(DATA_DIR, "cross_summary_bh_stratified.csv"),
            row.names = FALSE)

# ---- §7 PSEUDOREPLICATION ASSESSMENT (cortisol) --------------------------
pseudorep_table <- if (cort_re_label != "lm (no RE)") {
  pseudorep_compare(fixed_cort_pri, cort_re_formula, re_cort_cand,
                    data = cort_dat, term_name = "condition1")
} else NULL
if (is.null(pseudorep_table)) {
  # Try with fallback term name for non-sum-coded
  pseudorep_table <- if (cort_re_label != "lm (no RE)") {
    pseudorep_compare(fixed_cort_pri, cort_re_formula, re_cort_cand,
                      data = cort_dat, term_name = "conditiontreat")
  } else NULL
}
if (!is.null(pseudorep_table))
  write.csv(pseudorep_table, file.path(DATA_DIR, "pseudorep_check.csv"), row.names = FALSE)

# ---- §7b PSEUDOREPLICATION ASSESSMENT (monoamine per-cell models, M2) -----
# Now that (1|trial) competes with (1|date)/(1|tank) as an RE candidate for the
# primary per-cell models, run the same probe used for cortisol: refit each
# cell's winning model with every non-winning RE and flag beta/SE instability.
get_cell_data <- function(akey, ar) {
  mono_long %>%
    filter(analyte_key == akey, area == ar,
           is.finite(value), value > 0, !is.na(treatment)) %>%
    droplevels()
}
mono_cell_re_cand <- list("(1|date)" = "(1|date)", "(1|tank)" = "(1|tank)",
                          "(1|trial)" = "(1|trial)")
mono_pseudorep_table <- purrr::map_dfr(mono_cell_models, function(x) {
  if (grepl("^lm ", x$re_label)) return(NULL)  # no RE was fit; nothing to probe
  d <- get_cell_data(x$analyte_key, x$area)
  d$log_v <- apply_transform(d$value, list(type = x$trans_type, lambda = x$trans_lambda))
  d <- set_sum_contrasts(d, "treatment")
  # x$re_formula (from fit_one_mono_area's own select_re_aicc call), NOT a
  # re-derivation via x$re_label indexed into this function's local
  # mono_cell_re_cand -- see the cortisol fix above for why that silently
  # returns NULL whenever force_re mints a combo label this list doesn't
  # happen to spell the same way.
  win_term <- x$re_formula
  if (is.null(win_term)) return(NULL)
  out <- tryCatch(
    pseudorep_compare(log_v ~ treatment, win_term, mono_cell_re_cand,
                      data = d, term_name = "treatment1"),
    error = function(e) NULL)
  if (is.null(out)) return(NULL)
  out %>% mutate(analyte = x$analyte, area = x$area, .before = 1)
})
if (nrow(mono_pseudorep_table) > 0)
  write.csv(mono_pseudorep_table, file.path(DATA_DIR, "pseudorep_check_monoamines.csv"),
            row.names = FALSE)
message(sprintf("Monoamine pseudoreplication probe: %d cell x RE rows; %d flagged (|dSE| > 30%% or sign flip).",
                nrow(mono_pseudorep_table), sum(mono_pseudorep_table$flag, na.rm = TRUE)))

# ---- §4b ORDER / DAY-OF-TESTING EFFECT (M4/R2-5) --------------------------
# Trial order/day is tested as order_idx (chronological rank of testing date,
# 1..8) + treatment:order_idx, added to each model's existing winning formula.
#
# DESIGN LIMITATION (documented per M4's "test AND document" requirement):
# In this design, order_idx is NOT identifiable independently of tank or of
# fish_density:
#   - Between tanks: each tank occupies its own exclusive, non-overlapping
#     pair of testing days, so the GLOBAL day sequence (1..8) is a perfect
#     function of tank identity — order and tank cannot both be estimated
#     (aliased) when tank is also a random effect.
#   - Within a tank: the two testing days always run fish_density = {16,12}
#     on day 1 and {8,4} on day 2, and the first trial run each day always
#     has the higher of that day's two densities. So order_idx and
#     fish_density are related by an exact affine transform in this data
#     (fish_density = 20 - 4*order_idx_within_tank) — an order test and a
#     density test are mathematically the same test here and cannot be
#     told apart. Order effects reported below should be read as
#     "order-or-density", not as evidence of chronological drift per se.
cort_order_dat <- cort_dat
fit_cort_order <- tryCatch({
  if (cort_re_label != "lm (no RE)") {
    f <- as.formula(paste("log_cort ~ condition * order_idx +", cort_re_formula))
    suppressMessages(suppressWarnings(lmerTest::lmer(f, data = cort_order_dat, REML = TRUE)))
  } else {
    lm(log_cort ~ condition * order_idx, data = cort_order_dat)
  }
}, error = function(e) NULL)

extract_order_rows <- function(fit, outcome, model_label) {
  if (is.null(fit)) return(NULL)
  a <- tryCatch(as.data.frame(anova_kr(fit, type = 3)) %>%
                  tibble::rownames_to_column("term"),
                error = function(e) NULL)
  if (is.null(a)) return(NULL)
  a %>% filter(grepl("order_idx", term)) %>%
    mutate(outcome = outcome, model = model_label, .before = 1)
}

order_effect_rows <- list(
  extract_order_rows(fit_cort_order, "plasma_cortisol",
                     paste0("log_cort ~ condition*order_idx + ", cort_re_label))
)

mono_order_rows <- purrr::map(mono_cell_models, function(x) {
  d <- get_cell_data(x$analyte_key, x$area)
  d$log_v <- apply_transform(d$value, list(type = x$trans_type, lambda = x$trans_lambda))
  d <- set_sum_contrasts(d, "treatment")
  re_term <- x$re_formula  # see mono_pseudorep_table above for why not x$re_label-indexed
  fit_o <- tryCatch({
    if (!is.null(re_term)) {
      f <- as.formula(paste("log_v ~ treatment * order_idx +", re_term))
      suppressMessages(suppressWarnings(lmerTest::lmer(f, data = d, REML = TRUE)))
    } else {
      lm(log_v ~ treatment * order_idx, data = d)
    }
  }, error = function(e) NULL)
  extract_order_rows(fit_o, sprintf("%s [%s]", x$analyte, x$area),
                     paste0(x$model_formula, " + order_idx (+ interaction)"))
})

order_effect_summary <- dplyr::bind_rows(c(order_effect_rows, mono_order_rows))
write.csv(order_effect_summary, file.path(DATA_DIR, "order_effect_summary.csv"), row.names = FALSE)
message(sprintf("Order/day effect (M4): %d term-rows across %d models.",
                nrow(order_effect_summary),
                length(order_effect_rows) + length(mono_order_rows)))

# ---- §9 SECONDARY EXPLORATORY: SEX EFFECTS (NON-CONFIRMATORY) ------------
# 9.1 Cortisol: condition * sex
cort_sex_dat <- cort_full %>% filter(is.finite(log_cort), sex %in% c("F","M")) %>%
  droplevels() %>% set_sum_contrasts(c("condition","sex"))

fit_cort_sex <- tryCatch({
  if (cort_re_label != "lm (no RE)") {
    f <- as.formula(paste("log_cort ~ condition * sex +", cort_re_formula))
    suppressMessages(suppressWarnings(lmerTest::lmer(f, data = cort_sex_dat, REML = TRUE)))
  } else {
    lm(log_cort ~ condition * sex, data = cort_sex_dat)
  }
}, error = function(e) NULL)

sex_secondary_cortisol <- if (!is.null(fit_cort_sex)) {
  a <- as.data.frame(anova_kr(fit_cort_sex, type = 3)) %>%
    tibble::rownames_to_column("term") %>%
    mutate(outcome = "plasma_cortisol", model = "log_cort ~ condition*sex (secondary)",
           .before = 1)
  a
} else NULL
if (!is.null(sex_secondary_cortisol))
  write.csv(sex_secondary_cortisol,
            file.path(DATA_DIR, "sex_secondary_cortisol.csv"), row.names = FALSE)

cort_sex_emm <- if (!is.null(fit_cort_sex)) {
  tryCatch(as.data.frame(emmeans(fit_cort_sex, ~ condition * sex)),
           error = function(e) NULL)
} else NULL
if (!is.null(cort_sex_emm))
  write.csv(cort_sex_emm,
            file.path(DATA_DIR, "sex_secondary_cortisol_emmeans.csv"), row.names = FALSE)

# 9.2 Monoamines: per-analyte sex as additive covariate
fit_one_mono_sex <- function(akey) {
  d <- mono_long %>% filter(analyte_key == akey, is.finite(value), value > 0,
                            sex %in% c("F","M")) %>% droplevels()
  if (nrow(d) < 10) return(NULL)
  d$log_v <- log(d$value)
  d <- set_sum_contrasts(d, c("treatment","area","sex"))
  m <- tryCatch(suppressMessages(suppressWarnings(
    lmerTest::lmer(log_v ~ treatment * area + sex + (1|sample_id), data = d, REML = TRUE)
  )), error = function(e) NULL)
  if (is.null(m)) return(NULL)
  a <- as.data.frame(anova_kr(m, type = 3)) %>%
    tibble::rownames_to_column("term") %>%
    mutate(outcome = analyte_label[akey],
           model = "log_v ~ treatment*area + sex + (1|sample_id) (secondary)",
           .before = 1)
  a
}
sex_secondary_monoamines <- purrr::map_dfr(analyte_cols, fit_one_mono_sex)
if (nrow(sex_secondary_monoamines) > 0)
  write.csv(sex_secondary_monoamines,
            file.path(DATA_DIR, "sex_secondary_monoamines.csv"), row.names = FALSE)

# 9.3 Monoamines: per-cell (analyte × area) sex effect
fit_one_mono_sex_cell <- function(akey, ar) {
  d <- mono_long %>%
    filter(analyte_key == akey, area == ar,
           is.finite(value), value > 0,
           !is.na(treatment), sex %in% c("F", "M")) %>%
    droplevels()
  if (nrow(d) < 8 || length(unique(d$sex)) < 2) return(NULL)
  d$log_v <- log(d$value)
  d <- set_sum_contrasts(d, c("treatment", "sex"))

  fixed_form <- log_v ~ treatment + sex
  re_cand    <- list("(1|date)" = "(1|date)", "(1|tank)" = "(1|tank)")
  sel        <- select_re_aicc(fixed_form, re_cand, data = d, force_re = "(1|trial)")

  if (!is.null(sel$winner)) {
    fit           <- sel$winner
    re_label      <- sel$label
    model_formula <- paste("log_v ~ treatment + sex +", re_label)
  } else {
    fit           <- tryCatch(lm(log_v ~ treatment + sex, data = d),
                              error = function(e) NULL)
    re_label      <- "lm (no RE)"
    model_formula <- "log_v ~ treatment + sex [lm, no RE]"
  }
  if (is.null(fit)) return(NULL)

  # Extract sex term p-value
  if (inherits(fit, c("lmerMod", "lmerModLmerTest"))) {
    kr  <- tryCatch(as.data.frame(anova_kr(fit, type = 3)), error = function(e) NULL)
    if (is.null(kr)) return(NULL)
    sex_row <- kr[rownames(kr) == "sex", , drop = FALSE]
    if (nrow(sex_row) == 0) return(NULL)
    f_val  <- sex_row[["F value"]]
    df1    <- sex_row[["Num Df"]]
    df2    <- sex_row[["Den Df"]]
    p_val  <- sex_row[["Pr(>F)"]]
    method <- "KR"
  } else {
    s <- summary(fit)$coefficients
    sex_rows <- grep("^sex", rownames(s))
    if (length(sex_rows) == 0) return(NULL)
    aov_sex <- anova(fit)["sex", , drop = FALSE]
    f_val  <- aov_sex[["F value"]]
    df1    <- aov_sex[["Df"]]
    df2    <- aov_sex[["Df"]][length(aov_sex[["Df"]])]   # residual df from lm anova
    df2    <- df.residual(fit)
    p_val  <- aov_sex[["Pr(>F)"]]
    method <- "lm-F"
  }

  data.frame(
    area          = ar,
    analyte       = analyte_label[akey],
    n             = nrow(d),
    re_spec       = re_label,
    model         = model_formula,
    F_stat        = f_val,
    df1           = df1,
    df2           = df2,
    p_sex_raw     = p_val,
    method        = method,
    stringsAsFactors = FALSE
  )
}

sex_cell_combos <- expand.grid(akey = analyte_cols, area = expected_areas,
                               stringsAsFactors = FALSE)
sex_cell_results <- purrr::compact(purrr::map2(
  sex_cell_combos$akey, sex_cell_combos$area,
  function(ak, ar) tryCatch(
    fit_one_mono_sex_cell(ak, ar),
    error = function(e) {
      message(sprintf("fit_one_mono_sex_cell failed [%s x %s]: %s",
                      analyte_label[ak], ar, conditionMessage(e)))
      NULL
    })
))
sex_cell_df <- dplyr::bind_rows(sex_cell_results)

# No correction applied (see decision log); raw p only.
if (nrow(sex_cell_df) > 0) {
  sex_cell_df <- sex_cell_df %>%
    arrange(area, analyte)
  write.csv(sex_cell_df,
            file.path(DATA_DIR, "sex_cell_monoamines.csv"), row.names = FALSE)
  message(sprintf("Per-cell sex analysis: %d cells. Min raw p = %.4f (%s x %s)",
                  nrow(sex_cell_df),
                  min(sex_cell_df$p_sex_raw, na.rm = TRUE),
                  sex_cell_df$analyte[which.min(sex_cell_df$p_sex_raw)],
                  sex_cell_df$area[which.min(sex_cell_df$p_sex_raw)]))
}

# ---- §1 / §2 STUDY DESIGN + DECISION LOG (CSV exports) -------------------
re_reference <- data.frame(
  re   = c("tank","trial","plate","plate_date","sample_id","area"),
  role = c("housing tank",
           "arena trial = tank x removal-position/density (unit of replication for treatment; M2/R1-6; corrected 2026-08-08 from a raw trial column that mixed both treatment arms per group)",
           "cortisol assay plate","plate run date",
           "fish identity (within-fish for monoamines)",
           "brain region — fixed factor in monoamine model"),
  used_in = c("cortisol + monoamine RE candidate","cortisol + monoamine RE candidate (M2)",
              "cortisol RE candidate","diagnostic only",
              "monoamine within-fish RE","fixed factor (monoamines)")
)
write.csv(re_reference, file.path(DATA_DIR, "random_effects_reference.csv"), row.names = FALSE)

# Statistical Decisions Log (M3/R1-7): one row per analytical choice, each
# marked a priori (fixed before seeing treatment-effect results) vs data-driven
# (selected via a fitting criterion, e.g. AICc/Shapiro-Wilk), with rationale.
decision_log <- data.frame(
  step = paste0("2.", 1:15),
  decision = c(
    "Sum-to-zero contrasts on factors (Type III ANOVA interpretability).",
    "log()/Box-Cox transform on cortisol + monoamines (positive, right-skewed).",
    "Random-effect selection by AICc; winner refit with REML=TRUE.",
    "Type III KR F-test as primary inference (anti-conservative Wald χ² rejected).",
    "Satterthwaite ANOVA reported as sensitivity row.",
    "Tukey-adjusted pairwise EMM contrasts; multcomp::cld for letter display.",
    "No adjustment for multiple comparisons is applied across the 24 analyte x brain-region cells; raw (unadjusted) p-values are reported and interpreted at alpha = 0.05 per test (author decision, 2026-08-07). Family-based FDR/FWER correction (Benjamini-Hochberg/Yekutieli, Bonferroni, Sidak, Holm, Hochberg, Hommel; several candidate family structures) was evaluated but is not part of the reported analysis; that code is retained separately in 00_shared/correction_methods.R.",
    "Brain area = fixed within-fish factor; (1|sample_id) absorbs intra-fish correlation.",
    "Pseudoreplication probe: refit with non-winner REs; flag β/SE shifts > 30%.",
    "NE (noradrenaline) excluded from all analyses, tables and figures (author decision, 2026-08-04): 6 analytes x 4 areas = 24 cells.",
    "M2: (1|trial) added to monoamine RE candidate sets (stratified + per-cell) so fish sampled from the same video are not treated as independent replicates.",
    "M4: trial order/day added as a fixed covariate (order_idx) + treatment:order_idx interaction to cortisol + monoamine models, where identifiable from tank/density.",
    "M4 limitation: order_idx is a monotone (in fact affine) function of fish_density within each tank, and the between-tank order sequence is fully aliased with tank identity — order cannot be separated from density or tank in this design.",
    "M5: trial-level correlation between behavioural flow-preference and trial-aggregated monoamines computed where both are available for the same trial (see endocrine_behaviour_association.R); reported as association only, not causal/directional.",
    "Fixed random seed (20260706) and sessionInfo() recorded for reproducibility."
  ),
  basis = c(
    "a priori","data-driven (Box-Cox grid + Shapiro-Wilk)","data-driven (AICc)",
    "a priori","a priori (sensitivity)","a priori",
    "a priori (correction procedure evaluated then declined)",
    "a priori","a priori (probe design)",
    "a priori (analyte scope)",
    "a priori","a priori (test) / data-driven (identifiability depends on design)",
    "a priori (documented design constraint)","a priori","a priori"
  )
)
write.csv(decision_log, file.path(DATA_DIR, "decision_log.csv"), row.names = FALSE)

# Full 24-cell table for the Supplementary Materials, raw p only (no
# correction). The multi-procedure sensitivity comparison that used to be
# generated here has been moved to 00_shared/correction_methods.R
# (endocrine_correction_sensitivity()) and is no longer run by default.
local({
  full <- cell_cross_summary[order(match(cell_cross_summary$area, expected_areas),
                                   cell_cross_summary$analyte),
                             c("area","analyte","n","re_spec","F_stat","df1","df2","p_raw")]
  write.csv(full, file.path(DATA_DIR, "supplementary_full_24_cell_table.csv"),
            row.names = FALSE)
  message(sprintf("Full 24-cell table written (uncorrected): %d rows.", nrow(full)))
})

# ---- §8 DESCRIPTIVES — keep existing wide CV plots ------------------------
# Standard curve overlay
sc_df <- assay_long %>% filter(is_standard) %>%
  mutate(nominal = NOMINAL_STD[sample_raw]) %>% filter(!is.na(nominal) & !is.na(raw_OD))
p_curve <- ggplot(sc_df, aes(nominal, raw_OD, colour = plate)) +
  geom_point(alpha = 0.7) + stat_summary(fun = mean, geom = "line") +
  scale_x_log10() +
  labs(title = "Cortisol standard curve overlay across B3 plates",
       x = "Nominal pg/mL (log scale)", y = "Raw OD") + theme_bw()
ggsave(file.path(FIG_DIR, "fig_standard_curves.png"), p_curve, width = 7, height = 4.5, dpi = 200)
ggsave(file.path(FIG_DIR, "fig_standard_curves.pdf"), p_curve, width = 7, height = 4.5)

p_cv <- ggplot(intra_assay_cv %>% filter(is.finite(cv_intra_pct)),
               aes(plate, cv_intra_pct, fill = plate)) +
  geom_boxplot(alpha = 0.7) + geom_jitter(width = 0.15, alpha = 0.5) +
  geom_hline(yintercept = c(15, 20), linetype = "dashed", colour = "grey40") +
  labs(title = "Intra-assay CV by plate (B3 unknowns)",
       y = "CV (%)", x = "Plate") + theme_bw() + theme(legend.position = "none")
ggsave(file.path(FIG_DIR, "fig_intra_assay_cv.png"), p_cv, width = 6.5, height = 4.5, dpi = 200)
ggsave(file.path(FIG_DIR, "fig_intra_assay_cv.pdf"), p_cv, width = 6.5, height = 4.5)

# Cross-summary forest plot (8 outcomes, uncorrected raw p throughout).
# CI uses each row's own KR denominator df via qt() rather than a blanket
# z = 1.96 -- some cells have df as low as ~10-15, where the normal
# approximation is visibly too narrow relative to the t distribution.
plot_df <- cross_summary %>%
  mutate(.tcrit = ifelse(is.finite(df2), qt(0.975, df2), 1.96),
         label = paste0(outcome, " (", family, ")"),
         lower = beta - .tcrit * SE,
         upper = beta + .tcrit * SE) %>%
  dplyr::select(-.tcrit)
p_forest <- ggplot(plot_df, aes(x = beta, y = reorder(label, beta), color = family)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_pointrange(aes(xmin = lower, xmax = upper), size = 0.4) +
  labs(x = "treatment β (log scale) ± 95% CI", y = NULL,
       title = "Treatment effect across 8 outcomes (raw p, no correction applied)") +
  theme_bw() + theme(legend.position = "bottom")
ggsave(file.path(FIG_DIR, "emm_forest_treatment.png"), p_forest,
       width = 7.5, height = 5, dpi = 200)

# ---- 12. Cortisol × monoamine correlations (exploratory) -------------------
cort_mono_corr <- mono_long %>%
  filter(is.finite(value), is.finite(plasma_cortisol)) %>%
  group_by(area, analyte, sex) %>%
  summarise(n = n(),
            spearman_rho = suppressWarnings(cor(plasma_cortisol, value, method = "spearman")),
            .groups = "drop")
write.csv(cort_mono_corr, file.path(OUT, "cortisol_monoamine_spearman.csv"), row.names = FALSE)

# ---- 13. Reconciliation tests ----------------------------------------------
recon <- b3_cort %>%
  inner_join(assay_wide %>% select(sample_id, plate, mean_conc),
             by = c("sample_id","plate")) %>%
  mutate(derived_ngmL = mean_conc / 1000,
         abs_diff     = abs(plasma_cortisol - derived_ngmL),
         rel_diff_pct = abs_diff / plasma_cortisol * 100)
write.csv(recon, file.path(OUT, "concentration_reconciliation.csv"), row.names = FALSE)

# ---- 14b. Legacy-style publication plots ------------------------------------
# Ported from:
#   cort stats graph choice nov 2025.R    — jitter + SEM + mean line + N labels
#   monoamine B3 tests and graphs NOV 2025.R — per area × analyte jitter + SEM
# Style: Okabe-Ito palette, no grid, axis lines, conditional significance star
# Output: B3_Final_Report/
# -----------------------------------------------------------------------

FINAL_DIR <- file.path(OUT, "B3_Final_Report")   # honour the revised OUT tree
dir.create(FINAL_DIR, showWarnings = FALSE, recursive = TRUE)

# Okabe-Ito palette: Control = index 6, Exercise choice = index 5
pal_okabe <- as.character(cols4all::c4a("okabe"))
pal_treat  <- c("Control" = pal_okabe[6], "Exercise choice" = pal_okabe[5])

# Helper: extract F/chi-sq and p for a term from car::Anova or anova() result
extract_fp <- function(anv_obj, term_pattern = "condition|treatment") {
  tab <- as.data.frame(anv_obj)
  idx <- grep(term_pattern, rownames(tab), ignore.case = TRUE)[1]
  if (is.na(idx)) return(list(F = NA_real_, p = NA_real_, numDF = NA_real_, denDF = NA_real_))
  row <- tab[idx, , drop = FALSE]
  # F-test columns (lm or car::Anova ddf=KR/Satt: "F value" OR "F")
  has_F <- ("F value" %in% names(row)) || ("F" %in% names(row))
  if (has_F) {
    Fcol <- if ("F value" %in% names(row)) "F value" else "F"
    Fval  <- as.numeric(row[[Fcol]])
    pval  <- as.numeric(row[["Pr(>F)"]])
    numDF <- if ("NumDF" %in% names(row)) as.numeric(row[["NumDF"]])
             else if ("Df" %in% names(row)) as.numeric(row[["Df"]])
             else NA_real_
    denDF <- if ("DenDF" %in% names(row)) as.numeric(row[["DenDF"]])
             else if ("Df.res" %in% names(row)) as.numeric(row[["Df.res"]])
             else NA_real_
  # Chi-square columns (car::Anova on lmer default)
  } else if ("Chisq" %in% names(row)) {
    Fval  <- as.numeric(row[["Chisq"]])
    pval  <- as.numeric(row[["Pr(>Chisq)"]])
    numDF <- as.numeric(row[["Df"]])
    denDF <- NA_real_
  } else {
    return(list(F = NA_real_, p = NA_real_, numDF = NA_real_, denDF = NA_real_))
  }
  list(F = Fval, p = pval, numDF = numDF, denDF = denDF)
}

fmt_anova_cap <- function(fp) {
  if (length(fp$F) == 0 || is.na(fp$F)) return("")
  p_txt <- if (!is.na(fp$p) && fp$p < 0.001) "< 0.001" else fmt3(fp$p)
  if (!is.na(fp$denDF)) {
    # Partial eta-squared from the term's own F and df (F-tests only; the
    # chi-square branch of extract_fp() returns denDF = NA and is skipped).
    eta <- (fp$F * fp$numDF) / (fp$F * fp$numDF + fp$denDF)
    sprintf("F<sub>%d,%d</sub> = %s, p = %s, &eta;<sup>2</sup><sub>p</sub> = %s",
            # denDF rounded (2026-08-28): the manuscript now gives integer df and the
            # exact Kenward-Roger values live in the Supplementary Materials.
            round(fp$numDF), round(fp$denDF), fmt_Fstat(fp$F), p_txt, fmt_es3(eta))
  } else {
    sprintf("F<sub>%d</sub> = %s, p = %s", round(fp$numDF), fmt_Fstat(fp$F), p_txt)
  }
}

# Legacy theme shared between cortisol and monoamine plots
theme_legacy <- function(base_size = 14) {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      panel.grid    = element_blank(),
      axis.line     = element_line(colour = "black", linewidth = 1.0),
      axis.title.x  = element_blank(),
      axis.text.x   = element_text(size = base_size + 4, face = "bold"),
      axis.text.y   = element_text(size = base_size),
      axis.title.y  = element_text(size = base_size + 4, face = "bold",
                                   angle = 90, vjust = 0.5,
                                   margin = margin(r = 12)),
      plot.caption  = ggtext::element_markdown(size = 17, hjust = 0,
                                              face = "italic", lineheight = 1.3,
                                              margin = margin(t = 14)),
      legend.position = "none",
      plot.margin   = margin(10, 10, 30, 10)
    )
}

# ---- A. Cortisol: legacy jitter plot by condition ---------------------------
cort_leg_dat <- cort_sex %>%
  mutate(
    treatment_lbl = factor(
      dplyr::recode(as.character(condition),
                    "control" = "Control", "treat" = "Exercise choice",
                    .default  = as.character(condition)),
      levels = c("Control", "Exercise choice")
    )
  )

cort_leg_summ <- cort_leg_dat %>%
  group_by(treatment_lbl) %>%
  summarise(mean = mean(plasma_cortisol, na.rm = TRUE),
            sd   = sd(plasma_cortisol,   na.rm = TRUE),
            n    = sum(!is.na(plasma_cortisol)),
            .groups = "drop") %>%
  mutate(sem = ifelse(n > 0, sd / sqrt(n), NA_real_),
         y_N = 0)

max_y_c    <- max(cort_leg_dat$plasma_cortisol, na.rm = TRUE)
min_y_c    <- min(cort_leg_dat$plasma_cortisol, na.rm = TRUE)
rng_c      <- max(max_y_c - min_y_c, max_y_c * 0.1, 1)
bar_y_c    <- max_y_c + 0.10 * rng_c
cap_h_c    <- 0.025 * rng_c
star_y_c   <- bar_y_c + 0.04 * rng_c
y_max_c    <- max_y_c + 0.35 * rng_c
y_min_c    <- -0.18 * rng_c   # space below 0 for N= labels
cort_leg_summ$y_N <- -0.10 * rng_c

# Alias §4 KR ANOVA to legacy name expected by extract_fp()
anv <- if (exists("cort_anova_kr") && !is.null(cort_anova_kr)) {
  cort_anova_kr %>% tibble::column_to_rownames("term")
} else NULL
cort_fp  <- extract_fp(anv, "condition|treatment")
cort_cap <- fmt_anova_cap(cort_fp)

p_cort_leg <- ggplot(
  cort_leg_dat,
  aes(x = treatment_lbl, y = plasma_cortisol, color = treatment_lbl)
) +
  geom_jitter(width = 0.30, size = 4.5, alpha = 0.9) +
  geom_errorbar(
    data = cort_leg_summ,
    aes(x = treatment_lbl, y = mean, ymin = mean - sem, ymax = mean + sem),
    width = 0.15, linewidth = 1.4, colour = "black", inherit.aes = FALSE
  ) +
  geom_segment(
    data = cort_leg_summ,
    aes(x    = as.numeric(treatment_lbl) - 0.20,
        xend = as.numeric(treatment_lbl) + 0.20,
        y = mean, yend = mean),
    linewidth = 1.5, colour = "black", inherit.aes = FALSE
  ) +
  geom_point(
    data = cort_leg_summ,
    aes(x = treatment_lbl, y = mean),
    colour = "black", size = 2.5, inherit.aes = FALSE
  ) +
  geom_text(
    data = cort_leg_summ,
    # 2026-08-21 (MV): switched to lowercase n, then back to capital N later the
    # same day at the author's request -- the manuscript's convention is capital
    # N for EVERY sample size, group-level or fish-level, and the text was
    # changed to match. MV's comment asked for consistency, not for a particular
    # case; the author picked N.
    aes(x = treatment_lbl, y = y_N, label = paste0("N = ", n)),
    size = 5.5, inherit.aes = FALSE
  ) +
  scale_color_manual(values = pal_treat, guide = "none") +
  labs(x = NULL, y = "[Cortisol] (ng/ml)", caption = cort_cap) +
  scale_y_continuous(
    limits = c(y_min_c, y_max_c),
    breaks = function(x) seq(0, ceiling(x[2] / 10) * 10, by = 10),
    expand = c(0, 0)
  ) +
  theme_legacy(base_size = 14) +
  theme(
    aspect.ratio = 1.5,
    # Figure 13 styling: caption padding (14->40), 2x frame (margin),
    # 1.5x x-axis tick labels (14->21), 2x y-axis tick labels (14->28),
    # 2x y-axis title (18->36), 2x y-axis title spacing (r=12->r=24).
    # Caption size reduced from the original 2x (17->34) to 20 (2026-08-09):
    # even on this figure's 3740 px canvas, the one-line caption (F, p,
    # eta2p, no CI/wrap) overran the panel's left-aligned text box at 34pt
    # -- confirmed by direct visual inspection of the rendered PNG.
    plot.caption = ggtext::element_markdown(size = 20, hjust = 0, face = "italic",
                                             lineheight = 1.3,
                                             margin = margin(t = 40)),
    plot.margin  = margin(20, 20, 80, 20),
    axis.text.x  = element_text(size = 21),
    axis.text.y  = element_text(size = 28),
    axis.title.y = element_text(size = 36, face = "bold", angle = 90, vjust = 0.5,
                                 margin = margin(r = 24))
  )

# Conditional significance bar + star (p < 0.05 only)
if (!is.na(cort_fp$p) && cort_fp$p < 0.05) {
  p_cort_leg <- p_cort_leg +
    annotate("segment", x = 1, xend = 2, y = bar_y_c, yend = bar_y_c,
             linewidth = 0.7, colour = "black") +
    annotate("segment",
             x    = c(1, 2), xend = c(1, 2),
             y    = c(bar_y_c - cap_h_c, bar_y_c - cap_h_c),
             yend = c(bar_y_c, bar_y_c),
             linewidth = 0.7, colour = "black") +
    annotate("text", x = 1.5, y = star_y_c, label = "*",
             size = 12, fontface = "bold")
}

ggsave(file.path(FINAL_DIR, "Cort_jitterplot_legacy.png"),
       p_cort_leg, width = 7.48, height = 11.22, units = "in", dpi = 500)
message("Saved: Cort_jitterplot_legacy.png")

# ---- B. Monoamine: per-analyte legacy plots (4-area grid, SEM) --------------
analyte_y_lbl <- function(an) {
  switch(an,
    "5-HT"        = "[5-HT] (ng/mg)",
    "5-HIAA"      = "[5-HIAA] (ng/mg)",
    "5-HIAA/5-HT" = "[5-HIAA]/[5-HT]",
    "DA"          = "[DA] (ng/mg)",
    "DOPAC"       = "[DOPAC] (ng/mg)",
    "DOPAC/DA"    = "[DOPAC]/[DA]",
    "NE"          = "[NE] (ng/mg)",
    paste0("[", an, "]")
  )
}

# ANOVA caption lookup: per-cell (area-specific) models — primary
# cell_anova_caps built above from mono_cell_models; used for area-specific sig stars.
# Omnibus caps retained for any exploratory reference (see exploratory §5 note).
mono_anova_caps_omnibus <- purrr::map_dfr(mono_models, function(x) {
  fp <- extract_fp(x$anova_kr, "treatment|condition")
  tibble(analyte_key = x$analyte_key,
         analyte     = analyte_label[x$analyte_key],
         p_val       = fp$p,
         caption     = fmt_anova_cap(fp))
})

make_mono_area_plot <- function(d_ar, an_label, is_sig, cap_txt, show_cap, star_size = 9) {
  summ_ar <- d_ar %>%
    group_by(treatment_lbl) %>%
    summarise(mean = mean(value, na.rm = TRUE),
              sd   = sd(value, na.rm = TRUE),
              n    = sum(!is.na(value)),
              .groups = "drop") %>%
    mutate(sem = ifelse(n > 0, sd / sqrt(n), NA_real_)) %>%
    filter(n > 0)
  if (nrow(summ_ar) == 0) return(NULL)

  y_max_d <- max(c(d_ar$value, summ_ar$mean + summ_ar$sem), na.rm = TRUE)
  y_min_d <- min(c(d_ar$value, summ_ar$mean - summ_ar$sem), na.rm = TRUE)
  if (!is.finite(y_min_d)) y_min_d <- min(d_ar$value, na.rm = TRUE)
  if (!is.finite(y_max_d)) y_max_d <- max(d_ar$value, na.rm = TRUE)
  y_span  <- max(y_max_d - y_min_d, abs(y_max_d) * 0.01, 1e-6)
  y_N     <- y_min_d - 0.22 * y_span
  y_line  <- y_max_d + 0.10 * y_span
  cap_h   <- 0.025 * y_span
  y_star  <- y_line + 0.04 * y_span
  y_lo    <- y_min_d - 0.32 * y_span
  y_hi    <- y_max_d + 0.42 * y_span

  p_ar <- ggplot(d_ar, aes(x = treatment_lbl, y = value, color = treatment_lbl)) +
    geom_jitter(width = 0.25, size = 2.5) +
    geom_errorbar(
      data = summ_ar,
      aes(x = treatment_lbl, y = mean, ymin = mean - sem, ymax = mean + sem),
      width = 0.10, colour = "black", inherit.aes = FALSE
    ) +
    geom_point(
      data = summ_ar, aes(x = treatment_lbl, y = mean),
      colour = "black", size = 3.2, inherit.aes = FALSE
    ) +
    geom_text(
      data = summ_ar,
      # 2026-08-21: capital N -- see the note on the cortisol panel above. The
      # manuscript uses N for every sample size.
      aes(x = treatment_lbl, y = y_N, label = paste0("N = ", n)),
      inherit.aes = FALSE, vjust = 1, size = 5.2
    ) +
    scale_color_manual(values = pal_treat, guide = "none") +
    labs(x = NULL,
         y = analyte_y_lbl(an_label),
         caption = if (show_cap) cap_txt else "") +
    scale_y_continuous(limits = c(y_lo, y_hi),
                       breaks = ~ pretty(c(0, .x[2]))) +
    theme_legacy(base_size = 13) +
    theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
          # theme_legacy()'s plot.caption is sized for a single wide plot
          # (cortisol); in this 4-panel-per-row grid each panel is much
          # narrower, so the one-line caption (F, p, eta2p, no CI/wrap since
          # 2026-08-09) needs a smaller size to actually fit instead of
          # overrunning the panel edge. size=10 still clipped when F was
          # very small (4-sig-fig fmt_F makes e.g. 0.0009258 a long string,
          # confirmed by direct visual check on DA/VD) -- dropped further.
          plot.caption = ggtext::element_markdown(size = 8, hjust = 0,
                                                    face = "italic",
                                                    margin = margin(t = 10)))

  if (is_sig) {
    p_ar <- p_ar +
      annotate("segment", x = 1, xend = 2, y = y_line, yend = y_line,
               linewidth = 0.7, colour = "black") +
      annotate("segment",
               x    = c(1, 2), xend = c(1, 2),
               y    = c(y_line - cap_h, y_line - cap_h),
               yend = c(y_line, y_line),
               linewidth = 0.7, colour = "black") +
      annotate("text", x = 1.5, y = y_star, label = "*",
               size = star_size, fontface = "bold")
  }
  p_ar
}

make_mono_legacy_panel <- function(akey) {
  an_label <- analyte_label[akey]

  d <- mono_long %>%
    filter(analyte_key == akey, is.finite(value), !is.na(treatment)) %>%
    mutate(
      treatment_lbl = factor(
        dplyr::recode(as.character(treatment),
                      "control" = "Control", "treat" = "Exercise choice",
                      .default  = as.character(treatment)),
        levels = c("Control", "Exercise choice")
      ),
      area = factor(area, levels = expected_areas)
    )
  if (nrow(d) == 0) return(NULL)

  area_lvls  <- intersect(expected_areas, unique(as.character(d$area)))
  area_plots <- purrr::imap(
    stats::setNames(area_lvls, area_lvls),
    function(ar, i) {
      d_ar <- filter(d, area == ar)
      if (nrow(d_ar) == 0) return(NULL)
      # Area-specific caption from per-cell model; uncorrected, raw p only.
      cap_row <- cell_anova_caps %>% filter(analyte_key == akey, area == ar)
      cap_txt <- if (nrow(cap_row)) cap_row$caption[1] else ""
      is_sig  <- if (nrow(cap_row)) isTRUE(cap_row$is_sig_q[1]) else FALSE
      p <- make_mono_area_plot(
        d_ar, an_label, is_sig, cap_txt,
        show_cap = (ar == area_lvls[length(area_lvls)])  # caption on last panel only
      )
      p + ggtitle(ar)
    }
  )
  area_plots <- purrr::compact(area_plots)
  if (length(area_plots) == 0) return(NULL)

  patchwork::wrap_plots(area_plots, nrow = 1) +
    patchwork::plot_annotation(
      title  = an_label,
      theme  = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
    )
}

FIG_DIR <- file.path(OUT, "Figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

for (akey in analyte_cols) {
  p_leg <- tryCatch(make_mono_legacy_panel(akey),
                    error = function(e) {
                      message(sprintf("Legacy mono plot %s failed: %s", akey, conditionMessage(e)))
                      NULL
                    })
  if (is.null(p_leg)) next
  png_name <- sprintf("Mono_%s_jitterplot_legacy.png",
                      gsub("[^A-Za-z0-9]", "_", analyte_label[akey]))
  ggsave(file.path(FIG_DIR, png_name),
         p_leg, width = 14, height = 5.5, units = "in", dpi = 300)
  message(sprintf("Saved to Figures/: %s", png_name))
}

# ---- C. DM stacked composite plots → B3_Final_Report/ ----------------------
# Serotonin group: 5-HT / 5-HIAA / 5-HIAA:5-HT stacked (A/B/C)
# Catecholamine group: DA / DOPAC / DOPAC:DA stacked (A/B/C)
# One plot per analyte → DM area only → cowplot stack with panel labels

make_dm_single_plot <- function(akey) {
  an_label <- analyte_label[akey]
  # Uncorrected: significance follows raw p (no multiple-comparison adjustment).
  cap_row  <- cell_anova_caps %>% filter(analyte_key == akey, area == "DM")
  cap_txt  <- if (nrow(cap_row)) cap_row$caption[1] else ""
  is_sig   <- if (nrow(cap_row)) isTRUE(cap_row$is_sig_q[1]) else FALSE

  d <- mono_long %>%
    filter(analyte_key == akey, area == "DM", is.finite(value), !is.na(treatment)) %>%
    mutate(
      treatment_lbl = factor(
        dplyr::recode(as.character(treatment),
                      "control" = "Control", "treat" = "Exercise choice",
                      .default  = as.character(treatment)),
        levels = c("Control", "Exercise choice")
      )
    )
  if (nrow(d) == 0) return(NULL)
  p <- make_mono_area_plot(d, an_label, is_sig, cap_txt, show_cap = TRUE, star_size = 12)
  # Figure 13/14 styling (the DM stacks are the only monoamine figures that go
  # into the manuscript): 1.5x caption-to-plot padding (14 -> 21), 2x white
  # frame around plot (10/30/10 -> 20/60/20). Caption size reduced from the
  # original 2x (17 -> 34) to 20 (2026-08-09) and matched to the cortisol
  # figure: at 34pt the one-line caption (F, p, eta2p) ran off the 7.48-in
  # canvas and the effect size was clipped entirely -- confirmed by direct
  # visual inspection of the rendered PNG.
  p + ggplot2::theme(
    plot.caption = ggtext::element_markdown(size = 20, hjust = 0, face = "italic",
                                             lineheight = 1.3,
                                             margin = ggplot2::margin(t = 21)),
    plot.margin  = ggplot2::margin(20, 20, 60, 20)
  )
}

sero_keys <- c("ht_5", "hiaa_5", "hiaa_5_ratio")   # 5-HT, 5-HIAA, 5-HIAA/5-HT
cate_keys <- c("da",   "dopac",  "dopac_da_ratio")  # DA, DOPAC, DOPAC/DA

dm_sero_plots <- purrr::compact(purrr::map(sero_keys, make_dm_single_plot))
dm_cate_plots <- purrr::compact(purrr::map(cate_keys, make_dm_single_plot))

# Stack each group with A/B/C panel labels — portrait ratio 700:1330 from legacy
# 3740 px / 500 dpi = 7.48 in wide; height = 7.48 * (1330/700) = 14.22 in
DM_W <- 7.48
DM_H <- round(DM_W * (1330 / 700), 2)  # 14.22 in

stack_dm <- function(plot_list, fname) {
  if (length(plot_list) == 0) return(invisible(NULL))
  stk <- cowplot::plot_grid(
    plotlist       = plot_list,
    ncol           = 1,
    align          = "v",
    rel_heights    = rep(1.2, length(plot_list)),
    # 2026-08-21 (MV): two fixes in one place. The dot is gone -- every other
    # figure in the manuscript tags panels "A", "B", "C" without one, and this
    # was the only family that did not. And label_y dropped 1.05 -> 1.00: at
    # 1.05 the tag sits 5% of a panel height ABOVE its panel, which is fine for
    # panels B and C (they have the panel above to bleed into) but puts panel
    # A's tag off the top edge of the canvas, where it is clipped away entirely.
    # That is why every rendered version of Figures 8 and 9 showed "B." and
    # "C." but no "A" at all.
    labels         = LETTERS[seq_along(plot_list)],
    label_size     = 20,
    label_fontface = "bold",
    label_x        = 0.16,
    label_y        = 1.00,
    hjust          = 0,
    vjust          = 1.0
  )
  ggsave(file.path(FINAL_DIR, fname), stk,
         width = DM_W, height = DM_H, units = "in", dpi = 500)
  message(sprintf("Saved to B3_Final_Report/: %s", fname))
  stk
}

stk_sero <- stack_dm(dm_sero_plots, "DM_serotonin_group_SEM.png")
stk_cate <- stack_dm(dm_cate_plots, "DM_catecholamine_group_SEM.png")

# Mirror all three manu graphs to shared folder
ALL_MANU <- file.path(PROJECT_ROOT, "all_manu_graphs")
dir.create(ALL_MANU, showWarnings = FALSE, recursive = TRUE)

if (!is.null(stk_sero))
  ggsave(file.path(ALL_MANU, "DM_serotonin_group_SEM.png"),
         stk_sero, width = DM_W, height = DM_H, units = "in", dpi = 500)
if (!is.null(stk_cate))
  ggsave(file.path(ALL_MANU, "DM_catecholamine_group_SEM.png"),
         stk_cate, width = DM_W, height = DM_H, units = "in", dpi = 500)
ggsave(file.path(ALL_MANU, "Cort_jitterplot_legacy.png"),
       p_cort_leg, width = 7.48, height = 11.22, units = "in", dpi = 500)
message(sprintf("Saved 3 manu graphs → %s", ALL_MANU))

# Also save to dated folder manu_graphs_26.05.2026 (inside FINAL_DIR and ALL_MANU)
DATED_DIR_B3   <- file.path(FINAL_DIR, "manu_graphs_26.05.2026")
DATED_DIR_MANU <- file.path(ALL_MANU,  "manu_graphs_26.05.2026")
dir.create(DATED_DIR_B3,   showWarnings = FALSE, recursive = TRUE)
dir.create(DATED_DIR_MANU, showWarnings = FALSE, recursive = TRUE)
if (!is.null(stk_sero)) {
  ggsave(file.path(DATED_DIR_B3,   "DM_serotonin_group_SEM.png"),
         stk_sero, width = DM_W, height = DM_H, units = "in", dpi = 500)
  ggsave(file.path(DATED_DIR_MANU, "DM_serotonin_group_SEM.png"),
         stk_sero, width = DM_W, height = DM_H, units = "in", dpi = 500)
}
if (!is.null(stk_cate)) {
  ggsave(file.path(DATED_DIR_B3,   "DM_catecholamine_group_SEM.png"),
         stk_cate, width = DM_W, height = DM_H, units = "in", dpi = 500)
  ggsave(file.path(DATED_DIR_MANU, "DM_catecholamine_group_SEM.png"),
         stk_cate, width = DM_W, height = DM_H, units = "in", dpi = 500)
}
ggsave(file.path(DATED_DIR_B3,   "Cort_jitterplot_legacy.png"),
       p_cort_leg, width = 7.48, height = 11.22, units = "in", dpi = 500)
ggsave(file.path(DATED_DIR_MANU, "Cort_jitterplot_legacy.png"),
       p_cort_leg, width = 7.48, height = 11.22, units = "in", dpi = 500)
message(sprintf("Dated copies saved → %s", DATED_DIR_B3))

cat("\nDM composite plots → B3_Final_Report/\nAll-area mono plots → Figures/\n")

# ---- 14. Word report (template-style: §1 → §9) -----------------------------
mk_ft <- function(df, digits = 3) {
  ft <- flextable::flextable(df)
  ft <- flextable::colformat_double(ft, digits = digits)
  ft <- flextable::set_table_properties(ft, width = 1, layout = "autofit")
  ft
}

add_heading <- function(doc, text, level = 1) {
  body_add_par(doc, text, style = paste0("heading ", min(level, 3L)))
}

# Pre-compute CV sentence values (kept from prior version for §1 methods)
cv_x <- sum(intra_assay_summary$weighted_mean_cv * intra_assay_summary$n_samples, na.rm = TRUE) /
        sum(intra_assay_summary$n_samples, na.rm = TRUE)
if (nrow(inter_assay_cv_repeats) > 0 &&
    any(is.finite(inter_assay_cv_repeats$cv_inter_pct))) {
  cv_y <- mean(inter_assay_cv_repeats$cv_inter_pct, na.rm = TRUE)
} else {
  cv_y <- inter_assay_cv_controls$cv_inter_pct_OD[inter_assay_cv_controls$control == "B0"]
}
cv_sentence <- sprintf(
  paste0("Plasma cortisol concentrations were measured in duplicate using a commercial ELISA kit. ",
         "Mean intra-assay CV (duplicate wells per plate) was %.1f%%. ",
         "Inter-assay CV (internal control across plates) was %.1f%%."),
  cv_x, cv_y)

# Build document =============================================================
doc <- read_docx() %>%
  add_heading("B3 plasma cortisol & monoamines — analysis report", 1) %>%
  body_add_par(sprintf("Generated: %s", Sys.time()), style = "Normal") %>%
  body_add_par("Scope: B3 / batch 3 only. Mucus and MU* plates excluded.", style = "Normal")

# §1 Study Design + Data Summary --------------------------------------------
doc <- doc %>%
  add_heading("§1. Study Design and Data Summary", 1) %>%
  add_heading("1.1 Random effects reference", 2) %>%
  body_add_flextable(mk_ft(re_reference)) %>%
  add_heading("1.2 Glossary", 2) %>%
  body_add_par(paste0(
    "EMM = estimated marginal mean. EMM ± SE = estimated marginal mean ± standard error. ",
    "KR = Kenward-Roger F-test denominator-df approximation (primary inference). ",
    "Satt = Satterthwaite (sensitivity row). No adjustment for multiple comparisons is applied; ",
    "all p-values are raw. CLD = compact letter display. β = fixed-effect estimate on the log(concentration) scale ",
    "(treatment-coding is sum-to-zero, so β represents half the control–treat difference; ",
    "EMM contrasts on the response scale are reported as fold-changes when applicable)."),
    style = "Normal")

# §2 Statistical decision log -----------------------------------------------
doc <- doc %>%
  add_heading("§2. Statistical decision log", 1) %>%
  body_add_flextable(mk_ft(decision_log))

# §1 Methods (assay) — short prose to keep template feel
doc <- doc %>%
  add_heading("Assay quality (cortisol ELISA)", 1) %>%
  body_add_par(cv_sentence, style = "Normal") %>%
  add_heading("Per-plate standard curve quality (4PL)", 2) %>%
  body_add_flextable(mk_ft(plate_curve_quality)) %>%
  add_heading("Intra-assay CV per plate", 2) %>%
  body_add_flextable(mk_ft(intra_assay_summary)) %>%
  add_heading("Inter-assay CV — B0 / NSB controls", 2) %>%
  body_add_flextable(mk_ft(inter_assay_cv_controls)) %>%
  add_heading("Inter-assay CV — Standards", 2) %>%
  body_add_flextable(mk_ft(inter_assay_cv_standards)) %>%
  add_heading("Inter-assay CV — B3 samples repeated across plates", 2) %>%
  body_add_flextable(mk_ft(inter_assay_cv_repeats)) %>%
  body_add_gg(p_curve, width = 6.5, height = 4) %>%
  body_add_gg(p_cv,    width = 6.5, height = 4)

# §3 Sex-balance check ------------------------------------------------------
doc <- doc %>%
  add_heading("§3. Sex-balance check (treatment groups)", 1) %>%
  body_add_par(paste0(
    "Sex composition compared between treatment groups in two cohorts: full cortisol cohort ",
    "(n = ", nrow(cort_full), ") and monoamine subsample. Tests: overall χ²/Fisher; ",
    "Cochran-Mantel-Haenszel stratified by trial; GLMM I(sex=M)~condition+(1|trial)."),
    style = "Normal") %>%
  body_add_flextable(mk_ft(sex_balance))

# §4 Part 1: Treatment main effect — cortisol -----------------------------
doc <- doc %>%
  add_heading("§4. Treatment main effect — Part 1 (aggregated)", 1) %>%
  add_heading(sprintf("4.1 Cortisol [model = %s]", cort_model_type), 2) %>%
  add_heading("4.1.1 Normality / Levene", 3) %>%
  body_add_flextable(mk_ft(norm_lev_cort)) %>%
  add_heading("4.1.2 RE selection (AICc)", 3) %>%
  body_add_flextable(mk_ft(cort_re_aicc_table)) %>%
  add_heading("4.1.3 Model R²", 3) %>%
  body_add_flextable(mk_ft(cort_r2))
if (!is.null(cort_anova_kr))
  doc <- doc %>%
    add_heading("4.1.4 ANOVA (KR F-test, primary)", 3) %>%
    body_add_flextable(mk_ft(cort_anova_kr))
if (!is.null(cort_anova_satt))
  doc <- doc %>%
    add_heading("4.1.5 ANOVA (Satterthwaite, sensitivity)", 3) %>%
    body_add_flextable(mk_ft(cort_anova_satt))
if (!is.null(cort_emm_df)) {
  doc <- doc %>%
    add_heading("4.1.6 Estimated marginal means", 3) %>%
    body_add_flextable(mk_ft(cort_emm_df)) %>%
    add_heading("4.1.6 Pairwise (Tukey)", 3) %>%
    body_add_flextable(mk_ft(cort_pairs_df))
  if (!is.null(cort_cld_df))
    doc <- doc %>% add_heading("4.1.6 CLD", 3) %>% body_add_flextable(mk_ft(cort_cld_df))
}
if (file.exists(diag_cort_dharma_path))
  doc <- doc %>%
    add_heading("4.1.7 DHARMa diagnostic", 3) %>%
    body_add_img(diag_cort_dharma_path, width = 6.5, height = 3.5)
doc <- doc %>% body_add_gg(p_cort, width = 6.5, height = 4)

# §4 (cont.) Monoamines — PRIMARY: per-analyte × per-area, grouped by brain area ----------
doc <- doc %>%
  add_heading("§4 (cont.) — Monoamines: treatment main effect per analyte × brain area", 1) %>%
  body_add_par(paste0(
    "Each analyte × brain-area combination is modelled independently with treatment as the ",
    "sole fixed effect: log(concentration) ~ treatment + RE. Random-effect candidates: ",
    "(1|date), (1|tank), (1|trial) — trial nests fish sampled from the same video, the ",
    "unit of replication for the treatment manipulation (M2/R1-6). When all three are ",
    "singular, the model reduces to a fixed-effects lm ",
    "— this is flagged in the RE-AICc table. No adjustment for multiple comparisons is applied; ",
    "raw p is reported for each analyte x brain-area cell. ",
    "Descriptive statistics are on the raw (untransformed) concentration scale; ",
    "model estimates (β, EMM, contrasts) are on the log scale."),
    style = "Normal")

for (ar in expected_areas) {
  doc <- add_heading(doc, sprintf("Brain area: %s", ar), 2)

  # Area-level treatment-effect summary table first (uncorrected, raw p only)
  area_bh <- cell_cross_summary %>% filter(area == ar)
  if (nrow(area_bh) > 0) {
    doc <- doc %>%
      add_heading(sprintf("%s — treatment effect summary (%d monoamines)", ar, nrow(area_bh)), 3) %>%
      body_add_par(sprintf(
        "No adjustment for multiple comparisons applied; %d raw treatment p-values shown for %s.",
        nrow(area_bh), ar), style = "Normal") %>%
      body_add_flextable(mk_ft(area_bh %>%
        dplyr::select(analyte, n, re_spec, beta, SE, CI_lo_95, CI_hi_95,
                      F_stat, df1, df2, p_raw, sig_raw)))
  }

  area_cells <- Filter(function(x) x$area == ar, mono_cell_models)
  for (x in area_cells) {
    doc <- add_heading(doc, sprintf("%s — %s", ar, x$analyte), 3)

    # Descriptive statistics (raw scale)
    doc <- doc %>%
      add_heading("Descriptive statistics (raw scale, by treatment)", 4) %>%
      body_add_flextable(mk_ft(x$descr)) %>%
      body_add_par(sprintf(
        "Model: %s | n = %d | RE: %s", x$model_formula, x$n, x$re_label),
        style = "Normal")

    # RE selection
    doc <- doc %>%
      add_heading("RE selection (AICc)", 4) %>%
      body_add_flextable(mk_ft(x$re_aicc_table))

    # R²
    doc <- doc %>%
      add_heading("Model R²", 4) %>%
      body_add_flextable(mk_ft(x$r2))

    # KR ANOVA (primary)
    if (!is.null(x$anova_kr)) {
      doc <- doc %>%
        add_heading("ANOVA — Kenward-Roger F-test (primary)", 4) %>%
        body_add_flextable(mk_ft(x$anova_kr))
    }

    # Satterthwaite ANOVA (sensitivity)
    if (!is.null(x$anova_satt)) {
      doc <- doc %>%
        add_heading("ANOVA — Satterthwaite (sensitivity)", 4) %>%
        body_add_flextable(mk_ft(x$anova_satt))
    }

    # EMM + Tukey + CLD
    if (!is.null(x$emm)) {
      doc <- doc %>%
        add_heading("Estimated marginal means (log scale)", 4) %>%
        body_add_flextable(mk_ft(x$emm$emm))
      if (!is.null(x$emm$pairs)) {
        doc <- doc %>%
          add_heading("Pairwise contrasts — control vs. treat (Tukey)", 4) %>%
          body_add_flextable(mk_ft(x$emm$pairs))
      }
      if (!is.null(x$emm$cld)) {
        doc <- doc %>%
          add_heading("Compact letter display (CLD)", 4) %>%
          body_add_flextable(mk_ft(x$emm$cld))
      }
    }

    # Effect size
    if (!is.null(x$eff_d)) {
      eff_df <- tryCatch(as.data.frame(x$eff_d), error = function(e) NULL)
      if (!is.null(eff_df) && nrow(eff_df) > 0) {
        doc <- doc %>%
          add_heading("Effect size (Cohen's d, log scale)", 4) %>%
          body_add_flextable(mk_ft(eff_df))
      }
    }

    # Normality / Levene
    doc <- doc %>%
      add_heading("Normality / Levene check", 4) %>%
      body_add_flextable(mk_ft(x$norm_lev))

    # DHARMa
    if (!is.na(x$dharma_path) && file.exists(x$dharma_path)) {
      doc <- doc %>%
        add_heading("DHARMa residual diagnostic", 4) %>%
        body_add_img(x$dharma_path, width = 6.5, height = 3.5)
    }
  }
}

# §5 Exploratory: omnibus treatment × area (for context only) -----------------
doc <- doc %>%
  add_heading("§5. Exploratory: omnibus treatment × area models (reference only)", 1) %>%
  body_add_par(paste0(
    "EXPLORATORY / NOT CONFIRMATORY. These omnibus models (log_v ~ treatment * area + RE) ",
    "pool all four brain areas into a single model per analyte. They are retained here as ",
    "descriptive context for cross-area patterns. Primary inference uses the per-cell models ",
    "reported in §4. BH correction does NOT apply to these exploratory results."),
    style = "Normal") %>%
  add_heading("5.1 RE selection (omnibus, per analyte)", 2) %>%
  body_add_flextable(mk_ft(mono_re_aicc)) %>%
  add_heading("5.2 Model R² (omnibus, per analyte)", 2) %>%
  body_add_flextable(mk_ft(mono_r2)) %>%
  add_heading("5.3 KR ANOVA — treatment × area interaction (omnibus, per analyte)", 2) %>%
  body_add_flextable(mk_ft(mono_anova_kr)) %>%
  add_heading("5.4 EMM — treatment marginal over area (omnibus)", 2) %>%
  body_add_flextable(mk_ft(mono_emm_agg)) %>%
  add_heading("5.5 Stratified EMM — treatment within area (omnibus)", 2) %>%
  body_add_flextable(mk_ft(mono_emm_str)) %>%
  add_heading("5.6 Stratified pairwise contrasts (omnibus)", 2) %>%
  body_add_flextable(mk_ft(mono_pairs_str)) %>%
  add_heading("5.7 Stratified CLD per area (omnibus)", 2) %>%
  body_add_flextable(mk_ft(mono_cld_str))

# §6 Cross-summary BH -------------------------------------------------------
doc <- doc %>%
  add_heading("§6. Cross-analysis summary (BH-corrected; cortisol uncorrected)", 1) %>%
  add_heading("6.1 Per-cell BH summary — monoamines (primary; 4 families of 7)", 2) %>%
  body_add_par(paste0(
    "BH-FDR applied independently within each brain area over the 7 monoamine treatment ",
    "p-values (per-cell model: log_v ~ treatment + RE). Four separate BH families ",
    "(one per area). Cortisol is reported separately, uncorrected."),
    style = "Normal") %>%
  body_add_flextable(mk_ft(cell_cross_summary)) %>%
  add_heading("6.2 Cortisol (uncorrected) and monoamine omnibus aggregated (exploratory)", 2) %>%
  body_add_par(paste0(
    "Cortisol p-value is uncorrected (not part of the monoamine family). ",
    "Monoamine row uses the omnibus model marginal treatment p (exploratory reference only)."),
    style = "Normal") %>%
  body_add_flextable(mk_ft(cross_summary))
if (!is.null(stratified_summary))
  doc <- doc %>%
    add_heading("6.3 Omnibus stratified BH (treatment × area, exploratory)", 2) %>%
    body_add_flextable(mk_ft(stratified_summary))
if (file.exists(file.path(FIG_DIR, "emm_forest_treatment.png")))
  doc <- doc %>%
    add_heading("6.4 Treatment effect forest plot (omnibus betas, log scale)", 2) %>%
    body_add_img(file.path(FIG_DIR, "emm_forest_treatment.png"), width = 6.5, height = 4.5)

# §7 Pseudoreplication ------------------------------------------------------
doc <- doc %>%
  add_heading("§7. Pseudoreplication assessment (M2/R1-6)", 1) %>%
  body_add_par(paste0(
    "Cortisol winner model refit with each non-winning RE; treatment β and SE compared. ",
    "Flag = TRUE if SE shifts > 30% or β changes sign."), style = "Normal")
if (!is.null(pseudorep_table))
  doc <- doc %>% body_add_flextable(mk_ft(pseudorep_table))
doc <- doc %>%
  add_heading("7.1 Monoamine per-cell models (primary; Table 2 / Figures 11-12 source)", 2) %>%
  body_add_par(paste0(
    "Each cell's winning model (§4 primary) refit with every non-winning RE from ",
    "{(1|date), (1|tank), (1|trial)}; treatment β and SE compared across specifications. ",
    "`trial` (added for M2/R1-6) nests fish sampled together from the same video — the ",
    "unit of replication for the treatment manipulation. Flag = TRUE if SE shifts > 30% ",
    "or β changes sign; cells fit as lm (all RE candidates singular) are not probed."),
    style = "Normal")
if (nrow(mono_pseudorep_table) > 0)
  doc <- doc %>% body_add_flextable(mk_ft(mono_pseudorep_table))

# §7b Order/day-of-testing effect (M4/R2-5) ---------------------------------
doc <- doc %>%
  add_heading("§7b. Trial order / day-of-testing effect (M4/R2-5)", 1) %>%
  body_add_par(paste0(
    "Each model's winning fixed+random specification (from §4) is refit with ",
    "order_idx (chronological rank of testing date, 1-8) and its interaction with ",
    "treatment added: y ~ treatment * order_idx + RE. Rows below show the order_idx ",
    "main-effect and treatment:order_idx interaction terms only (KR F-test)."),
    style = "Normal") %>%
  body_add_par(paste0(
    "IMPORTANT LIMITATION: order_idx is not identifiable independently of tank or of ",
    "fish_density in this design. Between tanks, each tank occupies its own exclusive, ",
    "non-overlapping pair of testing days, so the day sequence is a perfect function of ",
    "tank identity (aliased with the tank random effect). Within a tank, the two testing ",
    "days always run fish_density = {16,12} then {8,4}, and the first trial each day always ",
    "has the higher of that day's two densities — order_idx and fish_density are related ",
    "by an exact affine transform, so an order test and a density test are mathematically ",
    "the same test here. The terms below should be read as \"order-or-density\", not as ",
    "independent evidence of chronological drift."), style = "Normal") %>%
  body_add_flextable(mk_ft(order_effect_summary))

# §8 Descriptive statistics -------------------------------------------------
doc <- doc %>%
  add_heading("§8. Descriptive statistics", 1) %>%
  add_heading("8.1 Cortisol (by treatment)", 2) %>%
  body_add_flextable(mk_ft(read.csv(file.path(DATA_DIR, "descriptive_cortisol.csv")))) %>%
  add_heading("8.2 Monoamines (by analyte × area × treatment)", 2) %>%
  body_add_flextable(mk_ft(mono_descr))

# §9 Secondary exploratory: sex effects -------------------------------------
doc <- doc %>%
  add_heading("§9. Secondary exploratory: sex effects (NON-CONFIRMATORY)", 1) %>%
  body_add_par(paste0(
    "Exploratory; not pre-registered; no FDR correction; hypothesis-generating only. ",
    "Sex enters the cortisol model as condition × sex interaction; for monoamines as an ",
    "additive covariate (small n=37 fish prevents stable interaction estimates)."),
    style = "Normal")
if (!is.null(sex_secondary_cortisol))
  doc <- doc %>%
    add_heading("9.1 Cortisol — condition × sex (KR ANOVA)", 2) %>%
    body_add_flextable(mk_ft(sex_secondary_cortisol))
if (!is.null(cort_sex_emm))
  doc <- doc %>%
    add_heading("9.1b Cortisol — EMM by condition × sex", 2) %>%
    body_add_flextable(mk_ft(cort_sex_emm))
if (nrow(sex_secondary_monoamines) > 0)
  doc <- doc %>%
    add_heading("9.2 Monoamines — sex covariate (KR ANOVA per analyte)", 2) %>%
    body_add_flextable(mk_ft(sex_secondary_monoamines))
if (exists("sex_cell_df") && nrow(sex_cell_df) > 0)
  doc <- doc %>%
    add_heading("9.3 Monoamines — per-cell sex effect (analyte x area)", 2) %>%
    body_add_flextable(mk_ft(sex_cell_df))

# Cortisol × monoamine correlations (exploratory) ---------------------------
doc <- doc %>%
  add_heading("Appendix A. Cortisol × monoamine correlations (exploratory)", 1) %>%
  body_add_flextable(mk_ft(cort_mono_corr))

# Audit appendix ------------------------------------------------------------
doc <- doc %>%
  add_heading("Appendix B. Audit", 1) %>%
  body_add_par("Sample-level matching audit between worksheet and assay files.", style = "Normal") %>%
  body_add_flextable(mk_ft(matching_audit)) %>%
  body_add_par("Replicate well counts per plate.", style = "Normal") %>%
  body_add_flextable(mk_ft(replicate_well_summary))

# Save report ---------------------------------------------------------------
report_path <- file.path(FINAL_DIR, "B3_cortisol_monoamines_report.docx")
tryCatch(
  print(doc, target = report_path),
  error = function(e) {
    ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
    report_path <<- file.path(FINAL_DIR,
                              sprintf("B3_cortisol_monoamines_report_%s.docx", ts))
    message(sprintf("Original docx locked — writing to %s", basename(report_path)))
    print(doc, target = report_path)
  }
)
message(sprintf("Saved: %s", basename(report_path)))

cat("\n=== Done ==========================================\n")
cat("Outputs in:\n  ", OUT, "\n", sep = "")
cat("Cortisol primary model: ", cort_model_type, "\n", sep = "")
cat("Monoamine RE summary:\n"); print(mono_re_summary)
cat("Cross-summary (BH on 7 monoamines; cortisol uncorrected):\n"); print(cross_summary)
cat("Sex-balance test (§3):\n"); print(sex_balance)

# ---- §15. Monoamine transformation info — standalone Word document ----------
# Source: mono_cell_models (primary per-analyte x per-area pipeline, §4).
# Transformation now selected adaptively per cell by select_transform():
#   prefer Box-Cox when Shapiro p >= 1.5x log; collapse to log() when lambda~0.
#   falls back to log() when Box-Cox unavailable.
# norm_lev carries: sw_raw, sw_log (Shapiro), levene_p_raw, levene_p_log (Levene).
# NOTE: when transforms differ across cells, beta/EMM are on different scales and
#   are not directly comparable without back-transformation.

get_cell_val <- function(nl, col) {
  if (col %in% names(nl) && length(nl[[col]]) > 0 && !is.na(nl[[col]][1]))
    sprintf("%.4f", nl[[col]][1])
  else "not available"
}

mono_transform_tbl <- purrr::map_dfr(mono_cell_models, function(x) {
  nl <- x$norm_lev
  post_lev_col <- if ("levene_p_log" %in% names(nl)) "levene_p_log" else "p"
  data.frame(
    Monoamine                      = x$analyte,
    Area                           = x$area,
    `Shapiro p (raw)`              = get_cell_val(nl, "sw_raw"),
    `Levene p (raw)`               = get_cell_val(nl, "levene_p_raw"),
    `Transformation selected`      = x$trans_label,
    `BC λ`                    = if (!is.na(x$trans_lambda)) sprintf("%.3f", x$trans_lambda)
                                     else "—",
    `Shapiro p (post-transform)`   = get_cell_val(nl, "sw_log"),
    `Levene p (post-transform)`    = get_cell_val(nl, post_lev_col),
    Rationale                      = x$trans_rationale,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
})

write.csv(mono_transform_tbl, file.path(DATA_DIR, "monoamine_transformation_info.csv"),
          row.names = FALSE)

n_bc <- sum(mono_transform_tbl[["Transformation selected"]] != "log(value)", na.rm = TRUE)
message(sprintf("Adaptive transform: %d cells use Box-Cox, %d use log().",
                n_bc, nrow(mono_transform_tbl) - n_bc))

ft_transform <- flextable::flextable(mono_transform_tbl)
ft_transform <- flextable::set_table_properties(ft_transform, width = 1, layout = "autofit")
ft_transform <- flextable::bold(ft_transform, part = "header")
ft_transform <- flextable::fontsize(ft_transform, size = 9, part = "all")
ft_transform <- flextable::bg(ft_transform, part = "header", bg = "#4472C4")
ft_transform <- flextable::color(ft_transform, part = "header", color = "white")

doc_transform <- officer::read_docx()
doc_transform <- officer::body_add_par(doc_transform,
  "Monoamine data-transformation summary (per analyte x brain area)", style = "heading 1")
doc_transform <- officer::body_add_par(doc_transform,
  paste0(
    "Per-cell (analyte x brain area) transformation summary for all primary-pipeline models. ",
    "Transformation is selected adaptively per cell: Box-Cox is preferred when its ",
    "Shapiro-Wilk p >= 1.5x the log() p and Box-Cox is available; otherwise log() is used. ",
    "Box-Cox with lambda near 0 (|lambda| < 0.1) is treated as log(). ",
    "Levene formula: raw: value ~ treatment; post-transform: trans(value) ~ treatment. ",
    "CAUTION: if transforms differ across cells, beta coefficients and EMM are on ",
    "different scales and are not directly comparable without back-transformation. ",
    "'not available' = test not computable for that cell."),
  style = "Normal")
doc_transform <- body_add_flextable(doc_transform, ft_transform)

transform_doc_path <- file.path(DATA_DIR, "monoamine_transformation_info.docx")
tryCatch(
  print(doc_transform, target = transform_doc_path),
  error = function(e) {
    message(sprintf("monoamine_transformation_info.docx write failed: %s", conditionMessage(e)))
  }
)
message(sprintf("Saved: %s", transform_doc_path))
n_missing <- sum(mono_transform_tbl == "not available")
message(sprintf("Transformation table: %d cells, %d 'not available' entries.",
                nrow(mono_transform_tbl), n_missing))


# =============================================================================
# ==== easy_scripts CSV export (endocrine) ====================================
# =============================================================================
# Writes a long-format CSV consumed by easy_scripts/endocrine/ mini-scripts.
# One row per sample × analyte (cortisol: NA area; monoamines: 4 areas).
# Pre-computes log_value and sqrt_value for the adaptive-transform models.
local({
  .easy_dir <- file.path(PROJECT_ROOT, "choice R pipeline/easy_scripts")
  dir.create(.easy_dir, showWarnings = FALSE, recursive = TRUE)

  # Cortisol long rows: one per sample, area = NA, analyte = "cort"
  .cort_rows <- cort_dat %>%
    dplyr::transmute(
      sample_id  = as.character(sample_id),
      tank       = as.character(tank),
      trial      = as.character(trial),
      treatment  = as.character(condition),
      sex        = as.character(sex),
      plate      = if ("plate" %in% names(.)) as.character(plate) else NA_character_,
      plate_date = if ("plate_date" %in% names(.)) as.character(plate_date) else NA_character_,
      area       = NA_character_,
      analyte    = "cort",
      value      = plasma_cortisol
    )

  # Monoamine long rows: one per sample × analyte × area
  .mono_rows <- mono_long %>%
    dplyr::transmute(
      sample_id  = as.character(sample_id),
      tank       = as.character(tank),
      trial      = as.character(trial),
      treatment  = as.character(treatment),
      sex        = as.character(sex),
      plate      = NA_character_,
      plate_date = NA_character_,
      area       = as.character(area),
      analyte    = analyte_key,
      value      = value
    )

  .endo <- dplyr::bind_rows(.cort_rows, .mono_rows) %>%
    dplyr::mutate(
      log_value  = ifelse(is.finite(value) & value > 0, log(value), NA_real_),
      sqrt_value = ifelse(is.finite(value) & value >= 0, sqrt(value), NA_real_)
    )

  .csv_path <- file.path(.easy_dir, "easy_scripts_endo_dataset.csv")
  readr::write_csv(.endo, .csv_path, na = "")
  message(sprintf("[easy_scripts] wrote %s: %d rows × %d cols",
                  basename(.csv_path), nrow(.endo), ncol(.endo)))
})
