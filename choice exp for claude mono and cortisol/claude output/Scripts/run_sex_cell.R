# Standalone: per-cell (analyte x area) sex effect analysis
# Loads monoamine_long.csv, runs fit_one_mono_sex_cell for all 28 combinations,
# exports sex_cell_monoamines.csv

suppressPackageStartupMessages({
  library(dplyr)
  library(lme4)
  library(lmerTest)
  library(car)
  library(MuMIn)
  library(purrr)
})

DATA_DIR <- "D:/choice exp for claude/claude output/Data"

# ---- helpers (copied from analysis_b3.R) ----------------------------------

set_sum_contrasts <- function(df, vars) {
  for (v in vars) if (v %in% names(df) && is.factor(df[[v]])) {
    contrasts(df[[v]]) <- contr.sum(nlevels(df[[v]]))
  }
  df
}

select_re_aicc <- function(fixed_formula, re_list, data) {
  fits <- list(); aics <- c(); singular <- c()
  for (lbl in names(re_list)) {
    re_term <- re_list[[lbl]]
    full_form <- update(fixed_formula, paste(". ~ . +", re_term))
    fit <- tryCatch(
      suppressMessages(suppressWarnings(
        lmerTest::lmer(full_form, data = data, REML = FALSE)
      )),
      error = function(e) NULL
    )
    if (is.null(fit)) { aics[lbl] <- NA_real_; singular[lbl] <- NA; next }
    aics[lbl]     <- tryCatch(MuMIn::AICc(fit), error = function(e) AIC(fit))
    singular[lbl] <- lme4::isSingular(fit, tol = 1e-4)
    fits[[lbl]]   <- fit
  }
  ok <- !is.na(aics) & !singular
  if (!any(ok)) {
    return(list(table = data.frame(re = names(re_list), AICc = aics, singular = singular,
                                   selected = FALSE), winner = NULL, label = NA_character_))
  }
  win_lbl  <- names(aics)[ok][which.min(aics[ok])]
  win_re   <- re_list[[win_lbl]]
  win_form <- update(fixed_formula, paste(". ~ . +", win_re))
  win_reml <- suppressMessages(suppressWarnings(
    lmerTest::lmer(win_form, data = data, REML = TRUE)
  ))
  list(winner = win_reml, label = win_lbl)
}

anova_kr <- function(model, type = 3) {
  if (inherits(model, "lmerMod") || inherits(model, "lmerModLmerTest")) {
    tryCatch(
      car::Anova(model, type = type, test = "F", ddf = "Kenward-Roger"),
      error = function(e) car::Anova(model, type = type, test = "F", ddf = "Satterthwaite")
    )
  } else {
    car::Anova(model, type = type)
  }
}

# ---- load data -------------------------------------------------------------

mono_long <- read.csv(file.path(DATA_DIR, "monoamine_long.csv"), stringsAsFactors = FALSE)

# Make factors
for (v in c("treatment", "area", "sex", "date", "tank")) {
  if (v %in% names(mono_long)) mono_long[[v]] <- as.factor(mono_long[[v]])
}

analyte_order <- c("5-HT", "5-HIAA", "DA", "DOPAC", "5-HIAA/5-HT", "DOPAC/DA", "NE")
expected_areas <- c("DM", "POA", "VV", "VD")

# analyte_key -> label map (from analyte column in data)
akey_map <- mono_long %>%
  distinct(analyte_key, analyte) %>%
  tibble::deframe()   # analyte_key -> analyte

analyte_cols <- names(akey_map)

# ---- per-cell sex function -------------------------------------------------

fit_one_mono_sex_cell <- function(akey, ar) {
  ana_label <- akey_map[akey]
  d <- mono_long %>%
    filter(analyte_key == akey, area == ar,
           is.finite(value), value > 0,
           !is.na(treatment), sex %in% c("F", "M")) %>%
    droplevels()
  if (nrow(d) < 8 || length(unique(as.character(d$sex))) < 2) {
    message(sprintf("  SKIP [%s x %s]: n=%d, sex levels=%d",
                    ana_label, ar, nrow(d), length(unique(as.character(d$sex)))))
    return(NULL)
  }
  d$log_v <- log(d$value)
  d <- set_sum_contrasts(d, c("treatment", "sex"))

  fixed_form <- log_v ~ treatment + sex
  re_cand    <- list("(1|date)" = "(1|date)", "(1|tank)" = "(1|tank)")
  sel        <- select_re_aicc(fixed_form, re_cand, data = d)

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

  # Extract sex term p-value via anova_kr (works for both lmer and lm)
  kr <- tryCatch(as.data.frame(anova_kr(fit, type = 3)), error = function(e) NULL)
  if (is.null(kr)) return(NULL)

  # car::Anova column names differ by model type
  p_col <- grep("^Pr", names(kr), value = TRUE)[1]
  f_col <- grep("^F", names(kr), value = TRUE)[1]
  df1_col <- grep("Df$|Num Df", names(kr), value = TRUE)[1]
  df2_col <- grep("Den Df|Res.Df", names(kr), value = TRUE)
  df2_col <- if (length(df2_col) > 0) df2_col[1] else NA

  sex_row <- kr[rownames(kr) == "sex", , drop = FALSE]
  if (nrow(sex_row) == 0) {
    # try partial match
    sex_idx <- grep("^sex$", rownames(kr))
    if (length(sex_idx) == 0) return(NULL)
    sex_row <- kr[sex_idx, , drop = FALSE]
  }

  f_val  <- if (!is.na(f_col))  sex_row[[f_col]]  else NA_real_
  df1    <- if (!is.na(df1_col)) sex_row[[df1_col]] else NA_real_
  df2    <- if (!is.na(df2_col) && !is.na(df2_col)) sex_row[[df2_col]] else NA_real_
  p_val  <- sex_row[[p_col]]

  # n per sex group
  sex_tab <- table(d$sex)

  data.frame(
    area          = ar,
    analyte       = ana_label,
    n_total       = nrow(d),
    n_F           = as.integer(sex_tab["F"]),
    n_M           = as.integer(sex_tab["M"]),
    re_spec       = re_label,
    model         = model_formula,
    F_stat        = f_val,
    df1           = df1,
    df2           = df2,
    p_sex_raw     = p_val,
    stringsAsFactors = FALSE
  )
}

# ---- run all 28 combinations -----------------------------------------------

combos <- expand.grid(akey = analyte_cols, area = expected_areas,
                      stringsAsFactors = FALSE)

message(sprintf("Running %d analyte x area sex models...", nrow(combos)))

results_list <- purrr::compact(purrr::map2(
  combos$akey, combos$area,
  function(ak, ar) tryCatch(
    fit_one_mono_sex_cell(ak, ar),
    error = function(e) {
      message(sprintf("  ERROR [%s x %s]: %s", akey_map[ak], ar, conditionMessage(e)))
      NULL
    })
))

sex_cell_df <- dplyr::bind_rows(results_list)

# BH correction per area
sex_cell_df <- sex_cell_df %>%
  group_by(area) %>%
  mutate(p_sex_BH = p.adjust(p_sex_raw, method = "BH")) %>%
  ungroup() %>%
  arrange(area, analyte)

# ---- report ----------------------------------------------------------------

cat("\n=== Per-cell sex effect results ===\n")
print(sex_cell_df[, c("area","analyte","n_total","re_spec","F_stat","df1","df2","p_sex_raw","p_sex_BH")],
      row.names = FALSE)

min_idx <- which.min(sex_cell_df$p_sex_raw)
cat(sprintf("\nMinimum raw p: %.6f  =>  %s x %s\n",
            sex_cell_df$p_sex_raw[min_idx],
            sex_cell_df$analyte[min_idx],
            sex_cell_df$area[min_idx]))
cat(sprintf("BH-adjusted p at minimum: %.6f\n", sex_cell_df$p_sex_BH[min_idx]))

# ---- export ----------------------------------------------------------------

out_path <- file.path(DATA_DIR, "sex_cell_monoamines.csv")
write.csv(sex_cell_df, out_path, row.names = FALSE)
cat(sprintf("\nSaved: %s\n", out_path))
