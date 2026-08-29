# 05_stats_validation/

Two-pronged statistical audit module.

## What it answers

1. **Are F-tests or chi-square tests actually being performed?**
   `validate_stat_types.R` scans every fitted indicator across both pipelines
   and labels the ANOVA mode actually used: F-KR / F-Satterthwaite / F-OLS /
   Wald-chisq (GLMM or beta-GLMM). Writes a per-indicator CSV +
   console summary.

2. **Are the post-hoc results actually time-specific?**
   `.run_ph_cld` in the main pipeline was patched to compute, in addition to
   the joint (all-pairs) Tukey, the **simple-effects** view for any
   interaction term:
   - `~ treatment | timepoint_f` → control vs treat **within each interval** (primary interest)
   - `~ timepoint_f | treatment` → interval pairs **within each treatment**
   For models using numeric `timepoint` (continuous RE formula), the same views are
   computed via `emmeans(..., at = list(timepoint = c(1, 2, 3)))`.
   Numeric stratifiers with ≤ 5 distinct observed values are treated as discrete.
   These are saved as `contrasts_<term>__simple_<focal>_within_<stratifier>.csv`
   in each indicator directory and are stored in
   `res$posthoc[[term]]$contrasts_simple`.

## Why simple-effects matter

The joint all-pairs Tukey over a 2 × 4 design produces **28 contrasts**.
Tukey adjustment for 28 contrasts is conservative enough that within-interval
differences visible in graphs often fail to reach significance even when the
omnibus interaction term is highly significant.

The simple-effects view (`~ treatment | timepoint_f`) yields **4 contrasts**
(one per interval), each adjusted by Tukey **within its stratum** — at most
one comparison per stratum when there are only two treatment levels. This
matches what the graphs visually display.

## Files

| File | Purpose |
|---|---|
| `validate_stat_types.R` | Scanner. Reads `STEP5_OUT/<label>/anova.csv` + `aicc_selection.csv`, plus live `fit_cort` / `mono_cell_models` for Pipeline B. Writes `STEP5_OUT/stat_type_validation.csv`. |
| `_run_validate.R`       | Convenience wrapper. Set `STEP5_OUT_TO_AUDIT` to a specific directory if not in-session. |
| `README.md`             | This file. |

## Usage

After both pipelines have been sourced in a single R session:

```r
source("D:/CHOICE R SCRIPTS/choice R pipeline/scripts/05_stats_validation/validate_stat_types.R")
```

Or stand-alone against any prior STEP5 directory:

```r
STEP5_OUT_TO_AUDIT <- "D:/CHOICE R SCRIPTS/choice R pipeline/STEP5_stats/STEP5_stats_YYYYMMDD_HHMMSS"
source("D:/CHOICE R SCRIPTS/choice R pipeline/scripts/05_stats_validation/validate_stat_types.R")
```

## CSV columns

| Column | Meaning |
|---|---|
| `label` | Indicator name (= `STEP5_OUT/<label>` subdir or `monoamine_cell__<cell_id>`) |
| `pipeline` | `"A"` (behaviour) or `"B"` (endocrine) |
| `aicc_fit_method` | `REML` (lmer LMM), `MLE` (glmmTMB beta), `OLS` (lm fallback), `NA` if not applicable |
| `re_winner` | Selected random-effects formula or `"1"` for OLS |
| `family_used` | `gaussian`, `beta`, `poisson`, `negbin`, etc. |
| `anova_mode` | The actual test used (the key output): `F-KR`, `F-OLS`, `Wald-chisq`, `Wald-chisq (beta-GLMM)`, `Wald-chisq (GLMM)`, `SKIPPED` |
| `has_df_denom` | `TRUE` iff denominator df was computed (i.e. an F-test was run, not a Wald-chisq) |
| `focal_term` | Highest-order term in the ANOVA (the interaction in 2-factor models) |
| `focal_statistic`, `focal_df1`, `focal_df2`, `focal_p` | Numeric stats for that term |
| `focal_stat_label` | Pretty string: `F(df1, df2) = X.XXX, p = X.XXXX` or `χ²(df) = ...` |

## Expected mode distribution (after a clean run)

- ~34 indicators **F-KR**: every `lmer` Gaussian (LMM) — collective movement,
  zone preference logit, Jacobs atanh, NND, polarisation, IID, hull, centroid
  speed, flux, etc.
- ~13 indicators **Wald-chisq (beta-GLMM)**: every `glmmTMB::beta_family`
  fit — Jacobs β-GLMM zones, zone main β-GLMM cells, zone sub CLR.
- 2 indicators **F-OLS**: `zone_main_desc` and `zone_sec_desc` (independent
  ANOVA paths).
- 0–2 **SKIPPED** when a response is near-constant (e.g. `prop_active` is
  often invariant in this dataset, hence the `[skipped] Response near-constant`
  note rather than a silent failure).
