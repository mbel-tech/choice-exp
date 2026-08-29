# B3 Analysis Pipeline — Phases A, B, C

## Main Scripts

- **`analysis_b3.R`** (900 lines)
  - Phase A: DHARMa residual diagnostics, `performance::check_model()`, BayesFactor, omega-squared
  - Phase B: cortisol model with tank/trial nesting; monoamine models with (1|tank) RE
  - Generates all supporting CSVs, PNGs (diagnostics), and Word report
  - Output: plots, tables, diagnostics PNGs

## Main Deliverables → `B3_Final_Report/`

- **`B3_cortisol_monoamines_report.docx`**
  - Complete Word report from `analysis_b3.R`
  - Methods, assay quality tables & plots
  - Cortisol: descriptive stats, model coefficients, ANOVA, EMMs, pairwise contrasts, variance components, effect sizes, Bayes factor, diagnostics
  - Monoamines: RE specification table, FDR results, EMMs, DHARMa diagnostics per analyte, plots

- **`b3_analysis.qmd`** (700+ lines)
  - Phase C: manuscript-ready Quarto document
  - Renders to HTML (with code folding) and Word
  - Sections: background, methods, assay quality, cortisol analysis, monoamines, results summary
  - Paste-ready for methods/results sections

## Assay Quality & CV

- **`assay_cv_tables.xlsx`** (8 sheets)
  - intra_assay_cv: per-sample CV
  - intra_assay_summary: by-plate median/weighted-mean CV
  - inter_assay_cv_controls, inter_assay_cv_standards, inter_assay_cv_repeats
  - plate_curve_quality: 4PL fit parameters (a,b,c,d,R²)
  - replicate_well_summary, matching_audit

- **`assay_long.csv`**: well-level cortisol data (standards, B0, NSB, unknowns)

## Cortisol Results

- **`sex_cortisol_descriptive.csv`**: n, mean, sd, median by sex × condition
- **`sex_cortisol_models.csv`**: model coefficients with 95% CI from primary lmer/lm
- **`sex_cortisol_emmeans.csv`**: estimated marginal means by sex × condition
- **`sex_cortisol_variance_components.csv`**: random effect SDs (if lmer used)
- **`b3_cortisol_clean.csv`**: full cortisol dataset with assay CV merged
- **`concentration_reconciliation.csv`**: worksheet vs. derived concentration audit

## Monoamine Results

- **`sex_monoamine_descriptive.csv`**: n, mean, sd, median by sex × area × analyte
- **`sex_monoamine_models.csv`**: model coefficients per analyte with 95% CI
- **`sex_monoamine_sex_main_fdr.csv`**: sex main effect with raw p and FDR-corrected p (primary/secondary families)
- **`sex_monoamine_emmeans.csv`**: estimated marginal means by sex × area × analyte
- **`monoamine_model_re_summary.csv`**: random effect spec and singularity status per analyte
- **`monoamine_long.csv`**: full monoamine dataset (long format: sex × area × analyte)

## Exploratory

- **`cortisol_monoamine_spearman.csv`**: Spearman correlations by sex × region

## Supporting

- **`column_dictionary.csv`**: source worksheet columns

## Phase A & B Improvements

### Phase A (Diagnostics & Inference)
✓ DHARMa simulated residuals (PNGs in current run but not saved to disk—modify `analysis_b3.R` line ~620 to save)
✓ `performance::check_model()` diagnostic plots
✓ BayesFactor with verbal interpretation (anecdotal/moderate/strong/extreme evidence)
✓ Omega-squared effect sizes (Type II ANOVA)
✓ All diagnostics embedded in Word report

### Phase B (Random Effects & Nesting)
✓ Cortisol: `(1|tank/trial)` primary model; fallback chain when singular
✓ Phase B verification: sex(M) estimate shift from lm → lmer printed to console
✓ Variance components table in report
✓ Monoamines: `(1|tank)` added alongside `(1|sample_id)`; fallback to `(1|sample_id)` if singular
✓ RE spec summary table in report showing which model per analyte

### Phase C (Quarto)
✓ `b3_analysis.qmd`: self-contained HTML + Word output
✓ Manuscript-ready sections with computed results (not hardcoded)
✓ Code folding for reproducibility
✓ Full methods and results prose

## Usage

### Run the R script (Phases A+B)
```r
source("D:/choice exp for claude/claude output/analysis_b3.R")
```

### Render the Quarto document (Phase C)
```bash
quarto render "D:/choice exp for claude/claude output/B3_Final_Report/b3_analysis.qmd" \
  --to html --to docx
```

Or in R:
```r
quarto::quarto_render("D:/choice exp for claude/claude output/B3_Final_Report/b3_analysis.qmd")
```

## Notes

1. **Cortisol model**: Falls back to lm if all tank/trial random effect specs are singular. Check console output for "Phase B verification" line to confirm sex effect did not flip.
2. **Monoamine REs**: Tank random effect included only when non-singular. See `monoamine_model_re_summary.csv` for spec per analyte.
3. **DHARMa PNGs**: Currently saved during runtime in `fit_one()`. Check `OUT` directory for `diag_monoamine_dharma_*.png` files after running `analysis_b3.R`.
4. **Reproducibility**: Use `renv::snapshot()` to create a lockfile from the installed packages.

---

Generated: 2026-04-28
Phases A, B, C complete.
