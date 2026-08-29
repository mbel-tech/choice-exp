# `easy_scripts/standalone/` — handwritten self-contained reproductions

The scripts in this folder reproduce two Pipeline-A models **without**
sourcing any project file. Only CRAN packages are loaded. Every choice
the pipeline makes silently (response transform, RE structure, post-hoc
adjustment, display-only letter remap) is spelled out as numbered code
with a comment explaining *why*.

## Why these two indicators?

Inter-individual distance (IID) and polarisation, both by interval,
produced "dataset / figure / post-hocs don't collimate" reports. The
auto-generated counterparts in `easy_scripts/by_timepoint/` faithfully
reproduce the pipeline by sourcing `_helpers.R`, but they're opaque —
the reader can't see why letters in the printed CLD differ from the
letters drawn on Figure_10, or why the figure dots don't fall on the
model's predicted line.

## What sits between raw data and Figure_10's letters

Three independent transformations stack up. Once you see them laid out
the "doesn't collimate" feeling goes away.

| layer | what the pipeline does | visible in standalone script? |
|---|---|---|
| 1. response transform | `sqrt(y)` for both indicators (auto-picked from `{identity, log1p, sqrt}` by Shapiro-Wilk + Levene) | Step 4 — prints the SW/Levene numbers; comment explains *why* sqrt |
| 2. continuous timepoint | `treatment * timepoint` (numeric) beats `treatment * timepoint_f` (3-level factor) on AICc, so the chosen model has `df_num = 1` for the interaction — the interaction is the *slope difference between treatments*, not the {tp2 vs tp1, tp3 vs tp1} pair you'd expect from a factor | Step 6 + Step 7 — formula and ANOVA table |
| 3. CLD display remap | `.mg_remap_cld()` walks cells in precedence order (control before exercise choice, intervals 1→2→3) and relabels letters in first-seen order. The saved CSV `cld_treatmentxtimepoint.csv` is the *raw* output; what you see on Figure_10 is *remapped*. For IID this remap is identity. For polarisation it changes most letters. | Step 9 — prints raw letters AND remapped letters side by side |

## The "where does each number come from" map

| artefact | what it shows |
|---|---|
| `easy_scripts_dataset.csv` (column `mean_iid_cm` / `mean_polarisation`) | raw per-trial-per-interval values |
| `STEP5_stats/STEP5_stats_<ts>/<indicator>_timepoint/anova.csv` | Type-III Kenward-Roger F-stats on the `sqrt(y) ~ treatment * timepoint + (1|RE)` model |
| `STEP5_stats/STEP5_stats_<ts>/<indicator>_timepoint/cld_treatmentxtimepoint.csv` | RAW emmeans/multcomp letters (before display remap) |
| `STEP5_stats/STEP5_stats_<ts>/<indicator>_timepoint/aicc_selection.csv` | AICc across all 12 RE × {factor, continuous} time candidates — the row with `selected = TRUE` is what got fitted |
| `STEP5_stats/STEP5_stats_<ts>/<indicator>_timepoint/normality_check.csv` | SW + Levene p-values that picked the transform |
| Figure_10 panel C (IID) / panel B (polarisation) | OBSERVED cell means (not emmeans) with the REMAPPED CLD letters; panel B additionally collapses interval 2 and forces interval 3 to split for cosmetics |

## Files

```
easy_scripts/standalone/
├─ README.md                                          ← this file
├─ by_timepoint/
│  ├─ mean_iid_cm_by_tp.R                             ← standalone IID
│  └─ mean_polarisation_by_tp.R                       ← standalone polarisation
└─ outputs/
   ├─ mean_iid_cm_by_tp.png                           ← created on first run
   └─ mean_polarisation_by_tp.png                     ← created on first run
```

## Running

In RStudio: open the script and source it.
From a terminal:

```bash
Rscript "easy_scripts/standalone/by_timepoint/mean_iid_cm_by_tp.R"
Rscript "easy_scripts/standalone/by_timepoint/mean_polarisation_by_tp.R"
```

Each script is around 130 lines. Run time is a few seconds.

## What they print, in order

1. Path to the CSV being loaded, and N rows / N per cell after filtering.
2. Why `sqrt` was the chosen transform (SW + Levene on raw vs `sqrt(y)`).
3. Why this RE wins (top-3 AICc shown inline; comment cites the full
   12-row pipeline table on disk).
4. Type-III Kenward-Roger ANOVA — matches `anova.csv` to ≥4 sig figs.
5. Tukey CLD on `emmeans` at integer intervals — RAW letters that match
   the saved `cld_treatmentxtimepoint.csv` exactly.
6. The same CLD table, *with an extra `.display` column* showing the
   Figure_10 remapped letters.
7. Observed cell means + SE — the values plotted as dots on Figure_10.
8. A minimal `ggplot` saved to `outputs/` showing dots with both letter
   rows (bold = raw, parenthesised = remapped).

## Important: continuous vs factor timepoint

The pipeline always **tries both** time encodings (factor and continuous)
and picks the lower-AICc fit. For these two indicators continuous wins.
That has a downstream consequence the user should know about:

- ANOVA: `treatment:timepoint` has `df_num = 1` (it's a slope-difference
  test), so a "no interval-level interaction" finding doesn't rule out
  *categorical* interval-specific effects you'd catch with `timepoint_f`.
- Predicted means: emmeans evaluates a *line*, so the predicted mean at
  interval 2 is exactly the midpoint of intervals 1 and 3 by
  construction. Observed cell means won't generally lie on that line.
- CLD: Tukey is run over only 6 model-implied means, so the letter
  pattern reflects the linear-fit gradient, not the raw cell pattern.

If you want the factor-time alternative, swap `timepoint` for
`timepoint_f` (after `d$timepoint_f <- factor(d$timepoint)`) in steps 6
and 8 and re-run. The pipeline's `aicc_selection.csv` shows what that
costs in AICc.

## Pass criteria (numbers to expect)

**IID by timepoint**
- ANOVA: treatment p = 0.740; timepoint p = 0.069;
  treatment:timepoint F(1, 30) = 6.096, p = 0.019.
- Raw CLD: ec/3 = a; ec/2 = b; ec/1 = c; ctrl/1 = ctrl/2 = ctrl/3 = abc.
- Display remap: identity (raw letters already in precedence order).

**Polarisation by timepoint**
- ANOVA: treatment F(1, 37) = 4.343, p = 0.044; timepoint p = 0.160;
  treatment:timepoint p = 0.337.
- Raw CLD: ec/3 = ab; ec/2 = a; ec/1 = abc; ctrl/3 = abc; ctrl/2 = bc;
  ctrl/1 = c.
- Display remap: ctrl/1 = a; ctrl/2 = ab; ctrl/3 = abc; ec/1 = abc;
  ec/2 = c; ec/3 = bc.

## Pipeline source pointers

- Model wrapper: `scripts/01_pipeline_analysis/activity_analysis_STATS_choice_exp.R`
  - `run_lmm_analysis()` — lines 1057–1342
  - `res_iid_tp <- ...` — line 2967
  - `res_pol_tp <- ...` — line 2942
- Figure assembly: same file, `.mg_make_line_with_tukey()` around L6079,
  `.fig9_C` (IID) at L6265, `.fig9_B` (polarisation) at L6253.
- CLD display remap: `.mg_remap_cld()` at lines 6056–6076.
