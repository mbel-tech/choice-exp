# Pipeline audit plan — systematic search for defects of the 2026-08-18 class

## Context

Three defects found on 2026-08-18 were not independent bugs. They share a shape:
**a value carried further along the chain than its name records, and nothing checks
the chain end-to-end.**

- Sub-zone proportions were area-normalised in STEP2 §6b, then normalised *again*
  in the log-ratio blocks, because `trial_occupancy_long$prop_time` does not record
  that it is already `prop_time_ac_*`.
- The Haldane floor `.eps_lr = 1e-4` was applied *after* that division, so it meant
  "0.42% of session time" for calm instead of the documented "~0.01%" — and clamped
  8 of 48 sessions, all in one treatment arm.
- The same four-line block sat in three files and was wrong in all three.

Two of these changed published statistics. None was caught by anything.

Scope: 52 R scripts, 39,381 lines (excluding `easy_scripts/` and backups).

## Why the existing validation did not catch it

The repo already has real validation infrastructure, and it should be reused — but
its axis is wrong for this defect class:

| existing | what it proves |
|---|---|
| `07_cross_validate/compare_data_identity.R` | workbook cells == easy-scripts CSV cells |
| `07_cross_validate/refit_equivalence.R` | same ANOVA from workbook-derived vs original data |
| `05_stats_validation/validate_stat_types.R` | which test type each indicator actually used |

All three check **consistency between representations of the same derived value**. A
value that is consistently wrong passes every one of them. `lr_high` was wrong in the
CSV, wrong in the workbook, and wrong on refit — perfect agreement, all layers green.

The missing axis is **correctness against the documented definition**, and the only
instrument that has it is independent recomputation from raw data. That is what
actually caught this one.

## Defect taxonomy

| class | observed instance | detectable by |
|---|---|---|
| **A** transform applied twice along a chain | area normalisation | Pass 3 |
| **B** guard/offset applied at the wrong point | `.eps_lr` after division | Pass 3 |
| **C** stale pinned artefact | caption caches, run pins | Pass 1 |
| **D** duplicated block drifts | `lr_*` block × 3 files | Pass 2 |
| **E** documentation contradicts code | Table 1 crossings row | Pass 4 |

## Pass 1 — provenance and freshness sweep (mechanical, fast)

Already partly evidenced. **13 live hardcoded run pins** remain, and most are worse
than stale — they use the *pre-reorganisation* path (`choice R pipeline/STEP5_stats/`)
rather than the real `output/STEP5_stats/`, so they point at directories that no
longer exist:

- `00_main/export_indicators_excel.R`, `02_utilities/diag_cells.R`,
  `04_reporting/build_reviewer_response_report.R`,
  `00_shared/effect_size_and_sample_size.R`,
  `05_stats_validation/_run_validate.R` + `validate_stat_types.R`
  → all pinned to `STEP5_stats_20260511_175522` (11 May, three months old)
- `00_shared/sensitivity_forced_trial_re.R`, `_v2.R`,
  `sensitivity_part2_loo_icc_equiv.R`, `sensitivity_multiverse.R`
  → pinned to `STEP5_stats_20260808_123821` + `STEP2b_output_20260807_214103`

Two of those matter directly: the **validation script audits a three-month-old run**,
and the **four sensitivity analyses that back methodological claims do not reflect
even the 17 August re-run**, let alone today's fix.

Steps:
1. Enumerate every pinned run and every one of the 126 absolute paths; for each record
   path-exists, run-age, and whether the script is live or superseded.
2. Replace live pins with the auto-detect block now standardised across
   `rebuild_figure_caption_table{,_interval,_cells}.R` and `run_manu_graphs_only.R`.
3. Add a freshness assertion to derived caches: each records the run it was built from,
   and its consumer stops if that run is not the current one. The pattern already
   exists — `bout_structure_analysis_choice_exp.R` §2 asserts its STEP1 input postdates
   the 2026-08-07 zone fix. Generalise that.

Output: `audit_pins.csv` — file, line, pinned run, exists, age_days, live, verdict.

## Pass 2 — duplicated-block drift (mechanical, fast)

The `lr_*` block existed verbatim in three files. Find the rest.

1. Normalise whitespace/comments and hash every sliding window of ≥8 statements across
   the 52 scripts; report clusters appearing in ≥2 files.
2. For each cluster, diff the copies. **Divergence is the signal** — identical copies are
   a maintenance smell, divergent copies are a live bug (one was fixed, the other wasn't).
3. Rank by whether the block computes a reported outcome.

Output: `audit_duplicates.csv` — block hash, files, line ranges, identical/divergent.

## Pass 3 — outcome provenance and independent recomputation (the core)

This is the only pass that catches A and B, and it is where the effort belongs.

For each of the 24 analyte × region cells' behavioural counterparts — concretely, every
outcome in manuscript Table 1 plus those in Table S10 — do two things:

**3a. Trace the chain.** Write the transform chain as one line, from the raw column in
`master_fish_by_frame.csv` to the number in the manuscript. Example, as it should read
after the fix:

```
prop_time_ac_high : n_in_high/n_total → weighted.mean(w=dt) → /10 → /Σ → log(·/calm) → LMM
                    STEP1              STEP2 §6            STEP2 §6b      STATS §2
```

Then assert, per chain: every transform appears **exactly once**; every guard, floor or
offset is applied at the point its documentation claims; every unit conversion happens
once and in the stated direction.

**3b. Recompute independently.** Write a standalone recomputation script that shares no
code with the pipeline, reads only `master_fish_by_frame.csv`, and derives each outcome
from its Table 1 definition. Compare to the pipeline's value per session.

Tolerance: exact for counts and durations, relative 1e-6 for floats. Any disagreement is
a finding — in either direction, since the recomputation can also be wrong, and the
disagreement is what forces the definition to be pinned down.

Order of work: the ~10 outcomes that appear in Results §3.1 first, since those are the
publishable claims. `zone_flux_per_session`, `switches_per_session`, the four collective
metrics and the bout metrics all have exact arithmetic definitions and are cheap to
recompute; the ALR family is already done.

Output: `audit_provenance.csv` (one row per transform step) and `audit_recompute.csv`
(one row per outcome × session, with `abs_diff`, `rel_diff`, `status`).

## Pass 4 — documentation reconciliation

Five places claim to define these outcomes, and they disagree:

- manuscript Table 1 (`Computation from the tracking data`)
- `00_shared/export_outcome_catalogue.R`
- `06_workbook/legend_descriptions.R`
- `06_workbook/indicator_column_map.R`
- `04_reporting/manifest_indicators.R`

Reconcile all five against the Pass 3 chain, one row per outcome × source. Already known
to be wrong: the catalogue calls `flux_timepoint` *"Count of zone crossings per session;
GLMM (Poisson/NB)"* — it is neither a count nor a Poisson GLMM, it is a continuous flux
on a Gaussian LMM with log1p. That wording is what the Table 1 row inherited.

Output: `audit_docs.csv` — outcome, source, claimed definition, actual, verdict.

## Rules for the audit

- **Find, don't fix.** Passes produce a ranked findings list. Fixes are a separate,
  approved pass — several findings will be judgment calls, not defects (the ALR(medium)
  random-effect question from today is exactly this shape).
- **Every finding carries its evidence**: file:line, the observed vs expected number, and
  the session or row where they diverge. No finding on reasoning alone.
- **State what a pass cannot see.** Pass 3 cannot validate STEP1's tracking or zone
  geometry — it takes `master_fish_by_frame.csv` as ground truth. That is a real limit
  and belongs in the report, not in a footnote.

## Verification of the audit itself

Before trusting a green result: take a scratch copy of the pipeline, re-introduce the
three defects already fixed today (double division, floor-after-division, stale caption
cache), and confirm each pass flags its own class. An audit that misses known-planted
bugs is worse than none — it manufactures confidence.

## Cost

Passes 1, 2 and 4 are mechanical and can be done in one working session. Pass 3 is the
bulk: roughly a day for the ~10 Results outcomes, more for the full Table S10 set.
Recommend Pass 1 + Pass 3-on-Results-outcomes first — that combination covers the two
classes that have already changed published numbers.
