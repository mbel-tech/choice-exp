# choice-exp

Analysis code for the study of behavioural, neural and physiological correlates
of swimming-exercise choice in fish.

This repository holds the code and the written record of the analytical
decisions behind it. It is a **curated subset** of a larger working tree: only
the scripts that back the manuscript's current claims are here. Superseded
versions, pipeline output, figures and the manuscript itself are deliberately
absent.

## Associated publication

Bellio, M., Alvarstein, H., Sneddon, L. U., Newberry, R. C., & Vindas, M. A.
(2026). Behavioural and neuroendocrine correlates of swimming exercise choice in
juvenile Atlantic salmon. *Hormones and Behavior*, 186, 106001.

> https://doi.org/10.1016/j.yhbeh.2026.106001

## Data availability

The dataset is deposited at Zenodo:

> https://doi.org/10.5281/zenodo.22162227

The DOI is reserved and the deposit is awaiting acceptance, so the link may not
resolve yet.

The data is not in this repository. What you can do without it splits in two:

### Reproducing the reported analyses — clone plus the DOI is enough

1. Download the record and unpack it to `zenodo_dataset/` in the repository
   root (or anywhere, and set `CHOICE_EXP_DATA_ROOT` to point there).
2. From the repository root:

```r
source("paths.R")
source("choice R pipeline/easy_scripts/by_timepoint/mean_nnd_cm_by_tp.R")
```

Each mini-script under `choice R pipeline/easy_scripts/` fits one outcome and
prints its model, ANOVA, post-hoc contrasts and figure. They read the deposited
CSVs through `easy_scripts/_data_access.R`, which translates the deposit's
published column names into the internal names the pipeline engines use.

To check the whole claim at once:

```r
source("choice R pipeline/scripts/07_cross_validate/verify_easy_scripts_from_deposit.R")
```

That refits each reported behavioural outcome from the deposit alone and checks
its treatment F and denominator df against the published value.

### Re-running the full pipeline — needs more than the DOI

The tracking stages need the raw idtracker.ai sessions (~14 GB), which are not
deposited: they are primary tracking output held in backup. Set
`CHOICE_EXP_DATA_DIR` to point at them. STEP1/STEP2 are the only stages that
read them; everything downstream works from derived indicators. The workbook
builders additionally expect
`Behavioural and neuroendocrine correlates datasets.xlsx` in `choice R pipeline/`
and `trial_summary_choice_exp.xlsx` in `choice R pipeline/data/`.

The scripts that built and verify the deposit are included:
`choice R pipeline/scripts/06_workbook/build_zenodo_deposit.R`,
`build_zenodo_docs.R`, and `verify_zenodo_deposit.R` — the last refits all
fourteen reported behavioural outcomes from the deposited CSV alone and checks
each treatment F/df against the published value.

## Layout

| Path | Contents |
|---|---|
| `choice R pipeline/scripts/` | The pipeline, in nine numbered stages |
| `choice R pipeline/easy_scripts/` | Single-outcome scripts, one model each |
| `choice exp for claude mono and cortisol/` | Cortisol and monoamine analysis |
| `handoff/` | How to read and extend this work — start with `HANDOFF.md` |
| `DECISIONS_LOG.md` | What changed, when, and why |
| `METHODS_CHANGES.md` | Drafted Materials & Methods text |
| `paths.R`, `config.R` | Root discovery; machine-specific locations |

## How paths work

No script hardcodes this project's drive letter. `paths.R` discovers
`PROJECT_ROOT` by walking up from the working directory to the
`.choice-exp-root` sentinel, or from the `CHOICE_EXP_ROOT` environment variable.
It fails loudly rather than guessing.

**Scripts do not source `paths.R` themselves** — they assume `PROJECT_ROOT` is
already defined, so source it first:

```r
source("paths.R")   # from the repository root
```

Machine-specific locations live in `config.R` and are overridable:

| Variable | Environment override | Default |
|---|---|---|
| `DATA_ROOT` | `CHOICE_EXP_DATA_ROOT` | `<PROJECT_ROOT>/zenodo_dataset` |
| `DATA_SOURCE` | `CHOICE_EXP_DATA_SOURCE` | `auto` (`deposit` / `internal`) |
| `PANDOC` | `CHOICE_EXP_PANDOC` | `D:/tools/pandoc-3.10.1/pandoc.exe` |
| `PIPELINE_DATA_DIR` | `CHOICE_EXP_DATA_DIR` | unset |

`DATA_ROOT` and `PIPELINE_DATA_DIR` differ by one word and mean opposite things.
`DATA_ROOT` is the ~200 kB Zenodo deposit that every *reported* analysis reads.
`PIPELINE_DATA_DIR` is ~14 GB of raw tracking sessions that only STEP1/STEP2
touch. If you are trying to reproduce a number from the paper, you want the
first one.

`DATA_SOURCE` chooses what the mini-scripts read when both are available.
`auto` prefers the deposit — so that the path a reader takes is the one the
author exercises too — and falls back to the pipeline-side
`easy_scripts_*.csv` when the deposit is absent. Every load announces its
source on one line.

## Decisions that govern interpretation

The full record is in `DECISIONS_LOG.md`. Three matter for reading the code:

**D1 — behavioural outcomes are reported at the trial × interval level only**
(N = 48 sessions). The trial-aggregated module is still present and still runs.
It is retained as a cross-check and is **not reported**.

**D6 — polarisation is excluded** from reported results, on identity-invariance
grounds. Code referencing polarisation nevertheless **remains in this
repository**: in the trial-level module, in `presentation_figures`, in a label
lookup, and in comments. It was not stripped — a previous removal left orphaned
references and killed a long run. Seeing polarisation code here does **not** mean
polarisation results are reported. They are not.

**Noradrenaline (NE) is excluded** throughout. The analyte × region grid is
24 cells, not 28.

## Terminology

The manuscript, the deposited workbook and this README all say **interval**. The
pipeline's internals say **`timepoint`**. They are the same thing.

There are two seams, and they are inverses of each other.

**Export** — `choice R pipeline/scripts/06_workbook/build_datasets_workbook.R`
renames internal names to published ones on the way out:

```r
names(d)[names(d) == "timepoint"] <- "interval"   # source header -> our name
```

**Import** — `choice R pipeline/easy_scripts/_data_access.R` reverses it on the
way back in, so the mini-scripts can read the deposit without changing a line of
their own vocabulary. The import seam re-derives the export seam's rename tables
at load and stops if they disagree, so the two cannot silently drift apart.

The rename is not only `interval`. `alr_flow` is internally `logit_flow`,
`mean_school_area_cm2` is `mean_hull_area_cm2`, and — the one that is easy to
miss — `crossings_per_session` is `zone_flux_per_session`, which is reported
outcome 5. `switches_per_session` exists under that name in *both* sources but
is a different, near-duplicate measure that backs no reported number.

Beware a false friend: `interval` in
`choice R pipeline/scripts/06_workbook/indicator_column_map.R` means something
else — it is a *level descriptor* (`"interval"` vs `"broadcast_mean"`) saying how
a column is aggregated, not the within-trial interval.

## Known rough edges

Inherited from the working tree and deliberately not silently "fixed", since
each would be a behavioural change beyond the scope of publishing the code:

- `choice R pipeline/scripts/02_utilities/diag_cells.R` carries a UTF-8 BOM and
  fails `parse()`. It fails identically in the working tree — pre-existing.
- `choice exp for claude mono and cortisol/claude output/Scripts/run_sex_cell.R`
  hardcodes `D:/choice exp for claude/claude output/Data`, a path from a
  different project root that does not exist here. Its sibling
  `endocrine_behaviour_association.R` resolves the same directory relative to its
  own output root; that is probably the intended form.
- `Figure 4 stuff script dataset image/choice exp occupancy through time UPDATED.R`
  reads an `.xlsx` from a OneDrive Desktop path specific to one machine.

## Provenance

Curated from the working tree at commit `3d79ce1ed4ea47255a4d37836ab826b4378a6772`. The absence of a file
here does not mean it never existed — it means it did not back a current claim.

Beyond replacing the hardcoded project root with `PROJECT_ROOT`, five categories
of change were made, all recorded here:

- a `source()` pointing at `choice R pipeline/activity_analysis_STATS_choice_exp.R`
  was repaired to `scripts/01_pipeline_analysis/`, the location the file moved to;
- three references to `analysis_b3.R` were redirected to the live
  `analysis_b3_REVISED.R`, which is the version this repository ships;
- `PANDOC` and `PIPELINE_DATA_DIR` were routed through `config.R` instead of
  being hardcoded, and the master pipeline's "not configured" message was updated
  to name `CHOICE_EXP_DATA_DIR` rather than a line to edit;
- a Quick-Start `setwd()` in `docs/pipeline_overview.qmd` pointing at a different
  project's root was replaced with `source("paths.R")`;
- the nineteen `easy_scripts` mini-scripts, plus `_template_mini_script.R`,
  `_generate_easy_scripts.R` and `_helpers.R`, were repointed at the Zenodo
  deposit through the new `easy_scripts/_data_access.R`. **The change is confined
  to how the data arrives**: every filter, factor level, model call, post-hoc and
  plot is untouched. `_helpers.R` additionally gained a graceful-degradation
  path, because the mini-scripts could not run from a clean clone at all before
  — the STATS engine loads pipeline output ~2,200 lines before the guard that
  was supposed to stop it, so a fresh checkout died before the helper functions
  were ever recovered.

Every script outside those five categories is byte-identical to the working tree
— verified by reversing the path substitution and comparing raw bytes, so line
endings and encoding are confirmed unchanged, not merely the text.

Two honest asterisks on the repointing:

- **`switches_per_session` is the one column whose value differs between the two
  sources**, by up to 0.014, because the deposit restores it to the integer count
  it always was (`build_datasets_workbook.R` rounds it on export). It backs no
  reported outcome — the data dictionary calls it a near-duplicate of
  `crossings_per_session` and explicitly not the reported measure. Every other
  mapped column agrees to ≤ 5e-12, which is CSV round-trip precision, not error.
- **`lr_medium_by_tp.R` does not reproduce the published `alr_medium` result**,
  and did not before this change either. It lets AICc choose the random effect
  freely and selects `(1 | tank)`; the manuscript reports `(1 | phys_trial_id)`,
  which `DECISIONS_LOG.md` D3 fixes *by design*. The mini-script does not
  implement that constraint, so it answers a slightly different question:
  F(1,39) = 8.93 rather than the reported F(1,14) = 5.28. Verified identical from
  both data sources, so it is not an artefact of reading the deposit. Recorded
  rather than quietly patched — changing a random-effect selection rule is an
  analysis decision, not a data-plumbing one.

## Verifying it yourself

| Script (under `choice R pipeline/scripts/07_cross_validate/`) | What it proves |
|---|---|
| `compare_deposit_to_internal.R` | The deposit and the pipeline CSVs are the same data, column by column, against pinned tolerances |
| `verify_easy_scripts_from_deposit.R` | Each reported outcome refits from the deposit alone to the published F and df |

The second is the one a stranger can run. It is also genuinely independent of
`scripts/06_workbook/verify_zenodo_deposit.R`, which checks the same numbers by
fitting hand-written models with hardcoded transforms; the mini-scripts get there
through AICc random-effect selection and an automatic transform search. Two
different routes, one deposited file, the same numbers.
