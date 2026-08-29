# choice-exp

Analysis code for the study of behavioural, neural and physiological correlates
of swimming-exercise choice in fish.

This repository holds the code and the written record of the analytical
decisions behind it. It is a **curated subset** of a larger working tree: only
the scripts that back the manuscript's current claims are here. Superseded
versions, pipeline output, figures and the manuscript itself are deliberately
absent.

## Data availability

The dataset is deposited at Zenodo:

> https://doi.org/10.5281/zenodo.22162227

The DOI is reserved and the deposit is awaiting acceptance, so the link may not
resolve yet.

**This repository does not run from a clean clone.** The input data is not
included — it lives at the DOI above. To run the pipeline:

1. Download the dataset record from Zenodo.
2. Place `Behavioural and neuroendocrine correlates datasets.xlsx` in
   `choice R pipeline/`.
3. Place `trial_summary_choice_exp.xlsx` in `choice R pipeline/data/`.
4. From the repository root, `source("paths.R")`, then run
   `choice R pipeline/scripts/00_main/00_master_pipeline_choice_exp.R`.

Raw idtracker.ai session data (~14 GB) is not deposited and is not needed for
the statistical pipeline; only the STEP1/STEP2 tracking stages read it. Set
`CHOICE_EXP_DATA_DIR` to point at it if you need those stages.

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
| `PANDOC` | `CHOICE_EXP_PANDOC` | `D:/tools/pandoc-3.10.1/pandoc.exe` |
| `PIPELINE_DATA_DIR` | `CHOICE_EXP_DATA_DIR` | unset |

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

`choice R pipeline/scripts/06_workbook/build_datasets_workbook.R` is the seam
where the translation happens, on export:

```r
names(d)[names(d) == "timepoint"] <- "interval"   # source header -> our name
```

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

Beyond replacing the hardcoded project root with `PROJECT_ROOT`, four changes
were made, all recorded here:

- a `source()` pointing at `choice R pipeline/activity_analysis_STATS_choice_exp.R`
  was repaired to `scripts/01_pipeline_analysis/`, the location the file moved to;
- three references to `analysis_b3.R` were redirected to the live
  `analysis_b3_REVISED.R`, which is the version this repository ships;
- `PANDOC` and `PIPELINE_DATA_DIR` were routed through `config.R` instead of
  being hardcoded, and the master pipeline's "not configured" message was updated
  to name `CHOICE_EXP_DATA_DIR` rather than a line to edit;
- a Quick-Start `setwd()` in `docs/pipeline_overview.qmd` pointing at a different
  project's root was replaced with `source("paths.R")`.

Every other script is byte-identical to the working tree — verified by reversing
the path substitution and comparing raw bytes, so line endings and encoding are
confirmed unchanged, not merely the text.
