# Choice R Pipeline — Folder Structure

## Organization

```
choice R pipeline/
├── scripts/                          # All R scripts
│   ├── 00_main/                      # Entry points
│   │   ├── 00_master_pipeline_choice_exp.R      # Full pipeline runner
│   │   └── run_manu_graphs_only.R               # Figure regeneration only
│   │
│   ├── 01_pipeline_analysis/         # Core analysis scripts (run in order)
│   │   ├── activity_analysis_STEP2_choice_exp.R
│   │   ├── activity_analysis_STEP3_choice_exp.R
│   │   ├── activity_analysis_STATS_choice_exp.R
│   │   ├── activity_analysis_STATS_trial_choice_exp.R
│   │   └── activity_analysis_GRAPHS_choice_exp.R
│   │
│   ├── 02_utilities/                 # Diagnostic & helper scripts
│   │   ├── idtracker_STEP1_choice_exp.R
│   │   ├── group_dynamics_STEP2b_choice_exp.R
│   │   ├── zone_switch_diagnostics_choice_exp.R
│   │   ├── diag_cells.R
│   │   └── aestetic stuff for new manu.R
│   │
│   └── 03_interactive_tools/         # Interactive Shiny/tuning apps
│       ├── style_tuner_app_choice_exp.R
│       └── skinny_graphs_tuner_choice_exp.R
│
├── output/                           # Pipeline output folders
│   ├── STEP1_output/                 # idTracker processing results
│   ├── STEP2_output/                 # Video-frame analysis
│   ├── STEP2b_output/                # Group dynamics
│   ├── STEP4_graphs/                 # Exploratory graphs
│   ├── STEP5_stats/                  # Statistical analysis & figures (main results)
│   ├── checked_sessions_choice_exp/  # Session validation data
│   └── zone_switch_diagnostics/      # Zone-switch diagnostic output
│
├── docs/                             # Documentation
│   ├── pipeline_overview.qmd         # Source documentation (Quarto)
│   ├── pipeline_overview.html        # Rendered pipeline overview
│   ├── pipeline_user_guide.docx      # User guide
│   └── Preference experiment results idtrackerai.docx
│
├── config/                           # Configuration & reference files
│   ├── zone_reference_choice_exp_FD.csv
│   └── zone_reference_choice_exp_FT.csv
│
├── reports/                          # Final reports
│   └── analysis_report_choice_exp.docx
│
├── data/                             # Input data files
│   └── trial_summary_choice_exp.xlsx
│
├── logs/                             # Script execution logs
│   ├── run_manu_graphs_only.log
│   └── run_manu_graphs_only_err.log
│
└── STEP5_stats/                      # All statistical analysis runs (timestamped)
    ├── STEP5_stats_20260511_175522/  # ← Latest (active)
    │   ├── manu_graphs/              #   Final publication figures
    │   ├── zone_main_cells_aggregated/
    │   ├── zone_sec_cells_aggregated/
    │   ├── *.csv files               #   Statistical tables & CLD
    │   └── analysis_report_choice_exp.docx
    ├── STEP5_stats_20260511_172318/
    ├── STEP5_stats_20260511_165402/
    └── ... (previous runs, archived)
```

## Quick Start

**Run full pipeline:**
```r
source("scripts/00_main/00_master_pipeline_choice_exp.R")
```

**Regenerate figures only (faster):**
```r
source("scripts/00_main/run_manu_graphs_only.R")
```

## Key Scripts

| Script | Purpose | Input | Output |
|--------|---------|-------|--------|
| `STEP2` | Frame-level video analysis | Raw video files | Spatial data per frame |
| `STEP3` | Zone classification | STEP2 output | Zone occupancy time series |
| `STATS` | Statistical models & main figures | STEP3/4 processed data | `STEP5_stats/` folder with figures, tables, Word report |
| `GRAPHS` | Publication-ready figure generation | Statistical results | PNG/PDF figures |

## Output Locations

- **Latest results folder:** `STEP5_stats/STEP5_stats_20260511_175522/` (see `run_manu_graphs_only.R` for current path)
  - Figures: `manu_graphs/` (Figure 7, 8, 9, S7, etc.)
  - Statistics: `*.csv` files per analysis (CLD tables, contrasts, etc.)
  - Word report: `analysis_report_choice_exp.docx`
- **Archived runs:** Previous timestamped folders in `STEP5_stats/` (kept for version history)

## Notes

- `STEP5_stats/` at root is in active use by the pipeline; will be moved to `output/` when not locked.
- All `.R` scripts in `scripts/` can be sourced independently if dependencies are satisfied.
- Zone reference CSVs define spatial layout per tank treatment.
